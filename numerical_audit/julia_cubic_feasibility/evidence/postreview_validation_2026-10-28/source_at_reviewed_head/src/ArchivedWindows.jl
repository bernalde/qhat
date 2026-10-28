module ArchivedWindows

using TOML
using ..StructureDP

export ArchiveData, load_archive, archive_factors, archived_integer_score,
       archive_full_order, enumeration, enumerate_archive, solve_archive_dp,
       validate_archive_solution, run_archive

_guard(g)=g===nothing ? nothing : g()
_cpu()=ccall(:clock,Int64,())/1.0e6 # Same Linux process-clock convention as StudyIO.
_big(x::Integer)=BigInt(x)
_big(x::AbstractString)=parse(BigInt,x)
_int(x::Integer)=Int(x)
_int(x::AbstractString)=parse(Int,x)
_row(x::AbstractString)=split(x,',')
_row(x)=x
_get(d,k,default=nothing)=haskey(d,k) ? d[k] : (haskey(d,Symbol(k)) ? d[Symbol(k)] : default)

"""The archived FIRST-ORDER JD window, retaining original zero-based labels.

Integer weights encode the archived supplied binary64 sinc-squared values.
They are never square-rooted or recomputed through a cubic-response kernel.
"""
struct ArchiveData
    case::String
    width::Int
    window::Vector{Int}
    warm_order::Vector{Int}
    outside_prefix::Vector{Int}
    outside_suffix::Vector{Int}
    pairs::Vector{Tuple{Int,Int,Int,BigInt}}
    integer_weights::Dict{Int,BigInt}
    objective_denominator::BigInt
    fixed_outside::Dict{Int,BigInt}
    internal_pairs::Dict{Tuple{Int,Int},Dict{Int,BigInt}}
    expected_minimum::Union{Nothing,BigInt}
    expected_minimizers::Union{Nothing,Int}
    expected_window::Vector{Int}
    metadata::Dict{String,Any}
end

function ArchiveData(raw::AbstractDict; weight_multiplier::Integer=1,
        max_width::Int=9,max_rows::Int=10_000,max_destinations::Int=256)
    weight_multiplier>0 || throw(ArgumentError("positive integer weight multiplier required"))
    window=_int.(_get(raw,"window"));warm=_int.(_get(raw,"warm_order"))
    prefix=_int.(_get(raw,"outside_prefix"));suffix=_int.(_get(raw,"outside_suffix"))
    width=_int(_get(raw,"width"))
    width==length(window) && 1<=width<=max_width || throw(ArgumentError("window width bound"))
    length(unique(window))==width || throw(ArgumentError("duplicate window label"))
    warm==[prefix;window;suffix] || throw(ArgumentError("outside/window inconsistency"))
    sort(warm)==collect(0:length(warm)-1) || throw(ArgumentError("archive full order must use original zero-based labels"))
    denominator=_big(_get(raw,"objective_denominator"))
    denominator>0 || throw(ArgumentError("positive objective denominator required"))
    pair_rows=_get(raw,"pairs");weight_rows=_get(raw,"integer_weights")
    length(pair_rows)<=max_rows || throw(ArgumentError("pair row cap exceeded"))
    length(weight_rows)<=max_destinations || throw(ArgumentError("destination cap exceeded"))
    weights=Dict{Int,BigInt}()
    for entry in weight_rows
        row=_row(entry);length(row)==2 || throw(ArgumentError("weight row arity"))
        destination=_int(row[1]);weight=_big(row[2])*weight_multiplier
        weight>=0 || throw(ArgumentError("negative archived weight"))
        haskey(weights,destination) && throw(ArgumentError("duplicate destination weight"))
        weights[destination]=weight
    end
    pairs=Tuple{Int,Int,Int,BigInt}[];destinations=Set{Int}()
    positions=Dict(v=>i for (i,v) in enumerate(warm));selected=Set(window)
    fixed=Dict{Int,BigInt}();internals=Dict{Tuple{Int,Int},Dict{Int,BigInt}}()
    for entry in pair_rows
        row=_row(entry);length(row)==4 || throw(ArgumentError("pair row arity"))
        i,j,d=_int(row[1]),_int(row[2]),_int(row[3]);b=_big(row[4])
        i<j && haskey(positions,i) && haskey(positions,j) || throw(ArgumentError("invalid canonical pair labels"))
        haskey(weights,d) || throw(ArgumentError("pair has no destination weight"))
        push!(pairs,(i,j,d,b));push!(destinations,d)
        if i in selected && j in selected
            vector=get!(internals,(i,j),Dict{Int,BigInt}())
            vector[d]=get(vector,d,big(0))+b
        else
            fixed[d]=get(fixed,d,big(0))+(positions[i]<positions[j] ? b : -b)
        end
    end
    Set(keys(weights))==destinations || throw(ArgumentError("destination weight coverage differs"))
    for vector in values(internals)
        filter!(kv->!iszero(last(kv)),vector)
    end
    filter!(kv->!isempty(last(kv)),internals)
    expected=_get(raw,"expected_minimum")
    expected=expected===nothing ? nothing : _big(expected)*weight_multiplier
    ties=_get(raw,"expected_minimizers")
    ties=ties===nothing ? nothing : _int(ties)
    expected_window=_int.(_get(raw,"expected_window",Int[]))
    isempty(expected_window) || sort(expected_window)==sort(window) || throw(ArgumentError("expected window is not the frozen label set"))
    metadata=Dict{String,Any}()
    for key in ("domain_sha256","certificate_sha256","archived_enumeration_cpu_seconds")
        value=_get(raw,key);value===nothing || (metadata[key]=value)
    end
    metadata["weight_multiplier"]=string(weight_multiplier)
    return ArchiveData(string(_get(raw,"case","archive")),width,window,warm,prefix,suffix,
        pairs,weights,denominator,fixed,internals,expected,ties,expected_window,metadata)
end

load_archive(path::AbstractString;kwargs...)=ArchiveData(TOML.parsefile(path);kwargs...)

function archive_full_order(data::ArchiveData, order)
    proposed=Int[x for x in order]
    if length(proposed)==data.width
        sort(proposed)==sort(data.window) || throw(ArgumentError("not a permutation of the frozen window"))
        return [data.outside_prefix;proposed;data.outside_suffix]
    end
    length(proposed)==length(data.warm_order) || throw(ArgumentError("invalid full-order length"))
    start=length(data.outside_prefix)+1;stop=start+data.width-1
    proposed[1:start-1]==data.outside_prefix && proposed[stop+1:end]==data.outside_suffix ||
        throw(ArgumentError("outside order or positions changed"))
    sort(proposed[start:stop])==sort(data.window) || throw(ArgumentError("window labels changed"))
    return proposed
end

"""Direct full-pair numerator, independent of the scalar-factor expansion."""
function archived_integer_score(data::ArchiveData, order)
    full=archive_full_order(data,order);positions=Dict(v=>i for (i,v) in enumerate(full))
    response=Dict{Int,BigInt}()
    for (i,j,d,b) in data.pairs
        response[d]=get(response,d,big(0))+(positions[i]<positions[j] ? b : -b)
    end
    result=big(0)
    for (d,weight) in data.integer_weights
        result+=weight*get(response,d,big(0))^2
    end
    return result
end

function _orders(scope)
    isempty(scope) && return [()]
    out=Tuple[]
    for v in scope,rest in _orders([x for x in scope if x!=v])
        push!(out,(v,rest...))
    end
    return out
end

_sign(order,e)=findfirst(==(e[1]),order)<findfirst(==(e[2]),order) ? 1 : -1
function _dot(data,a,b)
    result=big(0)
    for (d,value) in a
        result+=data.integer_weights[d]*value*get(b,d,big(0))
    end
    return result
end

"""Expand the archived coherent square into exact integer factors of arity <=4.

For b_d=c_d+sum_e sign_e v_ed, retain c², v_e², 2c v_e,
and 2v_e v_f. No fixed-outside interference or signed cross term is dropped.
"""
function archive_factors(data::ArchiveData;max_internal_pairs::Int=36,
        max_cross_terms::Int=630,max_factors::Int=256,max_entries::Int=10_000,
        guard=nothing)
    edges=sort!(collect(keys(data.internal_pairs)))
    length(edges)<=max_internal_pairs || throw(ArgumentError("internal pair cap exceeded"))
    length(edges)*(length(edges)-1)÷2<=max_cross_terms || throw(ArgumentError("factor cross-term cap exceeded"))
    grouped=Dict{Tuple,Dict{Tuple,BigInt}}();entry_count=0
    function add_term(scope,coefficient,value)
        iszero(coefficient) && return
        key=Tuple(sort!(unique(collect(scope))))
        if !haskey(grouped,key)
            length(grouped)<max_factors || throw(ArgumentError("factor count cap exceeded"))
            entries=Int(factorial(big(length(key))))
            entry_count+entries<=max_entries || throw(ArgumentError("factor entry cap exceeded"))
            table=Dict{Tuple,BigInt}()
            for order in _orders(collect(key))
                table[order]=big(0)
            end
            grouped[key]=table;entry_count+=entries
        end
        for order in keys(grouped[key])
            grouped[key][order]+=coefficient*value(order)
        end
    end
    constant=_dot(data,data.fixed_outside,data.fixed_outside)
    for e in edges
        _guard(guard);vector=data.internal_pairs[e]
        constant+=_dot(data,vector,vector)
        coefficient=2*_dot(data,data.fixed_outside,vector)
        add_term(e,coefficient,order->_sign(order,e))
    end
    add_term(Int[],constant,order->1)
    for i in eachindex(edges),j in 1:i-1
        _guard(guard);e,f=edges[i],edges[j]
        coefficient=2*_dot(data,data.internal_pairs[e],data.internal_pairs[f])
        add_term(vcat(collect(e),collect(f)),coefficient,order->_sign(order,e)*_sign(order,f))
    end
    factors=LocalFactor[]
    for scope in sort!(collect(keys(grouped));by=x->(length(x),x))
        push!(factors,LocalFactor(collect(scope),grouped[scope]))
    end
    return factors
end

function validate_archive_solution(data::ArchiveData,order,cost)
    full=archive_full_order(data,order);actual=archived_integer_score(data,full)
    integral=cost isa Integer || (cost isa Rational && denominator(cost)==1)
    encoded=integral ? BigInt(cost isa Rational ? numerator(cost) : cost) : nothing
    expected=data.expected_minimum
    return (valid=integral && encoded==actual && (expected===nothing || actual==expected),
            score_matches=integral && encoded==actual,
            expected_matches=expected===nothing ? nothing : actual==expected,
            numerator=actual,objective_denominator=data.objective_denominator,
            JD=actual//data.objective_denominator,full_order=full,
            scope="archived first-order JD; fixed outside; supplied exact weight ratios")
end

"""Fresh independent DFS baseline with exact integers and Lehmer coverage.

Only destinations appearing in a nonzero internal pair vector can vary within
the contiguous window. All other destinations are folded into an exact constant
sum(weight * fixed_outside^2). Active amplitudes retain their FULL outside term.
"""
function enumeration(data::ArchiveData;guard=nothing,max_candidates::Int=400_000,
        max_delta_entries::Int=1_000_000,max_bytes::Int=64*1024*1024)
    setup_cpu=_cpu();setup_wall=time_ns();m=data.width;domain=factorial(big(m))
    domain<=max_candidates || return (status=:candidate_limit,complete=false,evaluated=0,
        cost=nothing,order=Int[],reason="enumeration domain exceeds cap")
    labels=sort(data.window)
    active_destinations=Set{Int}()
    for vector in values(data.internal_pairs),destination in keys(vector)
        push!(active_destinations,destination)
    end
    destinations=sort!(collect(active_destinations));dcount=length(destinations)
    outside_constant=big(0)
    for (destination,weight) in data.integer_weights
        destination in active_destinations && continue
        outside_constant+=weight*get(data.fixed_outside,destination,big(0))^2
    end
    delta_entries=big(1<<m)*m*dcount
    amplitude_bits=maximum((ndigits(abs(b);base=2) for (_,_,_,b) in data.pairs);init=1)+
                   ndigits(max(1,length(data.pairs));base=2)+1
    estimated_bytes=(delta_entries+big(m)*m*dcount)*(96+cld(amplitude_bits,8))+cld(domain,8)
    if delta_entries>max_delta_entries || estimated_bytes>max_bytes
        return (status=:memory_limit,complete=false,evaluated=0,cost=nothing,order=Int[],
                reason="enumeration delta storage exceeds cap",estimated_bytes=estimated_bytes)
    end
    _guard(guard)
    positions=Dict(v=>i for (i,v) in enumerate(labels));destination_index=Dict(d=>i for (i,d) in enumerate(destinations))
    coupling=zeros(BigInt,m,m,dcount)
    for ((i,j),vector) in data.internal_pairs,(d,b) in vector
        a,bindex=positions[i],positions[j];k=destination_index[d]
        coupling[a,bindex,k]+=b;coupling[bindex,a,k]-=b
    end
    delta=zeros(BigInt,1<<m,m,dcount)
    for mask in 1:(1<<m)-1
        mask%64==0 && _guard(guard)
        earlier=trailing_zeros(mask)+1;previous=mask&(mask-1)
        for candidate in 1:m
            mask & (1<<(candidate-1)) == 0 || continue
            for k in 1:dcount
                delta[mask+1,candidate,k]=delta[previous+1,candidate,k]+coupling[earlier,candidate,k]
            end
        end
    end
    response=BigInt[get(data.fixed_outside,d,big(0)) for d in destinations]
    weights=BigInt[data.integer_weights[d] for d in destinations]
    coverage=falses(Int(domain));prefix=zeros(Int,m);facts=[Int(factorial(big(i))) for i in 0:m]
    best=nothing;winner=Int[];ties=0;evaluated=0
    setup_CPU_seconds=_cpu()-setup_cpu;setup_wall_seconds=(time_ns()-setup_wall)/1.0e9
    solve_cpu=_cpu();solve_wall=time_ns()
    function visit(depth,mask)
        if depth>m
            evaluated%256==0 && _guard(guard)
            rank=0
            for i in 1:m-1
                smaller=count(j->prefix[j]<prefix[i],i+1:m)
                rank+=smaller*facts[m-i+1]
            end
            coverage[rank+1] && error("duplicate Lehmer rank")
            coverage[rank+1]=true;evaluated+=1
            value=outside_constant
            for k in 1:dcount
                value+=weights[k]*response[k]^2
            end
            if best===nothing || value<best
                best=value;winner=Int[labels[i] for i in prefix];ties=1
            elseif value==best
                ties+=1
            end
            return
        end
        for candidate in 1:m
            bit=1<<(candidate-1);mask&bit==0 || continue
            prefix[depth]=candidate
            for k in 1:dcount
                response[k]+=delta[mask+1,candidate,k]
            end
            visit(depth+1,mask|bit)
            for k in 1:dcount
                response[k]-=delta[mask+1,candidate,k]
            end
        end
    end
    visit(1,0);_guard(guard)
    complete=evaluated==domain && all(coverage)
    complete || error("incomplete DFS coverage")
    solve_CPU_seconds=_cpu()-solve_cpu;solve_wall_seconds=(time_ns()-solve_wall)/1.0e9
    validation_cpu=_cpu();validation_wall=time_ns()
    validation=validate_archive_solution(data,winner,best)
    validation.score_matches || error("DFS/full-pair numerator mismatch")
    validation_CPU_seconds=_cpu()-validation_cpu
    validation_wall_seconds=(time_ns()-validation_wall)/1.0e9
    return (status=:complete,complete=true,evaluated=evaluated,cost=best,order=winner,
        minimizers=ties,coverage_lehmer_count=count(coverage),coverage_complete=all(coverage),
        validation=validation,setup_CPU_seconds=setup_CPU_seconds,
        setup_wall_seconds=setup_wall_seconds,solve_CPU_seconds=solve_CPU_seconds,
        solve_wall_seconds=solve_wall_seconds,
        validation_CPU_seconds=validation_CPU_seconds,validation_wall_seconds=validation_wall_seconds,
        estimated_bytes=estimated_bytes,active_destination_count=dcount,
        constant_destination_count=length(data.integer_weights)-dcount,
        original_destination_count=length(data.integer_weights),outside_constant=outside_constant,
        expected_minimizers_match=data.expected_minimizers===nothing ? nothing : ties==data.expected_minimizers)
end

enumerate_archive(args...;kwargs...)=enumeration(args...;kwargs...)

function solve_archive_dp(data::ArchiveData;guard=nothing,factors=nothing,
        max_states::Int=100_000,max_bag_states::Int=50_000,
        max_bytes::Int=128*1024*1024,kwargs...)
    setup_cpu=_cpu();setup_wall=time_ns()
    fs=factors===nothing ? archive_factors(data;guard=guard) : factors
    factor_count=length(fs);factor_entries=sum(length(f.values) for f in fs)
    factor_construction_CPU_seconds=_cpu()-setup_cpu
    factor_construction_wall_seconds=(time_ns()-setup_wall)/1.0e9
    factor_metadata=(factor_count=factor_count,factor_entries=factor_entries,
        factor_construction_CPU_seconds=factor_construction_CPU_seconds,
        factor_construction_wall_seconds=factor_construction_wall_seconds,
        factor_construction_performed=factors===nothing)
    graph_cpu=_cpu();graph_wall=time_ns()
    audit=audit_factors(data.window,fs;guard=guard)
    if !audit.complete
        return merge((status=:graph_incomplete,complete=false,cost=nothing,order=Int[],audit=audit),
            factor_metadata,(graph_decomposition_CPU_seconds=_cpu()-graph_cpu,
                graph_decomposition_wall_seconds=(time_ns()-graph_wall)/1.0e9,
                setup_CPU_seconds=_cpu()-setup_cpu,setup_wall_seconds=(time_ns()-setup_wall)/1.0e9))
    end
    decomposition=elimination_decomposition(data.window,audit.edges;guard=guard)
    if !decomposition.complete
        return merge((status=:decomposition_incomplete,complete=false,cost=nothing,
            order=Int[],audit=audit,decomposition=decomposition),factor_metadata,
            (graph_decomposition_CPU_seconds=_cpu()-graph_cpu,
             graph_decomposition_wall_seconds=(time_ns()-graph_wall)/1.0e9,
             setup_CPU_seconds=_cpu()-setup_cpu,setup_wall_seconds=(time_ns()-setup_wall)/1.0e9))
    end
    validate_decomposition(data.window,fs,decomposition.bags,decomposition.tree_edges)
    graph_decomposition_CPU_seconds=_cpu()-graph_cpu
    graph_decomposition_wall_seconds=(time_ns()-graph_wall)/1.0e9
    setup_CPU_seconds=_cpu()-setup_cpu;setup_wall_seconds=(time_ns()-setup_wall)/1.0e9
    solve_cpu=_cpu();solve_wall=time_ns()
    result=bag_order_dp(data.window,fs,decomposition.bags,decomposition.tree_edges;
        guard=guard,max_states=max_states,max_bag_states=max_bag_states,max_bytes=max_bytes,kwargs...)
    solve_CPU_seconds=_cpu()-solve_cpu;solve_wall_seconds=(time_ns()-solve_wall)/1.0e9
    validation_cpu=_cpu();validation_wall=time_ns()
    validation=result.complete ? validate_archive_solution(data,result.order,result.cost) : nothing
    validation===nothing || validation.score_matches || error("DP/full-pair numerator mismatch")
    validation_CPU_seconds=_cpu()-validation_cpu
    validation_wall_seconds=(time_ns()-validation_wall)/1.0e9
    return merge(result,factor_metadata,(validation=validation,graph_edges=length(audit.edges),
        actual_factor_clique_lower_bound=audit.clique_width_lower_bound,
        heuristic_width_upper_bound=decomposition.width_upper_bound,
        bag_count=length(decomposition.bags),bag_sizes=length.(decomposition.bags),
        setup_CPU_seconds=setup_CPU_seconds,setup_wall_seconds=setup_wall_seconds,
        graph_decomposition_CPU_seconds=graph_decomposition_CPU_seconds,
        graph_decomposition_wall_seconds=graph_decomposition_wall_seconds,
        solve_CPU_seconds=solve_CPU_seconds,solve_wall_seconds=solve_wall_seconds,
        validation_CPU_seconds=validation_CPU_seconds,validation_wall_seconds=validation_wall_seconds,
        minimizer_count_claim=false))
end

function run_archive(data::ArchiveData;method::Symbol=:dp,kwargs...)
    method==:dp && return solve_archive_dp(data;kwargs...)
    method==:enumeration && return enumeration(data;kwargs...)
    throw(ArgumentError("method must be :dp or :enumeration"))
end
run_archive(path::AbstractString;weight_multiplier::Integer=1,kwargs...)=
    run_archive(load_archive(path;weight_multiplier=weight_multiplier);kwargs...)

end # module
