module CubicDP

using ..WindowSearch: FlatWindow, exact_numerator
using ..StructureDP: LocalFactor, audit_factors, elimination_decomposition,
                     validate_decomposition, bag_order_dp

export cubic_factors, solve_cubic_dp

const Order = Tuple{Vararg{Int}}
const IntegerVector = Dict{Int,Tuple{BigInt,BigInt}}
const VALID_PAIR_COLUMNS = (1,2)
# Zero-based comparison codes 2 and 5 are impossible for a total order.
const VALID_TRIPLE_COLUMNS = (1,2,4,5,7,8)

struct AdmissionRefusal <: Exception
    status::Symbol
    reason::String
end

mutable struct Budget
    max_factor_entries::Int
    max_exact_vector_entries::Int
    max_bits::Int
    max_work::Int
    max_bytes::Int
    guard::Any
    work::Int
    vector_entries::Int
    factor_entries::Int
    dot_products::Int
    vector_bytes::BigInt
    estimated_peak_bytes::BigInt
    cost_bit_bound::Int
end

_cpu() = ccall(:clock,Int64,())/1.0e6
_bits(x::Integer) = iszero(x) ? 0 : ndigits(abs(x);base=2)
_ceil_log2(n::Integer) = n<=1 ? 0 : ndigits(BigInt(n)-1;base=2)
_refuse(status,reason) = throw(AdmissionRefusal(status,reason))
_guard(b::Budget) = b.guard===nothing ? nothing : b.guard()

function _work!(b::Budget,amount::Integer=1)
    amount>=0 || error("negative work charge")
    amount<=b.max_work-b.work || _refuse(:work_limit,"exact factor-construction work cap exceeded")
    previous=b.work;b.work+=Int(amount)
    (previous÷1024 != b.work÷1024) && _guard(b)
    nothing
end

function _memory!(b::Budget,extra::Integer=0)
    # Raw integer tables plus exact rational conversion/reduction and scratch.
    per_entry=256+4cld(max(1,b.cost_bit_bound),8)
    estimate=b.vector_bytes+BigInt(b.factor_entries)*per_entry+BigInt(extra)
    b.estimated_peak_bytes=max(b.estimated_peak_bytes,estimate)
    estimate<=b.max_bytes || _refuse(:memory_limit,"estimated exact-vector/factor storage exceeds cap")
end

function _checked(value::BigInt,b::Budget)
    _bits(value)<=b.max_bits || _refuse(:bit_limit,"exact integer cost exceeds bit cap")
    value
end

function _valid_columns(scope)
    length(scope)==2 ? VALID_PAIR_COLUMNS : VALID_TRIPLE_COLUMNS
end

function _validate_flat(flat::FlatWindow,b::Budget)
    k=length(flat.labels);d=length(flat.destinations)
    1<=k<=9 || throw(ArgumentError("cubic DP window width must be 1:9"))
    length(flat.baseline)==d || throw(ArgumentError("baseline response dimension differs"))
    supplied=big(2)*length(flat.baseline)
    supplied<=b.max_exact_vector_entries || _refuse(:vector_entry_limit,"supplied scalar vector-entry count exceeds cap")
    b.vector_bytes+=16big(k+d)+128big(length(flat.factors))
    _memory!(b)
    _work!(b,1+k+d)
    length(unique(flat.labels))==k || throw(ArgumentError("duplicate global window labels"))
    length(unique(flat.destinations))==d || throw(ArgumentError("duplicate response destinations"))
    # Admit the supplied representation before scanning or converting entries.
    for factor in flat.factors
        _work!(b,1+length(factor.scope)+length(factor.destinations))
        supplied+=2big(length(factor.values))
        supplied<=b.max_exact_vector_entries || _refuse(:vector_entry_limit,"supplied scalar vector-entry count exceeds cap")
        b.vector_bytes+=16big(length(factor.scope)+length(factor.destinations))
        _memory!(b)
        2<=length(factor.scope)<=3 || throw(ArgumentError("response factor arity must be two or three"))
        length(unique(factor.scope))==length(factor.scope) || throw(ArgumentError("duplicate response factor label"))
        all(i->1<=i<=k,factor.scope) || throw(ArgumentError("response factor uses nonlocal labels"))
        length(unique(factor.destinations))==length(factor.destinations) || throw(ArgumentError("duplicate factor destination"))
        all(i->1<=i<=d,factor.destinations) || throw(ArgumentError("factor destination outside baseline support"))
        size(factor.values)==(length(factor.destinations),length(factor.scope)==2 ? 2 : 8) ||
            throw(ArgumentError("invalid response factor matrix shape"))
    end
    for value in flat.baseline
        _work!(b)
        isfinite(value) || throw(ArgumentError("nonfinite frozen baseline"))
    end
    for factor in flat.factors,value in factor.values
        _work!(b)
        isfinite(value) || throw(ArgumentError("nonfinite frozen factor value"))
    end
    k,d
end

function _common_denominator(flat,b::Budget)
    denominator=big(1)
    function include_component(component)
        _work!(b)
        q=Rational{BigInt}(component)
        denominator=max(denominator,Base.denominator(q))
        _bits(denominator)<=b.max_bits || _refuse(:bit_limit,"common dyadic denominator exceeds bit cap")
    end
    for value in flat.baseline
        include_component(real(value));include_component(imag(value))
    end
    # Only the two/six realizable local-order variants enter the objective.
    for factor in flat.factors,column in _valid_columns(factor.scope),row in axes(factor.values,1)
        value=factor.values[row,column]
        include_component(real(value));include_component(imag(value))
    end
    denominator
end

function _integer_component(value::Float64,denominator::BigInt,b::Budget)
    q=Rational{BigInt}(value)
    quotient,remainder=divrem(denominator,Base.denominator(q))
    iszero(remainder) || error("Float64 denominator is not a divisor of the common dyadic scale")
    _checked(numerator(q)*quotient,b)
end

function _scaled_component_bits(value::Float64,denominator::BigInt)
    iszero(value) && return 0
    q=Rational{BigInt}(value)
    shift=_bits(denominator)-_bits(Base.denominator(q))
    shift>=0 || error("component denominator exceeds the admitted dyadic scale")
    _bits(numerator(q))+shift
end

function _integer_vector(destinations,values,denominator,b::Budget)
    result=IntegerVector()
    for (destination,value) in zip(destinations,values)
        _work!(b)
        iszero(value) && continue
        b.vector_entries+2<=b.max_exact_vector_entries || _refuse(:vector_entry_limit,"retained integer vector-entry cap exceeded")
        # Check conservative integer storage before allocating the BigInts.
        re_bits=_scaled_component_bits(real(value),denominator)
        im_bits=_scaled_component_bits(imag(value),denominator)
        bytes=128+cld(max(1,re_bits),8)+cld(max(1,im_bits),8)
        b.vector_entries+=2;b.vector_bytes+=bytes;_memory!(b)
        result[destination]=(_integer_component(real(value),denominator,b),
                             _integer_component(imag(value),denominator,b))
    end
    result
end

function _dot(left::IntegerVector,right::IntegerVector,b::Budget)
    b.dot_products+=1
    a,c=length(left)<=length(right) ? (left,right) : (right,left)
    _work!(b,1+4length(a))
    total=big(0)
    for (destination,(re,im)) in a
        other=get(c,destination,nothing)
        other===nothing && continue
        total+=re*other[1]+im*other[2]
    end
    _checked(total,b)
end

function _orders(labels)
    isempty(labels) && return [()]
    result=Order[]
    for label in labels,rest in _orders([x for x in labels if x!=label])
        push!(result,(label,rest...))
    end
    result
end

function _column(scope,order)
    a,b=scope[1],scope[2]
    at(x)=findfirst(==(x),order)
    code=Int(at(a)>at(b))
    length(scope)==3 && (code+=2Int(at(a)>at(scope[3]))+4Int(at(b)>at(scope[3])))
    code+1
end

function _add_term!(grouped,scope,value,b::Budget)
    key=Tuple(sort!(unique(collect(scope))))
    if !haskey(grouped,key)
        required=Int(factorial(big(length(key))))
        required<=b.max_factor_entries-b.factor_entries || _refuse(:factor_entry_limit,"scalar factor-entry cap exceeded")
        b.factor_entries+=required;_memory!(b)
        _work!(b,required)
        grouped[key]=Dict{Order,BigInt}(order=>big(0) for order in _orders(collect(key)))
    end
    table=grouped[key]
    _work!(b,length(table))
    for order in keys(table)
        table[order]=_checked(table[order]+value(order),b)
    end
end

function _stats(b::Budget;integerization_CPU_seconds=0.0,integerization_wall_seconds=0.0,
        factor_CPU_seconds=0.0,factor_wall_seconds=0.0,kwargs...)
    (;work=b.work,exact_vector_entries=b.vector_entries,
      scalar_entries_before_reduction=b.factor_entries,gram_dot_products=b.dot_products,
      estimated_peak_bytes=b.estimated_peak_bytes,cost_bit_bound=b.cost_bit_bound,
      integerization_CPU_seconds,integerization_wall_seconds,
      factor_CPU_seconds,factor_wall_seconds,kwargs...)
end

"""Expand EXACT frozen complex response DELTAS into signed scalar factors.

The supplied objective is Q=||R0+sum_S delta_S||^2. All factors use local labels
1:k. If D is the common dyadic vector denominator, factor costs sum to integer
N and Q=N/D^2. This does not re-evaluate transport, projection, or cubic columns.
Only valid pair/triple comparison columns are response variants.
"""
function cubic_factors(flat::FlatWindow;guard=nothing,max_factor_entries::Int=200_000,
        max_exact_vector_entries::Int=1_000_000,max_bits::Int=4096,
        max_work::Int=100_000_000,max_bytes::Int=128*1024*1024)
    all(x->x>=1,(max_factor_entries,max_exact_vector_entries,max_bits,max_work,max_bytes)) ||
        throw(ArgumentError("all factor admission caps must be positive"))
    budget=Budget(max_factor_entries,max_exact_vector_entries,max_bits,max_work,max_bytes,
        guard,0,0,0,0,big(0),big(0),max_bits)
    started_cpu=_cpu();started_wall=time_ns();phase=:integerization
    integer_cpu=0.0;integer_wall=0.0;D=nothing
    try
        _guard(budget)
        k,d=_validate_flat(flat,budget)
        D=_common_denominator(flat,budget)
        2*_bits(D)-1<=max_bits || _refuse(:bit_limit,"squared dyadic denominator exceeds bit cap")
        denominator=D*D
        baseline=_integer_vector(1:d,flat.baseline,D,budget)
        response_factors=NamedTuple[]
        for factor in flat.factors
            variants=Dict{Int,IntegerVector}()
            for column in _valid_columns(factor.scope)
                variants[column]=_integer_vector(factor.destinations,view(factor.values,:,column),D,budget)
            end
            all(isempty,values(variants)) && continue
            push!(response_factors,(scope=copy(factor.scope),variants=variants))
        end
        max_integer_bits=1
        for pair in values(baseline),component in pair
            max_integer_bits=max(max_integer_bits,_bits(component))
        end
        for factor in response_factors,vector in values(factor.variants),pair in values(vector),component in pair
            max_integer_bits=max(max_integer_bits,_bits(component))
        end
        count=length(response_factors)
        budget.cost_bit_bound=2max_integer_bits+_ceil_log2(max(1,2d))+
            _ceil_log2(max(big(1),2big(count)^2+2count+1))+4
        budget.cost_bit_bound<=max_bits || _refuse(:bit_limit,"conservative Gram/factor bit bound exceeds cap")
        _memory!(budget)
        integer_cpu=_cpu()-started_cpu;integer_wall=(time_ns()-started_wall)/1e9
        phase=:factor;factor_cpu_start=_cpu();factor_wall_start=time_ns()
        grouped=Dict{Order,Dict{Order,BigInt}}()
        baseline_norm=_dot(baseline,baseline,budget)
        iszero(baseline_norm) || _add_term!(grouped,Int[],order->baseline_norm,budget)
        for factor in response_factors
            _guard(budget)
            h=Dict{Int,BigInt}()
            for (column,vector) in factor.variants
                h[column]=_checked(2*_dot(baseline,vector,budget)+_dot(vector,vector,budget),budget)
            end
            all(iszero,values(h)) && continue
            _add_term!(grouped,factor.scope,order->h[_column(factor.scope,order)],budget)
        end
        for i in eachindex(response_factors),j in 1:i-1
            _guard(budget)
            left,right=response_factors[i],response_factors[j]
            # At most 6x6 scalar cached dots; no global dense Gram tensor.
            _memory!(budget,36*(96+cld(budget.cost_bit_bound,8)))
            gram=Dict{Tuple{Int,Int},BigInt}()
            for (a,va) in left.variants,(c,vc) in right.variants
                gram[a,c]=_checked(2*_dot(va,vc,budget),budget)
            end
            all(iszero,values(gram)) && continue
            scope=vcat(left.scope,right.scope)
            _add_term!(grouped,scope,
                order->gram[_column(left.scope,order),_column(right.scope,order)],budget)
        end
        factors=LocalFactor[];constant=big(0);zero_factors=0
        for scope in sort!(collect(keys(grouped));by=x->(length(x),x))
            _guard(budget)
            table=grouped[scope]
            if all(iszero,values(table))
                zero_factors+=1;continue
            end
            _work!(budget,1+length(table)*max(1,length(scope)))
            factor=LocalFactor(collect(scope),table;max_entries=720)
            if isempty(factor.scope)
                cost=only(values(factor.values))
                Base.denominator(cost)==1 || error("integer factor acquired a denominator")
                constant=_checked(constant+numerator(cost),budget)
            else
                push!(factors,factor)
            end
        end
        iszero(constant) || push!(factors,LocalFactor(Int[],Dict(()=>constant)))
        factor_cpu=_cpu()-factor_cpu_start;factor_wall=(time_ns()-factor_wall_start)/1e9
        stats=_stats(budget;integerization_CPU_seconds=integer_cpu,
            integerization_wall_seconds=integer_wall,factor_CPU_seconds=factor_cpu,
            factor_wall_seconds=factor_wall,input_response_factors=length(flat.factors),
            nonzero_response_factors=count,retained_factors=length(factors),
            retained_factor_entries=sum(length(f.values) for f in factors;init=0),
            exact_zero_groups_removed=zero_factors,constant_factor_value=constant,
            max_retained_scope=maximum((length(f.scope) for f in factors);init=0),
            local_labels=collect(1:k),valid_pair_variants=2,valid_triple_variants=6)
        return (status=:complete,complete=true,reason="",factors=factors,
            vector_denominator=D,denominator=denominator,stats=stats)
    catch exception
        exception isa AdmissionRefusal || rethrow()
        elapsed_cpu=_cpu()-started_cpu;elapsed_wall=(time_ns()-started_wall)/1e9
        stats=_stats(budget;integerization_CPU_seconds=phase==:integerization ? elapsed_cpu : integer_cpu,
            integerization_wall_seconds=phase==:integerization ? elapsed_wall : integer_wall,
            factor_CPU_seconds=phase==:factor ? elapsed_cpu-integer_cpu : 0.0,
            factor_wall_seconds=phase==:factor ? elapsed_wall-integer_wall : 0.0,
            stopped_phase=phase)
        return (status=exception.status,complete=false,reason=exception.reason,
            factors=LocalFactor[],vector_denominator=D,denominator=nothing,stats=stats)
    end
end

"""Exactly optimize the supplied cubic frozen-table Q with bag-order DP.

score is Rational{BigInt}, equal to WindowSearch.exact_numerator(flat,local_order).
order maps local labels through flat.labels exactly once. Division by any fixed
Strang r^4 belongs to the caller and is not performed by this adapter.
"""
function solve_cubic_dp(flat::FlatWindow;guard=nothing,max_factor_entries::Int=200_000,
        max_exact_vector_entries::Int=1_000_000,max_bits::Int=4096,
        max_work::Int=100_000_000,max_bytes::Int=128*1024*1024,
        max_states::Int=100_000,max_bag_states::Int=50_000,
        max_dp_work::Int=100_000_000,max_dp_bytes::Int=128*1024*1024)
    data=cubic_factors(flat;guard=guard,max_factor_entries=max_factor_entries,
        max_exact_vector_entries=max_exact_vector_entries,max_bits=max_bits,
        max_work=max_work,max_bytes=max_bytes)
    failure=(numerator=nothing,score=nothing,order=Int[],local_order=Int[])
    !data.complete && return merge(data,failure)
    labels=collect(1:length(flat.labels));graph_cpu=_cpu();graph_wall=time_ns()
    audit=audit_factors(labels,data.factors;max_entries=max_factor_entries,guard=guard)
    if !audit.complete
        return merge(data,failure,(status=:graph_incomplete,complete=false,
            reason="exact scalar-factor graph audit did not complete",audit=audit))
    end
    decomposition=elimination_decomposition(labels,audit.edges;guard=guard)
    if !decomposition.complete
        return merge(data,failure,(status=:decomposition_incomplete,complete=false,
            reason="graph decomposition did not complete",audit=audit,decomposition=decomposition))
    end
    validate_decomposition(labels,data.factors,decomposition.bags,decomposition.tree_edges)
    graph_CPU_seconds=_cpu()-graph_cpu;graph_wall_seconds=(time_ns()-graph_wall)/1e9
    solve_cpu=_cpu();solve_wall=time_ns()
    solved=bag_order_dp(labels,data.factors,decomposition.bags,decomposition.tree_edges;
        max_states=max_states,max_bag_states=max_bag_states,max_work=max_dp_work,
        max_factor_entries=max_factor_entries,max_bytes=max_dp_bytes,
        max_cost_bits=max_bits,guard=guard)
    solve_CPU_seconds=_cpu()-solve_cpu;solve_wall_seconds=(time_ns()-solve_wall)/1e9
    stats=merge(data.stats,(graph_decomposition_CPU_seconds=graph_CPU_seconds,
        graph_decomposition_wall_seconds=graph_wall_seconds,solve_CPU_seconds=solve_CPU_seconds,
        solve_wall_seconds=solve_wall_seconds,graph_edges=length(audit.edges),
        used_width=decomposition.width_upper_bound,bag_sizes=length.(decomposition.bags)))
    if !solved.complete
        return merge(data,failure,(status=solved.status,complete=false,reason=solved.reason,
            stats=stats,dp=solved))
    end
    witness_cpu=_cpu();witness_wall=time_ns()
    Base.denominator(solved.cost)==1 || error("cubic integer DP cost is nonintegral")
    N=BigInt(numerator(solved.cost));score=N//data.denominator
    N>=0 || error("squared response DP returned a negative exact objective")
    score==exact_numerator(flat,solved.order) || error("DP/scalar-factor frozen response identity failed")
    global_order=flat.labels[solved.order]
    stats=merge(stats,(witness_CPU_seconds=_cpu()-witness_cpu,
        witness_wall_seconds=(time_ns()-witness_wall)/1e9,dp_details=solved.details,
        witness_exact_score_matches=true))
    return (status=:complete,complete=true,reason="",numerator=N,
        vector_denominator=data.vector_denominator,denominator=data.denominator,
        score=score,order=global_order,local_order=solved.order,stats=stats,dp=solved,
        scope="Exact supplied frozen complex response DELTA tables; Q only, fixed r^4 external")
end

end # module
