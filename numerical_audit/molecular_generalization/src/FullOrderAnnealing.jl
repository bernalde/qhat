module FullOrderAnnealing

# Study-local kernels. The native trainer and conditioned-window backend stay
# unchanged; only the guarded scalar/Taylor/transport conventions are reused.
using SHA, Random
include(joinpath(@__DIR__, "..", "..", "..", "analysis", "julia_window", "ResponseCore.jl"))
using .ResponseCore: Unit, NumericalResolution

export Unit, FullOrderContext, full_column, full_score!, anneal, insertion,
       metropolis, temperature, validate_settings, Limit, noop

const Column = Dict{UInt64,Float64}
const Response = Dict{UInt64,ComplexF64}
const SETTINGS_DEFAULTS = Dict{String,Any}(
    "version"=>1, "calibrated"=>false, "calibration_sha256"=>"",
    "theta_start"=>0.0, "theta_end"=>0.0, "energy_floor"=>1e-30,
    "max_evaluations"=>1_000_000, "max_proposals"=>1_000_000,
    "max_actions"=>100_000_000, "cache_bytes"=>64*1024^2,
    "native_reserve_factor"=>2.0, "native_reserve_min"=>0.05,
    "cleanup_CPU_seconds"=>0.05)
noop(kind::Symbol=:loop, amount::Int=1) = nothing
cpu_seconds() = ccall(:clock, Int64, ()) / 1.0e6

struct Limit <: Exception
    status::String
    reason::String
end
Base.showerror(io::IO, e::Limit) = print(io, e.status, ": ", e.reason)
struct SliceExpired <: Exception end

function integer(x, name, lo, hi)
    x isa Integer && !(x isa Bool) && lo<=x<=hi || throw(ArgumentError("integer domain: "*name))
    Int(x)
end
function positive(x, name)
    x isa Real && !(x isa Bool) && isfinite(x) && x>0 || throw(ArgumentError("positive finite "*name))
    Float64(x)
end
function validate_settings(raw)
    raw isa AbstractDict && Set(keys(raw))==Set(keys(SETTINGS_DEFAULTS)) ||
        throw(ArgumentError("complete annealing settings schema required"))
    raw["version"] === 1 && raw["calibrated"] === true || throw(ArgumentError("calibrated settings required"))
    raw["calibration_sha256"] isa String && occursin(r"^[0-9a-f]{64}$",raw["calibration_sha256"]) ||
        throw(ArgumentError("calibration hash required"))
    hi=positive(raw["theta_start"],"theta_start");lo=positive(raw["theta_end"],"theta_end")
    hi>lo || throw(ArgumentError("temperature endpoints must decrease"))
    raw["energy_floor"] == 1e-30 || throw(ArgumentError("frozen energy floor"))
    integer(raw["max_evaluations"],"max_evaluations",4,1_000_000)
    integer(raw["max_proposals"],"max_proposals",1,1_000_000)
    integer(raw["max_actions"],"max_actions",1,100_000_000)
    integer(raw["cache_bytes"],"cache_bytes",0,64*1024^2)
    1<=positive(raw["native_reserve_factor"],"native_reserve_factor")<=10 || throw(ArgumentError("native reserve factor cap"))
    0.01<=positive(raw["native_reserve_min"],"native_reserve_min")<=10 || throw(ArgumentError("native reserve cap"))
    0.01<=positive(raw["cleanup_CPU_seconds"],"cleanup_CPU_seconds")<=10 || throw(ArgumentError("cleanup reserve cap"))
    Dict{String,Any}(raw)
end

mutable struct FullOrderContext
    units::Vector{Unit}
    centered::Vector{Unit}
    scalar_bases::Vector{Vector{UInt64}}
    hf::UInt64
    n::Int
    budget::Int
    time::Float64
    max_entries::Int
    max_cached_actions::Int
    cache_bytes::Int
    accounted_cache_bytes::Int
    actions::Dict{Tuple{Int,UInt64},Float64}
    kernels::Dict{UInt64,ComplexF64}
    exact::Column
    commuting::Bool
    buffers::Vector{Column}
    response::Response
    stats::Dict{String,Int}
end

"""Deterministic GF(2) row basis spanning every stored centered parity mask.

Input masks and coefficients are never reordered or rewritten. At most n<=14
masks are retained per unit: <=112 KiB of meaningful UInt64 masks over 1024
units, plus bounded vector capacity/containers. This static storage is separate from
the existing 1024-byte-per-cache-entry accounting and does not enlarge its cap.
"""
function _parity_basis(parities::AbstractVector{UInt64}, n::Int)
    1<=n<=14 || throw(ArgumentError("signature register cap"))
    pivots=zeros(UInt64,n)
    for mask in parities
        ResponseCore.check_budget()
        mask < UInt64(1)<<n || throw(ArgumentError("signature mask outside register"))
        row=mask
        # At most fourteen XOR/pivot steps between input-mask budget checks.
        for bit in n-1:-1:0
            iszero(row & (UInt64(1)<<bit)) && continue
            if iszero(pivots[bit+1])
                pivots[bit+1]=row
                break
            end
            row ⊻= pivots[bit+1]
        end
    end
    # Fixed descending pivot order; zero and dependent rows add no dimension.
    basis=UInt64[]
    for bit in n-1:-1:0
        !iszero(pivots[bit+1]) && push!(basis,pivots[bit+1])
    end
    basis
end

function _scalar_signature(basis::AbstractVector{UInt64}, state::UInt64)
    length(basis)<=14 || throw(ArgumentError("signature rank cap"))
    # The caller's unit-action guard covers this bounded, allocation-free loop;
    # no extra clock/callback per parity bit is needed.
    signature=UInt64(0)
    for (j,mask) in enumerate(basis)
        signature |= UInt64(isodd(count_ones(mask&state))) << (j-1)
    end
    signature
end

function _permutation(order, m)
    order isa AbstractVector && length(order)==m && all(x->x isa Integer && !(x isa Bool),order) ||
        throw(ArgumentError("complete integer order required"))
    out=Int.(order)
    sort(out)==collect(1:m) || throw(ArgumentError("order is not a permutation"))
    out
end

# Call only after the existing complete-permutation admission. Comparing the
# four vectors by value avoids compiling order-length-specific tuple hashing.
function _distinct_four_orders(orders::AbstractVector)
    length(orders)==4 || return false
    for a in 1:3, b in a+1:4
        orders[a]==orders[b] && return false
    end
    true
end

# Integer-vector equivalent of isless(Tuple(a), Tuple(b)), including equality
# and the prefix rule. Production orders have the same validated length.
function _lexless(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
    for j in 1:min(length(a),length(b))
        if a[j]!=b[j]
            return isless(a[j],b[j])
        end
    end
    length(a)<length(b)
end

function _units(units, n)
    units isa AbstractVector{Unit} && 1<=length(units)<=1024 || throw(ArgumentError("original unit cap"))
    words=0; seen=Set{Tuple{UInt64,UInt64}}()
    for u in units
        ResponseCore.check_budget()
        1<=length(u.coeff)<=8192 || throw(ArgumentError("original nonempty unit required"))
        words+=length(u.coeff);words<=8192 || throw(ArgumentError("original word cap"))
        u.flip < UInt64(1)<<n || throw(ArgumentError("flip outside register"))
        for j in eachindex(u.coeff)
            ResponseCore.check_budget()
            p=u.parity[j];ys=count_ones(u.flip&p)
            p < UInt64(1)<<n && iseven(ys) && u.signs[j]==(isodd(ys÷2) ? -1 : 1) ||
                throw(ArgumentError("invalid original real Pauli encoding"))
            !iszero(u.coeff[j]) && isfinite(u.coeff[j]) && !iszero(u.exact[j]) ||
                throw(ArgumentError("original coefficient must be finite and nonzero"))
            ndigits(abs(numerator(u.exact[j]));base=2)<=4096 && ndigits(denominator(u.exact[j]);base=2)<=4096 ||
                throw(ArgumentError("original coefficient bit cap"))
            Float64(u.exact[j])==u.coeff[j] || throw(ArgumentError("coefficient/exact metadata mismatch"))
            key=(u.flip,p);key in seen && throw(ArgumentError("duplicate global original Pauli word"));push!(seen,key)
            iseven(count_ones(u.flip&(p⊻u.parity[1]))) || throw(ArgumentError("noncommuting original unit"))
        end
    end
end

function _cache_admit!(c)
    ResponseCore.check_budget()
    if length(c.actions)+length(c.kernels)<c.max_cached_actions && c.accounted_cache_bytes+1024<=c.cache_bytes
        c.accounted_cache_bytes+=1024
        c.stats["cache_entries"]+=1
        c.stats["cache_accounted_bytes"]=c.accounted_cache_bytes
        return true
    end
    c.stats["cache_skipped"]+=1
    false
end
function _value(c, i, state)
    ResponseCore.check_budget(:unit_action)
    c.stats["unit_actions"]+=1
    # Equal signatures preserve every individual summand's parity sign, hence
    # the unchanged ordered floating sum AND its guarded exact fallback. Only
    # successful values enter the cache. Always use the same semantic key even
    # if admission limits change or the table fills; raw states can collide
    # numerically with packed signatures of different scalar values.
    key=(i,_scalar_signature(c.scalar_bases[i],state))
    if haskey(c.actions,key)
        c.stats["cache_hits"]+=1
        return c.actions[key]
    end
    c.stats["cache_misses"]+=1
    value=ResponseCore.unitvalue(c.centered[i],state;guard=true)
    _cache_admit!(c) && (c.actions[key]=value)
    value
end
function _apply!(out, c, i, input)
    out !== input || throw(ArgumentError("aliased unit action scratch"))
    empty!(out)
    for (state,value) in input
        ResponseCore.check_budget()
        ResponseCore._accumulate!(out,state⊻c.centered[i].flip,_value(c,i,state)*value,c.max_entries)
    end
    c.stats["peak_vector_entries"]=max(c.stats["peak_vector_entries"],length(out))
    out
end
function _hamiltonian(c, input)
    out=Column();scratch=Column()
    for i in eachindex(c.centered)
        ResponseCore.check_budget()
        _apply!(scratch,c,i,input)
        ResponseCore._add!(out,scratch,1.0,c.max_entries)
        c.stats["peak_vector_entries"]=max(c.stats["peak_vector_entries"],length(out),length(input))
        c.stats["peak_working_entries"]=max(c.stats["peak_working_entries"],length(out)+length(input)+length(scratch))
    end
    out
end
function _working!(c)
    ResponseCore.check_budget()
    total=length(c.exact)+length(c.response)+sum(length,c.buffers)
    total<=12c.max_entries || throw(Limit("support_limit","aggregate working entry cap"))
    c.stats["peak_working_entries"]=max(c.stats["peak_working_entries"],total)
    c.stats["peak_vector_entries"]=max(c.stats["peak_vector_entries"],length(c.exact),length(c.response),maximum(length,c.buffers))
end

"""One HF column, original labels and exact supplied coefficient metadata.

No window factors or pair/triple tables are constructed. Construction itself is
cooperative; the caller must also own an OS CPU/memory boundary for startup/JIT.
"""
function FullOrderContext(units::AbstractVector{Unit}, n::Integer, hf::Integer, budget::Integer;
        time::Real=1.0,max_entries::Int=50000,max_cached_actions::Int=100000,
        cache_bytes::Int=64*1024^2,check=noop)
    ResponseCore.with_budget(check) do
        nn=integer(n,"register",1,14);h=integer(hf,"HF occupation",0,(1<<nn)-1)
        b=integer(budget,"rotation budget",0,typemax(Int)÷2)
        t=positive(time,"training time");t==1.0 || throw(ArgumentError("training time must be one"))
        integer(max_entries,"max_entries",1,50000);integer(max_cached_actions,"max_cached_actions",0,100000)
        integer(cache_bytes,"cache_bytes",0,64*1024^2)
        _units(units,nn)
        b==30sum(length(u.coeff) for u in units) || throw(ArgumentError("training budget must be original 30U"))
        centered=ResponseCore.center_units(units,h)
        scalar_bases=[_parity_basis(u.parity,nn) for u in centered]
        commuting=ResponseCore._all_commuting(centered)
        stats=Dict(k=>0 for k in ("unit_actions","cache_hits","cache_misses","cache_skipped",
            "cache_entries","cache_accounted_bytes","kernel_hits","kernel_misses",
            "peak_vector_entries","peak_working_entries","full_scores"))
        c=FullOrderContext(collect(units),centered,scalar_bases,UInt64(h),nn,b,t,max_entries,max_cached_actions,
            cache_bytes,0,Dict{Tuple{Int,UInt64},Float64}(),Dict{UInt64,ComplexF64}(),Column(),commuting,
            [Column() for _ in 1:10],Response(),stats)
        if !commuting
            column=Column(c.hf=>1.0)
            for _ in 1:3
                ResponseCore.check_budget();column=_hamiltonian(c,column)
            end
            for state in keys(column)
                ResponseCore.check_budget();column[state]/=6.0
            end
            c.exact=column
        end
        _working!(c)
        c
    end
end

function _column!(c, order)
    pi=_permutation(order,length(c.units))
    for buffer in c.buffers;ResponseCore.check_budget();empty!(buffer);end
    empty!(c.response)
    c.commuting && return (c.buffers[4],0.0,pi)
    q=c.buffers[1:4];out=c.buffers[5:8]
    q[1][c.hf]=1.0
    # Every stage is evaluated, including both global endpoints. This is the
    # complete S3 Taylor recurrence, not a conditioned-window approximation.
    stages=Iterators.flatten((((i,0.5) for i in pi[1:end-1]),((last(pi),1.0),),
                             ((i,0.5) for i in reverse(pi[1:end-1]))))
    for (i,tau) in stages
        ResponseCore.check_budget()
        for degree in 1:4
            empty!(out[degree]);ResponseCore._add!(out[degree],q[degree],1.0,c.max_entries)
        end
        for degree in 0:2
            power=q[degree+1];factor=1.0
            for exponent in 1:3-degree
                target=c.buffers[isodd(exponent) ? 9 : 10]
                power=_apply!(target,c,i,power)
                factor*=tau/exponent
                ResponseCore._add!(out[degree+exponent+1],power,factor,c.max_entries)
            end
        end
        q,out=out,q
        _working!(c)
    end
    scale=ResponseCore._norm(q[4])+ResponseCore._norm(c.exact)
    iszero(scale) && throw(NumericalResolution("zero cubic scale is not structurally certified"))
    ResponseCore._add!(q[4],c.exact,-1.0,c.max_entries)
    ResponseCore._resolved(ResponseCore._norm(q[4]),scale,true)
    q[4],scale,pi
end

"""Return an owned full cubic column for independent identity checks."""
function full_column(c, order;check=noop)
    ResponseCore.with_budget(check) do
        column,scale,_=_column!(c,order)
        (column=copy(column),resolution_scale=scale)
    end
end
function _kernel(c,destination)
    if haskey(c.kernels,destination)
        c.stats["kernel_hits"]+=1
        return c.kernels[destination]
    end
    c.stats["kernel_misses"]+=1
    value=ResponseCore._kernel(ResponseCore._diagonal_gap(c.units,c.hf,destination),c.time)
    _cache_admit!(c) && (c.kernels[destination]=value)
    value
end
function full_score!(c, order;check=noop)
    ResponseCore.with_budget(check) do
        column,scale,pi=_column!(c,order)
        steps=ResponseCore.affordable_steps(c.units,pi,c.budget)
        for (destination,value) in column
            ResponseCore.check_budget()
            destination==c.hf && continue # HF orthogonal projection, not sector projection.
            ResponseCore._accumulate!(c.response,destination,_kernel(c,destination)*value,c.max_entries)
        end
        ResponseCore._resolved(ResponseCore._norm(c.response),abs(c.time)*scale,true)
        score=ResponseCore._score(c.response,steps)
        _working!(c);c.stats["full_scores"]+=1
        score
    end
end

function insertion(order, source::Int, destination::Int)
    1<=source<=length(order) && 1<=destination<=length(order) && source!=destination ||
        throw(ArgumentError("distinct insertion positions required"))
    out=copy(order);value=out[source];deleteat!(out,source);insert!(out,destination,value);out
end
function metropolis(delta::Real, theta::Real, uniform::Real)
    isfinite(delta) && isfinite(theta) && theta>0 && isfinite(uniform) && 0<uniform<1 ||
        throw(ArgumentError("finite Metropolis arguments required"))
    delta<=0 || log(uniform)<-delta/theta
end
function temperature(settings,s)
    isfinite(s) || throw(ArgumentError("finite cooling fraction required"))
    hi=settings["theta_start"];lo=settings["theta_end"]
    exp(log(hi)+(log(lo)-log(hi))*clamp(Float64(s),0.0,1.0))
end
function _seeds(material)
    material isa String && 0<ncodeunits(material)<=256 && isascii(material) || throw(ArgumentError("bounded ASCII seed material"))
    [foldl((a,b)->(a<<8)|UInt64(b),sha256(material*":molecular-generalization-sa-v1:chain"*string(k))[1:8];init=UInt64(0)) for k in 1:4]
end

"""Finite four-start whole-order insertion Metropolis search.

cpu_deadline is an ABSOLUTE owned-process CPU ceiling after serialization
reserve. A controlled clock can exercise slices without timing assertions.
Only complete scores update state; checkpoint snapshots own all mutable data.
"""
function anneal(c, starts, rawsettings;seed_material::String,cpu_deadline::Real,
        parent_value=nothing,clock=cpu_seconds,check=noop,checkpoint=(_)->nothing)
    settings=validate_settings(rawsettings)
    isfinite(cpu_deadline) && cpu_deadline>0 || throw(ArgumentError("finite CPU deadline required"))
    starts isa AbstractVector && length(starts)==4 || throw(ArgumentError("four starts required"))
    orders=[_permutation(p,length(c.units)) for p in starts]
    _distinct_four_orders(orders) || throw(ArgumentError("four distinct starts required"))
    seeds=_seeds(seed_material)
    state=Dict{String,Any}("status"=>"initializing","order"=>Int[],"score"=>0.0,
        "cubic_response"=>0.0,"steps"=>0,"start_scores_completed"=>0,"evaluations"=>0,
        "proposals"=>0,"accepted_moves"=>0,"uphill_attempts"=>0,"uphill_accepted"=>0,
        "chain_proposals"=>zeros(Int,4),"chain_accepted"=>zeros(Int,4),"chain_CPU_seconds"=>zeros(4),
        "seeds"=>string.(seeds),"stop_reason"=>"","kernel_stats"=>Dict{String,Int}())
    best=nothing;best_order=Int[];scores=Any[];normalizer=0.0
    function global_check(kind=:loop,amount=1)
        check(kind,amount)
        clock()<cpu_deadline || throw(Limit("cpu_limit","owned process search CPU exhausted"))
        c.stats["unit_actions"]<=settings["max_actions"] || throw(Limit("action_limit","sparse action ceiling"))
    end
    function publish()
        isnothing(best) && return
        state["order"]=copy(best_order);state["score"]=best.value
        state["cubic_response"]=best.cubic_response;state["steps"]=best.steps
        state["kernel_stats"]=copy(c.stats)
        checkpoint(deepcopy(state))
    end
    function candidate_score(order,callback)
        state["evaluations"]<settings["max_evaluations"] || throw(Limit("evaluation_limit","completed scoring ceiling"))
        score=full_score!(c,order;check=callback)
        state["evaluations"]+=1
        score
    end
    try
        for order in orders
            global_check()
            score=candidate_score(order,global_check);push!(scores,score)
            state["start_scores_completed"]+=1
            if isnothing(best) || score.value<best.value || score.value==best.value && _lexless(order,best_order)
                best=score;best_order=copy(order)
            end
            publish()
        end
        parent_J=isnothing(parent_value) ? scores[1].value : Float64(parent_value)
        isfinite(parent_J) && parent_J>=0 || throw(ArgumentError("finite parent energy normalization"))
        normalizer=max(parent_J,settings["energy_floor"])
        search_start=clock();slice=(cpu_deadline-search_start)/4
        slice>0 || throw(Limit("cpu_limit","setup exhausted search allowance"))
        for chain in 1:4
            began=clock();ending=search_start+chain*slice
            current=copy(orders[chain]);score=scores[chain];rng=Xoshiro(seeds[chain])
            local_check=(kind=:loop,amount=1)->begin
                global_check(kind,amount)
                clock()<ending || throw(SliceExpired())
            end
            try
                while true
                    local_check(:proposal)
                    state["proposals"]<settings["max_proposals"] || throw(Limit("proposal_limit","proposal ceiling"))
                    source=rand(rng,1:length(current));destination=rand(rng,1:length(current)-1)
                    destination>=source && (destination+=1)
                    proposal=insertion(current,source,destination)
                    theta=temperature(settings,(clock()-(ending-slice))/slice)
                    state["proposals"]+=1
                    actual=candidate_score(proposal,local_check)
                    state["chain_proposals"][chain]+=1
                    delta=actual.value/normalizer-score.value/normalizer
                    isfinite(delta) || throw(NumericalResolution("nonfinite normalized energy difference"))
                    uphill=delta>0;uphill && (state["uphill_attempts"]+=1)
                    uniform=rand(rng);while iszero(uniform);uniform=rand(rng);end
                    if metropolis(delta,theta,uniform)
                        current=proposal;score=actual
                        state["accepted_moves"]+=1;state["chain_accepted"][chain]+=1
                        uphill && (state["uphill_accepted"]+=1)
                    end
                    if actual.value<best.value || actual.value==best.value && _lexless(proposal,best_order)
                        best=actual;best_order=copy(proposal);state["status"]="checked_update";publish()
                    end
                end
            catch e
                e isa SliceExpired || rethrow()
            finally
                state["chain_CPU_seconds"][chain]=max(0.0,clock()-began)
            end
            state["status"]="slice_complete";publish()
        end
        state["status"]="complete"
    catch e
        if e isa Limit
            state["status"]=e.status;state["stop_reason"]=e.reason
        elseif e isa NumericalResolution
            state["status"]="numerical_limit";state["stop_reason"]=sprint(showerror,e)
        elseif e isa ArgumentError
            state["status"]="support_or_admission_limit";state["stop_reason"]=sprint(showerror,e)
        else
            rethrow()
        end
    end
    publish()
    state
end

end # module
