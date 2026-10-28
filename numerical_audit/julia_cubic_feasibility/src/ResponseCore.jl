module ResponseCore

export Unit, unitvalue, apply_unit, center_units, atom_column, cubic_column,
       diagonal_response, affordable_steps, score_order, build_window_model,
       window_response, window_score, WindowModel, ResponseFactor,
       AmplitudeCache, NumericalResolution

const RealColumn = Dict{UInt64,Float64}
const ResponseVector = Dict{UInt64,ComplexF64}
const OrderKey = Tuple{Vararg{Int}}
const RESOLUTION_RATIO = 2.0^-24

struct NumericalResolution <: Exception
    message::String
end
Base.showerror(io::IO, e::NumericalResolution) = print(io, e.message)

"""Equal-flip real Pauli unit. signs encode each word's real Y phase.

The caller authenticates Hermiticity, conservation, register size, and supplied
word counts. Empty units are allowed internally after scalar centering.
"""
struct Unit
    flip::UInt64
    parity::Vector{UInt64}
    signs::Vector{Int8}
    coeff::Vector{Float64}
    function Unit(flip::UInt64, parity::Vector{UInt64}, signs::Vector{Int8}, coeff::Vector{Float64})
        length(parity) == length(signs) == length(coeff) || throw(ArgumentError("unit array lengths differ"))
        all(s -> s == 1 || s == -1, signs) || throw(ArgumentError("real Pauli phases must be +/-1"))
        all(isfinite, coeff) || throw(ArgumentError("unit coefficients must be finite"))
        new(flip, copy(parity), copy(signs), copy(coeff))
    end
end

"""Scalar action: A|state> = unitvalue(A,state)|state xor A.flip>."""
function unitvalue(unit::Unit, state::UInt64)
    total = 0.0
    correction = 0.0
    @inbounds for j in eachindex(unit.coeff)
        sign = isodd(count_ones(state & unit.parity[j])) ? -unit.signs[j] : unit.signs[j]
        term = unit.coeff[j] * sign
        y = term - correction
        updated = total + y
        correction = (updated - total) - y
        total = updated
    end
    isfinite(total) || throw(NumericalResolution("nonfinite unit action"))
    total
end
unitvalue(unit::Unit, state::Integer) = unitvalue(unit, UInt64(state))

function _accumulate!(target::Dict{UInt64,T}, key::UInt64, value, maximum::Int) where T
    iszero(value) && return target
    updated = get(target, key, zero(T)) + value
    isfinite(updated) || throw(NumericalResolution("nonfinite response accumulation"))
    if iszero(updated)
        delete!(target, key)
    else
        !haskey(target, key) && length(target) >= maximum && throw(ArgumentError("response support cap exceeded"))
        target[key] = updated
    end
    target
end

function _add!(target::Dict{UInt64,T}, source, factor, maximum::Int) where T
    for (state, value) in source
        _accumulate!(target, state, factor * value, maximum)
    end
    target
end

function _norm(vector)
    result = 0.0
    for value in values(vector)
        result = hypot(result, abs(value))
    end
    isfinite(result) || throw(NumericalResolution("nonfinite response norm"))
    result
end

function _resolved(value::Float64, scale::Float64, guard::Bool)
    guard && scale > 0 && value <= RESOLUTION_RATIO * scale &&
        throw(NumericalResolution("response lies below the declared floating cancellation threshold"))
    nothing
end

function apply_unit(unit::Unit, vector::AbstractDict; max_entries::Int=50000)
    output = RealColumn()
    for (state0, value) in vector
        state = UInt64(state0)
        _accumulate!(output, state ⊻ unit.flip, unitvalue(unit, state) * value, max_entries)
    end
    output
end

"""Subtract each diagonal unit's HF scalar; preserve all labels.

These are response-only copies. Rotation accounting must use the ORIGINAL
units, since the numerical scalar term is not an implemented rotation.
"""
function center_units(units::AbstractVector{Unit}, hf::Integer)
    state = UInt64(hf)
    result = Unit[]
    for unit in units
        if iszero(unit.flip)
            shift = unitvalue(unit, state)
            parity, signs, coeff = copy(unit.parity), copy(unit.signs), copy(unit.coeff)
            if !iszero(shift)
                push!(parity, UInt64(0)); push!(signs, Int8(1)); push!(coeff, -shift)
            end
            push!(result, Unit(unit.flip, parity, signs, coeff))
        else
            push!(result, unit)
        end
    end
    result
end

function _validate_order(order, m::Int)
    result = Int[order...]
    sort(result) == collect(1:m) || throw(ArgumentError("order must contain every unit exactly once"))
    result
end

function _all_commuting(units)
    for a in eachindex(units), b in 1:a-1
        ua, ub = units[a], units[b]
        for i in eachindex(ua.coeff), j in eachindex(ub.coeff)
            (iszero(ua.coeff[i]) || iszero(ub.coeff[j])) && continue
            isodd(count_ones(ua.flip & ub.parity[j]) + count_ones(ub.flip & ua.parity[i])) && return false
        end
    end
    true
end

function _hamiltonian_action(units, vector, maximum)
    output = RealColumn()
    for unit in units
        _add!(output, apply_unit(unit, vector; max_entries=maximum), 1.0, maximum)
    end
    output
end

function _cubic_with_scale(units::AbstractVector{Unit}, order, hf::Integer;
                           max_entries::Int=50000, guard::Bool=true, center::Bool=true)
    permutation = _validate_order(order, length(units))
    isempty(permutation) && return (RealColumn(), 0.0)
    _all_commuting(units) && return (RealColumn(), 0.0)
    response_units = center ? center_units(units, hf) : units
    initial = RealColumn(UInt64(hf) => 1.0)
    exact = initial
    for _ in 1:3
        exact = _hamiltonian_action(response_units, exact, max_entries)
    end
    for state in keys(exact)
        exact[state] /= 6.0
    end
    q = [initial, RealColumn(), RealColumn(), RealColumn()]
    stages = Tuple{Int,Float64}[(i, 0.5) for i in permutation[1:end-1]]
    push!(stages, (last(permutation), 1.0))
    append!(stages, [(i, 0.5) for i in reverse(permutation[1:end-1])])
    for (index, tau) in stages
        output = [copy(column) for column in q]
        for degree in 0:2
            power = q[degree + 1]
            factor = 1.0
            for exponent in 1:3-degree
                power = apply_unit(response_units[index], power; max_entries=max_entries)
                factor *= tau / exponent
                _add!(output[degree + exponent + 1], power, factor, max_entries)
            end
        end
        q = output
    end
    scale = _norm(q[4]) + _norm(exact)
    guard && iszero(scale) && throw(NumericalResolution("zero cubic scale is not structurally certified"))
    _add!(q[4], exact, -1.0, max_entries)
    _resolved(_norm(q[4]), scale, guard)
    q[4], scale
end

"""Independent sparse Taylor recurrence through degree three, minus H^3/6."""
function cubic_column(units::AbstractVector{Unit}, order, hf::Integer; kwargs...)
    first(_cubic_with_scale(units, order, hf; kwargs...))
end

mutable struct AmplitudeCache
    units::Vector{Unit}
    amplitudes::Dict{Tuple{Int,UInt64},Float64}
    max_entries::Int
    hits::Int
    misses::Int
end
AmplitudeCache(units::AbstractVector{Unit}; max_entries::Int=100000) =
    AmplitudeCache(collect(units), Dict{Tuple{Int,UInt64},Float64}(), max_entries, 0, 0)

function _amplitude(cache::AmplitudeCache, index::Int, state::UInt64)
    key = (index, state)
    if haskey(cache.amplitudes, key)
        cache.hits += 1
        return cache.amplitudes[key]
    end
    cache.misses += 1
    value = unitvalue(cache.units[index], state)
    length(cache.amplitudes) < cache.max_entries && (cache.amplitudes[key] = value)
    value
end

function _product3(cache::AmplitudeCache, a::Int, b::Int, c::Int, state::UInt64)
    vc = _amplitude(cache, c, state)
    sb = state ⊻ cache.units[c].flip
    vb = _amplitude(cache, b, sb)
    sa = sb ⊻ cache.units[b].flip
    vc * vb * _amplitude(cache, a, sa)
end

function _atom_column(cache::AmplitudeCache, order::Tuple, hf::UInt64, maximum::Int)
    result = RealColumn()
    p(a,b,c) = _product3(cache, a, b, c, hf)
    if length(order) == 2
        a, b = order
        aa_b = p(a,a,b) - 2p(a,b,a) + p(b,a,a)
        b_ab = 2p(b,a,b) - p(b,b,a) - p(a,b,b)
        _accumulate!(result, hf ⊻ cache.units[b].flip, -aa_b / 24.0, maximum)
        _accumulate!(result, hf ⊻ cache.units[a].flip, -b_ab / 12.0, maximum)
    elseif length(order) == 3
        a, b, c = order
        # Distinct-label part of E(A,B+C)+E(B,C).
        bac = p(b,a,c) - p(b,c,a) - p(a,c,b) + p(c,a,b)
        cab = p(c,a,b) - p(c,b,a) - p(a,b,c) + p(b,a,c)
        destination = hf ⊻ cache.units[a].flip ⊻ cache.units[b].flip ⊻ cache.units[c].flip
        _accumulate!(result, destination, -(bac + cab) / 12.0, maximum)
    else
        throw(ArgumentError("an atom has two or three distinct labels"))
    end
    result
end

"""Cubic atom for an induced two/three-label order (not a full local defect).

To reuse a cache, construct it from center_units(units,hf), then pass center=false
and that cache. With a supplied cache its compiled units are authoritative and
`center` applies only when constructing a new cache. The unit-list cardinality
must match; callers must supply the cache compiled for their declared units.
Scalars cancel algebraically; centering improves arithmetic.
"""
function atom_column(units::AbstractVector{Unit}, order, hf::Integer;
                     center::Bool=true, cache::Union{Nothing,AmplitudeCache}=nothing,
                     max_entries::Int=50000)
    local_order = Tuple(Int[i for i in order])
    length(unique(local_order)) == length(local_order) || throw(ArgumentError("atom labels repeat"))
    all(i -> 1 <= i <= length(units), local_order) || throw(ArgumentError("atom label out of range"))
    isnothing(cache) || length(cache.units)==length(units) || throw(ArgumentError("atom cache unit count differs"))
    context = isnothing(cache) ? AmplitudeCache(center ? center_units(units,hf) : units) : cache
    _atom_column(context, local_order, UInt64(hf), max_entries)
end

function _diagonal_gap(units, source::UInt64, destination::UInt64)
    total = 0.0
    correction = 0.0
    for unit in units
        iszero(unit.flip) || continue
        for j in eachindex(unit.coeff)
            zs = isodd(count_ones(source & unit.parity[j])) ? -1 : 1
            zd = isodd(count_ones(destination & unit.parity[j])) ? -1 : 1
            term = unit.coeff[j] * unit.signs[j] * (zd-zs)
            y = term - correction
            updated = total + y
            correction = (updated-total)-y
            total = updated
        end
    end
    isfinite(total) || throw(NumericalResolution("nonfinite diagonal gap"))
    total
end

function _kernel(gap::Float64, time::Float64)
    x = time * gap / 2.0
    isfinite(x) || throw(NumericalResolution("nonfinite diagonal phase"))
    ratio = iszero(x) ? 1.0 : sin(x) / x
    value = time * ratio * cis(x)
    isfinite(value) || throw(NumericalResolution("nonfinite diagonal transport"))
    value
end

"""Apply HF diagonal transport and projection; never create a full state array."""
function diagonal_response(units::AbstractVector{Unit}, column::AbstractDict, hf::Integer;
                           time::Real=1.0, max_entries::Int=50000, guard::Bool=false,
                           resolution_scale::Real=0.0)
    t = Float64(time)
    isfinite(t) || throw(ArgumentError("time must be finite"))
    source = UInt64(hf)
    response = ResponseVector()
    for (destination0, value) in column
        destination = UInt64(destination0)
        destination == source && continue
        kernel = _kernel(_diagonal_gap(units, source, destination), t)
        _accumulate!(response, destination, kernel * value, max_entries)
    end
    _resolved(_norm(response), abs(t) * Float64(resolution_scale), guard)
    response
end

"""Affordable fused Strang steps using ORIGINAL unit word counts."""
function affordable_steps(units::AbstractVector{Unit}, order, budget::Integer)
    budget >= 0 || throw(ArgumentError("negative rotation budget"))
    permutation = _validate_order(order, length(units))
    lengths = [length(units[i].coeff) for i in permutation if !isempty(units[i].coeff)]
    isempty(lengths) && return 1
    total = sum(lengths)
    steps = length(lengths) == 1 ? (budget >= total ? 1 : 0) :
        max(0, fld(budget-first(lengths), 2total-first(lengths)-last(lengths)))
    steps >= 1 || throw(ArgumentError("rotation budget cannot afford a symmetric step"))
    steps
end

function _score(response, steps)
    norm_response = _norm(response)
    numerator = norm_response^2
    norm_response > 0 && iszero(numerator) && throw(NumericalResolution("squared response underflows"))
    denominator = Float64(steps)^4
    isfinite(denominator) || throw(NumericalResolution("step-count normalization overflows"))
    value = numerator / denominator
    numerator > 0 && iszero(value) && throw(NumericalResolution("normalized response score underflows"))
    isfinite(value) || throw(NumericalResolution("nonfinite response score"))
    (value=value, cubic_response=numerator, steps=steps)
end

function score_order(units::AbstractVector{Unit}, order, hf::Integer, budget::Integer;
                     time::Real=1.0, max_entries::Int=50000, guard::Bool=true)
    steps = affordable_steps(units, order, budget)
    column, scale = _cubic_with_scale(units, order, hf; max_entries=max_entries, guard=guard)
    response = diagonal_response(units, column, hf; time=time, max_entries=max_entries,
                                 guard=guard, resolution_scale=scale)
    _score(response, steps)
end

"""Frozen response DELTAS for all local permutations of one scope."""
struct ResponseFactor
    scope::Vector{Int}
    variants::Dict{OrderKey,ResponseVector}
end

struct WindowModel
    units::Vector{Unit}                 # ORIGINAL counts for budget accounting.
    baseline_order::Vector{Int}
    baseline_window::Vector{Int}
    window_positions::UnitRange{Int}
    hf::UInt64
    time::Float64
    baseline_column::RealColumn
    baseline_response::ResponseVector
    factors::Vector{ResponseFactor}
    max_entries::Int
    guard::Bool
    setup_stats::NamedTuple
end

_permutations2(a,b) = ((a,b),(b,a))
_permutations3(a,b,c) = ((a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a))

function _freeze_deltas(scope, raw, baseline_rank, maximum)
    baseline_key = Tuple(sort(collect(scope); by=i -> baseline_rank[i]))
    reference = raw[baseline_key]
    deltas = Dict{OrderKey,ResponseVector}()
    for (key, vector) in raw
        delta = ResponseVector()
        if key != baseline_key
            _add!(delta, vector, 1.0, maximum)
            _add!(delta, reference, -1.0, maximum)
        end
        deltas[key] = delta
    end
    ResponseFactor(sort(collect(scope)), deltas)
end

"""Condition a CONTIGUOUS window while retaining full outside interference.

Only atoms with >=2 window labels vary. Outside-third triple atoms are folded
into pair tables. The baseline is computed by the independent Taylor scorer.
Tables contain frozen Float64/ComplexF64 deltas, not exact kernel values.
Individual atom/DELTA construction has no interval or relative-error certificate.
The numerical resolution guard applies to the independent full Taylor baseline
and final coherent sum; it does not certify every small local delta. Changing
transport/aggregation order defines new supplied tables and needs separate
validation rather than inheriting this frozen study's exact optimum.
"""
function build_window_model(units::AbstractVector{Unit}, baseline_order, hf::Integer,
                            window_positions::AbstractUnitRange{<:Integer};
                            time::Real=1.0, max_entries::Int=50000,
                            max_cached_actions::Int=100000,
                            max_factor_entries::Int=1_000_000, guard::Bool=true)
    order = _validate_order(baseline_order, length(units))
    2 <= length(window_positions) <= 9 || throw(ArgumentError("pilot window width must be between two and nine"))
    max_factor_entries >= 1 || throw(ArgumentError("factor entry cap must be positive"))
    lo, hi = Int(first(window_positions)), Int(last(window_positions))
    1 <= lo <= hi <= length(order) || throw(ArgumentError("window positions outside order"))
    labels = order[lo:hi]
    baseline_rank = invperm(order)
    column, scale = _cubic_with_scale(units, order, hf; max_entries=max_entries, guard=guard)
    baseline = diagonal_response(units, column, hf; time=time, max_entries=max_entries,
                                 guard=guard, resolution_scale=scale)
    centered = center_units(units, hf)
    cache = AmplitudeCache(centered; max_entries=max_cached_actions)
    raw_atoms = 0
    factor_entries = 0
    factors = ResponseFactor[]
    function retain_factor!(factor)
        entries = sum(length(v) for v in values(factor.variants))
        factor_entries + entries <= max_factor_entries || throw(ArgumentError("conditioned factor storage cap exceeded"))
        factor_entries += entries
        push!(factors,factor)
    end
    outside = [i for i in order if !(i in labels)]
    function transported(local_order)
        raw_atoms += 1
        local_column = _atom_column(cache, local_order, UInt64(hf), max_entries)
        diagonal_response(units, local_column, hf; time=time, max_entries=max_entries)
    end
    for ia in 1:length(labels)-1, ib in ia+1:length(labels)
        a,b = labels[ia], labels[ib]
        raw = Dict{OrderKey,ResponseVector}()
        for pair_order in _permutations2(a,b)
            vector = transported(pair_order)
            for c in outside
                local_order = baseline_rank[c] < lo ? (c,pair_order...) : (pair_order...,c)
                _add!(vector, transported(local_order), 1.0, max_entries)
            end
            raw[pair_order] = vector
        end
        retain_factor!(_freeze_deltas((a,b), raw, baseline_rank, max_entries))
    end
    for ia in 1:length(labels)-2, ib in ia+1:length(labels)-1, ic in ib+1:length(labels)
        a,b,c = labels[ia],labels[ib],labels[ic]
        raw = Dict{OrderKey,ResponseVector}()
        for local_order in _permutations3(a,b,c)
            raw[local_order] = transported(local_order)
        end
        retain_factor!(_freeze_deltas((a,b,c), raw, baseline_rank, max_entries))
    end
    stats = (raw_atom_variants=raw_atoms, factors=length(factors),
             cached_actions=length(cache.amplitudes), cache_hits=cache.hits, cache_misses=cache.misses,
             delta_entries=sum(length(v) for factor in factors for v in values(factor.variants)))
    WindowModel(collect(units), order, labels, lo:hi, UInt64(hf), Float64(time),
                column, baseline, factors, max_entries, guard, stats)
end

function _window_order(model::WindowModel, window_order)
    local_order = Int[window_order...]
    sort(local_order) == sort(model.baseline_window) || throw(ArgumentError("window order has different labels"))
    order = copy(model.baseline_order)
    order[model.window_positions] = local_order
    order
end

function window_response(model::WindowModel, window_order; guard::Bool=model.guard)
    order = _window_order(model, window_order)
    rank = invperm(order)
    response = copy(model.baseline_response)
    scale = _norm(response)
    for factor in model.factors
        key = Tuple(sort(factor.scope; by=i -> rank[i]))
        delta = factor.variants[key]
        scale += _norm(delta)
        _add!(response, delta, 1.0, model.max_entries)
    end
    _resolved(_norm(response), scale, guard)
    response
end

function window_score(model::WindowModel, window_order, budget::Integer; guard::Bool=model.guard)
    order = _window_order(model, window_order)
    steps = affordable_steps(model.units, order, budget)
    _score(window_response(model, window_order; guard=guard), steps)
end

end # module
