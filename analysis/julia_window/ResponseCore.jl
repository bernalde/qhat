module ResponseCore

# Product local target v2; wire protocol v1. Historical sources remain immutable. This scoped
# callback is separate from floating-resolution `guard` and is used only in the
# single-threaded, single-request worker. It also meters cumulative atom/actions.
const BUDGET_CHECK = Ref{Function}((kind, amount) -> nothing)
check_budget(kind::Symbol=:loop, amount::Int=1) = BUDGET_CHECK[](kind, amount)
function with_budget(f::Function, callback::Function)
    previous = BUDGET_CHECK[]
    BUDGET_CHECK[] = callback
    try
        check_budget()
        f()
    finally
        BUDGET_CHECK[] = previous
    end
end

export Unit, unitvalue, apply_unit, center_units, atom_column, cubic_column,
       diagonal_response, affordable_steps, score_order, build_window_model,
       window_response, window_score, WindowModel, ResponseFactor,
       AmplitudeCache, NumericalResolution, with_budget, check_budget

const RealColumn = Dict{UInt64,Float64}
const ResponseVector = Dict{UInt64,ComplexF64}
const OrderKey = Tuple{Vararg{Int}}
const RESOLUTION_RATIO = 2.0^-24
const LOCAL_TARGET = "exact_real_deltas_before_transport_v2"
const EXACT_CACHE_BYTE_LIMIT = 64*1024^2
const ExactScalar = Rational{BigInt}
const ExactColumn = Dict{UInt64,ExactScalar}
const EXACT_ZERO = big(0)//big(1)

# These helpers bound every new local-arithmetic integer before products/adds.
_integer_bits(value::BigInt) = iszero(value) ? 0 : ndigits(abs(value);base=2)
function _exact_int_product(a::BigInt,b::BigInt)
    check_budget()
    max(_integer_bits(a),_integer_bits(b))<=16384 || throw(ArgumentError("local integer bit cap"))
    (iszero(a)||iszero(b)) && return big(0)
    a==1 && return b
    b==1 && return a
    a==-1 && return -b
    b==-1 && return -a
    _integer_bits(a)+_integer_bits(b)<=16384 || throw(ArgumentError("local prospective product bit cap"))
    a*b
end
function _exact_int_sum(a::BigInt,b::BigInt)
    check_budget()
    bits=max(_integer_bits(a),_integer_bits(b))
    bits<=16384 && (sign(a)!=sign(b)||bits<16384) || throw(ArgumentError("local prospective sum bit cap"))
    a+b
end
function _qadd(a::ExactScalar,b::ExactScalar)
    _check_rational_bits(a);_check_rational_bits(b)
    iszero(a) && return b
    iszero(b) && return a
    common=gcd(denominator(a),denominator(b));check_budget()
    left=div(denominator(b),common);right=div(denominator(a),common)
    num=_exact_int_sum(_exact_int_product(numerator(a),left),_exact_int_product(numerator(b),right))
    iszero(num) && return EXACT_ZERO
    den=_exact_int_product(denominator(a),left)
    _check_rational_bits(num//den)
end
function _qneg(a::ExactScalar)
    _check_rational_bits(a)
    -a
end
_qsub(a::ExactScalar,b::ExactScalar)=_qadd(a,_qneg(b))
function _qmul(a::ExactScalar,b::ExactScalar)
    _check_rational_bits(a);_check_rational_bits(b)
    (iszero(a)||iszero(b)) && return EXACT_ZERO
    first=gcd(abs(numerator(a)),denominator(b));check_budget()
    second=gcd(abs(numerator(b)),denominator(a));check_budget()
    num=_exact_int_product(div(numerator(a),first),div(numerator(b),second))
    den=_exact_int_product(div(denominator(a),second),div(denominator(b),first))
    _check_rational_bits(num//den)
end
_qscale(a::ExactScalar,n::Int,d::Int=1)=_qmul(a,big(n)//big(d))

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
    exact::Vector{Rational{BigInt}}
    function Unit(flip::UInt64, parity::Vector{UInt64}, signs::Vector{Int8}, coeff::Vector{Float64};
                  exact::Vector{Rational{BigInt}}=Rational{BigInt}.(coeff))
        check_budget()
        length(parity) == length(signs) == length(coeff) == length(exact) || throw(ArgumentError("unit array lengths differ"))
        all(s -> s == 1 || s == -1, signs) || throw(ArgumentError("real Pauli phases must be +/-1"))
        all(isfinite, coeff) || throw(ArgumentError("unit coefficients must be finite"))
        new(flip, copy(parity), copy(signs), copy(coeff), copy(exact))
    end
end

"""Exact rational comparison with bounded cross-products for rounding proofs."""
function _rounding_compare(a::Rational{BigInt}, b::Rational{BigInt})
    _check_rational_bits(a); _check_rational_bits(b)
    an, ad = numerator(a), denominator(a)
    bn, bd = numerator(b), denominator(b)
    check_budget()
    sign(an) != sign(bn) && return cmp(sign(an), sign(bn))
    iszero(an) && return 0 # Both numerators are zero after the sign check.
    # Divide positive common factors before comparing integer products. This
    # keeps dyadic midpoint proofs small without enlarging the intermediate cap.
    common_numerator = gcd(abs(an), abs(bn))
    common_denominator = gcd(ad, bd)
    check_budget()
    left_numerator = div(an, common_numerator)
    right_numerator = div(bn, common_numerator)
    left_denominator = div(bd, common_denominator)
    right_denominator = div(ad, common_denominator)
    ndigits(abs(left_numerator);base=2) + ndigits(left_denominator;base=2) <= 16384 &&
        ndigits(abs(right_numerator);base=2) + ndigits(right_denominator;base=2) <= 16384 ||
        throw(ArgumentError("rounding comparison bit cap"))
    check_budget()
    left = left_numerator * left_denominator
    right = right_numerator * right_denominator
    check_budget()
    cmp(left, right)
end

function _rounding_cell_contains(exact::Rational{BigInt}, candidate::Float64)
    check_budget()
    # Nonzero exact amplitudes rounding to zero are deliberately unsupported.
    isfinite(candidate) && !iszero(candidate) || return false
    below, above = prevfloat(candidate), nextfloat(candidate)
    isfinite(below) && isfinite(above) || return false
    center = _check_rational_bits(Rational{BigInt}(candidate))
    low = _check_rational_bits(Rational{BigInt}(below))
    high = _check_rational_bits(Rational{BigInt}(above))
    lower = _check_rational_bits((low + center) / 2)
    upper = _check_rational_bits((center + high) / 2)
    lower_comparison = _rounding_compare(exact, lower)
    upper_comparison = _rounding_compare(exact, upper)
    # The stored significand's low bit determines round-to-nearest ties-to-even
    # for both signs, subnormals and the normal/subnormal boundary.
    even = iseven(reinterpret(UInt64, candidate))
    (lower_comparison > 0 || (iszero(lower_comparison) && even)) &&
        (upper_comparison < 0 || (iszero(upper_comparison) && even))
end

"""Correctly rounded nonzero binary64 action, authenticated by exact midpoints.

Float64(exact) is only a proposal: an intermediate BigFloat conversion can
double-round near a midpoint. Check that proposal and at most its two neighbors;
return only a proved rounding cell. No precision or rounding mode is mutated.
"""
function _rounded_exact_action(exact::Rational{BigInt})
    _check_rational_bits(exact)
    iszero(exact) && return 0.0
    check_budget()
    candidate = try
        Float64(exact)
    catch error
        error isa Union{OverflowError,InexactError} || rethrow()
        throw(NumericalResolution("exact unit action conversion outside resolution"))
    end
    isfinite(candidate) || throw(NumericalResolution("nonfinite exact unit action conversion"))
    for value in (candidate, prevfloat(candidate), nextfloat(candidate))
        check_budget()
        _rounding_cell_contains(exact, value) && return value
    end
    throw(NumericalResolution("exact unit action has no certified supported nonzero binary64 rounding"))
end

"""Scalar action: A|state> = unitvalue(A,state)|state xor A.flip>."""
function unitvalue(unit::Unit, state::UInt64; guard::Bool=true)
    check_budget(:action)
    total = 0.0
    correction = 0.0
    scale = 0.0
    @inbounds for j in eachindex(unit.coeff)
        check_budget()
        sign = isodd(count_ones(state & unit.parity[j])) ? -unit.signs[j] : unit.signs[j]
        term = unit.coeff[j] * sign
        scale += abs(term)
        y = term - correction
        updated = total + y
        correction = (updated - total) - y
        total = updated
    end
    isfinite(total) || throw(NumericalResolution("nonfinite unit action"))
    if guard && scale>0 && abs(total)<=RESOLUTION_RATIO*scale
        exact=big(0)//big(1)
        for j in eachindex(unit.exact)
            check_budget()
            sign=isodd(count_ones(state&unit.parity[j])) ? -unit.signs[j] : unit.signs[j]
            exact=_qadd(exact,_qscale(unit.exact[j],Int(sign)))
        end
        # Centered diagonal and forbidden equal-flip actions commonly vanish
        # identically. Authenticate that identity instead of refusing roundoff.
        iszero(exact) && return 0.0
        isfinite(scale) || throw(NumericalResolution("nonfinite unit action scale"))
        # This replaces only the unresolved scalar action. Atom, DELTA, full
        # Taylor and global acceptance guards continue to apply independently.
        return _rounded_exact_action(exact)
    end
    _local_resolved(abs(total), scale, guard, "unit action")
    total
end
unitvalue(unit::Unit, state::Integer; kwargs...) = unitvalue(unit, UInt64(state); kwargs...)

function _accumulate!(target::Dict{UInt64,T}, key::UInt64, value, maximum::Int;
                      cap_message::String="response support cap exceeded") where T
    check_budget()
    iszero(value) && return target
    updated = get(target, key, zero(T)) + value
    isfinite(updated) || throw(NumericalResolution("nonfinite response accumulation"))
    if iszero(updated)
        delete!(target, key)
    else
        !haskey(target, key) && length(target) >= maximum && throw(ArgumentError(cap_message))
        target[key] = updated
    end
    target
end

function _add!(target::Dict{UInt64,T}, source, factor, maximum::Int) where T
    for (state, value) in source
        check_budget()
        _accumulate!(target, state, factor * value, maximum)
    end
    target
end

function _norm(vector)
    result = 0.0
    for value in values(vector)
        check_budget()
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

function apply_unit(unit::Unit, vector::AbstractDict; max_entries::Int=50000, guard::Bool=true)
    output = RealColumn()
    for (state0, value) in vector
        check_budget()
        state = UInt64(state0)
        _accumulate!(output, state ⊻ unit.flip, unitvalue(unit, state; guard=guard) * value, max_entries)
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
        check_budget()
        if iszero(unit.flip)
            parity=UInt64[]; signs=Int8[]; coeff=Float64[]; exact=Rational{BigInt}[]
            shift=big(0)//big(1)
            for j in eachindex(unit.coeff)
                check_budget()
                iszero(unit.parity[j]) && continue # Native response drops identity.
                push!(parity,unit.parity[j]);push!(signs,unit.signs[j])
                push!(coeff,unit.coeff[j]);push!(exact,unit.exact[j])
                sign=isodd(count_ones(state&unit.parity[j])) ? -unit.signs[j] : unit.signs[j]
                shift=_qadd(shift,_qscale(unit.exact[j],Int(sign)))
            end
            if !iszero(shift)
                scalar=Float64(shift)
                isfinite(scalar) && scalar!=0 || throw(NumericalResolution("centering conversion outside resolution"))
                push!(parity,UInt64(0));push!(signs,Int8(1));push!(coeff,-scalar);push!(exact,_qneg(shift))
            end
            push!(result,Unit(unit.flip,parity,signs,coeff;exact=exact))
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
        check_budget()
        ua, ub = units[a], units[b]
        for i in eachindex(ua.coeff), j in eachindex(ub.coeff)
            check_budget()
            (iszero(ua.coeff[i]) || iszero(ub.coeff[j])) && continue
            isodd(count_ones(ua.flip & ub.parity[j]) + count_ones(ub.flip & ua.parity[i])) && return false
        end
    end
    true
end

function _hamiltonian_action(units, vector, maximum; guard::Bool=true)
    output = RealColumn()
    for unit in units
        check_budget()
        _add!(output, apply_unit(unit, vector; max_entries=maximum, guard=guard), 1.0, maximum)
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
        check_budget()
        exact = _hamiltonian_action(response_units, exact, max_entries; guard=guard)
    end
    for state in keys(exact)
        check_budget()
        exact[state] /= 6.0
    end
    q = [initial, RealColumn(), RealColumn(), RealColumn()]
    stages = Tuple{Int,Float64}[(i, 0.5) for i in permutation[1:end-1]]
    push!(stages, (last(permutation), 1.0))
    append!(stages, [(i, 0.5) for i in reverse(permutation[1:end-1])])
    for (index, tau) in stages
        check_budget()
        output = [copy(column) for column in q]
        for degree in 0:2
            check_budget()
            power = q[degree + 1]
            factor = 1.0
            for exponent in 1:3-degree
                check_budget()
                power = apply_unit(response_units[index], power; max_entries=max_entries, guard=guard)
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

mutable struct AmplitudeCache{T}
    units::Vector{Unit}
    amplitudes::Dict{Tuple{Int,UInt64},T}
    max_entries::Int
    hits::Int
    misses::Int
    guard::Bool
    max_bytes::Int
    accounted_bytes::Int
    skipped::Int
    kernels::Dict{UInt64,ComplexF64}
    kernel_reference::Union{Nothing,UInt64}
    kernel_time::Union{Nothing,Float64}
end
function AmplitudeCache(units::AbstractVector{Unit}; max_entries::Int=100000, guard::Bool=true,
                        max_bytes::Int=EXACT_CACHE_BYTE_LIMIT)
    check_budget()
    0<=max_entries<=100000 && length(units)<=1024 || throw(ArgumentError("local cache count cap"))
    base=4096+128length(units) # Containers, units-vector copies, empty tables.
    base<=max_bytes<=EXACT_CACHE_BYTE_LIMIT || throw(ArgumentError("local cache byte allowance"))
    T=guard ? ExactScalar : Float64
    AmplitudeCache{T}(collect(units),Dict{Tuple{Int,UInt64},T}(),max_entries,0,0,guard,
        max_bytes,base,0,Dict{UInt64,ComplexF64}(),nothing,nothing)
end

function _cache_entry_bytes(value::ExactScalar)
    _check_rational_bits(value)
    # Julia's own summarysize(BigInt) uses alloc*sizeof(Base.GMP.Limb).
    # Charge allocated capacity, including retained excess after reduction,
    # twice;1024 extra bytes cover objects, keys, slots and resize coexistence.
    a=Int(getfield(numerator(value),:alloc));b=Int(getfield(denominator(value),:alloc))
    a>=0 && b>=0 || throw(ArgumentError("invalid BigInt allocation capacity"))
    1024+2*(a+b)*sizeof(Base.GMP.Limb)
end
_cache_entry_bytes(::Union{Float64,ComplexF64})=1024
function _cache_admit!(cache::AmplitudeCache,bytes::Int)
    check_budget()
    if length(cache.amplitudes)+length(cache.kernels)<cache.max_entries &&
            0<=bytes<=cache.max_bytes-cache.accounted_bytes
        cache.accounted_bytes+=bytes
        return true
    end
    cache.skipped+=1
    false # Recompute exactly next time; never thin an input or clip a value.
end

function _amplitude(cache::AmplitudeCache{Float64}, index::Int, state::UInt64)
    check_budget()
    key = (index, state)
    if haskey(cache.amplitudes, key)
        cache.hits += 1
        return cache.amplitudes[key]
    end
    cache.misses += 1
    value = unitvalue(cache.units[index], state; guard=cache.guard)
    _cache_admit!(cache,_cache_entry_bytes(value)) && (cache.amplitudes[key] = value)
    value
end

function _product3(cache::AmplitudeCache{Float64}, a::Int, b::Int, c::Int, state::UInt64)
    vc = _amplitude(cache, c, state)
    sb = state ⊻ cache.units[c].flip
    vb = _amplitude(cache, b, sb)
    sa = sb ⊻ cache.units[b].flip
    vc * vb * _amplitude(cache, a, sa)
end

function _atom_column(cache::AmplitudeCache{Float64}, order::Tuple, hf::UInt64, maximum::Int; guard::Bool=true)
    guard && throw(ArgumentError("guarded atoms require the exact-only cache"))
    check_budget(:atom)
    result = RealColumn()
    local_scale = 0.0
    function p(a,b,c)
        value = _product3(cache, a, b, c, hf)
        # Conservative local operand norm; exact zeros remain valid atoms.
        local_scale += abs(value)
        value
    end
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
    _local_resolved(_norm(result), local_scale / 6.0, guard, "local atom")
    result
end

function _amplitude(cache::AmplitudeCache{ExactScalar},index::Int,state::UInt64)
    check_budget()
    key=(index,state)
    if haskey(cache.amplitudes,key)
        cache.hits+=1
        return cache.amplitudes[key]
    end
    cache.misses+=1;check_budget(:action)
    unit=cache.units[index];value=EXACT_ZERO
    for j in eachindex(unit.exact)
        check_budget()
        sign=isodd(count_ones(state&unit.parity[j])) ? -unit.signs[j] : unit.signs[j]
        term=sign==1 ? unit.exact[j] : _qneg(unit.exact[j])
        value=_qadd(value,term)
    end
    _cache_admit!(cache,_cache_entry_bytes(value)) && (cache.amplitudes[key]=value)
    value
end
function _exact_product3(cache::AmplitudeCache{ExactScalar},a::Int,b::Int,c::Int,state::UInt64)
    vc=_amplitude(cache,c,state)
    iszero(vc) && return EXACT_ZERO
    sb=state⊻cache.units[c].flip
    vb=_amplitude(cache,b,sb)
    iszero(vb) && return EXACT_ZERO
    va=_amplitude(cache,a,sb⊻cache.units[b].flip)
    _qmul(_qmul(vc,vb),va)
end
function _exact_accumulate!(target::ExactColumn,state::UInt64,value::ExactScalar,maximum::Int)
    check_budget();_check_rational_bits(value)
    iszero(value) && return target
    updated=_qadd(get(target,state,EXACT_ZERO),value)
    if iszero(updated)
        delete!(target,state)
    else
        !haskey(target,state) && length(target)>=maximum && throw(ArgumentError("exact response support cap"))
        target[state]=updated
    end
    target
end
function _exact_add!(target::ExactColumn,source::ExactColumn,maximum::Int;subtract::Bool=false)
    for (state,value) in source
        check_budget()
        _exact_accumulate!(target,state,subtract ? _qneg(value) : value,maximum)
    end
    target
end
function _atom_column(cache::AmplitudeCache{ExactScalar},order::Tuple,hf::UInt64,maximum::Int;guard::Bool=true)
    guard || throw(ArgumentError("floating atoms require an explicit guard=false floating cache"))
    check_budget(:atom)
    length(unique(order))==length(order) && all(i->1<=i<=length(cache.units),order) ||
        throw(ArgumentError("exact atom labels must be distinct and in range"))
    result=ExactColumn()
    p(a,b,c)=_exact_product3(cache,a,b,c,hf)
    if length(order)==2
        a,b=order
        aa_b=_qadd(_qsub(p(a,a,b),_qscale(p(a,b,a),2)),p(b,a,a))
        b_ab=_qsub(_qsub(_qscale(p(b,a,b),2),p(b,b,a)),p(a,b,b))
        _exact_accumulate!(result,hf⊻cache.units[b].flip,_qscale(aa_b,-1,24),maximum)
        _exact_accumulate!(result,hf⊻cache.units[a].flip,_qscale(b_ab,-1,12),maximum)
    elseif length(order)==3
        a,b,c=order
        # Six distinct monomials, combining the two repeated ones exactly.
        positive=_qscale(_qadd(p(b,a,c),p(c,a,b)),2)
        negative=_qadd(_qadd(p(b,c,a),p(a,c,b)),_qadd(p(c,b,a),p(a,b,c)))
        value=_qscale(_qsub(positive,negative),-1,12)
        _exact_accumulate!(result,hf⊻cache.units[a].flip⊻cache.units[b].flip⊻cache.units[c].flip,value,maximum)
    else
        throw(ArgumentError("an atom has two or three distinct labels"))
    end
    result # Exact support includes HF until the final projection boundary.
end

function _certified_real_column(column::ExactColumn,maximum::Int)
    result=RealColumn()
    for (state,value) in column
        check_budget()
        _accumulate!(result,state,_rounded_exact_action(value),maximum)
    end
    result
end

function _local_kernel(cache::AmplitudeCache{ExactScalar},state::UInt64,hf::UInt64,time::Float64)
    check_budget()
    if isnothing(cache.kernel_reference)
        cache.kernel_reference=hf;cache.kernel_time=time
    else
        cache.kernel_reference==hf && cache.kernel_time==time || throw(ArgumentError("local kernel cache target differs"))
    end
    haskey(cache.kernels,state) && return cache.kernels[state]
    value=_kernel(_diagonal_gap(cache.units,hf,state),time)
    _cache_admit!(cache,_cache_entry_bytes(value)) && (cache.kernels[state]=value)
    value
end
function _transport_exact_delta(cache::AmplitudeCache{ExactScalar},column::ExactColumn,hf::UInt64,
                                time::Float64,maximum::Int;
                                cap_message::String="response support cap exceeded")
    isfinite(time) || throw(NumericalResolution("nonfinite exact DELTA time"))
    result=ResponseVector()
    for (state,value) in column
        check_budget()
        state==hf && continue # Project only after complete exact accumulation/subtraction.
        iszero(value) && continue
        real_value=_rounded_exact_action(value)
        transported=_local_kernel(cache,state,hf,time)*real_value
        isfinite(transported) || throw(NumericalResolution("nonfinite exact DELTA transport"))
        iszero(transported) && throw(NumericalResolution("nonzero exact DELTA erased by floating transport"))
        _accumulate!(result,state,transported,maximum;cap_message)
    end
    result
end

# This is a refusal heuristic, not an interval certificate. Computed zero atoms
# and deltas are admitted without claiming their mathematical exactness.
function _local_resolved(value, scale, guard, label)
    isfinite(value) && isfinite(scale) || throw(NumericalResolution("nonfinite $label scale"))
    guard && value > 0 && value <= RESOLUTION_RATIO * scale &&
        throw(NumericalResolution("$label lies below local floating resolution"))
    nothing
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
                     max_entries::Int=50000, guard::Bool=true)
    local_order = Tuple(Int[i for i in order])
    length(unique(local_order)) == length(local_order) || throw(ArgumentError("atom labels repeat"))
    all(i -> 1 <= i <= length(units), local_order) || throw(ArgumentError("atom label out of range"))
    isnothing(cache) || length(cache.units)==length(units) || throw(ArgumentError("atom cache unit count differs"))
    context = isnothing(cache) ? AmplitudeCache(center ? center_units(units,hf) : units; guard=guard) : cache
    column=_atom_column(context,local_order,UInt64(hf),max_entries;guard=guard)
    column isa ExactColumn ? _certified_real_column(column,max_entries) : column
end

function _check_rational_bits(value)
    check_budget()
    ndigits(abs(numerator(value));base=2)<=16384 && ndigits(denominator(value);base=2)<=16384 ||
        throw(ArgumentError("intermediate rational bit cap"))
    value
end

function _diagonal_gap(units, source::UInt64, destination::UInt64)
    gap=big(0)//big(1)
    for unit in units
        check_budget()
        iszero(unit.flip) || continue
        for j in eachindex(unit.coeff)
            check_budget()
            zs=isodd(count_ones(source&unit.parity[j])) ? -1 : 1
            zd=isodd(count_ones(destination&unit.parity[j])) ? -1 : 1
            gap=_qadd(gap,_qscale(unit.exact[j],Int((zd-zs)*unit.signs[j])))
        end
    end
    gap
end

function _kernel(gap::Rational, time::Float64)
    # Match the native exact dyadic-time/exact coefficient phase convention.
    _check_rational_bits(gap)
    phase=_qscale(_qmul(Rational{BigInt}(time),Rational{BigInt}(gap)),1,2)
    x=Float64(phase)
    isfinite(x) || throw(NumericalResolution("nonfinite diagonal phase"))
    ratio=iszero(x) ? 1.0 : sin(x)/x
    value=time*ratio*cis(x)
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
        check_budget()
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

function _freeze_deltas(scope, raw, baseline_rank, maximum; guard::Bool=true)
    baseline_key = Tuple(sort(collect(scope); by=i -> baseline_rank[i]))
    reference = raw[baseline_key]
    deltas = Dict{OrderKey,ResponseVector}()
    for (key, vector) in raw
        check_budget()
        delta = ResponseVector()
        if key != baseline_key
            _add!(delta, vector, 1.0, maximum)
            _add!(delta, reference, -1.0, maximum)
        end
        _local_resolved(_norm(delta), _norm(vector)+_norm(reference), guard, "DELTA")
        deltas[key] = delta
    end
    ResponseFactor(sort(collect(scope)), deltas)
end

"""Condition a contiguous interior window, preserving coherent exterior response.

Product v1 aggregates real pair/outside-triple columns BEFORE diagonal transport.
This defines new binary64 tables, distinct from the frozen research study.
Nonzero local atoms, aggregate columns and DELTAs have a 2^-24 cancellation
refusal relative to their operand scale. Computed zeros are admitted without
claiming mathematical exactness; this is not an
interval/relative-error certificate. Full Taylor acceptance remains mandatory.
"""
function _build_window_model_float(units::AbstractVector{Unit}, baseline_order, hf::Integer,
                            window_positions::AbstractUnitRange{<:Integer};
                            time::Real=1.0, max_entries::Int=50000,
                            max_cached_actions::Int=100000,
                            max_factor_entries::Int=250000, guard::Bool=false)
    guard && throw(ArgumentError("legacy floating builder requires guard=false"))
    order = _validate_order(baseline_order, length(units))
    2 <= length(window_positions) <= 9 || throw(ArgumentError("pilot window width must be between two and nine"))
    max_factor_entries >= 1 || throw(ArgumentError("factor entry cap must be positive"))
    lo, hi = Int(first(window_positions)), Int(last(window_positions))
    2 <= lo <= hi < length(order) || throw(ArgumentError("window must be strictly interior"))
    labels = order[lo:hi]
    baseline_rank = invperm(order)
    baseline_column, scale = _cubic_with_scale(units, order, hf; max_entries=max_entries, guard=guard)
    baseline = diagonal_response(units, baseline_column, hf; time=time, max_entries=max_entries,
                                 guard=guard, resolution_scale=scale)
    centered = center_units(units, hf)
    cache = AmplitudeCache(centered; max_entries=max_cached_actions, guard=guard)
    raw_atoms = 0
    factor_entries = 0
    factors = ResponseFactor[]
    function retain_factor!(factor)
        check_budget()
        entries = sum(length(v) for v in values(factor.variants))
        factor_entries + entries <= max_factor_entries || throw(ArgumentError("conditioned factor storage cap exceeded"))
        factor_entries += entries
        push!(factors,factor)
    end
    outside = [i for i in order if !(i in labels)]
    function atom(local_order)
        check_budget()
        raw_atoms += 1
        _atom_column(cache, local_order, UInt64(hf), max_entries; guard=guard)
    end
    transported(column) = diagonal_response(units, column, hf; time=time, max_entries=max_entries)
    for ia in 1:length(labels)-1, ib in ia+1:length(labels)
        check_budget()
        a,b = labels[ia], labels[ib]
        raw = Dict{OrderKey,ResponseVector}()
        for pair_order in _permutations2(a,b)
            check_budget()
            column = atom(pair_order)
            local_scale = _norm(column)
            for c in outside
                check_budget()
                local_order = baseline_rank[c] < lo ? (c,pair_order...) : (pair_order...,c)
                outside_column = atom(local_order)
                local_scale += _norm(outside_column)
                _add!(column, outside_column, 1.0, max_entries)
            end
            _local_resolved(_norm(column), local_scale, guard, "aggregated pair column")
            raw[pair_order] = transported(column)
        end
        retain_factor!(_freeze_deltas((a,b), raw, baseline_rank, max_entries; guard=guard))
    end
    for ia in 1:length(labels)-2, ib in ia+1:length(labels)-1, ic in ib+1:length(labels)
        check_budget()
        a,b,c = labels[ia],labels[ib],labels[ic]
        raw = Dict{OrderKey,ResponseVector}()
        for local_order in _permutations3(a,b,c)
            check_budget()
            raw[local_order] = transported(atom(local_order))
        end
        retain_factor!(_freeze_deltas((a,b,c), raw, baseline_rank, max_entries; guard=guard))
    end
    stats = (raw_atom_variants=raw_atoms, factors=length(factors),
             cached_actions=length(cache.amplitudes), cache_hits=cache.hits, cache_misses=cache.misses,
             aggregation="real_columns_before_transport_v1", local_resolution_ratio=RESOLUTION_RATIO,
             delta_entries=sum(length(v) for factor in factors for v in values(factor.variants)))
    WindowModel(collect(units), order, labels, lo:hi, UInt64(hf), Float64(time),
                baseline_column, baseline, factors, max_entries, guard, stats)
end

function _exact_pair_column(cache::AmplitudeCache{ExactScalar},pair_order::Tuple,outside,
                            baseline_rank,lo::Int,hf::UInt64,maximum::Int)
    column=_atom_column(cache,pair_order,hf,maximum)
    for c in outside
        check_budget()
        local_order=baseline_rank[c]<lo ? (c,pair_order...) : (pair_order...,c)
        _exact_add!(column,_atom_column(cache,local_order,hf,maximum),maximum)
    end
    column
end

"""Exact local differences anchored to the independent floating full Taylor.

Accumulate complete pair/exterior-triple columns and subtract their baseline
variant over original rationals BEFORE conversion/kernel. This new v2 frozen
target is not an exact global-response certificate. Full/global/native guards
remain mandatory. Explicit guard=false retains the legacy floating helper.
"""
function build_window_model(units::AbstractVector{Unit},baseline_order,hf::Integer,
                            window_positions::AbstractUnitRange{<:Integer};
                            time::Real=1.0,max_entries::Int=50000,
                            max_cached_actions::Int=100000,max_factor_entries::Int=250000,
                            guard::Bool=true)
    if !guard
        return _build_window_model_float(units,baseline_order,hf,window_positions;
            time,max_entries,max_cached_actions,max_factor_entries,guard=false)
    end
    order=_validate_order(baseline_order,length(units))
    2<=length(window_positions)<=9 || throw(ArgumentError("pilot window width must be between two and nine"))
    max_factor_entries>=1 || throw(ArgumentError("factor entry cap must be positive"))
    lo,hi=Int(first(window_positions)),Int(last(window_positions))
    2<=lo<=hi<length(order) || throw(ArgumentError("window must be strictly interior"))
    duration=Float64(time);isfinite(duration) || throw(ArgumentError("time must be finite"))
    labels=order[lo:hi];rank=invperm(order);source=UInt64(hf)
    baseline_column,scale=_cubic_with_scale(units,order,hf;max_entries,guard=true)
    baseline=diagonal_response(units,baseline_column,hf;time,max_entries,guard=true,resolution_scale=scale)
    cache=AmplitudeCache(center_units(units,hf);max_entries=max_cached_actions,guard=true)
    outside=[i for i in order if !(i in labels)]
    factors=ResponseFactor[];raw_atoms=0;factor_entries=0
    function freeze_exact_factor(scope,orders,column_for)
        baseline_key=Tuple(sort(collect(scope);by=i->rank[i]))
        reference=column_for(baseline_key)
        variants=Dict{OrderKey,ResponseVector}()
        for key in orders
            check_budget()
            vector=ResponseVector()
            if key!=baseline_key
                exact_delta=column_for(key)
                _exact_add!(exact_delta,reference,max_entries;subtract=true)
                remaining_factor=max_factor_entries-factor_entries
                remaining=min(max_entries,remaining_factor)
                cap_message=remaining_factor<max_entries ? "conditioned factor storage cap exceeded" : "response support cap exceeded"
                vector=_transport_exact_delta(cache,exact_delta,source,duration,remaining;cap_message)
            end
            factor_entries+=length(vector)
            variants[key]=vector
        end
        push!(factors,ResponseFactor(sort(collect(scope)),variants))
    end
    for ia in 1:length(labels)-1,ib in ia+1:length(labels)
        check_budget();a,b=labels[ia],labels[ib]
        function pair_column(key)
            raw_atoms+=1+length(outside)
            _exact_pair_column(cache,key,outside,rank,lo,source,max_entries)
        end
        freeze_exact_factor((a,b),_permutations2(a,b),pair_column)
    end
    for ia in 1:length(labels)-2,ib in ia+1:length(labels)-1,ic in ib+1:length(labels)
        check_budget();a,b,c=labels[ia],labels[ib],labels[ic]
        function triple_column(key)
            raw_atoms+=1
            _atom_column(cache,key,source,max_entries)
        end
        freeze_exact_factor((a,b,c),_permutations3(a,b,c),triple_column)
    end
    stats=(raw_atom_variants=raw_atoms,factors=length(factors),cached_actions=length(cache.amplitudes),
        cached_kernels=length(cache.kernels),cache_combined_entries=length(cache.amplitudes)+length(cache.kernels),
        cache_hits=cache.hits,cache_misses=cache.misses,cache_skipped=cache.skipped,
        exact_cache_accounted_bytes=cache.accounted_bytes,exact_cache_byte_limit=cache.max_bytes,
        cache_policy="allocated_limb_capacity_noinsert_v2",aggregation=LOCAL_TARGET,
        local_resolution_ratio=RESOLUTION_RATIO,delta_entries=factor_entries)
    WindowModel(collect(units),order,labels,lo:hi,source,duration,baseline_column,baseline,
        factors,max_entries,true,stats)
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
        check_budget()
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
