module PhysicalCheck

using LinearAlgebra
using ..ResponseCore: Unit, unitvalue, affordable_steps, NumericalResolution

export physical_compare, apply_word_rotation!, simulate_order, reference_sector,
       reference_state, state_leakage, SectorReference

struct SectorReference
    n::Int
    hf::UInt64
    basis::Vector{UInt64}
    eigenvalues::Vector{Float64}
    eigenvectors::Matrix{Float64}
    initial_coefficients::Vector{Float64}
    outside_coupling_max::Float64
    coupling_tolerance::Float64
    hermiticity_residual::Float64
end

function _validate(units::AbstractVector{Unit}, n::Integer, hf::Integer)
    1 <= n <= 12 || throw(ArgumentError("full-space circuit pilot supports 1:12 qubits"))
    0 <= hf < (1 << n) || throw(ArgumentError("HF occupation outside register"))
    isempty(units) && throw(ArgumentError("at least one physical unit required"))
    for unit in units
        isempty(unit.coeff) && throw(ArgumentError("physical units must retain nonempty ORIGINAL word lists"))
        unit.flip < (UInt64(1) << n) || throw(ArgumentError("flip outside register"))
        for parity in unit.parity
            parity < (UInt64(1) << n) || throw(ArgumentError("parity outside register"))
            iseven(count_ones(unit.flip & parity)) ||
                throw(ArgumentError("real-word phase does not define a Hermitian Pauli word"))
        end
        for i in eachindex(unit.parity), j in 1:i-1
            iseven(count_ones(unit.flip & (unit.parity[i] ⊻ unit.parity[j]))) ||
                throw(ArgumentError("physical unit words do not commute"))
        end
    end
    nothing
end

function _order(order, m)
    result=Int[order...]
    sort(result)==collect(1:m) || throw(ArgumentError("invalid complete physical order"))
    result
end

"""Apply exp(-i angle P) to the FULL state using distinct preallocated buffers.

The Pauli word is represented by flip, parity and its real Y phase sign.
Every destination is written exactly once. No particle-sector projection occurs.
"""
function apply_word_rotation!(output::Vector{ComplexF64}, input::Vector{ComplexF64},
                              flip::UInt64, parity::UInt64, sign::Integer, angle::Real)
    output !== input || throw(ArgumentError("rotation buffers must be distinct"))
    length(output)==length(input) || throw(ArgumentError("rotation buffer size mismatch"))
    sign in (-1,1) || throw(ArgumentError("Pauli sign must be +/-1"))
    flip < UInt64(length(input)) || throw(ArgumentError("rotation flip outside state"))
    iseven(count_ones(flip & parity)) || throw(ArgumentError("non-Hermitian real Pauli encoding"))
    theta=Float64(angle)
    isfinite(theta) || throw(NumericalResolution("nonfinite physical rotation angle"))
    sine,cosine=sincos(theta)
    @inbounds for at in eachindex(input)
        state=UInt64(at-1)
        destination=Int(state ⊻ flip)+1
        phase=isodd(count_ones(state & parity)) ? -sign : sign
        output[destination]=cosine*input[destination]-im*(sine*phase)*input[at]
    end
    output
end

function _schedule(order::Vector{Int}, steps::Int, formula::Int)
    schedule=Tuple{Int,Float64}[]
    function append_unit(index,weight)
        if !isempty(schedule) && last(schedule)[1]==index
            previous=pop!(schedule)
            push!(schedule,(index,previous[2]+weight))
        else
            push!(schedule,(index,weight))
        end
    end
    for _ in 1:steps
        if formula==1
            for index in order; append_unit(index,1.0); end
        else
            for index in order[1:end-1]; append_unit(index,0.5); end
            append_unit(last(order),1.0)
            for index in reverse(order[1:end-1]); append_unit(index,0.5); end
        end
    end
    schedule
end

function _step_count(units,order,budget,formula)
    budget>=0 || throw(ArgumentError("negative physical rotation budget"))
    total=sum(length(unit.coeff) for unit in units)
    total>0 || throw(ArgumentError("physical units have no words"))
    formula in (1,2) || throw(ArgumentError("formula must be 1 or 2"))
    if formula==2
        return Int(affordable_steps(units,order,budget))
    end
    steps=length(units)==1 ? (budget>=total ? 1 : 0) : fld(budget,total)
    steps>=1 || throw(ArgumentError("budget cannot afford a physical step"))
    Int(steps)
end

"""Ideal Pauli-rotation circuit on full 2^n space, with exact unit fusion.

Unit lengths are ORIGINAL implemented word counts. Only consecutive repetitions
of the same unit label are fused. Scalar words count as implemented rotations.
All norm_drift fields mean abs(sum(abs2,state)-1), the squared-norm residual.
"""
function simulate_order(units::AbstractVector{Unit}, n::Integer, hf::Integer,
                         order, time::Real, budget::Integer;
                         formula::Int=2, norm_tolerance::Real=1e-9,
                         buffers::Union{Nothing,Tuple}=nothing)
    _validate(units,n,hf)
    permutation=_order(order,length(units))
    isfinite(time) || throw(ArgumentError("time must be finite"))
    steps=_step_count(units,permutation,budget,formula)
    dimension=1 << n
    state,scratch=isnothing(buffers) ? (zeros(ComplexF64,dimension),zeros(ComplexF64,dimension)) : buffers
    length(state)==dimension && length(scratch)==dimension || throw(ArgumentError("full-state buffer dimension mismatch"))
    state !== scratch || throw(ArgumentError("physical buffers alias"))
    fill!(state,0.0im);fill!(scratch,0.0im);state[Int(hf)+1]=1.0+0.0im
    schedule=_schedule(permutation,steps,formula)
    rotations=0
    for (index,weight) in schedule
        unit=units[index]
        for j in eachindex(unit.coeff)
            angle=Float64(time)/steps*weight*unit.coeff[j]
            apply_word_rotation!(scratch,state,unit.flip,unit.parity[j],unit.signs[j],angle)
            state,scratch=scratch,state
            rotations+=1
        end
    end
    total=sum(length(unit.coeff) for unit in units)
    expected=length(units)==1 ? total : formula==1 ? steps*total :
        steps*(2total-length(units[first(permutation)].coeff)-length(units[last(permutation)].coeff))+
        length(units[first(permutation)].coeff)
    rotations==expected && rotations<=budget || error("physical fused rotation accounting differs")
    norm2=sum(abs2,state)
    drift=abs(norm2-1.0)
    isfinite(norm2) && drift<=norm_tolerance || throw(NumericalResolution("gross physical circuit norm drift"))
    (state=state,steps=steps,ideal_rotations=rotations,norm_drift=drift,
     formula=formula,budget=budget)
end

"""Dense numerical exact-H reference only within a checked fixed-N sector.

The circuit simulator does NOT use this restriction. Report the measured maximum
H coupling to the full complement and the declared numerical check threshold.
"""
function reference_sector(units::AbstractVector{Unit}, n::Integer, hf::Integer;
                          max_sector_dim::Int=2000, sector_tolerance::Real=1e-12)
    _validate(units,n,hf)
    sector_tolerance>=0 || throw(ArgumentError("negative sector check tolerance"))
    particles=count_ones(UInt64(hf))
    count=binomial(n,particles)
    count<=max_sector_dim || throw(ArgumentError("reference sector dimension exceeds cap"))
    basis=UInt64[state for state in UInt64(0):UInt64((1<<n)-1) if count_ones(state)==particles]
    index=Dict(state=>i for (i,state) in enumerate(basis))
    hamiltonian=zeros(Float64,length(basis),length(basis))
    outside_max=0.0
    coefficient_scale=sum(abs(c) for unit in units for c in unit.coeff)
    threshold=Float64(sector_tolerance)*max(1.0,coefficient_scale)
    for (column,state) in enumerate(basis)
        outside=Dict{UInt64,Float64}()
        for unit in units
            value=unitvalue(unit,state)
            destination=state ⊻ unit.flip
            if haskey(index,destination)
                hamiltonian[index[destination],column]+=value
            else
                outside[destination]=get(outside,destination,0.0)+value
            end
        end
        for value in values(outside)
            outside_max=max(outside_max,abs(value))
        end
    end
    outside_max<=threshold || throw(ArgumentError("Hamiltonian leaves the proposed fixed-N reference sector"))
    all(isfinite,hamiltonian) || throw(NumericalResolution("nonfinite reference Hamiltonian"))
    asymmetry=maximum(abs,hamiltonian-transpose(hamiltonian))
    asymmetry<=threshold || throw(NumericalResolution("reference Hamiltonian fails Hermiticity check"))
    # Symmetrization is bounded by the recorded antisymmetric residual.
    hamiltonian=(hamiltonian+transpose(hamiltonian))/2
    decomposition=eigen(Symmetric(hamiltonian))
    coefficients=collect(decomposition.vectors[index[UInt64(hf)],:])
    SectorReference(Int(n),UInt64(hf),basis,decomposition.values,decomposition.vectors,
                    coefficients,outside_max,threshold,asymmetry)
end

function reference_state(reference::SectorReference,time::Real;norm_tolerance::Real=1e-9)
    t=Float64(time)
    isfinite(t) || throw(ArgumentError("reference time must be finite"))
    weights=reference.initial_coefficients .* cis.(-t .* reference.eigenvalues)
    compact=reference.eigenvectors*weights
    full=zeros(ComplexF64,1 << reference.n)
    for (i,state) in enumerate(reference.basis)
        full[Int(state)+1]=compact[i]
    end
    drift=abs(sum(abs2,full)-1.0)
    isfinite(drift) && drift<=norm_tolerance || throw(NumericalResolution("gross reference norm drift"))
    (state=full,norm_drift=drift)
end

function _alpha_mask(n)
    result=UInt64(0)
    for site in 0:2:n-1;result|=UInt64(1)<<site;end
    result
end

"""Normalized full-state probabilities outside initial N and joint spin counts."""
function state_leakage(state::Vector{ComplexF64},n::Integer,hf::Integer;
                        alpha_mask::UInt64=_alpha_mask(n))
    length(state)==(1<<n) || throw(ArgumentError("leakage state dimension mismatch"))
    register=(UInt64(1)<<n)-1
    iszero(alpha_mask & ~register) || throw(ArgumentError("spin mask outside register"))
    beta_mask=register ⊻ alpha_mask
    particles=count_ones(UInt64(hf));alpha=count_ones(UInt64(hf)&alpha_mask)
    beta=count_ones(UInt64(hf)&beta_mask)
    total=0.0;number=0.0;spin=0.0
    for at in eachindex(state)
        occupation=UInt64(at-1);weight=abs2(state[at]);total+=weight
        count_ones(occupation)==particles || (number+=weight)
        (count_ones(occupation&alpha_mask)==alpha && count_ones(occupation&beta_mask)==beta) || (spin+=weight)
    end
    isfinite(total) && total>0 || throw(NumericalResolution("invalid state norm in leakage calculation"))
    (number_leakage=number/total,spin_leakage=spin/total)
end

function _projective_error(reference,state)
    refnorm2=sum(abs2,reference);norm2=sum(abs2,state)
    overlap=dot(reference,state)
    coefficient=overlap/refnorm2
    residual2=0.0
    for at in eachindex(state)
        residual2+=abs2(state[at]-coefficient*reference[at])
    end
    overlap_abs=abs(overlap)/sqrt(refnorm2*norm2)
    isfinite(overlap_abs) || throw(NumericalResolution("invalid physical overlap"))
    overlap_abs<=1+1e-10 || throw(NumericalResolution("physical overlap exceeds normalized bound"))
    infidelity=residual2/norm2
    isfinite(infidelity) && infidelity<=1+1e-10 || throw(NumericalResolution("invalid projective residual"))
    error=infidelity/(1+min(1.0,overlap_abs))
    (projective_error=error,infidelity=infidelity,overlap_abs=min(1.0,overlap_abs))
end

"""Compare ideal FULL-SPACE circuits with a checked exact-H sector reference.

Default spin assignment is OpenFermion's interleaved even-index alpha mask.
References are shared only within this call and time; every circuit is freshly
simulated from HF with two full-space buffers. No claim of hardware compilation
or interval certification accompanies these numerical state checks.
"""
function physical_compare(units::AbstractVector{Unit},n::Integer,hf::Integer,
                           orders::AbstractDict{String,<:AbstractVector{Int}};
                           times=[1.0],budget_multipliers=[10,30],formulas=[1,2],
                           alpha_mask::UInt64=_alpha_mask(n),max_sector_dim::Int=2000,
                           norm_tolerance::Real=1e-9,sector_tolerance::Real=1e-12)
    reference=reference_sector(units,n,hf;max_sector_dim=max_sector_dim,sector_tolerance=sector_tolerance)
    total=sum(length(unit.coeff) for unit in units)
    buffers=(zeros(ComplexF64,1<<n),zeros(ComplexF64,1<<n))
    records=NamedTuple[]
    for time in times
        exact=reference_state(reference,time;norm_tolerance=norm_tolerance)
        refleak=state_leakage(exact.state,n,hf;alpha_mask=alpha_mask)
        for name in sort!(collect(keys(orders))), multiplier in budget_multipliers, formula in formulas
            multiplier isa Integer && multiplier>=1 || throw(ArgumentError("budget multipliers must be positive integers"))
            budget=Int(multiplier)*total
            order=_order(orders[name],length(units))
            circuit=simulate_order(units,n,hf,order,time,budget;formula=Int(formula),
                                   norm_tolerance=norm_tolerance,buffers=buffers)
            metric=_projective_error(exact.state,circuit.state)
            leakage=state_leakage(circuit.state,n,hf;alpha_mask=alpha_mask)
            push!(records,(name=name,order=order,time=Float64(time),formula=Int(formula),
                budget_multiplier=Int(multiplier),rotation_budget=budget,steps=circuit.steps,
                ideal_rotations=circuit.ideal_rotations,projective_error=metric.projective_error,
                infidelity=metric.infidelity,overlap_abs=metric.overlap_abs,
                norm_drift=circuit.norm_drift,reference_norm_drift=exact.norm_drift,
                number_leakage=leakage.number_leakage,spin_leakage=leakage.spin_leakage,
                reference_number_leakage=refleak.number_leakage,
                reference_spin_leakage=refleak.spin_leakage,alpha_mask=string(alpha_mask),
                reference_sector_dimension=length(reference.basis),
                reference_outside_coupling_max=reference.outside_coupling_max,
                reference_coupling_tolerance=reference.coupling_tolerance,
                reference_hermiticity_residual=reference.hermiticity_residual))
        end
    end
    records
end

end # module
