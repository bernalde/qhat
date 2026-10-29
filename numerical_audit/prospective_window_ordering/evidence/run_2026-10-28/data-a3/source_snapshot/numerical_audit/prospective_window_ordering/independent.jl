module IndependentPhysical

# Independent raw-Pauli implementation. Do not import the primary physical or
# response kernels: only the exact model serialization is shared.
using LinearAlgebra, SparseArrays, TOML, SHA
export load_model, reference_action, unfused_state, rotate_word!, word_action,
       state_metrics, read_state, affordable, check_manifest, main,
       IndependentModel, PauliWord, Refusal

const MAX_FULL_DIM=16384
const MAX_WORDS=8192
const MAX_UNITS=1024
const MAX_MODEL_BYTES=8*1024^2
const MAX_SECTOR_DIM=2000
const MAX_KRYLOV=160
const WordKey=Tuple{Vararg{Tuple{Int,Char}}}
const ARMS=("canonical","magnitude","current","first","pauli","cubic","window_sparse","window_magnitude")
const STATE_TOLERANCE=1e-10
const METRIC_TOLERANCE=1e-12
struct Refusal <: Exception
    message::String
end
Base.showerror(io::IO,e::Refusal)=print(io,e.message)
require(condition,message)=condition || throw(Refusal(message))

struct PauliWord
    key::WordKey
    x::UInt32
    z::UInt32
    phase::ComplexF64
    numerator::BigInt
    coefficient::Float64
end
struct IndependentModel
    n::Int
    electrons::Int
    hf::UInt32
    words::Vector{PauliWord}       # Numeric structured Python tuple order.
    units::Vector{Vector{Int}}    # ORIGINAL conserving labels/counts.
    partition::String
    case::String
    denominator::BigInt
    input_sha256::String          # Exact target TOML bytes, not source JSON.
    source_json_sha256::String
end

function integer(value,name,lo,hi)
    require(value isa Integer && !(value isa Bool) && lo<=value<=hi,"integer domain: "*name)
    Int(value)
end
function integer_string(value,name;positive=false)
    require(value isa String && 1<=ncodeunits(value)<=1280,"integer spelling bound: "*name)
    require(occursin(r"^[+-]?[0-9]+$",value),"integer spelling: "*name)
    result=parse(BigInt,value)
    require(ndigits(abs(result);base=2)<=4096,"integer bit bound: "*name)
    require(!positive || result>0,"positive integer: "*name)
    result
end
function bounded_toml(path;max_bytes=MAX_MODEL_BYTES)
    require(isfile(path) && filesize(path)<=max_bytes,"TOML file admission bound")
    raw=read(path)
    require(length(raw)<=max_bytes,"TOML read bound")
    # Admit only shallow scientific tables/arrays before invoking the parser.
    quote_char=UInt8(0);escaped=false;depth=0
    for byte in raw
        if escaped
            escaped=false
        elseif quote_char==UInt8('"') && byte==UInt8('\\')
            escaped=true
        elseif quote_char!=0
            byte==quote_char && (quote_char=0)
        elseif byte in (UInt8('"'),UInt8('\''))
            quote_char=byte
        elseif byte==UInt8('[') || byte==UInt8('{')
            depth+=1;require(depth<=8,"TOML nesting admission bound")
        elseif byte==UInt8(']') || byte==UInt8('}')
            depth-=1;require(depth>=0,"malformed TOML nesting")
        end
    end
    require(!occursin("\"\"\"",String(copy(raw))) && !occursin("'''",String(copy(raw))),"multiline TOML strings are outside schema")
    (data=TOML.parse(String(copy(raw))),sha256=bytes2hex(sha256(raw)))
end
function compile_word(spelling,n,numerator,denominator)
    require(spelling isa String && ncodeunits(spelling)<=4n,"Pauli spelling admission")
    entries=Tuple{Int,Char}[];x=UInt32(0);z=UInt32(0);ys=0;previous=-1
    for component in (isempty(spelling) ? String[] : split(spelling,','))
        require(occursin(r"^(0|[1-9][0-9]*)[XYZ]$",component),"canonical Pauli component")
        q=parse(Int,component[1:end-1]);letter=component[end]
        require(previous<q<n,"numeric sorted unique Pauli sites")
        previous=q;push!(entries,(q,letter));bit=UInt32(1)<<q
        letter in ('X','Y') && (x|=bit)
        letter in ('Y','Z') && (z|=bit)
        ys+=letter=='Y'
    end
    coefficient=Float64(numerator//denominator)
    require(isfinite(coefficient) && coefficient!=0,"coefficient conversion refuses zero/overflow")
    phase=(ComplexF64(1),ComplexF64(im),ComplexF64(-1),ComplexF64(-im))[mod(ys,4)+1]
    PauliWord(Tuple(entries),x,z,phase,numerator,coefficient)
end

"""Exact integer coefficients of [N,H]/i or [Nalpha,H]/i.

N=(nI-sum Z)/2. At a flipped site, -[Z,X]/2=-iY and
-[Z,Y]/2=+iX. Toggling the canonical Z mask changes X<->Y.
The common original denominator cancels out of the zero-polynomial test.
"""
function number_commutator(words;alpha=false)
    polynomial=Dict{Tuple{UInt32,UInt32},BigInt}()
    for word in words, (site,letter) in word.key
        letter in ('X','Y') || continue
        alpha && isodd(site) && continue
        key=(word.x,word.z⊻(UInt32(1)<<site))
        value=get(polynomial,key,big(0))+(letter=='X' ? -word.numerator : word.numerator)
        require(ndigits(abs(value);base=2)<=8192,"exact conservation intermediate bit cap")
        iszero(value) ? delete!(polynomial,key) : (polynomial[key]=value)
    end
    polynomial
end
function load_model(path)
    parsed=bounded_toml(path);data=parsed.data
    n=integer(data["n"],"n",1,14);electrons=integer(data["electrons"],"electrons",0,n)
    require((1<<n)<=MAX_FULL_DIM,"full-state preallocation cap")
    case=data["case"]
    require(case isa String && 1<=ncodeunits(case)<=200,"case label bound")
    source=get(data,"source_json_sha256",nothing)
    require(source isa String && occursin(r"^[0-9a-f]{64}$",source),"explicit source JSON hash required")
    partition=get(data,"partition","conserving")
    require(partition in ("conserving","pauli"),"unknown input partition")
    denominator=integer_string(data["denominator"],"denominator";positive=true)
    blocks=data["units"]
    require(blocks isa Vector && 1<=length(blocks)<=(partition=="conserving" ? MAX_UNITS : MAX_WORDS),"unit cap before allocation")
    total=0
    for block in blocks
        require(block isa AbstractDict && haskey(block,"words") && haskey(block,"numerators"),"unit schema")
        require(block["words"] isa Vector && block["numerators"] isa Vector,"unit arrays")
        count=length(block["words"])
        require(1<=count==length(block["numerators"])<=MAX_WORDS,"nonempty original word arrays")
        partition=="pauli" && require(count==1,"singleton Pauli input partition")
        total+=count;require(total<=MAX_WORDS,"word cap before compilation")
    end
    raw=PauliWord[];raw_units=Vector{Vector{WordKey}}();seen=Set{WordKey}()
    for block in blocks
        keys=WordKey[];compiled=PauliWord[]
        for (spelling,value) in zip(block["words"],block["numerators"])
            numerator=integer_string(value,"numerator")
            require(numerator!=0,"zero original coefficient refused")
            word=compile_word(spelling,n,numerator,denominator)
            require(!(word.key in seen),"global original Pauli words must be unique")
            push!(seen,word.key);push!(raw,word);push!(keys,word.key);push!(compiled,word)
        end
        for i in eachindex(compiled),j in 1:i-1
            a,b=compiled[i],compiled[j]
            require(iseven(count_ones(a.x&b.z)+count_ones(a.z&b.x)),"unit words fail internal commutation")
        end
        if partition=="conserving"
            require(isempty(number_commutator(compiled)),"original conserving unit fails [N,A]=0")
            require(isempty(number_commutator(compiled;alpha=true)),"original conserving unit fails [Nalpha,A]=0")
        end
        push!(raw_units,keys)
    end
    require(isempty(number_commutator(raw)),"Hamiltonian fails exact [N,H]=0")
    require(isempty(number_commutator(raw;alpha=true)),"Hamiltonian fails exact [Nalpha,H]=0")
    # Structured tuples compare sites numerically: e.g. 2Z precedes 10Z.
    sort!(raw;by=word->word.key)
    label=Dict(word.key=>i for (i,word) in enumerate(raw))
    units=[Int[label[key] for key in keys] for keys in raw_units]
    order=data[haskey(data,"baseline_order") ? "baseline_order" : "cubic_order"]
    require(order isa Vector && all(i->i isa Integer && !(i isa Bool),order) && sort(order)==collect(1:length(units)),"complete original input order")
    IndependentModel(n,electrons,UInt32((1<<electrons)-1),raw,units,partition,case,denominator,parsed.sha256,source)
end

function partition_units(model,partition)
    if partition=="pauli"
        return [Int[i] for i in eachindex(model.words)]
    end
    require(partition=="conserving" && model.partition=="conserving","original conserving partition unavailable")
    model.units
end
function validate_order(order,m)
    require(order isa AbstractVector && length(order)==m,"complete order cardinality")
    require(all(i->i isa Integer && !(i isa Bool),order) && sort(order)==collect(1:m),"complete 1-based order labels")
    Int.(order)
end
phase_at(word,state)=isodd(count_ones(word.z&UInt32(state))) ? -word.phase : word.phase
function state_dimension(state)
    d=length(state)
    require(1<=d<=MAX_FULL_DIM && ispow2(d),"full-state dimension must be bounded positive power of two")
    require(all(isfinite,state),"nonfinite full-state input")
    d
end

"""Own raw-Pauli action, including complex Y phases, on a full state."""
function word_action(word::PauliWord,state)
    d=state_dimension(state)
    require(word.x<d && word.z<d,"Pauli masks outside full state")
    output=zeros(ComplexF64,d)
    for source in 0:d-1
        output[Int(UInt32(source)⊻word.x)+1]=phase_at(word,source)*state[source+1]
    end
    output
end

"""Independent in-place disjoint-pair rotation; never a sector projection."""
function rotate_word!(state::Vector{ComplexF64},word::PauliWord,angle::Real;checked=false)
    d=checked ? length(state) : state_dimension(state)
    require(1<=d<=MAX_FULL_DIM && ispow2(d) && word.x<d && word.z<d,"rotation bounds before indexing")
    theta=Float64(angle);require(isfinite(theta),"finite primitive angle")
    c=cos(theta);s=sin(theta)
    if iszero(word.x)
        @inbounds for basis in 0:d-1
            state[basis+1]*=c-im*s*phase_at(word,basis)
        end
    else
        anchor=1<<trailing_zeros(word.x)
        # Exactly one representative with anchor bit zero for each XOR pair.
        @inbounds for block in 0:2anchor:d-1, offset in 0:anchor-1
            a=block+offset;b=Int(UInt32(a)⊻word.x)
            left,right=state[a+1],state[b+1]
            pa,pb=phase_at(word,a),phase_at(word,b)
            state[a+1]=c*left-im*s*pb*right
            state[b+1]=c*right-im*s*pa*left
        end
    end
    state
end

function affordable(model,partition,order,budget;formula=2)
    units=partition_units(model,partition);permutation=validate_order(order,length(units))
    cap=integer(budget,"rotation budget",0,10^9)
    formula=integer(formula,"formula",1,2)
    lengths=length.(units);total=sum(lengths);m=length(units)
    steps=m==1 ? (cap>=total ? 1 : 0) : formula==1 ? fld(cap,total) :
        fld(cap-lengths[first(permutation)],2total-lengths[first(permutation)]-lengths[last(permutation)])
    require(steps>=1,"original-word budget cannot afford a step")
    rotations=m==1 ? total : formula==1 ? steps*total :
        steps*(2total-lengths[first(permutation)]-lengths[last(permutation)])+lengths[first(permutation)]
    require(rotations<=cap,"independent production-cost check")
    (steps=steps,ideal_rotations=rotations,total_words=total,word_lengths=lengths)
end

function unfused_state(model,partition,order,T,steps;formula=2,max_primitive_pulses=2_000_000,check_budget=()->nothing)
    units=partition_units(model,partition);permutation=validate_order(order,length(units))
    r=integer(steps,"unfused steps",1,1_000_000)
    formula=integer(formula,"unfused formula",1,2)
    require(!(T isa Bool),"Boolean physical time")
    t=Float64(T);require(isfinite(t),"finite physical time")
    pulses=BigInt(r)*sum(length,units)*(formula==2 ? 2 : 1)
    require(pulses<=max_primitive_pulses,"unfused primitive work cap before allocation")
    require(1<=model.n<=14 && (1<<model.n)<=MAX_FULL_DIM,"full circuit preallocation cap")
    state=zeros(ComplexF64,1<<model.n);state[Int(model.hf)+1]=1
    count=0
    for _ in 1:r
        directions=formula==1 ? (permutation,) : (permutation,reverse(permutation))
        weight=formula==1 ? 1.0 : .5
        for direction in directions,index in direction,id in units[index]
            check_budget()
            word=model.words[id]
            rotate_word!(state,word,(t/r)*weight*word.coefficient;checked=true)
            count+=1
        end
    end
    require(count==pulses,"literal auxiliary primitive count")
    drift=abs(sum(abs2,state)-1)
    require(all(isfinite,state) && drift<=1e-9,"unfused circuit norm check")
    (state=state,primitive_pulses=count,norm_drift=drift,
     pulse_scope="literal unfused auxiliary count; not production ideal resource count")
end

function compensated(sum,correction,value)
    updated=sum+value
    correction+=abs(sum)>=abs(value) ? (sum-updated)+value : (value-updated)+sum
    (updated,correction)
end
function sector_hamiltonian(model,max_sector_dim,check_budget)
    require(1<=model.n<=14 && (1<<model.n)<=MAX_FULL_DIM,"reference register cap BEFORE allocation")
    limit=integer(max_sector_dim,"reference sector cap",1,MAX_SECTOR_DIM)
    dimension=binomial(BigInt(model.n),model.electrons)
    require(dimension<=limit,"reference sector dimension cap BEFORE basis allocation")
    basis=UInt32[b for b in 0:(1<<model.n)-1 if count_ones(UInt32(b))==model.electrons]
    require(length(basis)==dimension,"independent fixed-number basis dimension")
    index=Dict(b=>i for (i,b) in enumerate(basis))
    rows=Int[];columns=Int[];values=ComplexF64[]
    scale=sum(abs(word.coefficient) for word in model.words)
    require(isfinite(scale),"Hamiltonian coefficient norm overflow")
    outside=0.0;spin_outside=0.0;alpha=sum(UInt32(1)<<site for site in 0:2:model.n-1;init=UInt32(0))
    for (column,state) in enumerate(basis)
        check_budget()
        amplitudes=Dict{UInt32,NTuple{4,Float64}}()
        for word in model.words
            check_budget()
            destination=state⊻word.x;value=word.coefficient*phase_at(word,state)
            sr,cr,si,ci=get(amplitudes,destination,(0.,0.,0.,0.))
            sr,cr=compensated(sr,cr,real(value));si,ci=compensated(si,ci,imag(value))
            amplitudes[destination]=(sr,cr,si,ci)
        end
        for (destination,(sr,cr,si,ci)) in amplitudes
            check_budget()
            value=ComplexF64(sr+cr,si+ci)
            if !haskey(index,destination)
                outside=max(outside,abs(value))
            else
                count_ones(destination&alpha)==count_ones(state&alpha) || (spin_outside=max(spin_outside,abs(value)))
                if !iszero(value)
                    require(length(values)<2_000_000,"sector sparse storage cap BEFORE append")
                    push!(rows,index[destination]);push!(columns,column);push!(values,value)
                end
            end
        end
    end
    threshold=1e-13*max(1.,scale)
    require(outside<=threshold && spin_outside<=threshold,"numeric action contradicts exact sector conservation")
    h=sparse(rows,columns,values,length(basis),length(basis))
    difference=h-adjoint(h);asymmetry=isempty(nonzeros(difference)) ? 0.0 : maximum(abs,nonzeros(difference))
    require(isfinite(asymmetry) && asymmetry<=threshold,"independent sector Hermiticity check")
    h=(h+adjoint(h))/2
    (h=h,basis=basis,index=index,outside=outside,spin_outside=spin_outside,asymmetry=asymmetry)
end

"""Independent full-reorthogonalized Hermitian Krylov exponential action.

Only the independently conservation-checked N sector is used for the reference;
full-space circuits never use it. Convergence diagnostics are numerical, not an
interval certificate. Scalar centering is undone in the returned global phase.
"""
function reference_action(model,T;max_sector_dim=2000,max_krylov=160,tolerance=1e-13,check_budget=()->nothing)
    require(1<=model.n<=14 && (1<<model.n)<=MAX_FULL_DIM,"reference full-state cap BEFORE allocation")
    require(0<=model.hf<(1<<model.n) && count_ones(model.hf)==model.electrons,"reference basis occupation binding")
    maximum=integer(max_krylov,"Krylov dimension cap",1,MAX_KRYLOV)
    tol=Float64(tolerance);require(isfinite(tol) && 0<tol<=1e-12,"reference convergence tolerance")
    require(!(T isa Bool),"Boolean reference time")
    t=Float64(T);require(isfinite(t),"finite reference time")
    sector=sector_hamiltonian(model,max_sector_dim,check_budget)
    d=length(sector.basis);kmax=min(maximum,d)
    origin=sector.index[model.hf];shift=real(sector.h[origin,origin])
    require(isfinite(shift),"finite reference scalar shift")
    vectors=zeros(ComplexF64,d,kmax+1);vectors[origin,1]=1
    projected=zeros(ComplexF64,kmax+1,kmax)
    work=zeros(ComplexF64,d);previous=nothing;agreement=0
    answer=nothing;used=0;indicator=Inf;change=Inf;beta=Inf;orthogonality=Inf
    for j in 1:kmax
        check_budget();mul!(work,sector.h,view(vectors,:,j))
        @inbounds for i in 1:d;work[i]-=shift*vectors[i,j];end
        # Two full modified Gram-Schmidt passes, independently of the primary
        # dense-sector diagonalization. Retain the complete projected matrix.
        for _ in 1:2,i in 1:j
            check_budget();coefficient=dot(view(vectors,:,i),work)
            projected[i,j]+=coefficient
            @inbounds for at in 1:d;work[at]-=coefficient*vectors[at,i];end
        end
        beta=norm(work)
        require(isfinite(beta),"nonfinite Krylov residual")
        projected[j+1,j]=beta
        small=Matrix(projected[1:j,1:j]);asymmetry=norm(small-adjoint(small))
        require(asymmetry<=1e-11*max(1.,norm(small)),"Krylov projected Hermiticity check")
        eig=eigen(Hermitian((small+adjoint(small))/2))
        coefficients=eig.vectors*(cis.(-t.*eig.values).*conj.(eig.vectors[1,:]))
        compact=view(vectors,:,1:j)*coefficients
        indicator=0.0
        for fraction in (0.0,.25,.5,.75,1.0)
            check_budget()
            witness=eig.vectors*(cis.((-t*fraction).*eig.values).*conj.(eig.vectors[1,:]))
            indicator=max(indicator,abs(t)*beta*abs(last(witness)))
        end
        change=isnothing(previous) ? Inf : norm(compact-previous)
        agreement=change<=tol/4 ? agreement+1 : 0
        breakdown=iszero(beta) || (beta<=32eps(Float64)*max(1.,norm(small)) && indicator<=tol/4)
        if breakdown || (j>=4 && agreement>=2 && indicator<=tol/4)
            answer=compact;used=j
            orthogonality=norm(adjoint(view(vectors,:,1:j))*view(vectors,:,1:j)-I)
            break
        end
        require(j<kmax,"Krylov cap reached without declared numerical convergence")
        require(beta>0,"unresolved zero Krylov norm")
        @inbounds for at in 1:d;vectors[at,j+1]=work[at]/beta;end
        previous=compact
    end
    require(answer!==nothing && orthogonality<=1e-10,"Krylov orthogonality/convergence refusal")
    answer .*= cis(-t*shift)
    drift=abs(sum(abs2,answer)-1)
    require(all(isfinite,answer) && drift<=1e-11,"independent reference norm check")
    full=zeros(ComplexF64,1<<model.n)
    for (i,basis) in enumerate(sector.basis);full[Int(basis)+1]=answer[i];end
    leakage=state_metrics(full,full,model.n,model.hf)
    require(leakage["number_leakage"]<=1e-12 && leakage["spin_leakage"]<=1e-12,"reference number/spin leakage check")
    diagnostic=Dict("method"=>"independent raw-Pauli sparse-sector fully reorthogonalized Hermitian Krylov",
        "sector_dimension"=>d,"krylov_dimension"=>used,"max_krylov"=>maximum,"tolerance"=>tol,
        "sector_nnz"=>nnz(sector.h),"original_words"=>length(model.words),
        "sampled_residual_indicator"=>indicator,"successive_action_change"=>(isfinite(change) ? change : "unavailable"),
        "scalar_shift"=>shift,"global_phase_restored"=>true,"norm_drift"=>drift,
        "orthogonality_residual"=>orthogonality,"sector_outside_max"=>sector.outside,
        "number_leakage"=>leakage["number_leakage"],"spin_leakage"=>leakage["spin_leakage"],
        "spin_outside_max"=>sector.spin_outside,"sector_hermiticity_residual"=>sector.asymmetry,
        "certificate"=>"numerical convergence/norm checks; no interval certificate")
    (state=full,diagnostic=diagnostic)
end

function state_metrics(reference,state,n,hf)
    require(length(reference)==length(state)==(1<<n) && 1<=n<=14,"metric full-state dimension")
    require(0<=hf<(1<<n),"metric occupation outside register")
    nr=sum(abs2,reference);ns=sum(abs2,state)
    require(all(isfinite,reference) && all(isfinite,state) && isfinite(nr) && isfinite(ns) && nr>0 && ns>0,"metric finite nonzero norms")
    require(abs(ns-1)<=1e-9 && abs(nr-1)<=1e-9,"gross saved-state norm drift")
    u=reference/sqrt(nr);v=state/sqrt(ns);overlap=dot(u,v)
    magnitude=abs(overlap);require(magnitude<=1+1e-10,"normalized overlap bound")
    # Stable orthogonal residual, retaining errors below the subtraction floor.
    residual=0.0
    for i in eachindex(v);residual+=abs2(v[i]-overlap*u[i]);end
    require(isfinite(residual) && residual<=1+1e-10,"projective residual bound")
    infidelity=max(0.,min(1.,residual));projective=infidelity/(1+min(1.,magnitude))
    alpha=UInt32(0)
    for site in 0:2:n-1;alpha|=UInt32(1)<<site;end
    beta=(UInt32(1)<<n)-1-alpha
    particles=count_ones(UInt32(hf));na=count_ones(UInt32(hf)&alpha);nb=count_ones(UInt32(hf)&beta)
    number=0.;spin=0.
    for basis in 0:length(state)-1
        weight=abs2(state[basis+1])/ns
        count_ones(UInt32(basis))==particles || (number+=weight)
        (count_ones(UInt32(basis)&alpha)==na && count_ones(UInt32(basis)&beta)==nb) || (spin+=weight)
    end
    Dict("projective_error"=>projective,"infidelity"=>infidelity,"norm_drift"=>abs(ns-1),
         "number_leakage"=>number,"spin_leakage"=>spin,"reference_norm_drift"=>abs(nr-1))
end
function phase_aligned_residual(reference,state)
    require(isfinite(norm(reference)) && isfinite(norm(state)) && norm(reference)>0 && norm(state)>0,"phase comparison requires finite nonzero states")
    a=reference/norm(reference);b=state/norm(state);overlap=dot(a,b)
    phase=iszero(overlap) ? one(overlap) : overlap/abs(overlap)
    norm(b-phase*a)
end
function safe_state_path(parent,relative)
    require(relative isa String && 1<=ncodeunits(relative)<=512 && !isabspath(relative),"relative bounded state path required")
    require(!any(component->component=="..",splitpath(relative)),"state path traversal refused")
    base=realpath(parent);candidate=realpath(joinpath(base,relative))
    require(!any(component->component=="..",splitpath(relpath(candidate,base))) && isfile(candidate),"state path/symlink escapes manifest directory")
    candidate
end
function read_state(parent,relative,digest,n)
    require(1<=n<=14 && (1<<n)<=MAX_FULL_DIM,"state dimension BEFORE allocation")
    require(digest isa String && occursin(r"^[0-9a-f]{64}$",digest),"state SHA-256 spelling")
    path=safe_state_path(parent,relative)
    require(filesize(path)==16*(1<<n),"little-endian state exact byte length")
    raw=read(path);require(length(raw)==16*(1<<n),"state byte length changed during read")
    require(bytes2hex(sha256(raw))==digest,"saved state SHA-256 mismatch")
    bits=reinterpret(UInt64,raw);state=Vector{ComplexF64}(undef,1<<n)
    for i in eachindex(state)
        state[i]=ComplexF64(reinterpret(Float64,ltoh(bits[2i-1])),reinterpret(Float64,ltoh(bits[2i])))
    end
    require(all(isfinite,state),"nonfinite saved state")
    state
end

function physical_time(value)
    require(!(value isa Bool),"Boolean physical time")
    t=value isa AbstractString ? parse(Float64,value) : Float64(value)
    require(isfinite(t) && t in (.4,1.6),"frozen physical times only")
    t
end
function expected_inventory(repetitions)
    Set((arm,repetition,formula,t,multiplier) for arm in ARMS for repetition in 1:repetitions
        for formula in (arm in ("canonical","current","first") ? (1,2) : (2,))
        for t in (.4,1.6) for multiplier in (8,24))
end
function row_identity(row,repetitions)
    arm=row["arm"];require(arm in ARMS,"unknown frozen arm")
    repetition=integer(row["repetition"],"repetition",1,repetitions)
    formula=integer(row["formula"],"formula",1,2)
    require(formula==2 || arm in ("canonical","current","first"),"frozen arm/formula policy")
    t=physical_time(row["time"]);multiplier=integer(row["budget_multiplier"],"budget multiplier",8,24)
    require(multiplier in (8,24),"frozen physical budget multipliers")
    (arm,repetition,formula,t,multiplier)
end
function exception_record(error)
    message=sprint(showerror,error)
    Dict("type"=>string(typeof(error)),"message"=>first(message,min(length(message),500)))
end

function check_manifest(model,manifest_path;check_budget=()->nothing)
    parsed=bounded_toml(manifest_path);manifest=parsed.data
    report=Dict{String,Any}("status"=>"incomplete","input_toml_sha256"=>model.input_sha256,
        "source_json_sha256"=>model.source_json_sha256,"manifest_sha256"=>parsed.sha256,
        "case"=>model.case,"n"=>model.n,"hf"=>string(model.hf),"rows"=>Any[],"failures"=>Any[],
        "reference_checks"=>Any[],"replay_scope"=>"only distinct T=1.6 / 24U circuits; other rows independently check saved-state metrics",
        "state_norm_tolerance"=>STATE_TOLERANCE,"absolute_metric_tolerance"=>METRIC_TOLERANCE,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "exact_polynomial_conservation"=>"[N,H]=[Nalpha,H]=0 and every original conserving unit checked over supplied integer coefficients",
        "checker_source_sha256"=>bytes2hex(sha256(read(@__FILE__))))
    try
        require(model.partition=="conserving","checker target must retain the original conserving partition")
        require(manifest["input_sha256"]==model.input_sha256,"manifest must authenticate target TOML bytes")
        require(manifest["case"]==model.case && integer(manifest["n"],"manifest n",1,14)==model.n && manifest["hf"]==string(model.hf),"manifest model/reference binding")
        mode=manifest["mode"];require(mode in ("main","preflight"),"explicit main/preflight mode")
        repetitions=mode=="main" ? 2 : 1
        expected=expected_inventory(repetitions)
        require(integer(manifest["expected_logical_rows"],"declared rows",44,88)==length(expected),"declared complete logical inventory size")
        rows=manifest["rows"];require(rows isa Vector && length(rows)<=88,"bounded logical rows")
        observed=Set{Tuple{String,Int,Int,Float64,Int}}()
        for row in rows
            require(row isa AbstractDict && length(row)<=64,"bounded row schema")
            identity=row_identity(row,repetitions)
            require(!(identity in observed),"duplicate frozen logical identity")
            push!(observed,identity)
        end
        report["expected_logical_rows"]=length(expected);report["provided_logical_rows"]=length(rows)
        complete_roster=observed==expected
        if !complete_roster
            push!(report["failures"],Dict("type"=>"IncompleteInventory","message"=>"missing frozen logical identities", "missing_count"=>length(setdiff(expected,observed))))
        end
        if get(manifest,"status",nothing)!="complete" || get(manifest,"physical_panel_complete",nothing)!==true
            push!(report["failures"],Dict("type"=>"IncompletePrimaryPanel","message"=>"primary panel is explicitly incomplete/censored"))
        end
        references=Dict{Float64,Any}()
        for t in (.4,1.6)
            check_budget()
            try
                references[t]=reference_action(model,t;check_budget=check_budget)
            catch error
                push!(report["failures"],merge(Dict("phase"=>"independent_reference","time"=>t),exception_record(error)))
            end
        end
        saved_reference_checks=Set{Tuple{Float64,String,String}}()
        replay_cache=Dict{Any,Any}();replay_failures=Dict{Any,Any}()
        parent=dirname(abspath(manifest_path));hard_identities=Set{Any}()
        for row in rows
            check_budget()
            entry=Dict{String,Any}("status"=>"incomplete","name"=>get(row,"name","unlabelled"),
                "arm"=>row["arm"],"repetition"=>row["repetition"],"formula"=>row["formula"],
                "time"=>row["time"],"budget_multiplier"=>row["budget_multiplier"])
            push!(report["rows"],entry)
            try
                identity=row_identity(row,repetitions);arm,repetition,formula,t,multiplier=identity
                require(row["name"] isa String && 1<=ncodeunits(row["name"])<=200,"bounded row name")
                partition=row["partition"]
                require(partition==(arm=="pauli" ? "pauli" : "conserving"),"arm/partition binding")
                units=partition_units(model,partition);order=validate_order(row["order"],length(units))
                budget=integer(row["rotation_budget"],"rotation budget",0,10^9)
                require(budget==multiplier*length(model.words),"ORIGINAL total word-count budget binding")
                plan=affordable(model,partition,order,budget;formula=formula)
                require(integer(row["steps"],"saved steps",1,1_000_000)==plan.steps &&
                    integer(row["ideal_rotations"],"saved ideal rotations",1,10^9)==plan.ideal_rotations,"independent original-word affordability/count mismatch")
                require(haskey(references,t),"independent reference incomplete")
                reference=read_state(parent,row["reference_file"],row["reference_sha256"],model.n)
                actual=read_state(parent,row["actual_file"],row["actual_sha256"],model.n)
                refidentity=(t,String(row["reference_file"]),String(row["reference_sha256"]))
                if !(refidentity in saved_reference_checks)
                    independent=references[t]
                    residual=norm(reference-independent.state)
                    aligned=phase_aligned_residual(independent.state,reference)
                    check=Dict("time"=>t,"reference_file"=>row["reference_file"],"reference_sha256"=>row["reference_sha256"],
                        "state_norm_residual"=>residual,"phase_aligned_residual"=>aligned,"diagnostic"=>independent.diagnostic,
                        "status"=>(residual<=STATE_TOLERANCE ? "complete" : "failed"))
                    push!(report["reference_checks"],check)
                    require(residual<=STATE_TOLERANCE,"independent Krylov/primary reference norm mismatch")
                    push!(saved_reference_checks,refidentity)
                end
                metrics=state_metrics(reference,actual,model.n,model.hf)
                for metric in ("projective_error","infidelity","norm_drift","number_leakage","spin_leakage")
                    require(row[metric] isa Real && !(row[metric] isa Bool),"numeric saved physical metric")
                    value=Float64(row[metric]);require(isfinite(value) && value>=0,"finite saved physical metric")
                    require(abs(value-metrics[metric])<=METRIC_TOLERANCE,"independent saved-state metric differs: "*metric)
                end
                entry["metrics"]=metrics;entry["saved_actual_sha256"]=row["actual_sha256"]
                entry["independent_unfused_replay"]=false
                if t==1.6 && multiplier==24
                    # Partition is part of the meaning of the original order IDs.
                    key=(model.input_sha256,partition,Tuple(order),formula,t,plan.steps)
                    push!(hard_identities,key)
                    if !haskey(replay_cache,key) && !haskey(replay_failures,key)
                        try
                            replay_cache[key]=unfused_state(model,partition,order,t,plan.steps;formula=formula,check_budget=check_budget)
                        catch error
                            replay_failures[key]=error
                        end
                    end
                    haskey(replay_failures,key) && throw(replay_failures[key])
                    replay=replay_cache[key]
                    residual=norm(replay.state-actual)
                    require(residual<=STATE_TOLERANCE,"independent unfused/main full-state norm mismatch")
                    replay_metrics=state_metrics(reference,replay.state,model.n,model.hf)
                    for metric in ("projective_error","infidelity","norm_drift","number_leakage","spin_leakage")
                        require(abs(replay_metrics[metric]-metrics[metric])<=METRIC_TOLERANCE,"unfused/main physical metric mismatch: "*metric)
                    end
                    entry["independent_unfused_replay"]=true
                    entry["unfused_state_norm_residual"]=residual
                    entry["unfused_phase_aligned_residual"]=phase_aligned_residual(actual,replay.state)
                    entry["auxiliary_unfused_primitive_pulses"]=replay.primitive_pulses
                end
                entry["status"]="complete"
            catch error
                entry["failure"]=exception_record(error)
            end
        end
        report["distinct_hard_circuit_identities"]=length(hard_identities)
        report["distinct_successful_unfused_replays"]=length(replay_cache)
        report["checked_logical_rows"]=count(row->row["status"]=="complete",report["rows"])
        if complete_roster && isempty(report["failures"]) && report["checked_logical_rows"]==length(expected)
            report["status"]="complete"
        end
    catch error
        push!(report["failures"],exception_record(error))
    end
    report["checked_logical_rows"]=count(row->get(row,"status",nothing)=="complete",report["rows"])
    report
end

function main(args=ARGS)
    length(args)==3 || throw(ArgumentError("usage: independent.jl INPUT_TOML RESULT_MANIFEST_TOML OUTPUT_CHECK_TOML"))
    input,manifest,output=args
    require(!ispath(output) && !ispath(output*".partial"),"refuse existing independent output")
    BLAS.set_num_threads(1)
    report=try
        require(Threads.nthreads()==1,"independent checker requires one Julia numerical thread")
        model=load_model(input)
        check_manifest(model,manifest)
    catch error
        Dict{String,Any}("status"=>"incomplete","rows"=>Any[],"failure"=>exception_record(error),
            "input_target"=>basename(input),"manifest_target"=>basename(manifest))
    end
    open(output*".partial","w") do io;TOML.print(io,report;sorted=true);end
    mv(output*".partial",output;force=false)
    report["status"]=="complete" ? 0 : 2
end

end # module
if abspath(PROGRAM_FILE)==@__FILE__
    exit(IndependentPhysical.main())
end
