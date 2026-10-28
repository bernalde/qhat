module ExactResponse
using TOML, SHA
using ..BoundsCommon
export ExactUnit, ExactModel, read_coefficients, freeze_kernel, load_exact_model,
       exact_full_response, project_full_atoms, endpoint_steps, maximum_steps,
       full_ratio, atom_integer, original_lengths, full_binding

struct ExactUnit
    flip::UInt64
    parity::Vector{UInt64}
    signs::Vector{Int8}
    numerators::Vector{BigInt}
end
struct ExactModel
    n::Int
    hf::UInt64
    units::Vector{ExactUnit}
    coefficient_denominator::BigInt
    kernel_denominator::BigInt
    kernel::Vector{Tuple{BigInt,BigInt}}
    input_sha256::String
    kernel_sha256::String
    time::Float64
end
original_lengths(m::ExactModel)=[length(u.numerators) for u in m.units]
function read_coefficients(path;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard)
    filesize(path)<=caps.max_input_bytes || throw(Refusal("input_limit","coefficient input exceeds byte cap"))
    raw=read(path);source=TOML.parse(String(copy(raw)));inputhash=bytes2hex(sha256(raw));n=Int(source["n"])
    1<=n<=caps.max_n && n<=20 || throw(Refusal("register_limit","register outside exact pilot"))
    blocks=source["units"];1<=length(blocks)<=caps.max_units || throw(Refusal("unit_limit","unit count outside cap"))
    DA=bigparse(source["denominator"],caps);DA>0 || throw(ArgumentError("coefficient denominator must be positive"))
    totalwords=sum(length(x["words"]) for x in blocks)
    totalwords<=caps.max_words || throw(Refusal("word_limit","original word count exceeds cap"))
    reserve!(b,totalwords;bytes=128length(blocks))
    units=ExactUnit[]
    for block in blocks
        words=block["words"];nums=block["numerators"]
        !isempty(words) && length(words)==length(nums) || throw(ArgumentError("invalid unit word inventory"))
        flip=nothing;parity=UInt64[];signs=Int8[];coeff=BigInt[]
        for (word,num) in zip(words,nums)
            tick!(b);length(word)<=4n || throw(ArgumentError("word text too long"))
            x=z=UInt64(0);ys=0;seen=Set{Int}()
            for entry in (isempty(word) ? String[] : split(word,','))
                q=parse(Int,entry[1:end-1]);p=entry[end]
                0<=q<n && !(q in seen) && p in ('X','Y','Z') || throw(ArgumentError("invalid Pauli word"))
                push!(seen,q);p in ('X','Y') && (x|=UInt64(1)<<q);p in ('Z','Y') && (z|=UInt64(1)<<q)
                ys+=p=='Y'
            end
            iseven(ys) || throw(ArgumentError("exact pilot requires real Hermitian words"))
            isnothing(flip) ? (flip=x) : (flip==x || throw(ArgumentError("unit has mixed flips")))
            a=bigparse(num,caps);!iszero(a) || throw(ArgumentError("zero ORIGINAL word coefficient"))
            push!(coeff,a);push!(parity,z);push!(signs,isodd(ys÷2) ? Int8(-1) : Int8(1))
        end
        for i in eachindex(parity),j in 1:i-1
            tick!(b)
            iseven(count_ones(flip & (parity[i]⊻parity[j]))) || throw(ArgumentError("noncommuting words in a unit"))
        end
        push!(units,ExactUnit(flip,parity,signs,coeff))
    end
    electrons=Int(source["electrons"]);0<=electrons<=n || throw(ArgumentError("invalid electron count"))
    declared=haskey(source,"hf") ? bigparse(source["hf"],caps) : (big(1)<<electrons)-1
    0<=declared<(big(1)<<n) || throw(ArgumentError("HF occupation outside register"))
    hf=UInt64(declared);count_ones(hf)==electrons || throw(ArgumentError("HF electron count differs"))
    (n=n,hf=hf,units=units,denominator=DA,source=source,input_sha256=inputhash)
end
function _unit_value(u::ExactUnit,state::UInt64,b)
    result=big(0)
    for (par,sign,a) in zip(u.parity,u.signs,u.numerators)
        tick!(b);s=isodd(count_ones(state&par)) ? -sign : sign
        result=addchecked(result,s==1 ? a : -a,b)
    end
    result
end
function freeze_kernel(input,output;time::Real=1.0,caps=Caps(),guard=()->nothing)
    ispath(output) && error("refuse to overwrite supplied kernel")
    data=read_coefficients(input;caps,guard);b=Budget(caps;guard);t=Float64(time)
    isfinite(t) || throw(ArgumentError("nonfinite time"))
    states=1<<data.n;states<=1<<caps.max_n || throw(Refusal("kernel_limit","kernel table outside cap"))
    reserve!(b,2states);floatvalues=ComplexF64[];DK=big(1)
    for occupation in 0:states-1
        target=UInt64(occupation);gap=big(0)
        for u in data.units
            iszero(u.flip) || continue
            gap=addchecked(gap,_unit_value(u,target,b)-_unit_value(u,data.hf,b),b)
        end
        x=t*Float64(gap//data.denominator)/2
        isfinite(x) || throw(ArgumentError("kernel phase outside Float64 range"))
        z=t*(iszero(x) ? 1.0 : sin(x)/x)*cis(x)
        isfinite(z) || throw(ArgumentError("nonfinite supplied kernel construction"))
        push!(floatvalues,z)
        DK=max(DK,denominator(Rational{BigInt}(real(z))),denominator(Rational{BigInt}(imag(z))))
        checked(DK,b)
    end
    realnums=BigInt[];imnums=BigInt[]
    for z in floatvalues
        a=Rational{BigInt}(real(z));c=Rational{BigInt}(imag(z))
        push!(realnums,mulchecked(numerator(a),div(DK,denominator(a)),b))
        push!(imnums,mulchecked(numerator(c),div(DK,denominator(c)),b))
    end
    record=Dict("schema"=>"supplied-dyadic-kernel-v1","input_sha256"=>data.input_sha256,
        "n"=>data.n,"hf"=>string(data.hf),"time"=>t,"denominator"=>DK,
        "real"=>realnums,"imag"=>imnums,
        "meaning"=>"Immutable supplied dyadic values; Float64 generation is not an exact-transcendental certificate.")
    write_record(output,record;caps);record
end
function load_exact_model(input,kernelpath;caps=Caps(),guard=()->nothing)
    data=read_coefficients(input;caps,guard);b=Budget(caps;guard)
    filesize(kernelpath)<=caps.max_input_bytes || throw(Refusal("input_limit","kernel input exceeds cap"))
    raw=read(kernelpath);k=TOML.parse(String(copy(raw)));kernelhash=bytes2hex(sha256(raw))
    k["schema"]=="supplied-dyadic-kernel-v1" && k["input_sha256"]==data.input_sha256 || throw(ArgumentError("kernel coefficient binding differs"))
    Int(k["n"])==data.n && bigparse(k["hf"],caps)==data.hf || throw(ArgumentError("kernel register/HF differs"))
    DK=bigparse(k["denominator"],caps)
    DK>0 && ispow2(DK) || throw(ArgumentError("supplied kernel denominator must be positive dyadic"))
    length(k["real"])==length(k["imag"])==1<<data.n || throw(ArgumentError("kernel is not a total register table"))
    reserve!(b,2length(k["real"]))
    kernel=Tuple{BigInt,BigInt}[(bigparse(a,caps),bigparse(c,caps)) for (a,c) in zip(k["real"],k["imag"])]
    t=Float64(k["time"]);isfinite(t) || throw(ArgumentError("invalid kernel provenance time"))
    DA=data.denominator
    E=mulchecked(big(48),mulchecked(mulchecked(DA,DA,b),DA,b),b);E=mulchecked(E,DK,b);mulchecked(E,E,b)
    ExactModel(data.n,data.hf,data.units,DA,DK,kernel,data.input_sha256,kernelhash,t)
end
function full_binding(m::ExactModel,targets)
    Dict{String,Any}("input_sha256"=>m.input_sha256,"kernel_sha256"=>m.kernel_sha256,
        "n"=>m.n,"hf"=>string(m.hf),"time"=>m.time,"targets"=>string.(targets),
        "original_lengths"=>original_lengths(m),"target"=>"exact supplied coefficients and immutable supplied dyadic kernel")
end

mutable struct ActionCache
    values::Dict{Tuple{Int,UInt64},BigInt}
end
ActionCache()=ActionCache(Dict{Tuple{Int,UInt64},BigInt}())
function _action(m,i,state,b,cache)
    key=(i,state)
    haskey(cache.values,key) && return cache.values[key]
    value=_unit_value(m.units[i],state,b)
    slot=256+4cld(b.caps.max_bits,8)
    # Caching is optional acceleration, never a reason to omit a unit action.
    # Reserve at most a quarter of the component allowance for cached actions.
    if length(cache.values)<b.caps.max_cache && (length(cache.values)+1)*slot<=b.caps.max_bytes÷4 && b.bytes+slot<=b.caps.max_bytes
        reserve!(b,1);cache.values[key]=value
    end
    value
end
function _add!(out,state,value,b)
    tick!(b);iszero(value) && return
    if !haskey(out,state)
        length(out)<b.caps.max_states || throw(Refusal("support_limit","exact sparse support exceeded cap"))
    end
    v=addchecked(get(out,state,big(0)),value,b)
    iszero(v) ? delete!(out,state) : (out[state]=v)
end
function _apply(m,i,v,b,cache)
    out=Dict{UInt64,BigInt}()
    for (state,a) in v
        tick!(b);_add!(out,state⊻m.units[i].flip,mulchecked(a,_action(m,i,state,b,cache),b),b)
    end
    # Six simultaneously live degree/action dictionaries plus model/cache.
    estimate=b.bytes+12length(out)*(128+cld(b.caps.max_bits,8))
    estimate<=b.caps.max_bytes || throw(Refusal("memory_limit","Taylor scratch forecast exceeded cap"))
    b.peak_bytes=max(b.peak_bytes,estimate);out
end
function _H(m,v,b,cache)
    out=Dict{UInt64,BigInt}()
    for i in eachindex(m.units),(state,a) in _apply(m,i,v,b,cache);_add!(out,state,a,b);end
    estimate=b.bytes+12length(out)*(128+cld(b.caps.max_bits,8))
    estimate<=b.caps.max_bytes || throw(Refusal("memory_limit","Hamiltonian scratch forecast exceeded cap"))
    b.peak_bytes=max(b.peak_bytes,estimate)
    out
end
function _response_scale(m,b)
    DA=m.coefficient_denominator
    mulchecked(mulchecked(big(48),mulchecked(mulchecked(DA,DA,b),DA,b),b),m.kernel_denominator,b)
end
function exact_full_response(m::ExactModel,order;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);p=Int[order...];sort(p)==collect(eachindex(m.units)) || throw(ArgumentError("not a full permutation"))
    # Model/kernel and action-cache residency are charged before Taylor work.
    reserve!(b,sum(original_lengths(m))+2length(m.kernel))
    cache=ActionCache();initial=Dict(m.hf=>big(1))
    exact=initial
    for _ in 1:3;exact=_H(m,exact,b,cache);end
    q=[initial,Dict{UInt64,BigInt}(),Dict{UInt64,BigInt}(),Dict{UInt64,BigInt}()]
    stages=[(i,1) for i in p[1:end-1]];push!(stages,(last(p),2));append!(stages,[(i,1) for i in reverse(p[1:end-1])])
    for (i,c) in stages
        b.guard();out=[copy(v) for v in q]
        for olddegree in 0:2
            power=q[olddegree+1]
            for exponent in 1:3-olddegree
                power=_apply(m,i,power,b,cache);j=olddegree+exponent
                multiplier=big(binomial(j,exponent)*c^exponent)
                for (state,a) in power;_add!(out[j+1],state,mulchecked(multiplier,a,b),b);end
            end
        end
        q=out
    end
    column=q[4]
    for (state,a) in exact;_add!(column,state,mulchecked(big(-8),a,b),b);end
    response=Dict{UInt64,Tuple{BigInt,BigInt}}();N=big(0)
    for (state,a) in column
        state==m.hf && continue
        re,im=m.kernel[Int(state)+1];x=mulchecked(a,re,b);y=mulchecked(a,im,b)
        iszero(x)&&iszero(y) && continue
        response[state]=(x,y);N=addchecked(N,addchecked(mulchecked(x,x,b),mulchecked(y,y,b),b),b)
    end
    E=_response_scale(m,b);Q=N//mulchecked(E,E,b)
    (Q=Q,numerator=N,denominator=E^2,response_denominator=E,response=response,column=column,
     order=p,work=b.work,estimated_peak_bytes=b.peak_bytes,cached_actions=length(cache.values))
end
function _product(m,a,c,d,b,cache)
    x=_action(m,d,m.hf,b,cache);state=m.hf⊻m.units[d].flip
    x=mulchecked(x,_action(m,c,state,b,cache),b);state⊻=m.units[c].flip
    mulchecked(x,_action(m,a,state,b,cache),b)
end
"""RAW cubic atom column at scale48DA^3; no kernel or HF projection."""
function atom_integer(m::ExactModel,order,b::Budget,cache=ActionCache())
    out=Dict{UInt64,BigInt}();p(a,c,d)=_product(m,a,c,d,b,cache)
    if length(order)==2
        a,c=order
        u=checked(p(a,a,c)-2p(a,c,a)+p(c,a,a),b)
        v=checked(2p(c,a,c)-p(c,c,a)-p(a,c,c),b)
        _add!(out,m.hf⊻m.units[c].flip,mulchecked(big(-2),u,b),b)
        _add!(out,m.hf⊻m.units[a].flip,mulchecked(big(-4),v,b),b)
    elseif length(order)==3
        a,c,d=order
        bac=checked(p(c,a,d)-p(c,d,a)-p(a,d,c)+p(d,a,c),b)
        cab=checked(p(d,a,c)-p(d,c,a)-p(a,c,d)+p(c,a,d),b)
        destination=m.hf⊻m.units[a].flip⊻m.units[c].flip⊻m.units[d].flip
        _add!(out,destination,mulchecked(big(-4),checked(bac+cab,b),b),b)
    else
        throw(ArgumentError("atom must have two/three labels"))
    end
    out
end
function project_full_atoms(m::ExactModel,targets;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);dest=sort!(UInt64[targets...]);s=length(dest)
    1<=s<=8 && length(unique(dest))==s || throw(ArgumentError("select one to eight distinct targets"))
    all(x->x!=m.hf && x<(UInt64(1)<<m.n),dest) || throw(ArgumentError("invalid projection target"))
    reserve!(b,sum(original_lengths(m))+2length(m.kernel))
    buckets=Dict{UInt64,Vector{Int}}()
    for (i,u) in enumerate(m.units);push!(get!(buckets,u.flip,Int[]),i);end
    flips=Set(d⊻m.hf for d in dest);scopes=Set{Tuple{Vararg{Int}}}();probes=0
    function retain!(scope)
        tick!(b)
        if !(scope in scopes)
            length(scopes)<caps.max_scopes || throw(Refusal("scope_limit","complete target scope inventory exceeds cap"))
            reserve!(b,0;bytes=128+16length(scope));push!(scopes,scope)
        end
    end
    countunits=length(m.units)
    for i in 1:countunits-1,j in i+1:countunits
        tick!(b)
        (m.units[i].flip in flips || m.units[j].flip in flips) && retain!((i,j))
    end
    for g in sort!(collect(flips)),i in 1:countunits-1,j in i+1:countunits
        probes+=1;probes<=caps.max_probes || throw(Refusal("probe_limit","complete target lookup exceeds cap"));tick!(b)
        needed=g⊻m.units[i].flip⊻m.units[j].flip
        for k in get(buckets,needed,Int[]);k>j && retain!((i,j,k));end
    end
    variants=sum(factorial(length(scope)) for scope in scopes;init=0)
    variants<=caps.max_variants || throw(Refusal("variant_limit","complete atom variants exceed cap"))
    reserve!(b,0;bytes=256variants+256length(scopes))
    factors=VectorFactor[];mapdest=Dict(d=>i for (i,d) in enumerate(dest));cache=ActionCache()
    for scope in sort!(collect(scopes))
        orders=permutations(collect(scope));table=Dict{Int,BigInt}[]
        for order in orders
            col=atom_integer(m,order,b,cache);v=Dict{Int,BigInt}()
            for (d,a) in col
                haskey(mapdest,d) || continue
                re,im=m.kernel[Int(d)+1];j=mapdest[d]
                x=mulchecked(a,re,b);y=mulchecked(a,im,b)
                if !iszero(x);reserve!(b,1);v[2j-1]=x;end
                if !iszero(y);reserve!(b,1);v[2j]=y;end
            end
            push!(table,v)
        end
        # Keep exact zero factors too: complete scope inventory is auditable.
        push!(factors,VectorFactor(collect(scope),orders,table))
    end
    E=_response_scale(m,b);bind=full_binding(m,dest)
    bind["coverage_scope_sha256"]=scope_digest(scopes);bind["coverage_scopes"]=length(scopes)
    response=IntegerResponse(collect(1:countunits),2s,E,Dict{Int,BigInt}(),factors,"full_order_projection",bind)
    validate_response(response,Budget(caps;guard))
    (response=response,targets=dest,scope_sha256=scope_digest(scopes),scopes=length(scopes),
     variants=variants,probes=probes,work=b.work,estimated_peak_bytes=b.peak_bytes)
end
function endpoint_steps(lengths,firstlabel,lastlabel,budget)
    B=BigInt(budget);B>=0 || throw(ArgumentError("negative budget"))
    bitlength(B)<=128 || throw(Refusal("budget_limit","ideal budget exceeds128-bit pilot cap"))
    !isempty(lengths) && all(x->x>0,lengths) || throw(ArgumentError("original units must be nonempty"))
    m=length(lengths);1<=firstlabel<=m && 1<=lastlabel<=m || throw(ArgumentError("bad endpoint"))
    m==1 && return B>=lengths[1] ? big(1) : big(0)
    firstlabel!=lastlabel || throw(ArgumentError("distinct endpoints required"))
    a=BigInt(lengths[firstlabel]);c=BigInt(lengths[lastlabel]);U=sum(BigInt.(lengths))
    max(big(0),fld(B-a,2U-a-c))
end
function maximum_steps(lengths,budget;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);m=length(lengths);m<=caps.max_units || throw(Refusal("unit_limit","endpoint pass exceeds cap"))
    !isempty(lengths) && all(x->x>0,lengths) || throw(ArgumentError("original unit counts must be positive"))
    B=BigInt(budget);B>=0 || throw(ArgumentError("negative budget"));U=sum(BigInt.(lengths))
    bitlength(B)<=128 || throw(Refusal("budget_limit","ideal budget exceeds128-bit pilot cap"))
    m==1 && return endpoint_steps(lengths,1,1,budget)
    result=big(0)
    for a in 1:m,c in 1:m
        a==c && continue;tick!(b)
        result=max(result,max(big(0),fld(B-lengths[a],2U-lengths[a]-lengths[c])))
    end
    result
end
function full_ratio(m,fullorder,Q,L,budget;caps=Caps(),guard=()->nothing)
    sort(Int[fullorder...])==collect(eachindex(m.units)) || throw(ArgumentError("not a full incumbent"))
    Q>=0 && L>=0 || throw(ArgumentError("negative norm bound"))
    lengths=original_lengths(m);r=endpoint_steps(lengths,first(fullorder),last(fullorder),budget)
    r>0 || throw(ArgumentError("incumbent cannot afford a step"))
    rmax=maximum_steps(lengths,budget;caps,guard);rmax>=r || error("endpoint bound inconsistent")
    LJ=L/rmax^4;UBJ=Q/r^4
    LJ<=UBJ || throw(ArgumentError("lower bound exceeds feasible same-model incumbent"))
    (L_J=LJ,UB_J=UBJ,ratio=LJ>0 ? UBJ/LJ : nothing,r_inc=r,r_max=rmax,
     positive_full_bound=LJ>0,desired_ratio_met=LJ>0 && UBJ<=2LJ)
end
end
