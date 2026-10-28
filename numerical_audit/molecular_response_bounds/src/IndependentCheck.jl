module IndependentCheck
using TOML, SHA
using ..BoundsCommon
export independent_bound, independent_full_check, independent_window_check, checker_projection

# This module deliberately does not import ExactResponse or ResponseBounds.
# Only representation, serialization and resource accounting are shared.
const Rat=Rational{BigInt}
const ZERO=big(0)//big(1)
struct CheckUnit
    flip::UInt64
    parity::Vector{UInt64}
    coefficient::Vector{Rat}  # Includes each real Y phase.
end
function _rcheck(q::Rat,b)
    checked(numerator(q),b);checked(denominator(q),b);q
end
function _plus(a::Rat,c::Rat,b)
    tick!(b);_rcheck(a+c,b)
end
function _times(a::Rat,c::Rat,b)
    tick!(b)
    bitlength(numerator(a))+bitlength(numerator(c))<=b.caps.max_bits &&
        bitlength(denominator(a))+bitlength(denominator(c))<=b.caps.max_bits || throw(Refusal("bit_limit","checker prospective rational product"))
    _rcheck(a*c,b)
end
function _checker_source(input,kernelpath,b)
    caps=b.caps
    filesize(input)<=caps.max_input_bytes && filesize(kernelpath)<=caps.max_input_bytes || throw(Refusal("input_limit","checker input exceeds cap"))
    inputbytes=read(input);kernelbytes=read(kernelpath)
    data=TOML.parse(String(copy(inputbytes)));k=TOML.parse(String(copy(kernelbytes)));n=Int(data["n"])
    1<=n<=caps.max_n || throw(Refusal("register_limit","checker register exceeds cap"))
    blocks=data["units"];1<=length(blocks)<=caps.max_units || throw(Refusal("unit_limit","checker units exceed cap"))
    den=bigparse(data["denominator"],caps);den>0 || error("nonpositive source denominator")
    words=sum(length(u["words"]) for u in blocks);words<=caps.max_words || throw(Refusal("word_limit","checker words exceed cap"))
    reserve!(b,words)
    units=CheckUnit[];lengths=Int[]
    for block in blocks
        strings=block["words"];nums=block["numerators"]
        length(strings)==length(nums)>0 || error("checker source unit cardinality")
        compiled=Tuple{UInt64,UInt64,Rat}[]
        for (word,raw) in zip(strings,nums)
            tick!(b);length(word)<=4n || error("checker Pauli string length")
            maskx=UInt64(0);maskz=UInt64(0);phase=0;used=Set{Int}()
            for token in (isempty(word) ? String[] : split(word,','))
                matchword=match(r"^([0-9]+)([XYZ])$",token);isnothing(matchword) && error("checker Pauli token")
                site=parse(Int,matchword.captures[1]);letter=only(matchword.captures[2])
                0<=site<n && !(site in used) || error("checker repeated/outside site");push!(used,site)
                letter!='Z' && (maskx⊻=UInt64(1)<<site)
                letter!='X' && (maskz⊻=UInt64(1)<<site)
                phase+=letter=='Y'
            end
            phase%2==0 || error("checker requires real words")
            numerator=bigparse(raw,caps);!iszero(numerator) || error("checker original word zero")
            signed=(phase%4==0 ? numerator : -numerator)//den
            push!(compiled,(maskx,maskz,signed))
        end
        all(t->t[1]==compiled[1][1],compiled) || error("checker equal-flip condition")
        for i in eachindex(compiled),j in 1:i-1
            tick!(b)
            iseven(count_ones(compiled[i][1]&compiled[j][2])+count_ones(compiled[j][1]&compiled[i][2])) || error("checker commuting-unit condition")
        end
        push!(units,CheckUnit(compiled[1][1],[t[2] for t in compiled],[t[3] for t in compiled]));push!(lengths,length(compiled))
    end
    electrons=Int(data["electrons"]);0<=electrons<=n || error("checker electron count")
    occupation=haskey(data,"hf") ? bigparse(data["hf"],caps) : big(2)^electrons-1
    0<=occupation<big(2)^n || error("checker HF outside register")
    hf=UInt64(occupation);count_ones(hf)==electrons || error("checker HF particle count")
    inputhash=bytes2hex(sha256(inputbytes))
    k["schema"]=="supplied-dyadic-kernel-v1" && k["input_sha256"]==inputhash || error("checker kernel/input authentication")
    Int(k["n"])==n && bigparse(k["hf"],caps)==hf || error("checker kernel HF/register")
    kd=bigparse(k["denominator"],caps);kd>0 && count_ones(kd)==1 || error("checker dyadic kernel denominator")
    length(k["real"])==length(k["imag"])==2^n || error("checker requires total kernel")
    reserve!(b,2length(k["real"]))
    kval=[(bigparse(x,caps)//kd,bigparse(y,caps)//kd) for (x,y) in zip(k["real"],k["imag"])]
    t=Float64(k["time"]);isfinite(t) || error("checker kernel provenance time")
    (n=n,hf=hf,units=units,lengths=lengths,DA=den,DK=kd,kernel=kval,time=t,
     inputhash=inputhash,kernelhash=bytes2hex(sha256(kernelbytes)))
end
function _value(u::CheckUnit,state,b)
    total=ZERO
    for j in eachindex(u.parity)
        tick!(b);v=u.coefficient[j]
        iseven(count_ones(state&u.parity[j])) || (v=-v)
        total=_plus(total,v,b)
    end
    total
end
function _put!(out,state,value,b)
    iszero(value) && return
    v=_plus(get(out,state,ZERO),value,b)
    if iszero(v);delete!(out,state)
    else
        haskey(out,state) || length(out)<b.caps.max_states || throw(Refusal("support_limit","checker support exceeded"))
        out[state]=v
    end
end
function _act(units,i,vector,b)
    out=Dict{UInt64,Rat}();u=units[i]
    for (state,coeff) in vector;_put!(out,state⊻u.flip,_times(_value(u,state,b),coeff,b),b);end
    estimate=b.bytes+12length(out)*(192+2cld(b.caps.max_bits,8))
    estimate<=b.caps.max_bytes || throw(Refusal("memory_limit","independent scratch forecast"))
    b.peak_bytes=max(b.peak_bytes,estimate);out
end
"""Unscaled rational polynomial action, separate from integer Taylor backend."""
function _polynomial_defect(data,order,b)
    sort(Int[order...])==collect(eachindex(data.units)) || error("checker full permutation")
    initial=Dict(data.hf=>big(1)//big(1))
    coefficients=[copy(initial),Dict{UInt64,Rat}(),Dict{UInt64,Rat}(),Dict{UInt64,Rat}()]
    schedule=[(x,big(1)//big(2)) for x in order[1:end-1]]
    push!(schedule,(last(order),big(1)//big(1)));append!(schedule,[(x,big(1)//big(2)) for x in reverse(order[1:end-1])])
    for (label,tau) in schedule
        b.guard();updated=[copy(v) for v in coefficients]
        for degree in 0:2
            vector=coefficients[degree+1];factor=big(1)//big(1)
            for p in 1:3-degree
                vector=_act(data.units,label,vector,b);factor=_times(factor,tau/p,b)
                for (state,a) in vector;_put!(updated[degree+p+1],state,_times(factor,a,b),b);end
            end
        end
        coefficients=updated
    end
    exact=initial
    for _ in 1:3
        total=Dict{UInt64,Rat}()
        for label in eachindex(data.units),(state,a) in _act(data.units,label,exact,b);_put!(total,state,a,b);end
        exact=total
    end
    for (state,a) in exact;_put!(coefficients[4],state,-a/6,b);end
    coefficients[4]
end
function _word_value(data,word,b)
    state=data.hf;amplitude=big(1)//big(1)
    for label in reverse(word)
        amplitude=_times(amplitude,_value(data.units[label],state,b),b);state⊻=data.units[label].flip
    end
    state,amplitude
end
function _literal_atom(data,order,b)
    a=order[1];c=order[2]
    if length(order)==2
        # Fully expanded polynomial rather than nested-commutator arithmetic.
        terms=[((a,a,c),-1//24),((a,c,a),1//12),((c,a,a),-1//24),
               ((c,a,c),-1//6),((c,c,a),1//12),((a,c,c),1//12)]
    else
        d=order[3]
        terms=[((c,a,d),-1//6),((c,d,a),1//12),((a,d,c),1//12),
               ((d,a,c),-1//6),((d,c,a),1//12),((a,c,d),1//12)]
    end
    out=Dict{UInt64,Rat}()
    for (word,scale) in terms
        state,v=_word_value(data,word,b);q=BigInt(numerator(scale))//BigInt(denominator(scale))
        _put!(out,state,_times(q,v,b),b)
    end
    out
end
function _variant_orders(scope)
    length(scope)==2 && return [[scope[1],scope[2]],[scope[2],scope[1]]]
    a,c,d=scope;[[a,c,d],[a,d,c],[c,a,d],[c,d,a],[d,a,c],[d,c,a]]
end
"""Complete reverse-endpoint XOR reconstruction, never producer scope lists."""
function checker_projection(data,targets,b)
    dest=sort(UInt64[targets...]);length(unique(dest))==length(dest) && 1<=length(dest)<=8 || error("checker target count")
    all(d->d!=data.hf && d<UInt64(2)^data.n,dest) || error("checker target outside register/HF")
    m=length(data.units);flips=Set(d⊻data.hf for d in dest)
    buckets=Dict{UInt64,Vector{Int}}()
    for label in 1:m;push!(get!(buckets,data.units[label].flip,Int[]),label);end
    scopes=Set{Tuple{Vararg{Int}}}();probes=0
    function keep!(scope)
        tick!(b)
        if !(scope in scopes)
            length(scopes)<b.caps.max_scopes || throw(Refusal("scope_limit","independent complete coverage cap"))
            reserve!(b,0;bytes=128+16length(scope));push!(scopes,scope)
        end
    end
    for right in 2:m,left in 1:right-1
        tick!(b)
        (data.units[left].flip in flips || data.units[right].flip in flips) && keep!((left,right))
    end
    for right in 3:m,middle in 2:right-1,g in flips
        probes+=1;probes<=b.caps.max_probes || throw(Refusal("probe_limit","independent target lookups capped"));tick!(b)
        need=g⊻data.units[middle].flip⊻data.units[right].flip
        for left in get(buckets,need,Int[]);left<middle && keep!((left,middle,right));end
    end
    sum(factorial(length(s)) for s in scopes;init=0)<=b.caps.max_variants || throw(Refusal("variant_limit","independent variants capped"))
    reserve!(b,0;bytes=256sum(factorial(length(s)) for s in scopes;init=0)+256length(scopes))
    3bitlength(data.DA)+bitlength(data.DK)+6<=b.caps.max_bits || throw(Refusal("bit_limit","independent prospective response scale"))
    E=checked(big(48)*data.DA^3*data.DK,b)
    2bitlength(E)<=b.caps.max_bits || throw(Refusal("bit_limit","independent prospective squared scale"))
    checked(E^2,b)
    dm=Dict(d=>i for (i,d) in enumerate(dest));factors=VectorFactor[]
    for scope in sort!(collect(scopes))
        orders=_variant_orders(collect(scope));vectors=Dict{Int,BigInt}[]
        for order in orders
            column=_literal_atom(data,order,b);v=Dict{Int,BigInt}()
            for (d,a) in column
                haskey(dm,d) || continue
                for component in 1:2
                    q=_times(_times(a,data.kernel[Int(d)+1][component],b),E//big(1),b)
                    denominator(q)==1 || error("independent48DA3DK does not clear atom denominator")
                    if !iszero(q);reserve!(b,1);v[2dm[d]-2+component]=checked(numerator(q),b);end
                end
            end
            push!(vectors,v)
        end
        push!(factors,VectorFactor(collect(scope),orders,vectors))
    end
    bind=Dict{String,Any}("input_sha256"=>data.inputhash,"kernel_sha256"=>data.kernelhash,
        "n"=>data.n,"hf"=>string(data.hf),"time"=>data.time,"targets"=>string.(dest),
        "original_lengths"=>data.lengths,"target"=>"exact supplied coefficients and immutable supplied dyadic kernel",
        "coverage_scope_sha256"=>scope_digest(scopes),"coverage_scopes"=>length(scopes))
    IntegerResponse(collect(1:m),2length(dest),E,Dict{Int,BigInt}(),factors,"full_order_projection",bind)
end

function _allowed(r,comparisons,mask,b)
    length(comparisons)<=3 || error("checker comparisons capped")
    m=length(r.labels);mapping=Dict(x=>i for (i,x) in enumerate(r.labels))
    edges=Set{Tuple{Int,Int}}();seen=Set{Tuple{Int,Int}}()
    for (i,p) in enumerate(comparisons)
        length(p)==2 && haskey(mapping,Int(p[1])) && haskey(mapping,Int(p[2])) || error("checker comparison labels")
        a,c=Int(p[1]),Int(p[2]);a!=c || error("checker self comparison")
        minmax(a,c) in seen && error("checker duplicate comparison");push!(seen,minmax(a,c))
        push!(edges,iszero(mask&(1<<(i-1))) ? (a,c) : (c,a))
    end
    # Independent BFS from each comparison endpoint (not producer closure).
    reach=Set{Tuple{Int,Int}}()
    nodes=Set(x for e in edges for x in e)
    for source in nodes
        stack=[c for (a,c) in edges if a==source];visited=Set{Int}()
        while !isempty(stack)
            v=pop!(stack);v in visited && continue;push!(visited,v);push!(reach,(source,v));tick!(b)
            append!(stack,[c for (a,c) in edges if a==v])
        end
    end
    any(e->e[1]==e[2],reach) && return true,Vector{Vector{Int}}()
    choices=Vector{Vector{Int}}()
    for f in r.factors
        good=Int[]
        for (j,order) in enumerate(f.orders)
            tick!(b,length(order)^2);positions=Dict(x=>i for (i,x) in enumerate(order))
            ok=true
            for (a,c) in reach
                haskey(positions,a) && haskey(positions,c) && positions[a]>=positions[c] && (ok=false)
            end
            ok && push!(good,j)
        end
        isempty(good) && error("checker acyclic branch lost all variants");push!(choices,good)
    end
    false,choices
end
function independent_bound(r::IntegerResponse,cert;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);validate_response(r,b)
    cert["schema"]=="response-bound-v1" && cert["complete"]===true && cert["response_sha256"]==response_digest(r) && cert["domain"]==r.domain || error("independent target/status binding")
    comparisons=cert["comparisons"];length(comparisons)<=3 || error("checker branch cap")
    records=cert["branches"];length(records)==2^length(comparisons) || error("independent incomplete cover")
    sort([Int(rec["assignment"]) for rec in records])==collect(0:2^length(comparisons)-1) || error("independent duplicate/missing assignment")
    L=nothing;D=r.denominator
    for rec in records
        cyclic,choices=_allowed(r,comparisons,Int(rec["assignment"]),b)
        cyclic==rec["cyclic"] || error("independent false cycle claim")
        if cyclic
            isempty(rec["witnesses"]) && rationalparse(rec["L_Q"],caps)==0 || error("independent cyclic branch claim")
            continue
        end
        # Direct rational coordinates, separate from integer interval helper.
        low=Rat[get(r.baseline,j,big(0))//D for j in 1:r.dimension];high=copy(low)
        for (f,indices) in zip(r.factors,choices)
            coords=Set(j for i in indices for j in keys(f.variants[i]))
            for j in coords
                tick!(b,length(indices));vals=Rat[get(f.variants[i],j,big(0))//D for i in indices]
                low[j]=_plus(low[j],minimum(vals),b);high[j]=_plus(high[j],maximum(vals),b)
            end
        end
        value=ZERO
        for (l,u) in zip(low,high)
            nearest=l>0 ? l : (u<0 ? u : ZERO);value=_plus(value,_times(nearest,nearest,b),b)
        end
        length(rec["witnesses"])<=8 || error("independent witness count cap")
        for raw in rec["witnesses"]
            length(raw)==r.dimension || error("independent witness dimension")
            W=BigInt[bigparse(v,caps) for v in raw]
            all(v->abs(v)<=typemax(Int64),W) && any(v->!iszero(v),W) || error("independent invalid witness")
            norm=ZERO;lower=ZERO;upper=ZERO
            for j in 1:r.dimension
                w=W[j]//big(1);norm=_plus(norm,_times(w,w,b),b)
                term=_times(w,get(r.baseline,j,big(0))//D,b);lower=_plus(lower,term,b)
            end
            upper=lower
            for (f,indices) in zip(r.factors,choices)
                vals=Rat[]
                for i in indices
                    v=ZERO
                    for (j,a) in f.variants[i];v=_plus(v,_times(W[j]//big(1),a//D,b),b);end
                    push!(vals,v)
                end
                lower=_plus(lower,minimum(vals),b);upper=_plus(upper,maximum(vals),b)
            end
            distance=lower>0 ? lower : (upper<0 ? -upper : ZERO)
            value=max(value,_rcheck(_times(distance,distance,b)/norm,b))
        end
        rationalparse(rec["L_Q"],caps)==value || error("independent branch arithmetic differs")
        L=isnothing(L) ? value : min(L,value)
    end
    !isnothing(L) && L==rationalparse(cert["L_Q"],caps) || error("independent bound differs")
    (L_Q=L,work=b.work,estimated_peak_bytes=b.peak_bytes)
end
function independent_window_check(r::IntegerResponse,record,input;source_path=nothing,caps=Caps(),guard=()->nothing)
    record["schema"]=="window-bound-record-v1" && record["complete"]===true &&
        record["domain"]=="fixed_window" && r.domain=="fixed_window" && record["positive_full_bound"]===false || error("independent window scope differs")
    record["response_sha256"]==response_digest(r) || error("independent window target digest")
    binding=record["binding"]
    expectedbinding=Dict(k=>v for (k,v) in r.binding if k!="destinations")
    clean(binding)==clean(expectedbinding) || error("independent window bindings differ")
    filesize(input)<=caps.max_input_bytes || throw(Refusal("input_limit","window source too large"))
    bytes=read(input);data=TOML.parse(String(copy(bytes)))
    bytes2hex(sha256(bytes))==binding["input_sha256"] || error("independent window source changed")
    !isnothing(source_path) && filehash(source_path)!=binding["response_source_sha256"] && error("independent frozen-response source changed")
    lengths=[length(u["words"]) for u in data["units"]]
    all(>(0),lengths) && lengths==Int.(binding["original_lengths"]) || error("independent original window word counts")
    parent=Int.(binding["parent_order"]);positions=Int.(binding["positions"])
    sort(parent)==collect(1:length(lengths)) && !isempty(positions) || error("independent parent order")
    positions==collect(first(positions):last(positions)) && first(positions)>1 && last(positions)<length(parent) &&
        parent[positions]==r.labels || error("independent interior-window mapping")
    B=bigparse(binding["ideal_budget"],caps);0<=B && bitlength(B)<=128 || error("independent window budget")
    U=sum(BigInt.(lengths));a=lengths[first(parent)];c=lengths[last(parent)]
    steps=max(big(0),fld(B-a,2U-a-c));steps>0 || error("independent window affordability")
    b=Budget(caps;guard);positions_by_label=Dict(x=>i for (i,x) in enumerate(r.labels))
    vector=Rat[get(r.baseline,j,big(0))//r.denominator for j in 1:r.dimension]
    for f in r.factors
        indices=[i for (i,order) in enumerate(f.orders) if all(positions_by_label[order[j]]<positions_by_label[order[j+1]] for j in 1:length(order)-1)]
        length(indices)==1 || error("independent window induced pattern")
        for (j,v) in f.variants[only(indices)];vector[j]=_plus(vector[j],v//r.denominator,b);end
    end
    Q=ZERO
    for v in vector;Q=_plus(Q,_times(v,v,b),b);end
    proof=independent_bound(r,record["certificate"];caps,guard);L=proof.L_Q
    Q==rationalparse(record["Q_inc"],caps) && L==rationalparse(record["L_Q"],caps) &&
        steps==bigparse(record["r_inc"],caps) && Q/steps^4==rationalparse(record["UB_J"],caps) &&
        L/steps^4==rationalparse(record["L_J"],caps) || error("independent fixed-window Q/J bridge")
    if L>0
        Q/L==rationalparse(record["ratio"],caps) || error("independent window ratio")
    else
        record["ratio"] in (nothing,"unavailable") || error("independent zero-window ratio")
    end
    (status="independently_verified_window",L_Q=L,UB_Q=Q,L_J=L/steps^4,UB_J=Q/steps^4,
        full_scope_verified=false,work=b.work+proof.work)
end
function independent_full_check(input,kernelpath,record;caps=Caps(),guard=()->nothing)
    record["schema"]=="molecular-bound-record-v1" || error("checker record schema")
    record["complete"]===true || error("refused/incomplete run is not a positive certificate")
    record["domain"]=="all complete orders; original unit partition" || error("checker full-domain tag differs")
    b=Budget(caps;guard);data=_checker_source(input,kernelpath,b)
    record["input_sha256"]==data.inputhash && record["kernel_sha256"]==data.kernelhash || error("checker input/kernel changed")
    targets=UInt64[parse(UInt64,string(x)) for x in record["targets"]]
    projection=checker_projection(data,targets,b)
    response_digest(projection)==record["projection_sha256"] || error("checker projected response differs")
    projection.binding["coverage_scope_sha256"]==record["coverage_scope_sha256"] || error("checker scope coverage differs")
    projection.binding["coverage_scopes"]==record["coverage_scopes"] || error("checker scope count differs")
    bound=independent_bound(projection,record["certificate"];caps,guard)
    order=Int.(record["full_order"]);column=_polynomial_defect(data,order,b)
    Q=ZERO;E=projection.denominator;N=big(0);energies=Dict{UInt64,Rat}()
    for (state,a) in column
        state==data.hf && continue
        energy=ZERO
        for component in 1:2
            value=_times(a,data.kernel[Int(state)+1][component],b)
            Q=_plus(Q,_times(value,value,b),b)
            energy=_plus(energy,_times(value,value,b),b)
            encoded=_times(value,E//big(1),b);denominator(encoded)==1 || error("checker48DA3DK Taylor normalization")
            N=addchecked(N,mulchecked(numerator(encoded),numerator(encoded),b),b)
        end
        !iszero(energy) && (energies[state]=energy)
    end
    count=Int(record["requested_target_count"]);1<=count<=8 || error("checker requested target count")
    expectedtargets=sort!(collect(keys(energies));by=d->(-energies[d],d))
    expectedtargets=isempty(expectedtargets) ? UInt64[data.hf==0 ? 1 : 0] : sort!(expectedtargets[1:min(count,length(expectedtargets))])
    sort(targets)==expectedtargets || error("checker deterministic target selection differs")
    Q==N//checked(E^2,b) || error("checker full norm scale")
    Q==rationalparse(record["Q_inc"],caps) && bigparse(record["UB_numerator"],caps)==N &&
        bigparse(record["UB_denominator"],caps)==E^2 && bigparse(record["response_denominator"],caps)==E || error("checker same-model UB differs")
    # Independently compute affordability, with the actual endpoints included.
    B=bigparse(record["ideal_budget"],caps);B>=0 || error("checker negative budget")
    bitlength(B)<=128 || throw(Refusal("budget_limit","independent ideal budget cap"))
    U=sum(BigInt.(data.lengths));m=length(data.lengths)
    r(a,c)=m==1 ? (B>=U ? big(1) : big(0)) : max(big(0),fld(B-data.lengths[a],2U-data.lengths[a]-data.lengths[c]))
    rmax=m==1 ? r(1,1) : maximum(r(a,c) for a in 1:m for c in 1:m if a!=c)
    rinc=r(first(order),last(order));rinc>0 || error("checker unaffordable incumbent")
    LJ=bound.L_Q/rmax^4;UBJ=Q/rinc^4
    LJ<=UBJ || error("checker lower exceeds feasible UB")
    ratio=LJ>0 ? UBJ/LJ : nothing
    rationalparse(record["L_Q"],caps)==bound.L_Q && rationalparse(record["L_J"],caps)==LJ && rationalparse(record["UB_J"],caps)==UBJ || error("checker J/response bridge differs")
    bigparse(record["r_inc"],caps)==rinc && bigparse(record["r_max"],caps)==rmax || error("checker endpoint normalization differs")
    if isnothing(ratio)
        record["ratio"] in (nothing,"unavailable") || error("checker zero floor has ratio")
    else
        rationalparse(record["ratio"],caps)==ratio || error("checker ratio differs")
    end
    record["positive_full_bound"]==(LJ>0) && record["desired_ratio_met"]==(LJ>0 && UBJ<=2LJ) || error("checker usefulness status differs")
    record["status"]==(LJ>0 ? "positive_full_bound" : "zero_full_bound") || error("checker outcome status differs")
    (status="independently_verified",L_Q=bound.L_Q,L_J=LJ,UB_J=UBJ,ratio=ratio,
     full_scope_verified=true,work=b.work+bound.work,estimated_peak_bytes=max(b.peak_bytes,bound.estimated_peak_bytes))
end
end
