module ResponseBounds
using ..BoundsCommon
export encode_window, bound_response, verify_bound, branch_closure, coordinate_floor,
       separator, allowed_variants, propose_directions, window_bound

function _sparse_float(values,indices,D,b)
    out=Dict{Int,BigInt}()
    for (i,x) in zip(indices,values)
        tick!(b);isfinite(x) || throw(ArgumentError("nonfinite window component"))
        iszero(x) && continue
        q=Rational{BigInt}(x);d=denominator(q)
        scale,rem=divrem(D,d);iszero(rem) || error("nondividing dyadic denominator")
        reserve!(b,1);out[i]=mulchecked(numerator(q),scale,b)
    end
    out
end
"""Encode supplied frozen Float vectors, without recomputing any transport.

FlatWindow is intentionally duck typed: no historical module is modified or
loaded by this backend. Binding must specify full exterior/endpoints for J.
"""
function encode_window(flat;binding=Dict{String,Any}(),caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);k=length(flat.labels);d=length(flat.destinations)
    1<=k<=9 || throw(ArgumentError("frozen window width must be 1:9"))
    length(unique(flat.labels))==k && length(unique(flat.destinations))==d || throw(ArgumentError("duplicate labels/destinations"))
    length(flat.baseline)==d || throw(ArgumentError("baseline dimension differs"))
    length(flat.factors)<=binomial(k,2)+binomial(k,3) || throw(Refusal("scope_limit","too many window factors"))
    d<=caps.max_states || throw(Refusal("dimension_limit","window support too large"))
    D=big(1)
    function include!(x)
        tick!(b);isfinite(x) || throw(ArgumentError("nonfinite supplied component"))
        q=Rational{BigInt}(x);D=max(D,denominator(q));checked(D,b)
    end
    for z in flat.baseline;include!(real(z));include!(imag(z));end
    seen=Set{Tuple{Vararg{Int}}}()
    for f in flat.factors
        tick!(b,1+length(f.scope)+length(f.destinations))
        a=length(f.scope);a in (2,3) || throw(ArgumentError("factor arity differs"))
        length(unique(f.scope))==a && all(x->1<=x<=k,f.scope) || throw(ArgumentError("invalid local scope"))
        key=Tuple(sort(f.scope));key in seen && throw(ArgumentError("duplicate local scope"));push!(seen,key)
        length(unique(f.destinations))==length(f.destinations) && all(x->1<=x<=d,f.destinations) || throw(ArgumentError("invalid factor destinations"))
        size(f.values)==(length(f.destinations),a==2 ? 2 : 8) || throw(ArgumentError("invalid table shape"))
        for col in (a==2 ? (1,2) : (1,2,4,5,7,8)),row in axes(f.values,1)
            z=f.values[row,col];include!(real(z));include!(imag(z))
        end
    end
    mulchecked(D,D,b)
    base=_sparse_float(Iterators.flatten(((real(z),imag(z)) for z in flat.baseline)),1:2d,D,b)
    factors=VectorFactor[]
    for f in flat.factors
        scope=sort!(Int[flat.labels[i] for i in f.scope]);orders=permutations(scope)
        variants=Dict{Int,BigInt}[]
        for order in orders
            pos=Dict(x=>i for (i,x) in enumerate(order))
            ids=flat.labels[f.scope];code=Int(pos[ids[1]]>pos[ids[2]])
            length(ids)==3 && (code+=2Int(pos[ids[1]]>pos[ids[3]])+4Int(pos[ids[2]]>pos[ids[3]]))
            vals=Iterators.flatten(((real(z),imag(z)) for z in f.values[:,code+1]))
            coords=Iterators.flatten(((2i-1,2i) for i in f.destinations))
            push!(variants,_sparse_float(vals,coords,D,b))
        end
        push!(factors,VectorFactor(scope,orders,variants))
    end
    sort!(factors;by=f->Tuple(f.scope))
    bind=Dict{String,Any}(string(k)=>v for (k,v) in binding)
    bind["destinations"]=string.(flat.destinations)
    r=IntegerResponse(Int.(flat.labels),2d,D,base,factors,"fixed_window",bind)
    validate_response(r,Budget(caps;guard));r
end

"""Exact transitive closure of at most three selected comparisons.
No cubic-in-molecule Floyd matrix is needed: only comparison endpoints occur.
"""
function branch_closure(labels,comparisons,assignment)
    length(comparisons)<=3 || throw(ArgumentError("at most three comparisons"))
    0<=assignment<(1<<length(comparisons)) || throw(ArgumentError("branch mask outside cover"))
    nodes=Set{Int}();edges=Set{Tuple{Int,Int}}();pairs=Set{Tuple{Int,Int}}()
    for (i,pair) in enumerate(comparisons)
        length(pair)==2 || throw(ArgumentError("comparison must have two labels"))
        a,c=Int(pair[1]),Int(pair[2]);a!=c && a in labels && c in labels || throw(ArgumentError("invalid branch comparison"))
        key=minmax(a,c);key in pairs && throw(ArgumentError("duplicate comparison"));push!(pairs,key)
        iszero(assignment & (1<<(i-1))) ? push!(edges,(a,c)) : push!(edges,(c,a))
        push!(nodes,a);push!(nodes,c)
    end
    reach=copy(edges)
    for middle in nodes,a in nodes,c in nodes
        (a,middle) in reach && (middle,c) in reach && push!(reach,(a,c))
    end
    any(a->(a,a) in reach,nodes),reach
end
function allowed_variants(r,closure,b)
    result=Vector{Vector{Int}}()
    for f in r.factors
        permitted=Int[]
        for (i,p) in enumerate(f.orders)
            tick!(b,length(p)^2)
            pos=Dict(a=>j for (j,a) in enumerate(p))
            all(e->!(e[1] in f.scope && e[2] in f.scope) || pos[e[1]]<pos[e[2]],closure) && push!(permitted,i)
        end
        isempty(permitted) && error("acyclic cover has an empty local variant set")
        push!(result,permitted)
    end
    result
end
function coordinate_floor(r,allowed,b)
    lo=copy(r.baseline);hi=copy(r.baseline)
    for (f,choices) in zip(r.factors,allowed)
        coords=Set(j for c in choices for j in keys(f.variants[c]))
        for j in coords
            tick!(b,length(choices))
            vals=[get(f.variants[c],j,big(0)) for c in choices]
            lo[j]=addchecked(get(lo,j,big(0)),minimum(vals),b)
            hi[j]=addchecked(get(hi,j,big(0)),maximum(vals),b)
        end
    end
    total=big(0)
    for j in union(keys(lo),keys(hi))
        l=get(lo,j,big(0));u=get(hi,j,big(0));g=l>0 ? l : (u<0 ? -u : big(0))
        total=addchecked(total,mulchecked(g,g,b),b)
    end
    total//mulchecked(r.denominator,r.denominator,b)
end
function _dot(W,v,b)
    x=big(0)
    for (j,a) in v;x=addchecked(x,mulchecked(W[j],a,b),b);end
    x
end
function separator(r,allowed,W,b)
    length(W)==r.dimension || throw(ArgumentError("witness dimension differs"))
    all(x->x isa Integer,W) || throw(ArgumentError("witness must be integer"))
    w=BigInt[checked(BigInt(x),b) for x in W]
    z=big(0)
    for x in w;z=addchecked(z,mulchecked(x,x,b),b);end
    z>0 || throw(ArgumentError("zero separator witness"))
    lo=_dot(w,r.baseline,b);hi=lo
    for (f,choices) in zip(r.factors,allowed)
        vals=[_dot(w,f.variants[c],b) for c in choices]
        lo=addchecked(lo,minimum(vals),b);hi=addchecked(hi,maximum(vals),b)
    end
    g=lo>0 ? lo : (hi<0 ? -hi : big(0))
    mulchecked(g,g,b)//mulchecked(mulchecked(r.denominator,r.denominator,b),z,b)
end

"""Floating proposal only. Every returned direction is independently exact-checked."""
function propose_directions(r,allowed,b;iterations=64,max_witnesses=8)
    0<=iterations<=64 && 0<=max_witnesses<=8 || throw(ArgumentError("proposal schedule exceeds frozen caps"))
    max_witnesses==0 && return Vector{BigInt}[]
    d=r.dimension;d==0 && return Vector{BigInt}[]
    base=zeros(Float64,d)
    for (j,v) in r.baseline;tick!(b);base[j]=Float64(v//r.denominator);end
    data=Vector{Vector{Dict{Int,Float64}}}()
    largest=maximum(abs,base;init=0.0)
    for (f,choices) in zip(r.factors,allowed)
        row=Dict{Int,Float64}[]
        for c in choices
            tick!(b,1+length(f.variants[c]))
            vec=Dict(j=>Float64(v//r.denominator) for (j,v) in f.variants[c])
            all(isfinite,values(vec)) || return Vector{BigInt}[]
            for v in values(vec);largest=max(largest,abs(v));end
            push!(row,vec)
        end
        push!(data,row)
    end
    isfinite(largest) && largest>0 || return Vector{BigInt}[]
    base./=largest
    for row in data,v in row,j in keys(v);v[j]/=largest;end
    x=copy(base)
    for row in data,(j,v) in first(row);x[j]+=v;end
    directions=Vector{BigInt}[]
    function retain!()
        length(directions)>=max_witnesses && return
        scale=maximum(abs,x;init=0.0)
        isfinite(scale) && scale>0 || return
        W=BigInt[round(Int64,(v/scale)*(2.0^30-1)) for v in x]
        any(v->!iszero(v),W) || return
        # Avoid constructing a thousands-coordinate Tuple/type for HF windows.
        tick!(b,length(directions)*d)
        any(old->old==W,directions) && return
        push!(directions,W)
    end
    retain!()
    for iteration in 1:iterations
        s=copy(base)
        for row in data
            tick!(b,1+sum(length(v) for v in row))
            costs=[sum(x[j]*a for (j,a) in v;init=0.0) for v in row]
            best=argmin(costs)
            for (j,a) in row[best];s[j]+=a;end
        end
        h=s-x;hh=sum(abs2,h);isfinite(hh) && hh>0 || break
        gamma=clamp(-sum(x.*h)/hh,0.0,1.0);x.+=gamma.*h
        iteration%8==0 && retain!()
    end
    retain!();directions
end
function bound_response(r::IntegerResponse;comparisons=Tuple{Int,Int}[],iterations=64,
        max_witnesses=8,caps=Caps(),guard=()->nothing)
    length(comparisons)<=3 || throw(ArgumentError("cover exceeds three-comparison cap"))
    b=Budget(caps;guard);validate_response(r,b)
    branches=Dict{String,Any}[];lower=nothing
    for assignment in 0:(1<<length(comparisons))-1
        tick!(b);cyclic,closure=branch_closure(r.labels,comparisons,assignment)
        if cyclic
            push!(branches,Dict("assignment"=>assignment,"cyclic"=>true,"witnesses"=>Vector{BigInt}[],"L_Q"=>big(0)//big(1)))
            continue
        end
        allowed=allowed_variants(r,closure,b);value=coordinate_floor(r,allowed,b)
        ws=propose_directions(r,allowed,b;iterations,max_witnesses)
        for W in ws;value=max(value,separator(r,allowed,W,b));end
        lower=isnothing(lower) ? value : min(lower,value)
        push!(branches,Dict("assignment"=>assignment,"cyclic"=>false,"witnesses"=>ws,"L_Q"=>value))
    end
    isnothing(lower) && error("comparison cover has no acyclic branch")
    Dict{String,Any}("schema"=>"response-bound-v1","response_sha256"=>response_digest(r),
        "domain"=>r.domain,"comparisons"=>[collect(p) for p in comparisons],"branches"=>branches,
        "L_Q"=>lower,"complete"=>true,"proposal_iterations"=>iterations,
        "work"=>b.work,"estimated_peak_bytes"=>b.peak_bytes)
end
function verify_bound(r::IntegerResponse,cert;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);validate_response(r,b)
    cert["schema"]=="response-bound-v1" && cert["response_sha256"]==response_digest(r) && cert["domain"]==r.domain ||
        throw(ArgumentError("certificate target binding differs"))
    cert["complete"]===true || throw(ArgumentError("incomplete certificate"))
    comparisons=cert["comparisons"];length(comparisons)<=3 || throw(ArgumentError("cover exceeds cap"))
    records=cert["branches"]
    length(records)==1<<length(comparisons) || throw(ArgumentError("incomplete branch cover"))
    Set(Int(x["assignment"]) for x in records)==Set(0:(1<<length(comparisons))-1) || throw(ArgumentError("branch cover differs"))
    result=nothing
    for rec in records
        cyclic,closure=branch_closure(r.labels,comparisons,Int(rec["assignment"]))
        cyclic==rec["cyclic"] || throw(ArgumentError("false cyclic-branch claim"))
        length(rec["witnesses"])<=8 || throw(ArgumentError("too many witnesses"))
        if cyclic
            isempty(rec["witnesses"]) && rationalparse(rec["L_Q"],caps)==0 || throw(ArgumentError("cyclic branch carries a claimed value"))
            continue
        end
        allowed=allowed_variants(r,closure,b);value=coordinate_floor(r,allowed,b)
        for raw in rec["witnesses"]
            length(raw)==r.dimension || throw(ArgumentError("witness dimension differs"))
            W=BigInt[bigparse(v,caps) for v in raw]
            all(v->abs(v)<=typemax(Int64),W) || throw(ArgumentError("witness component exceeds Int64 proposal contract"))
            value=max(value,separator(r,allowed,W,b))
        end
        value==rationalparse(rec["L_Q"],caps) || throw(ArgumentError("claimed branch lower bound differs"))
        result=isnothing(result) ? value : min(result,value)
    end
    !isnothing(result) && result==rationalparse(cert["L_Q"],caps) || throw(ArgumentError("claimed lower bound differs"))
    result
end
function window_bound(flat;binding=Dict{String,Any}(),caps=Caps(),guard=()->nothing,kwargs...)
    r=encode_window(flat;binding,caps,guard);cert=bound_response(r;caps,guard,kwargs...)
    (response=r,certificate=cert,lower=verify_bound(r,cert;caps,guard))
end
end
