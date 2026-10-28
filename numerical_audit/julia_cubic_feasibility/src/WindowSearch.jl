module WindowSearch
using ..ResponseCore: WindowModel
export enumerate_window, insertion_window, prepare_window, exact_numerator

struct FlatFactor
    scope::Vector{Int}
    destinations::Vector{Int}
    values::Matrix{ComplexF64}
end
struct FlatWindow
    labels::Vector{Int}
    destinations::Vector{UInt64}
    baseline::Vector{ComplexF64}
    factors::Vector{FlatFactor}
    error_bound::Float64
end

function pattern(scope,positions)
    a,b=scope[1:2]
    code=Int(positions[a]>positions[b])
    length(scope)==3 && (code+=2Int(positions[a]>positions[scope[3]])+4Int(positions[b]>positions[scope[3]]))
    code+1
end

function _validate_flat_for_screen(flat::FlatWindow)
    k=length(flat.labels);d=length(flat.destinations)
    1<=k<=9 || throw(ArgumentError("window width outside 1:9 cap"))
    length(flat.baseline)==d || throw(ArgumentError("baseline response dimension differs"))
    length(unique(flat.labels))==k || throw(ArgumentError("duplicate window labels"))
    length(unique(flat.destinations))==d || throw(ArgumentError("duplicate response destinations"))
    isfinite(flat.error_bound) && flat.error_bound>=0 || throw(ArgumentError("invalid roundoff-bound hint"))
    for factor in flat.factors
        arity=length(factor.scope)
        2<=arity<=3 || throw(ArgumentError("response factor arity must be two or three"))
        length(unique(factor.scope))==arity || throw(ArgumentError("duplicate response factor labels"))
        all(i->1<=i<=k,factor.scope) || throw(ArgumentError("response factor uses nonlocal labels"))
        length(unique(factor.destinations))==length(factor.destinations) || throw(ArgumentError("duplicate factor destinations"))
        all(i->1<=i<=d,factor.destinations) || throw(ArgumentError("factor destination outside baseline support"))
        size(factor.values)==(length(factor.destinations),arity==2 ? 2 : 8) ||
            throw(ArgumentError("response factor matrix shape differs"))
    end
    allvalues=Iterators.flatten((flat.baseline,(v for f in flat.factors for v in f.values)))
    all(v->isfinite(v) && (iszero(v)||1e-140<=abs(v)<=1e140),allvalues) ||
        throw(ArgumentError("frozen table magnitude outside certificate arithmetic range"))
    nothing
end

# Recompute from the current supplied arrays; a public FlatWindow may carry a
# manually supplied or stale hint, which cannot establish a certificate itself.
function _roundoff_bound(flat::FlatWindow)
    _validate_flat_for_screen(flat)
    baseline=flat.baseline;factors=flat.factors;destinations=flat.destinations
    envelope=setprecision(BigFloat,128) do
        setrounding(BigFloat,RoundUp) do
            re=BigFloat[abs(real(v)) for v in baseline]
            im=BigFloat[abs(imag(v)) for v in baseline]
            for f in factors, (i,at) in enumerate(f.destinations)
                re[at]+=maximum(abs(BigFloat(real(v))) for v in f.values[i,:])
                im[at]+=maximum(abs(BigFloat(imag(v))) for v in f.values[i,:])
            end
            h2=sum(x*x for x in re;init=BigFloat(0))+sum(x*x for x in im;init=BigFloat(0))
            a=BigFloat(length(factors))*BigFloat(eps(Float64))
            b=BigFloat(4length(destinations)+10)*BigFloat(eps(Float64))
            a<0.25 && b<0.25 || error("roundoff envelope count cap")
            # Include an absolute allowance for subnormal rounding as well.
            tiny=BigFloat(nextfloat(0.0))
            underflow=(length(factors)+4length(destinations)+10)*tiny*(1+sqrt(h2))
            ((2a+a*a)+b*(1+a)^2)*h2+underflow
        end
    end
    bound=nextfloat(Float64(envelope))
    isfinite(bound) || error("nonfinite roundoff envelope")
    bound
end

function prepare_window(model::WindowModel)
    labels=copy(model.baseline_window)
    2<=length(labels)<=9 || throw(ArgumentError("window width outside 2:9 cap"))
    length(unique(labels))==length(labels) || throw(ArgumentError("duplicate window labels"))
    labelmap=Dict(v=>i for (i,v) in enumerate(labels))
    for factor in model.factors
        arity=length(factor.scope)
        2<=arity<=3 || throw(ArgumentError("response factor arity must be two or three"))
        length(unique(factor.scope))==arity || throw(ArgumentError("duplicate response factor labels"))
        all(v->haskey(labelmap,v),factor.scope) || throw(ArgumentError("factor scope outside window"))
        length(factor.variants)==factorial(arity) || throw(ArgumentError("incomplete response order variants"))
        for key in keys(factor.variants)
            length(key)==arity && Set(key)==Set(factor.scope) || throw(ArgumentError("invalid response order variant"))
        end
    end
    destinations=sort!(collect(union(Set(keys(model.baseline_response)),
        Set(d for f in model.factors for v in values(f.variants) for d in keys(v)))))
    dmap=Dict(v=>i for (i,v) in enumerate(destinations))
    baseline=ComplexF64[get(model.baseline_response,d,0.0im) for d in destinations]
    factors=FlatFactor[]
    for factor in model.factors
        scope=Int[labelmap[v] for v in factor.scope]
        dest=sort!(collect(Set(dmap[d] for v in values(factor.variants) for d in keys(v))))
        table=zeros(ComplexF64,length(dest),length(scope)==2 ? 2 : 8)
        valid=length(scope)==2 ? (1,2) : (1,2,4,5,7,8)
        seen=Set{Int}()
        for (key,vector) in factor.variants
            pos=zeros(Int,length(labels))
            for (i,label) in enumerate(key);pos[labelmap[label]]=i;end
            column=pattern(scope,pos)
            column in valid && !(column in seen) || throw(ArgumentError("invalid or duplicate response comparison column"))
            push!(seen,column)
            for (i,at) in enumerate(dest);table[i,column]=get(vector,destinations[at],0.0im);end
        end
        length(seen)==length(valid) || throw(ArgumentError("incomplete response comparison columns"))
        push!(factors,FlatFactor(scope,dest,table))
    end
    flat=FlatWindow(labels,destinations,baseline,factors,0.0)
    FlatWindow(labels,destinations,baseline,factors,_roundoff_bound(flat))
end

function response!(buffer,flat::FlatWindow,order,positions)
    for (i,label) in enumerate(order); positions[label]=i; end
    copyto!(buffer,flat.baseline)
    for f in flat.factors
        column=pattern(f.scope,positions)
        @inbounds for (i,at) in enumerate(f.destinations)
            buffer[at]+=f.values[i,column]
        end
    end
    score=0.0
    @inbounds for v in buffer; score+=real(v)*real(v)+imag(v)*imag(v); end
    isfinite(score) && score>=0 || error("invalid frozen table score")
    score
end

function exact_numerator(flat::FlatWindow,order)
    positions=zeros(Int,length(flat.labels))
    for (i,label) in enumerate(order); positions[label]=i; end
    re=Rational{BigInt}[Rational{BigInt}(real(v)) for v in flat.baseline]
    im=Rational{BigInt}[Rational{BigInt}(imag(v)) for v in flat.baseline]
    for f in flat.factors
        column=pattern(f.scope,positions)
        for (i,at) in enumerate(f.destinations)
            re[at]+=Rational{BigInt}(real(f.values[i,column]))
            im[at]+=Rational{BigInt}(imag(f.values[i,column]))
        end
    end
    sum((x*x for x in re);init=big(0)//big(1))+sum((x*x for x in im);init=big(0)//big(1))
end

"""Heap enumeration visits k! permutations; stopping never becomes a certificate."""
function visit_orders(visit,order)
    function rec(k)
        k==1 && return visit(order)
        rec(k-1)
        for i in 1:k-1
            at=iseven(k) ? i : 1
            order[at],order[k]=order[k],order[at]
            rec(k-1)
        end
    end
    rec(length(order))
end

function enumerate_window(flat::FlatWindow; guard=()->nothing,max_exact_candidates=50_000)
    bound=max(flat.error_bound,_roundoff_bound(flat))
    k=length(flat.labels); total=Int(factorial(big(k)))
    total<=362880 || error("enumeration count cap")
    buffer=copy(flat.baseline); positions=zeros(Int,k)
    permutations=Matrix{Int16}(undef,k,total)
    scores=Vector{Float64}(undef,total)
    count=0; best=Inf; bestcol=1
    visit_orders(collect(1:k)) do order
        count+=1
        count%128==1 && guard()
        score=response!(buffer,flat,order,positions)
        scores[count]=score
        permutations[:,count]=order
        if score<best; best=score; bestcol=count; end
    end
    count==total || error("incomplete permutation coverage")
    guard()
    # Candidate intervals encompass exact evaluation of supplied dyadic tables.
    # Choose a preliminary exact incumbent; screen with outward rounded bounds.
    incumbent=exact_numerator(flat,Int.(permutations[:,bestcol]))
    incumbent_upper=nextfloat(Float64(incumbent))
    contenders=Int[i for i in 1:total if prevfloat(scores[i]-bound)<=incumbent_upper]
    length(contenders)<=max_exact_candidates || error("exact contender cap; no certificate")
    bestcol in contenders || error("roundoff envelope excluded its exact incumbent; no certificate")
    minimizers=0
    bestorder=Int.(permutations[:,bestcol])
    for (at,i) in enumerate(contenders)
        at%32==1 && guard()
        order=Int.(permutations[:,i]); value=exact_numerator(flat,order)
        if value<incumbent
            incumbent=value; bestorder=order; minimizers=1
        elseif value==incumbent
            minimizers+=1
            isless(Tuple(flat.labels[order]),Tuple(flat.labels[bestorder])) && (bestorder=order)
        end
    end
    minimizers>=1 || error("no exact minimizer retained; no certificate")
    # Every screened-out candidate has lower bound above the initial incumbent,
    # which is >= the final exact incumbent. Exhaustion plus exact ties certify
    # the restricted optimum of this frozen table model only.
    (status="complete",certificate="exact optimum of supplied frozen binary64 response DELTA tables",
     order=flat.labels[bestorder],numerator=incumbent,score=Float64(incumbent),
     permutations=count,exact_contenders=length(contenders),minimizers=minimizers,
     roundoff_bound=bound,destinations=length(flat.destinations),factors=length(flat.factors))
end

function insertion_window(flat::FlatWindow; guard=()->nothing,max_sweeps=10)
    bound=max(flat.error_bound,_roundoff_bound(flat))
    k=length(flat.labels); current=collect(1:k)
    positions=zeros(Int,k); buffer=copy(flat.baseline)
    score=response!(buffer,flat,current,positions); evaluations=1; moves=0
    for sweep in 1:max_sweeps
        best=copy(current); value=score
        for from in 1:k, to in 1:k
            from==to && continue
            guard()
            trial=copy(current); insert!(trial,to,popat!(trial,from))
            next=response!(buffer,flat,trial,positions); evaluations+=1
            if next<value; value=next; best=trial; end
        end
        value<score-bound*2 || break
        current=best; score=value; moves+=1
    end
    exact=exact_numerator(flat,current)
    (order=flat.labels[current],score=Float64(exact),numerator=exact,evaluations=evaluations,moves=moves)
end
end
