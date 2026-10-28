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

function prepare_window(model::WindowModel)
    labels=copy(model.baseline_window)
    2 <= length(labels) <= 9 || error("window width outside 2:9 cap")
    labelmap=Dict(v=>i for (i,v) in enumerate(labels))
    destinations=sort!(collect(union(Set(keys(model.baseline_response)),
        Set(d for f in model.factors for v in values(f.variants) for d in keys(v)))))
    dmap=Dict(v=>i for (i,v) in enumerate(destinations))
    baseline=ComplexF64[get(model.baseline_response,d,0.0im) for d in destinations]
    factors=FlatFactor[]
    for factor in model.factors
        scope=Int[labelmap[v] for v in factor.scope]
        dest=sort!(collect(Set(dmap[d] for v in values(factor.variants) for d in keys(v))))
        table=zeros(ComplexF64,length(dest),length(scope)==2 ? 2 : 8)
        for (key,vector) in factor.variants
            pos=zeros(Int,length(labels))
            for (i,label) in enumerate(key); pos[labelmap[label]]=i; end
            column=pattern(scope,pos)
            for (i,at) in enumerate(dest); table[i,column]=get(vector,destinations[at],0.0im); end
        end
        push!(factors,FlatFactor(scope,dest,table))
    end
    # Error relative to the EXACT dyadic values of these frozen DELTA tables.
    # No fastmath: each coordinate is formed by <=F additions. A generous
    # gamma bound also covers products/additions in the squared norm.
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
    allvalues=Iterators.flatten((baseline,(v for f in factors for v in f.values)))
    all(v->isfinite(v) && (iszero(v)||1e-140<=abs(v)<=1e140),allvalues) ||
        error("frozen table magnitude outside certificate arithmetic range")
    bound=nextfloat(Float64(envelope))
    isfinite(bound) || error("nonfinite roundoff envelope")
    FlatWindow(labels,destinations,baseline,factors,bound)
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
    contenders=Int[i for i in 1:total if prevfloat(scores[i]-flat.error_bound)<=incumbent_upper]
    length(contenders)<=max_exact_candidates || error("exact contender cap; no certificate")
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
    # Every screened-out candidate has lower bound above the initial incumbent,
    # which is >= the final exact incumbent. Exhaustion plus exact ties certify
    # the restricted optimum of this frozen table model only.
    (status="complete",certificate="exact optimum of supplied frozen binary64 response DELTA tables",
     order=flat.labels[bestorder],numerator=incumbent,score=Float64(incumbent),
     permutations=count,exact_contenders=length(contenders),minimizers=minimizers,
     roundoff_bound=flat.error_bound,destinations=length(flat.destinations),factors=length(flat.factors))
end

function insertion_window(flat::FlatWindow; guard=()->nothing,max_sweeps=10)
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
        value<score-flat.error_bound*2 || break
        current=best; score=value; moves+=1
    end
    exact=exact_numerator(flat,current)
    (order=flat.labels[current],score=Float64(exact),numerator=exact,evaluations=evaluations,moves=moves)
end
end
