module WindowSearch
using ..ResponseCore: WindowModel, check_budget
export insertion_window, prepare_window, exact_numerator

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
        check_budget()
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
    for value in allvalues
        check_budget()
        isfinite(value) && (iszero(value)||1e-140<=abs(value)<=1e140) ||
            throw(ArgumentError("frozen table magnitude outside exact arithmetic range"))
    end
    nothing
end

# Recompute from the current supplied arrays; a public FlatWindow may carry a
# manually supplied or stale hint, which cannot establish a certificate itself.
function _roundoff_bound(flat::FlatWindow)
    _validate_flat_for_screen(flat)
    baseline=flat.baseline;factors=flat.factors;destinations=flat.destinations
    envelope=setprecision(BigFloat,128) do
        setrounding(BigFloat,RoundUp) do
            re=BigFloat[];im=BigFloat[]
            for value in baseline
                check_budget()
                push!(re,abs(BigFloat(real(value))));push!(im,abs(BigFloat(imag(value))))
            end
            for f in factors, (i,at) in enumerate(f.destinations)
                check_budget()
                re[at]+=maximum(abs(BigFloat(real(v))) for v in f.values[i,:])
                im[at]+=maximum(abs(BigFloat(imag(v))) for v in f.values[i,:])
            end
            h2=BigFloat(0)
            for value in Iterators.flatten((re,im))
                check_budget()
                h2+=value*value
            end
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

function prepare_window(model::WindowModel; max_flat_entries::Int=1_000_000)
    check_budget()
    labels=copy(model.baseline_window)
    2<=length(labels)<=9 || throw(ArgumentError("window width outside 2:9 cap"))
    length(unique(labels))==length(labels) || throw(ArgumentError("duplicate window labels"))
    labelmap=Dict(v=>i for (i,v) in enumerate(labels))
    entries=length(model.baseline_response)
    entries<=model.max_entries || throw(ArgumentError("baseline support cap"))
    for factor in model.factors
        check_budget()
        arity=length(factor.scope)
        2<=arity<=3 || throw(ArgumentError("response factor arity must be two or three"))
        length(unique(factor.scope))==arity || throw(ArgumentError("duplicate response factor labels"))
        all(v->haskey(labelmap,v),factor.scope) || throw(ArgumentError("factor scope outside window"))
        length(factor.variants)==factorial(arity) || throw(ArgumentError("incomplete response order variants"))
        for key in keys(factor.variants)
            check_budget()
            entries+=length(factor.variants[key])
            entries<=250000+model.max_entries || throw(ArgumentError("factor entry cap before flatten"))
            length(key)==arity && Set(key)==Set(factor.scope) || throw(ArgumentError("invalid response order variant"))
        end
    end
    reached=Set{UInt64}()
    function admit_vector!(vector)
        check_budget()
        for destination in keys(vector)
            check_budget()
            if !(destination in reached)
                length(reached)<model.max_entries || throw(ArgumentError("flat destination cap"))
                push!(reached,destination)
            end
        end
    end
    admit_vector!(model.baseline_response)
    for factor in model.factors, vector in values(factor.variants)
        check_budget()
        admit_vector!(vector)
    end
    destinations=sort!(collect(reached));check_budget()
    dmap=Dict{UInt64,Int}();baseline=ComplexF64[]
    for (i,destination) in enumerate(destinations)
        check_budget()
        dmap[destination]=i
        push!(baseline,get(model.baseline_response,destination,0.0im))
    end
    factors=FlatFactor[]
    flat_entries=length(baseline)
    for factor in model.factors
        check_budget()
        scope=Int[labelmap[v] for v in factor.scope]
        local_destinations=Set{Int}()
        for vector in values(factor.variants), destination in keys(vector)
            check_budget()
            push!(local_destinations,dmap[destination])
        end
        dest=sort!(collect(local_destinations));check_budget()
        entries=length(dest)*(length(scope)==2 ? 2 : 8)
        flat_entries += entries
        flat_entries <= max_flat_entries || throw(ArgumentError("flat storage cap exceeded"))
        check_budget()
        table=zeros(ComplexF64,length(dest),length(scope)==2 ? 2 : 8)
        valid=length(scope)==2 ? (1,2) : (1,2,4,5,7,8)
        seen=Set{Int}()
        for (key,vector) in factor.variants
            check_budget()
            pos=zeros(Int,length(labels))
            for (i,label) in enumerate(key);pos[labelmap[label]]=i;end
            column=pattern(scope,pos)
            column in valid && !(column in seen) || throw(ArgumentError("invalid or duplicate response comparison column"))
            push!(seen,column)
            for (i,at) in enumerate(dest)
                check_budget()
                table[i,column]=get(vector,destinations[at],0.0im)
            end
        end
        length(seen)==length(valid) || throw(ArgumentError("incomplete response comparison columns"))
        push!(factors,FlatFactor(scope,dest,table))
    end
    flat=FlatWindow(labels,destinations,baseline,factors,0.0)
    FlatWindow(labels,destinations,baseline,factors,_roundoff_bound(flat))
end

function response!(buffer,flat::FlatWindow,order,positions)
    check_budget(:evaluation)
    for (i,label) in enumerate(order); positions[label]=i; end
    copyto!(buffer,flat.baseline)
    for f in flat.factors
        check_budget()
        column=pattern(f.scope,positions)
        @inbounds for (i,at) in enumerate(f.destinations)
            check_budget()
            buffer[at]+=f.values[i,column]
        end
    end
    score=0.0
    @inbounds for v in buffer
        check_budget()
        score+=real(v)*real(v)+imag(v)*imag(v)
    end
    isfinite(score) && score>=0 || error("invalid frozen table score")
    check_budget(:flat_candidate)
    score
end

function exact_numerator(flat::FlatWindow,order)
    check_budget(:evaluation)
    positions=zeros(Int,length(flat.labels))
    for (i,label) in enumerate(order); positions[label]=i; end
    re=Rational{BigInt}[];im=Rational{BigInt}[]
    for value in flat.baseline
        check_budget()
        push!(re,Rational{BigInt}(real(value)));push!(im,Rational{BigInt}(imag(value)))
    end
    for f in flat.factors
        check_budget()
        column=pattern(f.scope,positions)
        for (i,at) in enumerate(f.destinations)
            check_budget()
            re[at]+=Rational{BigInt}(real(f.values[i,column]))
            im[at]+=Rational{BigInt}(imag(f.values[i,column]))
        end
    end
    result=big(0)//big(1)
    for value in Iterators.flatten((re,im))
        check_budget()
        # Binary64 tables in [1e-140,1e140] and <=120 factors keep dyadic
        # squared sums comfortably inside this fixed intermediate bit cap.
        result+=value*value
        ndigits(abs(numerator(result));base=2)<=8192 && ndigits(denominator(result);base=2)<=8192 ||
            throw(ArgumentError("exact witness bit cap"))
    end
    check_budget(:exact_witness)
    result
end

function insertion_window(flat::FlatWindow; guard=()->nothing,max_sweeps=10)
    bound=max(flat.error_bound,_roundoff_bound(flat))
    k=length(flat.labels); current=collect(1:k)
    positions=zeros(Int,k); buffer=copy(flat.baseline)
    score=response!(buffer,flat,current,positions); evaluations=1; moves=0
    for sweep in 1:max_sweeps
        check_budget()
        best=copy(current); value=score
        for from in 1:k, to in 1:k
            check_budget()
            from==to && continue
            guard()
            trial=copy(current); insert!(trial,to,popat!(trial,from))
            next=response!(buffer,flat,trial,positions); evaluations+=1
            if next<value-bound*2
                value=next;best=trial
            elseif next<=value+bound*2
                # Overlapping floating screening intervals are resolved using
                # exact dyadic evaluation of these supplied tables only.
                guard()
                if exact_numerator(flat,trial)<exact_numerator(flat,best)
                    value=next;best=trial
                end
            end
        end
        best==current && break
        guard()
        exact_next=exact_numerator(flat,best)
        exact_current=exact_numerator(flat,current)
        exact_next < exact_current || break
        current=best; score=value; moves+=1
        check_budget(:insertion_move)
    end
    guard()
    exact=exact_numerator(flat,current)
    (order=flat.labels[current],score=Float64(exact),numerator=exact,evaluations=evaluations,moves=moves)
end
end
