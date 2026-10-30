using Test, LinearAlgebra, TOML, SHA
include(joinpath(@__DIR__,"..","anneal.jl"))
const W=MolecularAnnealingWorker
const A=W.FullOrderAnnealing
const R=A.ResponseCore

function fixtures()
    units=[A.Unit(UInt64(1),UInt64[0],Int8[1],[.75]),
           A.Unit(UInt64(0),UInt64[1],Int8[1],[-.5]),
           A.Unit(UInt64(0),UInt64[2],Int8[1],[.375]),
           A.Unit(UInt64(1),UInt64[2],Int8[1],[.25]),
           A.Unit(UInt64(2),UInt64[0,1],Int8[1,1],[.125,.1875])]
    x=[0. 1.;1. 0.];z=[1. 0.;0. -1.];id=Matrix{Float64}(I,2,2)
    matrices=[.75kron(id,x),-.5kron(id,z),.375kron(z,id),.25kron(z,x),
              .125kron(x,id)+.1875kron(x,z)]
    units,matrices
end
function permutations(values)
    isempty(values) && return [Int[]]
    [vcat([values[i]],tail) for i in eachindex(values) for tail in permutations(values[setdiff(eachindex(values),[i])])]
end
comm(a,b)=a*b-b*a
function nested_cubic(matrices,order)
    result=zeros(size(first(matrices)))
    for (position,index) in enumerate(order)
        suffix=sum((matrices[j] for j in order[position+1:end]);init=zeros(size(result)))
        a=matrices[index]
        result+=-comm(a,comm(a,suffix))/24+comm(suffix,comm(suffix,a))/12
    end
    result
end
function dense_s3(matrices,order)
    d=size(first(matrices),1);q=[Matrix{Float64}(I,d,d),zeros(d,d),zeros(d,d),zeros(d,d)]
    stages=vcat([(i,.5) for i in order[1:end-1]],[(last(order),1.)],[(i,.5) for i in reverse(order[1:end-1])])
    for (i,tau) in stages
        out=copy.(q)
        for degree in 0:2, exponent in 1:3-degree
            out[degree+exponent+1]+=(tau*matrices[i])^exponent/factorial(exponent)*q[degree+1]
        end
        q=out
    end
    q[4]-sum(matrices)^3/6
end
settings()=merge(copy(A.SETTINGS_DEFAULTS),Dict("calibrated"=>true,"calibration_sha256"=>repeat("a",64),
    "theta_start"=>.5,"theta_end"=>.005))
starts()=[[1,2,3,4,5],[2,1,3,4,5],[1,3,4,2,5],[1,4,2,3,5]]

@testset "Complete full-order S3, transport, projection and endpoint fusion" begin
    units,matrices=fixtures();budget=180
    context=A.FullOrderContext(units,2,0,budget)
    counts=Set{Int}()
    diagonal=diag(matrices[2]+matrices[3])
    for order in permutations(collect(1:5))
        expected=nested_cubic(matrices,order)
        @test dense_s3(matrices,order)≈expected atol=3e-14 rtol=1e-11
        column=A.full_column(context,order).column
        @test [get(column,UInt64(j),0.) for j in 0:3]≈expected[:,1] atol=3e-14 rtol=1e-11
        score=A.full_score!(context,order)
        response=0.
        for target in 2:4
            phase=(diagonal[target]-diagonal[1])/2
            kernel=(iszero(phase) ? 1. : sin(phase)/phase)*cis(phase)
            response+=abs2(kernel*expected[target,1])
        end
        # Literal repeated chronological stages count every boundary fusion.
        lengths=[length(u.coeff) for u in units]
        rotations(r)=begin
            pulses=Int[]
            for _ in 1:r
                for i in vcat(order[1:end-1],[last(order)],reverse(order[1:end-1]))
                    (isempty(pulses) || last(pulses)!=i) && push!(pulses,i)
                end
            end
            sum(lengths[i] for i in pulses)
        end
        @test rotations(score.steps)<=budget<rotations(score.steps+1)
        @test score.cubic_response≈response atol=1e-24 rtol=1e-9
        @test score.value≈response/score.steps^4 atol=1e-26 rtol=1e-9
        push!(counts,score.steps)
    end
    @test length(counts)>1
    @test A.insertion([1,2,3,4,5],1,5)==[2,3,4,5,1]
    @test A.insertion(A.insertion([1,2,3,4,5],1,5),5,1)==collect(1:5)
end

@testset "Caches preserve values, scratch ownership, exact scalar metadata" begin
    units,_=fixtures();cached=A.FullOrderContext(units,2,0,180)
    uncached=A.FullOrderContext(units,2,0,180;max_cached_actions=0,cache_bytes=0)
    order=[5,3,1,4,2]
    first=A.full_column(cached,order);saved=copy(first.column)
    for candidate in starts()
        a=A.full_score!(cached,candidate);b=A.full_score!(uncached,candidate)
        @test a==b
    end
    @test first.column==saved
    @test isempty(uncached.actions) && isempty(uncached.kernels)
    @test uncached.stats["cache_skipped"]>0 && cached.stats["cache_hits"]>0
    @test cached.accounted_cache_bytes<=cached.cache_bytes
    q=(big(2)^54+1)//big(2)^54
    exact=A.Unit(UInt64(0),UInt64[0,1],Int8[1,1],[1.,-1.];exact=[q,-big(1)//big(1)])
    @test R.unitvalue(exact,UInt64(0))==2.0^-54
    shifted=copy(units)
    shifted[2]=A.Unit(UInt64(0),UInt64[1,0],Int8[1,1],[-.5,1e8])
    shifted_context=A.FullOrderContext(shifted,2,0,210)
    shifted_column=A.full_column(shifted_context,order).column
    baseline_column=A.full_column(cached,order).column
    @test [get(shifted_column,UInt64(i),0.) for i in 0:3]≈[get(baseline_column,UInt64(i),0.) for i in 0:3]
    commuting=[A.Unit(UInt64(0),UInt64[1],Int8[1],[.5])]
    @test A.full_score!(A.FullOrderContext(commuting,1,0,30),[1]).value==0.
end

@testset "Finite annealing, complete incumbents, CPU slices and seed replay" begin
    @test A.metropolis(-1.,.1,.9)
    @test A.metropolis(0.,.1,.9)
    @test A.metropolis(1.,1.,.1)
    @test !A.metropolis(1.,1.,.9)
    @test_throws ArgumentError A.metropolis(1.,1.,0.)
    @test A.temperature(settings(),0.)≈.5
    @test A.temperature(settings(),1.)≈.005
    function run()
        units,_=fixtures();context=A.FullOrderContext(units,2,0,180)
        now=Ref(0.);saved=Any[]
        check=(kind=:loop,amount=1)->begin;kind==:proposal && (now[]+=1.);nothing;end
        state=A.anneal(context,starts(),settings();seed_material="tiny:1",cpu_deadline=40.,
            clock=()->now[],check=check,checkpoint=s->push!(saved,s))
        state,saved
    end
    state,saved=run();other,_=run()
    @test state==other
    @test state["start_scores_completed"]==4
    @test all(>(0),state["chain_proposals"])
    @test state["evaluations"]==4+sum(state["chain_proposals"])
    @test all(saved[i+1]["score"]<=saved[i]["score"] for i in 1:length(saved)-1)
    @test all(sort(s["order"])==collect(1:5) && s["steps"]>0 for s in saved)
    saved[1]["order"][1]=100
    @test sort(state["order"])==collect(1:5)
    @test_throws ArgumentError A.validate_settings(A.SETTINGS_DEFAULTS)
end

@testset "Refusal paths and isolated callbacks" begin
    units,_=fixtures()
    @test_throws ArgumentError A.FullOrderContext(units,2,0,179)
    @test_throws ArgumentError A.FullOrderContext(units,2,0,180;time=.5)
    @test_throws ArgumentError A.FullOrderContext(units,15,0,180)
    bad=copy(units);bad[1]=A.Unit(UInt64(1),UInt64[0],Int8[1],[.75];exact=[big(1)//big(2)])
    @test_throws ArgumentError A.FullOrderContext(bad,2,0,180)
    @test_throws ArgumentError A.FullOrderContext(vcat(units,[units[1]]),2,0,210)
    @test_throws ArgumentError A.FullOrderContext(units,2,0,180;max_entries=0)
    @test_throws ArgumentError A.FullOrderContext(units,2,0,180;max_entries=1)
    @test_throws ArgumentError A.Unit(UInt64(1),UInt64[0],Int8[1],[Inf])
    context=A.FullOrderContext(units,2,0,180)
    @test_throws ArgumentError A.full_score!(context,[1,2,3,4,4])
    @test_throws A.Limit A.full_score!(context,collect(1:5);check=(args...)->throw(A.Limit("cpu_limit","injected")))
    @test A.full_score!(context,collect(1:5)).value>0
    @test_throws ArgumentError A.metropolis(NaN,.1,.5)
    @test_throws ArgumentError A.anneal(context,fill(collect(1:5),4),settings();seed_material="tiny:1",cpu_deadline=1.)
    tiny=[A.Unit(UInt64(1),UInt64[0],Int8[1],[1e-200]),A.Unit(UInt64(0),UInt64[1],Int8[1],[1e-200])]
    tiny_context=A.FullOrderContext(tiny,1,0,60)
    @test_throws A.NumericalResolution A.full_score!(tiny_context,[1,2])
end

function wire_request(units)
    context=A.FullOrderContext(units,2,0,180);parent=A.full_score!(context,collect(1:5))
    rawunits=[Dict("flip"=>Int(u.flip),"parities"=>Int.(u.parity),"signs"=>Int.(u.signs),
        "coefficients"=>[string(numerator(q))*"/"*string(denominator(q)) for q in u.exact]) for u in units]
    Dict{String,Any}("version"=>1,"n"=>2,"reference"=>0,"time"=>"1/1","budget"=>180,
        "order"=>collect(0:4),"starts"=>[p.-1 for p in starts()],"word_lengths"=>[length(u.coeff) for u in units],
        "parent_score"=>parent.value,"parent_cubic_response"=>parent.cubic_response,"parent_steps"=>parent.steps,
        "max_entries"=>50000,"max_cached_actions"=>100000,"seed_material"=>"tiny:1",
        "seeds"=>string.(A._seeds("tiny:1")),"settings_sha256"=>repeat("b",64),
        "coefficients_sha256"=>W.coefficient_hash(rawunits),"settings"=>settings(),"units"=>rawunits,
        "source_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(W.ROOT,name)))) for name in W.SOURCE_NAMES))
end
@testset "Wire coefficients/source identity and atomic complete checkpoint" begin
    units,_=fixtures();data=wire_request(units)
    mktempdir() do directory
        requestpath=joinpath(directory,"request.toml")
        write_request()=open(requestpath,"w") do io;TOML.print(io,data;sorted=true);end
        write_request();request=W.read_request(requestpath)
        @test request.units[5].exact==units[5].exact
        state=Dict{String,Any}("order"=>collect(1:5),"score"=>data["parent_score"],
            "cubic_response"=>data["parent_cubic_response"],"steps"=>data["parent_steps"])
        checkpoint=joinpath(directory,"checkpoint.toml")
        W.atomic_checkpoint(checkpoint,request,state,Dict{String,Float64}())
        @test TOML.parsefile(checkpoint)["order"]==collect(0:4)
        @test !ispath(checkpoint*".tmp")
        bad=copy(state);bad["steps"]=0
        @test_throws ArgumentError W.atomic_checkpoint(checkpoint,request,bad,Dict{String,Float64}())
        bad=copy(state);bad["score"]=Inf
        @test_throws ArgumentError W.atomic_checkpoint(checkpoint,request,bad,Dict{String,Float64}())
        data["source_sha256"][first(W.SOURCE_NAMES)]=repeat("0",64)
        write_request();@test_throws ArgumentError W.read_request(requestpath)
        data=wire_request(units);data["coefficients_sha256"]=repeat("0",64)
        write_request();@test_throws ArgumentError W.read_request(requestpath)
    end
end
