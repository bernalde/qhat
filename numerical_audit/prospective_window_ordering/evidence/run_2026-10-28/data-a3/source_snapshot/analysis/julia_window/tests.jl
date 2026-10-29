# Tiny independent identities and interruption checks; no molecular generation.
using Test
const EXACT_ACTION_PROBE = length(ARGS)==2 && ARGS[1]=="--exact-action-probe"
isempty(ARGS) || EXACT_ACTION_PROBE || error("usage: tests.jl [--exact-action-probe CORE_PATH]")
include(EXACT_ACTION_PROBE ? abspath(ARGS[2]) : joinpath(@__DIR__,"ResponseCore.jl"))
using .ResponseCore

function cancellation_unit(delta::Rational{BigInt})
    exact=Rational{BigInt}[big(1)//big(1)+delta,-big(1)//big(1)]
    Unit(UInt64(1),UInt64[0,0],Int8[1,1],Float64.(exact);exact=exact)
end
if EXACT_ACTION_PROBE
    # Root may point this one test set at a preserved pre-repair Core. Its
    # expected old-code failure is a NumericalResolution, never a large job.
    @testset "Exact action repair probe: supplied nonzero cancellation and exact zero" begin
        delta=big(1)//(big(1)<<80)
        @test unitvalue(cancellation_unit(delta),0)==2.0^-80
        @test unitvalue(cancellation_unit(-delta),0)==-2.0^-80
        @test unitvalue(cancellation_unit(big(0)//big(1)),0)===0.0
    end
    exit()
end
include(joinpath(@__DIR__,"WindowSearch.jl"))
using .WindowSearch
struct InjectedDeadline <: Exception end

function fixture()
    u(f,p,s,c)=Unit(UInt64(f),UInt64.(p),Int8.(s),Float64.(c))
    [u(0,[0],[1],[.07]),u(3,[0,3],[1,-1],[.5,.5]),u(0,[1],[1],[.3]),
     u(6,[0,6],[1,-1],[.2,.2]),u(0,[6],[1],[.13]),u(0,[1,4],[1,1],[.11,-.19])]
end
function visit_orders(f, values)
    if length(values)<=1
        f(copy(values));return
    end
    for i in eachindex(values)
        head=values[i];rest=values[setdiff(eachindex(values),[i])]
        visit_orders(rest) do tail;f([head;tail]);end
    end
end

@testset "Product tables retain full exterior interference and original costs" begin
    units=fixture();parent=[6,4,2,5,3,1];hf=UInt64(2)
    for duration in (.37,1.6), budget in (17,32,80)
        model=build_window_model(units,parent,hf,2:5;time=duration)
        @test model.setup_stats.aggregation=="real_columns_before_transport_v1"
        flat=prepare_window(model)
        baseline_only=score_order(units[parent[2:5]],collect(1:4),hf,80;time=duration)
        @test !isapprox(baseline_only.cubic_response,score_order(units,parent,hf,budget;time=duration).cubic_response;rtol=1e-4)
        visit_orders(copy(parent[2:5])) do local_order
            order=copy(parent);order[2:5]=local_order
            direct=score_order(units,order,hf,budget;time=duration)
            conditioned=window_score(model,local_order,budget)
            @test conditioned.value≈direct.value rtol=1e-10 atol=1e-26
            lengths=length.(getfield.(units,:coeff));r=direct.steps
            @test r*(2sum(lengths)-lengths[first(order)]-lengths[last(order)])+lengths[first(order)]<=budget
            positions=Int[findfirst(==(label),flat.labels) for label in local_order]
            @test Float64(exact_numerator(flat,positions))≈direct.cubic_response rtol=1e-10 atol=1e-26
        end
        candidate=insertion_window(flat;max_sweeps=3)
        order=copy(parent);order[2:5]=candidate.order
        @test first(order)==first(parent) && last(order)==last(parent)
        @test score_order(units,order,hf,budget;time=duration).value<=score_order(units,parent,hf,budget;time=duration).value+1e-26
    end
end

@testset "Independent deadline callback reaches setup/action/transport/flatten/search" begin
    units=fixture();parent=[6,4,2,5,3,1];hf=UInt64(2)
    function interrupt_after(f, limit)
        calls=Ref(0)
        callback=(kind,amount)->begin
            calls[]+=1
            calls[]>=limit && throw(InjectedDeadline())
        end
        @test_throws InjectedDeadline with_budget(f,callback)
        @test calls[]==limit
    end
    interrupt_after(()->center_units(units,hf),2)
    interrupt_after(()->apply_unit(units[2],Dict(hf=>1.0)),2)
    interrupt_after(()->score_order(units,parent,hf,80;time=.37),32)
    interrupt_after(()->diagonal_response(units,Dict(UInt64(1)=>.1),hf;time=1.6),3)
    interrupt_after(()->build_window_model(units,parent,hf,2:5;time=.37),50)
    model=build_window_model(units,parent,hf,2:5;time=.37)
    interrupt_after(()->prepare_window(model),20)
    flat=prepare_window(model)
    interrupt_after(()->insertion_window(flat),25)
    interrupt_after(()->exact_numerator(flat,[1,2,3,4]),4)
    # Callback restoration, not a permanent mutation of the scientific guard.
    @test isfinite(score_order(units,parent,hf,80).value)
end

@testset "Local resolution semantics and preallocation/storage admission" begin
    @test ResponseCore._local_resolved(0.0,1.0,true,"structural zero")===nothing
    @test_throws NumericalResolution ResponseCore._local_resolved(1e-12,1.0,true,"nonzero cancellation")
    @test ResponseCore._local_resolved(1e-12,1.0,false,"guard disabled")===nothing
    units=fixture();parent=[6,4,2,5,3,1]
    @test_throws ArgumentError build_window_model(units,parent,2,1:4)
    @test_throws ArgumentError build_window_model(units,parent,2,2:5;max_factor_entries=1)
    model=build_window_model(units,parent,2,2:5)
    @test_throws ArgumentError prepare_window(model;max_flat_entries=1)
    huge=(big(1)<<16400)//big(1)
    @test_throws ArgumentError ResponseCore._check_rational_bits(huge)
end

@testset "Exact unit actions repair floating cancellation and preserve centered targets" begin
    delta=big(1)//(big(1)<<80)
    positive=cancellation_unit(delta);negative=cancellation_unit(-delta)
    # Both supplied Float64 sums vanish, while their exact source actions do not.
    @test positive.coeff==[1.0,-1.0] && negative.coeff==[1.0,-1.0]
    @test unitvalue(positive,0)==2.0^-80
    @test unitvalue(negative,0)==-2.0^-80
    @test unitvalue(positive,0;guard=false)==0.0 # Existing explicit opt-out.
    # A different rounded input has a nonzero floating residual but exact zero.
    large=big(10)^16
    zeroexact=Rational{BigInt}[large//big(1),-(large+1)//big(1),big(1)//big(1)]
    zerounit=Unit(UInt64(1),UInt64[0,0,0],Int8[1,1,1],Float64.(zeroexact);exact=zeroexact)
    @test unitvalue(zerounit,0;guard=false)==1.0
    @test unitvalue(zerounit,0)===0.0
    tinyexact=copy(zeroexact);tinyexact[2]+=delta
    tinyunit=Unit(UInt64(1),UInt64[0,0,0],Int8[1,1,1],Float64.(tinyexact);exact=tinyexact)
    @test unitvalue(tinyunit,0;guard=false)==1.0
    @test unitvalue(tinyunit,0)==2.0^-80
    original=Unit(UInt64(0),UInt64[1,2],Int8[1,1],[1.0,-1.0];exact=copy(positive.exact))
    centered=only(center_units([original],0))
    @test unitvalue(original,0)==2.0^-80
    @test unitvalue(centered,0)===0.0
    @test unitvalue(centered,3)==-2.0^-79
    @test original.exact==positive.exact && length(original.coeff)==2
    @test length(centered.coeff)==3 # Numerical scalar does not alter original counts.
    @test affordable_steps([original],[1],2)==1
    @test_throws ArgumentError affordable_steps([original],[1],1)
end

@testset "Exact midpoint cells authenticate ties and repair double-round proposals" begin
    q(x)=Rational{BigInt}(x)
    one=1.0;odd=nextfloat(one);even=nextfloat(odd)
    midpoint1=(q(one)+q(odd))/2
    midpoint2=(q(odd)+q(even))/2
    beyond256=big(1)//(big(1)<<400)
    # Expected values are adjacent binary64 literals, never Float64(rational).
    cases=[(midpoint1-beyond256,one),(midpoint1,one),(midpoint1+beyond256,odd),
           (midpoint2-beyond256,odd),(midpoint2,even),(midpoint2+beyond256,even)]
    for (exact,expected) in cases
        @test ResponseCore._rounded_exact_action(exact)===expected
        @test ResponseCore._rounded_exact_action(-exact)===-expected
        @test ResponseCore._rounding_cell_contains(exact,expected)
    end
    @test !ResponseCore._rounding_cell_contains(midpoint1,odd)
    @test !ResponseCore._rounding_cell_contains(midpoint2,odd)
    @test ResponseCore._rounded_exact_action(q(floatmin(Float64)))===floatmin(Float64)
end

@testset "Exact action rounding refuses unsupported range and bounded arithmetic" begin
    smallest=reinterpret(Float64,UInt64(1))
    second=reinterpret(Float64,UInt64(2))
    qsmall=big(1)//(big(1)<<1074)
    @test ResponseCore._rounded_exact_action(qsmall)===smallest
    @test ResponseCore._rounded_exact_action(-qsmall)===-smallest
    @test unitvalue(cancellation_unit(qsmall),0)===smallest
    @test ResponseCore._rounded_exact_action(3qsmall/2)===second # Subnormal tie-to-even.
    @test ResponseCore._rounded_exact_action(qsmall/2+big(1)//(big(1)<<1500))===smallest
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(qsmall/2)
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(qsmall/4)
    @test_throws NumericalResolution unitvalue(cancellation_unit(qsmall/4),0)
    @test_throws NumericalResolution ResponseCore._rounded_exact_action((big(1)<<1024)//big(1))
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(Rational{BigInt}(floatmax(Float64)))
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(-Rational{BigInt}(floatmax(Float64)))
    overcap=big(1)//(big(1)<<16400)
    @test_throws ArgumentError ResponseCore._rounded_exact_action(overcap)
    @test_throws ArgumentError unitvalue(cancellation_unit(overcap),0)
    # Candidate/midpoint components fit individually, but an odd denominator
    # requiring oversized proof cross-products still refuses rather than grows.
    nearcap=(big(1)<<16383)+1
    @test_throws ArgumentError ResponseCore._rounding_compare(nearcap//(nearcap+2),
        ((big(1)<<52)+1)//(big(1)<<52))
    overflowing=Unit(UInt64(1),UInt64[0,0],Int8[1,1],[floatmax(Float64),floatmax(Float64)])
    @test_throws NumericalResolution unitvalue(overflowing,0)
end

@testset "Exact fallback interruption and downstream global guards remain mandatory" begin
    exact=(big(1)//big(1))+(big(1)//(big(1)<<400))
    for limit in (2,4,8,14)
        calls=Ref(0)
        callback=(kind,amount)->begin
            calls[]+=1
            calls[]>=limit && throw(InjectedDeadline())
        end
        @test_throws InjectedDeadline with_budget(()->ResponseCore._rounded_exact_action(exact),callback)
        @test calls[]==limit
    end
    @test unitvalue(cancellation_unit(big(1)//(big(1)<<80)),0)==2.0^-80 # Callback restored.
    @test_throws NumericalResolution ResponseCore._resolved(0.0,1.0,true)
    @test_throws NumericalResolution ResponseCore._resolved(1e-30,1.0,true)
    @test_throws NumericalResolution diagonal_response(Unit[],Dict(UInt64(0)=>1.0),0;
        guard=true,resolution_scale=1.0) # Genuine HF projection is still guarded.
    zpi=Unit(UInt64(0),UInt64[1],Int8[1],[Float64(pi)])
    @test_throws NumericalResolution diagonal_response([zpi],Dict(UInt64(1)=>1.0),0;
        time=1.0,guard=true,resolution_scale=1.0) # Genuine kernel suppression.
    x=Unit(UInt64(1),UInt64[0],Int8[1],[1.0])
    ztiny=Unit(UInt64(0),UInt64[1],Int8[1],[2.0^-80])
    @test_throws NumericalResolution score_order([x,ztiny],[1,2],0,64;time=1.0)
    @test_throws NumericalResolution ResponseCore._score(Dict(UInt64(1)=>ComplexF64(reinterpret(Float64,UInt64(1)))),1)
    @test_throws NumericalResolution ResponseCore._score(Dict(UInt64(1)=>ComplexF64(1e-100)),big(1)<<180)
end

using TOML, SHA
include(joinpath(@__DIR__,"worker.jl"))
@testset "Explicit 1024-label/8192-word/20-bit reader admission" begin
    # Parsing only: there is no large-domain numerical search in this test.
    mktempdir() do directory
        path=joinpath(directory,"admission.toml")
        data=Dict{String,Any}("version"=>1,"n"=>20,"reference"=>17,"time"=>"37/100","budget"=>40000,
            "order"=>collect(1023:-1:0),"word_lengths"=>fill(8,1024),"parent_score"=>1.0,
            "width"=>2,"max_windows"=>1,"max_sweeps"=>1,"max_entries"=>50000,
            "max_cached_actions"=>50000,"max_factor_entries"=>250000,"max_flat_entries"=>1000000,
            "max_actions"=>2000000,"max_atoms"=>200000,"max_evaluations"=>100,
            "units"=>[Dict("flip"=>0,"parities"=>collect(0:7),"signs"=>fill(1,8),
                           "coefficients"=>fill("1/100",8)) for _ in 1:1024])
        open(path,"w") do io;TOML.print(io,data);end
        admitted=read_request(path,()->nothing)
        @test length(admitted.units)==1024 && sum(length(u.coeff) for u in admitted.units)==8192
        @test admitted.hf==UInt64(17)
        push!(data["units"],first(data["units"]))
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
        pop!(data["units"])
        push!(data["units"][end]["parities"],8)
        push!(data["units"][end]["signs"],1)
        push!(data["units"][end]["coefficients"],"1/100")
        data["word_lengths"][end]=9
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing) # 8193 words
        data["n"]=21
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
        data["n"]=20;data["version"]=true
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
    end
end
