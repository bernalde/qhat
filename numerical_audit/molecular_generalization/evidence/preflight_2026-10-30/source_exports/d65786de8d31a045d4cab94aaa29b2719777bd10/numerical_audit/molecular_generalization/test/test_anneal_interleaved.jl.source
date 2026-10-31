# Controlled CPU-clock scheduling tests; synthetic scores are not numerical evidence.
using Test, Random
length(ARGS)<=1 || error("at most one explicit kernel path")
include(isempty(ARGS) ? joinpath(@__DIR__,"..","src","FullOrderAnnealing.jl") : abspath(only(ARGS)))
const A=FullOrderAnnealing

mutable struct ControlledContext
    units::Vector{Int}
    stats::Dict{String,Int}
    now::Base.RefValue{Float64}
    calls::Int
    probes::Vector{Vector{Int}}
    costs::Vector{Float64}
    values::Vector{Float64}
end
function controlled_context(costs;values=ones(length(costs)))
    ControlledContext(collect(1:5),Dict("unit_actions"=>0,"full_scores"=>0),Ref(0.0),0,
        Vector{Vector{Int}}(),Float64.(costs),Float64.(values))
end
function A.full_score!(c::ControlledContext,order;check=A.noop)
    c.calls+=1
    check()
    value=1.0
    if c.calls>4
        push!(c.probes,copy(order));index=length(c.probes)
        c.now[]+=c.costs[min(index,length(c.costs))]
        c.stats["unit_actions"]+=1
        check(:unit_action)
        value=c.values[min(index,length(c.values))]
    end
    c.stats["full_scores"]+=1
    (value=value,cubic_response=value,steps=1)
end

# Observe only the new cooling fraction; delegate the unchanged formula through
# the original generic method. Never load this specialization in a real worker.
const COOLING_FRACTIONS=Float64[]
function A.temperature(settings::Dict{String,Any},fraction::Float64)
    push!(COOLING_FRACTIONS,fraction)
    invoke(A.temperature,Tuple{Any,Any},settings,fraction)
end
draft_settings(;hot=false)=merge(copy(A.SETTINGS_DEFAULTS),Dict("calibrated"=>true,
    "calibration_sha256"=>repeat("a",64),"theta_start"=>(hot ? 1e30 : .5),
    "theta_end"=>(hot ? 1e29 : .005)))
draft_starts()=[[1,2,3,4,5],[2,1,3,4,5],[1,3,4,2,5],[1,4,2,3,5]]

function controlled_run(c;deadline=100.0,stop_after=nothing,hot=false)
    empty!(COOLING_FRACTIONS)
    attempts=Ref(0);saved=Any[]
    check=(kind=:loop,amount=1)->begin
        if kind==:proposal
            attempts[]+=1
            if stop_after!==nothing && attempts[]>stop_after
                throw(A.Limit("cpu_limit","injected earlier whole-call cutoff"))
            end
        end
        nothing
    end
    state=A.anneal(c,draft_starts(),draft_settings(hot=hot);seed_material="interleaved-draft:1",
        cpu_deadline=deadline,clock=()->c.now[],check=check,checkpoint=s->push!(saved,s))
    state,saved,copy(COOLING_FRACTIONS)
end

function assert_counter_prefix(state)
    completed=sum(state["chain_proposals"])
    @test state["evaluations"]==state["start_scores_completed"]+completed
    @test completed<=state["proposals"]<=completed+4
    @test state["accepted_moves"]==sum(state["chain_accepted"])<=completed
    @test all(state["chain_accepted"].<=state["chain_proposals"])
    @test state["uphill_accepted"]<=min(state["uphill_attempts"],state["accepted_moves"])
    @test all(>=(0),state["chain_CPU_seconds"])
    @test sort(state["order"])==collect(1:5)
end

println("START interleaved controller: least actual CPU, ties and RNG ownership");flush(stdout)
@testset "least CPU and original-index ties survive an early global stop" begin
    context=controlled_context([4.,1.,1.,1.,1.])
    state,saved,cooling=controlled_run(context;stop_after=5)
    # First quantums: chains1,2,3,4. Next: chain2, tied for least spent.
    selected=[1,2,3,4,2]
    current=copy.(draft_starts());rngs=[Xoshiro(seed) for seed in A._seeds("interleaved-draft:1")]
    expected=Vector{Vector{Int}}()
    for chain in selected
        order=current[chain];rng=rngs[chain]
        source=rand(rng,1:length(order));destination=rand(rng,1:length(order)-1)
        destination>=source && (destination+=1)
        proposal=A.insertion(order,source,destination)
        uniform=rand(rng);while iszero(uniform);uniform=rand(rng);end
        current[chain]=proposal # Equal synthetic scores are accepted.
        push!(expected,copy(proposal))
    end
    @test context.probes==expected
    @test state["status"]=="cpu_limit"
    @test state["start_scores_completed"]==4
    @test state["chain_proposals"]==[1,2,1,1]
    @test state["chain_CPU_seconds"]==[4.,2.,1.,1.]
    @test cooling≈[0.,0.,0.,0.,1/25]
    # Durable first-round checkpoints contain all four completed chains before
    # the injected whole-call cutoff, without claiming a new best is necessary.
    @test any(snapshot->all(>(0),snapshot["chain_proposals"]),saved)
    foreach(assert_counter_prefix,saved)
    assert_counter_prefix(state)
    second=controlled_context([4.,1.,1.,1.,1.])
    repeated,_,_=controlled_run(second;stop_after=5)
    @test repeated==state && second.probes==context.probes
end

println("START interleaved controller: interrupted quantums are not retried");flush(stdout)
@testset "share expiry and final global interruption preserve partial counters" begin
    context=controlled_context([1.1])
    state,saved,_=controlled_run(context;deadline=4.) # Fixed per-chain share=1.
    @test state["status"]=="cpu_limit"
    @test length(context.probes)==4
    @test state["proposals"]==4
    @test state["chain_proposals"]==[0,0,0,0]
    @test state["evaluations"]==4 # Only the four starts completed.
    @test state["chain_CPU_seconds"]≈fill(1.1,4)
    # Cooperative checks can see an overrun, which is retained, not clipped or
    # turned into a completed proposal; real whole-call admission remains separate.
    @test sum(state["chain_CPU_seconds"])>4.
    foreach(assert_counter_prefix,saved)
    assert_counter_prefix(state)
end

println("START interleaved controller: uphill moves retain the global best");flush(stdout)
@testset "separate chain state and owned checkpoints preserve best score" begin
    context=controlled_context([1.,1.];values=[.5,2.])
    state,saved,_=controlled_run(context;stop_after=2,hot=true)
    @test state["accepted_moves"]==2
    @test state["uphill_attempts"]==state["uphill_accepted"]==1
    @test state["score"]==.5
    @test state["order"]==context.probes[1]
    @test all(saved[i+1]["score"]<=saved[i]["score"] for i in 1:length(saved)-1)
    @test state["seeds"]==string.(A._seeds("interleaved-draft:1"))
    saved[end]["order"][1]=999
    @test sort(state["order"])==collect(1:5)
    assert_counter_prefix(state)
end

# Before adoption, also run the existing real scorer/wire/native bridge tests,
# injected evaluation/action/support/numerical limits, callback restoration and
# checkpoint serialization CPU attribution. Then measure complete cold old-only
# B32/B64 constructors and require four scored starts, a completed proposal in
# every chain, native parity and strict complete-cost eligibility. These
# synthetic-controller tests cannot establish any of those performance gates.
