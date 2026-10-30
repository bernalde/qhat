"""Exact centered-parity cache identities and numerical refusal regressions."""
module ScalarSignatureCacheTests
using Test
length(ARGS)<=1 || error("at most one explicit kernel path")
include(isempty(ARGS) ? joinpath(@__DIR__,"..","src","FullOrderAnnealing.jl") : abspath(only(ARGS)))
const A=FullOrderAnnealing
const R=A.ResponseCore

context(units,n;hf=0,kwargs...)=A.FullOrderContext(units,n,hf,30sum(length(u.coeff) for u in units);kwargs...)
samebits(a::Float64,b::Float64)=reinterpret(UInt64,a)==reinterpret(UInt64,b)
phase(mask,state)=isodd(count_ones(mask&state))

@testset "Basis spans all masks, including the absolute anchor" begin
    fixtures=[(4,UInt64[],0),(4,UInt64[0,0],0),
              (4,UInt64[2],1),(4,UInt64[1,3,2,0],2),
              (4,UInt64[8,4,2,1],4)]
    for (n,masks,rank) in fixtures
        original_masks=copy(masks)
        basis=A._parity_basis(masks,n)
        @test masks==original_masks
        @test basis==A._parity_basis(masks,n)
        @test length(basis)==rank
        representatives=Dict{UInt64,UInt64}()
        for value in 0:(1<<n)-1
            state=UInt64(value);signature=A._scalar_signature(basis,state)
            @test signature<UInt64(1)<<rank
            if haskey(representatives,signature)
                other=representatives[signature]
                @test all(phase(mask,state)==phase(mask,other) for mask in masks)
            else
                representatives[signature]=state
            end
        end
        @test length(representatives)==1<<rank
    end
    # A relative-mask-only implementation would incorrectly collapse this unit.
    unit=A.Unit(UInt64(1),UInt64[2],Int8[1],[.75])
    c=context([unit],3)
    @test length(c.scalar_bases[1])==1
    @test A._scalar_signature(c.scalar_bases[1],UInt64(0))!=A._scalar_signature(c.scalar_bases[1],UInt64(2))
    @test A._value(c,1,UInt64(0))==.75
    @test A._value(c,1,UInt64(2))==-.75
    high=A._parity_basis(UInt64[8192,3,8195],14)
    @test length(high)==2
    for state in UInt64[0,1,2,3,8192,8193,8194,8195]
        for other in UInt64[0,1,2,3,8192,8193,8194,8195]
            if A._scalar_signature(high,state)==A._scalar_signature(high,other)
                @test all(phase(mask,state)==phase(mask,other) for mask in UInt64[8192,3,8195])
            end
        end
    end
    @test_throws ArgumentError A._parity_basis(UInt64[16],4)
    @test_throws ArgumentError A._parity_basis(UInt64[0],15)
    @test_throws ArgumentError A._scalar_signature(fill(UInt64(1),15),UInt64(0))
end

@testset "Centered diagonals and identities retain bitwise guarded actions" begin
    diagonal=A.Unit(UInt64(0),UInt64[1,2,3],Int8[1,1,1],[.5,.25,.125])
    identity=A.Unit(UInt64(0),UInt64[0],Int8[1],[.25])
    c=context([diagonal,identity],3;hf=1)
    @test length(c.scalar_bases[1])==2
    @test isempty(c.centered[2].coeff) && isempty(c.scalar_bases[2])
    for i in eachindex(c.centered), value in 0:7
        state=UInt64(value)
        @test samebits(A._value(c,i,state),R.unitvalue(c.centered[i],state;guard=true))
    end
    @test samebits(A._value(c,1,UInt64(1)),0.)
    @test length(c.actions)==5 # Four diagonal signatures plus centered identity.
    @test sum(length,c.scalar_bases)<=3length(c.units)
    @test c.accounted_cache_bytes==1024(length(c.actions)+length(c.kernels))
end

@testset "More than 64 terms pack basis bits, not term bits" begin
    exact=fill(big(1)//big(128),128)
    unit=A.Unit(UInt64(128),UInt64.(0:127),fill(Int8(1),128),fill(1/128,128);exact=exact)
    c=context([unit],8)
    @test length(c.scalar_bases[1])==7
    queries=Ref(0)
    callback=(kind,amount)->begin;kind==:unit_action && (queries[]+=amount);nothing;end
    R.with_budget(callback) do
        for value in 0:255
            state=UInt64(value)
            @test samebits(A._value(c,1,state),R.unitvalue(unit,state;guard=true))
        end
    end
    @test queries[]==c.stats["unit_actions"]==256
    @test c.stats["cache_misses"]==128
    @test c.stats["cache_hits"]==128
    @test length(c.actions)==128
end

@testset "Exact cancellation, certified zero and arithmetic refusal" begin
    exact=(big(2)^54+1)//big(2)^54
    unit=A.Unit(UInt64(2),UInt64[0,1],Int8[1,1],[1.,-1.];exact=[exact,-big(1)//big(1)])
    c=context([unit],2)
    for value in 0:3
        @test samebits(A._value(c,1,UInt64(value)),R.unitvalue(unit,UInt64(value);guard=true))
    end
    @test samebits(A._value(c,1,UInt64(0)),2.0^-54)
    zero_unit=A.Unit(UInt64(2),UInt64[0,1],Int8[1,1],[1.,-1.])
    zero_context=context([zero_unit],2)
    @test samebits(A._value(zero_context,1,UInt64(0)),0.)
    @test samebits(A._value(zero_context,1,UInt64(2)),0.)
    @test length(zero_context.actions)==1
    too_small=(big(2)^1200+1)//big(2)^1200
    refusing=A.Unit(UInt64(2),UInt64[0,1],Int8[1,1],[1.,-1.];exact=[too_small,-big(1)//big(1)])
    refused=context([refusing],2)
    for state in UInt64[0,2]
        @test_throws A.NumericalResolution R.unitvalue(refusing,state;guard=true)
        @test_throws A.NumericalResolution A._value(refused,1,state)
        @test isempty(refused.actions)
    end
    @test refused.stats["unit_actions"]==2 && refused.stats["cache_misses"]==2
    @test A._value(refused,1,UInt64(1))==2.
end

@testset "Disabled and saturated caches preserve values and counters" begin
    unit=A.Unit(UInt64(1),UInt64[2],Int8[1],[.75])
    for options in ((max_cached_actions=0,), (cache_bytes=0,), (cache_bytes=1023,))
        c=context([unit],3;options...)
        for value in 0:7
            @test samebits(A._value(c,1,UInt64(value)),R.unitvalue(unit,UInt64(value);guard=true))
        end
        @test isempty(c.actions) && isempty(c.kernels)
        @test c.stats["cache_misses"]==c.stats["cache_skipped"]==c.stats["unit_actions"]==8
    end
    saturated=context([unit],3;max_cached_actions=1,cache_bytes=1024)
    for state in UInt64[0,2,4,6]
        @test samebits(A._value(saturated,1,state),R.unitvalue(unit,state;guard=true))
    end
    @test length(saturated.actions)==1 && saturated.accounted_cache_bytes==1024
    @test saturated.stats["cache_hits"]==1 && saturated.stats["cache_misses"]==3
    @test saturated.stats["cache_skipped"]==2 && saturated.stats["unit_actions"]==4
end

@testset "Signature keys survive saturation and later admission-limit changes" begin
    # Centered Z on bit one is Z1-I: state2 has signature1/value-2, whereas
    # actual state1 has signature0/value0. A fallback to raw-state key1 would
    # incorrectly return the existing cached -2 after saturation/cap mutation.
    diagonal=A.Unit(UInt64(0),UInt64[2],Int8[1],[1.])
    for change in (:filled,:count,:bytes,:both)
        c=context([diagonal],3;max_cached_actions=1,cache_bytes=1024)
        @test A._value(c,1,UInt64(2))==-2.
        @test c.actions[(1,UInt64(1))]==-2.
        @test length(c.actions)==1 && c.accounted_cache_bytes==1024
        change in (:count,:both) && (c.max_cached_actions=0)
        change in (:bytes,:both) && (c.cache_bytes=0)
        @test samebits(A._value(c,1,UInt64(1)),0.)
        @test samebits(A._value(c,1,UInt64(5)),R.unitvalue(c.centered[1],UInt64(5);guard=true))
        @test A._value(c,1,UInt64(6))==-2.
        @test c.actions[(1,UInt64(1))]==-2.
        @test length(c.actions)==1 && c.accounted_cache_bytes==1024
        @test c.stats["unit_actions"]==4
        @test c.stats["cache_misses"]==3 && c.stats["cache_hits"]==1
        @test c.stats["cache_skipped"]==2
    end
end

@testset "Unit identities, actual XOR destinations and full-state transport" begin
    units=[A.Unit(UInt64(1),UInt64[0],Int8[1],[.75]),
           A.Unit(UInt64(2),UInt64[0],Int8[1],[.5])]
    separate=context(units,2)
    @test A._value(separate,1,UInt64(0))==.75
    @test A._value(separate,2,UInt64(0))==.5
    @test length(separate.actions)==2
    scalar=A.Unit(UInt64(1),UInt64[2],Int8[1],[.75])
    diagonal=A.Unit(UInt64(0),UInt64[4],Int8[1],[.25])
    c=context([scalar,diagonal],3)
    out=A.Column();input=A.Column(UInt64(0)=>1.,UInt64(4)=>2.)
    A._apply!(out,c,1,input)
    @test out==A.Column(UInt64(1)=>.75,UInt64(5)=>1.5)
    @test length(c.actions)==1
    @test A._kernel(c,UInt64(1))!=A._kernel(c,UInt64(5))
    @test Set(keys(c.kernels))==Set(UInt64[1,5])
end

@testset "Interrupted setup, outer action and scalar misses never cache a value" begin
    previous=R.BUDGET_CHECK[]
    loops=Ref(0)
    setup=(kind,amount)->begin
        kind==:loop && (loops[]+=amount)
        # One with_budget entry check, then one check per input mask. Stop on
        # the second mask, independently of its number of pivot operations.
        loops[]==3 && throw(A.Limit("cpu_limit","basis fixture"))
        nothing
    end
    @test_throws A.Limit R.with_budget(()->A._parity_basis(UInt64[1,2,3],3),setup)
    @test R.BUDGET_CHECK[]===previous
    unit=A.Unit(UInt64(1),UInt64[2],Int8[1],[.75])
    for interrupt in (:outer,:scalar)
        c=context([unit],3)
        callback=(kind,amount)->begin
            ((interrupt==:outer && kind==:unit_action) || (interrupt==:scalar && kind==:action)) &&
                throw(A.Limit("cpu_limit","miss fixture"))
            nothing
        end
        @test_throws A.Limit R.with_budget(()->A._value(c,1,UInt64(0)),callback)
        @test R.BUDGET_CHECK[]===previous
        @test isempty(c.actions) && c.stats["unit_actions"]==(interrupt==:outer ? 0 : 1)
        @test A._value(c,1,UInt64(0))==.75
    end
end

@testset "Full columns and scores agree with uncached and independent recurrence" begin
    units=[A.Unit(UInt64(1),UInt64[0],Int8[1],[.75]),
           A.Unit(UInt64(0),UInt64[1],Int8[1],[-.5]),
           A.Unit(UInt64(0),UInt64[2],Int8[1],[.375]),
           A.Unit(UInt64(1),UInt64[2],Int8[1],[.25]),
           A.Unit(UInt64(2),UInt64[0,1],Int8[1,1],[.125,.1875])]
    cached=context(units,2)
    uncached=context(units,2;max_cached_actions=0,cache_bytes=0)
    for order in ([1,2,3,4,5],[5,1,2,3,4],[2,1,3,4,5],[1,3,4,2,5])
        column=A.full_column(cached,order).column
        @test column==A.full_column(uncached,order).column
        reference=R.cubic_column(units,order,0;guard=true)
        @test [get(column,UInt64(i),0.) for i in 0:3]≈[get(reference,UInt64(i),0.) for i in 0:3] atol=3e-14 rtol=1e-10
        actual=A.full_score!(cached,order)
        @test actual==A.full_score!(uncached,order)
        independent=R.score_order(units,order,0,180;guard=true)
        @test actual.steps==independent.steps
        @test actual.cubic_response≈independent.cubic_response atol=1e-24 rtol=1e-8
        @test actual.value≈independent.value atol=1e-24 rtol=1e-8
    end
end

# The separate old-only adoption gate, not performed by this test file:
# retain the original 3144 bytes and fixed LiH requests; compare fixed-order
# cubic columns, steps/raw response/J against the unmodified scorer and native
# SparseResponse with the existing floor/refusal policy. Retain calibrated
# temperatures and every failed cost/producer identity. Then run fresh complete
# APIs at both 32/64 total CPU budgets, with startup/basis/JIT/polling/checks and
# parent prefix charged, and require four starts plus a scored proposal in each
# chain. Original/native/dense references stay unchanged; no new errors are read.
end # module
