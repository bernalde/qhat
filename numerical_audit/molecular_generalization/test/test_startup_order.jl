module StartupOrderTests
using Test
include(joinpath(@__DIR__,"..","anneal.jl"))
const W=MolecularAnnealingWorker
const A=W.FullOrderAnnealing

@testset "Four admitted orders compare by value" begin
    for size in (4,286)
        orders=[circshift(collect(1:size),shift) for shift in 0:3]
        saved=deepcopy(orders)
        @test A._distinct_four_orders(orders)
        @test orders==saved
        for a in 1:3, b in a+1:4
            duplicate=copy.(orders)
            duplicate[b]=copy(duplicate[a])
            @test duplicate[a]!==duplicate[b]
            @test duplicate[a]==duplicate[b]
            @test !A._distinct_four_orders(duplicate)
            # The tuple oracle is deliberately restricted to tiny orders.
            if size==4
                @test A._distinct_four_orders(duplicate)==(length(Set(Tuple(p) for p in duplicate))==4)
            end
        end
        @test orders==saved
    end
    @test !A._distinct_four_orders([collect(1:4) for _ in 1:3])
    @test !A._distinct_four_orders([collect(1:4) for _ in 1:5])
end

@testset "Integer-vector lexicographic semantics" begin
    tiny=[Int[],[1],[1,2],[1,3],[0,10],[1,2,0],[-1,3]]
    for a in tiny, b in tiny
        @test A._lexless(a,b)==isless(Tuple(a),Tuple(b))
    end
    order=collect(1:286)
    @test !A._lexless(order,copy(order))
    for position in (1,143,286)
        greater=copy(order);greater[position]+=1
        smaller=copy(order);smaller[position]-=1
        @test A._lexless(order,greater)
        @test !A._lexless(greater,order)
        @test A._lexless(smaller,order)
        @test !A._lexless(order,smaller)
    end
    prefix=order[1:end-1]
    @test A._lexless(prefix,order)
    @test !A._lexless(order,prefix)
    @test order==collect(1:286)
end

@testset "Observed worker and numerical module settings" begin
    # Exported in pinned Julia 1.12.6 include/julia/julia.h. The getter walks
    # parent modules for inherited settings; it reports actual module state,
    # not the independent command-line optimization/minimum clamps.
    module_level(m::Module)=ccall(:jl_get_module_optlevel,Cint,(Any,),m)
    @test module_level(W)==0
    @test module_level(A)==2
    @test module_level(A.ResponseCore)==2
end
end # module
