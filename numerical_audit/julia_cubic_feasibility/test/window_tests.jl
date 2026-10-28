using Test
if !isdefined(Main,:ResponseCore); include(joinpath(@__DIR__,"..","src","ResponseCore.jl")); end
if !isdefined(Main,:WindowSearch); include(joinpath(@__DIR__,"..","src","WindowSearch.jl")); end
using .ResponseCore, .WindowSearch

@testset "Exact frozen-table window certificate and enumeration" begin
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    z(a)=Unit(UInt64(0),UInt64[1],Int8[1],[a])
    units=[z(0.3),hop,z(0.5),z(0.2),hop]
    model=build_window_model(units,collect(1:5),UInt64(2),2:4;guard=false)
    flat=prepare_window(model)
    result=enumerate_window(flat)
    @test result.permutations==6
    @test sort(result.order)==[2,3,4]
    vals=[exact_numerator(flat,p) for p in ([1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1])]
    @test result.numerator==minimum(vals)
    @test result.minimizers==count(==(minimum(vals)),vals)
    insertion=insertion_window(flat)
    @test insertion.numerator>=result.numerator
    @test_throws ErrorException enumerate_window(flat;guard=()->error("deadline"))
end

@testset "Independent Heap coverage and global-label tie break" begin
    for k in 2:5
        seen=Set{Tuple{Vararg{Int}}}(); visits=Ref(0)
        WindowSearch.visit_orders(collect(1:k)) do order
            visits[]+=1;push!(seen,Tuple(order))
        end
        @test visits[]==factorial(k)
        @test length(seen)==factorial(k)
    end
    flat=WindowSearch.FlatWindow([20,10],UInt64[1],ComplexF64[1],WindowSearch.FlatFactor[],1e-14)
    result=enumerate_window(flat)
    @test result.order==[10,20]
    @test result.minimizers==2
end
