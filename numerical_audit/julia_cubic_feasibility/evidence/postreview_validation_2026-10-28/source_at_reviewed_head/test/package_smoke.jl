using CubicOrderingFeasibility, Test
const C=CubicOrderingFeasibility
@testset "Public Julia package loads from the exact project" begin
    @test normpath(Base.pathof(C))==normpath(joinpath(@__DIR__,"..","src","CubicOrderingFeasibility.jl"))
    factor=C.StructureDP.LocalFactor([1,2],Dict((1,2)=>-3,(2,1)=>7))
    solved=C.StructureDP.bag_order_dp([1,2],[factor],[[1,2]],Tuple{Int,Int}[])
    @test solved.complete && solved.cost==-3 && solved.order==[1,2]
    input=C.StudyIO.load_model(only(ARGS))
    @test input.case=="CO_r1p13" && input.n==12 && length(input.units)==195 && input.words==770
    @test C.ResponseCore.affordable_steps(input.units,input.order,30input.words)==15
end
println("Package path: ",Base.pathof(C))
