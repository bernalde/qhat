using Test, SHA
length(ARGS)==1 || error("pass exact ResponseCore source path")
include(abspath(ARGS[1]));using .ResponseCore
@testset "Authenticated exact nonzero action survives cancellation" begin
    d=big(2)^80
    u=Unit(UInt64(1),UInt64[0,0],Int8[1,-1],Float64[1,1];exact=[(d+1)//d,big(1)//big(1)])
    got=try unitvalue(u,UInt64(0)) catch e;e end
    @test !(got isa Exception) || got isa NumericalResolution
    @test got==2.0^-80
end
println("CORE_SHA256 ",bytes2hex(sha256(read(ARGS[1]))))
