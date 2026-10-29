using Test,SHA
length(ARGS)==1 || error("pass immutable ResponseCore source path")
include(abspath(ARGS[1]));using .ResponseCore
@testset "Exact commuting local atoms survive floating cancellation" begin
    t=big(1)//(big(1)<<27)
    ar=1+t;br=1-t
    a=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],fill(Float64(ar/2),2);exact=fill(ar/2,2))
    b=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],fill(Float64(br/2),2);exact=fill(br/2,2))
    @test ResponseCore._all_commuting([a,b])
    for order in ((1,2),(2,1))
        observed=try
            atom_column([a,b],collect(order),1)
        catch error
            error
        end
        @test !(observed isa Exception) || observed isa NumericalResolution
        @test observed isa AbstractDict && isempty(observed)
    end
end
println("CORE_SHA256 ",bytes2hex(sha256(read(ARGS[1]))))
