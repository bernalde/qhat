using Test
module WindowEntrypointSmoke end
function contains_parse_error(value)
    value isa Expr || return false
    value.head in (:error,:incomplete) || any(contains_parse_error,value.args)
end
@testset "Study entry points parse and can be loaded as libraries" begin
    root=joinpath(@__DIR__,"..")
    for file in ("archive_study.jl","window_study.jl","window_batch.jl","physical_study.jl")
        @test !contains_parse_error(Meta.parseall(read(joinpath(root,file),String)))
    end
    Base.include(WindowEntrypointSmoke,joinpath(root,"window_study.jl"))
    @test isdefined(WindowEntrypointSmoke,:main)
end
