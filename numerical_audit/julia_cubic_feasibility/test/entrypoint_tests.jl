using Test
module ArchiveEntrypointSmoke end
module WindowEntrypointSmoke end
module WindowBatchEntrypointSmoke end
module PhysicalEntrypointSmoke end
module CubicDPEntrypointSmoke end
module CubicCheckerEntrypointSmoke end
function contains_parse_error(value)
    value isa Expr || return false
    value.head in (:error,:incomplete) || any(contains_parse_error,value.args)
end
@testset "Study entry points parse and can be loaded as libraries" begin
    root=joinpath(@__DIR__,"..")
    entries=(("archive_study.jl",ArchiveEntrypointSmoke,:main),
        ("window_study.jl",WindowEntrypointSmoke,:main),
        ("window_batch.jl",WindowBatchEntrypointSmoke,:batch_main),
        ("physical_study.jl",PhysicalEntrypointSmoke,:main),
        ("cubic_dp_study.jl",CubicDPEntrypointSmoke,:main),
        ("independent_cubic_dp_check.jl",CubicCheckerEntrypointSmoke,:main))
    saved_args=copy(ARGS)
    try
        empty!(ARGS);push!(ARGS,"deliberately invalid library-import arguments")
        mktempdir() do directory
            cd(directory) do
                for (file,namespace,entry) in entries
                    @test !contains_parse_error(Meta.parseall(read(joinpath(root,file),String)))
                    Base.include(namespace,joinpath(root,file))
                    @test isdefined(namespace,entry)
                    # include creates methods in a newer world than this callback.
                    @test Base.invokelatest(hasmethod,getfield(namespace,entry),Tuple{Vector{String}})
                    @test ARGS==["deliberately invalid library-import arguments"]
                    @test isempty(readdir(directory))
                end
            end
        end
    finally
        empty!(ARGS);append!(ARGS,saved_args)
    end
end
