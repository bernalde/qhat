using Test
include(joinpath(@__DIR__, "..", "physical.jl"))

# Admission only: no model compilation, reference, optimization or propagation.
@testset "Single-line scientific TOML boundary" begin
    mktempdir() do folder
        path = joinpath(folder, "scientific.toml")
        text = "case = \"A\"\nword = '0Z'\nnumerators = [\"1\", \"-2\"]\n"
        text *= "# Comment quotes are ignored: " * repeat("'", 3) * " " * repeat("\"", 3) * "\n"
        text *= "escaped = \"" * repeat("\\\"", 3) * "\"\n"
        write(path, text)
        parsed = CostMatchedPhysical.bounded_toml(path)
        @test parsed.data["case"] == "A"
        @test parsed.data["word"] == "0Z"
        @test parsed.data["escaped"] == repeat("\"", 3)
        for quotechar in ("\"", "'"), value in ("A", "", "\nA\n")
            write(path, "case = " * repeat(quotechar, 3) * value * repeat(quotechar, 3) * "\n")
            @test_throws ArgumentError CostMatchedPhysical.bounded_toml(path)
        end
    end
end
