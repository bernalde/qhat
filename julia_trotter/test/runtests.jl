using Test, LinearAlgebra, SparseArrays
include("../src/operator_error_norm.jl")
include("../src/energy_scope.jl")

@testset "Operator error is a singular value, not a spectral radius" begin
    U = zeros(ComplexF64, 4, 4)
    V = zeros(ComplexF64, 4, 4)
    U[1:2, 1:2] = [0 1; 1 0]
    V[1:2, 1:2] = [0 -1; 1 0]
    U[3:4, 3:4] = im * Matrix{ComplexF64}(I, 2, 2)
    V[3:4, 3:4] = Matrix{ComplexF64}(I, 2, 2)
    @test U' * U ≈ Matrix{ComplexF64}(I, 4, 4)
    @test V' * V ≈ Matrix{ComplexF64}(I, 4, 4)
    E = U - V
    @test maximum(abs, eigvals(E)) ≈ sqrt(2)
    @test opnorm(E) ≈ 2
    @test operator_error_norm(sparse(E); dense_threshold=0) ≈ 2 rtol=1e-10
    @test operator_error_norm(E; dense_threshold=4) ≈ 2
    # Tiny errors must not be collapsed by a squared-singular-value threshold.
    @test operator_error_norm(sparse(1e-20 * E); dense_threshold=0) ≈ 2e-20 rtol=1e-10
    @test operator_error_norm(spzeros(ComplexF64, 4, 4); dense_threshold=0) == 0
    @test operator_error_norm(zeros(ComplexF64, 4, 4)) == 0
    normal = Diagonal(ComplexF64[1, 2im, -3, 0])
    @test operator_error_norm(sparse(normal); dense_threshold=0) ≈ opnorm(Matrix(normal)) rtol=1e-10
    @test_throws ArgumentError operator_error_norm(E; dense_threshold=-1)
    @test_throws ArgumentError operator_error_norm(E; maxiter=0)
    @test_throws ArgumentError operator_error_norm(fill(NaN, 4, 4))
end

@testset "HF mapping and energy scope are explicit" begin
    key = "fermion-to-qubit operator mapping"
    for mapping in ("jw", "JW", "jordan-wigner", "Jordan Wigner")
        @test require_jw_hf_seed(Dict(key => mapping)) === nothing
    end
    for meta in (Dict(key => "bk"), Dict(key => "bravyi-kitaev"), Dict{String,String}())
        @test_throws ArgumentError require_jw_hf_seed(meta)
    end
    for method in (:arpack, :arnoldi)
        io = IOBuffer()
        report_energy_scope(io, method)
        output = String(take!(io))
        @test occursin("full Fock space", output)
        @test occursin("sector provenance is not verified", output)
        if method == :arnoldi
            @test occursin("not a sector projection", output)
            @test occursin("global minimum is not guaranteed", output)
        end
    end
end
