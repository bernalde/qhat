using Test, LinearAlgebra, SparseArrays

module ReferenceBackend
include("../src/densesimulators.jl")
end

@testset "Scaled Chebyshev reference and opt-in metadata" begin
    H = ComplexF64[pi/2 + pi/8 pi/4; pi/4 pi/2 - pi/8]
    exact = exp(-im * H)
    legacy = ReferenceBackend.chebyshev_exponentiation(Hermitian(sparse(H)), 10, 1.0, pi)
    reference, info = ReferenceBackend.chebyshev_exponentiation(
        Hermitian(sparse(H)), 10, 1.0, pi; return_info=true)
    @test legacy isa AbstractMatrix
    @test reference ≈ legacy atol=2e-13 rtol=0
    # A numerical parity tolerance, not a certified error allowance.
    @test opnorm(Matrix(reference) - exact) < 2e-13
    @test isfinite(info.error) && 0 < info.error < 1e-12
    initial = ComplexF64[1, 0]
    state, _, state_info = ReferenceBackend.chebyshev_simulation(
        initial, 1.0, x -> H*x, pi, []; order=10, returnstate=true)
    @test norm(state - exact*initial) < 2e-13
    @test isfinite(state_info.error) && 0 < state_info.error < 1e-12
end

module GroundCLI
include("../trotter_ground_energy.jl")
end
module BondCLI
include("../bond_energy_curve.jl")
end

@testset "Reject unsupported seed before molecular work" begin
    @test_throws ArgumentError BondCLI.parse_cli(["--encoding", "bk", "--arnoldi"])
    @test BondCLI.parse_cli(["--encoding", "bk"]).method == :arpack
    @test BondCLI.parse_cli(["--encoding", "jw", "--arnoldi"]).method == :arnoldi
    # Deliberately omit numerical metadata: rejection must precede reconstruction.
    meta = Dict("fermion-to-qubit operator mapping" => "bravyi-kitaev")
    @test_throws ArgumentError GroundCLI.trotter_energy_arnoldi(meta, Dict{String,ComplexF64}(), 1)
    mktemp() do path, io
        write(io, "# fermion-to-qubit operator mapping = bravyi-kitaev\n")
        close(io)
        @test_throws ArgumentError GroundCLI.main(path; use_arnoldi=true)
        opts = BondCLI.parse_cli(["--arnoldi"])
        @test_throws ArgumentError BondCLI.parse_energies(path, opts)
    end
end
