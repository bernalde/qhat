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


module ExptsCLI
include("../trotter_expts.jl")
end

const TINY_JW_FIXTURE = joinpath(@__DIR__, "fixtures", "tiny_jw.dat")

function captured_stdout(f)
    mktemp() do path, io
        redirect_stdout(io) do
            f()
        end
        flush(io)
        seekstart(io)
        return read(io, String)
    end
end

function ground_cli_rows(text)
    rows = filter(line -> occursin(r"^(1|5|10),", line), split(text, '\n'))
    return Dict(parse(Int, split(row, ',')[1]) => parse(Float64, split(row, ',')[3])
                for row in rows)
end

@testset "Accept the real generator JW metadata convention" begin
    committed = joinpath(@__DIR__, "..", "He-He", "He-He_2.40_hgbs-5_as-002-002_jw.dat")
    meta, _ = GroundCLI.parse_hamiltonian_file(committed)
    @test meta["fermion-to-qubit operator mapping"] == "jordan-wigner"
    @test GroundCLI.require_jw_hf_seed(meta) === nothing
end

@testset "Accepted two-qubit energy CLI paths" begin
    meta, ham = GroundCLI.parse_hamiltonian_file(TINY_JW_FIXTURE)
    H = Matrix(GroundCLI.build_sparse_hamiltonian(ham, 2))
    # H is a sum of two independent noncommuting one-qubit Hamiltonians.
    expected = -(sqrt(.7^2 + .3^2) + sqrt(.2^2 + .15^2))
    @test minimum(eigvals(Hermitian(H))) ≈ expected atol=1e-13 rtol=0
    @test minimum(diff(sort(eigvals(Hermitian(H))))) > .4
    @test GroundCLI.metadata_ground_energy(meta) ≈ expected atol=1e-13 rtol=0
    # Old n=4/nev=3 code selected maxdim=5, violating ArnoldiMethod's maxdim≤n.
    nvals = GroundCLI.trotter_candidate_count(4)
    @test nvals == 3
    @test nvals <= GroundCLI.trotter_krylov_dimension(4, nvals) <= 4
    arpack_text = captured_stdout(() -> GroundCLI.main(TINY_JW_FIXTURE))
    arnoldi_text = captured_stdout(() -> GroundCLI.main(TINY_JW_FIXTURE; use_arnoldi=true))
    arpack = ground_cli_rows(arpack_text)
    arnoldi = ground_cli_rows(arnoldi_text)
    @test Set(keys(arpack)) == Set([1, 5, 10])
    @test Set(keys(arnoldi)) == Set([1, 5, 10])
    @test all(isfinite, values(arpack)) && all(isfinite, values(arnoldi))
    for nsteps in (1, 5, 10)
        @test arpack[nsteps] ≈ arnoldi[nsteps] atol=2e-8 rtol=0
    end
    @test abs(arpack[10] - expected) < abs(arpack[1] - expected)
    for args in (String[], ["--arnoldi"])
        result = BondCLI.parse_energies(TINY_JW_FIXTURE, BondCLI.parse_cli(args))
        @test result.exact ≈ expected atol=1e-13 rtol=0
        @test Set(keys(result.trotter)) == Set([1, 5])
        for nsteps in (1, 5)
            @test result.trotter[nsteps] ≈ arpack[nsteps] atol=2e-8 rtol=0
        end
    end
    # Also execute the real full-space Hamiltonian fallback, not just metadata.
    mktemp() do path, io
        text = read(TINY_JW_FIXTURE, String)
        write(io, join(filter(line -> !startswith(line, "# smallest eigenvalue"),
                              split(text, '\n')), "\n"))
        close(io)
        result = BondCLI.parse_energies(path, BondCLI.parse_cli(String[]))
        @test result.exact ≈ expected atol=2e-10 rtol=0
    end
end

@testset "Demonstration reference/operator call chain on two qubits" begin
    meta, ham = ExptsCLI.parse_hamiltonian_file(TINY_JW_FIXTURE)
    norm_info = ExptsCLI.normalize_hamiltonian(meta, ham)
    terms = ExptsCLI.build_hamiltonian_terms(ham, meta;
        normalize=true, scale_by_pi=true).H_terms
    initial = ExptsCLI.construct_hf_state(2, 1)
    H_full = ExptsCLI.build_sparse_hamiltonian(ham, 2)
    H_normalized = (Hermitian(H_full) - norm_info.shift * I) / norm_info.normalization
    H_qpe = (H_normalized + .5I) * pi
    # Execute the exact changed tuple-unpacking and H_qpe wrapper call shapes.
    state, observables, state_info = ExptsCLI.chebyshev_simulation(
        initial, ExptsCLI.T_TOTAL, x -> ExptsCLI.hamiltonian_matvec(x, terms),
        pi, []; order=ExptsCLI.CHEBY_ORDER, returnstate=true)
    exact_state = exp(-im * pi * Matrix(H_normalized)) * initial
    @test norm(state - exact_state) < 2e-12
    @test isempty(observables)
    @test isfinite(state_info.error) && state_info.error > 0
    reference, info = ExptsCLI.chebyshev_exponentiation(
        H_qpe, ExptsCLI.CHEBY_ORDER, 1.0, pi; return_info=true)
    @test opnorm(Matrix(reference) - exp(-im * Matrix(H_qpe))) < 2e-12
    @test isfinite(info.error) && info.error > 0
    trotter = ExptsCLI.reference_trotter_unitary(
        ham, norm_info.normalization, 2; numsteps=5, time=pi)
    error_operator = trotter - reference
    dense_oracle = opnorm(Matrix(error_operator))
    @test dense_oracle > 1e-8
    @test ExptsCLI.operator_error_norm(error_operator) ≈ dense_oracle atol=2e-12 rtol=1e-10
    @test ExptsCLI.operator_error_norm(sparse(error_operator); dense_threshold=0) ≈
          dense_oracle atol=2e-12 rtol=1e-10
end
