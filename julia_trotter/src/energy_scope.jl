# Scope/encoding checks shared by the energy CLIs and the HF-seeded backend.
function require_jw_hf_seed(meta::AbstractDict)
    raw = get(meta, "fermion-to-qubit operator mapping", "")
    mapping = replace(lowercase(strip(raw)), "_" => "-", " " => "-")
    mapping in ("jw", "jordan-wigner") || throw(ArgumentError(
        "HF-seeded Arnoldi requires explicit Jordan-Wigner metadata; " *
        "mapping '$raw' is unsupported. Use the full-space Arpack method for BK " *
        "or unknown encodings."))
    return nothing
end

function report_energy_scope(io::IO, method::Symbol)
    method in (:arpack, :arnoldi) || throw(ArgumentError("Unknown energy method: $method"))
    println(io, "WARNING: energies use full Fock space, without a fixed electron-number or spin-sector restriction. " *
                "Metadata energies are used as supplied; their sector provenance is not verified.")
    if method == :arnoldi
        println(io, "WARNING: the JW Hartree-Fock start is only a Krylov seed, not a sector projection. " *
                    "Its overlap/invariant subspace can affect the extrema found; agreement with a global minimum is not guaranteed.")
    end
    return nothing
end
