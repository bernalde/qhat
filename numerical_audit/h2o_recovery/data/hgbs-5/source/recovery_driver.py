#!/usr/bin/env python3
"""Reconstruct the 15 slide-26 H2O tensor inputs; do not run Trotter evolution.

Use QHAT's HF driver once per basis, retain its historical tensor policy, and
also export an explicitly unpruned spin-orbital expansion of the same integrals.
Completed outputs are immutable: resuming verifies provenance and every file.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.metadata
import inspect
import json
import logging
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace


GEOMETRY = [
    ("O", (0.0, 0.0, 0.0)),
    ("H", (0.757540797781, 0.0, 0.587079636589)),
    ("H", (-0.757540797781, 0.0, 0.587079636589)),
]
CASES = {
    "sto-6g": [(2, 2), (2, 4), (4, 4), (6, 4), (8, 4), (10, 4)],
    "hgbs-5": [(2, 2), (2, 4), (4, 4), (4, 6), (4, 8), (6, 8),
               (6, 10), (8, 10), (8, 12)],
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(canonical_json(value))
    temporary.replace(path)


def git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def config_path(root: Path, basis: str, occupied: int, vacant: int) -> Path:
    return (root / "hamiltonian_generator" / "polyatomic_library" / "H2O"
            / "s-1.00" / basis
            / f"H2O_s-1.00_{basis}_as-{occupied:03d}-{vacant:03d}_JW.config")


def verify_historical_config(path: Path, basis: str, occupied: int, vacant: int):
    """Read literal configuration fields without executing archived Python."""
    fields, geometry = {}, []
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name):
                fields[f"{target.value.id}.{target.attr}"] = ast.literal_eval(node.value)
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            call = node.value
            if (isinstance(call.func, ast.Attribute)
                    and isinstance(call.func.value, ast.Name)
                    and call.func.value.id == "hamiltonian"
                    and call.func.attr == "add_atom"):
                args = [ast.literal_eval(value) for value in call.args]
                geometry.append((args[0], tuple(args[1:])))
    expected = {
        "hamiltonian.basis": basis,
        "hamiltonian.num_active_occupied": occupied,
        "hamiltonian.num_active_vacant": vacant,
        "hamiltonian.f2q_mapping": "JW",
    }
    require(geometry == GEOMETRY, f"Geometry differs in archived config: {path}")
    for key, value in expected.items():
        require(fields.get(key) == value, f"Unexpected {key} in {path}")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qhat-root", type=Path, required=True)
    parser.add_argument("--historical-source", type=Path, required=True)
    parser.add_argument("--bases", nargs="+", choices=tuple(CASES), default=list(CASES))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-memory-mb", type=int, default=4096,
                        help="PySCF memory budget and full integral array size guard")
    parser.add_argument("--address-space-mb", type=int, default=12288,
                        help="Hard process virtual-memory limit; 0 disables")
    parser.add_argument("--plan-only", action="store_true",
                        help="Inspect basis sizes and provenance without running SCF")
    return parser.parse_args()


def main():
    args = parse_arguments()
    require(args.threads > 0 and args.max_memory_mb >= 256, "Invalid resource limits")
    require(args.address_space_mb >= 0, "Negative address-space limit")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = str(args.threads)
    if args.address_space_mb:
        requested = args.address_space_mb * 1024 * 1024
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        require(hard == resource.RLIM_INFINITY or requested <= hard,
                "Requested address-space limit exceeds current hard limit")
        resource.setrlimit(resource.RLIMIT_AS, (requested, hard))

    # Imports follow thread and memory controls. Import the explicitly selected
    # checkout instead of accidentally using another editable installation.
    qhat_root = args.qhat_root.resolve()
    require(qhat_root.name == "qhat", "QHAT checkout directory must be named qhat")
    sys.path.insert(0, str(qhat_root.parent))
    import numpy as np
    import openfermion.chem.molecular_data as molecular_data_module
    import openfermionpyscf._run_pyscf as ofpyscf_module
    import pyscf
    from openfermion import MolecularData
    from qhat.hamiltonian_generator import hamgen

    require(Path(hamgen.__file__).resolve().is_relative_to(qhat_root),
            "Imported QHAT outside requested checkout")
    pyscf.lib.num_threads(args.threads)
    historical_root = args.historical_source.resolve()
    config_records = {}
    for basis, pairs in CASES.items():
        for occupied, vacant in pairs:
            path = config_path(historical_root, basis, occupied, vacant)
            verify_historical_config(path, basis, occupied, vacant)
            config_records[str(path.relative_to(historical_root))] = digest(path)

    distributions = sorted(
        ((dist.metadata["Name"], dist.version)
         for dist in importlib.metadata.distributions()), key=lambda row: row[0].lower()
    )
    source_paths = {
        "recovery_driver": Path(__file__).resolve(),
        "qhat_hamgen": Path(hamgen.__file__).resolve(),
        "qhat_hamgen_types": qhat_root / "hamiltonian_generator" / "hamgen_types.py",
        "openfermion_molecular_data": Path(molecular_data_module.__file__).resolve(),
        "openfermionpyscf_driver": Path(ofpyscf_module.__file__).resolve(),
        "historical_qhat_hamgen": historical_root / "hamiltonian_generator" / "hamgen.py",
        "historical_qhat_hamgen_types": historical_root / "hamiltonian_generator" / "hamgen_types.py",
    }
    source_hashes = {key: digest(path) for key, path in source_paths.items()}
    archived_ast = ast.parse(source_paths["historical_qhat_hamgen"].read_text())
    archived_functions = {node.name: node for node in archived_ast.body
                          if isinstance(node, ast.FunctionDef)}
    identical_functions = {}
    for name in ("my_run_pyscf", "apply_active_space", "define_active_orbitals"):
        current_node = ast.parse(inspect.getsource(getattr(hamgen, name))).body[0]
        identical_functions[name] = ast.dump(current_node) == ast.dump(archived_functions[name])
        require(identical_functions[name], f"QHAT function {name} differs from archived source")
    provenance = {
        "schema_version": 1,
        "purpose": "reconstruction; original H2O tensors are unavailable",
        "geometry_angstrom": GEOMETRY,
        "charge": 0, "multiplicity": 1, "total_electrons": 10,
        "orbital_convention": "PySCF canonical RHF spatial orbitals; alpha/beta interleaved",
        "orbital_gauge": "Unmodified PySCF eigenvector signs; MO coefficients retained",
        "historical_qhat_commit": git_head(historical_root),
        "qhat_commit": git_head(qhat_root),
        "source_sha256": source_hashes,
        "historical_generator_function_ast_identical": identical_functions,
        "historical_configuration_sha256": config_records,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": distributions,
        "threads": args.threads,
        "max_memory_mb": args.max_memory_mb,
        "address_space_mb": args.address_space_mb,
        "historical_spinorb_integral_cutoff": float(molecular_data_module.EQ_TOLERANCE),
    }
    # Commit IDs remain evidence, but an unrelated/new commit with unchanged
    # generation source bytes must not prevent a safe resume.
    def compatible_provenance(value):
        return {key: item for key, item in value.items() if key != "qhat_commit"}

    fingerprint = hashlib.sha256(canonical_json(compatible_provenance(provenance)).encode()).hexdigest()

    plans = {}
    for basis in args.bases:
        probe = MolecularData(GEOMETRY, basis, 1, 0, filename="unused_probe")
        mol = ofpyscf_module.prepare_pyscf_molecule(probe)
        nao = int(mol.nao_nr())
        full_integral_bytes = 8 * nao ** 4
        require(full_integral_bytes < args.max_memory_mb * 1024 * 1024 / 3,
                f"{basis}: full four-index integral array needs {full_integral_bytes} bytes; "
                "exceeds one third of selected PySCF memory budget")
        plans[basis] = {
            "spatial_orbitals": nao, "full_spatial_eri_bytes": full_integral_bytes,
            "active_spin_orbital_counts": [a + b for a, b in CASES[basis]],
            "resolved_pyscf_basis": mol._basis,
            "pyscf_cartesian_basis": bool(mol.cart),
            "pyscf_geometry_unit": mol.unit,
        }
        require(int(mol.nelectron) == 10 and int(mol.spin) == 0,
                "Expected neutral singlet H2O")
    print(canonical_json({"provenance_sha256": fingerprint, "plan": plans}), flush=True)
    if args.plan_only:
        return

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        require(manifest.get("provenance_sha256") == fingerprint,
                "Existing output provenance differs; use a new output directory")
        require(compatible_provenance(manifest.get("provenance", {}))
                == json.loads(canonical_json(compatible_provenance(provenance))),
                "Existing provenance body does not match its fingerprint")
    else:
        require(not any(output.iterdir()), "Output exists without a manifest and is not empty")
        manifest = {"provenance_sha256": fingerprint, "provenance": provenance, "bases": {}}
        write_json(manifest_path, manifest)

    for basis, entry in manifest["bases"].items():
        basis_dir = output / basis
        actual_files = {str(path.relative_to(basis_dir)) for path in basis_dir.rglob("*")
                        if path.is_file()}
        optional = entry.get("optional_cache_sha256", {})
        require(actual_files - set(optional) == set(entry["files_sha256"]),
                f"Stored file set differs for completed basis {basis}")
        for relative, expected_hash in entry["files_sha256"].items():
            require(digest(basis_dir / relative) == expected_hash,
                    f"Stored output changed: {basis}/{relative}")
        for relative, expected_hash in optional.items():
            if (basis_dir / relative).exists():
                require(digest(basis_dir / relative) == expected_hash,
                        f"Stored optional cache changed: {basis}/{relative}")

    for basis in args.bases:
        if basis in manifest["bases"]:
            print(f"Verified complete basis {basis}; skipping SCF", flush=True)
            continue
        final_dir = output / basis
        partial = output / f".{basis}.partial"
        require(not partial.exists() and not final_dir.exists(),
                f"Incomplete/unregistered output exists for {basis}; inspect it before retrying")
        partial.mkdir()
        (partial / "configs").mkdir()
        (partial / "historical_policy").mkdir()
        (partial / "unpruned").mkdir()
        write_json(partial / "basis_definition.json", plans[basis])
        # Freeze the source files that govern this reconstruction, not just paths.
        (partial / "source").mkdir()
        for key, source in source_paths.items():
            shutil.copyfile(source, partial / "source" / f"{key}.py")
        log = logging.getLogger(f"h2o_recovery.{basis}")
        log.setLevel(logging.INFO)
        log_handler = logging.FileHandler(partial / "generation.log")
        log.addHandler(log_handler)
        state = SimpleNamespace(log=log.info, log_verbose=log.info)
        molecule = MolecularData(
            GEOMETRY, basis, 1, 0, filename=str(partial / "hf_molecule")
        )
        original_compute_scf = hamgen.compute_scf
        original_compute_integrals = hamgen.compute_integrals

        def bounded_scf(mol):
            result = original_compute_scf(mol)
            result.max_memory = args.max_memory_mb
            return result

        def checked_integrals(mol, scf):
            require(bool(scf.converged), "RHF did not converge; refusing to generate tensors")
            require(np.array_equal(np.asarray(scf.mo_occ),
                                   np.r_[np.full(5, 2.0), np.zeros(len(scf.mo_occ) - 5)]),
                    "Unexpected RHF occupied orbital pattern")
            return original_compute_integrals(mol, scf)

        hamgen.compute_scf = bounded_scf
        hamgen.compute_integrals = checked_integrals
        started = time.monotonic()
        old_cwd = Path.cwd()
        try:
            # MolecularData.save uses a temporary file in cwd. Keep that local.
            os.chdir(partial)
            hf = hamgen.my_run_pyscf(state, molecule)
        finally:
            os.chdir(old_cwd)
            hamgen.compute_scf = original_compute_scf
            hamgen.compute_integrals = original_compute_integrals
        hf.hf_time = time.monotonic() - started
        hf.basis = basis
        hf.separation = None
        scf = hf._pyscf_data["scf"]
        overlap = scf.get_ovlp()
        orthogonality_error = float(np.max(np.abs(scf.mo_coeff.T @ overlap @ scf.mo_coeff
                                                  - np.eye(hf.n_orbitals))))
        require(orthogonality_error < 1e-8, "Canonical orbitals are not S-orthonormal")
        np.savez_compressed(
            partial / "orbitals.npz", mo_coeff=scf.mo_coeff,
            mo_energy=scf.mo_energy, mo_occ=scf.mo_occ,
            ao_overlap=overlap, ao_core_hamiltonian=scf.get_hcore(),
        )
        summary = {
            "basis": basis, "hf_energy_hartree": hf.hf_energy,
            "nuclear_repulsion_hartree": hf.nuclear_repulsion,
            "scf_converged": bool(scf.converged),
            "scf_cycles": int(scf.cycles), "scf_conv_tol": float(scf.conv_tol),
            "scf_conv_tol_grad": scf.conv_tol_grad,
            "scf_max_cycle": int(scf.max_cycle),
            "hf_wall_seconds": hf.hf_time,
            "orbital_orthogonality_max_error": orthogonality_error,
            "spatial_orbitals": int(hf.n_orbitals),
            "historical_spinorb_integral_cutoff": float(molecular_data_module.EQ_TOLERANCE),
            "unpruned_policy": "No integral or spin-orbital coefficient cutoff in export",
            "cases": [],
        }

        for occupied, vacant in CASES[basis]:
            config = config_path(historical_root, basis, occupied, vacant)
            shutil.copyfile(config, partial / "configs" / config.name)
            stem = f"H2O_s-1.00_{basis}_as-{occupied:03d}-{vacant:03d}"
            config_text = (
                "# Reconstructed configuration. Coordinates are Angstrom.\n"
                "# This native QHAT replay writes its standard historical-policy tensors.\n"
                f"general.file_stub = 'H2O_s-1.00_{basis}'\n"
                f"general.output_directory = params['output_directory']\n"
                f"general.logfile = '{stem}.log'\n"
                "general.print_verbose()\n"
                + "".join(f"hamiltonian.add_atom({atom!r}, {xyz[0]!r}, {xyz[1]!r}, {xyz[2]!r})\n"
                          for atom, xyz in GEOMETRY)
                + f"hamiltonian.basis = {basis!r}\n"
                f"hamiltonian.num_active_occupied = {occupied}\n"
                f"hamiltonian.num_active_vacant = {vacant}\n"
                "hamiltonian.f2q_mapping = 'JW'\n"
            )
            (partial / "configs" / f"{stem}.relocated.py").write_text(config_text)
            state.config_hamiltonian = SimpleNamespace(
                num_active_occupied=occupied, num_active_vacant=vacant
            )
            historical = hamgen.apply_active_space(state, hf)
            frozen_spatial = (10 - occupied) // 2
            active_spatial = list(range(frozen_spatial,
                                        frozen_spatial + (occupied + vacant) // 2))
            core, spatial_one, spatial_two = hf.get_active_space_integrals(
                occupied_indices=list(range(frozen_spatial)), active_indices=active_spatial
            )
            n = occupied + vacant
            one = np.zeros((n, n))
            two = np.zeros((n, n, n, n))
            for sigma in (0, 1):
                one[sigma::2, sigma::2] = spatial_one
                for tau in (0, 1):
                    two[sigma::2, tau::2, tau::2, sigma::2] = spatial_two / 2
            constant = float(hf.nuclear_repulsion + core)
            hermiticity_one = float(np.max(np.abs(one - one.T)))
            hermiticity_two = float(np.max(np.abs(two - two.transpose(3, 2, 1, 0))))
            require(max(hermiticity_one, hermiticity_two) < 1e-9,
                    f"Tensor Hermiticity failed for {stem}")
            # Independent Slater-determinant contraction: <a†p a†q ar as>
            # = delta_ps delta_qr - delta_pr delta_qs for occupied p,q.
            def determinant_energy(h, v):
                return float(constant + sum(h[p, p] for p in range(occupied))
                             + sum(v[p, q, q, p] - v[p, q, p, q]
                                   for p in range(occupied) for q in range(occupied)))

            unpruned_hf = determinant_energy(one, two)
            require(abs(unpruned_hf - hf.hf_energy) < 1e-8,
                    f"Independent HF contraction failed for {stem}")
            cutoff = float(molecular_data_module.EQ_TOLERANCE)
            expected_one, expected_two = one.copy(), two.copy()
            expected_one[np.abs(expected_one) < cutoff] = 0.0
            expected_two[np.abs(expected_two) < cutoff / 2] = 0.0
            require(np.array_equal(historical.one_body_tensor, expected_one)
                    and np.array_equal(historical.two_body_tensor, expected_two)
                    and historical.constant == constant,
                    f"Historical QHAT spin-expansion policy differs for {stem}")
            delta_one = one - historical.one_body_tensor
            delta_two = two - historical.two_body_tensor
            truncation_bound = float(np.abs(delta_one).sum() + np.abs(delta_two).sum())
            historical_hf = determinant_energy(historical.one_body_tensor,
                                               historical.two_body_tensor)
            require(abs(historical_hf - unpruned_hf) <= truncation_bound + 1e-10,
                    f"Historical cutoff HF change exceeds coefficient-l1 bound for {stem}")
            np.savez_compressed(
                partial / "historical_policy" / f"{stem}.tensors.npz",
                constant=historical.constant,
                one_body=historical.one_body_tensor,
                two_body=historical.two_body_tensor,
            )
            np.savez_compressed(
                partial / "unpruned" / f"{stem}.tensors.npz",
                constant=constant, one_body=one, two_body=two,
            )
            summary["cases"].append({
                "case_id": stem, "active_electrons": occupied,
                "active_spin_orbitals": n, "active_spatial_orbitals": active_spatial,
                "frozen_doubly_occupied_spatial_orbitals": list(range(frozen_spatial)),
                "hf_computational_basis_index_big_endian": ((1 << occupied) - 1) << vacant,
                "hermiticity_one_max_error": hermiticity_one,
                "hermiticity_two_max_error": hermiticity_two,
                "unpruned_determinant_energy_hartree": unpruned_hf,
                "unpruned_hf_error": abs(unpruned_hf - hf.hf_energy),
                "historical_determinant_energy_hartree": historical_hf,
                "historical_hf_error": abs(historical_hf - hf.hf_energy),
                "matches_openfermion_cutoff_formula_bitwise": True,
                "cutoff_changed_one_entries": int(np.count_nonzero(delta_one)),
                "cutoff_changed_two_entries": int(np.count_nonzero(delta_two)),
                "cutoff_coefficient_l1_difference": truncation_bound,
                "cutoff_max_coefficient_difference": float(max(np.max(np.abs(delta_one)),
                                                                np.max(np.abs(delta_two)))),
            })
            print(f"Validated {stem}: {n} qubits, {occupied} electrons; "
                  f"HF residual {abs(unpruned_hf - hf.hf_energy):.3g}", flush=True)

        # Absolute HDF5 filename inside stored molecular metadata referred to the
        # partial directory during generation; loading by explicit filename works
        # after rename, and the portable tensor/orbital NPZ files have no paths.
        summary["max_rss_kb_process"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(partial / "summary.json", summary)
        log.removeHandler(log_handler)
        log_handler.close()
        optional_cache = {"hf_molecule.hdf5": digest(partial / "hf_molecule.hdf5")}
        files = {str(path.relative_to(partial)): digest(path)
                 for path in sorted(partial.rglob("*"))
                 if path.is_file() and str(path.relative_to(partial)) not in optional_cache}
        partial.rename(final_dir)
        manifest["bases"][basis] = {
            "files_sha256": files, "optional_cache_sha256": optional_cache,
            "case_count": len(CASES[basis]),
        }
        write_json(manifest_path, manifest)
        print(f"Completed and hashed {basis}: {len(CASES[basis])} cases", flush=True)
    print(f"Manifest: {manifest_path}", flush=True)


if __name__ == "__main__":
    main()
