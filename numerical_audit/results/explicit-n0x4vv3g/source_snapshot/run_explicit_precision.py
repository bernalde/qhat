"""Opt-in molecular ordering benchmark with explicit coefficient precision.

Use the raw frozen tensors, CAR-normalize without an intermediate cutoff,
consolidate final JW coefficients, then apply the stated final cutoff. Exact
and approximate evolution use precisely that same retained Pauli Hamiltonian.
This is a new, versioned experiment, not a rewrite of historical data.
"""
from __future__ import annotations

import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/qhat-audit-numba")

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile
import time

import numpy as np
from scipy.sparse.linalg import expm_multiply

from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis.fermionic_coefficients import (
    PRECISION_POLICY_VERSION, fermion_from_interaction,
)
from qhat.analysis.fermionic_ownership import build_molecular_ownership
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import (
    canonical, evolve, object_sha, physical_metrics, sha, write,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).with_name("fixtures")


def prepare(case, manifest, atol):
    entry = manifest[case]
    tensor = FIXTURES / case / entry["tensor"]
    if sha(tensor) != entry["files"][entry["tensor"]]:
        raise RuntimeError(f"Changed tensor: {case}")
    interaction, n = baseline.load_interaction_operator(tensor)
    owned = build_molecular_ownership(fermion_from_interaction(interaction), atol=atol)
    coefficients = owned.final_coefficients
    if any(complex(value).imag != 0 for value in coefficients.values()):
        raise ValueError("Unitary Pauli evolution requires real coefficients; no imaginary projection is applied by this runner")
    if owned.diagnostics["noncommuting_owned_block_count"]:
        raise ValueError("A purported commuting owned block contains noncommuting Pauli factors")
    keys = list(coefficients)
    schedules = {
        "fermionic_signed": list(owned.pauli_order),
        "jw_signed": baseline.signed_coefficient_lexicographic_order(keys, coefficients, n, atol),
        "jw_magnitude": baseline.magnitude_descending_lexicographic_order(keys, coefficients, n, atol),
    }
    if owned.fallback_pauli_keys:
        raise RuntimeError("Explicit molecular ownership has unowned Pauli terms")
    if any(set(schedule) != set(keys) or len(schedule) != len(keys) for schedule in schedules.values()):
        raise RuntimeError("A schedule is not a permutation of the common Hamiltonian")
    compiled = evolution.compile_ordered_terms(keys, coefficients, n, atol)
    initial_full = evolution.build_hartree_fock_state(n, entry["electrons"])
    hf_index = int(np.argmax(abs(initial_full)))
    states, reduced, _ = ref.affine_space(compiled, hf_index)
    H, parity = ref.assemble_hamiltonian(reduced, len(states))
    initial = np.zeros(len(states), complex)
    initial[0] = 1
    raw_index = {key: i for i, key in enumerate(keys)}
    terms = {name: [reduced[raw_index[key]] for key in schedule] for name, schedule in schedules.items()}
    number_mask = np.array([int(x).bit_count() == entry["electrons"] for x in states])
    spin_sector = set(map(int, evolution.spin_sector_basis_indices(n, entry["electrons"])))
    spin_mask = np.array([int(x) in spin_sector for x in states])
    old_setup = json.loads((FIXTURES / case / "setup.json").read_text())
    diagnostics = {
        "n_qubits": n, "n_electrons": entry["electrons"],
        "pauli_terms": len(keys), "closed_dimension": len(states),
        "tensor_sha256": sha(tensor), "coefficient_policy": PRECISION_POLICY_VERSION,
        "identity_free_jw_sha256": object_sha(canonical(coefficients)),
        "physical_order_sha256": {name: object_sha(
            [[key, complex(coefficients[key]).real.hex(), complex(coefficients[key]).imag.hex()]
             for key in order]) for name, order in schedules.items()},
        "historical_pauli_terms": old_setup["pauli_terms"],
        "historical_hamiltonian_sha256": old_setup["fingerprints"]["identity_free_jw_sha256"],
        "construction": owned.diagnostics,
    }
    return H, initial, terms, number_mask, spin_mask, diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--times", nargs="+", type=float, default=[1.])
    parser.add_argument("--steps", nargs="+", type=int, default=[100])
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).with_name("results"))
    args = parser.parse_args()
    if any(t <= 0 or not np.isfinite(t) for t in args.times) or any(r <= 0 for r in args.steps):
        parser.error("Times and step counts must be positive and finite")
    if not np.isfinite(args.atol) or args.atol < 0:
        parser.error("atol must be finite and nonnegative")
    fixture_manifest = json.loads((FIXTURES / "manifest.json").read_text())
    cases = args.cases or list(fixture_manifest)
    if set(cases) - fixture_manifest.keys():
        parser.error("Unknown frozen case")
    args.output_root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="explicit-", dir=args.output_root))
    modules = [
        ROOT / "analysis" / filename for filename in (
            "fermionic_coefficients.py", "fermionic_ownership.py", "state_metrics.py",
            "benchmark_b2_active_spaces_matrix_free.py", "benchmark_b2_signed_coefficient_baseline.py",
            "benchmark_L_sweep_trotter.py", "benchmark_b2_coloring_robustness.py",
        )
    ] + [Path(__file__), Path(ref.__file__), Path(__file__).with_name("run_frozen_audit.py")]
    hashes = {str(path.relative_to(ROOT)): sha(path) for path in modules}
    write(output / "manifest.json", {
        "cases": cases, "times": args.times, "steps": args.steps,
        "final_pauli_cutoff": args.atol, "intermediate_normal_order_cutoff": 0,
        "precision_policy": PRECISION_POLICY_VERSION, "source_hashes": hashes,
        "fixture_manifest_sha256": sha(FIXTURES / "manifest.json"),
        "versions": {p: importlib.metadata.version(p) for p in ["numpy", "scipy", "openfermion", "numba"]},
        "metric": "normalized_full_state_v2; no postselection",
        "budget": "Same retained Hamiltonian, initial state, time and Pauli exponential count across the three orderings within each case.",
        "historical_comparison": "Changes coefficient arithmetic and may change exact sorting ties. Not an isolated metric correction or identical-input historical reproduction.",
    })
    print("OUTPUT", output, flush=True)
    start = time.monotonic()
    rows = []
    try:
        for case in cases:
            H, initial, schedules, nm, sm, checks = prepare(case, fixture_manifest, args.atol)
            write(output / f"{case}.checks.json", checks)
            for total_time in args.times:
                exact = expm_multiply(-1j * total_time * H, initial, traceA=-1j * total_time * H.diagonal().sum())
                for name, terms in schedules.items():
                    for steps in args.steps:
                        approximate = evolve(initial, terms, total_time, steps)
                        row = {
                            "case": case, "ordering": name, "T": total_time, "steps": steps,
                            "pauli_terms": checks["pauli_terms"],
                            "nominal_exponential_count": steps * checks["pauli_terms"],
                            **physical_metrics(exact, approximate, nm, sm),
                        }
                        rows.append(row)
                        write(output / f"{case}.{name}.T{total_time:g}.r{steps}.json", row)
                        print(f"ROW {case} {name} T={total_time:g} r={steps} E={row['one_minus_overlap']:.9g}", flush=True)
        for rel, digest in hashes.items():
            if sha(ROOT / rel) != digest:
                raise RuntimeError(f"Source changed during run: {rel}")
        write(output / "results.json", rows)
        write(output / "completion.json", {"status": "complete", "rows": len(rows), "seconds": time.monotonic() - start})
    except Exception as exc:
        write(output / "failure.json", {"type": type(exc).__name__, "message": str(exc), "completed_rows": len(rows)})
        raise
    print("COMPLETE", len(rows), output, flush=True)


if __name__ == "__main__":
    main()
