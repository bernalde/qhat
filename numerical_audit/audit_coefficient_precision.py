"""Quantify implicit OpenFermion pruning on immutable molecular tensors.

No schedule or source input is changed. Both paths use the same final 1e-12
Pauli cutoff; the explicit path performs no intermediate coefficient cutoff.
The coefficient l1 difference bounds the operator-norm difference, including
on sectors not sampled by the numerical runs. It is not a rigorous interval
bound on floating-point rounding.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import time

from openfermion import get_fermion_operator, jordan_wigner

from qhat.analysis.benchmark_L_sweep_trotter import (
    clean_fermion_operator, load_interaction_operator,
)
from qhat.analysis.fermionic_coefficients import (
    PRECISION_POLICY_VERSION, fermion_from_interaction,
    jordan_wigner_molecular, normal_ordered_molecular,
)


def difference(left, right):
    """right minus left, summed without symbolic-operator arithmetic."""
    keys = sorted(left.keys() | right.keys())
    values = {key: right.get(key, 0) - left.get(key, 0) for key in keys}
    values = {key: value for key, value in values.items() if value != 0}
    return {
        "left_terms": len(left), "right_terms": len(right),
        "new_terms": len(right.keys() - left.keys()),
        "missing_terms": len(left.keys() - right.keys()),
        "changed_terms": len(values),
        "coefficient_l1_difference": math.fsum(abs(v) for v in values.values()),
        "max_coefficient_difference": max(map(abs, values.values()), default=0),
        "largest_differences": [
            {"word": key, "difference": [value.real, value.imag]}
            for key, value in sorted(values.items(), key=lambda item: -abs(item[1]))[:5]
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, default=Path(__file__).with_name("fixtures"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest_path = args.fixtures / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    rows = []
    for case, entry in manifest.items():
        start = time.monotonic()
        tensor = args.fixtures / case / entry["tensor"]
        digest = hashlib.sha256(tensor.read_bytes()).hexdigest()
        if digest != entry["files"][entry["tensor"]]:
            raise RuntimeError(f"Changed input: {case}")
        interaction, n = load_interaction_operator(tensor)
        legacy_raw = get_fermion_operator(interaction)
        raw = fermion_from_interaction(interaction)
        legacy_normal = clean_fermion_operator(legacy_raw, 1e-12)
        precise_normal = normal_ordered_molecular(raw, atol=0)
        legacy_jw = jordan_wigner(legacy_normal)
        legacy_jw.compress(abs_tol=1e-12)
        precise_jw = jordan_wigner_molecular(raw, atol=1e-12)
        normal_then_jw = jordan_wigner_molecular(precise_normal, atol=1e-12)
        row = {
            "case": case, "n_qubits": n, "tensor_sha256": digest,
            "raw_interaction_conversion": difference(legacy_raw.terms, raw.terms),
            "normal_ordering": difference(legacy_normal.terms, precise_normal.terms),
            "identity_free_jw": difference(
                {k: v for k, v in legacy_jw.terms.items() if k},
                {k: v for k, v in precise_jw.terms.items() if k}),
            "jw_direct_vs_exact_normal_then_jw": difference(precise_jw.terms, normal_then_jw.terms),
            "seconds": time.monotonic() - start,
        }
        rows.append(row)
        print(case, json.dumps(row["identity_free_jw"]), flush=True)
    result = {
        "precision_policy": PRECISION_POLICY_VERSION,
        "final_pauli_cutoff": 1e-12, "explicit_intermediate_cutoff": 0,
        "scope": "Numerical difference only; no assertion that any added terms improve accuracy or preserve exact molecular symmetries after cutoff.",
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "source_hashes": {
            str(path.relative_to(Path(__file__).resolve().parents[1])):
                hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [Path(__file__).resolve(),
                         Path(__file__).resolve().parents[1] / "analysis" / "fermionic_coefficients.py",
                         Path(__file__).resolve().parents[1] / "analysis" / "benchmark_L_sweep_trotter.py"]
        },
        "versions": {name: importlib.metadata.version(name) for name in ["numpy", "scipy", "openfermion"]},
        "rows": rows,
    }
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
