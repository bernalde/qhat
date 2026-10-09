"""Independent conventional-Pauli replay and dense transport quadrature.

Reads the frozen experiment; never imports QHAT or run.py.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import json
from pathlib import Path
import hashlib
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import eigh

ROOT = Path("/tmp/ordering-diagonal-response/results")
OUT = Path(__file__).with_name("diagonal_response_independent_checks.json")


def action(word, basis, n):
    index = {int(v): i for i, v in enumerate(basis)}
    destination, phase = [], []
    for v in basis:
        w, a = int(v), 1.0 + 0.0j
        for q, label in word:
            bit = 1 << (n - 1 - int(q))
            occupied = bool(int(v) & bit)
            if label == "Y":
                a *= -1j if occupied else 1j
            if label == "Z" and occupied:
                a *= -1
            if label in ("X", "Y"):
                w ^= bit
        destination.append(index[w])
        phase.append(a)
    return np.array(destination), np.array(phase)


def block_matrix(raw, den, basis, n):
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    columns = np.arange(len(basis))
    for word, numerator in raw:
        destination, phase = action(word, basis, n)
        out[destination, columns] += int(numerator) / den * phase
    return out


def stable_error(exact, approximate):
    exact = exact / np.linalg.norm(exact)
    approximate = approximate / np.linalg.norm(approximate)
    residual = approximate - exact * np.vdot(exact, approximate)
    infidelity = float(np.vdot(residual, residual).real)
    return infidelity / (1 + np.sqrt(max(0.0, 1 - infidelity)))


def transport_quadrature(H, K, initial):
    values, vectors = eigh(H)
    initial_eigen = vectors.conj().T @ initial
    K_eigen = vectors.conj().T @ K @ vectors
    def integrand(t):
        return vectors @ (np.exp(1j*t*values) * (
            K_eigen @ (np.exp(-1j*t*values)*initial_eigen)
        ))
    response, bound = quad_vec(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12)
    centered = response - initial*np.vdot(initial, response)
    return float(np.vdot(centered, centered).real), float(bound)


def main():
    diagnostics = []
    replay = []
    rows = json.loads((ROOT / "results.json").read_text())
    for case in json.loads((ROOT / "manifest.json").read_text())["cases"]:
        d = json.loads((ROOT / (case + ".design.json")).read_text())
        den, n = int(d["denominator"]), d["n_qubits"]
        basis = d["affine_basis_indices"]
        with np.load(ROOT / (case + ".states.npz")) as stored:
            initial, exact = stored["initial"], stored["exact"]
            H = block_matrix(d["coefficients"], den, basis, n)
            energies, eigenvectors = eigh(H)
            independent_exact = eigenvectors @ (
                np.exp(-1j*energies)*(eigenvectors.conj().T @ initial)
            )
            exact_distance = float(np.linalg.norm(independent_exact-exact))
            if exact_distance > 1e-11:
                raise RuntimeError("Independent Pauli reference mismatch")
            for row in (r for r in rows if r["case"] == case):
                order = d["orders"][row["ordering"]]
                rotations = []
                for block_index in order:
                    for word, numerator in d["blocks"]["class"][block_index]:
                        destination, phase = action(word, basis, n)
                        angle = int(numerator)/den/row["steps"]
                        rotations.append((destination, phase, np.cos(angle), np.sin(angle)))
                state = initial.copy()
                for _ in range(row["steps"]):
                    for destination, phase, cosine, sine in rotations:
                        pstate = np.empty_like(state)
                        pstate[destination] = phase*state
                        state = cosine*state - 1j*sine*pstate
                key = f"{row['ordering']}__r{row['steps']}__{row['role']}"
                distance = float(np.linalg.norm(state-stored[key]))
                metric_difference = abs(stable_error(independent_exact, state)-row["one_minus_overlap"])
                if distance > 2e-11 or metric_difference > 1e-13:
                    raise RuntimeError("Independent circuit replay mismatch")
                replay.append(dict(case=case, ordering=row["ordering"], steps=row["steps"],
                                   state_distance=distance, metric_difference=metric_difference,
                                   independent_exact_distance=exact_distance))
        if case not in ("CH2O_s1.00", "HCN_s1.25"):
            continue
        matrices = [block_matrix(raw, den, basis, n) for raw in d["blocks"]["class"]]
        expected = json.loads((ROOT / (case + ".transport_diagnostic.json")).read_text())["orders"]
        D = np.diag(np.diag(H))
        for name in ("class_lex", "class_magnitude", "class_diagonal_response"):
            order = d["orders"][name]
            suffix = sum((matrices[i] for i in order), np.zeros_like(H))
            B = np.zeros_like(H)
            for i in order:
                suffix = suffix - matrices[i]
                B += matrices[i] @ suffix - suffix @ matrices[i]
            # Chronological order: sign opposite to written order, irrelevant to variance.
            K = 0.5j*B
            JH, quadrature_error = transport_quadrature(H, K, initial)
            JD, diagonal_quadrature_error = transport_quadrature(D, K, initial)
            residual = max(abs(JH-expected[name]["J_H"]), abs(JD-expected[name]["J_D"]))
            if residual > 1e-11:
                raise RuntimeError("Independent dense transport mismatch")
            diagnostics.append(dict(case=case, ordering=name, J_H=JH, J_D=JD,
                                    maximum_difference=residual,
                                    quadrature_error=quadrature_error,
                                    diagonal_quadrature_error=diagonal_quadrature_error))
    report = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((ROOT/"run.py").read_bytes()).hexdigest(),
        "independent_circuit_rows": len(replay),
        "maximum_state_distance": max(r["state_distance"] for r in replay),
        "maximum_metric_difference": max(r["metric_difference"] for r in replay),
        "maximum_exact_reference_distance": max(r["independent_exact_distance"] for r in replay),
        "dense_transport_orders_checked": len(diagnostics),
        "maximum_transport_objective_difference": max(r["maximum_difference"] for r in diagnostics),
        "replay": replay,
        "transport": diagnostics,
    }
    OUT.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("replay","transport")}, indent=2))


if __name__ == "__main__":
    main()
