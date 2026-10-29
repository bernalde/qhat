"""
time_evolution.py -- Pauli-string Hamiltonian time evolution for QHAT.

Loads a Hamiltonian via the analysis configuration system and performs:
  exact     -- U = exp(-i H t)
  trotter1  -- 1st-order Lie-Trotter product formula
  trotter2  -- 2nd-order Suzuki-Trotter product formula

Place at the qhat-main/ root. Defaults to analysis/config.py if no config is given.

Hamiltonians can be supplied in two ways (mutually exclusive):
  config file  -- analysis .py config (fermionic tensors via hamiltonian)
  --pauli-file -- JSON file with Pauli strings directly (spin models, custom H)

Initial state notes:
  --psi0 zero  can be uninformative when every product-formula block preserves
               the vacuum. This does not follow from conservation by H alone.
               Use --psi0 hf for chemistry and report the chosen initial state.
  --psi0 hf    fermionic systems only; cannot be used with --pauli-file.

Usage
-----
    python time_evolution.py analysis/config.py --psi0 hf --n-electrons 3 --fidelity
    python time_evolution.py analysis/config.py --psi0 zero --fidelity
    python time_evolution.py --pauli-file xxz_cs2cocl4_10.json --psi0 zero --fidelity
    python time_evolution.py analysis/config.py --psi0 file --psi0-file my_state.npy --fidelity
"""

import argparse
import os as _os
import sys
import warnings
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import expm

from qhat.analysis.state_metrics import normalized_state, state_error_metrics

_ANALYSIS_DIR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "analysis")


def _load_state_from_config(config_file: str):
    from qhat.analysis.configuration import load_configuration

    # load_configuration() calls argparse internally, so temporarily isolate
    # sys.argv to just the config path before calling it.
    # Resolve any relative filename in config_hamiltonian to absolute so that
    # get_physical_hamiltonian() can find the file regardless of cwd.
    config_file  = _os.path.abspath(config_file)
    config_dir   = _os.path.dirname(config_file)
    saved_argv   = sys.argv[:]
    try:
        sys.argv = [sys.argv[0], config_file]
        state, _, _ = load_configuration()
    finally:
        sys.argv = saved_argv
    if getattr(state.config_hamiltonian, "filename", None) is not None:
        fname = state.config_hamiltonian.filename
        if not _os.path.isabs(fname):
            state.config_hamiltonian.filename = _os.path.join(config_dir, fname)
    return state


def load_pauli_file(path: str) -> Tuple[Dict[tuple, float], int]:
    import json
    with open(path) as fh:
        data = json.load(fh)
    n_qubits   = int(data["n_qubits"])
    pauli_dict: Dict[tuple, float] = {}
    for term in data["terms"]:
        key   = tuple((int(q), op) for q, op in term["ops"])
        coeff = float(term["coeff"])
        pauli_dict[key] = pauli_dict.get(key, 0.0) + coeff
    return pauli_dict, n_qubits


_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1],  [1,  0]], dtype=complex),
    "Y": np.array([[0, -1j],[1j, 0]], dtype=complex),
    "Z": np.array([[1, 0],  [0, -1]], dtype=complex),
}


def _pauli_string_to_matrix(term: tuple, n_qubits: int) -> np.ndarray:
    ops = {qubit: op for qubit, op in term}
    mat = np.array([[1.0 + 0j]])
    for q in range(n_qubits):
        mat = np.kron(mat, _PAULI.get(ops.get(q, "I")))
    return mat


def build_hamiltonian_matrix(pauli_dict: Dict[tuple, float], n_qubits: int) -> np.ndarray:
    H = np.zeros((2 ** n_qubits,) * 2, dtype=complex)
    for term, coeff in pauli_dict.items():
        H += coeff * _pauli_string_to_matrix(term, n_qubits)
    return H


def build_term_matrices(pauli_dict: Dict[tuple, float], n_qubits: int) -> List[np.ndarray]:
    return [coeff * _pauli_string_to_matrix(term, n_qubits)
            for term, coeff in pauli_dict.items()]


def exact_evolution(H: np.ndarray, t: float) -> np.ndarray:
    return expm(-1j * t * H)


def _validate_trotter_inputs(term_matrices, t, n_steps):
    if isinstance(n_steps, (bool, np.bool_)) or not isinstance(n_steps, (int, np.integer)) or n_steps <= 0:
        raise ValueError("n_steps must be a positive integer.")
    if not np.isfinite(t):
        raise ValueError("Evolution time must be finite.")
    if not term_matrices:
        raise ValueError("At least one term matrix is required (use a zero matrix for H=0).")


def trotter_1st_order(term_matrices: List[np.ndarray], t: float, n_steps: int) -> np.ndarray:
    """Apply the listed terms chronologically, with each new factor on the left."""
    _validate_trotter_inputs(term_matrices, t, n_steps)
    dt   = t / n_steps
    dim  = term_matrices[0].shape[0]
    step = np.eye(dim, dtype=complex)
    for Hk in term_matrices:
        step = expm(-1j * dt * Hk) @ step
    U = np.eye(dim, dtype=complex)
    for _ in range(n_steps):
        U = step @ U
    return U


def trotter_2nd_order(term_matrices: List[np.ndarray], t: float, n_steps: int) -> np.ndarray:
    """Apply a forward half sweep followed by the reverse half sweep."""
    _validate_trotter_inputs(term_matrices, t, n_steps)
    dt  = t / n_steps
    dim = term_matrices[0].shape[0]
    fwd = np.eye(dim, dtype=complex)
    for Hk in term_matrices:
        fwd = expm(-1j * (dt / 2) * Hk) @ fwd
    bwd = np.eye(dim, dtype=complex)
    for Hk in reversed(term_matrices):
        bwd = expm(-1j * (dt / 2) * Hk) @ bwd
    step = bwd @ fwd
    U = np.eye(dim, dtype=complex)
    for _ in range(n_steps):
        U = step @ U
    return U


def state_fidelity(U_exact: np.ndarray, U_approx: np.ndarray, psi0: np.ndarray) -> float:
    """Normalized squared-overlap fidelity; not the manuscript's 1-|overlap|.

    This remains a float-returning compatibility helper. For comparisons use
    state_error_metrics directly to retain stable small errors and norm checks.
    """
    psi0 = normalized_state(psi0, dimension=U_exact.shape[1])
    metrics = state_error_metrics(U_exact @ psi0, U_approx @ psi0)
    if max(abs(metrics["exact_norm"] - 1), abs(metrics["approx_norm"] - 1)) > 1e-8:
        warnings.warn("Evolution changes the state norm; fidelity alone cannot certify valid evolution.",
                      RuntimeWarning, stacklevel=2)
    return metrics["fidelity"]


_SEP = "=" * 65


def _print_unitary_summary(label, U, U_ref=None, psi0=None):
    print(f"\n{_SEP}\n  {label}\n{_SEP}")
    print(f"  Shape           : {U.shape}")
    print(f"  Unitarity check : ||U+U - I||_F = "
          f"{np.linalg.norm(U.conj().T @ U - np.eye(U.shape[0]), 'fro'):.3e}")
    if U_ref is not None and psi0 is not None:
        print(f"  State fidelity  : {state_fidelity(U_ref, U, psi0):.10f}")


def _print_fidelity_table(results, psi0, n_steps_1st, n_steps_2nd):
    U_exact = results.get("exact")
    if U_exact is None:
        print("\n[fidelity] Exact evolution not available.")
        return
    rows = []
    if "trotter1" in results:
        rows.append((f"1st-order Trotter (n={n_steps_1st})", results["trotter1"]))
    if "trotter2" in results:
        rows.append((f"2nd-order Trotter (n={n_steps_2nd})", results["trotter2"]))
    if not rows:
        print("\n[fidelity] No approximate methods to compare.")
        return
    print(f"\n{_SEP}\n  Fidelity summary  (state fidelity vs exact)\n{_SEP}")
    print("  F = squared overlap; E = 1 - |overlap| (manuscript metric).")
    print(f"  {'Method':<34} {'State F':>12} {'E':>12} {'1-F':>12} {'norm defect':>12}")
    psi0 = normalized_state(psi0, dimension=U_exact.shape[1])
    for label, U_approx in rows:
        metrics = state_error_metrics(U_exact @ psi0, U_approx @ psi0)
        norm_defect = max(abs(metrics['exact_norm'] - 1), abs(metrics['approx_norm'] - 1))
        print(f"  {label:<34} {metrics['fidelity']:>12.10f} "
              f"{metrics['one_minus_overlap']:>12.4e} {metrics['infidelity']:>12.4e} "
              f"{norm_defect:>12.4e}")
        if norm_defect > 1e-8:
            print("  [warning] Norm drift: overlap alone is not a valid propagation check.")
    print()


def run(pauli_dict, n_qubits, t, methods, n_steps_1st, n_steps_2nd,
        psi0, compute_fidelity, print_summaries):

    dim = 2 ** n_qubits
    if psi0 is None:
        psi0 = np.zeros(dim, dtype=complex)
        psi0[0] = 1.0
    psi0 = normalized_state(psi0, dimension=dim)
    if not np.isfinite(t):
        raise ValueError("Evolution time must be finite.")

    run_methods = set(methods)
    if compute_fidelity and "exact" not in run_methods:
        print("[info] --fidelity requested: adding exact evolution.")
        run_methods.add("exact")

    print(f"\nBuilding Hamiltonian for {n_qubits} qubits "
          f"({dim}x{dim} matrix, {len(pauli_dict)} Pauli terms) ...")
    H             = build_hamiltonian_matrix(pauli_dict, n_qubits)
    term_matrices = build_term_matrices(pauli_dict, n_qubits) or [np.zeros_like(H)]
    if not np.all(np.isfinite(H)) or not np.allclose(H, H.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError("The Hamiltonian must be finite and Hermitian.")
    print(f"Hamiltonian spectral norm : {np.linalg.norm(H, ord=2):.6f}")
    print(f"Evolution time t          : {t}")
    print(f"Methods requested         : {sorted(run_methods)}")

    results = {}

    if "exact" in run_methods:
        print("\nComputing exact evolution ...")
        results["exact"] = exact_evolution(H, t)
        if print_summaries:
            _print_unitary_summary("Exact evolution", results["exact"])

    if "trotter1" in run_methods:
        print(f"\nComputing 1st-order Trotter ({n_steps_1st} steps) ...")
        results["trotter1"] = trotter_1st_order(term_matrices, t, n_steps_1st)
        if print_summaries:
            _print_unitary_summary(f"1st-order Trotter (n={n_steps_1st})",
                                   results["trotter1"], results.get("exact"),
                                   psi0 if compute_fidelity else None)

    if "trotter2" in run_methods:
        print(f"\nComputing 2nd-order Trotter ({n_steps_2nd} steps) ...")
        results["trotter2"] = trotter_2nd_order(term_matrices, t, n_steps_2nd)
        if print_summaries:
            _print_unitary_summary(f"2nd-order Trotter (n={n_steps_2nd})",
                                   results["trotter2"], results.get("exact"),
                                   psi0 if compute_fidelity else None)

    if compute_fidelity:
        _print_fidelity_table(results, psi0, n_steps_1st, n_steps_2nd)

    return results


def make_zero_state(n_qubits: int) -> np.ndarray:
    psi0 = np.zeros(2 ** n_qubits, dtype=complex)
    psi0[0] = 1.0
    return psi0


def make_hartree_fock_state(n_qubits: int, n_electrons: int, encoding="JW") -> np.ndarray:
    """Encode a lowest-orbital occupation determinant, with qubit zero leftmost.

    This is a reference determinant, not a self-consistent HF calculation.
    BK bits are transported with the same binary code used by OpenFermion.
    Custom encodings require a supplied state instead of an assumed bit map.
    """
    if not isinstance(n_electrons, (int, np.integer)) or not (0 <= n_electrons <= n_qubits):
        raise ValueError(f"n_electrons={n_electrons} must be an integer in [0, {n_qubits}].")
    occupations = np.zeros(n_qubits, dtype=int)
    occupations[:n_electrons] = 1
    if encoding in ("JW", "jordan-wigner"):
        bits = occupations
    elif encoding in ("BK", "bravyi-kitaev"):
        from openfermion import bravyi_kitaev_code
        bits = np.asarray(bravyi_kitaev_code(n_qubits).encoder @ occupations).ravel() % 2
    else:
        raise ValueError(f"No reference-state transform for encoding {encoding!r}; use --psi0 file.")
    index = sum(int(bit) << (n_qubits - 1 - q) for q, bit in enumerate(bits))
    psi0 = np.zeros(2 ** n_qubits, dtype=complex)
    psi0[index] = 1.0
    return psi0


def make_file_state(path: str, n_qubits: int, key: str = "psi0") -> np.ndarray:
    dim = 2 ** n_qubits
    if path.endswith(".npz"):
        with np.load(path) as archive:
            if key not in archive:
                raise KeyError(f'Key "{key}" not found in "{path}". '
                               f'Available: {list(archive.keys())}. Use --psi0-key.')
            vec = archive[key].astype(complex)
    else:
        vec = np.load(path).astype(complex)
    return normalized_state(vec, dimension=dim)


def build_initial_state(args, n_qubits: int, using_pauli_file: bool, encoding="JW") -> Tuple[np.ndarray, str]:
    choice = args.psi0

    if choice == "hf" and using_pauli_file:
        raise ValueError("--psi0 hf is for fermionic systems only. "
                         "Use --psi0 file for spin models.")

    if choice == "zero":
        return make_zero_state(n_qubits), "|0...0> (computational basis)"

    if choice == "hf":
        if args.n_electrons is None:
            raise ValueError("--n-electrons is required when --psi0=hf.")
        psi0 = make_hartree_fock_state(n_qubits, args.n_electrons, encoding=encoding)
        idx  = int(np.argmax(np.abs(psi0)))
        return psi0, (f"Hartree-Fock reference ({args.n_electrons} electrons, "
                      f"orbitals 0-{args.n_electrons-1} occupied, {encoding} encoding, basis index {idx})")

    if choice == "file":
        if args.psi0_file is None:
            raise ValueError("--psi0-file is required when --psi0=file.")
        key   = getattr(args, "psi0_key", "psi0")
        psi0  = make_file_state(args.psi0_file, n_qubits, key)
        label = f'state loaded from "{args.psi0_file}"'
        if args.psi0_file.endswith(".npz"):
            label += f' (key="{key}")'
        return psi0, label

    raise ValueError(f'Unknown --psi0 choice: "{choice}".')


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Time-evolution of a QHAT Pauli-string Hamiltonian.\n\n"
            "--psi0 choices:\n"
            "  zero   |0...0> (default) -- can be an uninformative vacuum test\n"
            "  hf     Hartree-Fock reference (requires --n-electrons; fermionic only)\n"
            "  file   arbitrary state from .npy/.npz (requires --psi0-file)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("configuration_file", nargs="?", default=None,
                        help="Analysis config file (default: analysis/config.py)")
    parser.add_argument("--method", choices=["exact", "trotter1", "trotter2", "all"],
                        default="all", help="Evolution method (default: all)")
    parser.add_argument("--time", type=float, default=1.0, metavar="T",
                        help="Evolution parameter t/hbar in inverse Hamiltonian-energy units (default: 1.0)")
    parser.add_argument("--steps-1st", type=int, default=10, metavar="N", dest="steps_1st",
                        help="Trotter steps for 1st-order (default: 10)")
    parser.add_argument("--steps-2nd", type=int, default=5, metavar="N", dest="steps_2nd",
                        help="Trotter steps for 2nd-order (default: 5)")
    parser.add_argument("--psi0", choices=["zero", "hf", "file"],
                        default="zero", metavar="STATE", help="Initial state (default: zero)")
    parser.add_argument("--n-electrons", type=int, default=None, metavar="N",
                        dest="n_electrons", help="Electrons for HF state (--psi0=hf)")
    parser.add_argument("--psi0-file", type=str, default=None, metavar="PATH",
                        dest="psi0_file", help=".npy/.npz file for arbitrary initial state")
    parser.add_argument("--psi0-key", type=str, default="psi0", metavar="KEY",
                        dest="psi0_key", help='Array key in .npz archive (default: "psi0")')
    parser.add_argument("--fidelity", action="store_true",
                        help="Compute state fidelity vs exact (forces exact to run)")
    parser.add_argument("--no-print-unitaries", action="store_false", dest="print_summaries",
                        help="Suppress per-method unitary summary blocks")
    parser.add_argument("--pauli-file", type=str, default=None, metavar="PATH",
                        dest="pauli_file",
                        help="JSON file with Pauli strings (skips analysis config loading)")
    return parser.parse_args()


if __name__ == "__main__":

    from qhat.common.logging_utils import configure_logging
    configure_logging(level="warning")
    args    = _parse_args()
    methods = ["exact", "trotter1", "trotter2"] if args.method == "all" else [args.method]

    if args.pauli_file and args.configuration_file:
        raise SystemExit("error: --pauli-file and a config file are mutually exclusive.")

    if args.pauli_file:
        print(f'\nLoading Hamiltonian from Pauli file "{args.pauli_file}" ...')
        pauli_dict, n_qubits = load_pauli_file(args.pauli_file)
        using_pauli_file = True
        encoding = "JW"  # HF construction is rejected for bare Pauli input.
    else:
        from qhat.analysis.hamiltonian import get_physical_hamiltonian
        from openfermion import InteractionOperator
        config_file = args.configuration_file or _os.path.join(_ANALYSIS_DIR, "config.py")
        print(f'\nLoading Hamiltonian from config file "{config_file}" ...')
        state                = _load_state_from_config(config_file)
        physical_hamiltonian = get_physical_hamiltonian(state.config_hamiltonian)
        n_qubits         = physical_hamiltonian.num_qubits()
        pauli_dict       = physical_hamiltonian.get_all_pauli_strings(return_as="tuples")
        using_pauli_file = not isinstance(physical_hamiltonian.get_core_operator(), InteractionOperator)
        encoding = state.config_hamiltonian.fermion_to_qubit_transform

    print(f"  -> {n_qubits} qubits, {len(pauli_dict)} Pauli terms loaded.")

    psi0, psi0_label = build_initial_state(args, n_qubits, using_pauli_file, encoding=encoding)
    print(f"  -> Initial state : {psi0_label}")

    results = run(
        pauli_dict       = pauli_dict,
        n_qubits         = n_qubits,
        t                = args.time,
        methods          = methods,
        n_steps_1st      = args.steps_1st,
        n_steps_2nd      = args.steps_2nd,
        psi0             = psi0,
        compute_fidelity = args.fidelity,
        print_summaries  = args.print_summaries,
    )
