"""Independent small occupation-matrix checks; no QHAT imports."""

from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path
import random

import numpy as np

from pair_frobenius import pair_weight, signed_xz_commutator


def integer_matrix(n):
    rng = random.Random(7129 + n)
    dim = 1 << n
    H = [[0] * dim for _ in range(dim)]
    for u in range(dim):
        for v in range(u + 1):
            if u.bit_count() == v.bit_count():
                H[u][v] = H[v][u] = rng.randrange(-3, 4)
    return H


def coefficients(H, x):
    dim = len(H)
    return {
        z: c
        for z in range(dim)
        if (c := sum(
            Fraction((-1) ** ((z & v).bit_count() % 2) * H[v ^ x][v], dim)
            for v in range(dim)
        )) != 0
    }


def apply_weight(coeff, v):
    return sum(a * (-1) ** ((z & v).bit_count() % 2) for z, a in coeff.items())


def full_transport(H, K, T):
    energies, vectors = np.linalg.eigh(H)
    gaps = energies[:, None] - energies[None, :]
    kernel = T * np.exp(0.5j * T * gaps) * np.sinc(T * gaps / (2 * np.pi))
    return vectors @ (kernel * (vectors.conj().T @ K @ vectors)) @ vectors.conj().T


def main():
    exact_pair_checks = 0
    for n in (2, 3, 4):
        H = integer_matrix(n)
        classes = {x: coefficients(H, x) for x in range(1 << n)}
        for x, a in classes.items():
            for y, b in classes.items():
                mask, comm = signed_xz_commutator(x, a, y, b)
                for v in range(1 << n):
                    expected = H[(v ^ y) ^ x][v ^ y] * H[v ^ y][v] - (
                        H[(v ^ x) ^ y][v ^ x] * H[v ^ x][v]
                    )
                    actual = apply_weight(comm, v)
                    if actual != expected:
                        raise RuntimeError("Pauli commutator versus occupation matrices failed")
                for eta in range(n + 1):
                    expected = sum(
                        (
                            H[v ^ mask][v ^ y] * H[v ^ y][v]
                            - H[v ^ mask][v ^ x] * H[v ^ x][v]
                        ) ** 2
                        for v in range(1 << n) if v.bit_count() == eta
                    )
                    actual = pair_weight(n, eta, x, a, y, b)
                    if actual != expected:
                        raise RuntimeError("Krawtchouk weight versus exact occupation trace failed")
                    exact_pair_checks += 1

    n, eta, v, T = 4, 2, 3, 0.13
    H = np.asarray(integer_matrix(n), dtype=float)
    sector = np.array([u for u in range(1 << n) if u.bit_count() == eta])
    pos = list(sector).index(v)
    masks = [x for x in range(1 << n) if coefficients(H.astype(int).tolist(), x)]
    blocks = {
        x: np.array([[H[u, w] if u ^ w == x else 0.0 for w in sector] for u in sector])
        for x in masks
    }
    Hs = H[np.ix_(sector, sector)]
    D = np.diag(np.diag(Hs))
    gaps = np.diag(D) - D[pos, pos]
    rho = T**2 * np.sinc(T * gaps / (2 * np.pi))**2
    order = list(reversed(masks))
    B = sum(
        (blocks[x] @ blocks[y] - blocks[y] @ blocks[x])
        for i, x in enumerate(order) for y in order[i+1:]
    )
    K = -0.5j * B
    Qd = full_transport(D, K, T)
    j_direct = float(np.vdot(Qd[:, pos], Qd[:, pos]).real - abs(Qd[pos, pos])**2)
    j_formula = float(0.25 * np.dot(rho, B[:, pos]**2))
    maximum_residual = abs(j_direct - j_formula)
    for i in range(len(order) - 1):
        x, y = order[i:i+2]
        C = blocks[x] @ blocks[y] - blocks[y] @ blocks[x]
        Bnext = B - 2 * C
        predicted_delta = float(np.dot(rho, C[:, pos]**2 - B[:, pos] * C[:, pos]))
        actual_delta = float(0.25*np.dot(rho, Bnext[:, pos]**2)) - j_formula
        maximum_residual = max(maximum_residual, abs(predicted_delta - actual_delta))
    for s in range(1 << n):
        U = np.diag([(-1)**((s & int(u)).bit_count() % 2) for u in sector])
        transformed = {x: U @ blocks[x] @ U for x in masks}
        Bs = sum(
            transformed[x] @ transformed[y] - transformed[y] @ transformed[x]
            for i, x in enumerate(order) for y in order[i+1:]
        )
        maximum_residual = max(
            maximum_residual, abs(float(0.25*np.dot(rho, Bs[:, pos]**2))-j_formula)
        )
    Qh = full_transport(Hs, K, T)
    transport_error = float(np.linalg.norm(Qh-Qd, 2))
    transport_bound = float(T*T*np.linalg.norm(Hs-D, 2)*np.linalg.norm(K, 2))
    if maximum_residual > 1e-11 or transport_error > transport_bound + 1e-12:
        raise RuntimeError("Diagonal reference validation failed")
    report = {
        "exact_pair_sector_checks": exact_pair_checks,
        "dense_n": n,
        "dense_eta": eta,
        "dense_T": T,
        "dense_class_count": len(masks),
        "dense_objective": j_direct,
        "maximum_objective_update_gauge_residual": maximum_residual,
        "all_spin_orbital_sign_masks_checked": 1 << n,
        "exact_transport_difference_norm": transport_error,
        "conservative_transport_upper_bound": transport_bound,
        "scope": "small synthetic number-conserving real matrices; no molecular accuracy claim",
    }
    Path(__file__).with_name("prototype_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
