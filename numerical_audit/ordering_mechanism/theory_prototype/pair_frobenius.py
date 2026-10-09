"""Exact finite-sector trace prototype, independent of QHAT.

Coefficients use W(x,z)|v> = (-1)**popcount(z & v)|v ^ x>.
They are NOT conventional Hermitian tensor-Pauli coefficients when Y occurs.
All arithmetic is exact when the supplied coefficients are int/Fraction.
No Hilbert-space state enumeration occurs in the exported weight functions.
"""

from math import comb


def _binomial(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def fixed_weight_character(n, eta, k):
    """Sum (-1)**popcount(v&z) over |v|=eta, for |z|=k."""
    if not (0 <= eta <= n and 0 <= k <= n):
        raise ValueError("Need 0 <= eta,k <= n")
    return sum(
        (-1) ** j * _binomial(k, j) * _binomial(n - k, eta - j)
        for j in range(k + 1)
    )


def signed_xz_commutator(x, first, y, second):
    """Return (mask, coefficients) for [sum_z first[z]W(x,z), ...]."""
    output = {}
    for z, a in first.items():
        for w, b in second.items():
            phase = (-1) ** ((z & y).bit_count() % 2) - (-1) ** (
                (w & x).bit_count() % 2
            )
            if phase:
                key = z ^ w
                output[key] = output.get(key, 0) + phase * a * b
    return x ^ y, {z: c for z, c in output.items() if c != 0}


def sector_input_frobenius_sq(n, eta, coefficients):
    """Exact ||C P_eta||_F^2 for a common-flip-mask signed-XZ expansion.

    If C preserves eta, this is exactly its sector-restricted Frobenius norm.
    The flip mask does not enter this expression because all outputs for a
    fixed input coincide. Individual Pauli words need not preserve eta.
    """
    if not (0 <= eta <= n):
        raise ValueError("Need 0 <= eta <= n")
    if any(z < 0 or z >= 1 << n for z in coefficients):
        raise ValueError("phase mask out of range")
    chars = [fixed_weight_character(n, eta, k) for k in range(n + 1)]
    result = 0
    for z, a in coefficients.items():
        for w, b in coefficients.items():
            result += a.conjugate() * b * chars[(z ^ w).bit_count()]
    return result


def pair_weight(n, eta, x, first, y, second):
    """Squared sector Frobenius commutator weight for conserving inputs."""
    _, coefficients = signed_xz_commutator(x, first, y, second)
    return sector_input_frobenius_sq(n, eta, coefficients)
