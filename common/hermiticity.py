"""
Checks that Hamiltonian data are Hermitian and that derived physical quantities are real.

Hermitian Hamiltonians may be specified with complex integrals (e.g. relativistic chemistry), but
their Pauli coefficients, energies, evolution times, and phases are real. Library code (notably
OpenFermion) can nonetheless return those as complex numbers with round-off imaginary parts, which
then silently propagate. These helpers reject genuinely complex values and return real ones as
Python floats.
"""

import numbers

import numpy as np

DEFAULT_RTOL = 1e-12


def require_real_scalar(name, x, rtol=DEFAULT_RTOL):
    """Return `x` as a float, raising ValueError if it is not a finite real number.

    `x` may be complex with an imaginary part at most `rtol * |x|`, which is discarded.
    """
    if isinstance(x, np.ndarray) and x.ndim == 0:
        x = x[()]
    if not isinstance(x, numbers.Number) or isinstance(x, bool):
        raise TypeError(f"{name} must be a real number, got {type(x).__name__} {x!r}.")
    if not np.isfinite(x):
        raise ValueError(f"{name} must be finite, got {x}.")
    imag = abs(np.imag(x))
    if imag > rtol * abs(x):
        raise ValueError(f"{name} must be real, got {x} (|imag| / |value| = {imag / abs(x):.3g}, "
                         f"tolerance {rtol:.3g}).")
    return float(np.real(x))


def _require_finite(name, array):
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values.")


def _check_conjugate_symmetric(name, array, conjugate_transpose, description, rtol):
    scale = np.abs(array).max(initial=0.0)
    if scale == 0.0:
        return
    defect = np.abs(array - conjugate_transpose)
    worst = np.unravel_index(np.argmax(defect), defect.shape)
    ratio = defect[worst] / scale
    if ratio > rtol:
        raise ValueError(
            f"{name} is not Hermitian: {description} fails at index {tuple(int(i) for i in worst)} "
            f"(max |defect| / max |element| = {ratio:.3g}, tolerance {rtol:.3g}).")


def check_hermitian_one_body(name, h1, rtol=DEFAULT_RTOL):
    """Raise ValueError unless sum_pq h1[p,q] a+_p a_q is Hermitian, i.e. h1 = h1^dagger."""
    h1 = np.asarray(h1)
    _require_finite(name, h1)
    _check_conjugate_symmetric(name, h1, h1.conj().T, "h[p,q] = conj(h[q,p])", rtol)


def check_hermitian_two_body(name, h2, rtol=DEFAULT_RTOL):
    """Raise ValueError unless sum_pqrs h2[p,q,r,s] a+_p a+_q a_r a_s is Hermitian.

    The operator depends only on the part of h2 antisymmetric under p <-> q and under r <-> s, so
    that part is tested; a Hermitian operator stored with an unsymmetrized tensor is accepted.
    """
    h2 = np.asarray(h2)
    _require_finite(name, h2)
    a = (h2 - h2.transpose(1, 0, 2, 3) - h2.transpose(0, 1, 3, 2) + h2.transpose(1, 0, 3, 2)) / 4
    _check_conjugate_symmetric(name, a, a.transpose(3, 2, 1, 0).conj(),
                               "antisymmetrized h[p,q,r,s] = conj(h[s,r,q,p])", rtol)


def max_relative_imaginary(array):
    """Return max |imag(array)| / max |array| (0 for a zero or real-typed array)."""
    array = np.asarray(array)
    scale = np.abs(array).max(initial=0.0)
    if scale == 0.0 or not np.iscomplexobj(array):
        return 0.0
    return float(np.abs(array.imag).max() / scale)


def real_pauli_coefficients(terms, name="Pauli coefficients", rtol=DEFAULT_RTOL):
    """Return a copy of `terms` (Pauli key -> coefficient) with float coefficients.

    Raises ValueError if sum |Im c| > rtol * sum |c| over non-identity terms. The left side bounds
    the norm of the anti-Hermitian part of the operator; the identity term is excluded from the
    scale because a large energy offset would otherwise mask imaginary parts of the other terms.
    Identity keys are `()` (sparse) or all-'I' strings (dense).
    """
    def is_identity(key):
        return len(key) == 0 or (isinstance(key, str) and set(key) == {'I'})

    imag_weight = 0.0
    scale = 0.0
    for key, coefficient in terms.items():
        if not np.isfinite(coefficient):
            raise ValueError(f"{name}: coefficient of {key!r} is not finite ({coefficient}).")
        imag_weight += abs(np.imag(coefficient))
        if not is_identity(key):
            scale += abs(coefficient)
    if scale == 0.0:
        scale = sum(abs(c) for c in terms.values())
    if imag_weight > rtol * scale:
        worst = max(terms, key=lambda k: abs(np.imag(terms[k])))
        raise ValueError(
            f"{name} are not real, so the Hamiltonian is not Hermitian: sum |imag| / sum |coef| = "
            f"{imag_weight / scale:.3g} (tolerance {rtol:.3g}); largest imaginary part is "
            f"{terms[worst]} for {worst!r}.")
    return {key: float(np.real(coefficient)) for key, coefficient in terms.items()}


def real_pauli_file_coefficients(terms):
    """Per-term check used by the Pauli-string file loaders: |Im c| <= 1e-8 |c| for each term.

    Returns a copy of `terms` with float coefficients.
    """
    result = {}
    for pauli, coef in terms.items():
        coef = complex(coef)
        if abs(coef.imag) > abs(coef) * 1e-8:
            imag_ratio_percent = abs(coef.imag) / abs(coef) * 100
            raise ValueError(
                f"Hamiltonian must be Hermitian (real coefficients). "
                f"Found coefficient {coef} where imaginary part is "
                f"{imag_ratio_percent:.4g}% of magnitude (max allowed: 1e-6%).")
        result[pauli] = coef.real
    return result
