"""Pure-state error metrics with explicit normalization and phase conventions.

For normalized vectors e and a, the manuscript uses E = 1 - |<e,a>|,
whereas state infidelity is 1 - |<e,a>|**2 = E * (2 - E).  The phase-sensitive
vector distance measures a different observable and can rank orderings
oppositely. Comparisons must hold the Hamiltonian, input state, physical time,
and computational budget fixed.

Normalization makes the ray metrics meaningful for non-unit input vectors;
it does not certify unitary propagation. Always inspect the returned norms
and relative_norm_drift alongside those metrics.
"""

import numpy as np


def _state_and_norm(state, *, dimension=None):
    vector = np.asarray(state, dtype=complex)
    if vector.ndim != 1 or vector.size == 0:
        raise ValueError("State must be a nonempty one-dimensional vector.")
    if dimension is not None and vector.size != dimension:
        raise ValueError(f"State dimension {vector.size} does not match {dimension}.")
    if not np.all(np.isfinite(vector)):
        raise ValueError("State contains non-finite amplitudes.")
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm == 0:
        raise ValueError("State must have a finite, nonzero norm.")
    return vector, norm


def normalized_state(state, *, dimension=None):
    """Validate and normalize a state without changing the caller's array."""
    vector, norm = _state_and_norm(state, dimension=dimension)
    return vector / norm


def state_error_metrics(exact_final, approx_final, *, norm_tolerance=1e-8):
    """Compare final pure states, retaining norm-drift diagnostics.

    ``one_minus_overlap`` is evaluated as half the squared distance after
    aligning the global phase. This avoids subtracting an overlap close to
    one. ``infidelity`` uses E*(2-E) for the same reason. ``fidelity`` itself
    may round to one even when the stable infidelity is positive.

    ``norm_preserved`` only checks relative norm drift between these two
    outputs. It is not a test of unitarity or of correct dynamics. Non-finite
    and zero states are rejected instead of being normalized silently.
    """
    if not np.isfinite(norm_tolerance) or norm_tolerance < 0:
        raise ValueError("norm_tolerance must be finite and nonnegative.")
    e, exact_norm = _state_and_norm(exact_final)
    a, approx_norm = _state_and_norm(approx_final, dimension=e.size)
    e, a = e / exact_norm, a / approx_norm
    overlap = np.vdot(e, a)
    magnitude = float(abs(overlap))
    aligned_a = a if magnitude == 0 else a * (overlap.conjugate() / magnitude)
    phase_aligned_distance = float(np.linalg.norm(e - aligned_a))
    one_minus_overlap = float(np.clip(0.5 * phase_aligned_distance**2, 0.0, 1.0))
    infidelity = float(one_minus_overlap * (2.0 - one_minus_overlap))
    norm_drift = float(abs(approx_norm / exact_norm - 1.0))
    return {
        "overlap_magnitude": float(np.clip(magnitude, 0.0, 1.0)),
        "fidelity": float(np.clip(magnitude**2, 0.0, 1.0)),
        "one_minus_overlap": one_minus_overlap,
        "infidelity": infidelity,
        "phase_aligned_distance": phase_aligned_distance,
        "normalized_phase_sensitive_distance": float(np.linalg.norm(e - a)),
        "exact_norm": exact_norm,
        "approx_norm": approx_norm,
        "relative_norm_drift": norm_drift,
        "norm_preserved": bool(norm_drift <= norm_tolerance),
    }
