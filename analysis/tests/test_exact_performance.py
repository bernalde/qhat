"""Small correctness check; timing comparisons live in the opt-in benchmark."""
import pytest
from openfermion import QubitOperator

from qhat.analysis.trotter_coefficients_fast import trotter_error_estimator_fast


def test_auto_exact_matches_explicit_exact_for_weighted_terms():
    terms = [QubitOperator('X0', .7), QubitOperator('Y0', -.4),
             QubitOperator('Z0 Z1', 1.2), QubitOperator('X1', .3)]
    explicit = trotter_error_estimator_fast(terms, 60, mode='exact')
    automatic = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)
    assert automatic == pytest.approx(explicit, rel=1e-12, abs=1e-14)
    assert all(value > 0 for value in explicit)
