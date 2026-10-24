"""Small dense oracles for the chronological product formulas used in benchmarks.

These compare implemented circuits with independent matrix exponentials, not with
another expansion of the implementation's cached repeat pattern. In particular,
unitarity alone cannot detect evolving for the wrong total time.
"""

import numpy as np
import pytest
from scipy.linalg import expm

from qhat.analysis.ordering import reorder_paulis
from qhat.common.trotter_flattened import Trotterization


_PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.diag([1, -1]).astype(complex),
}


def matrix(label):
    result = np.ones((1, 1), dtype=complex)
    for char in label:
        result = np.kron(result, _PAULIS[char])
    return result


def chronological_product(factors):
    result = np.eye(len(factors[0]), dtype=complex)
    for factor in factors:
        result = factor @ result
    return result


def product_formula(matrices, order, dt):
    """Lie, Strang, and Suzuki definitions on whole Hamiltonian blocks."""
    if order == 1:
        return chronological_product([expm(-1j * dt * h) for h in matrices])
    if order == 2:
        halves = [expm(-0.5j * dt * h) for h in matrices]
        return chronological_product(halves + halves[::-1])
    assert order == 4
    p = 1 / (4 - 4 ** (1 / 3))
    outer = product_formula(matrices, 2, p * dt)
    center = product_formula(matrices, 2, (1 - 4 * p) * dt)
    return outer @ outer @ center @ outer @ outer


@pytest.mark.parametrize("label", ["Y", "I"])
@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("steps", [1, 3, 7])
@pytest.mark.parametrize("contraction", ["structured", "incremental", "qualtran"])
def test_single_term_preserves_total_time(label, order, steps, contraction):
    """Collapsing r steps must still implement exp(-i T h P / hbar)."""
    coefficient, time, hbar = -0.61, 0.23, 1.7
    actual = Trotterization.from_method(
        [(label, coefficient)], order, time, steps, hbar=hbar,
        tensor_contraction_method=contraction,
    ).tensor_contract()
    expected = expm(-1j * time * coefficient * matrix(label) / hbar)
    np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=0)


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("steps", [1, 4])
@pytest.mark.parametrize("combine", [False, True])
def test_commuting_hamiltonian_is_exact_including_identity(order, steps, combine):
    terms = (("XX", 0.31), ("YY", -0.23), ("ZZ", 0.47), ("II", 0.19))
    time = 0.38
    hamiltonian = sum(coefficient * matrix(label) for label, coefficient in terms)
    actual = Trotterization.from_method(
        terms, order, time, steps, combine_terms=combine,
    ).tensor_contract()
    np.testing.assert_allclose(actual, expm(-1j * time * hamiltonian), atol=2e-13, rtol=0)


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("steps", [1, 3])
@pytest.mark.parametrize("combine", [False, True])
@pytest.mark.parametrize("contraction", ["structured", "incremental", "qualtran"])
def test_noncommuting_schedule_matches_chronological_formula(order, steps, combine, contraction):
    # The repeated label is intentionally separated by noncommuting terms: its
    # two coefficients must not silently be overwritten or globally combined.
    terms = (("XI", 0.4), ("YX", -0.2), ("IZ", 0.7), ("XI", 0.1))
    matrices = [coefficient * matrix(label) for label, coefficient in terms]
    time, hbar = 0.31, 1.4
    expected = np.linalg.matrix_power(product_formula(matrices, order, time / (steps * hbar)), steps)
    actual = Trotterization.from_method(
        terms, order, time, steps, hbar=hbar, combine_terms=combine,
        tensor_contraction_method=contraction,
    ).tensor_contract()
    np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)


def test_first_order_sequence_is_application_order_not_matrix_order():
    time = 0.43
    terms = (("X", 0.7), ("Z", -0.4))
    ux, uz = [expm(-1j * time * coefficient * matrix(label)) for label, coefficient in terms]
    actual = Trotterization.from_method(terms, 1, time, 1).tensor_contract()
    np.testing.assert_allclose(actual, uz @ ux, atol=2e-14, rtol=0)
    assert np.linalg.norm(actual - ux @ uz) > 0.1


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("steps", [1, 3])
def test_flattening_contiguous_commuting_blocks_preserves_formula(order, steps):
    groups = (
        (("XX", 0.4), ("YY", -0.7), ("ZZ", 0.2)),
        (("ZI", 0.8), ("IZ", -0.1)),
        (("XY", 0.3), ("YX", -0.6)),
    )
    blocks = [sum(coefficient * matrix(label) for label, coefficient in group) for group in groups]
    terms = tuple(term for group in groups for term in group)
    time = 0.29
    expected = np.linalg.matrix_power(product_formula(blocks, order, time / steps), steps)
    actual = Trotterization.from_method(terms, order, time, steps).tensor_contract()
    np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)


@pytest.mark.parametrize("ordering", [None, "magnitude", "lexicographical", "group_evolve_xyz", "group_evolve_greedy"])
def test_reordering_preserves_the_same_hamiltonian(ordering):
    terms = {"ZI": -0.9, "II": 0.13, "YY": -0.4, "XI": 0.4, "XX": 0.8, "IZ": 0.0}
    reordered = reorder_paulis(terms, ordering)
    assert reordered == terms
    np.testing.assert_array_equal(
        sum(coefficient * matrix(label) for label, coefficient in reordered.items()),
        sum(coefficient * matrix(label) for label, coefficient in terms.items()),
    )


def test_xyz_does_not_silently_drop_mixed_molecular_paulis():
    with pytest.raises(Exception, match="at most one non-identity pauli type"):
        reorder_paulis({"XZY": 0.7, "YYI": -0.2}, "group_evolve_xyz")


@pytest.mark.parametrize("label", ["I", "Y", "IX", "YI", "XIY"])
def test_pauli_tensor_and_circuit_preserve_support_and_global_phase(label):
    from qualtran.simulation.tensor import bloq_to_dense
    from qhat.common.pauli_string_evolution import PauliStringEvolution

    evolution = PauliStringEvolution(label, coefficient=-0.61, time=0.23, hbar=1.7)
    expected = expm(1j * 0.61 * 0.23 / 1.7 * matrix(label))
    # Qualtran0.4 (required by pyLIQTR) predates the full_flatten keyword.
    # Keep both direct contraction and the explicit circuit decomposition.
    for bloq in (evolution, evolution.decompose_bloq()):
        actual = bloq_to_dense(bloq)
        np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=0)


@pytest.mark.parametrize("label", ["I", "Y", "IX"])
def test_controlled_pauli_evolution_preserves_relative_phase(label):
    from scipy.linalg import block_diag
    from qhat.common.pauli_string_evolution import PauliStringEvolution

    evolution = PauliStringEvolution(label, coefficient=0.61, time=0.23, hbar=1.7)
    target = expm(-1j * 0.61 * 0.23 / 1.7 * matrix(label))
    expected = block_diag(np.eye(len(target)), target)
    actual = evolution.controlled().tensor_contract()
    np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=0)
