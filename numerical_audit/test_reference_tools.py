"""Independent dense oracles for the frozen audit's reduced-space/BCH machinery.

Matrices below are constructed from literal Kronecker products, independently
of reference_tools and the benchmark's compiled Pauli-action kernels.
"""

import numpy as np
import pytest
from scipy import sparse
from scipy.integrate import quad_vec
from scipy.linalg import expm

from qhat.analysis.state_metrics import state_error_metrics
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import dynamic_prediction, physical_metrics


PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.diag([1., -1.]),
}


def dense_pauli(word):
    matrix = np.array([[1.]], dtype=complex)
    for char in word:
        matrix = np.kron(matrix, PAULI[char])
    return matrix


def compiled_terms(words):
    """Literal bit convention: qubit zero is the most significant bit."""
    terms = []
    for coefficient, word in words:
        x = sum((1 << (len(word) - 1 - q)) for q, char in enumerate(word) if char in 'XY')
        z = sum((1 << (len(word) - 1 - q)) for q, char in enumerate(word) if char in 'ZY')
        terms.append((coefficient, np.uint64(x), np.uint64(z), word.count('Y')))
    return terms


def dense_data(words):
    terms = compiled_terms(words)
    matrices = [coefficient * dense_pauli(word) for coefficient, word in words]
    hamiltonian = sum(matrices, np.zeros_like(matrices[0]))
    defect = sum((left @ right - right @ left
                  for i, left in enumerate(matrices) for right in matrices[i + 1:]),
                 np.zeros_like(hamiltonian)) / 2
    return terms, matrices, hamiltonian, defect


def dense_step(matrices, dt):
    result = np.eye(matrices[0].shape[0], dtype=complex)
    for matrix in matrices:
        result = expm(-1j * dt * matrix) @ result
    return result


def normalized(vector):
    vector = np.asarray(vector, dtype=complex)
    return vector / np.linalg.norm(vector)


# H conserves number, but this interleaved individual-Pauli schedule can leave
# the one-particle sector. Odd-Y summands also test imaginary hopping phases.
NUMBER_CONSERVING_WORDS = [
    (.4, 'XXI'), (-.3, 'IXX'), (.2, 'XYI'), (.4, 'YYI'),
    (-.3, 'IYY'), (-.2, 'YXI'), (.1, 'ZII'), (-.17, 'IZI'), (.06, 'IIZ'),
]


def test_affine_restriction_intertwines_each_odd_y_term_on_complex_states():
    words = [(.4, 'XIIX'), (-.7, 'ZYII'), (.3, 'YZZX'), (.2, 'IXYZ'), (.11, 'IIII')]
    terms, matrices, dense_h, _ = dense_data(words)
    states, reduced, basis = ref.affine_space(terms, hf_index=10)
    assert len(states) == 2 ** len(basis) < 16
    assert len(set(states)) == len(states)
    assert states[0] == 10
    embedding = np.eye(16, dtype=complex)[:, states]
    reduced_h, parity = ref.assemble_hamiltonian(reduced, len(states))
    np.testing.assert_allclose(dense_h @ embedding, embedding @ reduced_h.toarray(), atol=1e-14)
    rng = np.random.default_rng(31)
    psi = normalized(rng.normal(size=len(states)) + 1j * rng.normal(size=len(states)))
    for original, term in zip(matrices, reduced):
        reduced_action = term[0] * ref.apply_p(psi, term, parity)
        np.testing.assert_allclose(original @ (embedding @ psi), embedding @ reduced_action, atol=1e-14)


def test_affine_closure_keeps_trotter_leakage_outside_a_conserved_number_sector():
    terms, matrices, dense_h, _ = dense_data(NUMBER_CONSERVING_WORDS)
    states, reduced, _ = ref.affine_space(terms, hf_index=4)
    assert set(states) == {1, 2, 4, 7}  # three one-electron states AND |111>
    assert 7 in states
    number = np.diag([int(i).bit_count() for i in range(8)])
    np.testing.assert_allclose(number @ dense_h, dense_h @ number, atol=1e-14)
    psi_full = np.zeros(8, dtype=complex)
    psi_full[[1, 2, 4]] = normalized([1, 1j, 2])
    psi_small = psi_full[states]
    _, parity = ref.assemble_hamiltonian(reduced, len(states))
    # Exponentiate each reduced Pauli exactly using P^2=I; this explicitly
    # validates the affine embedding without the simulator's evolution kernel.
    approximate_small = psi_small.copy()
    for term in reduced:
        coefficient = term[0]
        approximate_small = (np.cos(.4 * coefficient) * approximate_small -
                             1j * np.sin(.4 * coefficient) * ref.apply_p(approximate_small, term, parity))
    approximate_full = dense_step(matrices, .4) @ psi_full
    embedded = np.zeros(8, dtype=complex)
    embedded[states] = approximate_small
    np.testing.assert_allclose(embedded, approximate_full, atol=2e-14)
    assert abs(approximate_full[7]) > .01
    exact_full = expm(-.4j * dense_h) @ psi_full
    assert abs(exact_full[7]) < 1e-14
    mask = np.array([int(i).bit_count() == 1 for i in range(8)])
    metrics = physical_metrics(exact_full, approximate_full, mask, mask)
    assert metrics['outside_sector_infidelity'] == pytest.approx(abs(approximate_full[7])**2)
    assert metrics['particle_number_leakage'] == pytest.approx(abs(approximate_full[7])**2)
    assert metrics['infidelity'] > metrics['conditional_sector_infidelity']


def test_leading_defect_has_the_chronological_commutator_sign():
    words = [(.4, 'X'), (-.2, 'Y'), (.7, 'Z'), (.13, 'I')]
    terms, matrices, dense_h, dense_defect = dense_data(words)
    hamiltonian, parity = ref.assemble_hamiltonian(terms, 2)
    psi = normalized([1 + 2j, .7 - .3j])
    expected = dense_defect @ psi
    actual = ref.leading_defect(terms, psi, hamiltonian, parity)
    np.testing.assert_allclose(actual, expected, atol=2e-14)
    # The imaginary expectation is a global-phase component which must be
    # projected out for fidelity, not counted as physical ray error.
    assert abs(np.vdot(psi, expected)) > .1
    assert np.linalg.norm(expected) > .2
    np.testing.assert_allclose(dense_defect.conj().T, -dense_defect, atol=1e-14)
    errors = []
    for dt in (.01, .005, .0025):
        observed = ((dense_step(matrices, dt) - expm(-1j * dt * dense_h)) @ psi) / dt**2
        errors.append(np.linalg.norm(observed - actual))
    assert 1.9 < errors[0] / errors[1] < 2.1
    assert 1.9 < errors[1] / errors[2] < 2.1
    reversed_defect = ref.leading_defect(list(reversed(terms)), psi, hamiltonian, parity)
    np.testing.assert_allclose(reversed_defect, -actual, atol=2e-14)


def test_defect_columns_match_dense_commutators_including_outside_sector_rows():
    terms, _, dense_h, dense_defect = dense_data(NUMBER_CONSERVING_WORDS)
    hamiltonian, parity = ref.assemble_hamiltonian(terms, 8)
    sector = np.array([1, 2, 4])
    columns = ref.defect_columns(terms, sector, hamiltonian, parity)
    np.testing.assert_allclose(columns, dense_defect[:, sector], atol=2e-14)
    assert np.linalg.norm(columns[7]) > .1
    coefficients = normalized([1, 1j, 2])
    initial = np.zeros(8, dtype=complex)
    initial[sector] = coefficients
    np.testing.assert_allclose(columns @ coefficients,
                               ref.leading_defect(terms, initial, hamiltonian, parity), atol=2e-14)
    np.testing.assert_allclose(hamiltonian.toarray(), dense_h, atol=1e-14)


@pytest.mark.parametrize('restricted', [False, True])
def test_dynamic_prediction_matches_independent_dense_dyson_integral_and_convergence(restricted):
    if restricted:
        words = NUMBER_CONSERVING_WORDS
        sector = np.array([1, 2, 4])
        initial = np.zeros(8, dtype=complex)
        initial[sector] = normalized([1, 1j, 2])
    else:
        words = [(.4, 'X'), (-.2, 'Y'), (.7, 'Z'), (.13, 'I')]
        sector = np.array([0, 1])
        initial = normalized([1 + 2j, .7 - .3j])
    terms, matrices, dense_h, dense_defect = dense_data(words)
    hamiltonian, parity = ref.assemble_hamiltonian(terms, len(initial))
    total_time = .73
    exact = expm(-1j * total_time * dense_h) @ initial
    eta, quadrature_error = quad_vec(
        lambda s: expm(-1j * (total_time - s) * dense_h) @ dense_defect @
                  expm(-1j * s * dense_h) @ initial,
        0., total_time, epsabs=1e-12, epsrel=1e-12,
    )
    assert quadrature_error < 1e-10
    assert abs(np.vdot(exact, eta)) > .01
    projected = eta - exact * np.vdot(exact, eta)
    expected = float(np.vdot(projected, projected).real)
    prediction = dynamic_prediction(hamiltonian, initial, terms, sector, parity, total_time)
    assert prediction['direct_column_residual'] < 2e-14
    assert prediction['integrated_projected_defect_norm2'] == pytest.approx(expected, rel=2e-12, abs=2e-14)
    # A norm-only coefficient is wrong when the integrated defect has a
    # nonzero phase expectation. This checks the variance/projection step.
    assert float(np.vdot(eta, eta).real) - expected > 1e-4
    if restricted:
        assert abs(projected[7]) > .02  # the prediction includes leakage
    relative_errors = []
    for steps in (32, 64, 128, 256):
        approximate = np.linalg.matrix_power(dense_step(matrices, total_time / steps), steps) @ initial
        actual = state_error_metrics(exact, approximate)['infidelity']
        predicted = (total_time / steps)**2 * expected
        relative_errors.append(abs(predicted / actual - 1.0))
    assert all(later < earlier for earlier, later in zip(relative_errors, relative_errors[1:]))
    assert relative_errors[-1] < .01
    assert relative_errors[0] / relative_errors[-1] > 6
