"""Small exact-matrix checks for the standalone comparison entrypoint."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest
from openfermion import FermionOperator, bravyi_kitaev, get_sparse_operator, jordan_wigner
from scipy.linalg import expm

from qhat import time_evolution as te
from qhat.analysis.state_metrics import state_error_metrics


X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1., -1.])


def test_fidelity_normalizes_its_public_input():
    assert te.state_fidelity(np.eye(2), np.eye(2), np.array([2., 0.])) == 1.


def test_fidelity_reports_nonunitary_norm_drift():
    with pytest.warns(RuntimeWarning, match='changes the state norm'):
        assert te.state_fidelity(np.eye(2), 2 * np.eye(2), np.array([1., 0.])) == 1.


@pytest.mark.parametrize('steps', [0, -1, 1.5, True])
@pytest.mark.parametrize('method', [te.trotter_1st_order, te.trotter_2nd_order])
def test_invalid_step_counts_are_rejected(steps, method):
    with pytest.raises(ValueError, match='positive integer'):
        method([X], 0.4, steps)


@pytest.mark.parametrize('method', [te.trotter_1st_order, te.trotter_2nd_order])
def test_single_term_evolution_is_exact_at_every_step_count(method):
    for steps in (1, 3, 7):
        np.testing.assert_allclose(method([0.7 * X], 0.4, steps), expm(-1j * 0.28 * X), atol=1e-14)


def test_chronological_order_is_explicit_and_noncommuting_reordering_matters():
    time = 0.3
    forward = te.trotter_1st_order([X, Z], time, 1)
    np.testing.assert_allclose(forward, expm(-1j * time * Z) @ expm(-1j * time * X), atol=1e-14)
    reverse = te.trotter_1st_order([Z, X], time, 1)
    assert np.linalg.norm(forward - reverse) > 0.1


def test_bk_reference_state_is_transported_with_the_hamiltonian():
    # The two-electron determinant has JW bits 1100 and BK bits 1000.
    jw = te.make_hartree_fock_state(4, 2, encoding='JW')
    bk = te.make_hartree_fock_state(4, 2, encoding='BK')
    assert np.argmax(abs(jw)) == 12
    assert np.argmax(abs(bk)) == 8
    f = FermionOperator('2^ 1', 0.7) + FermionOperator('1^ 2', 0.7)
    number = sum((FermionOperator(f'{q}^ {q}', 1.) for q in range(4)), FermionOperator())
    for transform, state in ((jordan_wigner, jw), (bravyi_kitaev, bk)):
        n = get_sparse_operator(transform(number), n_qubits=4).toarray()
        np.testing.assert_allclose(n @ state, 2 * state)
    hjw = get_sparse_operator(jordan_wigner(f), n_qubits=4).toarray()
    hbk = get_sparse_operator(bravyi_kitaev(f), n_qubits=4).toarray()
    ujw, ubk = expm(-1j * 0.3 * hjw), expm(-1j * 0.3 * hbk)
    assert state_error_metrics(jw, ujw @ jw)['one_minus_overlap'] == pytest.approx(
        state_error_metrics(bk, ubk @ bk)['one_minus_overlap'])
    # Reusing JW occupation bits under BK represents a different physical state.
    assert abs(state_error_metrics(jw, ubk @ jw)['one_minus_overlap'] -
               state_error_metrics(bk, ubk @ bk)['one_minus_overlap']) > 1e-3


def test_custom_encoding_requires_explicit_initial_state():
    with pytest.raises(ValueError, match='use --psi0 file'):
        te.make_hartree_fock_state(4, 2, encoding=lambda x: x)


def test_config_loader_unpacks_state_and_resolves_relative_file(tmp_path, monkeypatch):
    config = tmp_path / 'model.py'
    config.write_text('hamiltonian.load_second_quantization("molecule.npz", fermion_to_qubit_transform="BK")\n')
    monkeypatch.setattr(sys, 'argv', ['time_evolution.py', '--irrelevant-cli-argument'])
    original_argv = sys.argv[:]
    state = te._load_state_from_config(str(config))
    assert state.config_hamiltonian.filename == str(tmp_path / 'molecule.npz')
    assert state.config_hamiltonian.fermion_to_qubit_transform == 'BK'
    assert sys.argv == original_argv


@pytest.mark.parametrize('vector', [[0., 0.], [np.nan, 0.], [[1.], [0.]]])
def test_state_files_with_invalid_vectors_are_rejected(tmp_path, vector):
    path = tmp_path / 'psi.npy'
    np.save(path, vector)
    with pytest.raises(ValueError):
        te.make_file_state(str(path), 1)


def test_zero_hamiltonian_and_fidelity_output_are_well_defined(capsys):
    result = te.run({}, 1, 0.5, ['exact', 'trotter1', 'trotter2'], 3, 2,
                    np.array([1., 0.]), True, False)
    for matrix in result.values():
        np.testing.assert_allclose(matrix, np.eye(2))
    assert 'manuscript metric' in capsys.readouterr().out


@pytest.mark.parametrize('n_qubits', [3, 5, 6])
def test_bk_reference_matches_every_mode_occupation_for_non_power_of_two_sizes(n_qubits):
    # Testing every mapped occupation operator establishes the determinant's
    # encoding, rather than relying on a special small binary-code example.
    for n_electrons in sorted({0, 1, n_qubits // 2, n_qubits - 1, n_qubits}):
        state = te.make_hartree_fock_state(n_qubits, n_electrons, encoding='BK')
        assert np.linalg.norm(state) == 1
        for mode in range(n_qubits):
            occupation = FermionOperator(((mode, 1), (mode, 0)), 1.)
            encoded = bravyi_kitaev(occupation, n_qubits=n_qubits)
            matrix = get_sparse_operator(encoded, n_qubits=n_qubits)
            expected = state if mode < n_electrons else np.zeros_like(state)
            np.testing.assert_allclose(matrix @ state, expected, atol=1e-14)
