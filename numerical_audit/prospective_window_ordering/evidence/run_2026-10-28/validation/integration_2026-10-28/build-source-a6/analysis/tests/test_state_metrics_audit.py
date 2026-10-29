"""Regression cases for physically distinct ordering-comparison observables."""

import numpy as np
import pytest
from scipy.linalg import expm

from qhat.analysis.error_analysis import (
    _compute_eigenvalue_errors, _compute_state_dependent_errors, _save_error_results,
)
from qhat.analysis.operators import OperatorRepresentation
from qhat.analysis.state_metrics import normalized_state, state_error_metrics


def test_phase_sensitive_distance_can_reverse_a_state_accuracy_ranking():
    exact = np.array([1.0, 0.0])
    pure_phase = np.exp(1j) * exact
    physical_rotation = np.array([np.cos(0.1), np.sin(0.1)])
    phase = state_error_metrics(exact, pure_phase)
    rotation = state_error_metrics(exact, physical_rotation)
    assert phase['normalized_phase_sensitive_distance'] > rotation['normalized_phase_sensitive_distance']
    assert phase['one_minus_overlap'] < rotation['one_minus_overlap']
    assert phase['infidelity'] < 1e-28
    assert rotation['infidelity'] == pytest.approx(np.sin(0.1)**2)


def test_tiny_errors_survive_rounded_fidelity():
    theta = 1e-10
    metrics = state_error_metrics([1.0, 0.0], [np.cos(theta), np.sin(theta)])
    assert metrics['fidelity'] == 1.0  # ordinary subtraction would report zero
    assert metrics['infidelity'] == pytest.approx(theta**2, rel=1e-12, abs=0)
    assert metrics['one_minus_overlap'] == pytest.approx(theta**2 / 2, rel=1e-12, abs=0)


def test_rescaling_cannot_produce_fidelity_larger_than_one_and_drift_is_visible():
    metrics = state_error_metrics([2.0, 0.0], [3.0j, 0.0])
    assert metrics['fidelity'] == 1.0
    assert metrics['one_minus_overlap'] == 0.0
    assert metrics['exact_norm'] == 2.0
    assert metrics['approx_norm'] == 3.0
    assert metrics['relative_norm_drift'] == 0.5
    assert not metrics['norm_preserved']


@pytest.mark.parametrize('state', [[0, 0], [np.nan, 0], [np.inf, 0], [[1], [0]], []])
def test_invalid_state_is_rejected(state):
    with pytest.raises(ValueError):
        normalized_state(state)
    with pytest.raises(ValueError):
        state_error_metrics(state, [1, 0])


def test_mismatched_dimensions_are_rejected():
    with pytest.raises(ValueError, match='dimension'):
        state_error_metrics([1, 0], [1, 0, 0])


def test_orthogonal_states_have_the_expected_distinct_metrics():
    metrics = state_error_metrics([1, 0], [0, 1j])
    assert metrics['fidelity'] == 0
    assert metrics['infidelity'] == pytest.approx(1)
    assert metrics['one_minus_overlap'] == pytest.approx(1)
    assert metrics['phase_aligned_distance'] == pytest.approx(np.sqrt(2))


def test_legacy_state_fields_remain_phase_sensitive_and_metrics_are_saved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    np.save('state.npy', [2.0, 0.0])
    result = _compute_state_dependent_errors(np.eye(2), 1j * np.eye(2), 'state.npy')
    record = result['state_errors'][0]
    assert record['absolute_error'] == pytest.approx(2 * np.sqrt(2))
    assert record['relative_error'] == pytest.approx(np.sqrt(2))
    assert record['one_minus_overlap'] < 1e-28
    assert record['input_norm'] == record['exact_norm'] == record['approx_norm'] == 2
    output = _save_error_results(result)
    with np.load(output) as saved:
        assert saved['state_input_files'].tolist() == ['state.npy']
        assert saved['state_absolute_errors'][0] == record['absolute_error']
        assert saved['state_one_minus_overlap'][0] == record['one_minus_overlap']
        assert saved['state_infidelity'][0] == record['infidelity']


def test_equal_nonunitary_scaling_is_not_hidden_by_normalized_overlap(tmp_path, caplog):
    path = tmp_path / 'state.npy'
    np.save(path, [1.0, 0.0])
    result = _compute_state_dependent_errors(2 * np.eye(2), 2 * np.eye(2), [str(path)])
    record = result['state_errors'][0]
    assert record['fidelity'] == 1
    assert record['exact_norm_defect'] == record['approx_norm_defect'] == 1
    assert 'changes the state norm' in caplog.text


def _operator(energies, vectors=None, time=0.1):
    return OperatorRepresentation(
        data={'eigenvalues': np.asarray(energies),
              'eigenvectors': np.eye(len(energies)) if vectors is None else vectors},
        operator_type='hamiltonian', representation='eigendecomposition', tevol_hbar=time,
    )


def test_eigenvector_matching_is_bijective_even_when_greedy_matching_reuses_a_state():
    q, _ = np.linalg.qr(np.random.default_rng(4).normal(size=(4, 4)))
    assert len(set(np.argmax(abs(q), axis=1))) == 3
    result = _compute_eigenvalue_errors(_operator([1., 2., 3., 4.]),
                                       _operator([1.1, 2.1, 3.1, 4.1], q))['eigenenergy_errors']
    assert sorted(result['matched_approx_indices']) == [0, 1, 2, 3]
    assert result['matching'] == 'bijective_maximum_squared_overlap'


def test_energy_phase_alias_is_explicit_and_zero_reference_relative_error_is_undefined():
    result = _compute_eigenvalue_errors(_operator([0., 1.], time=1.),
                                       _operator([0., 1. + 2 * np.pi], time=1.))['eigenenergy_errors']
    np.testing.assert_allclose(result['absolute_errors'], [0., 0.], atol=1e-14)
    assert result['energy_alias_period'] == 2 * np.pi
    assert result['relative_error_defined'] == [False, True]
    assert np.isnan(result['relative_errors'][0])
    assert 'principal' in result['phase_branch']


@pytest.mark.parametrize('time', [0., np.nan, np.inf])
def test_energy_error_rejects_undefined_time(time):
    with pytest.raises(ValueError, match='finite, nonzero'):
        _compute_eigenvalue_errors(_operator([1., 2.], time=time), _operator([1., 2.], time=time))


def test_energy_error_rejects_mismatched_time():
    with pytest.raises(ValueError, match='same tevol_hbar'):
        _compute_eigenvalue_errors(_operator([1., 2.], time=0.1), _operator([1., 2.], time=0.2))


def test_energy_shift_restoration_preserves_the_reference_comparison():
    time, shift = 0.3, -8.
    h = np.diag([1., 2.])
    exact = OperatorRepresentation(h, 'hamiltonian', tevol_hbar=time)
    approximate = OperatorRepresentation(expm(-1j * time * (h + shift * np.eye(2))),
                                         'time_evolution', energy_shifted=True,
                                         tevol_hbar=time, energy_shift=shift)
    result = _compute_eigenvalue_errors(exact, approximate)['eigenenergy_errors']
    np.testing.assert_allclose(result['absolute_errors'], [0, 0], atol=2e-14)
