"""Regressions for the actual L-sweep molecular benchmark's metric contract."""

import csv
import io

import numpy as np
import pytest

from qhat.analysis import benchmark_b2_active_spaces_matrix_free as matrixfree
from qhat.analysis import benchmark_b2_coloring_robustness as robustness
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis import benchmark_fermionic_structure_ablation as ablation


def _compare(approx, **kwargs):
    return matrixfree.compare_states(np.array([1.0]), np.array([1]), np.asarray(approx),
                                     np.array([1, 2]), **kwargs)


def test_leakage_is_part_of_full_state_error_not_renormalized_away():
    leakage = 0.25
    result = _compare([np.sqrt(leakage), np.sqrt(1 - leakage), 0., 0.])
    assert result['particle_number_leakage'] == pytest.approx(leakage)
    assert result['spin_sector_leakage'] == pytest.approx(leakage)
    assert result['state_infidelity'] == pytest.approx(leakage)
    assert result['one_minus_overlap'] == pytest.approx(1 - np.sqrt(1 - leakage))
    # Renormalizing the projection alone would incorrectly return zero error.
    assert result['one_minus_overlap'] > 0.1


def test_number_and_spin_sector_leakage_are_distinguished():
    result = _compare([0., np.sqrt(0.75), 0.5, 0.])
    assert result['particle_number_leakage'] == 0
    assert result['spin_sector_leakage'] == pytest.approx(0.25)
    assert result['state_infidelity'] == pytest.approx(0.25)


def test_zero_projected_sector_weight_is_a_valid_orthogonal_full_state():
    result = _compare([1., 0., 0., 0.])
    assert result['state_overlap_abs'] == 0
    assert result['one_minus_overlap'] == pytest.approx(1)
    assert result['state_infidelity'] == pytest.approx(1)
    assert result['particle_number_leakage'] == 1


def test_global_phase_changes_only_phase_sensitive_distance():
    phase = _compare([0., 1j, 0., 0.])
    assert phase['state_vector_2norm_error'] == pytest.approx(np.sqrt(2))
    assert phase['one_minus_overlap'] == 0
    assert phase['state_infidelity'] == 0


def test_tiny_error_and_leakage_remain_positive_when_overlap_rounds_to_one():
    delta = 1e-10
    result = _compare([delta, 1., 0., 0.])
    assert result['state_overlap_abs'] == 1
    assert result['one_minus_overlap'] == pytest.approx(delta**2 / 2, rel=1e-12, abs=0)
    assert result['state_infidelity'] == pytest.approx(delta**2, rel=1e-12, abs=0)
    assert result['particle_number_leakage'] == pytest.approx(delta**2, rel=1e-12, abs=0)


@pytest.mark.parametrize('approx', [[0., 1.001, 0., 0.], [0., 0.99, 0., 0.]])
def test_norm_drift_is_not_converted_to_an_apparent_perfect_overlap(approx):
    with pytest.raises(ValueError, match='norm drift'):
        _compare(approx)


def test_reference_norm_drift_is_also_rejected():
    with pytest.raises(ValueError, match='norm drift'):
        matrixfree.compare_states(np.array([1.001]), np.array([1]), np.array([0., 1., 0., 0.]),
                                  np.array([1, 2]))


def test_roundoff_norm_defects_are_recorded_without_biasing_ray_error():
    result = _compare([0., 1. + 1e-12, 0., 0.])
    assert result['one_minus_overlap'] == 0
    assert result['approximate_state_norm_defect'] > 0
    assert result['state_norm_valid']
    assert result['state_metric_version'] == matrixfree.STATE_METRIC_VERSION


@pytest.mark.parametrize('indices', [np.array([1, 1]), np.array([4]), np.array([1.5])])
def test_invalid_sector_indices_are_rejected(indices):
    with pytest.raises(ValueError):
        matrixfree.compare_states(np.ones(indices.size), indices, np.array([0., 1., 0., 0.]),
                                  np.array([1, 2]))


@pytest.mark.parametrize('module', [matrixfree, robustness, baseline, ablation])
def test_all_metric_consumers_can_write_complete_records_without_duplicate_columns(module):
    result = _compare([0., 1., 0., 0.])
    assert len(module.FIELDNAMES) == len(set(module.FIELDNAMES))
    assert set(result) <= set(module.FIELDNAMES)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=module.FIELDNAMES)
    writer.writeheader()
    writer.writerow(result)
    stream.seek(0)
    restored = next(csv.DictReader(stream))
    assert restored['state_metric_version'] == matrixfree.STATE_METRIC_VERSION
    assert restored['one_minus_overlap'] == '0.0'


@pytest.mark.parametrize('loader', [matrixfree.load_resume_data, robustness.load_resume_data, ablation.load_completed])
def test_old_metric_csv_cannot_be_silently_resumed(loader, tmp_path):
    old = tmp_path / 'historical.csv'
    old.write_text('status,state_overlap_abs,state_infidelity\nsuccess,1.0,0.0\n')
    with pytest.raises(ValueError, match='schema|metric'):
        loader(old)
