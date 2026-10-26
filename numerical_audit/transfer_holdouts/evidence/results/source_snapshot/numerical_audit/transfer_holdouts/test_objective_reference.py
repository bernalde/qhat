"""Independent coherent-objective and exact legacy-decision checks."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
import objective_reference
import training
from cubic_response import Context
from test_states import toy
from test_cubic_response import model as interference_model, state as interference_state


@pytest.mark.parametrize('phase', [1, 1j])
def test_literal_nested_commutators_match_all_supported_cubic_columns_and_coherent_score(phase):
    model = toy()
    psi = np.zeros(16, complex)
    psi[[3, 6, 9, 12]] = [1, phase, -.2, .3j]
    psi /= np.linalg.norm(psi)
    reference = objective_reference.Reference(model, psi)
    producer = Context(model, psi)
    for order in [list(range(len(model['blocks']))), list(reversed(range(len(model['blocks']))))]:
        actual = producer.score(order, .83, 270)
        checked = reference.cubic_score(order, .83, 270)
        assert np.max(abs(np.asarray(actual['cubic_columns']) - checked['cubic_columns'])) < 2e-13
        assert actual['column_sources'] == checked['column_sources']
        for key in ('JD2', 'score', 'steps'):
            assert actual[key] == pytest.approx(checked[key], rel=2e-11, abs=2e-18)


@pytest.mark.parametrize('partition', ['conserving', 'pauli'])
def test_literal_first_order_complex_transport_matches_feature_objective(partition):
    model = toy()
    psi = np.zeros(16, complex)
    psi[[3, 6, 9, 12]] = [1, 1j, -.2, .3j]
    psi /= np.linalg.norm(psi)
    reference = objective_reference.Reference(model, psi)
    producer = training.FirstOrderContext(model, psi, partition)
    for order in [list(range(producer.m)), list(reversed(range(producer.m)))]:
        actual = producer.score(order, training.k.TRAIN_TIME, 270)
        checked = reference.first_score(order, training.k.TRAIN_TIME, 270, partition)
        for key in ('JD1', 'score', 'steps'):
            assert actual[key] == pytest.approx(checked[key], rel=2e-11, abs=2e-18)


def test_all_legacy_slot_decisions_match_independent_integer_bch_enumeration():
    model = toy()
    psi = np.zeros(16, complex)
    psi[3] = 1
    actual = training.current_construction(model)
    checked = objective_reference.Reference(model, psi).current_placement()
    for key in ('order', 'start_slot', 'slots', 'history', 'initial_slot_scores', 'exact_weight_denominator',
                'exact_score_numerator', 'evaluations', 'diagonal_priority', 'off_priority'):
        assert actual[key] == checked[key]
    assert actual['score'] == pytest.approx(checked['score'], rel=2e-11, abs=2e-18)


def test_literal_coherent_score_is_not_an_incoherent_state_mixture():
    model = interference_model()
    psi = interference_state()
    order = [2, 0, 4, 1, 3]
    coherent = objective_reference.Reference(model, psi).cubic_score(order, .83, 270)['JD2']
    mixture = 0.
    for index in np.flatnonzero(psi):
        basis = np.eye(len(psi), dtype=complex)[:, index]
        mixture += abs(psi[index])**2 * objective_reference.Reference(model, basis).cubic_score(order, .83, 270)['JD2']
    assert abs(coherent - mixture) > 1e-7
