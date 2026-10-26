"""Adversarial checks of independent retained-witness and complete-cohort audits."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys
import zipfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
import objective_reference
import panel
import training
from test_states import toy

_spec = importlib.util.spec_from_file_location('transfer_independent_verifier_tests', Path(__file__).with_name('verify.py'))
audit = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = audit
_spec.loader.exec_module(audit)


def counters():
    return dict(constructor_checks=0, main_constructors=0, warmup_constructors=0,
                nested_warm_constructors=0, objective_records=0, complete=0, censored=0)


def test_independently_rederived_complete_panel_matches_every_declared_identity():
    cases = [row['case'] for row in training.k.read(training.k.HERE / 'REGISTRY.json')['cases']]
    constructors, candidates, logical, times = audit.expected_panel(cases, False)
    expected = panel.build_plan(cases)
    assert (len(constructors), len(candidates), len(logical), len(times)) == (96, 1332, 1440, 3)
    assert constructors == expected['training'] and candidates == expected['circuits'] and logical == expected['logical']
    small = audit.expected_panel(['preflight_LiH_r1p841'], True)
    assert tuple(len(rows) for rows in small[:3]) == (16, 74, 80)


@pytest.fixture(scope='module')
def first_witness(tmp_path_factory):
    model, path = toy(), tmp_path_factory.mktemp('literal-first') / 'model.json'
    training.k.write(path, model)
    psi = np.eye(16, dtype=complex)[:, 3]
    record = training.first(path, training.k.sha(path), psi, 'conserving', reference_CPU_charge=.25)
    # Actual retained constructors must survive the same JSON round trip used
    # by the study, including every seed, sweep and accepted-score witness.
    record = json.loads(json.dumps(record, allow_nan=False))
    context = objective_reference.Reference(model, psi)
    return model, psi, context, context.current_placement(), record


def test_independent_audit_accepts_actual_retained_first_order_search(first_witness):
    _, _, context, current, record = first_witness
    checked = counters()
    audit.check_constructor(context, record, 'first_JD', current, .25, checked)
    assert checked['constructor_checks'] == 1 and checked['objective_records'] >= 7


def test_fabricated_seed_score_is_rejected_by_literal_operator(first_witness):
    _, _, context, current, original = first_witness
    record = copy.deepcopy(original)
    record['seeds'][0]['JD1'] += 1.
    with pytest.raises(ValueError, match='Independent numeric witness differs'):
        audit.check_constructor(context, record, 'first_JD', current, .25, counters())


def test_invalid_retained_insertion_cannot_be_accepted_as_a_search_witness(first_witness):
    _, _, context, current, original = first_witness
    record = copy.deepcopy(original)
    record['result']['starts'][0]['accepted'] = [dict(move=[0, 0], evaluations=0, record=record['seeds'][0])]
    with pytest.raises(ValueError, match='Invalid insertion witness'):
        audit.check_constructor(context, record, 'first_JD', current, .25, counters())


def test_wrong_reference_charge_is_rejected_even_when_search_scores_match(first_witness):
    _, _, context, current, original = first_witness
    with pytest.raises(ValueError, match='Independent numeric witness differs'):
        audit.check_constructor(context, original, 'first_JD', current, .5, counters())


def test_actual_cubic_warm_constructor_is_checked_and_reference_charged_once(tmp_path):
    model, path = toy(), tmp_path / 'model.json'
    training.k.write(path, model)
    psi = np.eye(16, dtype=complex)[:, 3]
    record = training.cubic(path, training.k.sha(path), psi, reference_CPU_charge=.25)
    record = json.loads(json.dumps(record, allow_nan=False))
    context = objective_reference.Reference(model, psi)
    checked = counters()
    audit.check_constructor(context, record, 'cubic_JD', context.current_placement(), .25, checked)
    assert checked['constructor_checks'] == 2 and checked['nested_warm_constructors'] == 1
    assert record['warm_start']['reference_preparation_CPU_charge'] == 0.
    assert checked['objective_records'] >= 14


@pytest.mark.parametrize('formula', [1, 2])
@pytest.mark.parametrize('lengths', [[2], [2, 1], [1, 2, 3], [3, 1, 2]])
def test_independent_expanded_affordability_respects_both_rotation_caps(formula, lengths):
    import compilation
    for cap in (10, 30):
        steps, rotations = audit.steps_and_rotations(lengths, formula, cap)
        assert steps == compilation.affordable_steps(lengths, cap, formula)
        assert rotations == compilation.rotation_count(lengths, steps, formula) <= cap * sum(lengths)


def test_independent_floor_loss_unavailable_and_nan_rules():
    complete = lambda error: dict(status='complete', error=error)
    assert audit.decision(complete(1e-14), complete(2e-14))['decision'] == 'floor_tie'
    assert audit.decision(complete(.1), complete(.2))['decision'] == 'a'
    assert audit.decision(complete(.2), complete(.1))['decision'] == 'b'
    assert audit.decision(dict(status='censored'), complete(.1))['decision'] == 'unavailable'
    with pytest.raises(ValueError, match='Independent finite comparison errors'):
        audit.decision(complete(float('nan')), complete(.1))


def test_independent_array_loader_rejects_oversized_shape_before_numpy_allocation(tmp_path, monkeypatch):
    raw, path = io.BytesIO(), tmp_path / 'vector.npz'
    np.lib.format.write_array_header_1_0(raw, dict(descr=np.dtype('complex128').str,
                                                fortran_order=False, shape=(1 << 40,)))
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('vector.npy', raw.getvalue())
    def forbidden(*args, **kwargs):
        pytest.fail('Independent allocation must follow the bounded NPY header')
    monkeypatch.setattr(audit.np, 'load', forbidden)
    with pytest.raises(ValueError, match='State array shape/dtype before allocation'):
        audit.array(path, 'vector', (16,), audit.sha(path))


def test_layout_coefficient_substitution_rejected_before_circuit_audit():
    model = toy()
    record = training.layout(model, list(range(len(model['blocks']))), 'conserving')
    record['schedule'][0][1] = '999'
    with pytest.raises(ValueError, match='Frozen layout differs'):
        audit.layout(model, record, 'conserving')


@pytest.mark.parametrize('status', ['complete', 'censored'])
def test_candidate_progress_keeps_verification_and_physical_status_separate(tmp_path, status):
    physical = dict(status=status)
    if status == 'complete':
        physical.update(error=.1, infidelity=.19, norm=1., number_leakage=0., joint_sector_leakage=0.)
    original = copy.deepcopy(physical)
    audit.write_candidate_progress(tmp_path, 'candidate', physical)
    progress = audit.read(tmp_path / 'progress/candidates/candidate.json')
    assert progress['status'] == 'passed'
    assert progress['physical'] == original
    assert physical == original
    other = dict(status='complete', error=.2)
    assert audit.decision(physical, other)['decision'] == ('a' if status == 'complete' else 'unavailable')
