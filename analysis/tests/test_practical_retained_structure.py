"""Audit retained practical layouts/dispatch without producer bookkeeping code."""

from collections import defaultdict, deque
import copy
import hashlib
import json
from pathlib import Path

import pytest


RESULTS = (Path(__file__).resolve().parents[2] / 'numerical_audit' /
           'practical_accuracy' / 'evidence' / 'results')


def read(path):
    return json.loads(path.read_text())


def canonical_terms(rows, n_qubits):
    coefficients = defaultdict(int)
    for raw_word, coefficient in rows:
        word = tuple((q, letter) for q, letter in raw_word)
        assert word == tuple(sorted(word))
        assert len({q for q, _ in word}) == len(word)
        assert all(type(q) is int and 0 <= q < n_qubits and
                   letter in ('X', 'Y', 'Z') for q, letter in word)
        coefficients[word] += int(coefficient)
    return [(word, str(value)) for word, value in sorted(coefficients.items()) if value]


def check_layout(model, record):
    assert record['formula'] == 2 and record['T'] == 2.11
    assert canonical_terms([term for block in model['blocks'] for term in block],
                           model['n']) == canonical_terms(model['coefficients'], model['n'])
    if record['partition'] == 'conserving':
        units = [canonical_terms(block, model['n']) for block in model['blocks']]
    else:
        assert record['partition'] == 'pauli'
        units = [[term] for term in canonical_terms(model['coefficients'], model['n'])]
    assert all(units)
    order = record['order']
    assert all(type(index) is int for index in order)
    assert sorted(order) == list(range(len(units)))
    expected = [[[list(pair) for pair in word], value]
                for index in order for word, value in units[index]]
    assert record['schedule'] == expected
    assert record['unit_lengths'] == [len(units[index]) for index in order]


def check_curve_model(curve, model_path):
    assert hashlib.sha256(model_path.read_bytes()).hexdigest() == curve['model_sha256']
    check_layout(read(model_path), curve['record'])


def test_all_81_retained_complete_unit_layouts_have_literal_independent_checks():
    curves = sorted((RESULTS / 'curves').glob('*.json'))
    assert len(curves) == 81
    models = {path.stem: path for path in (RESULTS / 'models').glob('*.json')}
    for path in curves:
        curve = read(path)
        check_curve_model(curve, models[curve['case']])


def test_retained_dispatch_matches_the_declared_cyclic_case_and_repetition_order():
    protocol = read(RESULTS / 'PROTOCOL_FROZEN.json')
    cases = ['Li2_r2p95', 'F2_r1p55', 'H2O_scale1p10']
    methods = ['first_JD', 'pauli_JD', 'cubic_JD', 'cubic_warm_JD']
    references = ['current', 'pauli_JD_large_sameT', 'cubic_insertion']
    assert protocol['original_cases'] == cases and protocol['methods'] == methods
    assert protocol['reference_arms'] == references
    assert protocol['repetitions'] == 3 and protocol['matched_total_CPU_budgets'] == [2., 8.]
    assert protocol['controller_policy_version'] == 2
    expected = []
    for case_index, case in enumerate(cases):
        cycle = deque(methods)
        cycle.rotate(-case_index)
        for repetition in range(3):
            if repetition:
                cycle.rotate(-1)
            for budget in (2, 8):
                for method in cycle:
                    expected.append(f'{case}.{method}.b{budget}.rep{repetition}')
        expected.extend(f'{case}.{method}.reference' for method in references)
    assert read(RESULTS / 'CONSTRUCTION_PLAN_FROZEN.json') == expected
    assert {path.stem for path in (RESULTS / 'curves').glob('*.json')} == set(expected)
    for index, name in enumerate(expected):
        curve = read(RESULTS / 'curves' / (name + '.json'))
        assert curve['dispatch_index'] == index and curve['curve'] == name


@pytest.mark.parametrize('corruption', ['missing_rotation', 'duplicate_unit'])
def test_literal_layout_audit_rejects_incomplete_or_duplicated_units(corruption):
    curve = read(sorted((RESULTS / 'curves').glob('*.json'))[0])
    model = read(RESULTS / 'models' / (curve['case'] + '.json'))
    broken = copy.deepcopy(curve['record'])
    if corruption == 'missing_rotation':
        broken['schedule'].pop()
    else:
        broken['order'][0] = broken['order'][1]
    with pytest.raises(AssertionError):
        check_layout(model, broken)


def test_literal_layout_audit_rejects_inconsistent_hamiltonian_blocks():
    curve = read(sorted((RESULTS / 'curves').glob('*.json'))[0])
    model = read(RESULTS / 'models' / (curve['case'] + '.json'))
    model['coefficients'][0][1] = str(int(model['coefficients'][0][1]) + 1)
    with pytest.raises(AssertionError):
        check_layout(model, curve['record'])


def test_curve_rejects_model_with_different_bytes(tmp_path):
    curve = read(sorted((RESULTS / 'curves').glob('*.json'))[0])
    original = RESULTS / 'models' / (curve['case'] + '.json')
    changed = tmp_path / 'model.json'
    changed.write_bytes(original.read_bytes() + b' ')
    with pytest.raises(AssertionError):
        check_curve_model(curve, changed)
