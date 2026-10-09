"""Independent dense checks for the fixed-partition ordering experiment."""
import json
from pathlib import Path

import numpy as np
import pytest
from openfermion import FermionOperator
from scipy.linalg import expm

from qhat.analysis.fermionic_ownership import build_molecular_ownership
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.numerical_audit import controlled_block_orders as experiment

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], complex),
    'Y': np.array([[0, -1j], [1j, 0]], complex),
    'Z': np.diag([1., -1.]),
}


def literal_pauli(key, n):
    labels = dict(key)
    result = np.ones((1, 1), complex)
    for q in range(n):
        result = np.kron(result, PAULI[labels.get(q, 'I')])
    return result


def molecular_fixture():
    raw = FermionOperator()
    raw.terms = {
        ((0, 1), (0, 0)): .7,
        ((1, 1), (1, 0)): -.23,
        ((2, 1), (2, 0)): .11,
        ((3, 1), (3, 0)): -.19,
        ((0, 1), (2, 0)): .41,
        ((2, 1), (0, 0)): .41,
        ((1, 1), (3, 0)): -.31,
        ((3, 1), (1, 0)): -.31,
        ((0, 1), (1, 1), (1, 0), (0, 0)): .17,
        ((2, 1), (3, 1), (3, 0), (2, 0)): .29,
    }
    return raw


def test_fixed_membership_and_within_block_order_are_preserved_for_every_treatment():
    owned = build_molecular_ownership(molecular_fixture(), atol=1e-12)
    original = [(block.parent_index, block.pauli_keys) for block in owned.owned_blocks]
    metadata, orders, seeds = experiment.plan_block_orders(owned, 'independent-toy')
    assert len(orders) == 24
    assert len(seeds) == 20
    assert experiment.flatten_blocks(owned, orders['block_signed']) == list(owned.pauli_order)
    for name, indices in orders.items():
        assert sorted(indices) == list(range(len(original)))
        flattened = experiment.flatten_blocks(owned, indices)
        assert set(flattened) == set(owned.final_coefficients)
        assert len(flattened) == len(owned.final_coefficients)
        offset = 0
        for i in indices:
            expected = list(original[i][1])
            assert flattened[offset:offset+len(expected)] == expected
            offset += len(expected)
    assert [(block.parent_index, block.pauli_keys) for block in owned.owned_blocks] == original
    assert orders['block_parent_magnitude'] != orders['block_signed']
    assert orders['block_reverse_signed'] == list(reversed(orders['block_signed']))
    for block in metadata:
        parent = owned.parents[block['parent_index']]
        assert block['signed_parent_weight'] == parent.operator.terms[sorted(parent.component_keys)[0]].real


def test_random_permutations_reproducible_and_case_specific():
    owned = build_molecular_ownership(molecular_fixture(), atol=1e-12)
    _, first, seeds = experiment.plan_block_orders(owned, 'toy-a')
    _, again, _ = experiment.plan_block_orders(owned, 'toy-a')
    _, other, _ = experiment.plan_block_orders(owned, 'toy-b')
    assert first == again
    random_names = list(seeds)
    assert all(first[name] != other[name] for name in random_names)
    # Independent uniform draws may repeat a permutation (especially with only
    # five toy blocks); rejecting duplicates would change the sampling law.
    assert len({tuple(first[name]) for name in random_names}) > 1
    assert len({tuple(spec['seed_sequence_entropy']) for spec in seeds.values()}) == 20


@pytest.mark.parametrize('schedule', experiment.BLOCK_DETERMINISTIC + ('block_random_2026100901', 'block_random_2026100920'))
def test_full_pauli_kernel_matches_literal_dense_completed_block_exponentials(schedule):
    owned = build_molecular_ownership(molecular_fixture(), atol=1e-12)
    _, orders, _ = experiment.plan_block_orders(owned, 'independent-toy')
    n, total_time, steps = 4, .61, 3
    number = np.diag([i.bit_count() for i in range(2**n)])
    matrices = []
    for block in owned.owned_blocks:
        summands = [owned.final_coefficients[key]*literal_pauli(key, n) for key in block.pauli_keys]
        matrix = sum(summands, np.zeros((2**n, 2**n), complex))
        np.testing.assert_allclose(matrix @ number, number @ matrix, atol=1e-14)
        for a in summands:
            for b in summands:
                np.testing.assert_allclose(a @ b, b @ a, atol=1e-14)
        matrices.append(matrix)
    initial = np.zeros(2**n, complex)
    initial[[12, 9, 6, 3]] = np.array([1., 1j, .3, -.2j])
    initial /= np.linalg.norm(initial)
    dense = initial.copy()
    for _ in range(steps):
        for i in orders[schedule]:
            dense = expm(-1j*total_time/steps*matrices[i]) @ dense
    keys = experiment.flatten_blocks(owned, orders[schedule])
    terms = evolution.compile_ordered_terms(keys, owned.final_coefficients, n, 1e-12)
    actual = experiment.evolve(initial, terms, total_time, steps)
    np.testing.assert_allclose(actual, dense, atol=2e-14)
    assert np.linalg.norm(actual[[i for i in range(2**n) if i.bit_count() != 2]]) < 1e-14


def test_changing_block_chronology_can_change_the_approximation_at_fixed_partition():
    owned = build_molecular_ownership(molecular_fixture(), atol=1e-12)
    _, orders, _ = experiment.plan_block_orders(owned, 'independent-toy')
    initial = np.zeros(16, complex)
    initial[12] = 1.
    states = []
    for name in ('block_signed', 'block_parent_magnitude'):
        keys = experiment.flatten_blocks(owned, orders[name])
        terms = evolution.compile_ordered_terms(keys, owned.final_coefficients, 4, 1e-12)
        states.append(experiment.evolve(initial, terms, .61, 3))
    assert np.linalg.norm(states[0]-states[1]) > 1e-3


def test_declared_campaign_contains_728_evolutions_and_preselected_anchors():
    manifest = json.loads((experiment.FIXTURES/'manifest.json').read_text())
    settings = sum(len(experiment.planned_settings(case)) for case in manifest)
    assert settings == 28
    assert settings*(4+20+2) == 728
    assert set(experiment.CONVERGENCE_CASES) == {'HCN_s1.00', 'li2', 'f2'}
    assert len({experiment.family(case) for case in manifest}) == 7


def test_random_rank_ties_and_floor_are_used_consistently():
    summary = experiment.random_comparison(1e-4, [5e-5, 1e-4, 1e-4+5e-13, 2e-4])
    assert (summary['random_better'], summary['random_tied'], summary['random_worse']) == (1, 2, 1)
    assert summary['signed_midrank_among_signed_and_random'] == 3
    assert summary['signed_midpercentile_against_random'] == .5
    assert summary['signed_vs_random_median'] == 'tie'
    assert experiment.random_comparison(1e-14, [2e-14, 3e-14])['signed_over_random_median'] is None


def test_real_fixture_constructs_ownership_only_once_and_saves_design_before_checks(tmp_path, monkeypatch):
    calls = []
    original = experiment.build_molecular_ownership
    def counted(*args, **kwargs):
        calls.append(kwargs)
        if len(calls) > 1:
            raise AssertionError('Ownership was recomputed for a competing priority.')
        return original(*args, **kwargs)
    monkeypatch.setattr(experiment, 'build_molecular_ownership', counted)
    manifest = json.loads((experiment.FIXTURES/'manifest.json').read_text())
    _, _, schedules, _, _, design, checks = experiment.prepare('HCN_s1.00', manifest, tmp_path)
    assert len(calls) == 1
    assert calls[0]['ordering'] == 'signed_ascending'
    assert len(schedules) == 26
    assert len(set(design['schedule_membership_sha256'].values())) == 1
    assert (tmp_path/'HCN_s1.00.design.json').exists()
    assert max(checks['full_space_one_step_residuals'].values()) < 1e-11
