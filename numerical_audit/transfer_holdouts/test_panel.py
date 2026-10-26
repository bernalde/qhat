"""Full declared transfer cohort and sharing/cost invariants without a solve."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pytest
import panel
import transfer_common as k


def test_complete_main_panel_keeps_all_states_times_arms_and_comparisons():
    cases = [row['case'] for row in k.read(k.HERE / 'REGISTRY.json')['cases']]
    plan = panel.build_plan(cases)
    assert (plan['expected_constructors'], plan['expected_unique_circuits'],
            plan['expected_logical_circuits'], plan['expected_reference_evolutions'],
            plan['expected_comparisons']) == (96, 1332, 1440, 90, 5040)
    cohort = {(row['case'], row['state'], row['T'], row['cap_multiple'], row['arm']) for row in plan['logical']}
    assert len(cohort) == 1440
    assert all((case, state, T, cap, arm) in cohort for case in cases for state in k.STATES
               for T in k.TIMES for cap in k.CAPS for arm in k.ARMS)
    assert set(plan['validation_times']).isdisjoint({k.TRAIN_TIME, 1.})


def test_shared_HF_roles_reuse_exact_work_without_collapsing_logical_observations():
    plan = panel.build_plan(['toy'])
    rows = {(row['state'], row['T'], row['cap_multiple'], row['arm']): row for row in plan['logical']}
    for T in k.TIMES:
        for cap in k.CAPS:
            for kind in ('first_JD', 'pauli_JD', 'cubic_JD'):
                hf = rows['HF', T, cap, 'hf_' + kind]
                same = rows['HF', T, cap, 'state_' + kind]
                assert hf['constructor'] == same['constructor'] and hf['candidate'] == same['candidate']
            assert rows['HF', T, cap, 'state_first_symmetric']['candidate'] != rows['HF', T, cap, 'state_first_JD']['candidate']
    ground = rows['ground', k.TIMES[0], k.CAPS[0], 'state_first_JD']
    frozen = rows['ground', k.TIMES[0], k.CAPS[0], 'hf_first_JD']
    assert ground['constructor'] != frozen['constructor'] and ground['candidate'] != frozen['candidate']
    assert ground['reference_charge_states'] == ['ground']
    assert frozen['reference_charge_states'] == ['HF', 'ground']


def test_first_order_order_is_shared_across_formula_arms_but_compilation_is_separate():
    plan = panel.build_plan(['toy'])
    rows = [row for row in plan['logical'] if row['state'] == 'coherent_imag' and
            row['T'] == k.TIMES[0] and row['cap_multiple'] == k.CAPS[0]]
    first = next(row for row in rows if row['arm'] == 'state_first_JD')
    symmetric = next(row for row in rows if row['arm'] == 'state_first_symmetric')
    assert first['constructor'] == symmetric['constructor']
    assert first['candidate'] != symmetric['candidate']
    assert all('evolution only' in row['compiled_resource_scope'] for row in rows)


def test_full_path_preflight_uses_all_five_states_and_controls_at_largest_time():
    plan = panel.build_plan(['preflight_LiH_r1p841'], preflight=True)
    assert (plan['expected_constructors'], plan['expected_unique_circuits'],
            plan['expected_logical_circuits'], plan['expected_reference_evolutions'],
            plan['expected_comparisons']) == (16, 74, 80, 5, 280)
    assert plan['validation_times'] == [max(k.TIMES)]
    assert plan['states'] == list(k.STATES) and plan['arms'] == list(k.ARMS)


def test_case_and_state_changes_prevent_unjustified_deduplication():
    plan = panel.build_plan(['case_a', 'case_b'])
    candidates = {row['candidate'] for row in plan['circuits']}
    assert len(candidates) == 444
    assert all(row['candidate'] in candidates for row in plan['logical'])
    assert len({row['constructor'] for row in plan['training']}) == 32


@pytest.mark.parametrize('cases', [[], ['same', 'same'], ['../escape'], ['case'] * 7])
def test_invalid_case_inventory_refused_before_plan_allocation(cases):
    with pytest.raises(ValueError):
        panel.build_plan(cases)


def test_floor_loss_tie_and_unavailable_observations_are_retained():
    complete = lambda error: dict(status='complete', error=error)
    assert panel.compare(complete(1e-14), complete(2e-14))['decision'] == 'floor_tie'
    assert not panel.compare(complete(1e-14), complete(2e-14))['informative']
    assert panel.compare(complete(1e-5), complete(1e-5 + 1e-13))['decision'] == 'numerical_tie'
    assert panel.compare(complete(1e-5), complete(2e-5))['decision'] == 'a'
    assert panel.compare(complete(2e-5), complete(1e-5))['decision'] == 'b'
    missing = panel.compare(dict(status='censored'), complete(.1))
    assert missing['decision'] == 'unavailable' and missing['a_error'] is None and missing['b_error'] == .1


def test_nonfinite_error_cannot_be_classified_as_a_win():
    with pytest.raises(ValueError, match='Invalid comparison errors'):
        panel.compare(dict(status='complete', error=float('nan')), dict(status='complete', error=.1))
