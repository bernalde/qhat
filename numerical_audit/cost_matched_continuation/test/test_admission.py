"""Cheap metadata regressions; no backend or scientific-runtime dependency."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


a = load('continuation_execution_admission_test', HERE/'admission.py')
analysis = load('continuation_primary_gate_test', HERE/'analyze.py')
PARENT = list(range(12))
CHANGED = [0, 2, 1]+list(range(3, 12))
EXPECTED_CHILD_STOPS = ('child_cpu_limit', 'child_wall_limit', 'child_evaluation_limit',
                        'child_action_limit', 'child_atom_limit')


def native(status='complete', seeds=4, evaluations=3000):
    return dict(order=list(PARENT), score=dict(value=1e-9, cubic_response=1e-5, steps=10),
                status=status, evaluations=evaluations, stats=dict(seed_evaluations=seeds))


def window(status='complete', child='complete', changed=False):
    result = native()
    result['order'] = list(CHANGED if changed else PARENT)
    stats = dict(method='window', parent_method='magnitude', parent_evaluations=1,
                 window_status=status, native_parent_revalidation_evaluations=1,
                 native_verification_attempts=1 if changed else 0,
                 native_verifications=1 if changed else 0)
    if child is not None:
        stats.update(backend_runtime_assets_verified=True, child_outcome=child,
            child_returncode=0 if child == 'complete' else -9,
            child_checkpoint_status='complete' if child == 'complete' else 'checked_update',
            child_evaluation_scope='checkpoint_prefix', child_completed_evaluations=10,
            child_allocated_evaluations=100, child_checked_updates=1 if changed else 0,
            child_counters=dict(evaluation=10, full_score=2 if changed else 1,
                flat_candidate=7, exact_witness=1, action=100, atom=10, insertion_move=1))
        if not changed and child != 'complete':
            stats['child_checkpoint_status'] = 'parent_checked'
    result['stats'] = stats
    return result


class ExecutionAdmissionTests(unittest.TestCase):
    def admit(self, result, arm='window', **kwargs):
        before = copy.deepcopy(result)
        outcome = a.execution_admission(arm, result, parent_order=PARENT, **kwargs)
        self.assertEqual(result, before)
        self.assertEqual(set(outcome), {'version', 'admitted', 'category', 'reason',
            'window_status', 'child_outcome', 'child_refusal'})
        json.dumps(outcome, allow_nan=False)
        return outcome

    def test_native_parent_declared_limit_and_strong_full_control(self):
        cases = [
            ('parent', native('complete'), True), ('parent', native('cpu_limit'), True),
            ('parent', native('evaluation_limit', evaluations=2999), True),
            ('parent', native('evaluation_limit', evaluations=3000), True),
            ('parent', native('evaluation_limit', evaluations=2998), False),
            ('parent', native('evaluation_limit', evaluations=3001), False),
            ('parent', native('support_limit'), False), ('parent', native('numerical_limit'), False),
            ('parent', native('unknown_error'), False), ('full', native('complete'), True),
            ('full', native('cpu_limit'), True), ('full', native('complete', seeds=3), False),
            ('full', native('cpu_limit', seeds=True), False), ('full', native('evaluation_limit'), False),
            ('full', native('support_limit'), False), ('full', native('numerical_limit'), False),
            ('full', native('unknown_error'), False)]
        for arm, result, expected in cases:
            with self.subTest(arm=arm, status=result['status'], stats=result['stats'], evaluations=result['evaluations']):
                self.assertEqual(self.admit(result, arm)['admitted'], expected)

    def test_completed_and_interrupted_checked_windows_remain_admitted(self):
        cases = [('complete_no_change', window()),
                 ('complete_improvement', window('verified_improvement', changed=True))]
        no_gain = window('native_no_resolved_improvement')
        no_gain['stats'].update(native_verification_attempts=1, native_verifications=1,
                                child_checked_updates=1)
        no_gain['stats']['child_counters']['full_score'] = 2
        cases.append(('native_checked_no_gain', no_gain))
        for child in EXPECTED_CHILD_STOPS:
            result = window('verified_checkpoint_after_'+child, child, changed=True)
            if child in ('child_evaluation_limit', 'child_action_limit', 'child_atom_limit'):
                result['stats'].update(child_returncode=2, child_refusal=child.removeprefix('child_'))
            cases.append((child, result))
        interrupted_no_gain = window('native_no_resolved_improvement', 'child_cpu_limit')
        interrupted_no_gain['stats'].update(native_verification_attempts=1, native_verifications=1,
            child_checked_updates=1, child_checkpoint_status='checked_update')
        interrupted_no_gain['stats']['child_counters']['full_score'] = 2
        cases.append(('interrupted_checked_no_gain', interrupted_no_gain))
        for label, result in cases:
            with self.subTest(case=label):
                self.assertTrue(self.admit(result)['admitted'])

    def test_checked_parent_budget_and_exact_capacity_stops(self):
        cases = [(status, window(status, None)) for status in ('whole_call_limit', 'evaluation_limit',
            'insufficient_evaluation_budget', 'verification_reserve_unavailable')]
        for child in EXPECTED_CHILD_STOPS:
            result = window(child, child)
            for key in ('child_counters', 'child_completed_evaluations', 'child_checkpoint_status'):
                result['stats'].pop(key)
            cases.append((child+'_without_checkpoint', result))
        for reason in ('reachable vector exceeds max_entries', 'cached reference columns exceed max_entries',
                       'reachable unit column exceeds max_entries', 'combined cubic columns exceed max_entries',
                       'updated cubic columns exceed max_entries', 'window rational exceeds 4096-bit admission cap',
                       'window request exceeds 8 MiB'):
            result = window('support_limit', None)
            result['stats']['refusal_reason'] = reason
            cases.append((reason, result))
        for label, result in cases:
            with self.subTest(case=label):
                self.assertTrue(self.admit(result)['admitted'])

    def test_scored_parents_do_not_hide_window_failures(self):
        cases = [(status, window(status, None)) for status in sorted(a.RUNTIME_FAILURES | a.DOMAIN_REFUSALS)]
        cases += [(status, window(status)) for status in sorted(a.REPLY_FAILURES)]
        cases += [('launch', window('child_launch_or_transport_failure', None)),
                  ('numerical', window('numerical_limit', None)),
                  ('unknown', window('future_status')), ('no_status', window(None, None))]
        for child in ('child_failure', 'child_worker_failure', 'child_numerical_limit',
                      'child_admission_or_support_limit', 'owned_child_cpu_unavailable', 'future_child_outcome'):
            cases.append(('masked_'+child, window('native_no_resolved_improvement', child)))
            cases.append(('prefixed_'+child, window('verified_checkpoint_after_'+child, child, changed=True)))
        for reason in ('installed Julia worker assets are missing', 'backend source manifest missing or unreadable',
                       'malformed backend source manifest', 'backend runtime asset byte cap: worker.jl',
                       'backend runtime asset SHA-256 mismatch: worker.jl', 'unrecognized support cap'):
            result = window('support_limit', None)
            result['stats']['refusal_reason'] = reason
            cases.append((reason, result))
        for label, result in cases:
            with self.subTest(case=label):
                self.assertFalse(self.admit(result)['admitted'])

    def test_checkpoint_metadata_and_status_combinations_fail_closed(self):
        original = window('verified_checkpoint_after_child_cpu_limit', 'child_cpu_limit', changed=True)
        changes = dict(backend_runtime_assets_verified=False, child_checkpoint_status=None,
            child_counters=None, child_completed_evaluations=None, child_checked_updates=0,
            native_parent_revalidation_evaluations=True, native_verifications=0,
            native_verification_attempts=0, child_outcome=None, child_refusal='worker_failure',
            child_evaluation_scope='unvalidated', parent_evaluations=0, parent_method='sparse')
        for field, value in changes.items():
            with self.subTest(field=field):
                result = copy.deepcopy(original)
                result['stats'][field] = value
                self.assertFalse(self.admit(result)['admitted'])
        contradictory = [window('complete', 'child_cpu_limit'), window('child_cpu_limit', 'child_cpu_limit', changed=True),
            window('verified_improvement', 'child_cpu_limit', changed=True)]
        result = window('verified_improvement', changed=True)
        result['stats']['child_returncode'] = -9
        contradictory.append(result)
        for result in contradictory:
            with self.subTest(status=result['stats']['window_status']):
                self.assertFalse(self.admit(result)['admitted'])
        self.assertFalse(a.execution_admission('window', original, parent_order=None)['admitted'])

    def test_refusals_and_invalid_incumbents_are_json_safe_and_ineligible(self):
        self.assertFalse(self.admit(window(), refusal=dict(type='OSError', message='retained'))['admitted'])
        self.assertFalse(a.execution_admission('window', None, parent_order=PARENT)['admitted'])
        for value in (None, float('inf'), float('nan'), True, -1., 10**1000):
            result = window()
            if value is None:
                result['score'] = None
            else:
                result['score']['value'] = value
            outcome = a.execution_admission('window', result, parent_order=PARENT)
            self.assertFalse(outcome['admitted'])
            json.dumps(outcome, allow_nan=False)

    def test_primary_gate_preserves_all_window_losses_and_valid_checkpoints(self):
        config = dict(total_CPU_budgets=[32, 64], physical=dict(times=[.65, 1.35],
            rotation_cap_multipliers=[8, 24], ambiguity=1e-12, ratio_floor=1e-12,
            primary=dict(total_CPU_budget=64, time=1.35, rotation_cap_multiplier=8),
            repetitions_scope='Deterministic timing repetitions; no confidence intervals'))
        cases, constructors, physical = ['A', 'B'], {}, {}
        for case in cases:
            for rep in (1, 2):
                for arm in ('full', 'window'):
                    for budget in (32, 64):
                        constructors[case, rep, arm, budget] = dict(analysis_eligible=True)
                        for duration in (.65, 1.35):
                            for cap in (8, 24):
                                physical[case, arm, rep, budget, duration, cap] = dict(
                                    projective_error=.01 if arm == 'full' else .02)
        verified = set(physical)
        original = analysis.pair_grid(cases, 2, config, physical, constructors, verified)
        self.assertEqual(analysis.primary_decision(original, cases, 2, config, True)['decision'], 'full_advantage')
        for case in cases:
            for rep in (1, 2):
                constructors[case, rep, 'window', 64]['analysis_eligible'] = self.admit(
                    window('julia_unavailable', None))['admitted']
        rejected = analysis.pair_grid(cases, 2, config, physical, constructors, verified)
        self.assertEqual(len(rejected), 32)
        self.assertEqual([(row['difference'], row['outcome']) for row in rejected],
                         [(row['difference'], row['outcome']) for row in original])
        self.assertTrue(all(row['outcome'] == 'loss' for row in rejected))
        self.assertEqual(analysis.primary_decision(rejected, cases, 2, config, False)['decision'], 'study_incomplete')
        for case in cases:
            for rep in (1, 2):
                constructors[case, rep, 'window', 64]['analysis_eligible'] = self.admit(
                    window('verified_checkpoint_after_child_cpu_limit', 'child_cpu_limit', changed=True))['admitted']
        restored = analysis.pair_grid(cases, 2, config, physical, constructors, verified)
        self.assertEqual(analysis.primary_decision(restored, cases, 2, config, True)['decision'], 'full_advantage')


if __name__ == '__main__':
    unittest.main()
