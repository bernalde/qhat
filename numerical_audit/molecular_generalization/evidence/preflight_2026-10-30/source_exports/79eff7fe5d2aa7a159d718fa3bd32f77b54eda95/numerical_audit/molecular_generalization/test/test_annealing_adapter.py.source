"""Tiny native parity and synthetic owned-child accounting; no Julia launch.

Run with the pinned scientific Python under the coordinator's OS resource guard.
Injected child witnesses/costs below are test fixtures, never study evidence.
"""
import copy
from fractions import Fraction
import importlib.machinery
import importlib.util
import math
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
if 'qhat' not in sys.modules:
    spec = importlib.machinery.ModuleSpec('qhat', loader=None, is_package=True)
    spec.submodule_search_locations = [str(ROOT)]
    sys.modules['qhat'] = importlib.util.module_from_spec(spec)
spec = importlib.util.spec_from_file_location('molecular_annealing_adapter_test', HERE/'annealing.py')
a = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = a
spec.loader.exec_module(a)
native = a._native()


def settings():
    return dict(calibrated=True, calibration_sha256='a'*64, theta_start=.5, theta_end=.005)


def write_toml(path, data):
    lines = [key+'='+a._toml(value) for key, value in data.items() if not isinstance(value, dict)]
    for name, mapping in data.items():
        if isinstance(mapping, dict):
            lines.append('['+name+']')
            lines.extend(a._toml(key)+'='+a._toml(value) for key, value in mapping.items())
    path.write_text('\n'.join(lines)+'\n', encoding='ascii')


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.blocks = ({((0, 'X'),): Fraction(3, 4)},
                       {((0, 'Z'),): Fraction(-1, 2)}, {(): Fraction(1, 8)})
        self.context = native.SparseResponse(self.blocks, 1, 0, time=1.)
        left, right = (0, 1, 2), (1, 0, 2)
        pairs = [(self.context.score(p, 90).value, p) for p in (left, right)]
        self.parent = max(pairs)[1]
        self.improved = min(pairs)[1]
        self.score = self.context.score(self.parent, 90)
        rest = [p for p in ((0, 1, 2), (1, 0, 2), (2, 0, 1), (2, 1, 0), (0, 2, 1)) if p != self.parent]
        self.starts = [self.parent]+rest[:3]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.worker_calls = 0
        self.last_reply = None

    def call(self, **changes):
        arguments = dict(parent_order=self.parent, starts=self.starts, time=1., rotation_budget=90,
            parent_score=self.score, settings=settings(), seed_material='tiny:1', cpu_seconds=10.,
            wall_seconds=50., checkpoint_path=self.directory/'checkpoint.toml')
        arguments.update(changes)
        return a.train_annealing(self.blocks, 1, 0, **arguments)

    def fake_worker(self, executable, request, checkpoint, child_cap, wall_cap, check, stats):
        self.worker_calls += 1
        self.assertGreater(child_cap, 0)
        self.assertLess(child_cap, 10.)
        check(0.)
        data = a.tomllib.loads(request.read_text())
        score = self.context.score(self.improved, 90)
        reply = dict(version=1, request_sha256=a.hashlib.sha256(request.read_bytes()).hexdigest(),
            settings_sha256=data['settings_sha256'], coefficients_sha256=data['coefficients_sha256'],
            source_sha256=data['source_sha256'], word_lengths=[1, 1, 1], time=1., budget=90,
            julia_version='test-fixture-no-runtime', rng='Random.Xoshiro',
            score_scope='full_order_transported_projected_cubic_J_not_error_bound',
            phases=dict(request_CPU_seconds=.01, context_CPU_seconds=.02), status='cpu_limit',
            order=list(self.improved), score=score.value, cubic_response=score.cubic_response, steps=score.steps,
            start_scores_completed=4, evaluations=8, proposals=4, accepted_moves=1,
            uphill_attempts=1, uphill_accepted=0, chain_proposals=[1, 1, 1, 1],
            chain_accepted=[1, 0, 0, 0], chain_CPU_seconds=[.01]*4, seeds=data['seeds'], stop_reason='fixture',
            kernel_stats=dict(unit_actions=20, cache_hits=10, cache_misses=10, cache_skipped=0,
                cache_entries=11, cache_accounted_bytes=11264, kernel_hits=1, kernel_misses=1,
                peak_vector_entries=2, peak_working_entries=10, full_scores=8))
        self.last_reply = copy.deepcopy(reply)
        write_toml(checkpoint, reply)
        stats.update(child_process_CPU_seconds=.2, child_peak_RSS_KiB=42, child_returncode=-9,
                     child_reaped=True, child_wall_seconds=.8)
        return 'child_cpu_limit'

    def environment(self, worker=None):
        return patch.multiple(a, _run_worker=worker or self.fake_worker)

    def test_improving_interrupted_witness_is_native_checked_and_cpu_charged(self):
        previous = native._EXTERNAL_BUDGET_CHECK.get()
        with self.environment(), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            result = self.call()
        self.assertIs(native._EXTERNAL_BUDGET_CHECK.get(), previous)
        self.assertTrue(result.eligible, result.refusal)
        self.assertEqual(result.order, self.improved)
        self.assertEqual(result.native_promotion, 'native_verified_improvement')
        self.assertEqual(result.score, self.context.score(self.improved, 90))
        self.assertLess(result.score.value, self.score.value)
        self.assertEqual(self.worker_calls, 1)
        self.assertEqual(result.stats['child_process_CPU_seconds'], .2)
        self.assertAlmostEqual(result.CPU_seconds, result.stats['parent_process_CPU_seconds']+.2)
        self.assertFalse(result.stats['fallback'])
        self.assertEqual(result.raw_julia_candidate['native_score']['steps'], result.score.steps)
        request = a.tomllib.loads((self.directory/'checkpoint.toml.request.toml').read_text())
        self.assertEqual(request['units'][0]['coefficients'], ['3/4'])
        self.assertEqual(request['seeds'], a._seeds('tiny:1'))

    def test_no_fresh_output_retry(self):
        with self.environment(), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            first = self.call()
            second = self.call()
        self.assertTrue(first.eligible)
        self.assertFalse(second.eligible)
        self.assertEqual(second.refusal['code'], 'existing_output')
        self.assertEqual(self.worker_calls, 1)

    def test_uncalibrated_and_wrong_objective_refuse_before_child(self):
        with self.environment(), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            result = self.call(settings={})
            wrong_time = self.call(time=.5)
            wrong_budget = self.call(rotation_budget=89)
        self.assertEqual(result.refusal['code'], 'uncalibrated_settings')
        self.assertFalse(wrong_time.eligible)
        self.assertFalse(wrong_budget.eligible)
        self.assertEqual(self.worker_calls, 0)

    def test_controlled_startup_expiry_retains_distinct_parent_fallback(self):
        readings = iter([0., 10.])
        clock = lambda: next(readings, 10.)
        with patch.object(a._time, 'process_time', side_effect=clock), self.environment():
            result = self.call(cpu_seconds=1.)
        self.assertFalse(result.eligible)
        self.assertEqual(result.order, self.parent)
        self.assertEqual(result.score, self.score)
        self.assertTrue(result.stats['fallback'])
        self.assertFalse(result.stats['native_parent_revalidated'])
        self.assertIsNone(result.raw_julia_candidate)
        self.assertEqual(self.worker_calls, 0)

    def test_reply_mismatch_and_incomplete_chains_are_ineligible(self):
        original = self.fake_worker
        def corrupt(*args):
            outcome = original(*args)
            reply = copy.deepcopy(self.last_reply)
            reply['coefficients_sha256'] = '0'*64
            write_toml(args[2], reply)
            return outcome
        with self.environment(corrupt), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            result = self.call()
        self.assertFalse(result.eligible)
        self.assertEqual(result.refusal['code'], 'invalid_worker_reply')
        self.assertEqual(result.order, self.parent)

    def test_valid_witness_survives_method_failure_with_raw_score(self):
        original = self.fake_worker
        def refused(*args):
            original(*args)
            reply = copy.deepcopy(self.last_reply)
            reply['status'] = 'numerical_limit'
            write_toml(args[2], reply)
            return 'child_numerical_limit'
        with self.environment(refused), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            result = self.call()
        self.assertFalse(result.eligible)
        self.assertEqual(result.refusal['code'], 'ineligible_method_execution')
        self.assertEqual(result.order, self.improved)
        self.assertIsNotNone(result.raw_julia_candidate['native_score'])

    def test_missing_julia_is_checked_parent_fallback(self):
        with patch.object(a.shutil, 'which', return_value=None), self.environment():
            result = self.call()
        self.assertFalse(result.eligible)
        self.assertEqual(result.refusal['code'], 'julia_unavailable')
        self.assertTrue(result.stats['native_parent_revalidated'])
        self.assertEqual(result.order, self.parent)
        self.assertEqual(self.worker_calls, 0)

    def test_source_identity_model_and_nonfinite_guards(self):
        with patch.object(native, '__file__', '/stale/install/analysis/sparse_response_ordering.py'):
            with self.assertRaisesRegex(a.Refusal, 'exact checkout'):
                a._native()
        with self.assertRaises(a.Refusal):
            a.validate_settings(dict(settings(), theta_start=math.inf))
        with self.assertRaises(a.Refusal):
            a._score_dict(dict(value=1., cubic_response=1., steps=True))
        with self.assertRaises(a.Refusal):
            a._seeds('x'*257)
        with self.assertRaises(a.Refusal):
            a._compiled_units(native, ({((0, 'Y'),): Fraction(1)},), 1, lambda: None)
        with self.assertRaises(a.Refusal):
            a._compiled_units(native, ({((0, 'X'),): Fraction(0)},), 1, lambda: None)
        with self.assertRaises(a.Refusal):
            a._compiled_units(native, self.blocks, 15, lambda: None)
        result = a.train_annealing(self.blocks, 1, {0: 1.}, parent_order=self.parent,
            starts=self.starts, time=1., rotation_budget=90, parent_score=self.score,
            settings=settings(), seed_material='tiny:1', cpu_seconds=1., wall_seconds=5.)
        self.assertFalse(result.eligible)
        self.assertIn('coherent', result.refusal['message'])

    def test_exact_fraction_wire_survives_float_rounding(self):
        q = Fraction((1 << 54)+1, 1 << 54)
        _, units = a._compiled_units(native, ({((0, 'Z'),): q},), 1, lambda: None)
        self.assertEqual(float(q), 1.)
        self.assertEqual(units[0]['coefficients'], [str(q.numerator)+'/'+str(q.denominator)])

    def test_reply_byte_and_nesting_refusals(self):
        path = self.directory/'bad.toml'
        path.write_bytes(b'x'*(64*1024+1))
        with self.assertRaises(a.Refusal):
            a._bounded_reply(path)
        path.write_text('a=[[[[1]]]]\n')
        with self.assertRaises(a.Refusal):
            a._bounded_reply(path)


class OwnedChildAccountingTests(unittest.TestCase):
    def exercise(self, *, limit):
        stats = {}
        pid = 876543
        process = SimpleNamespace(pid=pid, returncode=None)
        usage = SimpleNamespace(ru_utime=.1, ru_stime=.2, ru_maxrss=42)
        waits = [(0, 0, None), (pid, 9 if limit else 0, usage)]
        seen = []
        with patch.object(a.subprocess, 'Popen', return_value=process) as launch, \
             patch.object(a.os, 'wait4', side_effect=waits) as wait, \
             patch.object(a.os, 'killpg') as kill, \
             patch.object(a, '_proc_cpu', return_value=.2 if limit else .02), \
             patch.object(a._time, 'sleep'), patch.object(a._time, 'monotonic', return_value=0.):
            result = a._run_worker('/fixture/julia', Path('/fixture/request'), Path('/fixture/reply'),
                .1 if limit else 1., 5., lambda cpu: seen.append(cpu), stats)
        self.assertEqual(wait.call_args_list[0].args, (pid, a.os.WNOHANG))
        self.assertEqual(wait.call_args_list[-1].args, (pid, 0 if limit else a.os.WNOHANG))
        self.assertEqual(kill.call_count, int(limit))
        self.assertTrue(stats['child_reaped'])
        self.assertAlmostEqual(stats['child_process_CPU_seconds'], .3)
        argv = launch.call_args.args[0]
        self.assertIn('--startup-file=no', argv)
        self.assertIn('--threads=1', argv)
        self.assertEqual(launch.call_args.kwargs['env']['JULIA_PKG_PRECOMPILE_AUTO'], '0')
        self.assertEqual(result, 'child_cpu_limit' if limit else 'complete')
        self.assertEqual(len(seen), 0 if limit else 1)

    def test_exact_owned_wait4_user_system_and_environment(self):
        self.exercise(limit=False)

    def test_cpu_expiry_kills_and_reaps_only_owned_process_group(self):
        self.exercise(limit=True)


if __name__ == '__main__':
    unittest.main()
