"""Old-only calibration contracts; synthetic witnesses are not study results.

These tests do not spawn Julia, run chemistry, or evaluate physical errors.
The actual capped calibration CLI independently scores every fixed native/full/
archived witness; these fixtures test fitting and refusal/accounting boundaries.
"""
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('molecular_calibration_contract_test', HERE/'calibration.py')
c = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = c
spec.loader.exec_module(c)


class CalibrationContracts(unittest.TestCase):
    def test_retained_old_inputs_preserve_legacy_identity_and_rationals(self):
        for declared in c.CASES:
            model = c.read_model(c.ROOT/declared['path'], declared)
            self.assertEqual(model['data']['source_model_sha256'], declared['legacy_source_model_sha256'])
            self.assertNotIn('source_json_sha256', model['data'])
            self.assertEqual(model['words'], declared['words'])
            self.assertEqual(len(model['blocks']), declared['units'])
            self.assertEqual(model['reference'], 63)
            self.assertEqual(hashlib.sha256(model['raw']).hexdigest(), declared['sha256'])
            self.assertTrue(all(isinstance(value, Fraction) for block in model['blocks'] for value in block.values()))
            wrong = dict(declared, legacy_source_model_sha256='0'*64)
            with self.assertRaisesRegex(c.Refusal, 'legacy'):
                c.read_model(c.ROOT/declared['path'], wrong)

    def test_exactly_32_predeclared_full_insertions_with_duplicate_observations(self):
        parent = [0, 1, 2, 3, 4]
        plan = c.probe_plan('CO_r1p13', parent)
        self.assertEqual(plan, c.probe_plan('CO_r1p13', parent))
        self.assertEqual(len(plan['pairs']), 32)
        self.assertEqual(len(plan['orders']), 33)
        material = 'CO_r1p13:molecular-generalization-sa-v1:calibration-v1'
        expected_seed = int.from_bytes(hashlib.sha256(material.encode()).digest()[:8], 'big')
        self.assertEqual(plan['seed'], str(expected_seed))
        for pair, order in zip(plan['pairs'], plan['orders'][1:]):
            source, destination = pair
            self.assertNotEqual(source, destination)
            trial = list(parent);trial.insert(destination, trial.pop(source))
            self.assertEqual(order, trial)
            inverse = list(order);inverse.insert(source, inverse.pop(destination))
            self.assertEqual(inverse, parent)
        # Every draw has the same changed order on two labels. None disappear.
        duplicate = c.probe_plan('CO_r1p13', [0, 1])
        self.assertEqual(len(duplicate['orders'][1:]), 32)
        self.assertEqual(duplicate['orders'][1:], [[1, 0]]*32)
        with self.assertRaises(c.Refusal):
            c.probe_plan('new_holdout_family', parent)

    def test_parity_tolerance_and_floors_keep_missing_evidence_ineligible(self):
        self.assertTrue(c.comparison(1., 1.+5e-9)['verified'])
        self.assertFalse(c.comparison(1., 1.+5e-7)['verified'])
        floor = c.comparison(1e-25, 0.)
        self.assertTrue(floor['within_tolerance'])
        self.assertTrue(floor['numerical_floor'])
        self.assertFalse(floor['verified'])
        for bad in (math.inf, math.nan, -1., True):
            with self.assertRaises(c.Refusal):
                c.comparison(bad, 1.)

    def test_temperature_fit_uses_every_probe_but_only_positive_native_deltas(self):
        rows = [dict(parent_J=2., native_J=value, verified=True) for value in (4., 8., 1., 2.)]
        before = copy.deepcopy(rows)
        result = c.fit_temperature(rows, expected_count=4)
        self.assertEqual(rows, before)
        self.assertEqual(result['positive_delta_count'], 2)
        self.assertEqual(result['pooled_positive_median'], 2.)
        self.assertAlmostEqual(math.exp(-2./result['theta_start']), .8)
        self.assertAlmostEqual(math.exp(-2./result['theta_end']), .01)
        censored = copy.deepcopy(rows);censored[2]['verified'] = False
        with self.assertRaises(c.Refusal):
            c.fit_temperature(censored, expected_count=4)
        with self.assertRaises(c.Refusal):
            c.fit_temperature(rows[:3], expected_count=4)
        with self.assertRaisesRegex(c.Refusal, 'Insufficient'):
            c.fit_temperature([dict(parent_J=2., native_J=1., verified=True)]*4, expected_count=4)

    def test_full_fit_rejects_duplicate_model_or_probe_identity(self):
        rows = [dict(case=declared['case'], index=index, parent_J=1., native_J=2., verified=True)
                for declared in c.CASES for index in range(1, 33)]
        self.assertEqual(c.fit_temperature(rows)['positive_delta_count'], 64)
        forged = copy.deepcopy(rows);forged[-1] = copy.deepcopy(forged[0])
        with self.assertRaisesRegex(c.Refusal, 'roster'):
            c.fit_temperature(forged)
        forged = copy.deepcopy(rows);forged[1]['parent_J'] = 1.1
        with self.assertRaisesRegex(c.Refusal, 'normalization'):
            c.fit_temperature(forged)

    def row_fixture(self, path, *, full_J=1., archived_J=1., model='a'*64):
        legacy = c.CASES[0]['legacy_source_model_sha256']
        path.write_text('\n'.join([
            'version=1', 'case="CO_r1p13"', 'index=1', 'order=[0,1]', 'status="complete"',
            'legacy_source_model_sha256='+json.dumps(legacy), 'model_sha256='+json.dumps(model),
            'request_sha256="'+'b'*64+'"', 'julia_version="fixture_not_a_runtime_result"',
            '[native_claim]', 'value=1.0', 'cubic_response=16.0', 'steps=2',
            '[full_order_score]', 'value='+repr(full_J), 'cubic_response='+repr(16*full_J), 'steps=2',
            '[archived_score]', 'value='+repr(archived_J), 'cubic_response='+repr(16*archived_J), 'steps=2',
            '[source_sha256]', '"source.py"="'+'c'*64+'"', '']), encoding='ascii')

    def admit(self, path):
        return c._admit_child_row(path, 'CO_r1p13', 1, [0, 1],
            dict(score=dict(value=1., cubic_response=16., steps=2)), 'a'*64, 'b'*64, {'source.py': 'c'*64})

    def test_every_probe_checks_native_new_and_archived_scores_and_source_joins(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'row001.toml'
            self.row_fixture(path)
            self.assertTrue(self.admit(path)['verified'])
            self.row_fixture(path, archived_J=1.01)
            bad = self.admit(path)
            self.assertFalse(bad['verified'])
            self.assertEqual(bad['independent_python_parity']['archived_score']['value']['status'], 'mismatch')
            self.row_fixture(path, full_J=1.+1e-6)
            self.assertFalse(self.admit(path)['verified'])
            self.row_fixture(path, model='0'*64)
            with self.assertRaisesRegex(c.Refusal, 'join'):
                self.admit(path)
            self.row_fixture(path)
            path.write_text(path.read_text().replace('index=1', 'index=true'))
            with self.assertRaises(c.Refusal):
                self.admit(path)

    def test_floor_censor_is_preserved_not_turned_into_a_verified_probe(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'floor.toml'
            self.row_fixture(path, full_J=0., archived_J=0.)
            raw = path.read_text().replace('value=1.0', 'value=0.0').replace('cubic_response=16.0', 'cubic_response=0.0')
            path.write_text(raw)
            result = c._admit_child_row(path, 'CO_r1p13', 1, [0, 1],
                dict(score=dict(value=0., cubic_response=0., steps=2)), 'a'*64, 'b'*64, {'source.py': 'c'*64})
            self.assertFalse(result['verified'])
            self.assertTrue(result['independent_python_parity']['archived_score']['value']['numerical_floor'])

    def test_timeout_inventory_keeps_all_declared_probes_and_native_losses(self):
        spec = c.CASES[0]
        plan = c.probe_plan(spec['case'], [0, 1])
        declared = dict(spec, pairs=plan['pairs'], orders=plan['orders'])
        report = dict(probes=[], cases=[dict(case=spec['case'], native_rows=[
            dict(index=0, score=dict(value=2.)), dict(index=1, score=dict(value=3.)),
            dict(index=2, score=dict(value=1.))])])
        c.retain_probe_inventory(report, [declared])
        self.assertEqual(len(report['probes']), 32)
        self.assertEqual(report['probes'][0]['normalized_delta'], .5)
        self.assertEqual(report['probes'][1]['normalized_delta'], -.5)
        self.assertIsNone(report['probes'][2]['native_J'])
        self.assertTrue(all(row['verified'] is False for row in report['probes']))
        self.assertEqual(report['probes'][2]['status'], 'unprocessed')

    def test_fresh_output_and_refused_setup_leave_no_fabricated_settings(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)/'new'
            with patch.object(c, 'adapter', side_effect=c.Refusal('injected setup refusal')), \
                 patch.object(c.subprocess, 'Popen') as launch:
                report = c.run_calibration(output, 100., 500.)
            self.assertFalse(report['complete'])
            self.assertFalse((output/'ANNEALING_SETTINGS.json').exists())
            self.assertFalse(report['physical_errors_read'])
            self.assertEqual(json.loads((output/'CALIBRATION.json').read_text())['status'], 'incomplete')
            launch.assert_not_called()
            with self.assertRaisesRegex(c.Refusal, 'Fresh'):
                c.run_calibration(output, 100., 500.)
            with self.assertRaises(c.Refusal):
                c.run_calibration(Path(temporary)/'bad', math.inf, 500.)

    def test_exact_owned_child_cost_includes_fresh_runtime_start_and_teardown(self):
        with tempfile.TemporaryDirectory() as temporary:
            pid = 876543
            usage = SimpleNamespace(ru_utime=.2, ru_stime=.1, ru_maxrss=42)
            process = SimpleNamespace(pid=pid, returncode=None)
            bridge = SimpleNamespace(_proc_cpu=lambda owned: .01)
            waits = [(0, 0, None), (pid, 0, usage)]
            with patch.object(c.subprocess, 'Popen', return_value=process) as launch, \
                 patch.object(c.os, 'wait4', side_effect=waits) as wait, \
                 patch.object(c._time, 'sleep'), patch.object(c._time, 'monotonic', return_value=0.):
                result = c.run_child('/fixture/julia', Path('/fixture/model'), Path('/fixture/request'),
                    Path(temporary)/'fresh', 10., 50., lambda live: None, bridge)
            self.assertAlmostEqual(result['raw_CPU_seconds'], .3)
            self.assertTrue(result['reaped'])
            self.assertEqual(result['stop'], 'complete')
            self.assertTrue(all(call.args[0] == pid for call in wait.call_args_list))
            argv = launch.call_args.args[0]
            self.assertIn('--startup-file=no', argv)
            self.assertIn('--threads=1', argv)
            self.assertIn(str(c.HERE/'calibration.jl'), argv)
            self.assertEqual(launch.call_args.kwargs['env']['JULIA_PKG_PRECOMPILE_AUTO'], '0')


if __name__=='__main__':
    unittest.main()
