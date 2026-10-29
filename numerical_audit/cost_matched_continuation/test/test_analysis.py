"""Small identity/arithmetic fixtures; no numerical optimizer or molecule jobs."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('continuation_analysis_under_test', HERE/'analyze.py')
a = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(a)


def costs(cpu, wall=None):
    return dict(raw_CPU_seconds=cpu, charged_CPU_seconds=cpu, wall_seconds=4*cpu if wall is None else wall)


def record(arm, budget=0, prefix=1., continuation=0.):
    return dict(case='A', repetition=1, arm=arm, total_CPU_budget=budget, analysis_eligible=True, guard_exit_code=0,
                costs=dict(prefix=costs(prefix), continuation=costs(continuation), total=costs(prefix+continuation)))


def write_json(root, name, value):
    path = root/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, allow_nan=False))
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalysisArithmeticTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((HERE/'PROTOCOL.json').read_text())
        self.cases = ['A', 'B']
        self.physical, self.constructors = {}, {}
        for case in self.cases:
            for rep in (1, 2):
                for arm in ('full', 'window'):
                    for budget in (32, 64):
                        self.constructors[case, rep, arm, budget] = dict(analysis_eligible=True)
                        for duration in (.65, 1.35):
                            for cap in (8, 24):
                                self.physical[case, arm, rep, budget, duration, cap] = dict(
                                    projective_error=.02 if arm == 'full' else .01)
        self.verified = set(self.physical)

    def pairs(self):
        return a.pair_grid(self.cases, 2, self.config, self.physical, self.constructors, self.verified)

    def test_complete_pair_grid_and_primary_rule(self):
        pairs = self.pairs()
        self.assertEqual(len(pairs), 32)
        self.assertTrue(all(row['outcome'] == 'win' for row in pairs))
        primary = a.primary_decision(pairs, self.cases, 2, self.config, True)
        self.assertEqual(primary['decision'], 'window_advantage')
        self.assertEqual(primary['per_case_median_paired_difference'], dict(A=-.01, B=-.01))
        self.assertEqual(len(primary['rows']), 4)
        self.assertIsNone(primary['confidence_interval'])

    def test_primary_repetition_loss_vetoes_negative_case_median(self):
        self.physical['A', 'window', 1, 64, 1.35, 8]['projective_error'] = .021
        primary = a.primary_decision(self.pairs(), self.cases, 2, self.config, True)
        self.assertLess(primary['per_case_median_paired_difference']['A'], -1e-12)
        self.assertEqual(primary['decision'], 'mixed_or_unresolved')

    def test_symmetric_full_advantage(self):
        for key, row in self.physical.items():
            row['projective_error'] = .01 if key[1] == 'full' else .02
        self.assertEqual(a.primary_decision(self.pairs(), self.cases, 2, self.config, True)['decision'], 'full_advantage')

    def test_unverified_refusal_retains_loss_and_prevents_advantage(self):
        key = ('A', 'window', 1, 64, 1.35, 8)
        self.physical[key]['projective_error'] = .03
        self.constructors['A', 1, 'window', 64]['analysis_eligible'] = False
        self.verified.remove(key)
        pairs = self.pairs()
        failed = next(row for row in pairs if (row['case'], row['repetition'], row['total_CPU_budget'],
            row['time'], row['rotation_cap_multiplier']) == ('A', 1, 64, 1.35, 8))
        self.assertEqual(failed['outcome'], 'loss')
        self.assertAlmostEqual(failed['difference'], .01)
        self.assertFalse(failed['eligible'])
        self.assertEqual(a.primary_decision(pairs, self.cases, 2, self.config, False)['decision'], 'study_incomplete')

    def test_missing_pair_stays_on_grid_and_cannot_be_selected_away(self):
        del self.physical['A', 'window', 1, 64, 1.35, 8]
        pairs = self.pairs()
        self.assertEqual(len(pairs), 32)
        unavailable = [row for row in pairs if row['outcome'] == 'unavailable']
        self.assertEqual(len(unavailable), 1)
        self.assertIsNone(unavailable[0]['difference'])
        self.assertEqual(a.primary_decision(pairs, self.cases, 2, self.config, True)['decision'], 'study_incomplete')

    def test_floor_suppresses_ratio_without_suppressing_paired_error(self):
        for arm, value in [('full', 0.), ('window', 5e-13)]:
            self.physical['A', arm, 1, 64, 1.35, 8]['projective_error'] = value
        row = next(row for row in self.pairs() if row['case'] == 'A' and row['repetition'] == 1 and
            row['total_CPU_budget'] == 64 and row['time'] == 1.35 and row['rotation_cap_multiplier'] == 8)
        self.assertEqual(row['difference'], 5e-13)
        self.assertEqual(row['outcome'], 'tie')
        self.assertTrue(row['full_at_or_below_floor'] and row['window_at_or_below_floor'])
        self.assertIsNone(row['window_to_full_error_ratio'])

    def test_whole_budget_and_prefix_reconciliation(self):
        parent, child = record('parent'), record('full', 32, continuation=31.)
        a.cost_check(parent, parent, self.config)
        a.cost_check(child, parent, self.config)
        over = copy.deepcopy(child)
        over['costs']['continuation']['charged_CPU_seconds'] += 1e-10
        over['costs']['total']['charged_CPU_seconds'] += 1e-10
        with self.assertRaisesRegex(ValueError, 'Strict whole-construction CPU'):
            a.cost_check(over, parent, self.config)
        wrong = copy.deepcopy(child)
        wrong['costs']['prefix']['wall_seconds'] += .1
        with self.assertRaisesRegex(ValueError, 'Shared parent prefix'):
            a.cost_check(wrong, parent, self.config)
        wall = record('window', 32, continuation=1.)
        wall['costs']['continuation']['wall_seconds'] = 157.
        wall['costs']['total']['wall_seconds'] = 161.
        with self.assertRaisesRegex(ValueError, 'wall ceiling'):
            a.cost_check(wall, parent, self.config)

    def test_shared_parent_charges_are_not_multiplied_by_physical_rows(self):
        records = {('A', 1, 'parent', 0): record('parent', prefix=3.)}
        for arm in ('full', 'window'):
            for budget in (32, 64):
                records['A', 1, arm, budget] = record(arm, budget, prefix=3., continuation=10.)
        result = a.resource_totals(records)['charged_CPU_seconds']
        self.assertEqual(result['actual_shared_constructor_compute'], 43.)
        self.assertEqual(result['attributed_continuation_comparisons'], 52.)
        self.assertEqual(result['actual_parent_compute'], 3.)

    def test_duplicate_join_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate identity'):
            a.unique([dict(name='same'), dict(name='same')], lambda row: row['name'])

    def test_unknown_raw_cpu_retains_known_lower_bound_and_conservative_totals(self):
        parent = record('parent', prefix=3.)
        parent['costs']['total']['raw_CPU_seconds'] = None
        parent['costs']['prefix']['raw_CPU_seconds'] = None
        records = {('A', 1, 'parent', 0): parent}
        for arm in ('full', 'window'):
            for budget in (32, 64):
                row = record(arm, budget, prefix=3., continuation=10.)
                row['costs']['prefix']['raw_CPU_seconds'] = None
                row['costs']['total']['raw_CPU_seconds'] = None
                records['A', 1, arm, budget] = row
        summary = a.resource_totals(records)
        raw = summary['raw_CPU_seconds']
        self.assertIsNone(raw['actual_shared_constructor_compute'])
        self.assertEqual(raw['known_actual_shared_constructor_compute'], 40.)
        self.assertEqual(raw['actual_shared_constructor_compute_unknown_components'], 1)
        self.assertEqual(raw['attributed_continuation_comparisons_unknown_components'], 4)
        self.assertEqual(summary['charged_CPU_seconds']['actual_shared_constructor_compute'], 43.)


class ArchivedEvidenceTests(unittest.TestCase):
    def worker_fixture(self, root, arm='full', budget=32):
        score = dict(value=1e-9, cubic_response=.001, steps=15)
        result = dict(order=[0, 1, 2, 3, 4], score=score, accepted_moves=0, status='cpu_limit')
        identity = dict(case='A', repetition=1, arm=arm, total_CPU_budget=budget,
            model_json_sha256='a'*64, model_toml_sha256='b'*64, protocol_sha256='c'*64,
            source_freeze_sha256='d'*64, parent_sha256=None if arm == 'parent' else 'e'*64)
        worker = dict(identity, schema='cost-matched-worker-v1', result=copy.deepcopy(result), refusal=None)
        path = 'workers/'+arm+str(budget)+'.json'
        digest = write_json(root, path, worker)
        constructor = dict(identity, controller_fallback=False, worker_output_path=path,
            worker_output_sha256=digest, result=copy.deepcopy(result),
            order=list(result['order']), score=copy.deepcopy(score))
        return constructor, worker

    def test_worker_file_binds_native_result_and_all_arm_identities(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for arm, budget in [('parent', 0), ('full', 32), ('window', 64)]:
                with self.subTest(arm=arm):
                    constructor, _ = self.worker_fixture(root, arm, budget)
                    audit = a.Audit(root)
                    self.assertTrue(audit.worker(constructor)['verified'])
                    self.assertEqual(audit.inputs[constructor['worker_output_path']], constructor['worker_output_sha256'])

    def test_worker_header_identity_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructor, _ = self.worker_fixture(root)
            replacements = dict(case='B', repetition=2, arm='window', total_CPU_budget=64,
                model_json_sha256='0'*64, model_toml_sha256='1'*64, protocol_sha256='2'*64,
                source_freeze_sha256='3'*64, parent_sha256='4'*64)
            for field, value in replacements.items():
                with self.subTest(field=field):
                    changed = copy.deepcopy(constructor)
                    changed[field] = value
                    with self.assertRaisesRegex(ValueError, 'Worker/canonical identity'):
                        a.Audit(root).worker(changed)

    def test_worker_selected_order_score_and_original_result_cannot_change(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructor, _ = self.worker_fixture(root)
            order = copy.deepcopy(constructor)
            order['order'] = [1, 0, 2, 3, 4]
            score = copy.deepcopy(constructor)
            score['score']['value'] *= 2
            native = copy.deepcopy(constructor)
            native['result']['accepted_moves'] = 99
            for label, changed in [('order', order), ('score', score), ('native', native)]:
                with self.subTest(change=label):
                    with self.assertRaises(ValueError):
                        a.Audit(root).worker(changed)

    def test_worker_digest_mismatch_is_rejected_before_metadata_join(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructor, worker = self.worker_fixture(root)
            worker['result']['score']['value'] *= 2
            write_json(root, constructor['worker_output_path'], worker)
            with self.assertRaisesRegex(ValueError, 'Input digest differs'):
                a.Audit(root).worker(constructor)

    def test_missing_unscored_and_refused_workers_keep_fallback_reasons(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructor, worker = self.worker_fixture(root)
            constructor['controller_fallback'] = True
            worker['result']['score'] = None
            constructor['result'] = copy.deepcopy(worker['result'])
            constructor['worker_output_sha256'] = write_json(root, constructor['worker_output_path'], worker)
            unscored = a.Audit(root).worker(constructor)
            self.assertFalse(unscored['verified'])
            self.assertEqual(unscored['reason'], 'unscored_worker_result')
            worker['result'] = constructor['result'] = None
            worker['refusal'] = dict(type='SearchDeadline', message='retained worker refusal')
            constructor['worker_output_sha256'] = write_json(root, constructor['worker_output_path'], worker)
            refused = a.Audit(root).worker(constructor)
            self.assertFalse(refused['verified'])
            self.assertEqual(refused['reason'], 'worker_refusal')
            self.assertEqual(refused['refusal'], worker['refusal'])
            (root/constructor['worker_output_path']).unlink()
            constructor['worker_output_sha256'] = None
            missing = a.Audit(root).worker(constructor)
            self.assertFalse(missing['verified'])
            self.assertEqual(missing['reason'], 'missing_worker_output')
            constructor['controller_fallback'] = False
            with self.assertRaisesRegex(ValueError, 'Missing worker digest'):
                a.Audit(root).worker(constructor)

    def test_digest_failure_keeps_actual_digest_and_never_uses_current_source(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            digest = write_json(root, 'artifact.json', dict(retained='yes'))
            audit = a.Audit(root)
            with self.assertRaisesRegex(ValueError, 'digest differs'):
                audit.load('artifact.json', '0'*64)
            self.assertEqual(audit.inputs['artifact.json'], digest)
            with self.assertRaisesRegex(ValueError, 'outside retained run'):
                audit.load('../outside.json')
            copied = root/'source_snapshot/old/dependency.py'
            copied.parent.mkdir(parents=True)
            copied.write_text('archived = True\n')
            frozen_hash = hashlib.sha256(copied.read_bytes()).hexdigest()
            audit.source_hashes({'old/dependency.py': frozen_hash}, {'sources': {'old/dependency.py': frozen_hash}})

    def test_new_output_preserves_existing_file_and_partial(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'ANALYSIS.json'
            a.write_new(path, dict(value=1))
            with self.assertRaisesRegex(ValueError, 'Preserve existing'):
                a.write_new(path, dict(value=2))
            self.assertEqual(json.loads(path.read_text()), dict(value=1))
            other = Path(folder)/'OTHER.json'
            other.with_name(other.name+'.partial').write_text('retained')
            with self.assertRaisesRegex(ValueError, 'Preserve existing'):
                a.write_new(other, dict(value=3))

    def ledger_fixture(self, root, extra_failure=False):
        guard = root/'campaign/resource_guard.py'
        guard.parent.mkdir(parents=True, exist_ok=True)
        guard.write_text('# Archived guard fixture; not executed.\n')
        guard_sha = hashlib.sha256(guard.read_bytes()).hexdigest()
        dispatcher = root/'campaign/run_guarded_job.py'
        dispatcher.write_text('# Archived dispatcher fixture; not executed.\n')
        dispatcher_sha = hashlib.sha256(dispatcher.read_bytes()).hexdigest()
        pools_sha = write_json(root, 'campaign/RESOURCE_POOLS.json', dict(
            phases={'constructors': dict(CPU_seconds=100., wall_seconds=500.)},
            guard_path='campaign/resource_guard.py', guard_sha256=guard_sha))
        settings = dict(CPUQuota='25%', Nice='19', CPUWeight='5', IOWeight='10',
                        MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0')
        jobs, references = [], {}
        labels = ['parent', 'full32', 'full64', 'window32', 'window64']
        if extra_failure:
            labels.append('failed-attempt')
        for label in labels:
            cpu = 1. if label == 'parent' else 10.
            reference = dict(label=label, final_path='campaign/jobs/'+label+'.FINAL.json',
                receipt_path='campaign/logs/'+label+'.resources.json')
            reference['receipt_sha256'] = write_json(root, reference['receipt_path'], dict(settings=settings,
                exit_code=2 if label == 'failed-attempt' else 0, elapsed=cpu*4, final_accounting=dict(CPU_seconds=cpu)))
            reference['final_sha256'] = write_json(root, reference['final_path'], dict(label=label, phase='constructors',
                receipt_sha256=reference['receipt_sha256'], charged_CPU_seconds=cpu,
                absolute_guarded_wall_seconds=cpu*4, exit_code=2 if label == 'failed-attempt' else 0))
            plan_path = 'campaign/jobs/'+label+'.PLAN.json'
            plan_sha = write_json(root, plan_path, dict(label=label, phase='constructors', pools_sha256=pools_sha,
                dispatcher_sha256=dispatcher_sha, attempt_CPU_seconds=cpu+1, attempt_wall_seconds=cpu*4+5))
            jobs.append(dict(reference, phase='constructors', plan_path=plan_path, plan_sha256=plan_sha))
            references[label] = dict(costs(cpu), **reference)
        write_json(root, 'campaign/LEDGER_CUTOFF.json', dict(schema='cost-matched-ledger-v1', cutoff_UTC='2026-10-29T00:00:00Z',
            scope='Scientific jobs through cutoff; later analysis/reporting/review/packaging excluded',
            resource_pools=dict(path='campaign/RESOURCE_POOLS.json', sha256=pools_sha), jobs=jobs))
        constructors = {('A', 1, 'parent', 0): record('parent')}
        constructors['A', 1, 'parent', 0]['costs']['prefix'] = references['parent']
        for arm in ('full', 'window'):
            for budget in (32, 64):
                row = record(arm, budget, continuation=10.)
                row['costs']['prefix'] = references['parent']
                row['costs']['continuation'] = references[arm+str(budget)]
                constructors['A', 1, arm, budget] = row
        return constructors

    def test_ledger_receipts_in_logs_and_extra_receipt_refusal(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root)
            audit = a.Audit(root)
            result = audit.ledger(constructors)
            self.assertEqual(audit.issues, [])
            self.assertEqual(result['reconciled_jobs'], 5)
            write_json(root, 'campaign/logs/unindexed.resources.json', dict(retained='yes'))
            audit = a.Audit(root)
            audit.ledger(constructors)
            self.assertTrue(any('Unindexed/missing archived job' in issue['message'] for issue in audit.issues))

    def test_ledger_joins_once_per_actual_job_and_retains_failures(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root, extra_failure=True)
            audit = a.Audit(root)
            result = audit.ledger(constructors)
            self.assertEqual(result['declared_jobs'], 6)
            self.assertEqual(result['reconciled_jobs'], 6)
            self.assertEqual(result['charged_CPU_seconds'], 51.)
            self.assertEqual(len(result['phases']['constructors']['jobs']), 6)
            self.assertTrue(any('unattributed attempts' in issue['message'] for issue in audit.issues))
            self.assertTrue(result['later_analysis_reporting_review_packaging_excluded'])

    def test_ledger_cost_mutation_is_not_trusted(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root)
            constructors['A', 1, 'full', 32]['costs']['continuation']['charged_CPU_seconds'] = 9.
            audit = a.Audit(root)
            audit.ledger(constructors)
            self.assertTrue(any('differ from OS evidence' in issue['message'] for issue in audit.issues))

    def test_rejected_overshoot_remains_in_measured_phase_totals(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root)
            cutoff = json.loads((root/'campaign/LEDGER_CUTOFF.json').read_text())
            entry = cutoff['jobs'][0]
            final = json.loads((root/entry['final_path']).read_text())
            final['charged_CPU_seconds'] = 3.  # Parent reserved only two CPU seconds.
            entry['final_sha256'] = write_json(root, entry['final_path'], final)
            write_json(root, 'campaign/LEDGER_CUTOFF.json', cutoff)
            audit = a.Audit(root)
            result = audit.ledger(constructors)
            self.assertEqual(result['charged_CPU_seconds'], 43.)
            self.assertEqual(len(result['jobs']), 5)
            self.assertFalse(result['jobs'][0]['reconciled'])
            self.assertTrue(any('exceeds reserved' in issue['message'] for issue in audit.issues))

    def test_dispatcher_bytes_are_bound_to_each_plan_without_dropping_charges(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root)
            (root/'campaign/run_guarded_job.py').write_text('# Replaced dispatcher.\n')
            audit = a.Audit(root)
            result = audit.ledger(constructors)
            self.assertEqual(result['charged_CPU_seconds'], 41.)
            self.assertEqual(result['declared_jobs'], 5)
            self.assertEqual(result['reconciled_jobs'], 0)
            self.assertTrue(any('Archived source digest differs: campaign/run_guarded_job.py' in issue['message']
                                for issue in audit.issues))


if __name__ == '__main__':
    unittest.main()
