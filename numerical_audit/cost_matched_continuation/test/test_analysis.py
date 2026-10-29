"""Small identity/arithmetic fixtures; no numerical optimizer or molecule jobs."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

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

    def ledger_fixture(self, root, extra_failure=False, bound=False):
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
            command = ['python', 'fixture.py']
            cap, wall_cap = cpu+1, cpu*4+5
            if bound:
                arm = 'parent' if label == 'parent' else 'full' if label.startswith('full') else 'window'
                budget = 0 if arm == 'parent' else int(label[-2:])
                cap, wall_cap = (8., 40) if arm == 'parent' else (budget-1., 5*budget-4)
                command = ['python', '/source/numerical_audit/cost_matched_continuation/construct.py', '--worker',
                    '--root', '/original/data', '--case', 'A', '--repetition', '1', '--arm', arm,
                    '--budget', str(budget), '--raw', '/original/data/workers/'+label+'.json',
                    '--cpu-cap', str(cap), '--wall-cap', str(wall_cap)]
            reference = dict(label=label, final_path='campaign/jobs/'+label+'.FINAL.json',
                receipt_path='campaign/logs/'+label+'.resources.json')
            reference['receipt_sha256'] = write_json(root, reference['receipt_path'], dict(settings=settings,
                command=command, exit_code=2 if label == 'failed-attempt' else 0, elapsed=cpu*4,
                budget=dict(cpu_seconds=cap, wall_seconds=wall_cap), final_accounting=dict(CPU_seconds=cpu)))
            reference['final_sha256'] = write_json(root, reference['final_path'], dict(label=label, phase='constructors',
                receipt_sha256=reference['receipt_sha256'], charged_CPU_seconds=cpu,
                final_service_query_exit=0, final_service_snapshot=dict(CPUUsageNSec=str(int(cpu*1e9)), MemoryPeak='1048576'),
                absolute_guarded_wall_seconds=cpu*4, exit_code=2 if label == 'failed-attempt' else 0))
            plan_path = 'campaign/jobs/'+label+'.PLAN.json'
            plan_sha = write_json(root, plan_path, dict(label=label, phase='constructors', pools_sha256=pools_sha,
                command=command, cwd='/source', dispatcher_sha256=dispatcher_sha,
                attempt_CPU_seconds=cap, attempt_wall_seconds=wall_cap))
            log_path = 'campaign/logs/'+label+'.log'
            (root/log_path).write_text('Memory peak: 1.0M\n')
            jobs.append(dict(reference, phase='constructors', plan_path=plan_path, plan_sha256=plan_sha,
                log_path=log_path, log_sha256=hashlib.sha256((root/log_path).read_bytes()).hexdigest()))
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
        if bound:
            for row in constructors.values():
                label = 'parent' if row['arm'] == 'parent' else row['arm']+str(row['total_CPU_budget'])
                cap, wall = (8., 40) if label == 'parent' else (row['total_CPU_budget']-1., 5*row['total_CPU_budget']-4)
                row['worker_output_path'] = 'workers/'+label+'.json'
                row['worker_verification'] = dict(worker_data=dict(guarded_CPU_cap=cap, guarded_wall_cap=wall,
                    process_CPU_before_output=.9 if label == 'parent' else 9.99))
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

    def test_swapped_full_window_receipts_reject_owners_but_retain_all_charges(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            constructors = self.ledger_fixture(root, bound=True)
            full, window = constructors['A', 1, 'full', 32], constructors['A', 1, 'window', 32]
            full['costs']['continuation'], window['costs']['continuation'] = window['costs']['continuation'], full['costs']['continuation']
            audit = a.Audit(root)
            result = audit.ledger(constructors, json.loads((HERE/'PROTOCOL.json').read_text()), 'preflight', ['A'])
            self.assertEqual(result['charged_CPU_seconds'], 41.)
            self.assertEqual(result['reconciled_jobs'], 5)
            self.assertFalse(full['analysis_eligible'])
            self.assertFalse(window['analysis_eligible'])
            self.assertTrue(any('PLAN constructor identity differs' in issue['message'] for issue in audit.issues))


class ConsumerBoundaryTests(unittest.TestCase):
    def graph_fixture(self):
        order = [0, 1, 2, 3, 4]
        score = dict(value=1e-9, cubic_response=1e-5, steps=10)
        restart = dict(seed=int.from_bytes(hashlib.sha256(b'A:1:cost-matched-v1').digest()[:8], 'big'),
            interior_perturbation_counts=[0, 1, 2, 4], orders=[order, [0, 2, 1, 3, 4], [0, 3, 1, 2, 4], [0, 2, 3, 1, 4]])
        result = dict(order=order, score=score, status='complete', stats=dict(seed_evaluations=3))
        worker = dict(result=result, refusal=None, starts=None, restarts=copy.deepcopy(restart))
        parent = dict(case='A', repetition=1, arm='parent', total_CPU_budget=0, eligible=True,
            controller_fallback=False, guard_exit_code=0, order=order, score=score, restarts=restart,
            worker_verification=dict(verified=True, worker_data=worker))
        return parent

    def test_finite_unique_json_retains_digest_on_parse_rejection(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for raw in (b'{"key":1,"key":2}', b'{"key":NaN}', b'{"key":Infinity}', b'{"key":1e999}'):
                with self.subTest(raw=raw):
                    (root/'bad.json').write_bytes(raw)
                    audit = a.Audit(root)
                    with self.assertRaises(ValueError):
                        audit.load('bad.json')
                    self.assertEqual(audit.inputs['bad.json'], hashlib.sha256(raw).hexdigest())
            with self.assertRaisesRegex(ValueError, 'nesting'):
                a.parse_json(b'['*65+b'0'+b']'*65)

    def test_ineligible_unknown_cpu_and_overshoot_are_retained(self):
        config = json.loads((HERE/'PROTOCOL.json').read_text())
        parent, child = record('parent'), record('window', 32, continuation=40.)
        child['eligible'] = False
        child['costs']['continuation']['raw_CPU_seconds'] = None
        child['costs']['total']['raw_CPU_seconds'] = None
        a.cost_check(child, parent, config)
        summary = a.resource_totals({('A', 1, 'parent', 0): parent, ('A', 1, 'window', 32): child})
        self.assertIsNone(summary['raw_CPU_seconds']['actual_shared_constructor_compute'])
        self.assertEqual(summary['charged_CPU_seconds']['actual_shared_constructor_compute'], 41.)

    def test_terminal_cpu_rejects_poll_only_measurement_and_preserves_coarse_clock(self):
        receipt = dict(final_accounting=dict(CPU_seconds=.01, process_CPU_timing_available=False),
                       samples=[dict(service='CPUUsageNSec=10000000\n')])
        unknown = a.terminal_accounting(receipt, {}, '', 32.)
        self.assertFalse(unknown['terminal_CPU_measurement_available'])
        self.assertIsNone(unknown['raw_terminal_CPU_seconds'])
        self.assertEqual(unknown['minimum_charged_CPU_seconds'], 32.)
        receipt['final_accounting']['process_CPU_timing_available'] = True
        receipt['process_measurements'] = dict(user_cpu_seconds=1., system_cpu_seconds=.02)
        known = a.terminal_accounting(receipt, {}, '', 32.)
        self.assertEqual(known['raw_terminal_CPU_seconds'], 1.02)
        self.assertEqual(known['minimum_charged_CPU_seconds'], 1.04)
        footer = a.terminal_accounting(receipt, {}, 'CPU time consumed: 1min 2.125s\n', 64.)
        self.assertEqual(footer['raw_terminal_CPU_seconds'], 62.125)
        self.assertEqual(footer['minimum_charged_CPU_seconds'], 62.126)

    def test_preflight_cgroup_peak_uses_samples_and_rounded_footer_not_rss(self):
        receipt = dict(samples=[dict(service='MemoryPeak=406847488\n')],
                       process_measurements=dict(max_rss_kib=557232))
        peak = a.cgroup_peak(receipt, {}, 'Memory peak: 388.0M\n')
        self.assertGreater(peak, 406847488)
        self.assertLess(peak, 448*1024**2)
        self.assertGreater(a.cgroup_peak(receipt, {}, 'Memory peak: 450M\n'), 448*1024**2)

    def test_independent_execution_gate_rejects_masked_failure_and_keeps_checkpoint(self):
        parent = list(range(12))
        result = dict(order=[0, 2, 1]+list(range(3, 12)), score=dict(value=1e-9, cubic_response=1e-5, steps=10),
            status='complete', stats=dict(parent_method='magnitude', parent_evaluations=1,
                window_status='verified_checkpoint_after_child_cpu_limit', child_outcome='child_cpu_limit',
                child_returncode=-9, backend_runtime_assets_verified=True, native_parent_revalidation_evaluations=1,
                native_verifications=1, native_verification_attempts=1, child_completed_evaluations=10,
                child_allocated_evaluations=100, child_evaluation_scope='checkpoint_prefix', child_checked_updates=1,
                child_checkpoint_status='checked_update', child_counters=dict(evaluation=10, full_score=2,
                    flat_candidate=7, exact_witness=1, action=100, atom=10, insertion_move=1)))
        self.assertTrue(a.derived_execution('window', result, parent)['admitted'])
        for child in ('child_worker_failure', 'child_numerical_limit', 'child_admission_or_support_limit'):
            with self.subTest(child=child):
                bad = copy.deepcopy(result)
                bad['order'] = parent
                bad['stats'].update(window_status='native_no_resolved_improvement', child_outcome=child)
                self.assertFalse(a.derived_execution('window', bad, parent)['admitted'])
        complete = dict(order=parent, score=result['score'], status='complete', stats=dict(seed_evaluations=4))
        self.assertTrue(a.derived_execution('full', complete)['admitted'])
        complete['stats']['seed_evaluations'] = 3
        self.assertFalse(a.derived_execution('full', complete)['admitted'])

    def test_plan_grammar_binding_and_downward_rounding(self):
        config = json.loads((HERE/'PROTOCOL.json').read_text())
        parent, owner = record('parent'), record('full', 32, continuation=10.)
        owner['worker_output_path'] = 'workers/full32.json'
        command = ['python', '/source/numerical_audit/cost_matched_continuation/construct.py', '--worker',
            '--root', '/old/data', '--case', 'A', '--repetition', '1', '--arm', 'full', '--budget', '32',
            '--raw', '/old/data/workers/full32.json', '--cpu-cap', '31.0', '--wall-cap', '156']
        plan = dict(command=command, phase='constructors', cwd='/source', attempt_CPU_seconds=31., attempt_wall_seconds=156)
        receipt = dict(command=list(command), budget=dict(cpu_seconds=31., wall_seconds=156),
                       final_accounting=dict(CPU_seconds=9.98))
        worker = dict(guarded_CPU_cap=31., guarded_wall_cap=156, process_CPU_before_output=9.99)
        final = dict(charged_CPU_seconds=10.)
        self.assertEqual(a.bind_constructor(plan, receipt, final, owner, parent, worker, config), '/old/data')
        for suffix in (['--arm', 'full'], ['--unknown', 'value']):
            with self.subTest(suffix=suffix):
                bad = copy.deepcopy(plan)
                bad['command'] += suffix
                with self.assertRaisesRegex(ValueError, 'Duplicate/unknown'):
                    a.bind_constructor(bad, receipt, final, owner, parent, worker, config)
        for cpu, expected in ((31.-1e-13, True), (31.-1e-9, False), (31.+1e-13, False)):
            with self.subTest(cpu=cpu):
                altered = copy.deepcopy(plan)
                altered['command'][altered['command'].index('--cpu-cap')+1] = str(cpu)
                altered['attempt_CPU_seconds'] = cpu
                rec = dict(command=altered['command'], budget=dict(cpu_seconds=cpu, wall_seconds=156))
                witness = dict(worker, guarded_CPU_cap=cpu)
                if expected:
                    a.bind_constructor(altered, rec, final, owner, parent, witness, config)
                else:
                    with self.assertRaisesRegex(ValueError, 'remaining allowance'):
                        a.bind_constructor(altered, rec, final, owner, parent, witness, config)

    def test_scientific_job_guard_failure_filters_only_the_affected_case(self):
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            audit.original_roots = {'/old/data'}
            def add(label, phase, command, exit_code=0):
                audit.reconciled[label] = dict(plan=dict(phase=phase, command=command, cwd='/source'),
                    final=dict(exit_code=exit_code), receipt=dict(command=command))
            add('inputs', 'inputs', ['python', '/source/numerical_audit/cost_matched_continuation/preflight.py', '--output', '/old/data'])
            add('panels', 'packaging', ['python', '/source/numerical_audit/cost_matched_continuation/panels.py', '--root', '/old/data'])
            for case in ('A', 'B'):
                model, panel = '/old/data/inputs/'+case+'/model.toml', '/old/data/panels/'+case+'/candidates.toml'
                for phase, script, outputs in (('physics', 'physical.jl', ['/old/data/physical/'+case]),
                    ('objective', 'objective.jl', ['/old/data/objective/'+case]),
                    ('independent', 'verify.jl', ['/old/data/physical/'+case, '/old/data/independent/'+case])):
                    add(case+phase, phase, ['julia', '--startup-file=no', '--threads=1',
                        '/source/numerical_audit/cost_matched_continuation/'+script, model, panel]+outputs,
                        2 if case == 'A' and phase == 'objective' else 0)
            result = audit.scientific_jobs('preflight', ['A', 'B'])
            self.assertTrue(result['global_pass'])
            self.assertEqual(result['verified_cases'], ['B'])
            self.assertTrue(any('Missing/ambiguous successful scientific guard' in issue['message'] for issue in audit.issues))

    def test_restart_witness_and_full_start_graph_must_match(self):
        parent = self.graph_fixture()
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            self.assertTrue(audit.constructor_graph(parent, parent)['admitted'])
            changed = copy.deepcopy(parent)
            changed['restarts']['orders'][1] = [0, 3, 2, 1, 4]
            with self.assertRaisesRegex(ValueError, 'restart witness differs'):
                audit.constructor_graph(changed, changed)
            child = copy.deepcopy(parent)
            child.update(arm='full', total_CPU_budget=32)
            child['worker_verification']['worker_data'].update(starts=copy.deepcopy(parent['restarts']['orders']))
            child['worker_verification']['worker_data']['result']['stats']['seed_evaluations'] = 4
            self.assertTrue(audit.constructor_graph(child, parent)['admitted'])
            child['worker_verification']['worker_data']['starts'].reverse()
            with self.assertRaisesRegex(ValueError, 'starts differ'):
                audit.constructor_graph(child, parent)

    def test_fallback_binds_parent_and_parent_refusal_is_retained(self):
        parent = self.graph_fixture()
        child = copy.deepcopy(parent)
        child.update(arm='window', total_CPU_budget=32, eligible=False, controller_fallback=True)
        child['worker_verification'] = dict(verified=False, worker_data=dict(result=None, refusal=None, starts=[parent['order']]))
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            self.assertFalse(audit.constructor_graph(child, parent)['admitted'])
            child['order'] = [0, 2, 1, 3, 4]
            with self.assertRaisesRegex(ValueError, 'Fallback differs'):
                audit.constructor_graph(child, parent)
            refused_parent = copy.deepcopy(parent)
            refused_parent.update(eligible=False, controller_fallback=True, score=None)
            refused_parent['worker_verification'] = dict(verified=False, worker_data=None)
            self.assertFalse(audit.constructor_graph(refused_parent, refused_parent)['admitted'])

    def test_absent_worker_uses_explicit_refusal_without_losing_fallback_costs(self):
        expected = dict(version=1, admitted=False, category='worker_refusal', reason='outer_worker_refusal',
                        window_status=None, child_outcome=None, child_refusal=None)
        config = json.loads((HERE/'PROTOCOL.json').read_text())
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            parent = self.graph_fixture()
            parent['costs'] = record('parent')['costs']
            for arm in ('parent', 'full', 'window'):
                with self.subTest(arm=arm):
                    row = record(arm, 0 if arm == 'parent' else 32, continuation=0. if arm == 'parent' else 34.)
                    row.update(eligible=False, controller_fallback=True, guard_exit_code=124, result=None,
                        order=list(parent['order']), score=None if arm == 'parent' else copy.deepcopy(parent['score']),
                        restarts=None, worker_output_path='workers/missing-'+arm+'.json', worker_output_sha256=None,
                        execution_admission=copy.deepcopy(expected), analysis_eligible=False)
                    row['costs']['prefix' if arm == 'parent' else 'continuation']['raw_CPU_seconds'] = None
                    row['costs']['total']['raw_CPU_seconds'] = None
                    preserved_costs = copy.deepcopy(row['costs'])
                    row['worker_verification'] = audit.worker(row)
                    self.assertEqual(audit.constructor_graph(row, row if arm == 'parent' else parent), expected)
                    a.cost_check(row, row if arm == 'parent' else parent, config)
                    self.assertEqual(row['costs'], preserved_costs)
                    self.assertFalse(row['derived_execution_admission']['admitted'])
                    if arm == 'window':
                        physical = {('A', 'full', 1, 32, 1.35, 8): dict(projective_error=.01),
                                    ('A', 'window', 1, 32, 1.35, 8): dict(projective_error=.02)}
                        constructors = {('A', 1, 'full', 32): dict(analysis_eligible=True), ('A', 1, 'window', 32): row}
                        paired = a.pair_grid(['A'], 1, config, physical, constructors, set(physical))
                        retained = next(value for value in paired if value['total_CPU_budget'] == 32 and
                                        value['time'] == 1.35 and value['rotation_cap_multiplier'] == 8)
                        self.assertEqual(retained['outcome'], 'loss')
                        self.assertAlmostEqual(retained['difference'], .01)
                        self.assertFalse(retained['eligible'])

    def test_scored_parent_missing_measured_restarts_is_a_classified_failure(self):
        expected = dict(version=1, admitted=False, category='missing_parent_preparation',
            reason='missing_measured_restart_plan', window_status=None, child_outcome=None, child_refusal=None)
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            for refusal in (None, dict(type='ResourceLimit', message='Preparation interrupted')):
                with self.subTest(refusal=refusal):
                    parent = self.graph_fixture()
                    parent.update(eligible=False, restarts=None, execution_admission=copy.deepcopy(expected))
                    parent['worker_verification']['worker_data'].update(restarts=None, refusal=refusal)
                    selected = copy.deepcopy((parent['order'], parent['score']))
                    self.assertEqual(audit.constructor_graph(parent, parent), expected)
                    self.assertEqual((parent['order'], parent['score']), selected)
                    parent['eligible'] = True
                    with self.assertRaisesRegex(ValueError, 'Forged eligible'):
                        audit.constructor_graph(parent, parent)

    def test_parent_refusal_allows_unstarted_child_witness_and_never_admits_it(self):
        expected = dict(version=1, admitted=False, category='ineligible_parent', reason='frozen_parent_not_admitted',
                        window_status=None, child_outcome=None, child_refusal=None)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            audit = a.Audit(root)
            for parent_scored in (True, False):
                parent = self.graph_fixture()
                parent.update(eligible=False, restarts=None)
                parent['worker_verification']['worker_data']['restarts'] = None
                if not parent_scored:
                    parent.update(controller_fallback=True, score=None)
                    parent['worker_verification']['verified'] = False
                    parent['worker_verification']['worker_data'].update(result=None,
                        refusal=dict(type='ValueError', message='Parent construction refused'))
                self.assertFalse(audit.constructor_graph(parent, parent)['admitted'])
                for arm in ('full', 'window'):
                    with self.subTest(parent_scored=parent_scored, arm=arm):
                        identity = dict(case='A', repetition=1, arm=arm, total_CPU_budget=32,
                            model_json_sha256='a'*64, model_toml_sha256='b'*64, protocol_sha256='c'*64,
                            source_freeze_sha256='d'*64, parent_sha256='e'*64)
                        raw = dict(identity, schema='cost-matched-worker-v1', result=None, starts=None, restarts=None,
                            refusal=dict(type='ValueError', message='Frozen parent is ineligible; dependent construction refused'))
                        path = 'workers/refused-'+arm+'-'+str(parent_scored)+'.json'
                        digest = write_json(root, path, raw)
                        child = dict(identity, eligible=False, controller_fallback=True, guard_exit_code=0,
                            result=None, order=list(parent['order']), score=copy.deepcopy(parent['score']),
                            worker_output_path=path, worker_output_sha256=digest, execution_admission=copy.deepcopy(expected))
                        child['worker_verification'] = audit.worker(child)
                        self.assertEqual(audit.constructor_graph(child, parent), expected)
                        self.assertFalse(child['worker_verification']['verified'])
                        self.assertIsNone(child['worker_verification']['worker_data']['starts'])
            healthy = self.graph_fixture()
            child = copy.deepcopy(healthy)
            child.update(arm='full', total_CPU_budget=32, eligible=False)
            child['worker_verification']['worker_data'].update(starts=None)
            child['worker_verification']['worker_data']['result']['stats']['seed_evaluations'] = 4
            with self.assertRaisesRegex(ValueError, 'Absent continuation starts'):
                audit.constructor_graph(child, healthy)

    def test_forged_eligible_scored_parent_window_failure_is_rejected(self):
        parent = self.graph_fixture()
        child = copy.deepcopy(parent)
        child.update(arm='window', total_CPU_budget=64)
        raw = child['worker_verification']['worker_data']
        raw['starts'] = [parent['order']]
        raw['result']['stats'] = dict(window_status='julia_unavailable', parent_method='magnitude', parent_evaluations=1)
        with tempfile.TemporaryDirectory() as folder:
            audit = a.Audit(folder)
            with self.assertRaisesRegex(ValueError, 'Forged eligible'):
                audit.constructor_graph(child, parent)
            child['eligible'] = False
            derived = audit.constructor_graph(child, parent)
            self.assertFalse(derived['admitted'])
            self.assertEqual(child['order'], parent['order'])
            self.assertEqual(child['score'], parent['score'])
            child['execution_admission'] = dict(derived, admitted=True)
            with self.assertRaisesRegex(ValueError, 'execution admission differs'):
                audit.constructor_graph(child, parent)

    def test_main_preflight_proof_binds_consumer_delta_and_recomputes_gate(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            consumer = 'numerical_audit/cost_matched_continuation/analyze.py'
            source_bytes = b'# Archived consumer repair fixture.\n'
            digest = hashlib.sha256(source_bytes).hexdigest()
            pref = dict(mode='preflight', sources={consumer: '0'*64, 'kernel.py': '1'*64},
                configuration=dict(training_time=1.), registry=dict(cases=['A', 'B']), environment=dict(python='fixture'))
            pref_sha = write_json(root, 'preflight/INPUT_PROTOCOL_FROZEN.json', pref)
            main = copy.deepcopy(pref)
            main['mode'] = 'main'
            main['sources'][consumer] = digest
            source_path = 'preflight/validation/analyze.'+digest+'.py'
            (root/source_path).parent.mkdir(parents=True, exist_ok=True)
            (root/source_path).write_bytes(source_bytes)
            saved = dict(study_complete=True, mode='preflight', source_freeze_sha256=pref_sha,
                analyzer_sha256=digest, primary=dict(decision='preflight_only'), memory_gate=dict(passed=True, ceiling_MiB=448))
            analysis_path = 'preflight/validation/ANALYSIS.'+digest+'.json'
            analysis_sha = write_json(root, analysis_path, saved)
            main['preflight'] = dict(path='preflight', analysis_path=analysis_path, analysis_sha256=analysis_sha,
                source_freeze_sha256=pref_sha, consumer_source_path=source_path, consumer_sha256=digest,
                consumer_only_source_delta={consumer: dict(original='0'*64, replayed=digest)})
            replay = dict(saved, input_digests={}, scientific_cutoff_ledger=dict(
                scientific_output_guards=dict(original_data_root='/old/pref')))
            with mock.patch.object(a, 'analyze', return_value=replay) as recursive:
                audit = a.Audit(root)
                self.assertTrue(a.preflight_proof(audit, main)['verified'])
                recursive.assert_called_once_with(root/'preflight', preflight_only=True)
                self.assertEqual(audit.preflight_original_root, '/old/pref')
            changed = copy.deepcopy(main)
            changed['sources']['kernel.py'] = '2'*64
            with self.assertRaisesRegex(ValueError, 'source delta'):
                a.preflight_proof(a.Audit(root), changed)
            failed = copy.deepcopy(replay)
            failed['study_complete'] = False
            with mock.patch.object(a, 'analyze', return_value=failed):
                with self.assertRaisesRegex(ValueError, 'does not admit'):
                    a.preflight_proof(a.Audit(root), main)

    def test_incomplete_analyze_requires_input_freeze_without_erasing_metrics(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            config = json.loads((HERE/'PROTOCOL.json').read_text())
            config['preflight']['case'] = 'A'
            protocol = 'numerical_audit/cost_matched_continuation/PROTOCOL.json'
            protocol_sha = write_json(root, 'source_snapshot/'+protocol, config)
            amendment = 'numerical_audit/cost_matched_continuation/PREFLIGHT_RESOURCE_AMENDMENT.json'
            amendment_sha = write_json(root, 'source_snapshot/'+amendment, json.loads((HERE/'PREFLIGHT_RESOURCE_AMENDMENT.json').read_text()))
            source_sha = write_json(root, 'INPUT_PROTOCOL_FROZEN.json', dict(mode='preflight',
                configuration=config, configuration_sha256=protocol_sha,
                sources={protocol: protocol_sha, amendment: amendment_sha},
                registry=dict(cases=[dict(case='X'), dict(case='Y')]), environment={}))
            header = dict(mode='preflight', cases=['A'], repetitions=1, protocol_sha256=protocol_sha,
                source_freeze_sha256=source_sha, protocol_path='source_snapshot/'+protocol,
                source_freeze_path='INPUT_PROTOCOL_FROZEN.json', entries=[])
            parents_sha = write_json(root, 'PARENTS_FROZEN.json', header)
            constructors_sha = write_json(root, 'CONSTRUCTORS_FROZEN.json', header)
            write_json(root, 'PANELS_FROZEN.json', dict(header, input_freeze_sha256='0'*64,
                parents_freeze_sha256=parents_sha, constructors_freeze_sha256=constructors_sha))
            rows = []
            for arm in ('full', 'window'):
                for budget in (32, 64):
                    for duration in (.65, 1.35):
                        for cap in (8, 24):
                            values = dict(arm=arm, repetition=1, cpu_budget=budget, time=duration,
                                budget_multiplier=cap, name=arm+str(budget)+str(duration)+str(cap),
                                projective_error=.02 if arm == 'window' else .01)
                            rows.append('[[rows]]\n'+'\n'.join(key+'='+json.dumps(value) for key, value in values.items()))
            target = root/'physical/A/PHYSICAL.toml'
            target.parent.mkdir(parents=True)
            target.write_text('\n'.join(rows)+'\n')
            result = a.analyze(root, preflight_only=True)
            self.assertFalse(result['study_complete'])
            self.assertEqual(result['primary']['decision'], 'study_incomplete')
            self.assertEqual(len(result['physical_rows']), 16)
            self.assertEqual(len(result['paired_rows']), 8)
            self.assertTrue(all(row['outcome'] == 'loss' and not row['eligible'] for row in result['paired_rows']))
            self.assertTrue(any(issue['scope'] == 'input freeze' for issue in result['issues']))
            self.assertFalse(result['resources']['inventory_complete'])


if __name__ == '__main__':
    unittest.main()
