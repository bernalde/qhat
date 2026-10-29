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
                receipt_path='campaign/jobs/'+label+'.receipt.json')
            reference['receipt_sha256'] = write_json(root, reference['receipt_path'], dict(settings=settings,
                exit_code=2 if label == 'failed-attempt' else 0, elapsed=cpu*4, final_accounting=dict(CPU_seconds=cpu)))
            reference['final_sha256'] = write_json(root, reference['final_path'], dict(label=label, phase='constructors',
                receipt_sha256=reference['receipt_sha256'], charged_CPU_seconds=cpu,
                absolute_guarded_wall_seconds=cpu*4, exit_code=2 if label == 'failed-attempt' else 0))
            plan_path = 'campaign/jobs/'+label+'.PLAN.json'
            plan_sha = write_json(root, plan_path, dict(label=label, phase='constructors', pools_sha256=pools_sha,
                attempt_CPU_seconds=cpu+1, attempt_wall_seconds=cpu*4+5))
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


if __name__ == '__main__':
    unittest.main()
