"""Standard-library metadata fixtures; no guard or dispatcher execution."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1]/'ledger.py'
spec = importlib.util.spec_from_file_location('continuation_ledger_metadata', SOURCE)
ledger = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ledger)


class LedgerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.run = Path(self.temporary.name)/'run'
        self.data = self.run/'data'
        for path in (self.data, self.run/'jobs', self.run/'logs', self.run/'campaign'):
            path.mkdir(parents=True, exist_ok=True)
        self.put('campaign/resource_guard.py', b'guard archive, never executed\n')
        self.put('campaign/run_guarded_job.py', b'original dispatcher, never executed\n')
        self.pools = dict(schema='cost-matched-resource-pools-v1',
            guard=dict(path='campaign/resource_guard.py', sha256=self.hash('campaign/resource_guard.py')),
            phases=dict(constructors=dict(CPU_seconds=100, wall_seconds=500)))
        self.put_json('RESOURCE_POOLS.json', self.pools)

    def put(self, name, raw):
        path = self.run/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)

    def put_json(self, name, value):
        self.put(name, (json.dumps(value, indent=2)+'\n').encode())

    def hash(self, name):
        return hashlib.sha256((self.run/name).read_bytes()).hexdigest()

    def job(self, label='job-a1', versioned=False, exit_code=0):
        source = 'campaign/run_guarded_job.py'
        if versioned:
            raw = b'new dispatcher version, never executed\n'
            source = 'campaign/dispatcher_versions/'+hashlib.sha256(raw).hexdigest()+'.py'
            self.put(source, raw)
        plan = dict(label=label, phase='constructors', pools_sha256=self.hash('RESOURCE_POOLS.json'),
                    dispatcher_sha256=self.hash(source), command=['python', 'worker.py'])
        if versioned:
            plan['dispatcher_archive_path'] = source
        self.put_json('jobs/'+label+'.PLAN.json', plan)
        self.put_json('logs/'+label+'.resources.json', dict(exit_code=exit_code, elapsed=1.))
        self.put('logs/'+label+'.log', b'failed or successful job output retained\n')
        self.put_json('jobs/'+label+'.FINAL.json', dict(label=label, phase='constructors', exit_code=exit_code,
            receipt_sha256=self.hash('logs/'+label+'.resources.json'), charged_CPU_seconds=2.))
        return plan

    def test_canonical_layout_preserves_every_job_and_source_byte(self):
        self.job(); self.job('failed-a1', versioned=True, exit_code=7)
        result = ledger.archive(self.data, self.run)
        self.assertEqual(len(result['jobs']), 2)
        self.assertEqual(result['data_root_relative_to_run_root'], 'data')
        self.assertIn('cumulative RUN_ROOT', result['disk_accounting_scope'])
        for entry in result['jobs']:
            for kind in ('plan', 'final', 'receipt', 'log'):
                archived = self.data/entry[kind+'_path']
                source = self.run/entry[kind+'_path'].removeprefix('campaign/')
                self.assertEqual(archived.read_bytes(), source.read_bytes())
                self.assertEqual(hashlib.sha256(archived.read_bytes()).hexdigest(), entry[kind+'_sha256'])
            self.assertEqual((self.data/entry['dispatcher_path']).read_bytes(),
                             (self.run/entry['dispatcher_path']).read_bytes())
        self.assertEqual((self.data/'campaign/resource_guard.py').read_bytes(),
                         (self.run/'campaign/resource_guard.py').read_bytes())
        self.assertEqual(json.loads((self.data/'campaign/jobs/failed-a1.FINAL.json').read_text())['exit_code'], 7)

    def test_identical_existing_evidence_reused_and_cutoff_preserved(self):
        self.job()
        target = self.data/'campaign/jobs/job-a1.PLAN.json'; target.parent.mkdir(parents=True)
        target.write_bytes((self.run/'jobs/job-a1.PLAN.json').read_bytes())
        ledger.archive(self.data, self.run)
        before = (self.data/'campaign/LEDGER_CUTOFF.json').read_bytes()
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertEqual((self.data/'campaign/LEDGER_CUTOFF.json').read_bytes(), before)

    def test_each_orphan_artifact_refuses_before_archive_writes(self):
        self.job()
        for name in ('jobs/orphan.FINAL.json', 'logs/orphan.resources.json', 'logs/orphan.log'):
            with self.subTest(name=name):
                self.put(name, b'{}\n')
                with self.assertRaises(ValueError):
                    ledger.archive(self.data, self.run)
                self.assertFalse((self.data/'campaign').exists())
                (self.run/name).unlink()

    def test_plan_filename_and_ascii_label_authenticated(self):
        self.job()
        path = self.run/'jobs/job-a1.PLAN.json'
        plan = json.loads(path.read_text()); plan['label'] = 'other-a1'
        path.write_text(json.dumps(plan))
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertFalse((self.data/'campaign').exists())
        for name in ('étiquette', '../escape', 'a'*129, '-leading'):
            self.assertIsNone(ledger.LABEL.fullmatch(name))

    def test_missing_final_blocks_cutoff(self):
        self.job(); (self.run/'jobs/job-a1.FINAL.json').unlink()
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertFalse((self.data/'campaign').exists())

    def test_guard_schema_and_external_data_root_rejected(self):
        self.job(); self.pools['guard']['path'] = 'resource_guard.py'
        self.put_json('RESOURCE_POOLS.json', self.pools)
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        outside = Path(self.temporary.name)/'outside'; outside.mkdir()
        with self.assertRaises(ValueError):
            ledger.archive(outside, self.run)

    def test_external_source_symlink_and_stale_archive_preserved(self):
        self.job()
        path = self.run/'logs/job-a1.log'; original = path.read_bytes(); path.unlink()
        external = Path(self.temporary.name)/'outside.log'; external.write_bytes(original)
        path.symlink_to(external)
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertFalse((self.data/'campaign').exists())
        path.unlink(); path.write_bytes(original)
        stale = self.data/'campaign/logs/job-a1.log'; stale.parent.mkdir(parents=True)
        stale.write_bytes(b'stale original archive\n')
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertEqual(stale.read_bytes(), b'stale original archive\n')
        self.assertFalse((self.data/'campaign/LEDGER_CUTOFF.json').exists())

    def test_dispatcher_path_bound_to_digest(self):
        self.job(versioned=True)
        path = self.run/'jobs/job-a1.PLAN.json'
        plan = json.loads(path.read_text()); plan['dispatcher_archive_path'] = '../outside.py'
        path.write_text(json.dumps(plan))
        with self.assertRaises(ValueError):
            ledger.archive(self.data, self.run)
        self.assertFalse((self.data/'campaign').exists())


if __name__ == '__main__':
    unittest.main()
