"""Synthetic producer contracts; no molecules, native optimization or services.

All costs and clocks below are injected fixture values, not scientific evidence.
"""
import copy
from dataclasses import asdict, dataclass
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


k = load('continuation_contract_common', HERE/'common.py')
admission = load('continuation_contract_admission', HERE/'admission.py')
provenance = load('continuation_contract_provenance', HERE/'provenance.py')
with patch.dict(sys.modules, {'common': k, 'admission': admission, 'provenance': provenance}):
    producer = load('continuation_contract_producer', HERE/'construct.py')
    panels = load('continuation_contract_panels', HERE/'panels.py')


@dataclass
class FakeScore:
    value: float = 1e-6
    cubic_response: float = 1e-4
    steps: int = 3


@dataclass
class FakeResult:
    order: list
    score: FakeScore
    CPU_seconds: float
    evaluations: int
    accepted_moves: int
    status: str
    stats: dict


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, allow_nan=False))


class ProducerContracts(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)/'data'; self.root.mkdir()
        self.campaign = Path(temporary.name)/'run'
        self.config = json.loads((HERE/'PROTOCOL.json').read_text())
        self.info = dict(protocol_sha256='a'*64, source_freeze_sha256='b'*64)
        self.model = dict(case='A', n=14, electrons=4, blocks=[{} for _ in range(6)], coefficients=[0]*6)
        write_json(self.root/'inputs/A/model.json', self.model)
        (self.root/'inputs/A/model.toml').write_text('case = "synthetic fixture A"\n')
        self.parent = dict(case='A', repetition=1, arm='parent', order=list(range(6)),
            score=dict(value=1e-6, cubic_response=1e-4, steps=3), eligible=True,
            costs=dict(prefix=dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.),
                       total=dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.)))
        self.parent['restarts'] = producer.restarts(self.parent['order'], 'A', 1)

    def owning_job(self, arm, raw):
        budget = 0 if arm == 'parent' else 32
        cpu, wall = (8., 40) if arm == 'parent' else (30., 150)
        command = [sys.executable, str(HERE/'construct.py'), '--worker', '--root', str(self.root),
            '--case', 'A', '--repetition', '1', '--arm', arm, '--budget', str(budget),
            '--raw', str(raw), '--cpu-cap', str(cpu), '--wall-cap', str(wall)]
        plan = dict(phase='constructors', label='fixture-a1', cwd=str(k.ROOT), command=command,
                    attempt_CPU_seconds=cpu, attempt_wall_seconds=wall)
        receipt = dict(command=list(command), budget=dict(cpu_seconds=cpu, wall_seconds=wall))
        owner = dict(case='A', repetition=1, arm=arm, total_CPU_budget=budget,
                     worker_output_path=str(raw.relative_to(self.root)))
        return plan, receipt, owner

    def finalize(self, arm, raw, cost, exit_code):
        plan, receipt, _ = self.owning_job(arm, raw)
        write_json(self.campaign/'jobs/fixture-a1.PLAN.json', plan)
        write_json(self.campaign/'logs/fixture-a1.resources.json', receipt)
        with patch.object(producer, 'identity', return_value=(self.info, self.config, self.model)), \
                patch.object(producer, 'parent_record', return_value=(self.parent, 'c'*64)), \
                patch.object(producer, 'measured_cost', return_value=(cost, exit_code)):
            return producer.finalize(self.root, 'A', 1, arm, 0 if arm == 'parent' else 32,
                                     raw, self.campaign, 'fixture-a1')

    def parent_worker(self, raw, trained, restart_failure=None):
        with patch.object(producer, 'identity', return_value=(self.info, self.config, self.model)), \
                patch.object(k.prior, 'blocks', return_value=tuple({(): 1} for _ in range(6))), \
                patch.object(k, 'bootstrap', return_value=SimpleNamespace(train_ordering=lambda *a, **kw: trained)), \
                patch.object(producer, 'time', SimpleNamespace(process_time=lambda: .1, monotonic=lambda: 1.)), \
                patch.object(producer, 'BOOT_WALL', 0.), \
                patch.object(producer, 'restarts', wraps=producer.restarts, side_effect=restart_failure):
            producer.worker(self.root, 'A', 1, 'parent', 0, 8., 40, raw, None)
        return k.read(raw)

    def child_worker(self, arm):
        return dict(schema='cost-matched-worker-v1', case='A', repetition=1, arm=arm, total_CPU_budget=32,
            parent_sha256='c'*64, **self.info,
            model_json_sha256=k.sha(self.root/'inputs/A/model.json'),
            model_toml_sha256=k.sha(self.root/'inputs/A/model.toml'),
            result=asdict(FakeResult(list(self.parent['order']), FakeScore(), .1, 4, 0, 'complete', dict(seed_evaluations=4))),
            refusal=None, starts=copy.deepcopy(self.parent['restarts']['orders'] if arm == 'full' else [self.parent['order']]),
            guarded_CPU_cap=30., guarded_wall_cap=150, process_CPU_before_output=1.)

    def test_parent_restart_work_precedes_raw_write_and_finalizer_copies_without_rng(self):
        events, clock = [], SimpleNamespace(value=0.)
        def process_time():
            clock.value += .1
            events.append('clock')
            return clock.value
        def train(*args, **kwargs):
            events.append('train')
            self.assertIsNone(kwargs['starts'])
            return FakeResult(list(range(6)), FakeScore(), .1, 1, 0, 'complete', dict(seed_evaluations=1))
        original_restarts, original_write = producer.restarts, k.write_new
        def restarts(*args):
            events.append('restarts')
            return original_restarts(*args)
        def store(path, value):
            events.append('raw_write')
            self.assertIsNotNone(value['restarts'])
            return original_write(path, value)
        raw = self.root/'workers/parent.json'
        with patch.object(producer, 'identity', return_value=(self.info, self.config, self.model)), \
                patch.object(k.prior, 'blocks', return_value=tuple({(): 1} for _ in range(6))), \
                patch.object(k, 'bootstrap', return_value=SimpleNamespace(train_ordering=train)), \
                patch.object(producer, 'time', SimpleNamespace(process_time=process_time, monotonic=lambda: 1.)), \
                patch.object(producer, 'BOOT_WALL', 0.), patch.object(producer, 'restarts', side_effect=restarts), \
                patch.object(k, 'write_new', side_effect=store):
            producer.worker(self.root, 'A', 1, 'parent', 0, 8., 40, raw, None)
        worker = k.read(raw)
        self.assertIsNone(worker['refusal'])
        self.assertLess(events.index('train'), events.index('restarts'))
        self.assertLess(events.index('restarts'), events.index('raw_write'))
        self.assertIn('clock', events[events.index('restarts')+1:events.index('raw_write')])
        self.assertIsInstance(worker['restart_generation_CPU_seconds'], float)
        self.assertEqual(worker['restarts']['orders'][0], worker['result']['order'])
        self.assertEqual(len({tuple(order) for order in worker['restarts']['orders']}), 4)
        cost = dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.)
        with patch.object(producer, 'restarts', side_effect=AssertionError('finalizer cannot regenerate starts')), \
                patch.object(producer.random, 'Random', side_effect=AssertionError('finalizer cannot call RNG')):
            parent = self.finalize('parent', raw, cost, 0)
        self.assertTrue(parent['eligible'])
        self.assertEqual(parent['restarts'], worker['restarts'])
        self.assertEqual(k.read(producer.canonical_path(self.root, 'A', 1, 'parent', 0))['restarts'], worker['restarts'])

    def test_invalid_candidates_leave_complete_refusals_without_promoting_result(self):
        invalid = [FakeResult(list(range(6)), FakeScore(value=value), .1, 1, 0, 'complete', {})
                   for value in (float('nan'), float('inf'), -1.)]
        invalid += [FakeResult(list(range(6)), FakeScore(cubic_response=float('nan')), .1, 1, 0, 'complete', {}),
                    FakeResult(list(range(6)), FakeScore(steps=float('nan')), .1, 1, 0, 'complete', {}),
                    FakeResult([0, 0, 2, 3, 4, 5], FakeScore(), .1, 1, 0, 'complete', {})]
        for index, trained in enumerate(invalid):
            with self.subTest(index=index):
                raw = self.root/'workers'/('invalid-'+str(index)+'.json')
                worker = self.parent_worker(raw, trained)
                self.assertIsNone(worker['result'])
                self.assertIsNotNone(worker['refusal'])
                self.assertIsNone(worker['restarts'])
                self.assertFalse(raw.with_name(raw.name+'.partial').exists())

    def test_restart_preparation_failure_retains_the_valid_scored_parent(self):
        raw = self.root/'workers/preparation-refusal.json'
        trained = FakeResult([0, 2, 1, 3, 4, 5], FakeScore(), .1, 1, 0, 'complete', {})
        worker = self.parent_worker(raw, trained, ValueError('synthetic preparation refusal'))
        self.assertEqual(worker['result'], asdict(trained))
        self.assertIsNone(worker['restarts'])
        self.assertIn('synthetic preparation refusal', worker['refusal']['message'])
        with patch.object(producer, 'restarts', side_effect=AssertionError('no replacement starts')):
            record = self.finalize('parent', raw, dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.), 0)
        self.assertEqual(record['order'], trained.order)
        self.assertEqual(record['score'], asdict(trained.score))
        self.assertFalse(record['controller_fallback'])
        self.assertFalse(record['eligible'])
        self.assertEqual(record['execution_admission']['category'], 'missing_parent_preparation')
        self.assertEqual(panels.worker_witness(panels.Artifacts(self.root), record, record), worker)

    def test_finalizer_binds_nonfallback_full_and_window_starts(self):
        for arm in ('full', 'window'):
            with self.subTest(arm=arm):
                raw = self.root/'workers'/(arm+'.json')
                worker = self.child_worker(arm)
                worker['starts'] = []
                write_json(raw, worker)
                original = raw.read_bytes()
                with self.assertRaisesRegex(ValueError, 'Continuation input starts'):
                    self.finalize(arm, raw, dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.), 0)
                self.assertEqual(raw.read_bytes(), original)
                self.assertFalse(producer.canonical_path(self.root, 'A', 1, arm, 32).exists())
        raw = self.root/'workers/valid-full.json'
        write_json(raw, self.child_worker('full'))
        self.assertTrue(self.finalize('full', raw, dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.), 0)['eligible'])

    def test_panel_checks_restart_seed_counts_and_interior_permutations_without_rng(self):
        with patch.object(producer.random, 'Random', side_effect=AssertionError('verification cannot select starts')):
            panels.parent_preparation(self.parent)
        changes = [lambda p: p.update(seed=p['seed']+1),
                   lambda p: p.update(interior_perturbation_counts=[0, 1, 2, 3]),
                   lambda p: p['orders'].__setitem__(1, list(p['orders'][0])),
                   lambda p: p['orders'].__setitem__(1, [1, 0, 2, 3, 4, 5]),
                   lambda p: p['orders'].__setitem__(1, [0, 1, 1, 3, 4, 5])]
        for index, change in enumerate(changes):
            with self.subTest(index=index):
                parent = copy.deepcopy(self.parent)
                change(parent['restarts'])
                with self.assertRaises(ValueError):
                    panels.parent_preparation(parent)

    def test_worker_owner_receipt_destination_and_remaining_caps_are_bound(self):
        raw = self.root/'workers/full.json'
        plan, receipt, owner = self.owning_job('full', raw)
        worker = dict(guarded_CPU_cap=30., guarded_wall_cap=150)
        self.assertTrue(provenance.bind_worker(plan, receipt, owner, self.parent, worker, self.config)['verified'])
        def argument(record, key, value):
            index = record['command'].index(key)+1
            record['command'][index] = value
        changes = (
            ('swapped_receipt_arm', lambda p, r, w: argument(r, '--arm', 'window')),
            ('swapped_owning_arm', lambda p, r, w: (argument(p, '--arm', 'window'), argument(r, '--arm', 'window'))),
            ('wrong_raw_destination', lambda p, r, w: (argument(p, '--raw', str(self.root/'workers/window.json')),
                                                     argument(r, '--raw', str(self.root/'workers/window.json')))),
            ('reduced_remaining_cpu', lambda p, r, w: (argument(p, '--cpu-cap', '29'), argument(r, '--cpu-cap', '29'))),
            ('receipt_cpu_cap', lambda p, r, w: r['budget'].update(cpu_seconds=31.)),
            ('worker_cpu_cap', lambda p, r, w: w.update(guarded_CPU_cap=29.)),
            ('wrong_phase', lambda p, r, w: p.update(phase='physics')))
        for name, change in changes:
            with self.subTest(name=name):
                changed_plan, changed_receipt, changed_worker = copy.deepcopy((plan, receipt, worker))
                change(changed_plan, changed_receipt, changed_worker)
                with self.assertRaises(ValueError):
                    provenance.bind_worker(changed_plan, changed_receipt, owner, self.parent, changed_worker, self.config)

    def test_parent_identity_is_bound_without_requiring_parent_eligibility(self):
        raw = self.root/'workers/full.json'
        plan, receipt, owner = self.owning_job('full', raw)
        for field, value in (('case', 'B'), ('repetition', 2), ('arm', 'full')):
            with self.subTest(field=field):
                parent = dict(self.parent, **{field: value})
                with self.assertRaisesRegex(ValueError, 'Frozen parent identity'):
                    provenance.bind_worker(plan, receipt, owner, parent, None, self.config)
        refused = dict(self.parent, eligible=False)
        self.assertTrue(provenance.bind_worker(plan, receipt, owner, refused, None, self.config)['verified'])
        with self.assertRaisesRegex(ValueError, 'Frozen parent identity'):
            provenance.bind_worker(plan, receipt, owner, None, None, self.config)

    def test_legacy_cpu_needs_terminal_process_or_successful_postexit_evidence(self):
        for partial in (0., 3.5):
            with self.subTest(partial=partial):
                receipt = dict(final_accounting=dict(CPU_seconds=partial, process_CPU_timing_available=False),
                               samples=[dict(service='CPUUsageNSec=0')])
                self.assertIsNone(provenance.raw_cpu({}, receipt))
        receipt = dict(final_accounting=dict(CPU_seconds=1.25, process_CPU_timing_available=True),
                       process_measurements=dict(user_cpu_seconds=1., system_cpu_seconds=.2))
        self.assertEqual(provenance.raw_cpu({}, receipt), 1.25)
        receipt['process_measurements']['system_cpu_seconds'] = None
        self.assertIsNone(provenance.raw_cpu({}, receipt))
        receipt = dict(final_accounting=dict(CPU_seconds=0., process_CPU_timing_available=True),
                       process_measurements=dict(user_cpu_seconds=0., system_cpu_seconds=0.))
        self.assertEqual(provenance.raw_cpu({}, receipt), 0.)
        receipt['final_accounting']['process_CPU_timing_available'] = False
        final = dict(final_service_query_exit=0, final_service_snapshot=dict(CPUUsageNSec='1234000000'))
        self.assertEqual(provenance.raw_cpu(final, receipt), 1.234)
        self.assertIsNone(provenance.raw_cpu(dict(final, final_service_query_exit=1), receipt))
        self.assertIsNone(provenance.raw_cpu(dict(final, final_service_snapshot=dict(CPUUsageNSec='٠')), receipt))
        current = dict(terminal_CPU_measurement_available=False, raw_terminal_CPU_seconds=None)
        self.assertIsNone(provenance.raw_cpu(current, receipt))
        with self.assertRaises(ValueError):
            provenance.raw_cpu(dict(current, raw_terminal_CPU_seconds=0.), receipt)

    def test_artifact_null_pins_reject_but_unindexed_reads_and_legacy_plan_omission_remain(self):
        raw = self.root/'workers/full.json'
        plan, receipt, owner = self.owning_job('full', raw)
        receipt.update(exit_code=0, elapsed=8., final_accounting=dict(CPU_seconds=1., process_CPU_timing_available=True),
                       process_measurements=dict(user_cpu_seconds=.8, system_cpu_seconds=.2))
        prefix = 'campaign/'
        plan_path = prefix+'jobs/fixture-a1.PLAN.json'
        receipt_path = prefix+'logs/fixture-a1.resources.json'
        final_path = prefix+'jobs/fixture-a1.FINAL.json'
        write_json(self.root/plan_path, plan)
        write_json(self.root/receipt_path, receipt)
        final = dict(phase='constructors', label='fixture-a1', exit_code=0, receipt_sha256=k.sha(self.root/receipt_path),
                     charged_CPU_seconds=2., absolute_guarded_wall_seconds=8.)
        write_json(self.root/final_path, final)
        cost = dict(label='fixture-a1', final_path=final_path, final_sha256=k.sha(self.root/final_path),
            receipt_path=receipt_path, receipt_sha256=k.sha(self.root/receipt_path),
            raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=8.)
        artifacts = panels.Artifacts(self.root)
        self.assertEqual(artifacts.load(plan_path), plan)
        self.assertEqual(artifacts.load(plan_path, k.sha(self.root/plan_path)), plan)
        with self.assertRaisesRegex(ValueError, 'SHA-256 spelling'):
            artifacts.load(plan_path, None)
        with self.assertRaisesRegex(ValueError, 'SHA-256 spelling'):
            panels.authenticate_model(artifacts, 'A', {'model.json': None,
                'model.toml': k.sha(self.root/'inputs/A/model.toml')})
        owner['guard_exit_code'] = 0
        panels.cost_evidence(artifacts, cost, owner, self.parent, None, self.config)
        for fields in (dict(final_sha256=None), dict(receipt_sha256=None), dict(plan_path=plan_path),
                       dict(plan_path=plan_path, plan_sha256=None)):
            with self.subTest(fields=fields):
                with self.assertRaises(ValueError):
                    panels.cost_evidence(artifacts, dict(cost, **fields), owner, self.parent, None, self.config)

    def test_unknown_terminal_cpu_retains_full_charge_and_ineligible_parent_fallback(self):
        raw = self.root/'workers/missing-full.json'
        cost = dict(raw_CPU_seconds=None, charged_CPU_seconds=30., wall_seconds=20.)
        record = self.finalize('full', raw, cost, 1)
        self.assertIsNone(record['costs']['continuation']['raw_CPU_seconds'])
        self.assertIsNone(record['costs']['total']['raw_CPU_seconds'])
        self.assertEqual(record['costs']['continuation']['charged_CPU_seconds'], 30.)
        self.assertEqual(record['costs']['total']['charged_CPU_seconds'], 32.)
        self.assertFalse(record['eligible'])
        self.assertTrue(record['controller_fallback'])
        self.assertEqual(record['order'], self.parent['order'])
        self.assertEqual(record['score'], self.parent['score'])
        self.assertEqual(record['parent_sha256'], 'c'*64)
        self.assertIsNone(record['worker_output_sha256'])
        self.assertTrue(record['cost_binding']['verified'])
        self.assertIsNone(panels.worker_witness(panels.Artifacts(self.root), record, self.parent))

    def test_source_inventory_covers_consumer_and_new_producers_but_not_runtime_subdirectories(self):
        folder = self.root/'source'; folder.mkdir()
        for name in ('analyze.py', 'construct.py', 'new_producer.py', 'PROTOCOL.json'):
            (folder/name).write_text('{}\n' if name.endswith('.json') else '# synthetic source\n')
        (folder/'runtime').mkdir(); (folder/'runtime/constructor.json').write_text('{}\n')
        with patch.object(k, 'HERE', folder), patch.object(k, 'ROOT', self.root), \
                patch.object(k.prior, 'source_paths', return_value=[]):
            inventory = k.source_map()
        self.assertIn('source/analyze.py', inventory)
        self.assertIn('source/new_producer.py', inventory)
        self.assertNotIn('source/runtime/constructor.json', inventory)


if __name__ == '__main__':
    unittest.main()
