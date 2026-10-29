"""Synthetic producer contracts; no molecules, native optimization or services.

All costs and clocks below are injected fixture values, not scientific evidence.
"""
import copy
from dataclasses import dataclass
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
        self.parent = dict(order=list(range(6)), score=dict(value=1e-6, cubic_response=1e-4, steps=3), eligible=True,
            costs=dict(prefix=dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.),
                       total=dict(raw_CPU_seconds=1., charged_CPU_seconds=2., wall_seconds=10.)))

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
