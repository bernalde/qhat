"""Synthetic metadata contracts only; these fixtures are not scored evidence."""
import copy
import hashlib
import importlib
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
p = importlib.import_module('numerical_audit.molecular_generalization.panels')
a = importlib.import_module('numerical_audit.molecular_generalization.analyze')
CASE = 'LiH_r1p487_631g_14m4e'


def put(root, relative, value):
    path = root/relative; path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False)+'\n')
    return p.k.sha(path)


def starts(order):
    seed = int.from_bytes(hashlib.sha256(f'{CASE}:1:cost-matched-v1'.encode()).digest()[:8], 'big')
    rng = random.Random(seed); orders = [list(order)]
    for count in (1, 2, 4):
        for _ in range(100):
            trial = list(order)
            for _ in range(count):
                source, target = rng.sample(range(1, len(order)-1), 2)
                trial.insert(target, trial.pop(source))
            if trial not in orders:
                orders.append(trial); break
    return dict(seed=seed, orders=orders, interior_perturbation_counts=[0, 1, 2, 4])


def bundle(root):
    config = json.loads((p.k.HERE/'PROTOCOL.json').read_text())
    config.update(main_execution_allowed=False, protocol_state='old-case_preflight_development')
    protocol = 'source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json'
    guard = 'source_snapshot/numerical_audit/molecular_generalization/tools/resource_guard.py'
    protocol_sha = put(root, protocol, config)
    (root/guard).parent.mkdir(parents=True, exist_ok=True); (root/guard).write_text('# synthetic guard source\n')
    freeze = dict(schema='molecular-source-freeze-v1', mode='preflight', configuration=config,
        sources={protocol.removeprefix('source_snapshot/'): protocol_sha,
                 guard.removeprefix('source_snapshot/'): p.k.sha(root/guard)}, registry=dict(cases=[]),
        producer_root='/original/source', producer_data_root='/original/data')
    source_sha = put(root, 'INPUT_PROTOCOL_FROZEN.json', freeze)
    info = dict(mode='preflight', cases=[CASE], repetitions=1, protocol_sha256=protocol_sha,
        source_freeze_sha256=source_sha, protocol_path=protocol, source_freeze_path='INPUT_PROTOCOL_FROZEN.json')
    terms = [[[[q, 'Z']], str(q+1)] for q in range(5)]
    model = dict(case=CASE, n=14, electrons=4, denominator='100', coefficients=terms,
                 blocks=[[term] for term in terms])
    model_sha = put(root, f'inputs/{CASE}/model.json', model)
    target = dict(case=CASE, n=14, electrons=4, denominator='100', partition='conserving',
        source_json_sha256=model_sha, baseline_order=[1, 2, 3, 4, 5],
        units=[dict(words=[str(q)+'Z'], numerators=[str(q+1)]) for q in range(5)])
    p.write_manifest(root/f'inputs/{CASE}/model.toml', target, ('units',))
    model_toml_sha = p.k.sha(root/f'inputs/{CASE}/model.toml')
    put(root, 'INPUTS_FROZEN.json', dict(info, status='complete', refusals={},
        successes={CASE: {'model.json': model_sha, 'model.toml': model_toml_sha}}))
    order = list(range(5)); prepared = starts(order); parent = None; parent_sha = None; entries = []
    for arm, budget in p.IDENTITIES:
        selected = list(reversed(order)) if arm == 'magnitude' else list(order)
        score = dict(value=1e-6, cubic_response=1e-6*18**4, steps=18)
        stats = dict(seed_evaluations=4)
        result = dict(order=selected, score=score, status='complete', evaluations=4, accepted_moves=0, stats=stats)
        if arm in ('canonical', 'magnitude'):
            result['stats'] = dict(seed_evaluations=1)
        elif arm == 'window':
            result['stats'] = dict(window_status='whole_call_limit', parent_method='magnitude', parent_evaluations=1,
                raw_candidate_inventory_complete=True)
            result['status'] = 'cpu_limit'
        elif arm == 'annealing':
            raw = dict(order=selected, score=score, native_score=score,
                       witness=dict(start_scores_completed=4, chain_proposals=[1, 1, 1, 1]))
            result.update(eligible=True, raw_julia_candidate=raw)
            result['stats'] = dict(four_chain_coverage=True, native_parent_revalidated=True, native_verifications=1)
        identity = dict(case=CASE, repetition=1, arm=arm, cpu_budget=budget)
        base = dict(identity, n=14, electrons=4, model_json_sha256=model_sha, model_toml_sha256=model_toml_sha,
            protocol_sha256=protocol_sha, source_freeze_sha256=source_sha, parent_sha256=parent_sha if arm in p.k.ARMS else None)
        own = dict(raw_CPU_seconds=1., charged_CPU_seconds=1.011, wall_seconds=2.)
        zero = {key: 0. for key in own}
        prefix = own if arm == 'parent' else parent['costs']['total'] if arm in p.k.ARMS else zero
        continuation = zero if arm == 'parent' else own
        costs = dict(prefix=prefix, continuation=continuation,
                     total={key: prefix[key]+continuation[key] for key in own})
        label = f'fixture-{arm}-{budget}'; raw_path = f'workers/{label}.json'
        cpu_cap = 8 if arm not in p.k.ARMS else budget-prefix['charged_CPU_seconds']
        wall_cap = 40 if arm not in p.k.ARMS else int(5*budget-prefix['wall_seconds'])
        supplied = None if arm not in p.k.ARMS else prepared['orders'] if arm in ('full', 'annealing') else [order]
        worker = dict(base, schema='molecular-worker-v1', result=result, refusal=None, starts=supplied,
            restarts=prepared if arm == 'parent' else None, observed_raw_witnesses=[],
            process_CPU_before_output=.5, guarded_CPU_cap=cpu_cap, guarded_wall_cap=wall_cap,
            physical_error_read=False, response_is_error_bound=False)
        worker_sha = put(root, raw_path, worker)
        command = ['/pinned/python', '/original/source/numerical_audit/molecular_generalization/construct.py',
            '--data', '/original/data', '--worker', '--case', CASE, '--repetition', '1', '--arm', arm,
            '--budget', str(budget), '--raw', '/original/data/'+raw_path, '--cpu-cap', str(cpu_cap), '--wall-cap', str(wall_cap)]
        plan = dict(schema='molecular-job-plan-v1', phase='constructors', label=label, command=command,
            identity=identity, expected_raw=raw_path, cpu_cap=cpu_cap, wall_cap=wall_cap, disk_mib=512,
            source_freeze_sha256=source_sha, guard_sha256=p.k.sha(root/guard), producer_root='/original/source',
            producer_data_root='/original/data')
        receipt = dict(command=command, budget=dict(cpu_seconds=cpu_cap, wall_seconds=wall_cap, disk_mib=512),
            process_measurements=dict(user_cpu_seconds=.7, system_cpu_seconds=.3), elapsed=2., samples=[],
            settings=dict(CPUQuota='25%', CPUWeight='5', IOWeight='10', Nice='19', MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0'),
            environment={key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
            exit_code=0, termination=None, guard_error=None, cleanup_errors=[])
        plan_path = f'jobs/{label}.PLAN.json'; receipt_path = f'logs/{label}.resources.json'
        record = dict(base, schema='molecular-constructor-v1', order=selected, score=score['value'],
            cubic_response=score['cubic_response'], steps=score['steps'], score_kind=p.k.SCORE_KIND,
            result=result, eligible=True, controller_fallback=False, refusal=None, costs=costs,
            own_plan=dict(path=plan_path, sha256=put(root, plan_path, plan)),
            own_receipt=dict(path=receipt_path, sha256=put(root, receipt_path, receipt)),
            worker_output_path=raw_path, worker_output_sha256=worker_sha,
            restarts=prepared if arm == 'parent' else None, physical_error_read=False, response_is_error_bound=False)
        relative = f'constructors/{CASE}.rep1.{arm}.B{budget}.json'; sha = put(root, relative, record)
        entries.append(dict(identity, path=relative, sha256=sha))
        if arm == 'parent':
            parent, parent_sha = record, sha
    put(root, 'PARENTS_FROZEN.json', dict(info, entries=[entry for entry in entries if entry['arm'] == 'parent']))
    put(root, 'CONSTRUCTORS_FROZEN.json', dict(info, entries=[entry for entry in entries if entry['arm'] != 'parent']))
    return p.context(root, preflight=True)


class MetadataTests(unittest.TestCase):
    def test_real_byte_bindings_roster_raw_claims_and_relocated_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); context = bundle(root)
            inventory, panel, witnesses = p.validate_case(*context, CASE, emitted=False)
            self.assertEqual(len(panel['candidates']), 36)
            self.assertEqual(len(panel['constructor_artifacts']), 9)
            self.assertEqual(len(witnesses['witnesses']), 2)
            self.assertTrue(all('parent_sha256' not in row for row in panel['candidates'] if row['arm'] in p.k.CHEAP))
            p.write_manifest(root/f'panels/{CASE}.toml', panel, ('constructor_artifacts', 'candidates'))
            p.write_manifest(root/f'witnesses/{CASE}.toml', witnesses, ('constructor_artifacts', 'witnesses'))
            report = p.validate(root, preflight=True)
            self.assertTrue(report['raw_inventory_independently_confirmed'])
            self.assertEqual(report['cases_metadata'][0]['constructor_count'], 9)
            # Missing physics produces an incomplete analysis, never synthesized errors.
            analysis = a.analyze(root, preflight=True)
            self.assertFalse(analysis['complete'])
            self.assertTrue(all(row['difference'] is None for row in analysis['named_pairs']))
            with self.assertRaises(ValueError):
                p.validate(root) # Explicit old-only mode is required.
            bad = copy.deepcopy(witnesses); bad['witnesses'].pop(); bad['expected_witnesses'] -= 1
            path = root/f'witnesses/{CASE}.toml'; path.unlink()
            p.write_manifest(path, bad, ('constructor_artifacts', 'witnesses'))
            with self.assertRaises(ValueError):
                p.validate(root, preflight=True)

    def test_empty_witnesses_explicit_freshness_and_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); root.joinpath('real').write_text('x')
            artifacts = p.Artifacts(root)
            with self.assertRaises(ValueError):
                artifacts.raw('../escape')
            root.joinpath('alias').symlink_to(root/'real')
            with self.assertRaises(ValueError):
                artifacts.raw('alias')
            target = root/'empty.toml'; manifest = dict(expected_witnesses=0, witnesses=[])
            p.write_manifest(target, manifest, ('witnesses',))
            self.assertEqual(artifacts.toml('empty.toml'), manifest)
            with self.assertRaises(ValueError):
                p.write_manifest(target, manifest, ('witnesses',))
            root.joinpath('duplicate.json').write_text('{"score": 1, "score": 2}')
            with self.assertRaises(ValueError):
                artifacts.json('duplicate.json')

    def test_unknown_cpu_and_sampled_cpu_are_not_zero(self):
        receipt = dict(process_measurements={}, elapsed=3., samples=[])
        self.assertIsNone(p.measured_receipt(receipt)['charged_CPU_seconds'])
        receipt.update(process_measurements=dict(user_cpu_seconds=.2, system_cpu_seconds=.3),
                       samples=[dict(service='CPUUsageNSec=900000000\n')])
        self.assertEqual(p.measured_receipt(receipt)['charged_CPU_seconds'], .9)
        self.assertEqual(p.measured_receipt(receipt)['raw_CPU_seconds'], .5)

    def test_two_named_primary_rules_keep_opposite_repetitions_and_floor_differences(self):
        cases = ['a', 'b']; physical = {}; eligible = set()
        for case in cases:
            for rep in (1, 2):
                for arm in ('full', 'window', 'annealing'):
                    key = (case, arm, rep, 64, 1.35, 8)
                    physical[key] = dict(projective_error=1e-8 if arm == 'full' else 1e-9, numerical_floor=False)
                    eligible.add(key)
        pairs = a.paired_grid(cases, 2, physical, eligible)
        self.assertEqual(len(pairs), 64)
        result = a.primary_decision(pairs, cases, complete=True)
        self.assertEqual(result['window_minus_full']['decision'], 'candidate_advantage')
        physical[('a', 'window', 1, 64, 1.35, 8)]['projective_error'] = 1.1e-8
        pairs = a.paired_grid(cases, 2, physical, eligible)
        self.assertEqual(a.primary_decision(pairs, cases, complete=True)['window_minus_full']['decision'], 'mixed_or_unresolved')
        self.assertEqual(a.primary_decision(pairs, cases, complete=True)['annealing_minus_full']['decision'], 'candidate_advantage')
        eligible.remove(('a', 'window', 1, 64, 1.35, 8))
        pairs = a.paired_grid(cases, 2, physical, eligible)
        self.assertEqual(a.primary_decision(pairs, cases, complete=True)['window_minus_full']['decision'], 'incomplete')
        self.assertIsNotNone(next(row['difference'] for row in pairs if row['question'] == 'window_minus_full' and row['time'] == 1.35 and row['cpu_budget'] == 64))
        self.assertFalse(a.comparison(0., 0.)['verified'])

    def test_cost_accounting_parent_once_and_unknown_inventory(self):
        own = dict(raw_CPU_seconds=1., charged_CPU_seconds=1.1, wall_seconds=2.)
        zero = {key: 0. for key in own}; records = {}
        for arm, budget in p.IDENTITIES:
            prefix = own if arm == 'parent' or arm in p.k.ARMS else zero
            continuation = zero if arm == 'parent' else own
            records[(arm, budget)] = dict(case='x', repetition=1, arm=arm, cpu_budget=budget,
                costs=dict(prefix=prefix, continuation=continuation, total={key: prefix[key]+continuation[key] for key in own}), analysis_eligible=True)
        totals = a.cost_totals({'x': dict(records=records)}, expected_count=9)
        self.assertAlmostEqual(totals['actual']['charged_CPU_seconds'], 9*1.1)
        self.assertAlmostEqual(totals['comparison_attributed']['charged_CPU_seconds'], 14*1.1)
        totals = a.cost_totals({}, expected_count=9)
        self.assertIsNone(totals['actual']['charged_CPU_seconds'])
        self.assertFalse(totals['inventory_complete'])


if __name__ == '__main__':
    unittest.main()
