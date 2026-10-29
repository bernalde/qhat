"""Construct without physics; bind the result to complete owning-cgroup costs.

The worker and post-process finalizer are separate: only the latter can read
the actual process/child teardown receipt. A failed worker never gets a new B.
"""
from __future__ import annotations

import time
BOOT_WALL = time.monotonic()
import argparse
from dataclasses import asdict
import hashlib
import math
import random
import shutil
import common as k


def header(root):
    frozen = k.frozen_protocol(root)
    config = k.validate_protocol(frozen['configuration'])
    inputs = k.read(root/'INPUTS_FROZEN.json')
    k.require(inputs['status'] == 'complete' and not inputs['refusals'], 'Complete input freeze required')
    mode = frozen.get('mode', 'main')
    k.require(mode in ('main', 'preflight'), 'Known campaign mode')
    cases = [config['preflight']['case']] if mode == 'preflight' else [v['case'] for v in frozen['registry']['cases']]
    k.require(set(inputs['successes']) == set(cases), 'Exact declared successful input grid')
    result = dict(mode=mode, repetitions=1 if mode == 'preflight' else config['repetitions'], cases=cases,
        protocol_sha256=k.sha(k.HERE/'PROTOCOL.json'), source_freeze_sha256=k.sha(root/'INPUT_PROTOCOL_FROZEN.json'),
        protocol_path='source_snapshot/numerical_audit/cost_matched_continuation/PROTOCOL.json',
        source_freeze_path='INPUT_PROTOCOL_FROZEN.json')
    k.require(inputs['protocol_sha256'] == result['protocol_sha256'] and
              inputs['source_freeze_sha256'] == result['source_freeze_sha256'], 'Input freeze identity')
    return result, config, inputs


def identity(root, case, repetition, arm, budget):
    info, config, inputs = header(root)
    k.require(case in info['cases'] and 1 <= repetition <= info['repetitions'], 'Frozen case/repetition')
    k.require((arm == 'parent' and budget == 0) or
              (arm in k.ARMS and budget in k.BUDGETS), 'Frozen arm and total budget')
    folder = root/'inputs'/case
    for name in ('model.json', 'model.toml'):
        k.require(k.sha(folder/name) == inputs['successes'][case][name], 'Frozen input bytes changed')
    model = k.read(folder/'model.json')
    k.require(model['case'] == case, 'Input case mismatch')
    return info, config, model


def restarts(order, case, repetition):
    """Four fixed modest perturbations, chosen before either continuation."""
    material = (case+':'+str(repetition)+':cost-matched-v1').encode('utf-8')
    seed = int.from_bytes(hashlib.sha256(material).digest()[:8], 'big')
    rng = random.Random(seed)
    k.require(len(order) >= 5, 'Interior perturbation requires at least five units')
    orders = [list(order)]
    for count in (1, 2, 4):
        for _ in range(100):
            trial = list(order)
            for _ in range(count):
                source, target = rng.sample(range(1, len(order)-1), 2)
                trial.insert(target, trial.pop(source))
            if trial not in orders:
                orders.append(trial)
                break
        else:
            raise ValueError('Could not construct four distinct frozen starts')
    return dict(seed=seed, orders=orders, interior_perturbation_counts=[0, 1, 2, 4])


def canonical_path(root, case, repetition, arm, budget):
    return (root/'parents'/case/('rep'+str(repetition)+'.json') if arm == 'parent' else
            root/'orders'/case/('rep'+str(repetition))/(arm+'_B'+str(budget)+'.json'))


def parent_record(root, case, repetition):
    index = k.read(root/'PARENTS_FROZEN.json')
    info, _, _ = header(root)
    k.require(all(index[key] == info[key] for key in info), 'Parent index header identity')
    expected = {(name, rep) for name in info['cases'] for rep in range(1, info['repetitions']+1)}
    k.require({(v['case'], v['repetition']) for v in index['entries']} == expected and
              len(index['entries']) == len(expected), 'Freeze every parent before any continuation')
    rows = [v for v in index['entries'] if (v['case'], v['repetition']) == (case, repetition)]
    path = canonical_path(root, case, repetition, 'parent', 0)
    k.require(len(rows) == 1 and rows[0]['path'] == str(path.relative_to(root)) and
              k.sha(path) == rows[0]['sha256'], 'Frozen parent bytes changed')
    value = k.read(path)
    k.require(value['eligible'] and value['score'] is not None, 'Parent must be checked and within its cap')
    return value, rows[0]['sha256']


def worker(root, case, repetition, arm, budget, cpu_cap, wall_cap, output, julia):
    info, config, model = identity(root, case, repetition, arm, budget)
    k.require(not output.exists(), 'Fresh worker output')
    parent, parent_sha = (None, None) if arm == 'parent' else parent_record(root, case, repetition)
    prefix = 0. if parent is None else parent['costs']['total']['charged_CPU_seconds']
    prefix_wall = 0. if parent is None else parent['costs']['total']['wall_seconds']
    expected_cpu = config['parent']['CPU_seconds'] if parent is None else budget-prefix
    expected_wall = config['parent']['wall_seconds'] if parent is None else math.floor(5*budget-prefix_wall)
    k.require(cpu_cap == expected_cpu and wall_cap == expected_wall, 'Guard/API complete remaining allowance differs')
    units = k.prior.blocks(model)
    ordering = k.bootstrap()
    starts = None if parent is None else (parent['restarts']['orders'] if arm == 'full' else [parent['order']])
    available = cpu_cap-time.process_time()-config['accounting']['CPU_cleanup_reserve']
    wall_available = wall_cap-(time.monotonic()-BOOT_WALL)-config['accounting']['wall_cleanup_reserve']
    result = None
    refusal = None
    try:
        k.require(available > 0 and wall_available > 0, 'Startup exhausted the continuation allowance')
        options = {}
        if arm == 'window':
            from qhat.analysis.window_ordering import WindowOptions
            fields = {key: value for key, value in config['window'].items() if key in
                ('width', 'max_windows', 'max_sweeps', 'parent_method', 'parent_fraction',
                 'max_factor_entries', 'max_flat_entries', 'max_actions', 'max_atoms')}
            options = dict(method='window', window_options=WindowOptions(**fields,
                julia_executable=julia, wall_seconds=wall_available))
        limits = config['parent'] if arm == 'parent' else config['native']
        trained = ordering.train_ordering(units, model['n'], (1 << model['electrons'])-1,
            time=config['training_time'], rotation_budget=config['training_rotation_multiplier']*len(model['coefficients']),
            cpu_seconds=available, max_evaluations=limits['max_evaluations'], max_moves=limits['max_moves'],
            max_entries=config['native']['max_entries'], max_cached_actions=config['native']['max_cached_actions'],
            starts=starts, **options)
        result = asdict(trained)
        result['order'] = list(result['order'])
        k.validate_order(result['order'], len(units))
        if result['score'] is not None:
            k.require(all(math.isfinite(result['score'][key]) and result['score'][key] >= 0
                for key in ('value', 'cubic_response')), 'Finite fully scored result')
    except Exception as exc:
        refusal = dict(type=type(exc).__name__, message=str(exc))
    k.write_new(output, dict(schema='cost-matched-worker-v1', case=case, repetition=repetition,
        arm=arm, total_CPU_budget=budget, protocol_sha256=info['protocol_sha256'],
        source_freeze_sha256=info['source_freeze_sha256'], parent_sha256=parent_sha,
        model_json_sha256=k.sha(root/'inputs'/case/'model.json'),
        model_toml_sha256=k.sha(root/'inputs'/case/'model.toml'),
        result=result, refusal=refusal, starts=starts, API_CPU_allowance=available,
        API_wall_allowance=wall_available, guarded_CPU_cap=cpu_cap, guarded_wall_cap=wall_cap,
        process_CPU_before_output=time.process_time(), wall_before_output=time.monotonic()-BOOT_WALL,
        physical_error_read=False, response_is_error_bound=False))
    print('CONSTRUCTED', case, repetition, arm, budget, None if result is None else result['status'], flush=True)


def copy_evidence(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        k.require(k.sha(source) == k.sha(target), 'Existing copied evidence differs')
    else:
        with source.open('rb') as stream, target.open('xb') as destination:
            shutil.copyfileobj(stream, destination)


def measured_cost(root, campaign, label):
    final_name, receipt_name = 'jobs/'+label+'.FINAL.json', 'logs/'+label+'.resources.json'
    for name in ('jobs/'+label+'.PLAN.json', final_name, receipt_name):
        copy_evidence(campaign/name, root/'campaign'/name)
    final = k.read(campaign/final_name)
    receipt = k.read(campaign/receipt_name)
    k.require(final['phase'] == 'constructors' and final['label'] == label and
        final['receipt_sha256'] == k.sha(campaign/receipt_name), 'Constructor receipt identity')
    return dict(raw_CPU_seconds=receipt.get('final_accounting', {}).get('CPU_seconds'),
        charged_CPU_seconds=final['charged_CPU_seconds'], wall_seconds=final['absolute_guarded_wall_seconds'],
        label=label, final_path='campaign/'+final_name, final_sha256=k.sha(campaign/final_name),
        receipt_path='campaign/'+receipt_name, receipt_sha256=k.sha(campaign/receipt_name)), final['exit_code']


def finalize(root, case, repetition, arm, budget, raw, campaign, label):
    info, config, model = identity(root, case, repetition, arm, budget)
    path = canonical_path(root, case, repetition, arm, budget)
    k.require(not path.exists(), 'One canonical construction per identity; preserve attempts')
    parent, parent_sha = (None, None) if arm == 'parent' else parent_record(root, case, repetition)
    cost, exit_code = measured_cost(root, campaign, label)
    record = k.read(raw) if raw.exists() else None
    if record is not None:
        k.require(all(record[key] == value for key, value in dict(case=case, repetition=repetition,
            arm=arm, total_CPU_budget=budget, parent_sha256=parent_sha,
            protocol_sha256=info['protocol_sha256'], source_freeze_sha256=info['source_freeze_sha256'],
            model_json_sha256=k.sha(root/'inputs'/case/'model.json'),
            model_toml_sha256=k.sha(root/'inputs'/case/'model.toml')).items()), 'Worker result identity changed')
    result = None if record is None else record['result']
    fallback = result is None or result.get('score') is None
    if fallback:
        order = list(range(len(model['blocks']))) if parent is None else parent['order']
        score = None if parent is None else parent['score']
    else:
        order, score = result['order'], result['score']
    k.validate_order(order, len(model['blocks']))
    fields = ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds')
    zero = {key: 0. for key in fields}
    prefix, continuation = (cost, zero) if parent is None else (parent['costs']['prefix'], cost)
    total = {key: None if prefix[key] is None or continuation[key] is None else prefix[key]+continuation[key]
             for key in fields}
    ceiling = config['parent']['CPU_seconds'] if parent is None else budget
    wall_ceiling = config['parent']['wall_seconds'] if parent is None else 5*budget
    eligible = (exit_code == 0 and not fallback and score is not None and
                total['raw_CPU_seconds'] is not None and total['raw_CPU_seconds'] <= ceiling and
                total['charged_CPU_seconds'] <= ceiling and total['wall_seconds'] <= wall_ceiling)
    value = dict(schema='cost-matched-constructor-v1', case=case, repetition=repetition, arm=arm,
        total_CPU_budget=budget, order=order, score=score, result=result,
        parent_sha256=parent_sha, protocol_sha256=info['protocol_sha256'],
        source_freeze_sha256=info['source_freeze_sha256'],
        model_json_sha256=k.sha(root/'inputs'/case/'model.json'),
        model_toml_sha256=k.sha(root/'inputs'/case/'model.toml'),
        costs=dict(prefix=prefix, continuation=continuation, total=total),
        eligible=bool(eligible), guard_exit_code=exit_code, controller_fallback=bool(fallback),
        worker_output_path=str(raw.relative_to(root)), worker_output_sha256=k.sha(raw) if raw.exists() else None,
        physical_error_read=False, response_is_error_bound=False)
    if arm == 'parent':
        value['restarts'] = restarts(order, case, repetition)
    k.write_new(path, value)
    return value


def freeze(root, parents):
    info, _, _ = header(root)
    if not parents:
        k.require((root/'PARENTS_FROZEN.json').exists(), 'Parents frozen before continuation freeze')
    entries = []
    for case in info['cases']:
        for rep in range(1, info['repetitions']+1):
            choices = [('parent', 0)] if parents else [(arm, budget) for arm in k.ARMS for budget in k.BUDGETS]
            for arm, budget in choices:
                path = canonical_path(root, case, rep, arm, budget)
                record = k.read(path)
                k.require((record['case'], record['repetition'], record['arm'], record['total_CPU_budget']) ==
                          (case, rep, arm, budget), 'Canonical constructor identity')
                entry = dict(case=case, repetition=rep, path=str(path.relative_to(root)), sha256=k.sha(path))
                if not parents:
                    entry.update(arm=arm, total_CPU_budget=budget)
                entries.append(entry)
    name = 'PARENTS_FROZEN.json' if parents else 'CONSTRUCTORS_FROZEN.json'
    k.write_new(root/name, dict(info, schema='cost-matched-'+('parent' if parents else 'constructor')+'-freeze-v1', entries=entries))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=k.Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--worker', action='store_true')
    group.add_argument('--finalize', action='store_true')
    group.add_argument('--freeze', choices=('parents', 'continuations'))
    parser.add_argument('--case')
    parser.add_argument('--repetition', type=int)
    parser.add_argument('--arm', choices=('parent',)+k.ARMS)
    parser.add_argument('--budget', type=int, default=0)
    parser.add_argument('--cpu-cap', type=float)
    parser.add_argument('--wall-cap', type=int)
    parser.add_argument('--raw', type=k.Path)
    parser.add_argument('--campaign', type=k.Path)
    parser.add_argument('--label')
    parser.add_argument('--julia-executable')
    args = parser.parse_args()
    if args.freeze:
        freeze(args.root, args.freeze == 'parents')
    elif args.worker:
        worker(args.root, args.case, args.repetition, args.arm, args.budget,
               args.cpu_cap, args.wall_cap, args.raw, args.julia_executable)
    else:
        finalize(args.root, args.case, args.repetition, args.arm, args.budget,
                 args.raw, args.campaign, args.label)
