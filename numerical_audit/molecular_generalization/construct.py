"""Fresh constructor workers and metadata-only post-teardown finalization.

The OS receipt owns complete startup, all children, checking, serialization and
teardown costs. A missed scientific cap retains an ineligible scored incumbent.
"""
from __future__ import annotations
import time
BOOT_WALL = time.monotonic()
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import hashlib
import math
import os
from pathlib import Path
import random
import common as k

admission = k.load_module('_molecular_execution_admission',
    k.ROOT/'numerical_audit/cost_matched_continuation/admission.py')


def identity(data, case, rep, arm, budget, *, producer=False):
    h, config = k.header(data, producer=producer)
    k.require(case in h['cases'] and (rep, arm, budget) in k.identities(h['repetitions']),
              'Frozen constructor identity')
    inputs = k.read(data/'INPUTS_FROZEN.json')
    k.require(inputs['status'] == 'complete' and not inputs['refusals'] and
              set(inputs['successes']) == set(h['cases']), 'Complete input grid required')
    for field in ('protocol_sha256', 'source_freeze_sha256', 'mode'):
        k.require(inputs[field] == h[field], 'Frozen input/source join')
    folder = data/'inputs'/case
    for name, digest in inputs['successes'][case].items():
        k.require(k.sha(folder/name) == digest, 'Frozen input bytes changed')
    model = k.read(folder/'model.json')
    k.require(model['case'] == case, 'Model identity')
    return h, config, model


def restarts(order, case, rep):
    """Same four frozen modest interior perturbations as retained full control."""
    seed = int.from_bytes(hashlib.sha256(f'{case}:{rep}:cost-matched-v1'.encode()).digest()[:8], 'big')
    rng, orders = random.Random(seed), [list(order)]
    k.require(len(order) >= 5, 'Four distinct perturbations need at least five units')
    for count in (1, 2, 4):
        for _ in range(100):
            candidate = list(order)
            for _ in range(count):
                source, target = rng.sample(range(1, len(order)-1), 2)
                candidate.insert(target, candidate.pop(source))
            if candidate not in orders:
                orders.append(candidate)
                break
        else:
            raise ValueError('Could not construct distinct frozen starts')
    return dict(seed=seed, orders=orders, interior_perturbation_counts=[0, 1, 2, 4])


def parent(data, case, rep):
    index = k.read(data/'PARENTS_FROZEN.json')
    h, _ = k.header(data)
    expected = {(c, r) for c in h['cases'] for r in range(1, h['repetitions']+1)}
    k.require(len(index['entries']) == len(expected) and
              {(e['case'], e['repetition']) for e in index['entries']} == expected, 'Complete parent grid')
    path = k.constructor_path(data, case, rep, 'parent', 0)
    entry = next(e for e in index['entries'] if (e['case'], e['repetition']) == (case, rep))
    k.require(entry['path'] == str(path.relative_to(data)) and entry['sha256'] == k.sha(path),
              'Immutable measured parent')
    value = k.read(path)
    k.require(all(value[f] == h[f] for f in ('protocol_sha256', 'source_freeze_sha256')), 'Parent/source identity')
    return value, entry['sha256']


def allowances(config, arm, budget, measured_parent=None):
    if arm == 'parent':
        return config['parent']['CPU_seconds'], config['parent']['wall_seconds']
    if arm in k.CHEAP:
        return config['cheap']['CPU_ceiling'], config['cheap']['wall_ceiling']
    cost = measured_parent['costs']['total']
    k.require(cost['charged_CPU_seconds'] is not None and cost['wall_seconds'] is not None,
              'Measured parent costs required')
    return budget-cost['charged_CPU_seconds'], math.floor(5*budget-cost['wall_seconds'])


@contextmanager
def observe_window(output, artifact_base=None):
    """Observe authenticated replies and actual native checks, without new scores.

    This study-only observer restores both functions and changes no argument or
    return value. Original numerical source bytes and selection logic are intact.
    """
    from qhat.analysis import window_ordering as window
    original = window._reply
    restored = []
    def observed(*args, **kwargs):
        candidate, reply = original(*args, **kwargs)
        context = args[2]
        raw = Path(args[0]).read_bytes()
        score = dict(value=float(reply['score']), cubic_response=float(reply['cubic_response']),
                     steps=int(reply['steps']))
        witness = dict(order=list(candidate), score=score, witness=reply,
                       checkpoint_sha256=hashlib.sha256(raw).hexdigest(),
                       native_refusal='No subsequent native score observed')
        if artifact_base is not None:
            checkpoint = Path(str(artifact_base)+'.checkpoint.toml')
            request = Path(str(artifact_base)+'.request.toml')
            request_raw = Path(args[0]).with_name('request.toml').read_bytes()
            k.require(len(raw) <= 64*1024 and len(request_raw) <= 8*1024**2,
                      'Bounded retained owned window wire artifacts')
            with checkpoint.open('xb') as stream:
                stream.write(raw)
            with request.open('xb') as stream:
                stream.write(request_raw)
            witness['retained_checkpoint_name'] = checkpoint.name
            witness['retained_request_name'] = request.name
            witness['request_sha256'] = hashlib.sha256(request_raw).hexdigest()
        output.append(witness)
        method = context.score
        def observed_score(order, *score_args, **score_kwargs):
            try:
                value = method(order, *score_args, **score_kwargs)
            except Exception as exc:
                if tuple(order) == tuple(candidate):
                    witness['native_refusal'] = dict(type=type(exc).__name__, message=str(exc))
                raise
            if tuple(order) == tuple(candidate):
                witness['native_score'] = asdict(value)
                witness.pop('native_refusal', None)
            return value
        context.score = observed_score
        restored.append((context, method))
        return candidate, reply
    window._reply = observed
    try:
        yield
    finally:
        window._reply = original
        for context, method in restored:
            context.score = method


def worker(data, case, rep, arm, budget, cpu_cap, wall_cap, raw, julia):
    h, config, model = identity(data, case, rep, arm, budget, producer=True)
    k.require(not raw.exists(), 'Fresh one-attempt constructor output')
    p, parent_sha = parent(data, case, rep) if arm in k.ARMS else (None, None)
    expected_cpu, expected_wall = allowances(config, arm, budget, p)
    k.require(math.isfinite(cpu_cap) and cpu_cap > 0 and abs(cpu_cap-expected_cpu) <= 1e-12
              and wall_cap == expected_wall, 'Remaining whole-call allowance mismatch')
    os.environ['QHAT_JULIA_EXECUTABLE'] = julia
    result = refusal = starts = restart_plan = None
    witnesses = []
    available = wall_available = None
    try:
        k.require(p is None or p['eligible'], 'Dependent arm refused: parent ineligible')
        ordering, units = k.bootstrap(), k.blocks(model)
        starts = None if p is None else (p['restarts']['orders'] if arm in ('full', 'annealing') else [p['order']])
        available = cpu_cap-time.process_time()-config['accounting']['CPU_cleanup_reserve']
        wall_available = wall_cap-(time.monotonic()-BOOT_WALL)-config['accounting']['wall_cleanup_reserve']
        k.require(available > 0 and wall_available > 0, 'Startup exhausted allowance')
        rotation_cap = 30*len(model['coefficients'])
        if arm in k.CHEAP:
            order = list(range(len(units)))
            if arm == 'magnitude':
                order.sort(key=lambda i: (-sum(abs(c) for c in units[i].values()), i))
            # One native scored fixed order, no insertion search.
            trained = ordering.train_ordering(units, model['n'], (1 << model['electrons'])-1,
                time=1., rotation_budget=rotation_cap, cpu_seconds=available,
                max_evaluations=2, max_moves=0, max_entries=50000, max_cached_actions=100000,
                starts=[order])
            result = asdict(trained)
        elif arm == 'annealing':
            import annealing
            a = config['annealing']
            settings = dict(calibrated=True, calibration_sha256=a['calibration_sha256'],
                theta_start=a['temperature_start'], theta_end=a['temperature_end'],
                max_evaluations=a['max_evaluations'], max_proposals=a['max_proposals'],
                max_actions=a['max_actions'], cache_bytes=a['max_cache_bytes'])
            checkpoint = raw.with_suffix('.annealing.toml')
            trained = annealing.train_annealing(units, model['n'], (1 << model['electrons'])-1,
                parent_order=p['order'], starts=starts, time=1., rotation_budget=rotation_cap,
                parent_score=dict(value=p['score'], cubic_response=p['cubic_response'], steps=p['steps']),
                settings=settings, seed_material=f'{case}:{rep}', cpu_seconds=available,
                wall_seconds=wall_available, max_entries=50000, max_cached_actions=100000,
                checkpoint_path=checkpoint)
            result = asdict(trained)
        else:
            options = {}
            if arm == 'window':
                from qhat.analysis.window_ordering import WindowOptions
                fields = {key: value for key, value in config['window'].items() if key in (
                    'width', 'max_windows', 'max_sweeps', 'parent_method', 'parent_fraction',
                    'max_factor_entries', 'max_flat_entries', 'max_actions', 'max_atoms')}
                options = dict(method='window', window_options=WindowOptions(**fields,
                    julia_executable=julia, wall_seconds=wall_available))
            limits = config['parent'] if arm == 'parent' else config['native']
            def train():
                return ordering.train_ordering(units, model['n'], (1 << model['electrons'])-1,
                    time=1., rotation_budget=rotation_cap, cpu_seconds=available,
                    max_evaluations=limits['max_evaluations'], max_moves=limits['max_moves'],
                    max_entries=50000, max_cached_actions=100000, starts=starts, **options)
            if arm == 'window':
                with observe_window(witnesses, raw.with_suffix('.window')):
                    trained = train()
            else:
                trained = train()
            result = asdict(trained)
        result['order'] = k.validate_order(result['order'], len(units))
        if result['score'] is not None:
            k.finite_score(result['score'])
        if arm == 'window':
            k.require(len(witnesses) <= 1, 'One owned window child')
            if witnesses:
                witness = witnesses[0]
                if witness['order'] == result['order']:
                    witness['native_score'] = result['score']
                    witness.pop('native_refusal', None)
                result['stats']['raw_julia_candidate'] = witness
            result['stats']['raw_candidate_inventory_complete'] = True
        if arm == 'parent' and result['score'] is not None:
            restart_plan = restarts(result['order'], case, rep)
    except Exception as exc:
        refusal = dict(type=type(exc).__name__, message=str(exc))
    k.write_new(raw, dict(schema='molecular-worker-v1', case=case, repetition=rep, arm=arm, cpu_budget=budget,
        n=model['n'], electrons=model['electrons'], protocol_sha256=h['protocol_sha256'],
        source_freeze_sha256=h['source_freeze_sha256'], parent_sha256=parent_sha,
        model_json_sha256=k.sha(data/'inputs'/case/'model.json'),
        model_toml_sha256=k.sha(data/'inputs'/case/'model.toml'), result=result, refusal=refusal,
        starts=starts, restarts=restart_plan, observed_raw_witnesses=witnesses,
        API_CPU_allowance=available, API_wall_allowance=wall_available,
        guarded_CPU_cap=cpu_cap, guarded_wall_cap=wall_cap,
        process_CPU_before_output=time.process_time(), wall_before_output=time.monotonic()-BOOT_WALL,
        physical_error_read=False, response_is_error_bound=False))
    print('CONSTRUCTED', case, rep, arm, budget, None if result is None else result['status'], flush=True)


def cost_from_receipt(receipt):
    measurements = receipt.get('process_measurements', {})
    values = [measurements.get(k) for k in ('user_cpu_seconds', 'system_cpu_seconds')]
    raw = sum(values) if all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in values) else None
    sampled = []
    for sample in receipt.get('samples', []):
        fields = dict(line.split('=', 1) for line in sample.get('service', '').splitlines() if '=' in line)
        value = fields.get('CPUUsageNSec', '')
        if value.isdigit():
            sampled.append(int(value)/1e9)
    # GNU time rounds two fields to .01s. Keep a conservative .011s margin;
    # sampled cgroup values may be lower bounds and are never terminal proof.
    charged = None if raw is None else max([raw+.011]+sampled)
    wall = receipt.get('elapsed')
    k.require(type(wall) in (int, float) and math.isfinite(wall) and wall >= 0, 'Measured absolute guarded wall')
    return dict(raw_CPU_seconds=raw, charged_CPU_seconds=charged, wall_seconds=wall)


def finalize(data, case, rep, arm, budget, raw, plan_path, receipt_path):
    h, config, model = identity(data, case, rep, arm, budget)
    target = k.constructor_path(data, case, rep, arm, budget)
    k.require(not target.exists(), 'No fresh-budget retry for same constructor identity')
    plan, receipt = k.read(plan_path), k.read(receipt_path)
    owner = dict(case=case, repetition=rep, arm=arm, cpu_budget=budget)
    k.require(plan['identity'] == owner and plan['phase'] == 'constructors'
              and plan['source_freeze_sha256'] == h['source_freeze_sha256'], 'Owning plan identity')
    k.require(receipt['command'] == plan['command'] and receipt['budget'] == dict(
        cpu_seconds=plan['cpu_cap'], wall_seconds=plan['wall_cap'], disk_mib=plan['disk_mib']), 'Plan/receipt join')
    k.require(plan['expected_raw'] == str(raw.relative_to(data)), 'Owning raw result path')
    record = k.read(raw) if raw.is_file() else None
    p, psha = parent(data, case, rep) if arm in k.ARMS else (None, None)
    expected_cpu, expected_wall = allowances(config, arm, budget, p)
    k.require(abs(plan['cpu_cap']-expected_cpu) <= 1e-12 and plan['wall_cap'] == expected_wall, 'Complete allowance')
    if record is not None:
        expected = dict(owner, n=model['n'], electrons=model['electrons'], parent_sha256=psha,
            protocol_sha256=h['protocol_sha256'], source_freeze_sha256=h['source_freeze_sha256'],
            model_json_sha256=k.sha(data/'inputs'/case/'model.json'), model_toml_sha256=k.sha(data/'inputs'/case/'model.toml'))
        k.require(all(record[key] == value for key, value in expected.items()), 'Raw constructor identity')
    result = None if record is None else record['result']
    fallback = result is None or result.get('score') is None
    if fallback:
        order = list(range(len(model['blocks']))) if p is None else p['order']
        score = None if p is None else dict(value=p['score'], cubic_response=p['cubic_response'], steps=p['steps'])
    else:
        order, score = result['order'], k.finite_score(result['score'])
    k.validate_order(order, len(model['blocks']))
    own_cost = cost_from_receipt(receipt)
    zero = {key: 0. for key in own_cost}
    if arm == 'parent':
        prefix, continuation = own_cost, zero
    elif arm in k.CHEAP:
        prefix, continuation = zero, own_cost
    else:
        prefix, continuation = p['costs']['total'], own_cost
    total = {key: None if prefix[key] is None or continuation[key] is None else prefix[key]+continuation[key]
             for key in own_cost}
    refusal = {'type': 'missing_worker'} if record is None else record['refusal']
    if arm == 'annealing':
        admitted = bool(result is not None and result.get('eligible') is True and refusal is None)
        execution = dict(admitted=admitted, category='annealing', reason=None if admitted else 'Incomplete/refused annealing')
    elif arm in k.CHEAP:
        admitted = not fallback and refusal is None and result['status'] == 'complete' and result['order'] == order
        execution = dict(admitted=bool(admitted), category='fixed_order', reason=None if admitted else 'Fixed order not checked')
    else:
        execution = admission.execution_admission(arm, result, None if p is None else p['order'], refusal)
    rest = None if record is None else record['restarts']
    if arm == 'parent':
        execution['admitted'] = bool(execution['admitted'] and rest is not None and
            len(rest['orders']) == 4 and rest['orders'][0] == order and len({tuple(o) for o in rest['orders']}) == 4)
    if p is not None:
        saved_starts = (p.get('restarts') or {}).get('orders')
        expected_starts = saved_starts if arm in ('full', 'annealing') else [p['order']]
        execution['admitted'] = bool(execution['admitted'] and p['eligible'] and record is not None
            and expected_starts is not None and record['starts'] == expected_starts)
    ceiling = config['parent']['CPU_seconds'] if arm == 'parent' else config['cheap']['CPU_ceiling'] if arm in k.CHEAP else budget
    wall_cap = config['parent']['wall_seconds'] if arm == 'parent' else config['cheap']['wall_ceiling'] if arm in k.CHEAP else 5*budget
    eligible = (not fallback and score is not None and execution['admitted'] and receipt['exit_code'] == 0
        and receipt.get('termination') is None and receipt.get('guard_error') is None and not receipt.get('cleanup_errors')
        and total['raw_CPU_seconds'] is not None and total['charged_CPU_seconds'] is not None
        and total['charged_CPU_seconds'] <= ceiling and total['wall_seconds'] <= wall_cap)
    if record is not None and own_cost['charged_CPU_seconds'] is not None:
        k.require(own_cost['charged_CPU_seconds'] >= record['process_CPU_before_output'], 'Receipt below worker CPU')
    def reference(path):
        k.require(path.resolve().is_relative_to(data.resolve()), 'Owning evidence contained inDATA')
        return dict(path=str(path.relative_to(data)), sha256=k.sha(path))
    value = dict(schema='molecular-constructor-v1', **owner, n=model['n'], electrons=model['electrons'],
        model_json_sha256=k.sha(data/'inputs'/case/'model.json'), model_toml_sha256=k.sha(data/'inputs'/case/'model.toml'),
        protocol_sha256=h['protocol_sha256'], source_freeze_sha256=h['source_freeze_sha256'],
        parent_sha256=psha, order=order, score=None if score is None else score['value'],
        cubic_response=None if score is None else score['cubic_response'], steps=None if score is None else score['steps'],
        score_kind=k.SCORE_KIND, eligible=bool(eligible), costs=dict(prefix=prefix, continuation=continuation, total=total),
        own_plan=reference(plan_path), own_receipt=reference(receipt_path),
        worker_output_path=str(raw.relative_to(data)), worker_output_sha256=k.sha(raw) if raw.is_file() else None,
        result=result, refusal=refusal, controller_fallback=bool(fallback), execution_admission=execution,
        restarts=rest, physical_error_read=False, response_is_error_bound=False)
    k.write_new(target, value)
    return value


def freeze(data, parents_only):
    h, _ = k.header(data)
    entries = []
    for case in h['cases']:
        for rep, arm, budget in k.identities(h['repetitions']):
            if parents_only != (arm == 'parent'):
                continue
            path = k.constructor_path(data, case, rep, arm, budget)
            value = k.read(path)
            k.require((value['case'], value['repetition'], value['arm'], value['cpu_budget']) ==
                      (case, rep, arm, budget), 'Freeze exact constructor identity')
            entries.append(dict(case=case, repetition=rep, arm=arm, cpu_budget=budget,
                path=str(path.relative_to(data)), sha256=k.sha(path)))
    k.write_new(data/('PARENTS_FROZEN.json' if parents_only else 'CONSTRUCTORS_FROZEN.json'), dict(h, entries=entries))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--worker', action='store_true')
    group.add_argument('--finalize', action='store_true')
    group.add_argument('--freeze', choices=('parents', 'others'))
    parser.add_argument('--case'); parser.add_argument('--repetition', type=int)
    parser.add_argument('--arm', choices=('parent',)+k.ARMS+k.CHEAP)
    parser.add_argument('--budget', type=int, default=0)
    parser.add_argument('--cpu-cap', type=float); parser.add_argument('--wall-cap', type=int)
    parser.add_argument('--raw', type=Path); parser.add_argument('--plan', type=Path)
    parser.add_argument('--receipt', type=Path); parser.add_argument('--julia-executable')
    args = parser.parse_args()
    if args.freeze:
        freeze(args.data, args.freeze == 'parents')
    elif args.worker:
        worker(args.data, args.case, args.repetition, args.arm, args.budget,
            args.cpu_cap, args.wall_cap, args.raw, args.julia_executable)
    else:
        finalize(args.data, args.case, args.repetition, args.arm, args.budget, args.raw, args.plan, args.receipt)
