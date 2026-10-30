"""Single-worker dispatcher with immutable plans and cumulative phase charges.

This is bookkeeping outside the scientific constructor. It never evaluates J,
chooses an order, or promotes a candidate after the owning service exits.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import subprocess
import sys
import time
import common as k
import construct

GUARD = k.HERE/'tools/resource_guard.py'


def initialize(data, limits):
    k.header(data, producer=True)
    k.require(not (data/'PHASE_LIMITS.json').exists(), 'Immutable phase pool declaration')
    k.require(isinstance(limits, dict) and limits, 'Named finite phase pools')
    for name, value in limits.items():
        k.require(name in ('inputs', 'constructors', 'panels', 'physical', 'independent', 'objective', 'analysis', 'package'),
                  'Known resource phase')
        k.require(set(value) == {'CPU_seconds', 'wall_seconds', 'disk_MiB'} and
                  all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in value.values()),
                  'Finite positive phase limits')
    k.write_new(data/'PHASE_LIMITS.json', dict(schema='molecular-phase-pools-v1', phases=limits,
        source_freeze_sha256=k.sha(data/'INPUT_PROTOCOL_FROZEN.json')))


def remaining(data, phase):
    declaration = k.read(data/'PHASE_LIMITS.json')
    k.require(declaration['source_freeze_sha256'] == k.sha(data/'INPUT_PROTOCOL_FROZEN.json'), 'Phase/source freeze')
    limits = declaration['phases'][phase]
    cpu, wall = 0., 0.
    for path in (data/'jobs').glob('*.PLAN.json'):
        plan = k.read(path)
        if plan['phase'] != phase:
            continue
        final = path.with_name(path.name.replace('.PLAN.json', '.FINAL.json'))
        k.require(final.is_file(), 'Unfinished owning job blocks next dispatch')
        record = k.read(final)
        k.require(record['plan_sha256'] == k.sha(path), 'Completed job/plan hash')
        # Unknown costs spend the entire declared attempt; never free retries.
        cpu += plan['cpu_cap'] if record['charged_CPU_seconds'] is None else record['charged_CPU_seconds']
        wall += plan['wall_cap'] if record['wall_seconds'] is None else record['wall_seconds']
    return dict(CPU_seconds=limits['CPU_seconds']-cpu, wall_seconds=limits['wall_seconds']-wall,
                disk_MiB=limits['disk_MiB'])


def dispatch(data, label, phase, command, cpu, wall, disk, *, identity=None, expected_raw=None):
    h, _ = k.header(data, producer=True)
    k.require(label.isascii() and label.replace('-', '').isalnum(), 'Bounded owning label')
    k.require(all(isinstance(a, str) and '\x00' not in a for a in command), 'Literal argv list')
    allowance = remaining(data, phase)
    k.require(0 < cpu <= allowance['CPU_seconds'] and 0 < wall <= allowance['wall_seconds']
              and 0 < disk <= allowance['disk_MiB'], 'Attempt exceeds remaining phase pool')
    jobs, logs = data/'jobs', data/'logs'
    jobs.mkdir(exist_ok=True); logs.mkdir(exist_ok=True)
    plan_path, final_path = jobs/(label+'.PLAN.json'), jobs/(label+'.FINAL.json')
    receipt_path = logs/(label+'.resources.json')
    k.require(not any(p.exists() for p in (plan_path, final_path, receipt_path, logs/(label+'.log'))),
              'Preserve attempted job identity and all costs')
    plan = dict(schema='molecular-job-plan-v1', label=label, phase=phase, command=command,
        cpu_cap=cpu, wall_cap=wall, disk_mib=disk, source_freeze_sha256=h['source_freeze_sha256'],
        guard_sha256=k.sha(GUARD), producer_root=str(k.ROOT), producer_data_root=str(data.resolve()),
        identity=identity, expected_raw=expected_raw)
    k.write_new(plan_path, plan)
    started = time.monotonic()
    invocation = [sys.executable, str(GUARD), '--label', label, '--logs', str(logs),
        '--cpu-seconds', repr(cpu), '--wall-seconds', str(wall), '--disk-mib', str(disk),
        '--output', str(data), '--']+command
    status = subprocess.run(invocation, cwd=k.ROOT, check=False).returncode
    elapsed = time.monotonic()-started
    receipt = k.read(receipt_path) if receipt_path.is_file() else None
    costs = dict(raw_CPU_seconds=None, charged_CPU_seconds=None, wall_seconds=None)
    if receipt is not None:
        k.require(receipt['command'] == command, 'Actual guarded command changed')
        costs = construct.cost_from_receipt(receipt)
    k.write_new(final_path, dict(schema='molecular-job-final-v1', label=label, phase=phase,
        plan_sha256=k.sha(plan_path), receipt_path=str(receipt_path.relative_to(data)),
        receipt_sha256=k.sha(receipt_path) if receipt is not None else None,
        guard_supervisor_exit_code=status, dispatch_wall_seconds=elapsed, **costs))
    k.require(receipt is not None, 'Missing owning service receipt; retained failed attempt')
    return plan_path, receipt_path


def constructors(data):
    h, config = k.header(data, producer=True)
    julia = k.frozen(data)['runtime']['julia']['executable']
    # No daemon, warm shared worker, or inherited free JIT process.
    def run_one(case, rep, arm, budget):
        target = k.constructor_path(data, case, rep, arm, budget)
        k.require(not target.exists(), 'Constructor already attempted/finalized')
        p, _ = construct.parent(data, case, rep) if arm in k.ARMS else (None, None)
        cpu, wall = construct.allowances(config, arm, budget, p)
        label = f'mg-{case.lower().replace("_", "-")}-r{rep}-{arm}-b{budget}'
        raw = data/'workers'/(label+'.json')
        raw.parent.mkdir(exist_ok=True)
        command = [sys.executable, str(k.HERE/'construct.py'), '--data', str(data), '--worker',
            '--case', case, '--repetition', str(rep), '--arm', arm, '--budget', str(budget),
            '--cpu-cap', repr(cpu), '--wall-cap', str(wall), '--raw', str(raw), '--julia-executable', julia]
        plan, receipt = dispatch(data, label, 'constructors', command, cpu, wall, 512,
            identity=dict(case=case, repetition=rep, arm=arm, cpu_budget=budget),
            expected_raw=str(raw.relative_to(data)))
        value = construct.finalize(data, case, rep, arm, budget, raw, plan, receipt)
        print('FINALIZED', case, rep, arm, budget, 'eligible', value['eligible'], flush=True)
    for case in h['cases']:
        for rep in range(1, h['repetitions']+1):
            run_one(case, rep, 'parent', 0)
    construct.freeze(data, True)
    for case in h['cases']:
        for rep, arm, budget in k.identities(h['repetitions']):
            if arm != 'parent':
                run_one(case, rep, arm, budget)
    construct.freeze(data, False)


def validation(data, phase, cpu, wall):
    h, _ = k.header(data, producer=True)
    julia = k.frozen(data)['runtime']['julia']['executable']
    base = [julia, '--startup-file=no', '--history-file=no', '--threads=1', '--gcthreads=1',
            '--heap-size-hint=256M', '--project='+str(k.ROOT/'analysis/julia_window')]
    for case in h['cases']:
        model, panel = data/'inputs'/case/'model.toml', data/'panels'/(case+'.toml')
        folder = data/phase/case
        folder.parent.mkdir(parents=True, exist_ok=True)
        command = base+[str(k.HERE/('verify.jl' if phase == 'independent' else phase+'.jl')),
                        str(model), str(panel)]
        if phase == 'independent':
            command += [str(data/'physical'/case), str(folder)]
        else:
            command.append(str(folder))
            if phase == 'objective':
                command.append(str(data/'witnesses'/(case+'.toml')))
        dispatch(data, 'mg-'+phase+'-'+case.lower().replace('_', '-'), phase, command, cpu, wall, 512,
                 identity=dict(case=case))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--initialize', type=Path, help='JSON map of finite phase limits')
    group.add_argument('--constructors', action='store_true')
    group.add_argument('--validate', choices=('physical', 'independent', 'objective'))
    parser.add_argument('--cpu-cap', type=float); parser.add_argument('--wall-cap', type=int)
    args = parser.parse_args(); args.data = args.data.resolve()
    if args.initialize:
        initialize(args.data, k.read(args.initialize))
    elif args.constructors:
        constructors(args.data)
    else:
        validation(args.data, args.validate, args.cpu_cap, args.wall_cap)
