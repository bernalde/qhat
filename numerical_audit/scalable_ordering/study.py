"""Freeze all orders, then compile and independently validate the complete panel."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from collections import Counter
import itertools
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
import common as k
import training
import compilation
import compiler
from reference import Reference, metrics


def authenticate_inputs(inputs):
    frozen = k.read(inputs / 'INPUTS_FROZEN.json')
    k.require(not frozen['refusals'] and frozen['expected_cases'] == 4 and
              frozen['protocol_sha256'] == k.sha(inputs / 'INPUT_PROTOCOL_FROZEN.json'), 'Incomplete/changed chemistry')
    protocol = k.read(inputs / 'INPUT_PROTOCOL_FROZEN.json')
    k.require(protocol['registry'] == k.read(k.HERE / 'REGISTRY.json'), 'Registry differs from frozen input panel')
    for name, digest in protocol['sources'].items():
        k.require(k.sha(inputs / 'source_snapshot' / name) == digest, 'Chemistry source snapshot changed')
    for case, digest in frozen['successes'].items():
        path = inputs / 'inputs' / case
        model = k.read(path / 'model.json')
        k.require(k.sha(path / 'model.json') == digest and
                  k.sha(path / 'tensors.npz') == model['tensors_sha256'] and
                  k.sha(path / 'orbitals.npz') == model['orbitals_sha256'], 'Changed chemistry bytes')
    return frozen


def prepare(inputs, output):
    k.require(not output.exists(), 'Fresh study directory required; preserve old attempts')
    head = k.committed_sources()
    frozen = authenticate_inputs(inputs)
    sources = k.source_map()
    output.mkdir(parents=True)
    shutil.copytree(inputs, output / 'chemistry')
    for name in sources:
        target = output / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT / name, target)
    (output / 'constructors').mkdir()
    constructors = {}
    started, cpu = time.monotonic(), time.process_time()
    for index, (case, digest) in enumerate(frozen['successes'].items()):
        path = output / 'chemistry' / 'inputs' / case / 'model.json'
        # Rotate deterministic dispatch to distribute process/cache order effects.
        arms = k.ARMS[index % len(k.ARMS):] + k.ARMS[:index % len(k.ARMS)]
        for arm in arms:
            record = training.construct(path, digest, arm)
            relative = f'constructors/{case}.{arm}.json'
            k.write(output / relative, record)
            constructors[relative] = k.sha(output / relative)
            print('FROZEN_ORDER', case, arm, record['status'], record['evaluations'], flush=True)
    k.require(len(constructors) == 24, 'Freeze every arm before any validation error')
    k.write(output / 'ALL_ORDERS_FROZEN.json', dict(source_commit=head, sources=sources,
        environment=k.environment(), input_manifest_sha256=k.sha(output / 'chemistry/INPUTS_FROZEN.json'),
        models=frozen['successes'], constructors=constructors, arms=k.ARMS,
        train_time=k.TRAIN_TIME, times=k.TIMES, caps=k.CAPS, train_CPU_seconds=k.TRAIN_CPU_SECONDS,
        expected_circuits=96, expected_references=8,
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - started,
        note='All orders frozen; no validation reference, compiled state or error has been evaluated'))


def bind(output):
    frozen = k.read(output / 'ALL_ORDERS_FROZEN.json')
    k.require(frozen['sources'] == k.source_map() and frozen['environment'] == k.environment(),
              'Execution source or environment changed after order freeze')
    k.require(frozen['times'] == list(k.TIMES) and frozen['caps'] == list(k.CAPS) and
              frozen['arms'] == list(k.ARMS) and len(frozen['constructors']) == 24, 'Frozen panel differs')
    authenticate_inputs(output / 'chemistry')
    for name, digest in frozen['constructors'].items():
        k.require(k.sha(output / name) == digest, 'Frozen constructor changed')
    return frozen


def run(output):
    frozen = bind(output)
    for name in ('traces', 'rows', 'states', 'references', 'objectives'):
        k.require(not (output / name).exists(), 'Fresh run required; preserve every partial attempt')
        (output / name).mkdir()
    references, rows = {}, []
    for case, digest in frozen['models'].items():
        model = k.read(output / 'chemistry/inputs' / case / 'model.json')
        context = Reference(model)
        for arm in k.ARMS:
            record = k.read(output / f'constructors/{case}.{arm}.json')
            objective = context.objective(record, k.TRAIN_TIME, 30,
                                          kind='first' if arm in ('first', 'pauli') else 'cubic')
            if record['score'] is not None:
                k.require(math.isclose(record['score'], objective['value'], rel_tol=1e-7, abs_tol=1e-20),
                          'Independent final frozen training objective differs')
            k.write(output / f'objectives/{case}.{arm}.json', objective)
        for ti, T in enumerate(k.TIMES):
            cpu, wall = time.process_time(), time.monotonic()
            exact = context.exact_action(T)
            refname = f'references/{case}.t{ti}.npz'
            np.savez_compressed(output / refname, exact=exact)
            references[refname] = dict(sha256=k.sha(output / refname), case=case, T=T,
                CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall)
            for arm in k.ARMS:
                record = k.read(output / f'constructors/{case}.{arm}.json')
                for cap in k.CAPS:
                    label = f'{case}.{arm}.t{ti}.c{cap}'
                    trace = output / f'traces/{label}.ndjson.gz'
                    identity = dict(case=case, arm=arm, model_sha256=digest,
                                    constructor_sha256=frozen['constructors'][f'constructors/{case}.{arm}.json'])
                    result = compilation.compile_candidate(model, record, T, cap, trace, identity)
                    row = dict(candidate=label, case=case, arm=arm, T=T, cap_multiple=cap,
                        ordering_CPU_seconds=record['CPU_seconds'], ordering_status=record['status'],
                        ordering_evaluations=record['evaluations'], trace_file=str(trace.relative_to(output)),
                        reference_file=refname, **result)
                    if result['status'] == 'complete':
                        cpu = time.process_time()
                        actual, costs = context.audit_trace(trace, result['identity'], compiler.VERSION)
                        expected_costs = dict(costs)
                        expected_costs.pop('actual_Rz_rotations')
                        k.require(expected_costs == result['compiled_resources'], 'Independent DAG costs differ')
                        ideal = context.ideal_formula(record, result['steps'], T)
                        overlap = np.vdot(ideal, actual)
                        k.require(abs(overlap) > 0 and
                                  float(np.max(np.abs(actual - ideal * overlap / abs(overlap)))) < 2e-9,
                                  'Actual primitive circuit differs from independent ideal formula')
                        row.update(metrics(exact, actual, model))
                        k.require(abs(row['norm'] - 1) < 1e-8, 'Compiled state norm drift')
                        if record['partition'] == 'conserving':
                            k.require(row['joint_sector_leakage'] < 1e-16, 'Conserving-unit sector leakage')
                        row['validation_CPU_seconds'] = time.process_time() - cpu
                        state = output / f'states/{label}.npz'
                        np.savez_compressed(state, actual=actual)
                        row.update(state_file=str(state.relative_to(output)), state_sha256=k.sha(state))
                    k.write(output / f'rows/{label}.json', row)
                    rows.append(label)
                    print('VALIDATED', label, result['status'], flush=True)
        del context
    k.require(len(rows) == 96 and len(references) == 8, 'Complete declared validation panel')
    k.write(output / 'RUN_COMPLETE.json', dict(rows=rows, references=references,
        order_freeze_sha256=k.sha(output / 'ALL_ORDERS_FROZEN.json'), source_map=k.source_map()))


def summarize(output):
    bind(output)
    complete = k.read(output / 'RUN_COMPLETE.json')
    rows = [k.read(output / f'rows/{name}.json') for name in complete['rows']]
    k.require(len(rows) == 96 and len({r['candidate'] for r in rows}) == 96, 'Row count/identity')
    pairs = []
    for case, T, cap in itertools.product(k.read(output / 'ALL_ORDERS_FROZEN.json')['models'], k.TIMES, k.CAPS):
        panel = {r['arm']: r for r in rows if (r['case'], r['T'], r['cap_multiple']) == (case, T, cap)}
        k.require(set(panel) == set(k.ARMS), 'Complete within-setting arms')
        for left, right in itertools.combinations(k.ARMS, 2):
            a, b = panel[left], panel[right]
            record = dict(case=case, T=T, cap_multiple=cap, left=left, right=right,
                          status='complete' if a['status'] == b['status'] == 'complete' else 'censored')
            if record['status'] == 'complete':
                delta = a['error'] - b['error']
                record.update(error_left=a['error'], error_right=b['error'], error_difference=delta,
                    error_decision=left if delta < -1e-12 else right if delta > 1e-12 else 'tie',
                    CX_left=a['compiled_resources']['two_qubit_gates'], CX_right=b['compiled_resources']['two_qubit_gates'],
                    depth_left=a['compiled_resources']['unit_duration_ASAP_depth'],
                    depth_right=b['compiled_resources']['unit_duration_ASAP_depth'],
                    construction_CPU_left=a['ordering_CPU_seconds'] + a['CPU_seconds'],
                    construction_CPU_right=b['ordering_CPU_seconds'] + b['CPU_seconds'])
            pairs.append(record)
    k.require(len(pairs) == 240, 'Retain all paired comparisons, not selected winners')
    groups = {}
    for control in k.ARMS[:-1]:
        group = [r for r in pairs if {r['left'], r['right']} == {'cubic', control}]
        groups[control] = dict(Counter(r.get('error_decision', 'censored') for r in group))
    k.write(output / 'PAIRS.json', pairs)
    summary = dict(circuits=len(rows), references=len(complete['references']),
        statuses=dict(Counter(r['status'] for r in rows)), cubic_error_comparisons=groups,
        all_pair_count=len(pairs), metric='1-abs(overlap), numerical error comparison with absolute ambiguity1e-12',
        limits=['Four12-modecases/threefamilies;dependenttimes/caps/geometries.',
            'IdealprototypecompilerandHFpreparation;nohardwareorSCF/basiserrorclaim.',
            'MatchedsearchCPUallowancesarenotmatchedconsumedCPUorconstructioncharges.',
            'Exactreferencevalidationisoffline;responseheuristicdoesnotcertifyfiniteerror.',
            'Full-spacevalidationisbounded;thetraineraloneavoids2**nallocation.'])
    k.write(output / 'SUMMARY.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('prepare', 'run', 'summarize'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inputs', type=Path)
    args = parser.parse_args()
    if args.phase == 'prepare':
        k.require(args.inputs is not None, 'Preparation needs authenticated chemistry')
        prepare(args.inputs, args.output)
    elif args.phase == 'run':
        run(args.output)
    else:
        print(json.dumps(summarize(args.output), indent=2))
