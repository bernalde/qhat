"""Resource-only representative test on a 7/5-scaled twelve-mode input."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
import numpy as np
import common as k
import compiler
import compilation
from reference import Reference, metrics
import training
from study import authenticate_inputs


def scaled_model(model):
    altered = json.loads(json.dumps(model))
    altered['case'] = 'preflight_' + model['case'] + '_scaled7over5'
    altered['denominator'] = str(5 * int(model['denominator']))
    altered['coefficients'] = [[word, str(7 * int(c))] for word, c in altered['coefficients']]
    altered['blocks'] = [[[word, str(7 * int(c))] for word, c in block] for block in altered['blocks']]
    for parent in altered['fermionic_parents']:
        parent['components'] = [[word, str(7 * int(c))] for word, c in parent['components']]
        parent['removed_scalar'] = str(7 * int(parent['removed_scalar']))
    altered['removed_scalar_exact'] = str(Fraction(model['removed_scalar_exact']) * Fraction(7, 5))
    altered['preflight'] = dict(scale='7/5', original_case=model['case'],
        original_HF_fields_are_provenance_only=True,
        scope='Altered supplied Hamiltonian for mechanics/resources; no main-case validation claim')
    return altered


def run(inputs, output):
    k.require(not output.exists(), 'Fresh preflight output required')
    head = k.committed_sources()
    frozen = authenticate_inputs(inputs)
    candidates = []
    for case, digest in frozen['successes'].items():
        model = k.read(inputs / 'inputs' / case / 'model.json')
        candidates.append((-len(model['coefficients']), case, model, digest))
    _, original_case, original, original_digest = min(candidates, key=lambda x: x[:2])
    model = scaled_model(original)
    output.mkdir(parents=True)
    k.write(output / 'model.json', model)
    digest = k.sha(output / 'model.json')
    orders = {}
    cpu, wall = time.process_time(), time.monotonic()
    for arm in k.ARMS:
        record = training.construct(output / 'model.json', digest, arm)
        k.write(output / f'order-{arm}.json', record)
        orders[arm] = k.sha(output / f'order-{arm}.json')
    k.write(output / 'PREFLIGHT_ORDERS_FROZEN.json', dict(source_commit=head, sources=k.source_map(),
        environment=k.environment(), selection='Largest Pauli-word count; lexical case tie; no error used',
        original_case=original_case, original_model_sha256=original_digest, model_sha256=digest,
        scale='7/5', orders=orders, T=1.5, cap_multiple=30))
    reference = Reference(model)
    for arm in k.ARMS:
        record = k.read(output / f'order-{arm}.json')
        objective = reference.objective(record, k.TRAIN_TIME, 30,
                                        kind='first' if arm in ('first', 'pauli') else 'cubic')
        if record['score'] is not None:
            k.require(abs(record['score'] - objective['value']) <=
                      1e-20 + 1e-7 * max(record['score'], objective['value']),
                      'Preflight independent frozen training objective')
        k.write(output / f'objective-{arm}.json', objective)
    exact = reference.exact_action(1.5)
    independent = reference.series_action(1.5)
    k.require(float(np.max(abs(exact - independent))) < 2e-9, 'Preflight reference parity')
    np.savez_compressed(output / 'reference.npz', exact=exact)
    rows = []
    for arm in k.ARMS:
        record = k.read(output / f'order-{arm}.json')
        trace = output / f'{arm}.ndjson.gz'
        result = compilation.compile_candidate(model, record, 1.5, 30, trace,
                                               dict(case=model['case'], arm=arm, mode='preflight'))
        row = dict(arm=arm, ordering_CPU_seconds=record['CPU_seconds'], ordering_status=record['status'],
                   ordering_evaluations=record['evaluations'], **result)
        if result['status'] == 'complete':
            actual, counts = reference.audit_trace(trace, result['identity'], compiler.VERSION)
            counts.pop('actual_Rz_rotations')
            k.require(counts == result['compiled_resources'], 'Preflight independent DAG parity')
            ideal = reference.ideal_formula(record, result['steps'], 1.5)
            overlap = np.vdot(ideal, actual)
            k.require(abs(overlap) > 0 and float(np.max(abs(actual - ideal * overlap / abs(overlap)))) < 2e-9,
                      'Preflight actual primitive/ideal parity')
            row.update(metrics(exact, actual, model))
        rows.append(row)
        k.write(output / f'result-{arm}.json', row)
    complete = all(r['status'] == 'complete' for r in rows)
    result = dict(status='passed' if complete else 'censored', mode='altered-input preflight',
        original_case=original_case, model_sha256=digest, arms=rows,
        max_trace_bytes=max(r['trace_bytes'] for r in rows),
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
        main_cases_evaluated=False, main_errors_exposed=False,
        note='Resource feasibility only; do not tune main protocol using these accuracy values')
    k.write(output / 'PREFLIGHT_COMPLETE.json', result)
    k.require(complete, 'Preflight contains retained incomplete circuits')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.inputs, args.output)
    print('PREFLIGHT_PASS', result['original_case'], len(result['arms']), flush=True)
