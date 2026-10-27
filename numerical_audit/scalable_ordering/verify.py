"""Re-audit every frozen candidate, actual primitive trace and paired outcome."""
import argparse
import itertools
from pathlib import Path
import time
import numpy as np
import common as k
import compiler
from reference import Reference, metrics


def verify(output, destination):
    k.require(not destination.exists(), 'Fresh independent output required')
    frozen = k.read(output / 'ALL_ORDERS_FROZEN.json')
    k.require(frozen['sources'] == k.source_map() and frozen['environment'] == k.environment(),
              'Current verification sources/environment differ from frozen producer')
    k.require(frozen['arms'] == list(k.ARMS) and frozen['times'] == list(k.TIMES) and
              frozen['caps'] == list(k.CAPS) and len(frozen['models']) == 4, 'Complete frozen panel')
    expected = {f'{case}.{arm}.t{ti}.c{cap}' for case, arm, ti, cap in
                itertools.product(frozen['models'], k.ARMS, range(len(k.TIMES)), k.CAPS)}
    actual = {p.stem for p in (output / 'rows').glob('*.json')}
    k.require(actual == expected and len(actual) == 96, 'Exact candidate inventory')
    destination.mkdir(parents=True)
    checks, references = [], []
    for case, digest in frozen['models'].items():
        path = output / 'chemistry/inputs' / case / 'model.json'
        k.require(k.sha(path) == digest, 'Independent model authentication')
        model = k.read(path)
        reference = Reference(model)
        for arm in k.ARMS:
            record_path = output / f'constructors/{case}.{arm}.json'
            k.require(k.sha(record_path) == frozen['constructors'][str(record_path.relative_to(output))],
                      'Independent objective constructor identity')
            record = k.read(record_path)
            objective = reference.objective(record, k.TRAIN_TIME, 30,
                                            kind='first' if arm in ('first', 'pauli') else 'cubic')
            original_objective = k.read(output / f'objectives/{case}.{arm}.json')
            k.require(objective['kind'] == original_objective['kind'] and
                      objective['steps'] == original_objective['steps'] and
                      abs(objective['value'] - original_objective['value']) <=
                      1e-20 + 1e-7 * max(objective['value'], original_objective['value']),
                      'Independent saved objective')
            if record['score'] is not None:
                k.require(abs(objective['value'] - record['score']) <=
                          1e-20 + 1e-7 * max(objective['value'], record['score']),
                          'Independent frozen trained objective')
        for ti, T in enumerate(k.TIMES):
            cpu = time.process_time()
            exact = reference.series_action(T)
            saved_path = output / f'references/{case}.t{ti}.npz'
            saved = k.load_vector(saved_path, 'exact', model['n'])
            difference = float(np.max(np.abs(saved - exact)))
            k.require(difference < 2e-9, 'Independent Taylor and exponential references differ')
            references.append(dict(case=case, T=T, max_difference=difference,
                                   CPU_seconds=time.process_time() - cpu))
            for arm, cap in itertools.product(k.ARMS, k.CAPS):
                label = f'{case}.{arm}.t{ti}.c{cap}'
                row = k.read(output / f'rows/{label}.json')
                record_path = output / f'constructors/{case}.{arm}.json'
                k.require(k.sha(record_path) == frozen['constructors'][str(record_path.relative_to(output))],
                          'Independent constructor authentication')
                record = k.read(record_path)
                k.require((row['candidate'], row['case'], row['arm'], row['T'], row['cap_multiple']) ==
                          (label, case, arm, T, cap), 'Row logical identity')
                trace = output / f'traces/{label}.ndjson.gz'
                k.require(row['trace_file'] == str(trace.relative_to(output)) and
                          row['trace_sha256'] == k.sha(trace), 'Actual trace hash/path')
                lengths = record['unit_lengths']
                steps = 1 if len(lengths) == 1 else (cap * sum(lengths) - lengths[0]) // (
                    2 * sum(lengths) - lengths[0] - lengths[-1])
                k.require(row['steps'] == steps and row['identity'] == dict(case=case, arm=arm,
                    model_sha256=digest, constructor_sha256=frozen['constructors'][str(record_path.relative_to(output))],
                    T=T, steps=steps, cap_multiple=cap, formula=2, preparation='HF_prefix',
                    compilation_status=row['status']), 'Independent circuit identity/step choice')
                if row['status'] != 'complete':
                    k.require(row['status'] == 'censored' and row['compiled_resources'] is None,
                              'Retained censor metadata')
                    checks.append(dict(candidate=label, status='censored'))
                    continue
                state, counts = reference.audit_trace(trace, row['identity'], compiler.VERSION)
                counts.pop('actual_Rz_rotations')
                k.require(counts == row['compiled_resources'], 'Independent actual circuit resources')
                ideal = reference.ideal_formula(record, steps, T)
                overlap = np.vdot(ideal, state)
                k.require(abs(overlap) > 0 and float(np.max(np.abs(state - ideal * overlap / abs(overlap)))) < 2e-9,
                          'Independent compiled/ideal state parity')
                values = metrics(exact, state, model)
                for key in ('error', 'squared_infidelity', 'norm', 'number_leakage', 'joint_sector_leakage'):
                    k.require(abs(values[key] - row[key]) < 2e-9, 'Independent physical metric: ' + key)
                saved_state = output / f'states/{label}.npz'
                k.require(row['state_file'] == str(saved_state.relative_to(output)) and
                          row['state_sha256'] == k.sha(saved_state), 'Independent state identity')
                original = k.load_vector(saved_state, 'actual', model['n'])
                k.require(original.shape == state.shape and float(np.max(np.abs(original - state))) < 2e-11,
                          'Independent saved state')
                checks.append(dict(candidate=label, status='complete', error=values['error'], resources=counts))
        del reference
    k.require(len(checks) == 96 and len(references) == 8, 'Complete independent cohort')
    pairs = []
    for case, ti, cap in itertools.product(frozen['models'], range(len(k.TIMES)), k.CAPS):
        panel = {arm: next(r for r in checks if r['candidate'] == f'{case}.{arm}.t{ti}.c{cap}') for arm in k.ARMS}
        for left, right in itertools.combinations(k.ARMS, 2):
            a, b = panel[left], panel[right]
            if a['status'] == b['status'] == 'complete':
                delta = a['error'] - b['error']
                decision = left if delta < -1e-12 else right if delta > 1e-12 else 'tie'
            else:
                decision = 'censored'
            pairs.append(dict(case=case, T=k.TIMES[ti], cap_multiple=cap, left=left, right=right,
                              independent_error_decision=decision))
    original = k.read(output / 'PAIRS.json')
    k.require(len(original) == len(pairs) == 240, 'Independent complete pair count')
    by_key = {(r['case'], r['T'], r['cap_multiple'], r['left'], r['right']): r for r in original}
    k.require(len(by_key) == 240, 'No duplicate producer pairs')
    for r in pairs:
        other = by_key[r['case'], r['T'], r['cap_multiple'], r['left'], r['right']]
        k.require(other.get('error_decision', 'censored') == r['independent_error_decision'],
                  'Independent pair decision')
    result = dict(status='passed', candidates=checks, references=references, pairs=pairs,
        source_map=k.source_map(), environment=k.environment(), source_snapshot_executed=False,
        scope='Actual traces/DAGs/ideal states and separate Taylor references; no Lean or interval certificate')
    k.write(destination / 'INDEPENDENT_VERIFIED.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.output, args.destination)
    print('INDEPENDENT_PASS', len(result['candidates']), len(result['pairs']), flush=True)
