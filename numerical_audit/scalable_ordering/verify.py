"""Re-audit every frozen candidate, actual primitive trace and paired outcome."""
import argparse
from collections import Counter
import itertools
import math
from pathlib import Path
import time
import numpy as np
import common as k
import compiler
from reference import Reference, metrics


def validate_derived_data(frozen, rows, constructors, checks, references, saved_pairs, summary):
    """Check saved report arithmetic against authenticated, verified records.

    This pure derived-data gate does no propagation or file I/O. Its caller
    authenticates constructors and verifies each row's state/resources first.
    Saved error values must exactly match their producer rows; independent
    state errors must also give the same winner/tie/censor decision.
    """
    def nonnegative(value, label):
        k.require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
                  'Independent finite derived value: ' + label)

    expected_rows = {(case, T, cap, arm): f'{case}.{arm}.t{ti}.c{cap}'
        for case, (ti, T), cap, arm in itertools.product(
            frozen['models'], enumerate(k.TIMES), k.CAPS, k.ARMS)}
    by_row = {(r['case'], r['T'], r['cap_multiple'], r['arm']): r for r in rows}
    by_check = {r['candidate']: r for r in checks}
    expected_constructors = set(itertools.product(frozen['models'], k.ARMS))
    expected_references = set(itertools.product(frozen['models'], k.TIMES))
    k.require(len(rows) == len(by_row) == len(checks) == len(by_check) == len(expected_rows) == 96 and
              set(by_row) == set(expected_rows) and set(by_check) == set(expected_rows.values()),
              'Independent derived row inventory')
    k.require(set(constructors) == expected_constructors and len(constructors) == 24 and
              len(references) == len(expected_references) == 8 and
              {(r['case'], r['T']) for r in references} == expected_references,
              'Independent derived constructor/reference inventory')
    for key, row in by_row.items():
        case, T, cap, arm = key
        check, constructor = by_check[expected_rows[key]], constructors[case, arm]
        k.require(row['candidate'] == expected_rows[key] and row['status'] in ('complete', 'censored') and
                  check['status'] == row['status'], 'Independent derived row identity/status')
        nonnegative(constructor['CPU_seconds'], 'constructor CPU')
        nonnegative(row['CPU_seconds'], 'compilation CPU')
        nonnegative(row['ordering_CPU_seconds'], 'row ordering CPU')
        k.require(row['ordering_CPU_seconds'] == constructor['CPU_seconds'] and
                  row['ordering_status'] == constructor['status'] and
                  row['ordering_evaluations'] == constructor['evaluations'],
                  'Independent row ordering cost/status')
        if row['status'] == 'complete':
            nonnegative(row['error'], 'producer state error')
            nonnegative(check['error'], 'independent state error')
            k.require(row['compiled_resources'] == check['resources'],
                      'Independent verified resource fields')
            k.require(abs(row['error'] - check['error']) < 2e-9, 'Independent verified error field')
        else:
            k.require(row['compiled_resources'] is None, 'Independent censored resource fields')

    expected_pairs = {(case, T, cap, left, right) for case, T, cap in
        itertools.product(frozen['models'], k.TIMES, k.CAPS)
        for left, right in itertools.combinations(k.ARMS, 2)}
    by_pair = {(r['case'], r['T'], r['cap_multiple'], r['left'], r['right']): r for r in saved_pairs}
    k.require(len(saved_pairs) == len(by_pair) == len(expected_pairs) == 240 and
              set(by_pair) == expected_pairs, 'Independent complete pair inventory')
    metric_fields = ('error_left', 'error_right', 'error_difference', 'CX_left', 'CX_right',
                     'depth_left', 'depth_right', 'construction_CPU_left', 'construction_CPU_right')
    independent_pairs = []
    groups = {control: Counter() for control in k.ARMS[:-1]}
    def winner(left, right, error_left, error_right):
        delta = error_left - error_right
        return left if delta < -1e-12 else right if delta > 1e-12 else 'tie'

    for case, T, cap in itertools.product(frozen['models'], k.TIMES, k.CAPS):
        for left, right in itertools.combinations(k.ARMS, 2):
            a, b = by_row[case, T, cap, left], by_row[case, T, cap, right]
            ca, cb = by_check[a['candidate']], by_check[b['candidate']]
            pair = by_pair[case, T, cap, left, right]
            complete = a['status'] == b['status'] == 'complete'
            k.require(pair['status'] == ('complete' if complete else 'censored'),
                      'Independent saved pair status')
            if complete:
                expected = dict(error_left=a['error'], error_right=b['error'],
                    error_difference=a['error'] - b['error'],
                    CX_left=ca['resources']['two_qubit_gates'], CX_right=cb['resources']['two_qubit_gates'],
                    depth_left=ca['resources']['unit_duration_ASAP_depth'],
                    depth_right=cb['resources']['unit_duration_ASAP_depth'],
                    construction_CPU_left=constructors[case, left]['CPU_seconds'] + a['CPU_seconds'],
                    construction_CPU_right=constructors[case, right]['CPU_seconds'] + b['CPU_seconds'])
                for field in metric_fields:
                    k.require(type(pair[field]) is type(expected[field]) and pair[field] == expected[field],
                              'Independent saved pair field: ' + field)
                decision = winner(left, right, ca['error'], cb['error'])
                k.require(pair['error_decision'] == winner(left, right, a['error'], b['error']) == decision,
                          'Independent pair decision')
            else:
                k.require(not any(field in pair for field in metric_fields + ('error_decision',)),
                          'Censored pair must omit completed metrics')
                decision = 'censored'
            independent_pairs.append(dict(case=case, T=T, cap_multiple=cap, left=left, right=right,
                                          independent_error_decision=decision))
            if right == 'cubic':
                groups[left][decision] += 1

    counts = dict(circuits=len(rows), references=len(references), all_pair_count=len(saved_pairs))
    for field, expected in counts.items():
        k.require(type(summary[field]) is int and summary[field] == expected,
                  'Independent summary count: ' + field)
    k.require(isinstance(summary['statuses'], dict) and
              all(type(count) is int and count >= 0 for count in summary['statuses'].values()) and
              summary['statuses'] == dict(Counter(r['status'] for r in rows)),
              'Independent summary statuses')
    k.require(isinstance(summary['cubic_error_comparisons'], dict) and
              all(isinstance(group, dict) and all(type(count) is int and count >= 0 for count in group.values())
                  for group in summary['cubic_error_comparisons'].values()) and
              summary['cubic_error_comparisons'] == {control: dict(group) for control, group in groups.items()},
              'Independent summary cubic comparisons')
    k.require(summary['metric'] == '1-abs(overlap), numerical error comparison with absolute ambiguity1e-12',
              'Independent summary metric/ambiguity')
    return dict(status='passed', **counts,
        checked_pair_fields=list(metric_fields) + ['status', 'error_decision'],
        checked_summary_fields=list(counts) + ['statuses', 'cubic_error_comparisons', 'metric']), independent_pairs


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
    checks, references, verified_rows, constructors = [], [], [], {}
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
            constructors[case, arm] = {field: record[field] for field in ('CPU_seconds', 'status', 'evaluations')}
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
                verified_rows.append(row)
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
    derived, pairs = validate_derived_data(frozen, verified_rows, constructors, checks, references,
                                         k.read(output / 'PAIRS.json'), k.read(output / 'SUMMARY.json'))
    result = dict(status='passed', candidates=checks, references=references, pairs=pairs,
        derived_data_validation=derived,
        derived_input_sha256={name: k.sha(output / name) for name in ('PAIRS.json', 'SUMMARY.json')},
        source_map=k.source_map(), environment=k.environment(), source_snapshot_executed=False,
        scope='Actual traces/DAGs/ideal states, separate Taylor references and saved pair/summary arithmetic; '
              'no Lean or interval certificate')
    k.write(destination / 'INDEPENDENT_VERIFIED.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.output, args.destination)
    print('INDEPENDENT_PASS', len(result['candidates']), len(result['pairs']), flush=True)
