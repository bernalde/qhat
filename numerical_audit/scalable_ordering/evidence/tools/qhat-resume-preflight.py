#!/usr/bin/env python3
"""Continue an interrupted preflight using exactly its frozen model and orders.

Run only after the original owned worker has exited, under the task's external
resource guard and a prospectively allocated continuation allowance. No
constructor or objective search is rerun. Existing completed circuit rows are
authenticated and reused, and orphan traces are retained before replacement.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import uuid

sys.dont_write_bytecode = True


PROJECT_ROOT = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')
ARMS = ('canonical', 'magnitude', 'current', 'first', 'pauli', 'cubic')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, label, *, positive=False):
    require(type(value) in (float, int) and math.isfinite(value) and
            (value > 0 if positive else value >= 0), 'Invalid finite number: ' + label)
    return value


def bounded_integer(value, lower, upper, label):
    require(type(value) is int and lower <= value <= upper, 'Invalid bounded integer: ' + label)
    return value


def artifact(path, root=None):
    """Streaming authentication; no entire trace is allocated."""
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'Expected owned regular artifact: ' + str(path))
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return dict(path=str(path.relative_to(root)) if root is not None else str(path),
                bytes=path.stat().st_size, sha256=value.hexdigest())


class Journal:
    def __init__(self, path, output, started_cpu, started_wall):
        self.path, self.output = path, output
        self.cpu, self.wall = started_cpu, started_wall
        require(not path.exists() and not path.with_suffix(path.suffix + '.partial').exists(),
                'Fresh continuation record and partial path required')
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x', encoding='utf-8') as stream:
            stream.write('{}\n')
        self.data = dict(schema='qhat-preflight-continuation-v1', status='initializing',
            started_UTC=datetime.now(timezone.utc).isoformat(), output=str(output),
            renderer=artifact(Path(__file__).resolve()), reused_artifacts={}, retained_interrupted=[],
            arm_events=[], reference_check=None, same_frozen_orders=True,
            new_training=False, main_cases_evaluated=False, main_errors_exposed=False,
            old_attempt_costs_excluded=True,
            cost_scope='Own continuation process CPU/wall only, starting after stdlib imports. '
                'Original-attempt charges remain in the external task ledger and reused row metadata; '
                'they are not added to this continuation stage. The external guard is authoritative '
                'for subprocess CPU and total phase allowances.',
            resource_policy='Requires existing one-worker, quarter-core, one-thread, 512 MiB/no-swap '
                'guard and headroom rules; this script grants no additional allowance.')
        self.save()

    def save(self):
        self.data['own_stage_CPU_seconds'] = time.process_time() - self.cpu
        self.data['own_stage_wall_seconds'] = time.monotonic() - self.wall
        temporary = self.path.with_suffix(self.path.suffix + '.partial')
        with temporary.open('w', encoding='utf-8') as stream:
            json.dump(self.data, stream, indent=2, allow_nan=False)
            stream.write('\n')
        temporary.replace(self.path)

    def bind(self, path, *, label=None):
        value = artifact(path, self.output)
        key = label if label is not None else value['path']
        previous = self.data['reused_artifacts'].get(key)
        require(previous is None or previous == value, 'Reused artifact changed: ' + key)
        self.data['reused_artifacts'][key] = value
        return value


def validate_resources(resources, compiler):
    require(isinstance(resources, dict), 'Missing retained circuit resources')
    for key in ('gates', 'single_qubit_gates', 'two_qubit_gates',
                'unit_duration_ASAP_depth', 'CX_dependency_depth'):
        bounded_integer(resources[key], 0, compiler.MAX_STORED_GATES, 'resources/' + key)
    require(resources['gates'] == resources['single_qubit_gates'] + resources['two_qubit_gates'] and
            resources['CX_dependency_depth'] <= resources['unit_duration_ASAP_depth'] <= resources['gates'],
            'Retained resource counts/depth differ')


def authenticate_row(row, arm, model, order, trace, journal, k, compiler, compilation, *, new=False):
    require(row['arm'] == arm and row['status'] in ('complete', 'censored') and
            row['compiler'] == compiler.VERSION, 'Retained row arm/status/compiler differs')
    require(row['ordering_CPU_seconds'] == order['CPU_seconds'] and
            row['ordering_status'] == order['status'] and row['ordering_evaluations'] == order['evaluations'],
            'Retained row ordering charge/status differs')
    finite(row['CPU_seconds'], 'retained compilation CPU')
    finite(row['wall_seconds'], 'retained compilation wall')
    lengths = order['unit_lengths']
    steps = compilation.affordable_steps(lengths, 30)
    expected = compilation.rotation_count(lengths, steps)
    identity = dict(case=model['case'], arm=arm, mode='preflight', T=1.5, steps=steps,
                    cap_multiple=30, formula=2, preparation='HF_prefix', compilation_status=row['status'])
    require(row['identity'] == identity and row['steps'] == steps and row['cap_multiple'] == 30 and
            row['rotation_cap'] == 30 * sum(lengths) and row['expected_ideal_rotations'] == expected and
            row['quantum_preparation_circuit_cost'] == 'HF prefix included',
            'Retained circuit logical identity/rotation cap differs')
    bounded_integer(row['generated'], 0, compiler.MAX_GENERATED_GATES, 'generated primitive gates')
    bounded_integer(row['ideal_rotations'], 0, expected, 'retained ideal rotations')
    # The validator's existing line and gate bounds imply this conservative
    # compressed-file bound, without introducing a smaller circuit gate cap.
    maximum_trace_bytes = (compiler.MAX_STORED_GATES + 1) * 16384 + 1024**2
    bounded_integer(row['trace_bytes'], 1, maximum_trace_bytes, 'retained compressed trace bytes')
    require(trace.is_file() and not trace.is_symlink() and trace.stat().st_size == row['trace_bytes'],
            'Retained trace is missing or has changed size')
    saved_trace = journal.bind(trace, label='new/' + trace.name if new else None)
    require(saved_trace['sha256'] == row['trace_sha256'], 'Retained trace SHA-256 differs')
    with gzip.open(trace, 'rt') as stream:
        header = stream.readline(16385)
        require(0 < len(header) <= 16384, 'Retained trace header bound')
        header = json.loads(header)
    require(header['compiler'] == compiler.VERSION and header['n'] == model['n'] and
            header['identity'] == identity, 'Retained trace header identity differs')
    if row['status'] == 'complete':
        require(row['ideal_rotations'] == expected and row['partial_compiled_resources'] is None,
                'Retained complete circuit is partial')
        validate_resources(row['compiled_resources'], compiler)
        for key in ('error', 'squared_infidelity', 'number_leakage', 'joint_sector_leakage'):
            require(finite(row[key], 'retained ' + key) <= 1 + 1e-8, 'Retained state metric exceeds range')
        finite(row['norm'], 'retained norm', positive=True)
        require(abs(row['norm'] - 1) < 1e-8, 'Preflight state norm drift')
        if order['partition'] == 'conserving':
            require(row['joint_sector_leakage'] < 1e-16, 'Preflight conserving-unit sector leakage')
    else:
        require(row['compiled_resources'] is None, 'Retained censor claims complete resources')
        validate_resources(row['partial_compiled_resources'], compiler)
    return row


def archive_owned(path, archive_root, journal):
    """Retain an interrupted owned artifact before using its original path."""
    if not path.exists():
        return
    require(path.parent == journal.output and not path.is_symlink() and path.is_file(),
            'Only owned direct preflight artifacts may be archived')
    before = artifact(path, journal.output)
    archive_root.mkdir(parents=True, exist_ok=True)
    destination = archive_root / path.name
    require(not destination.exists(), 'Never overwrite retained interrupted artifacts')
    event = dict(original=before, retained_path=str(destination.relative_to(journal.output)), status='planned')
    journal.data['retained_interrupted'].append(event)
    journal.save()
    path.rename(destination)
    after = artifact(destination, journal.output)
    require(after['sha256'] == before['sha256'] and after['bytes'] == before['bytes'],
            'Interrupted artifact changed during retention')
    event.update(status='retained', retained=after)
    journal.save()


def resume(args, started_cpu, started_wall):
    output = args.output.resolve(strict=True)
    project = args.project_root.resolve(strict=True)
    require(output.is_dir() and not output.is_relative_to(project), 'Preflight output must be outside Q')
    record_path = args.record.resolve()
    require(not record_path.is_relative_to(project), 'Continuation record must be outside Q')
    require(not (output / 'PREFLIGHT_COMPLETE.json').exists(), 'Preflight already has a completion artifact')
    # Prevent the caller's record path from claiming an existing or future
    # scientific artifact; its filename remains an explicit continuation label.
    require(record_path.name not in {'model.json', 'reference.npz', 'PREFLIGHT_ORDERS_FROZEN.json',
            'PREFLIGHT_COMPLETE.json'} and not record_path.name.startswith(('order-', 'result-', 'objective-')),
            'Continuation record name collides with preflight artifacts')
    source_dir = project / 'numerical_audit/scalable_ordering'
    sys.path.insert(0, str(source_dir))
    import common as k
    import compiler
    import compilation
    import reference as reference_module
    import numpy as np
    for module in (k, compiler, compilation, reference_module):
        require(Path(module.__file__).resolve().parent == source_dir, 'Wrong continuation source import')
    journal = Journal(record_path, output, started_cpu, started_wall)
    try:
        frozen_path = output / 'PREFLIGHT_ORDERS_FROZEN.json'
        frozen = k.read(frozen_path)
        original_frozen = journal.bind(frozen_path)
        journal.data.update(original_frozen_sha256=original_frozen['sha256'],
                            frozen_source_commit=frozen['source_commit'])
        require(frozen['sources'] == k.source_map() and frozen['environment'] == k.environment(),
                'Frozen preflight sources/environment differ from current execution')
        require(k.committed_sources() == frozen['source_commit'], 'Current committed head differs from frozen head')
        require(tuple(k.ARMS) == ARMS and set(frozen['orders']) == set(ARMS) and
                frozen['scale'] == '7/5' and frozen['T'] == 1.5 and frozen['cap_multiple'] == 30,
                'Frozen preflight panel differs')
        model_path = output / 'model.json'
        model = k.read(model_path)
        require(journal.bind(model_path)['sha256'] == frozen['model_sha256'], 'Frozen preflight model hash differs')
        require(model['n'] == 12 and model['electrons'] == 6 and
                model['case'] == 'preflight_' + frozen['original_case'] + '_scaled7over5' and
                model['preflight']['scale'] == '7/5' and
                model['preflight']['original_case'] == frozen['original_case'], 'Preflight model identity differs')
        k.blocks(model)
        orders, existing = {}, {}
        for arm in ARMS:
            path = output / f'order-{arm}.json'
            require(journal.bind(path)['sha256'] == frozen['orders'][arm], 'Frozen order hash differs: ' + arm)
            order = k.read(path)
            require(order['arm'] == arm and order['model_sha256'] == frozen['model_sha256'] and
                    order['time'] == k.TRAIN_TIME and order['response_is_error_bound'] is False,
                    'Frozen order logical identity differs: ' + arm)
            compilation.validate_record(model, order)
            finite(order['CPU_seconds'], 'ordering CPU')
            bounded_integer(order['evaluations'], 0, 3000, 'ordering evaluations')
            bounded_integer(order['accepted_moves'], 0, 15, 'ordering accepted moves')
            objective_path = output / f'objective-{arm}.json'
            objective = k.read(objective_path)
            journal.bind(objective_path)
            expected_kind = 'first' if arm in ('first', 'pauli') else 'cubic'
            steps = compilation.affordable_steps(order['unit_lengths'], 30)
            require(objective['kind'] == expected_kind and objective['steps'] == steps,
                    'Saved independent objective identity differs: ' + arm)
            finite(objective['value'], 'saved objective')
            finite(objective['raw_response'], 'saved raw response')
            expected_score = objective['raw_response'] if expected_kind == 'first' else objective['raw_response'] / steps**4
            require(abs(objective['value'] - expected_score) <= 1e-20 + 1e-12 * max(objective['value'], expected_score),
                    'Saved objective scaling differs')
            if order['score'] is not None:
                finite(order['score'], 'stored ordering score')
                require(abs(order['score'] - objective['value']) <=
                        1e-20 + 1e-7 * max(order['score'], objective['value']), 'Saved objective/score parity differs')
            orders[arm] = order
            result_path = output / f'result-{arm}.json'
            if result_path.exists():
                row = k.read(result_path)
                journal.bind(result_path)
                existing[arm] = authenticate_row(row, arm, model, order, output / f'{arm}.ndjson.gz',
                                                  journal, k, compiler, compilation)
                journal.data['arm_events'].append(dict(arm=arm, action='reused', status=row['status'],
                    result_sha256=journal.data['reused_artifacts'][result_path.name]['sha256'],
                    repeats_completed_circuit=False))
                journal.save()
        require(all(row['status'] == 'complete' for row in existing.values()),
                'A retained compiler censor cannot become complete under the same frozen bounds; no rerun performed')
        reference_path = output / 'reference.npz'
        journal.bind(reference_path)
        exact = k.load_vector(reference_path, 'exact', model['n'])
        journal.data.update(status='reference_parity', reused_arms=list(existing),
                            missing_arms=[arm for arm in ARMS if arm not in existing])
        journal.save()
        cpu, wall = time.process_time(), time.monotonic()
        reference = reference_module.Reference(model)
        independent = reference.series_action(1.5)
        difference = float(np.max(abs(exact - independent)))
        require(math.isfinite(difference) and difference <= 2e-9, 'Saved preflight reference/independent series parity')
        journal.data['reference_check'] = dict(status='passed', max_difference=difference,
            CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
            saved_reference_sha256=journal.data['reused_artifacts']['reference.npz']['sha256'],
            method='One fresh independent series_action; saved reference loaded through bounded load_vector')
        journal.save()

        archive_root = output / 'retained-interrupted' / (record_path.stem + '-' + uuid.uuid4().hex)
        rows, continued = [], []
        for arm in ARMS:
            if arm in existing:
                rows.append(existing[arm])
                continue
            journal.data.update(status='continuing_arm', active_arm=arm)
            journal.save()
            trace = output / f'{arm}.ndjson.gz'
            result_path = output / f'result-{arm}.json'
            require(not result_path.exists(), 'Missing result appeared during continuation: ' + arm)
            archive_owned(trace, archive_root, journal)
            archive_owned(result_path.with_suffix('.json.partial'), archive_root, journal)
            cpu, wall = time.process_time(), time.monotonic()
            order = orders[arm]
            result = compilation.compile_candidate(model, order, 1.5, 30, trace,
                                                    dict(case=model['case'], arm=arm, mode='preflight'))
            row = dict(arm=arm, ordering_CPU_seconds=order['CPU_seconds'], ordering_status=order['status'],
                       ordering_evaluations=order['evaluations'], **result)
            if result['status'] == 'complete':
                actual, counts = reference.audit_trace(trace, result['identity'], compiler.VERSION)
                counts.pop('actual_Rz_rotations')
                require(counts == result['compiled_resources'], 'Continuation independent DAG parity')
                ideal = reference.ideal_formula(order, result['steps'], 1.5)
                overlap = np.vdot(ideal, actual)
                require(abs(overlap) > 0 and
                        float(np.max(abs(actual - ideal * overlap / abs(overlap)))) < 2e-9,
                        'Continuation actual primitive/ideal parity')
                row.update(reference_module.metrics(exact, actual, model))
            require(not result_path.exists(), 'Never overwrite a result during continuation')
            k.write(result_path, row)
            authenticate_row(row, arm, model, order, trace, journal, k, compiler, compilation, new=True)
            rows.append(row)
            continued.append(arm)
            journal.bind(result_path, label='new/' + result_path.name)
            journal.data['arm_events'].append(dict(arm=arm, action='completed_in_continuation', status=row['status'],
                result_sha256=journal.data['reused_artifacts']['new/' + result_path.name]['sha256'],
                own_arm_CPU_seconds=time.process_time() - cpu, own_arm_wall_seconds=time.monotonic() - wall,
                repeats_completed_circuit=False))
            journal.save()
            print('PREFLIGHT_CONTINUED', arm, row['status'], flush=True)
            require(row['status'] == 'complete', 'Continuation retained a compiler censor; no completion claim')

        require(len(rows) == 6 and {row['arm'] for row in rows} == set(ARMS) and
                all(row['status'] == 'complete' for row in rows), 'All six complete preflight arms required')
        require(k.sha(frozen_path) == original_frozen['sha256'] and
                frozen['sources'] == k.source_map() and frozen['environment'] == k.environment(),
                'Frozen identity changed during continuation')
        require(k.sha(model_path) == frozen['model_sha256'] and
                all(k.sha(output / f'order-{arm}.json') == frozen['orders'][arm] for arm in ARMS),
                'Frozen model/order bytes changed during continuation')
        archive_owned(output / 'PREFLIGHT_COMPLETE.json.partial', archive_root, journal)
        require(not (output / 'PREFLIGHT_COMPLETE.json').exists(), 'Never overwrite a preflight completion')
        result = dict(status='passed', mode='altered-input preflight', original_case=frozen['original_case'],
            model_sha256=frozen['model_sha256'], arms=rows,
            max_trace_bytes=max(row['trace_bytes'] for row in rows),
            CPU_seconds=time.process_time() - started_cpu, wall_seconds=time.monotonic() - started_wall,
            main_cases_evaluated=False, main_errors_exposed=False,
            same_frozen_orders=True, new_training=False, continuation=True,
            original_frozen_sha256=original_frozen['sha256'], source_commit=frozen['source_commit'],
            continuation_record=str(record_path), continuation_script_sha256=journal.data['renderer']['sha256'],
            reused_arms=list(existing), continued_arms=continued,
            old_attempt_costs_excluded=True, cost_scope=journal.data['cost_scope'],
            note='Resource feasibility only; same frozen preflight model/orders and source. '
                 'Original attempt costs remain retained separately. Do not tune the main protocol '
                 'using these altered-input accuracy values.')
        k.write(output / 'PREFLIGHT_COMPLETE.json', result)
        journal.bind(output / 'PREFLIGHT_COMPLETE.json', label='new/PREFLIGHT_COMPLETE.json')
        journal.data.update(status='passed', active_arm=None, continued_arms=continued,
                            complete_arms=6, main_errors_exposed=False)
        journal.save()
        print('PREFLIGHT_CONTINUATION_PASS', frozen['original_case'], len(rows), flush=True)
        return result
    except BaseException as error:
        journal.data.update(status='failed', exception=dict(type=type(error).__name__, message=str(error)))
        journal.save()
        raise


def main():
    cpu, wall = time.process_time(), time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True, help='Fresh continuation journal JSON')
    parser.add_argument('--project-root', type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    resume(args, cpu, wall)


if __name__ == '__main__':
    main()
