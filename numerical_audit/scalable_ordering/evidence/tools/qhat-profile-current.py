#!/usr/bin/env python3
"""Paired serial current()-only profiling; never evaluates physical errors."""
import argparse
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time
import traceback


COMMON_PATH = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree/'
                   'numerical_audit/scalable_ordering/common.py')
MAX_SOURCE_BYTES = 4 * 1024**2
MAX_MODEL_BYTES = 32 * 1024**2
DISPATCH = (('old', 'new'), ('new', 'old'), ('old', 'new'))


def authenticated_bytes(path, maximum):
    path = path.resolve(strict=True)
    with path.open('rb') as stream:
        raw = stream.read(maximum + 1)
    if len(raw) > maximum:
        raise ValueError('input exceeds byte bound: ' + str(path))
    return raw, dict(path=str(path), bytes=len(raw),
                     sha256=hashlib.sha256(raw).hexdigest())


def save_owned(path, record):
    """Atomically update only the output exclusively created by this run."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent,
                                         prefix=path.name + '.', suffix='.partial',
                                         delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load_authenticated_module(label, path, raw):
    """Use importlib identity while executing the exact hashed source bytes."""
    name = '_qhat_current_profile_' + label
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError('cannot create source module: ' + str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    # Bypass timestamp-based bytecode caches: the authenticated bytes execute.
    exec(compile(raw, str(path), 'exec'), module.__dict__)
    if not callable(getattr(module, 'current', None)):
        raise ValueError('source has no current() function: ' + str(path))
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old', required=True, type=Path)
    parser.add_argument('--new', required=True, type=Path)
    parser.add_argument('--model', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    # Own the requested fresh output before performing imports or any calls.
    with output.open('x') as stream:
        json.dump({'status': 'initializing'}, stream)
        stream.write('\n')
    record = dict(status='initializing',
        measured_boundary='current(model, deadline=inf) call only',
        separate_costs='Source/model reading, parsing, imports, validation and JSON writes are outside each measured call.',
        scope='Initializer performance and exact returned-order equality only; no full-training or accuracy claim.',
        physical_error_files_read=False, dispatch=[list(pair) for pair in DISPATCH],
        attempts=[], pairs=[])
    current_attempt = None
    try:
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                     'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
            os.environ[name] = '1'
        before_cpu, before_wall = time.process_time_ns(), time.perf_counter_ns()
        old_raw, old_identity = authenticated_bytes(args.old, MAX_SOURCE_BYTES)
        new_raw, new_identity = authenticated_bytes(args.new, MAX_SOURCE_BYTES)
        model_raw, model_identity = authenticated_bytes(args.model, MAX_MODEL_BYTES)
        common_raw, common_identity = authenticated_bytes(COMMON_PATH, MAX_SOURCE_BYTES)
        model = json.loads(model_raw)
        if not isinstance(model, dict) or not isinstance(model.get('blocks'), list):
            raise ValueError('model must contain a blocks list')
        expected_order = list(range(len(model['blocks'])))
        record.update(sources={'old': old_identity, 'new': new_identity,
                               'common': common_identity}, model=model_identity,
                      n_qubits=model.get('n'), unit_count=len(expected_order),
                      input_preparation=dict(CPU_seconds=(time.process_time_ns()-before_cpu)/1e9,
                                             wall_seconds=(time.perf_counter_ns()-before_wall)/1e9))
        save_owned(output, record)

        before_cpu, before_wall = time.process_time_ns(), time.perf_counter_ns()
        sys.path.insert(0, str(COMMON_PATH.parent))
        common = importlib.import_module('common')
        if Path(common.__file__).resolve() != COMMON_PATH.resolve():
            raise ValueError('wrong common module import')
        # The unchanged common dependency is also bound to retained bytes.
        if COMMON_PATH.read_bytes() != common_raw:
            raise ValueError('common source changed during import')
        modules = {
            'old': load_authenticated_module('old', args.old.resolve(), old_raw),
            'new': load_authenticated_module('new', args.new.resolve(), new_raw),
        }
        record['imports'] = dict(CPU_seconds=(time.process_time_ns()-before_cpu)/1e9,
                                 wall_seconds=(time.perf_counter_ns()-before_wall)/1e9,
                                 python=sys.version, executable=sys.executable)
        record['status'] = 'running'
        save_owned(output, record)

        reference_order = None
        for repetition, pair in enumerate(DISPATCH, start=1):
            completed = {}
            for label in pair:
                # Parse outside the measured boundary; each call receives fresh,
                # identical authenticated data rather than another call's object.
                call_model = json.loads(model_raw)
                current_attempt = dict(repetition=repetition, source=label,
                                       dispatch_position=len(completed)+1, status='running')
                record['attempts'].append(current_attempt)
                save_owned(output, record)
                wall_start = time.perf_counter_ns()
                cpu_start = time.process_time_ns()
                try:
                    order = modules[label].current(call_model, deadline=math.inf)
                finally:
                    cpu_elapsed = time.process_time_ns() - cpu_start
                    wall_elapsed = time.perf_counter_ns() - wall_start
                    current_attempt.update(CPU_seconds=cpu_elapsed/1e9,
                                           wall_seconds=wall_elapsed/1e9)
                if type(order) not in (list, tuple) or any(type(i) is not int for i in order):
                    raise ValueError('current() returned an invalid order type')
                order = list(order)
                current_attempt['order'] = order
                if sorted(order) != expected_order:
                    raise ValueError('current() returned an incomplete permutation')
                if reference_order is None:
                    reference_order = list(order)
                elif order != reference_order:
                    raise ValueError('exact returned-order mismatch across source/repetition')
                current_attempt.update(status='complete', exact_order_equal=True)
                completed[label] = current_attempt
                save_owned(output, record)
                current_attempt = None
            if completed['new']['CPU_seconds'] <= 0 or completed['new']['wall_seconds'] <= 0:
                raise ValueError('timer resolution cannot support a paired ratio')
            record['pairs'].append(dict(repetition=repetition, dispatch=list(pair),
                exact_order_equal=completed['old']['order'] == completed['new']['order'],
                CPU_old_over_new=completed['old']['CPU_seconds']/completed['new']['CPU_seconds'],
                wall_old_over_new=completed['old']['wall_seconds']/completed['new']['wall_seconds']))
            save_owned(output, record)

        # Authenticate paths again; executed modules used retained source bytes.
        for identity in (*record['sources'].values(), record['model']):
            if hashlib.sha256(Path(identity['path']).read_bytes()).hexdigest() != identity['sha256']:
                raise ValueError('authenticated input changed during campaign: ' + identity['path'])
        record.update(status='passed', exact_order_equal_across_all_attempts=True,
                      returned_order=reference_order,
                      summary=dict(repetitions=3, calls=6,
                          ratio_definition='Paired old/new ratios computed per repetition before taking the median; values above1 mean new is faster.',
                          median_paired_CPU_old_over_new=statistics.median(
                              pair['CPU_old_over_new'] for pair in record['pairs']),
                          median_paired_wall_old_over_new=statistics.median(
                              pair['wall_old_over_new'] for pair in record['pairs']),
                          median_CPU_seconds={label: statistics.median(
                              a['CPU_seconds'] for a in record['attempts'] if a['source'] == label)
                              for label in ('old', 'new')},
                          median_wall_seconds={label: statistics.median(
                              a['wall_seconds'] for a in record['attempts'] if a['source'] == label)
                              for label in ('old', 'new')}))
        save_owned(output, record)
    except BaseException as error:
        if current_attempt is not None:
            current_attempt.update(status='failed', error_type=type(error).__name__,
                                   error=str(error))
        record.update(status='failed', error_type=type(error).__name__, error=str(error),
                      traceback=traceback.format_exc())
        save_owned(output, record)
        raise
    print('CURRENT_PROFILE_PASS', output, flush=True)


if __name__ == '__main__':
    main()
