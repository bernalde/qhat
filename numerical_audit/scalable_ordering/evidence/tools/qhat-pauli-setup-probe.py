#!/usr/bin/env python3
"""Resource-only 32-CPU-second Pauli construction on the frozen altered CH4.

Imports the same reference/SciPy module as study.py before training, but never
allocates Reference or propagates a state. Calls training.construct exactly once.
All previous inputs remain read-only; use the coordinator's external guard.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.dont_write_bytecode = True
PROJECT_ROOT = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def artifact(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'Expected regular artifact: ' + str(path))
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=value.hexdigest())


def save(output, report, cpu, wall):
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    report.update(CPU_seconds_since_entry=time.process_time() - cpu,
        wall_seconds_since_entry=time.monotonic() - wall,
        whole_outer_CPU_seconds_before_metadata_write=time.process_time() + children.ru_utime + children.ru_stime,
        max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path = output / 'PROBE.json'
    temporary = path.with_suffix('.json.partial')
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    temporary.replace(path)


def main():
    cpu, wall = time.process_time(), time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Fresh directory outside Q and old preflight')
    parser.add_argument('--project-root', type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    model_path, project, output = args.model.resolve(strict=True), args.project_root.resolve(strict=True), args.output.resolve()
    require(not output.exists(), 'Fresh output directory required')
    for root in (project, model_path.parent):
        require(not output.is_relative_to(root) and not root.is_relative_to(output), 'Probe output must be outside inputs/Q')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='qhat-pauli-setup-resource-probe-v1', status='importing_reference',
        started_UTC=datetime.now(timezone.utc).isoformat(), probe=artifact(Path(__file__).resolve()),
        arm='pauli', expected_nonzero_pair_entries=588360, nominal_constructor_CPU_seconds=32.,
        reference_allocated=False, physical_errors_read=False, main_errors_exposed=False,
        exact_training_construct_calls=0, completed_training_construct_calls=0,
        RSS_units='KiB on Linux; lifetime process high-water mark',
        cost_scope='Whole outer CPU includes waited child CPU and startup through the current metadata '
            'write point. Final metadata writing/interpreter shutdown and the parent guard are excluded; '
            'the external task ledger owns the complete charge. Constructor includes its public cleanup charge.')
    save(output, report, cpu, wall)
    try:
        source_dir = project / 'numerical_audit/scalable_ordering'
        sys.path.insert(0, str(source_dir))
        # Match study.py's loaded numerical environment before timing construction.
        import reference as reference_module
        import common as k
        import training
        for module in (reference_module, k, training):
            require(Path(module.__file__).resolve().parent == source_dir, 'Wrong project-module import')
        sources, environment, head = k.source_map(), k.environment(), k.committed_sources()
        require(k.TRAIN_CPU_SECONDS == 32. and training.MAX_STORED_PAIRS == 750000,
                'Commit the declared 32-second/one-dictionary Pauli setup first')
        frozen_path = model_path.parent / 'PREFLIGHT_ORDERS_FROZEN.json'
        frozen = k.read(frozen_path)
        model = k.read(model_path)
        model_binding, frozen_binding = artifact(model_path), artifact(frozen_path)
        require(model_binding['sha256'] == frozen['model_sha256'] and
                frozen['original_case'] == 'CH4_r1p09' and frozen['scale'] == '7/5' and
                frozen['environment'] == environment and model['n'] == 12 and model['electrons'] == 6 and
                model['case'] == 'preflight_CH4_r1p09_scaled7over5' and len(model['coefficients']) == 1722,
                'Frozen altered-CH4 model identity differs')
        require(model['preflight']['scale'] == '7/5', 'Expected retained 7/5 altered input')
        report.update(status='constructing', source_commit=head, source_sha256=sources, environment=environment,
            model=model_binding, input_freeze=frozen_binding, original_frozen_source_commit=frozen['source_commit'],
            imported_reference_before_training=True, import_and_authentication_CPU_seconds=time.process_time() - cpu,
            exact_training_construct_calls=1,
            pair_storage='Single encoded-pair dictionary with alias, bounded at 750000; external 512 MiB unchanged')
        save(output, report, cpu, wall)
        call_cpu, call_wall = time.process_time(), time.monotonic()
        record = training.construct(model_path, model_binding['sha256'], 'pauli')
        caller_cpu, caller_wall = time.process_time() - call_cpu, time.monotonic() - call_wall
        record_path = output / 'order-pauli.json'
        k.write(record_path, record)
        report.update(status='checking_completed_seed', order_record=artifact(record_path),
            completed_training_construct_calls=1,
            standalone_call_CPU_seconds=caller_cpu, standalone_call_wall_seconds=caller_wall,
            constructor_CPU_seconds=record['CPU_seconds'], constructor_wall_seconds=record['wall_seconds'],
            ordering_status=record['status'], evaluations=record['evaluations'], accepted_moves=record['accepted_moves'],
            score_completed=record['score'] is not None, stats=record.get('stats', {}))
        save(output, report, cpu, wall)
        require(record['arm'] == 'pauli' and record['partition'] == 'pauli' and
                record['model_sha256'] == model_binding['sha256'] and
                sorted(record['order']) == list(range(1722)) and record['unit_lengths'] == [1] * 1722,
                'Returned Pauli layout differs from the frozen target')
        stats = record.get('stats', {})
        require(stats.get('pair_setup_completed') is True and stats.get('pair_entries') == 588360,
                'Pauli setup did not complete with the independently counted support')
        require(type(stats.get('seed_evaluations')) is int and stats['seed_evaluations'] >= 1 and
                record['score'] is not None and math.isfinite(record['score']) and record['score'] >= 0,
                'No completed Pauli seed was demonstrated')
        require(record['evaluations'] == stats['seed_evaluations'] + stats['candidate_evaluations'],
                'Seed/candidate accounting differs')
        require(k.sha(model_path) == model_binding['sha256'] and k.sha(frozen_path) == frozen_binding['sha256'] and
                sources == k.source_map() and environment == k.environment() and
                subprocess.check_output(['git', '-C', str(project), 'rev-parse', 'HEAD'], text=True).strip() == head,
                'Probe model/freeze/source/environment changed during construction')
        report.update(status='passed', demonstrated_completed_seed=True, reference_allocated=False,
                      physical_errors_read=False, main_errors_exposed=False)
        save(output, report, cpu, wall)
        print(json.dumps(dict(status='passed', ordering_status=record['status'],
            pair_entries=stats['pair_entries'], seed_evaluations=stats['seed_evaluations'],
            candidate_evaluations=stats['candidate_evaluations'], score_completed=True,
            constructor_CPU_seconds=record['CPU_seconds'], max_RSS_KiB=report['max_RSS_KiB']),
            sort_keys=True), flush=True)
    except BaseException as error:
        report.update(status='failed', exception=dict(type=type(error).__name__, message=str(error)))
        save(output, report, cpu, wall)
        raise


if __name__ == '__main__':
    main()
