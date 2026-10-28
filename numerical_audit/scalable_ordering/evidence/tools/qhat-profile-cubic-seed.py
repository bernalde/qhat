#!/usr/bin/env python3
"""Resource-only calibration of one authenticated sparse cubic seed."""
import argparse
import gc
import hashlib
import importlib
import json
import math
from pathlib import Path
import resource
import sys
import tempfile
import os
import time
import traceback
import weakref


ROOT = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')
STUDY = ROOT / 'numerical_audit/scalable_ordering'
MAX_INPUT_BYTES = 32 * 1024**2
MAX_ENTRIES = 50000
MAX_CACHED_ACTIONS = 50000
INITIALIZER_CPU_NS = 700_000_000


def read_authenticated(path):
    path = path.resolve(strict=True)
    with path.open('rb') as stream:
        raw = stream.read(MAX_INPUT_BYTES + 1)
    if len(raw) > MAX_INPUT_BYTES:
        raise ValueError('input exceeds byte bound: ' + str(path))
    return raw, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw))


def save_owned(path, record):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent,
                prefix=path.name + '.', suffix='.partial', delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def peak_rss():
    return dict(process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                scope='Linux process lifetime peak; not an isolated phase allocation')


def run(args, record, output):
    phase_start = None

    def start(name):
        nonlocal phase_start
        record['current_phase'] = name
        record['phases'][name] = dict(status='running', memory_before=peak_rss())
        save_owned(output, record)
        phase_start = (time.process_time_ns(), time.perf_counter_ns())

    def finish(name, *, stopped=None, **extra):
        stopped = stopped or (time.process_time_ns(), time.perf_counter_ns())
        cpu_ns = stopped[0] - phase_start[0]
        wall_ns = stopped[1] - phase_start[1]
        record['phases'][name].update(status='complete', CPU_nanoseconds=cpu_ns,
            wall_nanoseconds=wall_ns, CPU_seconds=cpu_ns/1e9, wall_seconds=wall_ns/1e9,
            memory_after=peak_rss(), **extra)
        save_owned(output, record)

    context = None
    try:
        start('input_bootstrap')
        model_raw, model_identity = read_authenticated(args.model)
        order_raw, order_identity = read_authenticated(args.order)
        model, frozen_order = json.loads(model_raw), json.loads(order_raw)
        sys.path.insert(0, str(STUDY))
        k = importlib.import_module('common')
        if Path(k.__file__).resolve() != STUDY / 'common.py':
            raise ValueError('wrong study common module')
        sources = k.source_map()
        ordering = k.bootstrap()
        planner = importlib.import_module('qhat.analysis.strang_schedule')
        if (Path(ordering.__file__).resolve() != ROOT / 'analysis/sparse_response_ordering.py'
                or Path(planner.__file__).resolve() != ROOT / 'analysis/strang_schedule.py'):
            raise ValueError('wrong core/planner module import')
        if k.source_map() != sources:
            raise ValueError('scientific source changed during bootstrap')
        environment = k.environment()
        blocks = k.blocks(model)
        n, electrons = model['n'], model['electrons']
        if type(electrons) is not int or not 0 <= electrons <= n:
            raise ValueError('invalid HF electron count')
        order = frozen_order['order']
        if (not isinstance(order, list) or any(type(i) is not int for i in order)
                or sorted(order) != list(range(len(blocks)))):
            raise ValueError('order is not a complete unit permutation')
        if (frozen_order.get('model_sha256') != model_identity['sha256']
                or frozen_order.get('arm') != 'current'
                or frozen_order.get('partition') != 'conserving'
                or frozen_order.get('formula') != 2):
            raise ValueError('current order/model identity differs')
        raw_blocks = [k.decode(block) for block in model['blocks']]
        expected_schedule = [[[[q, p] for q, p in word], str(c)]
                             for i in order for word, c in sorted(raw_blocks[i].items())]
        if (frozen_order.get('schedule') != expected_schedule
                or frozen_order.get('unit_lengths') != [len(raw_blocks[i]) for i in order]):
            raise ValueError('frozen current layout differs from supplied Hamiltonian')
        budget, hf = 30 * len(model['coefficients']), (1 << electrons) - 1
        record.update(model=model_identity, order=order_identity, source_map=sources,
            profiler_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            environment=environment, module_paths=dict(common=str(Path(k.__file__).resolve()),
                core=str(Path(ordering.__file__).resolve()), planner=str(Path(planner.__file__).resolve())),
            logical_input=dict(case=model.get('case'), n_qubits=n, electrons=electrons,
                Pauli_words=len(model['coefficients']), units=len(blocks), reference=hf,
                source_order_status=frozen_order.get('status'), time=1.0, rotation_budget=budget),
            support_bounds=dict(max_entries=MAX_ENTRIES, max_cached_actions=MAX_CACHED_ACTIONS,
                supplied_reference_components=1,
                scope='Explicit reachable-vector/cache entry caps; not byte-size estimates. RAM/CPU/wall enforcement belongs to the external OS guard.'))
        # Release parsing-only copies before context work; their cleanup belongs
        # to input/setup rather than being omitted from the outer measurement.
        del raw_blocks, expected_schedule, model_raw, order_raw, frozen_order, model
        finish('input_bootstrap')

        start('context_construction')
        context = ordering.SparseResponse(blocks, n, hf, time=1.0,
            max_entries=MAX_ENTRIES, max_cached_actions=MAX_CACHED_ACTIONS, deadline=math.inf)
        stopped = (time.process_time_ns(), time.perf_counter_ns())
        finish('context_construction', stopped=stopped, stats=dict(context.stats))

        start('first_seed_score')
        score = context.score(order, budget, deadline=math.inf)
        stopped = (time.process_time_ns(), time.perf_counter_ns())
        record['score'] = dict(status='complete', value=score.value,
            cubic_response=score.cubic_response, steps=score.steps,
            response_is_error_bound=False,
            scope='Cubic ranking proxy only; no physical validation error was evaluated.')
        record['support_stats'] = dict(context.stats)
        record['retained_support'] = dict(action_cache_entries=context.cache_entries,
            action_cache_columns=len(context.cache),
            reference_column_entries=sum(len(column) for column in context.reference_cubic.values()),
            last_cubic_column_entries=0 if context.last_columns is None else
                sum(len(column) for column in context.last_columns.values()))
        finish('first_seed_score', stopped=stopped)

        start('complete_context_cleanup')
        context_weak = weakref.ref(context)
        del context, blocks, score
        collected = gc.collect()
        if context_weak() is not None:
            raise RuntimeError('context remains live after cleanup')
        finish('complete_context_cleanup', context_collected=True,
               cyclic_objects_collected=collected, explicit_gc_included=True)

        # Apply the prospective resource-only rule to measured integer clock
        # values. This proposes a cap; it changes no source/global/protocol.
        context_ns = record['phases']['context_construction']['CPU_nanoseconds']
        score_ns = record['phases']['first_seed_score']['CPU_nanoseconds']
        required_ns = 4 * (context_ns + score_ns + INITIALIZER_CPU_NS)
        candidate = 2
        while candidate * 1_000_000_000 < required_ns:
            candidate *= 2
        record['resource_calibration'] = dict(
            rule='Smallest power of two >=4*(context CPU + first seed CPU +0.7 initializer CPU), minimum2, maximum32CPU seconds.',
            initializer_CPU_seconds=.7, required_nanoseconds=required_ns,
            uncapped_next_power_of_two_CPU_seconds=candidate,
            status='resource_question' if candidate > 32 else 'within_declared_range',
            proposed_shared_searched_arm_CPU_seconds=None if candidate > 32 else candidate,
            applies_to=['first', 'pauli', 'cubic'], source_or_protocol_changed=False,
            scope='Resource-only candidate allowance; not proof that three seeds or a search move will finish. Actual complete constructor CPU remains separately charged.')
        if k.source_map() != sources or k.environment() != environment:
            raise ValueError('source/environment changed during profile')
        for identity in (model_identity, order_identity):
            if hashlib.sha256(Path(identity['path']).read_bytes()).hexdigest() != identity['sha256']:
                raise ValueError('model/order file changed during profile')
        record.update(status='passed', current_phase='complete')
    except BaseException as error:
        record.update(status='failed', error_type=type(error).__name__, error=str(error),
                      traceback=traceback.format_exc())
        if phase_start is not None:
            phase = record['phases'][record['current_phase']]
            if phase['status'] == 'running':
                phase.update(status='failed', CPU_seconds=(time.process_time_ns()-phase_start[0])/1e9,
                             wall_seconds=(time.perf_counter_ns()-phase_start[1])/1e9,
                             memory_after=peak_rss())
        save_owned(output, record)
        raise


def main():
    outer_cpu, outer_wall = time.process_time_ns(), time.perf_counter_ns()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True, type=Path)
    parser.add_argument('--order', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    record = dict(status='initializing', current_phase='not_started', phases={},
        physical_validation_errors_read=False, other_orders_scored=False,
        internal_deadline='infinity; root external OS guard owns the finite allowance',
        scope='One current-order cubic seed only; imports/input, context, score and explicit cleanup are separate.')
    with output.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    try:
        run(args, record, output)
    finally:
        # run() has returned/unwound: its local object cleanup is included.
        record['outer'] = dict(CPU_seconds=(time.process_time_ns()-outer_cpu)/1e9,
            wall_seconds=(time.perf_counter_ns()-outer_wall)/1e9,
            process_CPU_seconds_at_completion=time.process_time(), memory=peak_rss(),
            boundary='main entry through run() return/unwind, including input/imports, phase writes and local cleanup; final JSON serialization and process shutdown excluded. External guard retains whole-process accounting.')
        save_owned(output, record)
    print('CUBIC_SEED_PROFILE', record['status'], output, flush=True)


if __name__ == '__main__':
    main()
