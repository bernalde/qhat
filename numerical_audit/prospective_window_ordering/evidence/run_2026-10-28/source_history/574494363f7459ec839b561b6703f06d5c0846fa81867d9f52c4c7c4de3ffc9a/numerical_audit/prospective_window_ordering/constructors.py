"""One fresh, bounded normal-API or unchanged-control construction.

The external owning cgroup is authoritative for complete CPU, including child
startup/JIT, final serialization and process teardown. No physical error is read.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
import time
import study_common as k


def construct(model_path, digest, arm, output, *, repetition, source_freeze, julia_executable=None):
    wall = time.monotonic()
    k.require(arm in k.ARMS and not output.exists(), 'Known arm and fresh output required')
    k.require(type(repetition) is int and 1 <= repetition <= k.REPETITIONS, 'Frozen repetition identity')
    k.require(k.sha(model_path) == digest, 'Frozen model bytes changed')
    frozen = k.read(source_freeze)
    k.require(frozen['sources'] == k.source_map(), 'Scientific source differs from constructor freeze')
    source_freeze_sha256 = k.sha(source_freeze)
    model = k.read(model_path)
    units = k.blocks(model)
    ordering = k.bootstrap()
    canonical = list(range(len(units)))
    magnitude = sorted(canonical, key=lambda index:
                       (-sum(abs(complex(value)) for value in units[index].values()), index))
    partition = 'pauli' if arm == 'pauli' else 'conserving'
    physical_units = k.blocks(model, pauli=partition == 'pauli')
    result = dict(order=list(range(len(physical_units))), score=None, evaluations=0,
                  accepted_moves=0, status='complete', partition=partition)
    # Reserve final order/layout serialization inside the SAME absolute cap.
    # process_time is lifetime process CPU here, so Python startup/import work
    # already consumes the allowance. Julia child CPU is additionally charged
    # by the normal bridge and, authoritatively, the owning cgroup.
    deadline = k.TRAIN_CPU_SECONDS - 0.5
    initializer_cpu = None
    score_kind = None
    try:
        if arm == 'canonical':
            result['order'] = canonical
        elif arm == 'magnitude':
            result['order'] = magnitude
        else:
            control = None if arm == 'window_magnitude' else k.retained_training()
            if time.process_time() >= deadline:
                raise ordering.SearchDeadline('startup exhausted the construction allowance')
            began = time.process_time()
            if arm == 'window_magnitude':
                placed = magnitude
            else:
                placed = control.current(model, deadline)
            initializer_cpu = time.process_time() - began
            if partition == 'conserving':
                result['order'] = list(placed)
            if arm in ('first', 'pauli'):
                result = control.first_search(model, arm, placed, deadline)
                score_kind = 'retained_first_order_response'
            elif arm in ('cubic', 'window_sparse', 'window_magnitude'):
                ordering = k.bootstrap()
                remaining = deadline - time.process_time()
                if remaining <= 0:
                    result['status'] = 'cpu_limit'
                else:
                    options = {}
                    if arm.startswith('window_'):
                        from qhat.analysis.window_ordering import WindowOptions
                        parent = 'sparse' if arm == 'window_sparse' else 'magnitude'
                        options = dict(method='window', window_options=WindowOptions(
                            parent_method=parent, width=8, max_windows=1, max_sweeps=10,
                            julia_executable=julia_executable,
                            wall_seconds=max(0.001, k.TRAIN_WALL_SECONDS - (time.monotonic() - wall))))
                    starts = [placed] if arm == 'window_magnitude' else [placed, list(reversed(placed)), canonical]
                    trained = ordering.train_ordering(units, model['n'], (1 << model['electrons']) - 1,
                        time=k.TRAIN_TIME, rotation_budget=k.TRAIN_CAP * len(model['coefficients']),
                        cpu_seconds=remaining, max_evaluations=3000, max_moves=5,
                        max_entries=50000, max_cached_actions=100000,
                        starts=starts, **options)
                    result = asdict(trained)
                    result['score_components'] = result['score']
                    result['score'] = None if trained.score is None else trained.score.value
                    result['partition'] = partition
                    score_kind = 'diagonal_reference_cubic_response'
            else:
                k.require(arm == 'current', 'Unhandled frozen arm')
    except Exception as error:
        # A complete permutation was established before expensive construction.
        # It remains an explicitly unscored fallback if setup cannot finish.
        recognized = isinstance(error, (ValueError, ordering.ResourceLimit, ordering.SearchDeadline))
        if type(error).__name__ == 'Deadline':
            recognized = True
        result['status'] = ('cpu_limit' if type(error).__name__ in ('Deadline', 'SearchDeadline')
                            else 'control_or_backend_refusal' if recognized else 'unexpected_error')
        result['reason'] = dict(type=type(error).__name__, message=str(error))
    order = list(result['order'])
    k.require(all(type(index) is int for index in order) and
              sorted(order) == list(range(len(physical_units))), 'Complete original-label witness')
    from qhat.analysis.strang_schedule import plan_strang_schedule
    lengths = [len(physical_units[index]) for index in order]
    plan = plan_strang_schedule(lengths, k.TRAIN_CAP * len(model['coefficients']))
    k.require(plan.feasible, 'Returned initializer cannot afford the training Strang step')
    if isinstance(result.get('score'), (int, float)):
        k.require(math.isfinite(result['score']) and result['score'] >= 0, 'Finite nonnegative control score')
    # The model, order, unit lengths and hash specify the original word layout;
    # do not duplicate every word/primitive trace in every repetition file.
    layout = [[[[site, letter] for site, letter in word], str(value)]
              for index in order for word, value in sorted(physical_units[index].items())]
    record = dict(schema='prospective-window-constructor-v1', artifact_complete=True,
        arm=arm, repetition=repetition, model_sha256=digest, case=model['case'], partition=partition,
        source_freeze_sha256=source_freeze_sha256,
        training_time=k.TRAIN_TIME, training_rotation_budget=k.TRAIN_CAP * len(model['coefficients']),
        constructor=result, score_kind=score_kind, original_word_count=len(model['coefficients']),
        unit_lengths=lengths, pulse_layout_sha256=hashlib.sha256(
            json.dumps(layout, separators=(',', ':')).encode()).hexdigest(),
        initializer_CPU_seconds=initializer_cpu, process_CPU_before_output=time.process_time(),
        wall_seconds_before_output=time.monotonic() - wall,
        total_CPU_ceiling=k.TRAIN_CPU_SECONDS,
        complete_cost_source='Owning process-group receipt, including child/startup/JIT/serialization/teardown',
        response_is_error_bound=False, physical_error_read=False)
    k.write_new(output, record)
    print('FROZE_CONSTRUCTOR', arm, result['status'], flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=k.Path)
    parser.add_argument('sha256')
    parser.add_argument('arm', choices=k.ARMS)
    parser.add_argument('output', type=k.Path)
    parser.add_argument('--repetition', type=int, choices=range(1, k.REPETITIONS + 1), required=True)
    parser.add_argument('--julia-executable')
    parser.add_argument('--source-freeze', type=k.Path, required=True)
    args = parser.parse_args()
    record = construct(args.model, args.sha256, args.arm, args.output,
                       repetition=args.repetition, source_freeze=args.source_freeze,
                       julia_executable=args.julia_executable)
    if record['constructor']['status'] == 'unexpected_error':
        raise SystemExit(2)
