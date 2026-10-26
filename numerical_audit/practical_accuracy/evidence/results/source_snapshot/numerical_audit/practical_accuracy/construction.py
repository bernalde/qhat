"""Caller-through-return ordering and exact primitive-circuit construction.

No physical error is available to either constructor. Nominal order work gets
half the total cap; actual eligibility uses the full measured construction sum.
"""
from fractions import Fraction
import time
from common import (METHODS, TIME, GRID, require, read, sha, cubic, compiler, placement)

CLEANUP_RESERVATION = .02
EVALUATION_CAP = 3000
MAX_MOVES = 5


def _search(model, partition, starts, deadline, phases, key, clock):
    began = clock()
    context = (cubic.Context(model) if partition == 'cubic' else
               cubic.FirstOrderContext(model, TIME, partition))
    phases[key + '_context_CPU_seconds'] = clock() - began
    seeds = []
    result = None
    status = 'complete seeds; no search allowance'
    began = clock()
    try:
        for order in starts:
            seeds.append(context.score(order, TIME, 30 * len(model['coefficients']),
                                       deadline - CLEANUP_RESERVATION, clock))
    except cubic.Interrupted as interrupted:
        status = 'incomplete seed evaluation'
        partial_stages = interrupted.completed_stages
    else:
        partial_stages = None
    phases[key + '_seed_CPU_seconds'] = clock() - began
    began = clock()
    remaining = deadline - clock() - CLEANUP_RESERVATION
    if len(seeds) == 3 and remaining > .025:
        result = cubic.insertion(context, seeds, TIME, 30 * len(model['coefficients']),
                                 min(8., remaining), EVALUATION_CAP, MAX_MOVES, clock)
        best = result['record']
        status = 'bounded search returned'
    else:
        best = min(seeds, key=lambda r: (r['score'], tuple(r['order']))) if seeds else None
    phases[key + '_search_CPU_seconds'] = clock() - began
    began = clock()
    del context
    phases[key + '_cleanup_CPU_seconds'] = clock() - began
    return dict(best=best, seeds=seeds, search=result, status=status,
                interrupted_seed_stages=partial_stages)


def _ordering(model_path, model_digest, method, deadline, clock):
    phases = {}
    began = clock()
    require(sha(model_path) == model_digest, 'Constructor input digest differs')
    model = read(model_path)
    blocks = cubic.validate_model(model)
    phases['model_load_validation_CPU_seconds'] = clock() - began
    began = clock()
    diagonal = sorted((i for i, b in enumerate(blocks) if cubic.engine.descriptor(b)[0] == 0),
                      key=lambda i: (-abs(next(iter(blocks[i].values()))), next(iter(blocks[i]))))
    off = sorted((i for i, b in enumerate(blocks) if cubic.engine.descriptor(b)[0]),
                 key=lambda i: (-sum(abs(c) for c in blocks[i].values()), cubic.engine.descriptor(blocks[i])[0]))
    current = placement.current_placement(blocks, diagonal, off,
                                          (1 << model['electrons']) - 1, int(model['denominator']))
    current_order = current['order']
    phases['base_and_placement_CPU_seconds'] = clock() - began
    began = clock()
    partition = 'pauli' if method == 'pauli_JD' else 'conserving'
    if partition == 'pauli':
        index = {w: i for i, w in enumerate(sorted(cubic.decode(model['coefficients'])))}
        initial = [index[w] for i in current_order for w in sorted(blocks[i])]
    else:
        initial = list(current_order)
    starts = [initial, list(reversed(initial)), list(range(len(initial)))]
    phases['initializer_CPU_seconds'] = clock() - began
    warm = None
    if method == 'cubic_warm_JD' and clock() + .05 < deadline:
        # Warm-start work is inside the same caller timer and consumes 45% of
        # the remaining nominal ordering allowance, including its cleanup.
        warm_deadline = clock() + .45 * (deadline - clock())
        warm = _search(model, 'conserving', starts, warm_deadline, phases, 'warm', clock)
        if warm['best'] is not None:
            order = warm['best']['order']
            starts = [order, list(reversed(order)), list(current_order)]
    objective = 'cubic' if method.startswith('cubic') else partition
    if clock() + .025 < deadline:
        result = _search(model, objective, starts, deadline, phases, 'objective', clock)
        chosen = result['best']['order'] if result['best'] is not None else initial
    else:
        result = dict(best=None, seeds=[], search=None, status='setup exhausted nominal ordering allowance')
        chosen = warm['best']['order'] if warm and warm['best'] is not None else initial
    # Circuit layout conversion and cleanup belong to the ordering constructor.
    began = clock()
    if partition == 'pauli':
        layout = [{w: c} for w, c in sorted(cubic.decode(model['coefficients']).items())]
    else:
        layout = blocks
    record = dict(arm=method, formula=2, T=TIME, partition=partition, order=list(chosen),
                  unit_lengths=[len(layout[i]) for i in chosen],
                  schedule=[[[[q, p] for q, p in w], str(c)]
                            for i in chosen for w, c in sorted(layout[i].items())])
    del blocks, layout, model, current
    phases['layout_and_cleanup_CPU_seconds'] = clock() - began
    return dict(record=record, phases=phases, result=result, warm_start=warm,
                seed_policy='current, reverse current, canonical; cubic warm uses charged first-order result')


def construct_order(model_path, model_digest, method, total_budget, *, clock=time.process_time):
    require(method in METHODS and total_budget in (2., 8.), 'Undeclared construction cell')
    began = clock()
    wall = time.monotonic()
    result = _ordering(model_path, model_digest, method, began + total_budget / 2, clock)
    elapsed = clock() - began
    result.update(ordering_CPU_seconds=elapsed, ordering_wall_seconds=time.monotonic() - wall,
                  nominal_ordering_CPU_seconds=total_budget / 2,
                  nominal_ordering_overrun_seconds=max(0., elapsed - total_budget / 2),
                  accounting_overhead_CPU_seconds=elapsed - sum(result['phases'].values()),
                  timing='caller through full ordering return and context cleanup')
    return result


def fused_units(lengths, steps):
    require(type(steps) is int and steps in GRID and lengths and all(type(k) is int and k > 0 for k in lengths),
            'Declared positive fine-grid circuit required')
    last, duration = None, Fraction(0)
    m = len(lengths)
    for _ in range(steps):
        sequence = ([(i, Fraction(1, 2 * steps)) for i in range(m - 1)] +
                    [(m - 1, Fraction(1, steps))] +
                    [(i, Fraction(1, 2 * steps)) for i in reversed(range(m - 1))])
        for unit, tau in sequence:
            if unit == last:
                duration += tau
            else:
                if last is not None:
                    yield last, duration
                last, duration = unit, tau
    yield last, duration


def validate_record(model, record):
    cubic.validate_model(model)
    require(record['formula'] == 2 and record['T'] == TIME, 'Formula/time changed')
    require(record['partition'] in ('conserving', 'pauli'), 'Unknown complete-unit partition')
    blocks = [cubic.decode(b) for b in model['blocks']]
    if record['partition'] == 'pauli':
        blocks = [{w: c} for w, c in sorted(cubic.decode(model['coefficients']).items())]
    order = record['order']
    cubic.permutation(order, len(blocks))
    expected = [[[[q, p] for q, p in w], str(c)] for i in order for w, c in sorted(blocks[i].items())]
    require(record['schedule'] == expected and record['unit_lengths'] == [len(blocks[i]) for i in order],
            'Frozen layout does not match complete ordered target units')
    # Reuse the independent parity-of-Pauli-support commutation guard before fusion.
    from common import compiled
    compiled.validate_model(model)


def eligibility(ordering_cpu, compilation_cpu, total_budget, complete):
    require(ordering_cpu >= 0 and compilation_cpu >= 0, 'Negative construction charge')
    charge = ordering_cpu + compilation_cpu
    return dict(total_construction_CPU_seconds=charge,
                total_CPU_budget=total_budget,
                CPU_overrun_seconds=max(0., charge - total_budget) if total_budget is not None else 0.,
                within_total_CPU_budget=charge <= total_budget if total_budget is not None else None,
                eligible=bool(complete and (total_budget is None or charge <= total_budget)))


def _compile(model_path, model_digest, record, steps, trace_path, identity, remaining, clock):
    began = clock()
    require(sha(model_path) == model_digest, 'Compilation input digest differs')
    model = read(model_path)
    validate_record(model, record)
    machine = compiler.Compiler(model['n'])
    rotations = 0
    status, reason = 'complete', 'full exact-policy circuit'
    # Initialization, validation, lowering, counting and trace writing are all
    # charged. Deadline interruption checks bound lowering; later overrun is
    # retained and makes a completed circuit ineligible.
    deadline = began + remaining if remaining is not None else float('inf')
    try:
        if clock() >= deadline:
            raise compiler.CompilationLimit('Total CPU remainder exhausted before HF preparation')
        machine.prepare_hf(model['electrons'])
        offsets = [0]
        for size in record['unit_lengths']:
            offsets.append(offsets[-1] + size)
        duration = Fraction(*TIME.as_integer_ratio())
        denominator = int(model['denominator'])
        for unit, tau in fused_units(record['unit_lengths'], steps):
            for at in range(offsets[unit], offsets[unit + 1]):
                if clock() >= deadline:
                    raise compiler.CompilationLimit('Total CPU remainder exhausted during lowering')
                word, c = record['schedule'][at]
                machine.rotation(word, Fraction(int(c), denominator) * duration * tau)
                rotations += 1
    except compiler.CompilationLimit as error:
        status, reason = 'censored', str(error)
    lowering = clock() - began
    began = clock()
    resources = compiler.resources(model['n'], machine.gates)
    counting = clock() - began
    trace_identity = dict(identity, compilation_status=status)
    measurements = {}
    digest = compiler.write_trace(trace_path, model['n'], machine.gates, trace_identity,
                                  measurements=measurements)
    result = dict(status=status, reason=reason, identity=trace_identity,
                  compiled_resources=resources if status == 'complete' else None,
                  partial_compiled_resources=resources if status != 'complete' else None,
                  generated=machine.generated, completed_ideal_rotations=rotations,
                  lowering_CPU_seconds=lowering, counting_CPU_seconds=counting,
                  trace_sha256=digest, trace_bytes=trace_path.stat().st_size, **measurements)
    del machine, model
    return result


def construct_candidate(model_path, model_digest, record, steps, trace_path, identity,
                        ordering_CPU_seconds, total_budget, *, clock=time.process_time):
    began = clock()
    wall = time.monotonic()
    remaining = None if total_budget is None else max(0., total_budget - ordering_CPU_seconds)
    result = _compile(model_path, model_digest, record, steps, trace_path, identity, remaining, clock)
    elapsed = clock() - began
    result.update(compilation_CPU_seconds=elapsed, compilation_wall_seconds=time.monotonic() - wall,
                  ordering_CPU_charge=ordering_CPU_seconds,
                  **eligibility(ordering_CPU_seconds, elapsed, total_budget, result['status'] == 'complete'),
                  timing='caller through validation/lowering/counting/serialization/hash/cleanup return')
    return result
