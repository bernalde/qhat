"""Actual twelve-mode line circuits from authenticated frozen layouts."""
from fractions import Fraction
import math
import time
import common as k
import compiler


def rotation_count(lengths, steps):
    k.require(lengths and all(type(x) is int and x > 0 for x in lengths) and
              type(steps) is int and 1 <= steps <= 128, 'Rotation accounting domain')
    if len(lengths) == 1:
        return lengths[0]
    return steps * (2 * sum(lengths) - lengths[0] - lengths[-1]) + lengths[0]


def affordable_steps(lengths, multiple):
    k.require(type(multiple) is int and multiple in k.CAPS, 'Declared rotation budget')
    cap = multiple * sum(lengths)
    if len(lengths) == 1:
        return 1
    steps = (cap - lengths[0]) // (2 * sum(lengths) - lengths[0] - lengths[-1])
    k.require(1 <= steps <= 128 and rotation_count(lengths, steps) <= cap and
              rotation_count(lengths, steps + 1) > cap, 'Maximal affordable symmetric step count')
    return steps


def fused_units(lengths, steps):
    rotation_count(lengths, steps)
    last, duration = None, Fraction(0)
    for _ in range(steps):
        sequence = ([(i, Fraction(1, 2 * steps)) for i in range(len(lengths) - 1)] +
                    [(len(lengths) - 1, Fraction(1, steps))] +
                    [(i, Fraction(1, 2 * steps)) for i in reversed(range(len(lengths) - 1))])
        for unit, tau in sequence:
            if unit == last:
                duration += tau
            else:
                if last is not None:
                    yield last, duration
                last, duration = unit, tau
    yield last, duration


def validate_record(model, record):
    k.blocks(model)
    raw = ([{w: c} for w, c in sorted(k.decode(model['coefficients']).items())]
           if record['partition'] == 'pauli' else [k.decode(b) for b in model['blocks']])
    k.require(record['partition'] in ('conserving', 'pauli') and record['formula'] == 2 and
              sorted(record['order']) == list(range(len(raw))), 'Frozen order/partition')
    expected = [[[[q, p] for q, p in word], str(c)] for i in record['order'] for word, c in sorted(raw[i].items())]
    k.require(record['schedule'] == expected and record['unit_lengths'] == [len(raw[i]) for i in record['order']],
              'Frozen layout differs from exact target')
    for block in raw:
        words = [dict(w) for w in block]
        for at, a in enumerate(words):
            for b in words[at + 1:]:
                k.require(sum(a[q] != b[q] for q in a.keys() & b.keys()) % 2 == 0,
                          'Unit boundary fusion requires internal commutation')


def compile_candidate(model, record, T, multiple, trace, identity):
    cpu, wall = time.process_time(), time.monotonic()
    validate_record(model, record)
    k.require(type(T) in (float, int) and math.isfinite(T) and T > 0, 'Finite evolution time')
    steps = affordable_steps(record['unit_lengths'], multiple)
    machine = compiler.Compiler(model['n'])
    machine.prepare_hf(model['electrons'])
    offsets = [0]
    for count in record['unit_lengths']:
        offsets.append(offsets[-1] + count)
    duration, denominator = Fraction(T), int(model['denominator'])
    status, reason, rotations = 'complete', 'complete declared circuit', 0
    try:
        for unit, tau in fused_units(record['unit_lengths'], steps):
            for at in range(offsets[unit], offsets[unit + 1]):
                word, coefficient = record['schedule'][at]
                machine.rotation(word, Fraction(int(coefficient), denominator) * duration * tau)
                rotations += 1
    except compiler.CompilationLimit as error:
        status, reason = 'censored', str(error)
    resources = compiler.resources(model['n'], machine.gates)
    trace_identity = dict(identity, T=T, steps=steps, cap_multiple=multiple, formula=2,
                          preparation='HF_prefix', compilation_status=status)
    measurements = {}
    digest = compiler.write_trace(trace, model['n'], machine.gates, trace_identity, measurements=measurements)
    expected = rotation_count(record['unit_lengths'], steps)
    k.require(status != 'complete' or expected == rotations, 'Actual fused rotation count')
    result = dict(status=status, reason=reason, identity=trace_identity, steps=steps, cap_multiple=multiple,
        rotation_cap=multiple * sum(record['unit_lengths']), ideal_rotations=rotations,
        expected_ideal_rotations=expected, generated=machine.generated,
        compiled_resources=resources if status == 'complete' else None,
        partial_compiled_resources=resources if status != 'complete' else None,
        trace_sha256=digest, trace_bytes=trace.stat().st_size,
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
        quantum_preparation_circuit_cost='HF prefix included', compiler=compiler.VERSION, **measurements)
    del machine
    return result
