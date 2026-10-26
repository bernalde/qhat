"""Transfer circuits with explicit external-state preparation and complete costs.

Search is separate. This module accepts a frozen complete-unit order and uses
the unchanged reviewed exact-angle line compiler. No validation error influences
the chosen order, formula, rotation cap or affordable step count.
"""
from fractions import Fraction
from functools import lru_cache
import math
from pathlib import Path
import time
import numpy as np
import transfer_common as k
from cubic_response import prior, validate

compiler = k.load('transfer_exact_angle_compiler', k.AUDIT / 'compiled_accuracy/compiler.py')


def validate_record(model, record):
    blocks = validate(model)
    k.require(len(model['coefficients']) <= 128 and len(blocks) <= 128,
              'Inherited compiler/reference term domain exceeded')
    k.require(type(record['formula']) is int and record['formula'] in (1, 2),
              'First or symmetric second-order formula required')
    k.require(record['partition'] in ('conserving', 'pauli'), 'Unknown partition')
    if record['partition'] == 'pauli':
        blocks = [{w: c} for w, c in sorted(prior.decode(model['coefficients']).items())]
    order = record['order']
    k.require(isinstance(order, list) and all(type(i) is int for i in order), 'Integral order required')
    prior.permutation(order, len(blocks))
    schedule = [[[[q, p] for q, p in word], str(c)]
                for i in order for word, c in sorted(blocks[i].items())]
    k.require(record['unit_lengths'] == [len(blocks[i]) for i in order] and
              record['schedule'] == schedule, 'Layout differs from frozen complete target units')
    # Equal-unit boundary fusion requires commutation within each complete unit.
    for block in blocks:
        words = [dict(word) for word in block]
        for at, left in enumerate(words):
            for right in words[at + 1:]:
                k.require(sum(left[q] != right[q] for q in left.keys() & right.keys()) % 2 == 0,
                          'Noncommuting words cannot be fused as one unit')
    return blocks


def rotation_count(lengths, steps, formula):
    k.require(isinstance(lengths, list) and lengths and len(lengths) <= 128 and
              all(type(size) is int and size > 0 for size in lengths) and sum(lengths) <= 128,
              'Bounded complete unit lengths required')
    k.require(type(steps) is int and 1 <= steps <= 128 and type(formula) is int and formula in (1, 2),
              'Bounded formula/step domain required')
    if formula == 1:
        return steps * sum(lengths)
    if len(lengths) == 1:
        return lengths[0]
    return steps * (2 * sum(lengths) - lengths[0] - lengths[-1]) + lengths[0]


def affordable_steps(lengths, cap_multiple, formula):
    k.require(type(cap_multiple) is int and cap_multiple in k.CAPS, 'Declared 10U or 30U cap required')
    rotation_count(lengths, 1, formula)
    cap = cap_multiple * sum(lengths)
    steps = k.engine.feasible_steps(lengths, cap, formula == 2)
    k.require(type(steps) is int and 1 <= steps <= 128 and
              rotation_count(lengths, steps, formula) <= cap, 'No bounded affordable circuit')
    if formula == 1 or len(lengths) > 1:
        k.require(rotation_count(lengths, steps + 1, formula) > cap, 'Step count is not maximal')
    return steps


def fused_units(lengths, steps, formula):
    rotation_count(lengths, steps, formula)
    if formula == 1:
        for _ in range(steps):
            for unit in range(len(lengths)):
                yield unit, Fraction(1, steps)
        return
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


def _compile(model_path, digest, record, T, cap_multiple, state_label, trace_path, identity, clock):
    cpu = clock()
    model_path, trace_path = Path(model_path), Path(trace_path)
    k.require(state_label in k.STATES and type(T) in (float, int) and math.isfinite(T) and T > 0,
              'Declared state and positive finite evolution time required')
    k.require(model_path.stat().st_size <= 4 * 1024**2 and k.sha(model_path) == digest,
              'Model input size/hash differs')
    model = k.read(model_path)
    validate_record(model, record)
    steps = affordable_steps(record['unit_lengths'], cap_multiple, record['formula'])
    preparation = 'HF_prefix' if state_label == 'HF' else 'external_state_not_compiled'
    reserved = {'state', 'T', 'formula', 'steps', 'cap_multiple', 'preparation', 'compilation_status'}
    k.require(isinstance(identity, dict) and not reserved.intersection(identity), 'Reserved trace identity key')
    trace_identity = dict(identity, state=state_label, T=T, formula=record['formula'], steps=steps,
                          cap_multiple=cap_multiple, preparation=preparation)
    phases = {'input_validation_CPU_seconds': clock() - cpu}
    cpu = clock()
    machine = compiler.Compiler(model['n'])
    status, reason, rotations = 'complete', 'full declared circuit', 0
    offsets = np.cumsum([0] + record['unit_lengths'])
    duration, denominator = Fraction(*float(T).as_integer_ratio()), int(model['denominator'])
    try:
        if state_label == 'HF':
            machine.prepare_hf(model['electrons'])
        for unit, tau in fused_units(record['unit_lengths'], steps, record['formula']):
            for at in range(offsets[unit], offsets[unit + 1]):
                word, coefficient = record['schedule'][at]
                machine.rotation(word, Fraction(int(coefficient), denominator) * duration * tau)
                rotations += 1
    except compiler.CompilationLimit as error:
        status, reason = 'censored', str(error)
    phases['lowering_CPU_seconds'] = clock() - cpu
    cpu = clock()
    resources = compiler.resources(model['n'], machine.gates)
    phases['counting_CPU_seconds'] = clock() - cpu
    trace_identity['compilation_status'] = status
    measurements = {}
    trace_digest = compiler.write_trace(trace_path, model['n'], machine.gates, trace_identity,
                                        measurements=measurements)
    expected = rotation_count(record['unit_lengths'], steps, record['formula'])
    k.require(status != 'complete' or rotations == expected, 'Actual fused rotation count differs')
    result = dict(status=status, reason=reason, identity=trace_identity, steps=steps,
                  rotation_cap=cap_multiple * sum(record['unit_lengths']),
                  ideal_rotations=rotations, expected_ideal_rotations=expected,
                  compiled_resources=resources if status == 'complete' else None,
                  partial_compiled_resources=resources if status != 'complete' else None,
                  generated=machine.generated, trace_sha256=trace_digest, trace_bytes=trace_path.stat().st_size,
                  phases=phases, **measurements,
                  quantum_preparation_circuit_cost='HF prefix included' if state_label == 'HF' else 'unknown',
                  compiled_resource_scope='HF preparation plus evolution' if state_label == 'HF' else
                                          'evolution only; correlated/coherent preparation cost unknown')
    del machine, model, offsets
    return result


def construct_candidate(model_path, digest, record, T, cap_multiple, state_label, trace_path, identity,
                        *, clock=time.process_time):
    cpu, wall = clock(), time.monotonic()
    result = _compile(model_path, digest, record, T, cap_multiple, state_label, trace_path, identity, clock)
    result.update(compilation_CPU_seconds=clock() - cpu, compilation_wall_seconds=time.monotonic() - wall,
                  timing='caller through input validation,lowering,counting,serialization,hash and cleanup return')
    return result


def state_vector(model, psi):
    # Register validation precedes exponentiation and any state-sized allocation.
    validate(model)
    k.require(np.finfo(np.longdouble).nmant >= 63, 'Extended-precision replay requires 63 mantissa bits')
    k.require(isinstance(psi, np.ndarray) and psi.shape == (1 << model['n'],),
              'Reference shape before allocation')
    value = np.asarray(psi, dtype=np.clongdouble)
    k.require(np.all(np.isfinite(value)) and
              abs(np.vdot(value, value) - 1) < 1e-10 and 0 < np.count_nonzero(value) <= 9,
              'Declared normalized full-space reference/support required')
    return value.copy()


class Replay:
    """Primitive action on the supplied full vector; no sector projection."""
    def __init__(self, model, psi):
        self.initial = state_vector(model, psi)
        self.model, self.n = model, model['n']
        basis = np.arange(len(self.initial))
        self.low = {q: basis[(basis & (1 << q)) == 0] for q in range(self.n)}
        self.high = {q: self.low[q] ^ (1 << q) for q in range(self.n)}
        self.cx = {(q, t): basis ^ (((basis >> q) & 1) << t)
                   for q in range(self.n) for t in range(self.n) if abs(q - t) == 1}

    @lru_cache(maxsize=128)
    def phase(self, numerator, denominator):
        angle = np.longdouble(numerator) / np.longdouble(denominator) / 2
        k.require(bool(np.isfinite(angle)), 'Nonfinite primitive angle')
        return np.clongdouble(np.cos(angle) - 1j * np.sin(angle))

    def replay(self, path, identity):
        k.require(Path(path).stat().st_size <= 32 * 1024**2, 'Trace byte bound')
        iterator = compiler.read_trace(path)
        header = next(iterator)
        k.require(header['n'] == self.n and header['identity'] == identity and
                  identity['compilation_status'] == 'complete', 'Complete trace identity differs')
        if identity['preparation'] == 'HF_prefix':
            hf = (1 << self.model['electrons']) - 1
            k.require(identity['state'] == 'HF' and abs(self.initial[hf] - 1) < 1e-14 and
                      np.count_nonzero(self.initial) == 1, 'HF prefix/reference differs')
            state = np.zeros_like(self.initial)
            state[0] = 1
        else:
            k.require(identity['preparation'] == 'external_state_not_compiled' and
                      identity['state'] in k.STATES[1:], 'External preparation declaration differs')
            state = self.initial.copy()
        for gate in iterator:
            if gate.kind == 'cx':
                state = state[self.cx[gate.a, gate.b]]
            elif gate.kind == 'h':
                low, high = state[self.low[gate.a]].copy(), state[self.high[gate.a]].copy()
                state[self.low[gate.a]] = (low + high) / np.sqrt(np.longdouble(2))
                state[self.high[gate.a]] = (low - high) / np.sqrt(np.longdouble(2))
            elif gate.kind == 's':
                state[self.high[gate.a]] *= 1j
            elif gate.kind == 'sdg':
                state[self.high[gate.a]] *= -1j
            else:
                phase = self.phase(gate.angle.numerator, gate.angle.denominator)
                state[self.low[gate.a]] *= phase
                state[self.high[gate.a]] *= phase.conjugate()
        return state


def metrics(exact, actual, model):
    norm = float(np.linalg.norm(actual))
    k.require(math.isfinite(norm) and norm > 0, 'Invalid evolved-state norm')
    a, z = actual / norm, exact / np.linalg.norm(exact)
    inner = np.vdot(z, a)
    residual = a - inner * z
    loss = float(np.vdot(residual, residual).real)
    n, hf = model['n'], (1 << model['electrons']) - 1
    sector = [sum((hf >> q) & 1 for q in range(parity, n, 2)) for parity in (0, 1)]
    number = np.array([v.bit_count() != model['electrons'] for v in range(len(a))])
    joint = np.array([[sum((v >> q) & 1 for q in range(parity, n, 2)) for parity in (0, 1)] != sector
                      for v in range(len(a))])
    return dict(error=loss / (1 + float(abs(inner))), infidelity=loss, norm=norm,
                number_leakage=float(np.vdot(a[number], a[number]).real),
                joint_sector_leakage=float(np.vdot(a[joint], a[joint]).real))
