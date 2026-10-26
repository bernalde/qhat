"""Frozen transfer-panel identities and exact sharing; no search or propagation.

Logical observations keep every declared arm. Actual work is shared only when
model, training state/method and evaluation/formula/time/cap identities agree.
"""
import itertools
import math
import re
import transfer_common as k

METHODS = ('first_JD', 'pauli_JD', 'cubic_JD')
COMPARISON_TOLERANCE = 2e-12


def arm_spec(arm, evaluation_state):
    k.require(arm in k.ARMS and evaluation_state in k.STATES, 'Undeclared arm/state')
    method = {
        'hf_current': 'current', 'hf_first_JD': 'first_JD',
        'hf_pauli_JD': 'pauli_JD', 'hf_cubic_JD': 'cubic_JD',
        'state_first_JD': 'first_JD', 'state_first_symmetric': 'first_JD',
        'state_pauli_JD': 'pauli_JD', 'state_cubic_JD': 'cubic_JD',
    }[arm]
    return dict(training_state='HF' if arm.startswith('hf_') else evaluation_state,
                method=method, partition='pauli' if method == 'pauli_JD' else 'conserving',
                formula=1 if arm in ('hf_current', 'hf_first_JD', 'state_first_JD') else 2)


def constructor_id(case, state, method):
    k.require(re.fullmatch('[A-Za-z0-9_]+', case) is not None and state in k.STATES and
              method in ('current',) + METHODS, 'Unsafe or undeclared constructor identity')
    k.require(method != 'current' or state == 'HF', 'Legacy current placement is HF-defined')
    return case + '.' + state + '.' + method


def reference_charge_states(training_state, evaluation_state):
    """Standalone reference charges, with the same archived reference used once."""
    k.require(training_state in k.STATES and evaluation_state in k.STATES, 'Unknown reference')
    return [state for state in k.STATES if state in (training_state, evaluation_state)]


def build_plan(cases, *, preflight=False):
    k.require(isinstance(cases, list) and 1 <= len(cases) <= 6 and
              all(isinstance(case, str) for case in cases) and len(set(cases)) == len(cases),
              'Bounded unique case inventory required')
    k.require(all(isinstance(case, str) and re.fullmatch('[A-Za-z0-9_]+', case) for case in cases),
              'Unsafe case identifier')
    k.require(type(preflight) is bool and (not preflight or len(cases) == 1), 'One-case full-path preflight')
    times = (max(k.TIMES),) if preflight else k.TIMES
    training, circuits, logical, seen = [], [], [], set()
    for ci, case in enumerate(cases):
        training.append(dict(constructor=constructor_id(case, 'HF', 'current'), case=case,
                             state='HF', method='current', training_time=1., training_rotation_cap_units=None))
        for si, state in enumerate(k.STATES):
            offset = (ci + si) % len(METHODS)
            for method in METHODS[offset:] + METHODS[:offset]:
                training.append(dict(constructor=constructor_id(case, state, method), case=case,
                                     state=state, method=method, training_time=k.TRAIN_TIME,
                                     training_rotation_cap_units=30))
            for ti, T in enumerate(times):
                for cap_index, cap in enumerate(k.CAPS):
                    offset = (ci + si + ti + cap_index) % len(k.ARMS)
                    for arm in k.ARMS[offset:] + k.ARMS[:offset]:
                        specification = arm_spec(arm, state)
                        trained = constructor_id(case, specification['training_state'], specification['method'])
                        candidate = (trained + '.eval.' + state + '.f' + str(specification['formula']) +
                                     '.t' + str(ti) + '.cap' + str(cap))
                        identity = dict(candidate=candidate, case=case, state=state, time_index=ti, T=T,
                                        cap_multiple=cap, constructor=trained, **specification)
                        if candidate not in seen:
                            circuits.append(identity)
                            seen.add(candidate)
                        logical.append(dict(case=case, state=state, time_index=ti, T=T, cap_multiple=cap,
                                            arm=arm, candidate=candidate, constructor=trained,
                                            reference_charge_states=reference_charge_states(
                                                specification['training_state'], state),
                                            compiled_resource_scope='HF preparation plus evolution' if state == 'HF'
                                                else 'evolution only; quantum preparation cost unknown'))
    constructor_names = {row['constructor'] for row in training}
    k.require(len(training) == len(constructor_names) == 16 * len(cases), 'Constructor sharing differs')
    k.require(all(row['constructor'] in constructor_names for row in logical), 'Missing constructor')
    k.require(len(logical) == len(cases) * len(k.STATES) * len(times) * len(k.CAPS) * len(k.ARMS),
              'Logical panel incomplete')
    # HF roles share three identical formula arms. The first-order-trained
    # symmetric arm still needs its own circuit, making five distinct HF arms.
    k.require(len(circuits) == len(cases) * (5 + 4 * 8) * len(times) * len(k.CAPS),
              'Exact circuit sharing differs')
    pairs = [dict(a=a, b=b) for a, b in itertools.combinations(k.ARMS, 2)]
    return dict(schema=1, preflight=preflight, cases=cases, states=list(k.STATES),
                training_time=k.TRAIN_TIME, legacy_current_time=1., validation_times=list(times),
                cap_multiples=list(k.CAPS), arms=list(k.ARMS), training=training, circuits=circuits,
                logical=logical, pairs=pairs, expected_constructors=len(training),
                expected_unique_circuits=len(circuits), expected_logical_circuits=len(logical),
                expected_reference_evolutions=len(cases) * len(k.STATES) * len(times),
                expected_comparisons=len(cases) * len(k.STATES) * len(times) * len(k.CAPS) * len(pairs),
                error_floor=k.FLOOR, comparison_absolute_tolerance=COMPARISON_TOLERANCE,
                sharing='Only exact model/training state/method and evaluation state/formula/time/cap identity',
                dispatch='State/case-rotated method training;state/case/time/cap-rotated circuit arms',
                no_validation_errors_used_for_construction=True)


def compare(a, b):
    """Retain unavailable and numerical-floor pairs, including their raw errors."""
    k.require(a['status'] in ('complete', 'censored') and b['status'] in ('complete', 'censored'),
              'Unknown circuit status')
    if a['status'] != 'complete' or b['status'] != 'complete':
        return dict(decision='unavailable', informative=False,
                    a_status=a['status'], b_status=b['status'],
                    a_error=a.get('error'), b_error=b.get('error'))
    x, y = a['error'], b['error']
    k.require(all(type(v) in (float, int) and math.isfinite(v) and 0 <= v <= 1 + 1e-12 for v in (x, y)),
              'Invalid comparison errors')
    both_floor = max(x, y) <= k.FLOOR
    decision = ('floor_tie' if both_floor else 'numerical_tie' if abs(x - y) <= COMPARISON_TOLERANCE
                else 'a' if x < y else 'b')
    return dict(decision=decision, informative=not both_floor, a_status='complete', b_status='complete',
                a_error=x, b_error=y, a_minus_b=x - y, both_at_floor=both_floor)
