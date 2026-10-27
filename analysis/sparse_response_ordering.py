"""Bounded ordering from a sparse symmetric-step cubic response.

Qubit zero is the least significant occupation bit. Input units are sums of
commuting Hermitian Pauli words; coefficients are real, in the same energy units
as inverse evolution time (hbar=1). No array indexed by the Hilbert space, exact
time propagation, or validation error is used during training.

The response is a diagonal-reference approximation, not a finite-error bound.
``conservative_steps`` separately bounds the ideal product formula by a Taylor
tail and unitary telescoping. It makes no claim about SCF, synthesis or hardware.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from fractions import Fraction
import math
import time as _time
from typing import Mapping, Sequence

from qhat.analysis.strang_schedule import plan_strang_schedule

Word = tuple[tuple[int, str], ...]
Vector = dict[int, complex]


class ResourceLimit(RuntimeError):
    """The declared reachable-support or cache resource domain was exceeded."""


class SearchDeadline(RuntimeError):
    """A candidate was interrupted; it must not be reported as a full score."""


def _check(deadline):
    if _time.process_time() >= deadline:
        raise SearchDeadline('ordering CPU deadline reached')


def _integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f'{name} is outside its integer domain')
    return value


def _finite(value, name, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a finite real number')
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real number') from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f'{name} must be finite' + (' and positive' if positive else ''))
    return result


def _compile(blocks, n_qubits, deadline=float('inf')):
    _integer(n_qubits, 'n_qubits', 1, 4096)
    if len(blocks) > 4096:
        raise ResourceLimit('too many units')
    result = []
    total = 0
    for block in blocks:
        _check(deadline)
        if not isinstance(block, Mapping):
            raise ValueError('every unit must be a Pauli-word mapping')
        total += len(block)
        if total > 50000:
            raise ResourceLimit('too many Pauli terms')
        words = []
        for word, coefficient in block.items():
            _check(deadline)
            if isinstance(coefficient, complex):
                if coefficient.imag != 0:
                    raise ValueError('Pauli coefficients must be real for Hermiticity')
                coefficient = coefficient.real
            c = _finite(coefficient, 'coefficient')
            flip = parity = ys = 0
            seen = set()
            for q, letter in word:
                _integer(q, 'Pauli index', 0, n_qubits - 1)
                if q in seen or letter not in ('X', 'Y', 'Z'):
                    raise ValueError('Pauli word has a repeated index or invalid letter')
                seen.add(q)
                if letter in ('X', 'Y'):
                    flip |= 1 << q
                if letter in ('Z', 'Y'):
                    parity |= 1 << q
                ys += letter == 'Y'
            if not c and coefficient != 0:
                raise ValueError('nonzero coefficient underflows response arithmetic')
            if c:
                words.append((flip, parity, (1., 1j, -1., -1j)[ys % 4] * c, Fraction(coefficient)))
        # A unit is exponentiated as a whole, or as commuting Pauli rotations.
        # Rejecting a noncommuting block protects both the formula and its cost.
        for i, (x, z, _, _) in enumerate(words):
            _check(deadline)
            for y, w, _, _ in words[:i]:
                if ((x & w).bit_count() + (z & y).bit_count()) % 2:
                    raise ValueError('terms within a unit must commute')
        result.append(tuple(words))
    return tuple(result)


def _reference(reference, n_qubits, maximum):
    if type(reference) is int:
        reference = {reference: 1.}
    if not isinstance(reference, Mapping) or not 0 < len(reference) <= 16:
        raise ValueError('reference must have 1..16 explicitly supplied components')
    result = {}
    for key, value in reference.items():
        if type(key) is not int or key < 0 or key.bit_length() > n_qubits:
            raise ValueError('reference occupation is outside the register')
        value = complex(value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError('reference amplitudes must be finite')
        if value:
            result[key] = value
    norm = sum(x.real * x.real + x.imag * x.imag for x in result.values())
    if not math.isfinite(norm) or abs(norm - 1.) > 1e-12:
        raise ValueError('reference must be normalized')
    if len(result) > maximum:
        raise ResourceLimit('reference exceeds max_entries')
    return result


def molecular_units(terms, n_qubits):
    """Group a real molecular Pauli expansion by flip mask, keeping Z terms.

    ``terms`` may be an OpenFermion ``QubitOperator.terms`` mapping. The
    resulting deterministic units retain every supplied coefficient. Molecular
    equal-flip terms must commute; a mixed noncommuting group is rejected rather
    than being silently exponentiated as a product. No optional backend import
    is needed. The identity, if supplied, is retained as a singleton unit.
    """
    _integer(n_qubits, 'n_qubits', 1, 4096)
    if not isinstance(terms, Mapping):
        raise ValueError('terms must be a Pauli-word mapping')
    if len(terms) > 50000:
        raise ResourceLimit('too many Pauli terms')
    diagonal, groups = [], {}
    for word, coefficient in terms.items():
        # Validate each word independently before using its bit mask. Cross-unit
        # noncommutation is expected and must not be rejected at this stage.
        unit = _compile([{word: coefficient}], n_qubits)[0]
        if not unit:
            continue
        mask = unit[0][0]
        normalized = tuple(sorted(word))
        if not mask:
            diagonal.append({normalized: coefficient})
        else:
            group = groups.setdefault(mask, {})
            group[normalized] = group.get(normalized, 0) + coefficient
    diagonal.sort(key=lambda unit: next(iter(unit)))
    output = diagonal + [dict(sorted(groups[k].items())) for k in sorted(groups)]
    _compile(output, n_qubits)
    return tuple(output)


def _permutation(order, size):
    order = tuple(order)
    if len(order) != size or any(type(i) is not int for i in order) or set(order) != set(range(size)):
        raise ValueError('order must be a complete unit permutation')
    return order


def _stages(order):
    if not order:
        return ()
    return tuple((i, .5) for i in order[:-1]) + ((order[-1], 1.),) + tuple((i, .5) for i in reversed(order[:-1]))


def _add(target, source, factor, maximum, stats=None):
    for key, value in source.items():
        _accumulate(target, key, factor * value, maximum, stats)


def _accumulate(target, key, value, maximum, stats=None):
    value = target.get(key, 0j) + value
    if value:
        if key not in target and len(target) >= maximum:
            raise ResourceLimit('reachable vector exceeds max_entries')
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError('nonfinite response arithmetic')
        target[key] = value
        if stats is not None:
            stats['peak_vector_entries'] = max(stats['peak_vector_entries'], len(target))
    else:
        target.pop(key, None)


@dataclass(frozen=True)
class Score:
    value: float
    cubic_response: float
    steps: int


class SparseResponse:
    """Sparse third-degree columns and coherent diagonal-reference transport.

    Storage follows reached configurations, capped explicitly. A bounded LRU
    caches a unit's action on one occupation bit string, not its full matrix.
    The default at most 16 reference components covers HF and sparse coherent
    references without implying efficient quantum state preparation.
    """
    def __init__(self, blocks, n_qubits, reference, *, time=1., max_entries=50000,
                 max_cached_actions=50000, deadline=float('inf')):
        self.time = _finite(time, 'time', positive=True)
        self.max_entries = _integer(max_entries, 'max_entries', 1, 1000000)
        self.max_cached_actions = _integer(max_cached_actions, 'max_cached_actions', 0, 1000000)
        self.units = _compile(blocks, n_qubits, deadline)
        self.n_qubits = n_qubits
        self.reference = _reference(reference, n_qubits, self.max_entries)
        self.lengths = tuple(len(u) for u in self.units)
        self.active_units = tuple(i for i, unit in enumerate(self.units) if unit)
        self.diagonal_words = tuple((parity, phase) for unit in self.units
                                    for flip, parity, phase, _ in unit if not flip)
        self.cache = OrderedDict()
        self.cache_entries = 0
        self.stats = dict(action_evaluations=0, cache_hits=0, peak_vector_entries=len(self.reference),
                          peak_cached_actions=0, reference_column_entries=0)
        self.reference_cubic = {}
        self.last_order = None
        self.last_columns = None

    def _reference_column(self, source, deadline):
        if source not in self.reference_cubic:
            column = {source: 1. + 0j}
            for _ in range(3):
                column = self.hamiltonian_action(column, deadline)
            column = {k: v / 6 for k, v in column.items()}
            if self.stats['reference_column_entries'] + len(column) > self.max_entries:
                raise ResourceLimit('cached reference columns exceed max_entries')
            self.reference_cubic[source] = column
            self.stats['reference_column_entries'] += len(column)
        return self.reference_cubic[source]

    def _unit_column(self, index, source, deadline):
        _check(deadline)
        key = index, source
        if key in self.cache:
            self.stats['cache_hits'] += 1
            value = self.cache.pop(key)
            self.cache[key] = value
            return value
        output = {}
        for flip, parity, phase, _ in self.units[index]:
            destination = source ^ flip
            value = phase * (-1 if (source & parity).bit_count() % 2 else 1)
            _accumulate(output, destination, value, self.max_entries, self.stats)
        result = tuple(output.items())
        self.stats['action_evaluations'] += 1
        cost = max(1, len(result))
        if cost <= self.max_cached_actions:
            while self.cache and self.cache_entries + cost > self.max_cached_actions:
                _, removed = self.cache.popitem(last=False)
                self.cache_entries -= max(1, len(removed))
            self.cache[key] = result
            self.cache_entries += cost
            self.stats['peak_cached_actions'] = max(self.stats['peak_cached_actions'], self.cache_entries)
        return result

    def apply(self, unit, vector, deadline=float('inf')):
        output = {}
        for source, amplitude in vector.items():
            _check(deadline)
            for destination, value in self._unit_column(unit, source, deadline):
                _accumulate(output, destination, value * amplitude, self.max_entries, self.stats)
        self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(output))
        return output

    def hamiltonian_action(self, vector, deadline=float('inf')):
        output = {}
        for unit in range(len(self.units)):
            _add(output, self.apply(unit, vector, deadline), 1., self.max_entries, self.stats)
        self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(output))
        return output

    def diagonal_energy(self, source, deadline=float('inf')):
        result = 0j
        for parity, phase in self.diagonal_words:
            _check(deadline)
            result += phase * (-1 if (source & parity).bit_count() % 2 else 1)
        if abs(result.imag) > 1e-12 or not math.isfinite(result.real):
            raise ValueError('invalid real diagonal energy')
        return result.real

    def cubic_column(self, order, source, deadline=float('inf')):
        """Coefficient of h^3 in S2(h)-exp(hH); common physical i is stripped."""
        order = _permutation(order, len(self.units))
        if source not in self.reference:
            raise ValueError('source must belong to the supplied reference support')
        _check(deadline)
        if len(self.active_units) <= 1:
            return {}
        reference = self._reference_column(source, deadline)
        q = [{source: 1. + 0j}, {}, {}, {}]
        for unit, tau in _stages(order):
            _check(deadline)
            out = [dict(a) for a in q]
            powers = list(q)
            factor = 1.
            for degree in range(1, 4):
                factor *= tau / degree
                for k in range(4 - degree):
                    powers[k] = self.apply(unit, powers[k], deadline)
                    _add(out[k + degree], powers[k], factor, self.max_entries, self.stats)
            q = out
            self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], *(len(a) for a in q))
        _add(q[3], reference, -1., self.max_entries, self.stats)
        self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(q[3]))
        return q[3]

    def affordable_steps(self, order, rotation_budget):
        order = _permutation(order, len(self.units))
        _integer(rotation_budget, 'rotation_budget')
        lengths = tuple(self.lengths[i] for i in order if self.lengths[i])
        if not lengths:
            return 1
        plan = plan_strang_schedule(lengths, rotation_budget)
        if not plan.feasible:
            raise ValueError('rotation budget cannot afford one symmetric step')
        return plan.steps

    def columns(self, order, deadline=float('inf')):
        order = _permutation(order, len(self.units))
        if order == self.last_order and self.last_columns is not None:
            _check(deadline)
            return {source: dict(column) for source, column in self.last_columns.items()}
        columns = {source: self.cubic_column(order, source, deadline)
                   for source in self.reference}
        if sum(map(len, columns.values())) > self.max_entries:
            raise ResourceLimit('combined cubic columns exceed max_entries')
        return columns

    def adjacent_columns(self, order, columns, position, deadline=float('inf')):
        """Swap adjacent A,B: delta C=[(A+B)/8+suffix/4,[A,B]].

        This is an exact algebraic update of the cubic coefficient, evaluated
        numerically on each reference column before diagonal transport. Accepted
        search moves are checked again through the full Taylor recurrence.
        """
        order = _permutation(order, len(self.units))
        _integer(position, 'adjacent position', 0, len(order) - 2)
        a, b = order[position:position + 2]
        suffix = order[position + 2:]
        def commutator(vector):
            left = self.apply(a, self.apply(b, vector, deadline), deadline)
            _add(left, self.apply(b, self.apply(a, vector, deadline), deadline), -1., self.max_entries, self.stats)
            return left
        def linear(vector):
            out = {}
            for unit, factor in ((a, .125), (b, .125)):
                _add(out, self.apply(unit, vector, deadline), factor, self.max_entries, self.stats)
            for unit in suffix:
                _add(out, self.apply(unit, vector, deadline), .25, self.max_entries, self.stats)
            return out
        result = {}
        for source in self.reference:
            initial = {source: 1. + 0j}
            delta = linear(commutator(initial))
            _add(delta, commutator(linear(initial)), -1., self.max_entries, self.stats)
            column = dict(columns[source])
            _add(column, delta, 1., self.max_entries, self.stats)
            result[source] = column
            self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(delta), len(column))
        if sum(map(len, result.values())) > self.max_entries:
            raise ResourceLimit('updated cubic columns exceed max_entries')
        new_order = list(order)
        new_order[position:position + 2] = [b, a]
        return tuple(new_order), result

    def score_columns(self, order, rotation_budget, columns, deadline=float('inf')):
        _check(deadline)
        order = _permutation(order, len(self.units))
        steps = self.affordable_steps(order, rotation_budget)
        if set(columns) != set(self.reference):
            raise ValueError('cubic columns must cover the reference support')
        action = {}
        energies = {}
        for source, amplitude in self.reference.items():
            if source not in energies:
                energies[source] = self.diagonal_energy(source, deadline)
            es = energies[source]
            for destination, coefficient in columns[source].items():
                _check(deadline)
                if destination not in energies:
                    if len(energies) >= self.max_entries:
                        raise ResourceLimit('energy cache exceeds max_entries')
                    energies[destination] = self.diagonal_energy(destination, deadline)
                x = self.time * (energies[destination] - es) / 2
                if not math.isfinite(x):
                    raise ValueError('nonfinite diagonal response phase')
                kernel = self.time * (math.sin(x) / x if x else 1.) * complex(math.cos(x), math.sin(x))
                _add(action, {destination: coefficient * amplitude * kernel}, 1., self.max_entries, self.stats)
                self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(action))
        projection = sum(a.conjugate() * action.get(k, 0j) for k, a in self.reference.items())
        _add(action, self.reference, -projection, self.max_entries, self.stats)
        self.stats['peak_vector_entries'] = max(self.stats['peak_vector_entries'], len(action))
        response = sum(x.real * x.real + x.imag * x.imag for x in action.values())
        value = response / steps ** 4
        if not math.isfinite(value) or not math.isfinite(response):
            raise ValueError('nonfinite cubic response score')
        return Score(value=value, cubic_response=response, steps=steps)

    def score(self, order, rotation_budget, deadline=float('inf')):
        _check(deadline)
        order = _permutation(order, len(self.units))
        self.affordable_steps(order, rotation_budget)
        columns = self.columns(order, deadline)
        result = self.score_columns(order, rotation_budget, columns, deadline)
        self.last_order, self.last_columns = order, columns
        return result


@dataclass(frozen=True)
class TrainingResult:
    order: tuple[int, ...]
    score: Score | None
    CPU_seconds: float
    evaluations: int
    accepted_moves: int
    status: str
    stats: dict


def _incremental_search(context, scored, budget, deadline, max_evaluations, max_moves,
                        best_order, best_score, evaluations, reserve):
    """Visit insertion destinations by adjacent updates; retain rescored gains."""
    moves = 0
    status = 'complete'
    start = _time.process_time()
    slice_seconds = max(0., deadline - start) / max(1, len(scored))
    context.stats['incremental_candidates'] = 0
    context.stats['full_move_rescores'] = 0
    for seed_index, (order, score) in enumerate(scored):
        end = min(deadline, start + (seed_index + 1) * slice_seconds)
        current, current_score = order, score
        for _ in range(max_moves):
            proposal, proposal_score = current, current_score
            try:
                columns = context.columns(current, end - reserve)
                for source in range(len(current)):
                    for direction in (-1, 1):
                        trial, trial_columns = current, columns
                        at = source
                        while 0 <= at + direction < len(current):
                            _check(end - reserve)
                            if evaluations >= max_evaluations - 1:
                                raise SearchDeadline('evaluation limit reached')
                            position = at - 1 if direction < 0 else at
                            trial, trial_columns = context.adjacent_columns(trial, trial_columns, position, end - reserve)
                            try:
                                context.affordable_steps(trial, budget)
                            except ValueError:
                                at += direction
                                continue
                            trial_score = context.score_columns(trial, budget, trial_columns, end - reserve)
                            evaluations += 1
                            context.stats['incremental_candidates'] += 1
                            if (trial_score.value, trial) < (proposal_score.value, proposal):
                                proposal, proposal_score = trial, trial_score
                            at += direction
            except SearchDeadline:
                status = 'evaluation_limit' if evaluations >= max_evaluations - 1 else 'cpu_limit'
            except ResourceLimit:
                status = 'support_limit'
            threshold = max(1e-30, abs(current_score.value) * 1e-12)
            if proposal_score.value >= current_score.value - threshold:
                break
            try:
                if evaluations >= max_evaluations:
                    status = 'evaluation_limit'
                    break
                actual = context.score(proposal, budget, end)
                evaluations += 1
                context.stats['full_move_rescores'] += 1
            except (SearchDeadline, ResourceLimit):
                status = 'cpu_limit' if _time.process_time() >= end else 'support_limit'
                break
            if actual.value >= current_score.value - threshold:
                break
            current, current_score = proposal, actual
            moves += 1
            if (actual.value, proposal) < (best_score.value, best_order):
                best_order, best_score = proposal, actual
            if status != 'complete':
                break
        if _time.process_time() >= deadline or evaluations >= max_evaluations:
            break
    return best_order, best_score, evaluations, moves, status


def train_ordering(blocks, n_qubits, reference, *, time=1., rotation_budget,
                   cpu_seconds=2., max_evaluations=3000, max_moves=5,
                   max_entries=50000, max_cached_actions=50000, starts=None):
    """Return the best fully scored candidate under a cooperative CPU allowance.

    If context/seed work exhausts the cap, a valid initializer is returned with
    score=None. Deadline/support failures never become a zero or partial score.
    A score is a response heuristic, not an actual finite-state-error certificate.
    Callers needing strict total construction accounting must time the complete
    call, including input loading and object cleanup, as the benchmark does.
    """
    began = _time.process_time()
    allowance = _finite(cpu_seconds, 'cpu_seconds', positive=True)
    _integer(max_evaluations, 'max_evaluations', 0, 1000000)
    _integer(max_moves, 'max_moves', 0, 1000)
    _integer(rotation_budget, 'rotation_budget')
    deadline = began + allowance
    # Materialize and validate once; even an interrupted search returns a valid
    # permutation of the caller's units. Validation cannot be bypassed by timeout.
    _integer(max_entries, 'max_entries', 1, 1000000)
    _integer(max_cached_actions, 'max_cached_actions', 0, 1000000)
    _finite(time, 'time', positive=True)
    reference = _reference(reference, _integer(n_qubits, 'n_qubits', 1, 4096), max_entries)
    compiled = _compile(blocks, n_qubits, deadline)
    canonical = tuple(range(len(compiled)))
    magnitude = tuple(sorted(canonical, key=lambda i: (-sum(abs(c) for _, _, c, _ in compiled[i]), i)))
    seeds = [magnitude, tuple(reversed(magnitude)), canonical] if starts is None else list(starts)
    if not seeds:
        raise ValueError('at least one seed is required')
    seeds = list(dict.fromkeys(_permutation(order, len(compiled)) for order in seeds))
    def affordable(order):
        lengths = [len(compiled[i]) for i in order if compiled[i]]
        return not lengths or plan_strang_schedule(lengths, rotation_budget).feasible
    seeds = [order for order in seeds if affordable(order)]
    if not seeds:
        fallback = tuple(sorted(canonical, key=lambda i: (len(compiled[i]), i)))
        if not affordable(fallback):
            raise ValueError('rotation budget cannot afford one symmetric step')
        seeds = [fallback]
    best_order, best_score = seeds[0], None
    evaluations = moves = 0
    context = None
    status = 'complete'
    try:
        context = SparseResponse(blocks, n_qubits, reference, time=time, max_entries=max_entries,
                                 max_cached_actions=max_cached_actions, deadline=deadline)
        scored = []
        longest_seed = 0.
        for order in seeds:
            if evaluations >= max_evaluations:
                status = 'evaluation_limit'
                break
            seed_began = _time.process_time()
            score = context.score(order, rotation_budget, deadline)
            longest_seed = max(longest_seed, _time.process_time() - seed_began)
            evaluations += 1
            scored.append((order, score))
            if best_score is None or (score.value, order) < (best_score.value, best_order):
                best_order, best_score = order, score
        if len(compiled) > 8 and scored:
            best_order, best_score, evaluations, moves, status = _incremental_search(
                context, scored, rotation_budget, deadline, max_evaluations, max_moves,
                best_order, best_score, evaluations, max(.02, longest_seed * 1.25))
            stats = dict(context.stats)
            stats['response_is_error_bound'] = False
            del context
            return TrainingResult(best_order, best_score, _time.process_time() - began,
                                  evaluations, moves, status, stats)
        for order, score in scored:
            current, current_score = order, score
            for _ in range(max_moves):
                candidate, candidate_score = current, current_score
                for source in range(len(current)):
                    for destination in range(len(current)):
                        if source == destination:
                            continue
                        if evaluations >= max_evaluations:
                            raise SearchDeadline('evaluation limit reached')
                        trial = list(current)
                        trial.insert(destination, trial.pop(source))
                        trial = tuple(trial)
                        _check(deadline)
                        if not affordable(trial):
                            continue
                        trial_score = context.score(trial, rotation_budget, deadline)
                        evaluations += 1
                        if (trial_score.value, trial) < (candidate_score.value, candidate):
                            candidate, candidate_score = trial, trial_score
                        # A completed candidate remains usable if a later sweep
                        # is interrupted, without inventing an accepted move.
                        if (trial_score.value, trial) < (best_score.value, best_order):
                            best_order, best_score = trial, trial_score
                threshold = max(1e-30, abs(current_score.value) * 1e-12)
                if candidate_score.value >= current_score.value - threshold:
                    break
                current, current_score = candidate, candidate_score
                moves += 1
    except SearchDeadline:
        status = 'evaluation_limit' if evaluations >= max_evaluations else 'cpu_limit'
    except ResourceLimit:
        status = 'support_limit'
    stats = dict(context.stats) if context is not None else {}
    stats['response_is_error_bound'] = False
    del context
    return TrainingResult(best_order, best_score, _time.process_time() - began,
                          evaluations, moves, status, stats)


@dataclass(frozen=True)
class StepCertificate:
    steps: int
    operator_error_bound: Fraction
    metric: str
    error_bound: Fraction


def conservative_steps(blocks, n_qubits, *, time, error_tolerance,
                       metric='operator', max_steps=10**9):
    """Exact-rational upper bound for an ideal symmetric product formula.

    Let a=|T| sum|c| and r>=2a. Both degree>=3 Taylor tails are bounded
    by x^3/[6(1-x/4)], x=a/r<=1/2. Unitary telescoping therefore gives
    ||S2(T/r)^r-exp(-iTH)|| <= (2/3)*a^3/r^2. The loose 2/3 constant is
    deliberate. For normalized states, 1-|overlap| <= delta^2/2 and
    1-|overlap|^2 <= delta^2. Arithmetic here is rational, not a floating
    interval assertion. Input floats denote their exact supplied dyadics.
    """
    _finite(time, 'time')
    _finite(error_tolerance, 'error_tolerance', positive=True)
    _integer(max_steps, 'max_steps', 1)
    if metric not in ('operator', 'infidelity', 'squared_infidelity'):
        raise ValueError('unknown error metric')
    units = _compile(blocks, n_qubits)
    a = abs(Fraction(time)) * sum((abs(c) for unit in units for _, _, _, c in unit), Fraction(0))
    target = Fraction(error_tolerance)
    if not a or sum(bool(unit) for unit in units) <= 1:
        return StepCertificate(1, Fraction(0), metric, Fraction(0))
    lower = max(1, (2 * a.numerator + a.denominator - 1) // a.denominator)
    def bounds(r):
        delta = min(Fraction(2), Fraction(2, 3) * a**3 / r**2)
        error = delta if metric == 'operator' else delta**2 / (2 if metric == 'infidelity' else 1)
        return delta, error
    if lower > max_steps or bounds(max_steps)[1] > target:
        raise ResourceLimit('conservative target requires more than max_steps')
    lo, hi = lower, max_steps
    while lo < hi:
        middle = (lo + hi) // 2
        if bounds(middle)[1] <= target:
            hi = middle
        else:
            lo = middle + 1
    delta, error = bounds(lo)
    return StepCertificate(lo, delta, metric, error)
