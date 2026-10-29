"""Analytic regressions for scalar invariance and unresolved cubic responses."""
import copy
import inspect
import itertools
import math

import numpy as np
import pytest

from qhat.analysis import sparse_response_ordering as response


_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.diag([1., -1.]).astype(complex)
_I = np.eye(2, dtype=complex)


def _literal_blocks(blocks, n):
    matrices = []
    for block in blocks:
        matrix = np.zeros((2**n, 2**n), dtype=complex)
        for word, coefficient in block.items():
            labels = dict(word)
            term = np.ones((1, 1), dtype=complex)
            for q in reversed(range(n)):
                term = np.kron(term, {'X': _X, 'Y': _Y, 'Z': _Z}.get(labels.get(q), _I))
            matrix += coefficient * term
        matrices.append(matrix)
    return matrices


def _comm(a, b):
    return a @ b - b @ a


def _nested_cubic(matrices, order):
    """Independent commutator identity, not a difference of Taylor columns."""
    cubic = np.zeros_like(matrices[0])
    for at, index in enumerate(order[:-1]):
        a = matrices[index]
        suffix = sum((matrices[j] for j in order[at + 1:]), np.zeros_like(a))
        cubic += -_comm(a, _comm(a, suffix)) / 24 + _comm(suffix, _comm(suffix, a)) / 12
    return cubic


def _transported_response(matrices, cubic, reference, time):
    psi = np.zeros(matrices[0].shape[0], dtype=complex)
    for source, amplitude in reference.items():
        psi[source] = amplitude
    diagonal = np.diag(sum(matrices)).real
    gaps = diagonal[:, None] - diagonal[None, :]
    kernel = time * np.sinc(time * gaps / (2*np.pi)) * np.exp(.5j * time * gaps)
    action = (kernel * cubic) @ psi
    action -= psi * np.vdot(psi, action)
    return float(np.vdot(action, action).real)


def _columns(cubic, support):
    return {source: {destination: value for destination, value in enumerate(cubic[:, source]) if value}
            for source in support}


def _assert_columns(actual, expected, support):
    for source in support:
        column = np.array([actual[source].get(destination, 0j)
                           for destination in range(expected.shape[0])])
        np.testing.assert_allclose(column, expected[:, source], rtol=2e-12, atol=3e-13)


def test_large_identity_shift_preserves_analytic_cubic_and_public_unit_actions():
    scalar = 2**40
    blocks = [{(): scalar, ((0, 'X'),): 1.}, {((0, 'Z'),): 1.}]
    context = response.SparseResponse(blocks, 1, 0, time=.7)
    # [X,[X,Z]]=4Z and [Z,[Z,X]]=4X, irrespective of scalar I.
    column = context.cubic_column((0, 1), 0)
    assert column.get(0, 0) == pytest.approx(-1/6, rel=2e-12)
    assert column.get(1, 0) == pytest.approx(1/3, rel=2e-12)
    assert context.lengths == (2, 1)
    # Response centering must not change the public full-H application API.
    assert context.apply(0, {0: 1.}) == {0: scalar + 0j, 1: 1. + 0j}
    assert context.apply(1, {0: 1., 1: 2j}) == {0: 1. + 0j, 1: -2j}
    assert context.hamiltonian_action({0: 1.}) == {0: scalar + 1. + 0j, 1: 1. + 0j}
    score = context.score((0, 1), 30)
    assert score.steps == 9  # Cost is 3*r+2, including the supplied identity.
    assert score.cubic_response == pytest.approx(math.sin(.7)**2 / 9, rel=2e-12)
    assert score.value == pytest.approx(math.sin(.7)**2 / (9 * 9**4), rel=2e-12)


def _shifted_case():
    base = [{((0, 'X'),): .5, ((1, 'Z'),): .25},
            {((0, 'Z'),): .75},
            {((0, 'Y'), (1, 'X')): .125}]
    shifts = (2**60, -(2**58), 3 * 2**57)
    shifted = [dict(block) for block in base]
    for block, scalar in zip(shifted, shifts):
        block[()] = scalar
    reference = {0: math.sqrt(3)/2, 3: .5j}
    return base, shifted, reference


def test_multiple_unit_identity_shifts_preserve_columns_and_coherent_transport():
    base, shifted, reference = _shifted_case()
    matrices = _literal_blocks(base, 2)  # Never form huge scalar dense products.
    context = response.SparseResponse(shifted, 2, reference, time=.73)
    assert context.lengths == (3, 2, 2)
    for order in itertools.permutations(range(3)):
        expected = _nested_cubic(matrices, order)
        actual = {source: context.cubic_column(order, source) for source in reference}
        _assert_columns(actual, expected, reference)
        expected_response = _transported_response(matrices, expected, reference, .73)
        # Supply independent columns to isolate diagonal-gap stability as well
        # as checking the complete computed-column scoring path.
        transported = context.score_columns(order, 80, _columns(expected, reference))
        complete = context.score(order, 80)
        assert transported.cubic_response == pytest.approx(expected_response, rel=2e-12, abs=3e-13)
        assert complete.cubic_response == pytest.approx(expected_response, rel=2e-12, abs=3e-13)
        assert complete.value == pytest.approx(expected_response / complete.steps**4,
                                              rel=2e-12, abs=1e-16)


def test_adjacent_updates_are_scalar_invariant_for_coherent_reference_columns():
    base, shifted, reference = _shifted_case()
    matrices = _literal_blocks(base, 2)
    context = response.SparseResponse(shifted, 2, reference, time=.73)
    for order in itertools.permutations(range(3)):
        columns = _columns(_nested_cubic(matrices, order), reference)
        saved = copy.deepcopy(columns)
        for position in (0, 1):
            new_order, updated = context.adjacent_columns(order, columns, position)
            expected_order = list(order)
            expected_order[position:position + 2] = reversed(expected_order[position:position + 2])
            assert new_order == tuple(expected_order)
            _assert_columns(updated, _nested_cubic(matrices, new_order), reference)
            assert columns == saved


def test_scalar_only_unit_keeps_rotation_cost_without_changing_the_defect():
    blocks = [{(): 2**60}, {((0, 'X'),): 1.}, {((0, 'Z'),): 1.}]
    context = response.SparseResponse(blocks, 1, 0)
    scalar_free = [np.zeros((2, 2), dtype=complex), _X, _Z]
    assert context.lengths == (1, 1, 1)
    for order in itertools.permutations(range(3)):
        column = context.cubic_column(order, 0)
        _assert_columns({0: column}, _nested_cubic(scalar_free, order), [0])
        assert context.affordable_steps(order, 30) == 7


def test_scalar_free_dynamic_range_is_resolved_or_explicitly_refused():
    magnitude = 2**40
    time = 1/magnitude
    context = response.SparseResponse(
        [{((0, 'X'),): magnitude}, {((0, 'Z'),): 1.}], 1, 0, time=time)
    numerical_resolution = getattr(response, 'NumericalResolution', ValueError)
    try:
        column = context.cubic_column((0, 1), 0)
        score = context.score((0, 1), 30)
    except numerical_resolution:
        return
    # Exact nested commutators give -M^2 Z/6 + M X/3. At T=1/M,
    # the perpendicular response is sinc(T)^2/9, approximately1/9.
    assert column.get(0, 0) == pytest.approx(-magnitude**2 / 6, rel=2e-12)
    assert column.get(1, 0) == pytest.approx(magnitude / 3, rel=2e-12)
    expected = float(np.sinc(time/np.pi)**2 / 9)
    assert score.cubic_response == pytest.approx(expected, rel=2e-12)
    assert score.value == pytest.approx(expected / score.steps**4, rel=2e-12)
    assert score.value > 0


def test_unresolved_first_seed_never_becomes_a_completed_zero_training_score():
    magnitude = 2**40
    result = response.train_ordering(
        [{((0, 'X'),): magnitude}, {((0, 'Z'),): 1.}], 1, 0,
        time=1/magnitude, rotation_budget=30, cpu_seconds=5.,
        max_evaluations=1, max_moves=0, starts=[(0, 1)])
    if result.score is None:
        assert result.status == 'numerical_limit'
        assert result.evaluations == 0
    else:
        assert result.evaluations == 1
        assert result.score.cubic_response == pytest.approx(1/9, rel=2e-12)
        assert result.score.value > 0


def test_numerical_refusal_retains_the_earlier_fully_scored_seed(monkeypatch):
    numerical_resolution = getattr(response, 'NumericalResolution', ValueError)
    original = response.SparseResponse.score
    calls = []

    def first_complete_then_refuse(self, order, rotation_budget, deadline=math.inf):
        calls.append(tuple(order))
        if len(calls) > 1:
            raise numerical_resolution('deliberate unresolved later candidate')
        return original(self, order, rotation_budget, deadline)

    monkeypatch.setattr(response.SparseResponse, 'score', first_complete_then_refuse)
    result = response.train_ordering(
        [{((0, 'X'),): 1.}, {((0, 'Z'),): 1.}], 1, 0,
        time=.7, rotation_budget=30, cpu_seconds=5., max_moves=1,
        starts=[(0, 1), (1, 0)])
    assert result.status == 'numerical_limit'
    assert result.evaluations == 1
    assert result.order == (0, 1)
    assert result.score is not None
    assert result.score.cubic_response == pytest.approx(math.sin(.7)**2 / 9, rel=2e-12)


def test_shared_rescore_reserve_can_fund_a_checked_best_seed_improvement(monkeypatch):
    """Three slices individually shorter than a rescore must still be useful."""
    class Clock:
        now = 0.

        def process_time(self):
            return self.now

        def spend(self, amount, deadline):
            if self.now >= deadline or self.now + amount >= deadline:
                self.now = max(self.now, deadline)
                raise response.SearchDeadline('deterministic phase deadline')
            self.now += amount

    clock = Clock()
    monkeypatch.setattr(response._time, 'process_time', clock.process_time)
    best_seed = tuple(reversed(range(9)))
    bad_seed = tuple(range(9))
    other_seed = tuple(range(1, 9)) + (0,)
    improvement = (best_seed[1], best_seed[0]) + best_seed[2:]
    snapshot = {0: {1: 1. + 0j}}
    saved_snapshot = copy.deepcopy(snapshot)

    class Context:
        max_entries = 50000
        reference = {0: 1.}
        last_order = other_seed
        last_columns = {0: {1: 3. + 0j}}

        def __init__(self):
            self.stats = {}
            self.visited = []
            self.full_rescores = []
            self.column_requests = []

        def columns(self, order, deadline=math.inf):
            self.column_requests.append(order)
            clock.spend(4., deadline)
            return copy.deepcopy(snapshot)

        def adjacent_columns(self, order, columns, position, deadline=math.inf):
            clock.spend(.03, deadline)
            self.visited.append(order)
            trial = list(order)
            trial[position:position + 2] = reversed(trial[position:position + 2])
            return tuple(trial), copy.deepcopy(columns)

        def affordable_steps(self, order, budget):
            return 1

        def score_columns(self, order, budget, columns, deadline=math.inf):
            clock.spend(.02, deadline)
            value = 1. if order == improvement else 10.
            return response.Score(value, value, 1)

        def score(self, order, budget, deadline=math.inf):
            self.full_rescores.append((order, deadline))
            clock.spend(4., deadline)
            value = 1. if order == improvement else 10.
            return response.Score(value, value, 1)

    context = Context()
    best_score = response.Score(10., 10., 1)
    scored = [(bad_seed, response.Score(30., 30., 1)),
              (best_seed, best_score), (other_seed, response.Score(20., 20., 1))]
    # Keep the baseline reproduction behavioral: before the optional snapshot
    # parameter exists, failure should expose unused search time, not TypeError.
    options = ({'seed_columns': {best_seed: snapshot}}
               if 'seed_columns' in inspect.signature(response._incremental_search).parameters else {})
    order, score, evaluations, moves, status = response._incremental_search(
        context, scored, 30, 9., 3000, 1, best_seed, best_score, 3, 4.,
        **options)
    assert context.visited and context.visited[0] == best_seed
    assert best_seed not in context.column_requests
    assert context.full_rescores == [(improvement, 9.)]
    assert order == improvement
    assert score.value == 1.
    assert evaluations > 3 and moves == 1
    assert clock.now <= 9.
    assert snapshot == saved_snapshot
