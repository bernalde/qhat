"""Independent small dense oracles for the bounded sparse ordering path.

The cubic oracle uses nested commutators, rather than the implementation's
Taylor-vector recurrence. Only the explicit large-register test exceeds two
qubits; it compares a four-configuration reachable subspace with a dense model.
"""
import itertools
import math
from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import expm

from qhat.analysis.sparse_response_ordering import (
    ResourceLimit,
    SearchDeadline,
    SparseResponse,
    conservative_steps,
    molecular_units,
    train_ordering,
)


_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.diag([1., -1.]).astype(complex),
}
_XZ = [{((0, "X"),): 1.}, {((0, "Z"),): 1.}]
_MIXED = [
    {((0, "X"),): .375},
    {((0, "Z"),): -.5, ((1, "Z"),): .25},
    {((1, "Y"),): .625},
]


def test_later_restart_can_accept_multiple_moves_after_an_earlier_slice_expires(monkeypatch):
    """A prior seed's stop flag must not truncate later local searches."""
    from types import SimpleNamespace
    from qhat.analysis import sparse_response_ordering as native

    early = (3, 2, 1, 0)
    start = (0, 1, 2, 3)
    first_gain = (1, 0, 2, 3)
    second_gain = (1, 0, 3, 2)
    values = {early: 1., start: 10., first_gain: 8., second_gain: .5}

    class Context:
        stats = {}

        def columns(self, order, deadline):
            return order

        def adjacent_columns(self, order, columns, position, deadline):
            trial = list(order)
            trial[position], trial[position+1] = trial[position+1], trial[position]
            return tuple(trial), tuple(trial)

        def affordable_steps(self, order, budget):
            return 1

        def score_columns(self, order, budget, columns, deadline):
            return native.Score(values.get(order, 100.), 0., 1)

        def score(self, order, budget, deadline):
            return self.score_columns(order, budget, order, deadline)

    def check(deadline):
        if deadline < 20.:
            raise native.SearchDeadline('first seed exhausted its own slice')

    monkeypatch.setattr(native, '_time', SimpleNamespace(process_time=lambda: 0.))
    monkeypatch.setattr(native, '_check', check)
    best, score, _, moves, status = native._incremental_search(
        Context(), [(early, native.Score(1., 0., 1)), (start, native.Score(10., 0., 1))],
        100, 20., 1000, 5, early, native.Score(1., 0., 1), 0, 0.)
    assert best == second_gain
    assert score.value == .5
    assert moves == 2
    assert status == 'cpu_limit'  # Preserve the earlier truncation in metadata.


def _dense_word(word, n_qubits):
    # Deliberately do not use sparse bit masks or the implementation's action.
    assert n_qubits <= 6, "dense oracle is restricted to small registers"
    letters = dict(word)
    result = np.ones((1, 1), dtype=complex)
    for qubit in reversed(range(n_qubits)):
        result = np.kron(result, _PAULI[letters.get(qubit, "I")])
    return result


def _dense_blocks(blocks, n_qubits):
    return [sum((complex(c) * _dense_word(word, n_qubits)
                 for word, c in block.items()),
                np.zeros((2**n_qubits, 2**n_qubits), dtype=complex))
            for block in blocks]


def _dense_reference(reference, n_qubits):
    vector = np.zeros(2**n_qubits, dtype=complex)
    if isinstance(reference, int):
        vector[reference] = 1.
    else:
        for occupation, amplitude in reference.items():
            vector[occupation] = amplitude
    return vector


def _commutator(a, b):
    return a @ b - b @ a


def _dense_cubic(matrices, order):
    """Coefficient of h^3 in the positive-H symmetric formula minus exp(hH)."""
    result = np.zeros_like(matrices[0])
    for position, unit in enumerate(order[:-1]):
        a = matrices[unit]
        suffix = sum((matrices[index] for index in order[position + 1:]),
                     np.zeros_like(a))
        result -= _commutator(a, _commutator(a, suffix)) / 24
        result += _commutator(suffix, _commutator(suffix, a)) / 12
    return result


def _fused_rotation_cost(blocks, order, steps):
    """Count actual adjacent fusions by enumerating the repeated stage list."""
    if not order:
        return 0
    stages = list(order[:-1]) + [order[-1]] + list(reversed(order[:-1]))
    expanded = stages * steps
    fused = [unit for position, unit in enumerate(expanded)
             if position == 0 or unit != expanded[position - 1]]
    return sum(sum(c != 0 for c in blocks[unit].values()) for unit in fused)


def _dense_steps(blocks, order, budget):
    if not any(c != 0 for block in blocks for c in block.values()):
        return 1
    if len(order) == 1:
        assert _fused_rotation_cost(blocks, order, 1) <= budget
        return 1
    steps = 0
    while _fused_rotation_cost(blocks, order, steps + 1) <= budget:
        steps += 1
    assert steps > 0
    return steps


def _dense_proxy(blocks, n_qubits, reference, order, time, budget):
    matrices = _dense_blocks(blocks, n_qubits)
    cubic = _dense_cubic(matrices, order)
    hamiltonian = sum(matrices)
    energies = np.diag(hamiltonian).real
    gaps = energies[:, None] - energies[None, :]
    kernel = (time * np.sinc(time * gaps / (2 * np.pi))
              * np.exp(.5j * time * gaps))
    psi = _dense_reference(reference, n_qubits)
    transported = (kernel * cubic) @ psi
    perpendicular = transported - psi * np.vdot(psi, transported)
    response = float(np.vdot(perpendicular, perpendicular).real)
    steps = _dense_steps(blocks, order, budget)
    return cubic, response, steps, response / steps**4


def _dense_symmetric(matrices, order, time):
    stages = ([(unit, .5) for unit in order[:-1]]
              + [(order[-1], 1.)]
              + [(unit, .5) for unit in reversed(order[:-1])])
    result = np.eye(matrices[0].shape[0], dtype=complex)
    for unit, scale in stages:
        result = expm(-1j * time * scale * matrices[unit]) @ result
    return result


@pytest.mark.parametrize("letter, expected", [
    ("X", {0: -1/6, 1: 1/3}),
    ("Y", {0: -1/6, 1: 1j/3}),
])
def test_hand_derived_cubic_column_and_complex_y_sign(letter, expected):
    blocks = [{((0, letter),): 1.}, {((0, "Z"),): 1.}]
    context = SparseResponse(blocks, 1, 0)
    column = context.cubic_column((0, 1), 0)
    assert set(column) == set(expected)
    for destination, coefficient in expected.items():
        assert column[destination] == pytest.approx(coefficient, abs=2e-15)
    reverse = context.cubic_column((1, 0), 0)
    assert reverse[0] == pytest.approx(1/3, abs=2e-15)
    assert reverse[1] == pytest.approx(-1j/6 if letter == "Y" else -1/6,
                                       abs=2e-15)


@pytest.mark.parametrize("reference", [
    0,
    {0: 1 / math.sqrt(3), 1: 1j / math.sqrt(3), 3: 1 / math.sqrt(3)},
])
@pytest.mark.parametrize("cache_entries", [0, 2, 128])
def test_all_permutations_match_commutator_and_coherent_dense_oracles(reference,
                                                                    cache_entries):
    context = SparseResponse(_MIXED, 2, reference, time=.73,
                             max_cached_actions=cache_entries)
    support = [reference] if isinstance(reference, int) else list(reference)
    for order in itertools.permutations(range(3)):
        cubic, response, steps, value = _dense_proxy(
            _MIXED, 2, reference, order, .73, 40)
        for source in support:
            actual = context.cubic_column(order, source)
            vector = np.array([actual.get(k, 0j) for k in range(4)])
            np.testing.assert_allclose(vector, cubic[:, source], rtol=1e-11,
                                       atol=2e-13)
        score = context.score(order, 40)
        assert score.steps == steps
        assert score.cubic_response == pytest.approx(response, rel=1e-11, abs=2e-13)
        assert score.value == pytest.approx(value, rel=1e-11, abs=2e-16)
    assert context.stats["peak_cached_actions"] <= cache_entries
    if not cache_entries:
        assert not context.cache


def test_complex_reference_columns_interfere_before_projection():
    plus = {0: math.sqrt(3)/2, 1: .5j}
    minus = {0: math.sqrt(3)/2, 1: -.5j}
    responses = []
    for reference in (plus, minus):
        expected = _dense_proxy(_XZ, 1, reference, (0, 1), .6, 30)[1]
        actual = SparseResponse(_XZ, 1, reference, time=.6).score((0, 1), 30)
        assert actual.cubic_response == pytest.approx(expected, rel=1e-12)
        responses.append(actual.cubic_response)
    # Identical occupation probabilities are insufficient to predict response.
    assert abs(responses[0] - responses[1]) > 1e-3
    phase = np.exp(.43j)
    phased = {key: phase * value for key, value in plus.items()}
    assert SparseResponse(_XZ, 1, phased, time=.6).score((0, 1), 30).value == \
        pytest.approx(SparseResponse(_XZ, 1, plus, time=.6).score((0, 1), 30).value)


def test_physical_exponential_has_positive_i_times_the_cubic_column():
    matrices = _dense_blocks(_MIXED, 2)
    order = (2, 0, 1)
    h = .03
    cubic = _dense_cubic(matrices, order)
    difference = (_dense_symmetric(matrices, order, h)
                  - expm(-1j * h * sum(matrices)) - 1j * h**3 * cubic)
    coefficient_norm = sum(abs(c) for block in _MIXED for c in block.values())
    # Independent absolute bound for the two degree>=4 Taylor tails.
    z = h * coefficient_norm
    assert np.linalg.norm(difference, 2) <= math.exp(z) * z**4 / 12 + 1e-13


def test_137_qubit_response_only_visits_four_configurations():
    high = 136
    spectator = (1 << 128) | (1 << 42)
    blocks = [{((0, "Y"),): .3}, {((high, "X"),): -.2},
              {((0, "Z"), (high, "Z")): .7}]
    small = [{((0, "Y"),): .3}, {((1, "X"),): -.2},
             {((0, "Z"), (1, "Z")): .7}]
    context = SparseResponse(blocks, 137, spectator, time=.8,
                             max_entries=128, max_cached_actions=8)
    order = (2, 0, 1)
    expected = _dense_proxy(small, 2, 0, order, .8, 40)
    actual = context.score(order, 40)
    assert actual.cubic_response == pytest.approx(expected[1], rel=1e-12)
    assert actual.value == pytest.approx(expected[3], rel=1e-12)
    reachable = {spectator, spectator ^ 1, spectator ^ (1 << high),
                 spectator ^ 1 ^ (1 << high)}
    assert set(context.cubic_column(order, spectator)) <= reachable
    assert context.stats["peak_vector_entries"] <= 4
    assert context.stats["peak_cached_actions"] <= 8


def test_molecular_adapter_preserves_coefficients_and_conserving_units():
    terms = {
        (): .125 + 0j,
        ((0, "Z"),): -.5 + 0j,
        ((1, "Z"),): .25 + 0j,
        ((0, "X"), (1, "X")): .125 + 0j,
        ((0, "Y"), (1, "Y")): .125 + 0j,
    }
    units = molecular_units(terms, 2)
    assert len(units) == 4
    assert sorted(map(len, units)) == [1, 1, 1, 2]
    matrices = _dense_blocks(units, 2)
    np.testing.assert_allclose(sum(matrices), _dense_blocks([terms], 2)[0], atol=0)
    number = np.diag([0, 1, 1, 2])
    for matrix in matrices:
        np.testing.assert_allclose(_commutator(matrix, number), 0, atol=0)
    order = tuple(reversed(range(len(units))))
    expected = _dense_proxy(units, 2, 1, order, .7, 40)
    score = SparseResponse(units, 2, 1, time=.7).score(order, 40)
    assert score.cubic_response == pytest.approx(expected[1], abs=2e-14)
    assert score.value == pytest.approx(expected[3], abs=2e-17)


def test_adapter_accepts_hermitian_odd_y_but_rejects_noncommuting_flip_group():
    assert molecular_units({((0, "Y"),): .5 + 0j}, 1) == ({((0, "Y"),): .5 + 0j},)
    with pytest.raises(ValueError, match="real"):
        molecular_units({((0, "Y"),): .5 + .1j}, 1)
    with pytest.raises(ValueError, match="commute"):
        molecular_units({((0, "X"),): 1., ((0, "Y"),): 1.}, 1)


@pytest.mark.parametrize("blocks,n_qubits", [
    (_XZ, 0), (_XZ, -1), (_XZ, True), (_XZ, 1.5),
    ([{((-1, "X"),): 1.}], 1),
    ([{((1, "X"),): 1.}], 1),
    ([{((0, "X"), (0, "Y")): 1.}], 1),
    ([{((0, "Q"),): 1.}], 1),
    ([{((0, "X"),): 1j}], 1),
    ([{((0, "X"),): float("nan")}], 1),
    ([{((0, "X"),): float("inf")}], 1),
    ([{((0, "X"),): True}], 1),
    ([{((0, "X"),): 1., ((0, "Z"),): 1.}], 1),
])
def test_invalid_pauli_domains_are_rejected(blocks, n_qubits):
    with pytest.raises(ValueError):
        SparseResponse(blocks, n_qubits, 0)


@pytest.mark.parametrize("reference", [
    -1, 2, True, {}, {0: .5}, {0: float("nan")}, {2: 1.},
    {0: 1., 1: 1.}, "0",
])
def test_invalid_reference_is_not_accepted(reference):
    with pytest.raises(ValueError):
        SparseResponse(_XZ, 1, reference)


def test_reference_support_limit_and_working_support_refuse_truncation():
    with pytest.raises(ValueError, match="16"):
        SparseResponse([], 5, {key: 1/math.sqrt(17) for key in range(17)})
    with pytest.raises(ResourceLimit):
        SparseResponse([], 1, {0: 1/math.sqrt(2), 1: 1j/math.sqrt(2)}, max_entries=1)
    with pytest.raises(ResourceLimit):
        context = SparseResponse([{((0, "X"),): 1., ((1, "X"),): 1.},
                                  {((0, "Z"), (1, "Z")): .5}], 2, 0,
                                 max_entries=1)
        context.score((0, 1), 10)


def test_single_commuting_unit_needs_no_expanding_cubic_column():
    context = SparseResponse([{((0, "X"),): 1., ((1, "X"),): 1.}], 2, 0,
                             max_entries=1)
    assert context.score((0,), 10).value == 0


@pytest.mark.parametrize("order", [(0, 0), (0,), (0, 2), (False, 1), (0, 1, 2)])
def test_invalid_permutation_is_rejected(order):
    context = SparseResponse(_XZ, 1, 0)
    with pytest.raises(ValueError, match="permutation"):
        context.score(order, 20)


@pytest.mark.parametrize("kwargs", [
    {"max_entries": 0}, {"max_entries": True},
    {"max_cached_actions": -1}, {"time": 0},
    {"time": -1}, {"time": float("inf")},
])
def test_invalid_context_resource_and_time_settings(kwargs):
    with pytest.raises(ValueError):
        SparseResponse(_XZ, 1, 0, **kwargs)


@pytest.mark.parametrize("budget", [0, 2, -1, True, 3.5])
def test_unaffordable_or_invalid_rotation_budget_is_rejected(budget):
    with pytest.raises(ValueError):
        SparseResponse(_XZ, 1, 0).score((0, 1), budget)


def test_expired_deadlines_do_not_become_zero_response():
    with pytest.raises(SearchDeadline):
        SparseResponse(_XZ, 1, 0, deadline=-math.inf)
    context = SparseResponse(_XZ, 1, 0)
    with pytest.raises(SearchDeadline):
        context.score((0, 1), 30, deadline=-math.inf)
    with pytest.raises(SearchDeadline):
        context.cubic_column((0, 1), 0, deadline=-math.inf)
    with pytest.raises(ValueError, match="support"):
        context.cubic_column((0, 1), 1)


@pytest.mark.parametrize("allowance", [0, -1, float("inf"), float("nan"), True])
def test_training_requires_a_positive_finite_cpu_allowance(allowance):
    with pytest.raises(ValueError):
        train_ordering(_XZ, 1, 0, rotation_budget=30, cpu_seconds=allowance)


def test_zero_evaluations_returns_no_invented_score():
    result = train_ordering(_XZ, 1, 0, rotation_budget=30, cpu_seconds=5.,
                            max_evaluations=0, starts=[(0, 1)])
    assert result.order == (0, 1)
    assert result.score is None
    assert result.evaluations == result.accepted_moves == 0
    assert result.status == "evaluation_limit"
    assert result.stats["response_is_error_bound"] is False


def test_completed_improvement_survives_an_interrupted_insertion_sweep():
    result = train_ordering(_XZ, 1, 0, time=.6, rotation_budget=30,
                            cpu_seconds=5., max_evaluations=2,
                            max_moves=5, starts=[(0, 1)])
    expected = _dense_proxy(_XZ, 1, 0, (1, 0), .6, 30)
    assert result.order == (1, 0)
    assert result.score is not None
    assert result.score.value == pytest.approx(expected[3], rel=1e-12)
    assert result.evaluations == 2
    assert result.accepted_moves == 0  # The full insertion sweep did not finish.
    assert result.status == "evaluation_limit"


def test_completed_insertion_sweep_records_a_real_improving_move():
    result = train_ordering(_XZ, 1, 0, time=.6, rotation_budget=30,
                            cpu_seconds=5., max_evaluations=8,
                            max_moves=2, starts=[(0, 1)])
    initial = _dense_proxy(_XZ, 1, 0, (0, 1), .6, 30)[3]
    assert result.order == (1, 0)
    assert result.score.value < initial
    assert result.accepted_moves == 1
    assert result.status == "complete"


@pytest.mark.parametrize("metric,tolerance", [
    ("operator", Fraction(1, 100)),
    ("infidelity", Fraction(1, 10000)),
    ("squared_infidelity", Fraction(1, 10000)),
])
def test_conservative_certificate_bounds_independent_physical_evolution(metric,
                                                                      tolerance):
    blocks = [{((0, "X"),): Fraction(3, 10)},
              {((0, "Z"),): Fraction(1, 5)}]
    time = Fraction(7, 10)
    certificate = conservative_steps(blocks, 1, time=time,
                                     error_tolerance=tolerance, metric=metric)
    assert isinstance(certificate.operator_error_bound, Fraction)
    assert isinstance(certificate.error_bound, Fraction)
    assert certificate.metric == metric
    assert 0 <= certificate.error_bound <= tolerance
    matrices = _dense_blocks(blocks, 1)
    actual = np.linalg.matrix_power(
        _dense_symmetric(matrices, (0, 1), float(time)/certificate.steps),
        certificate.steps)
    exact = expm(-1j * float(time) * sum(matrices))
    operator_error = np.linalg.norm(actual - exact, 2)
    assert operator_error <= float(certificate.operator_error_bound) + 2e-13
    psi = np.array([math.sqrt(3)/2, .5j])
    overlap = abs(np.vdot(exact @ psi, actual @ psi))
    observed = {"operator": operator_error,
                "infidelity": 1-overlap,
                "squared_infidelity": 1-overlap**2}[metric]
    assert observed <= float(certificate.error_bound) + 2e-13
    if metric != "operator":
        denominator = 2 if metric == "infidelity" else 1
        assert certificate.error_bound == certificate.operator_error_bound**2 / denominator


@pytest.mark.parametrize("blocks,time,n_qubits", [
    ([], 5., 1), ([{}], 1., 1), (_XZ, 0., 1),
    ([{((0, "Y"),): .4, ((1, "X"),): -.2}], -.7, 2),
])
def test_exact_zero_and_single_commuting_unit_certificates(blocks, time, n_qubits):
    certificate = conservative_steps(blocks, n_qubits, time=time,
                                     error_tolerance=Fraction(1, 10**30), max_steps=1)
    assert certificate.steps == 1
    assert certificate.operator_error_bound == certificate.error_bound == 0


def test_step_certificate_refuses_an_unaffordable_target():
    with pytest.raises(ResourceLimit):
        conservative_steps(_XZ, 1, time=1., error_tolerance=Fraction(1, 10**12),
                           max_steps=1)


@pytest.mark.parametrize("kwargs", [
    {"error_tolerance": 0}, {"error_tolerance": -1},
    {"error_tolerance": float("nan")}, {"error_tolerance": float("inf")},
    {"metric": "response"}, {"max_steps": 0}, {"max_steps": True},
    {"time": float("inf")},
])
def test_invalid_certificate_settings_are_rejected(kwargs):
    settings = dict(time=1., error_tolerance=.01)
    settings.update(kwargs)
    with pytest.raises(ValueError):
        conservative_steps(_XZ, 1, **settings)


class _IncrementalBudgetContext:
    """Two literal objective values; exercise search accounting, no physics mock."""
    def __init__(self, *, rescore_deadline=False, slice_once=False):
        from qhat.analysis import sparse_response_ordering as native
        self.stats = {}
        self.values = {(0, 1): native.Score(1., 1., 1),
                       (1, 0): native.Score(2., 2., 1)}
        self.columns_calls = []
        self.rescore_deadline = rescore_deadline
        self.slice_once = slice_once

    def columns(self, order, deadline):
        self.columns_calls.append(tuple(order))
        if self.slice_once and len(self.columns_calls) == 1:
            raise SearchDeadline('first seed local slice expired')
        return tuple(order)

    def adjacent_columns(self, order, columns, position, deadline):
        result = list(order)
        result[position], result[position + 1] = result[position + 1], result[position]
        return tuple(result), tuple(result)

    def affordable_steps(self, order, budget):
        return 1

    def score_columns(self, order, budget, columns, deadline):
        return self.values[tuple(order)]

    def score(self, order, budget, deadline):
        if self.rescore_deadline:
            raise SearchDeadline('external whole-call deadline')
        return self.values[tuple(order)]


def test_incremental_full_rescore_deadline_is_not_a_support_refusal():
    from qhat.analysis import sparse_response_ordering as native
    context = _IncrementalBudgetContext(rescore_deadline=True)
    initial = (1, 0)
    result = native._incremental_search(context, [(initial, context.values[initial])],
        10, native._time.process_time() + 10., 8, 1,
        initial, context.values[initial], 1, .01)
    assert result[0] == initial and result[1] == context.values[initial]
    assert result[4] == 'cpu_limit'
    assert context.stats['full_move_rescore_attempts'] == 1
    assert context.stats['full_move_rescores'] == 0


def test_incremental_reserved_last_evaluation_does_not_rebuild_other_seeds():
    from qhat.analysis import sparse_response_ordering as native
    context = _IncrementalBudgetContext()
    first, second = (0, 1), (1, 0)
    result = native._incremental_search(context,
        [(first, context.values[first]), (second, context.values[second])],
        10, native._time.process_time() + 10., 3, 1,
        first, context.values[first], 2, .01, seed_columns={first: first})
    assert context.columns_calls == []
    assert result[2:5] == (2, 0, 'evaluation_limit')
    assert context.stats['incremental_candidates'] == 0
    assert context.stats['full_move_rescore_attempts'] == 0


@pytest.mark.parametrize('seed_count,max_moves,expected', [
    (2, 1, 'evaluation_limit'), (1, 1, 'complete'), (2, 0, 'complete'),
])
def test_complete_sweep_at_evaluation_boundary_reports_only_skipped_work(
        monkeypatch, seed_count, max_moves, expected):
    from qhat.analysis import sparse_response_ordering as native
    context = _IncrementalBudgetContext()
    monkeypatch.setattr(native, '_check', lambda deadline: None)
    monkeypatch.setattr(native._time, 'process_time', lambda: 0.)
    first, second = (0, 1), (1, 0)
    seeds = [first, second][:seed_count]
    candidates = 2 if max_moves else 0
    result = native._incremental_search(context,
        [(order, context.values[order]) for order in seeds], 10, 100.,
        seed_count+candidates+1, max_moves, first, context.values[first], seed_count, 0.)
    assert result[:2] == (first, context.values[first])
    assert result[2:] == (seed_count+candidates, 0, expected)
    assert context.columns_calls == ([first] if max_moves else [])
    assert context.stats['incremental_candidates'] == candidates
    assert context.stats['full_move_rescore_attempts'] == 0


@pytest.mark.parametrize('exception,status', [
    ('NumericalResolution', 'numerical_limit'), ('ResourceLimit', 'support_limit'),
])
def test_evaluation_boundary_preserves_full_rescore_failure(monkeypatch, exception, status):
    from qhat.analysis import sparse_response_ordering as native
    context = _IncrementalBudgetContext()
    monkeypatch.setattr(native, '_check', lambda deadline: None)
    monkeypatch.setattr(native._time, 'process_time', lambda: 0.)
    seeds = [(0, 1, 2), (1, 0, 2), (2, 1, 0)]
    context.values = {order: native.Score(float(index), float(index), 1)
                      for index, order in enumerate(seeds, 1)}
    context.score_columns = lambda order, *args: context.values.get(order, native.Score(4., 4., 1))
    def refused_rescore(*args):
        raise getattr(native, exception)('literal boundary refusal')
    context.score = refused_rescore
    result = native._incremental_search(context,
        [(order, context.values[order]) for order in seeds], 10, 100., 16, 1,
        seeds[0], context.values[seeds[0]], 3, 0.)
    assert result[:2] == (seeds[0], context.values[seeds[0]])
    assert result[2:] == (15, 0, status)
    assert context.columns_calls == seeds[:2]
    assert context.stats['full_move_rescore_attempts'] == 1


def test_local_seed_slice_expiry_preserves_later_global_budget_allocation():
    from qhat.analysis import sparse_response_ordering as native
    context = _IncrementalBudgetContext(slice_once=True)
    first, second = (0, 1), (1, 0)
    result = native._incremental_search(context,
        [(first, context.values[first]), (second, context.values[second])],
        10, native._time.process_time() + 10., 8, 1,
        first, context.values[first], 2, .01)
    assert context.columns_calls == [first, second]
    assert result[0] == first and result[3] == 1
    assert context.stats['full_move_rescores'] == 1


def test_interrupted_admission_raises_before_returning_an_unvalidated_initializer():
    from qhat.analysis import sparse_response_ordering as native
    def exhausted():
        raise SearchDeadline('whole-call admission deadline')
    token = native._EXTERNAL_BUDGET_CHECK.set(exhausted)
    try:
        with pytest.raises(SearchDeadline):
            native.train_ordering(_XZ, 1, 0, rotation_budget=30,
                                  cpu_seconds=5., max_evaluations=0)
    finally:
        native._EXTERNAL_BUDGET_CHECK.reset(token)
