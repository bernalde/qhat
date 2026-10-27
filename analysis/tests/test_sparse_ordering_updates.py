"""Adjacent-response updates checked against full independent matrix algebra."""
import itertools
import numpy as np

from qhat.analysis.sparse_response_ordering import SparseResponse, train_ordering


def comm(a, b):
    return a @ b - b @ a


def dense_cubic(units, order):
    result = np.zeros_like(units[0])
    for at, i in enumerate(order):
        rest = sum((units[j] for j in order[at + 1:]), np.zeros_like(result))
        result += -comm(units[i], comm(units[i], rest)) / 24 + comm(rest, comm(rest, units[i])) / 12
    return result


def test_adjacent_swaps_with_nonzero_suffix_and_reverse_insertion():
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    units = [x * .75, y * .5, z * -.625]
    blocks = [{((0, letter),): coefficient}
              for letter, coefficient in [('X', .75), ('Y', .5), ('Z', -.625)]]
    reference = {0: 1 / np.sqrt(2), 1: 1j / np.sqrt(2)}
    context = SparseResponse(blocks, 1, reference, time=.9)
    for order in itertools.permutations(range(3)):
        columns = context.columns(order)
        for position in (0, 1):
            changed, updated = context.adjacent_columns(order, columns, position)
            expected = dense_cubic(units, changed)
            for source in reference:
                actual = np.array([updated[source].get(i, 0) for i in range(2)])
                np.testing.assert_allclose(actual, expected[:, source], atol=2e-15)
            assert abs(context.score_columns(changed, 90, updated).value - context.score(changed, 90).value) < 1e-16
            restored, reversed_columns = context.adjacent_columns(changed, updated, position)
            assert restored == order
            for source in reference:
                for target in range(2):
                    assert abs(reversed_columns[source].get(target, 0) - columns[source].get(target, 0)) < 2e-15
        trial, updated = order, columns
        for position in (0, 1, 0):
            trial, updated = context.adjacent_columns(trial, updated, position)
            for source in reference:
                actual = np.array([updated[source].get(i, 0) for i in range(2)])
                np.testing.assert_allclose(actual, dense_cubic(units, trial)[:, source], atol=3e-15)


def test_larger_training_returns_fully_rescored_incumbent():
    blocks = [{((q % 2, 'X' if q % 3 else 'Z'),): (q + 1) / 16} for q in range(9)]
    result = train_ordering(blocks, 2, 0, time=.7, rotation_budget=300,
                            cpu_seconds=2., max_evaluations=35, max_moves=2)
    assert sorted(result.order) == list(range(9))
    assert result.score is not None
    full = SparseResponse(blocks, 2, 0, time=.7).score(result.order, 300)
    assert result.score == full
    assert result.evaluations <= 35
    assert result.stats['incremental_candidates'] > 0


def test_many_y_phases_are_exact_and_column_cache_is_not_exposed():
    count = 137
    word = tuple((q, 'Y') for q in range(count))
    context = SparseResponse([{word: 1.}], count, 0)
    full = (1 << count) - 1
    assert context.apply(0, {0: 1.}) == {full: 1j}
    assert context.apply(0, {full: 1.}) == {0: -1j}
    small = SparseResponse([{((0, 'X'),): 1.}, {((0, 'Z'),): 1.}], 1, 0)
    expected = small.score((0, 1), 30)
    columns = small.columns((0, 1));columns[0].clear()
    assert small.score((0, 1), 30) == expected


def test_huge_finite_reference_is_rejected_as_invalid_norm():
    import pytest
    with pytest.raises(ValueError, match='normalized'):
        SparseResponse([{((0, 'X'),): 1.}], 1, {0: 1e300})
