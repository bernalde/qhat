"""Independent matrix oracles for Claude's legacy operator/coefficient findings."""
from itertools import product
import numpy as np
import pytest
from openfermion import QubitOperator, get_sparse_operator
from qhat.analysis.matrix_operations import PauliStringOperator, pauli_dict_to_matrix
from qhat.analysis.trotter_coefficients_fast import trotter_error_estimator_fast


@pytest.mark.parametrize('n', [1, 2, 3])
def test_matrix_free_all_pauli_basis_columns_and_adjoint(n):
    rng = np.random.default_rng(1947 + n)
    for letters in product('IXYZ', repeat=n):
        word = ''.join(letters)
        dense = pauli_dict_to_matrix({word: 1 - 0.25j}, n)
        operator = PauliStringOperator({word: 1 - 0.25j}, n)
        for column in np.eye(2**n, dtype=complex):
            np.testing.assert_allclose(operator @ column, dense @ column, atol=1e-14)
        vector = rng.normal(size=2**n) + 1j*rng.normal(size=2**n)
        np.testing.assert_allclose(operator @ vector, dense @ vector, atol=1e-14)
        np.testing.assert_allclose(operator.H @ vector, dense.conj().T @ vector, atol=1e-14)


def dense_coefficients(terms):
    matrices = [get_sparse_operator(term, n_qubits=2).toarray() for term in terms]
    def comm(a, b):
        return a @ b - b @ a
    c1 = sum(np.linalg.norm(comm(a, b), 2) for k, a in enumerate(matrices)
             for b in matrices[k+1:]) / 2
    c21 = sum(np.linalg.norm(comm(a, comm(b, h)), 2)
              for k, h in enumerate(matrices) for a in matrices[k+1:]
              for b in matrices[k+1:])
    c22 = sum(np.linalg.norm(comm(h, comm(h, b)), 2)
              for k, h in enumerate(matrices) for b in matrices[k+1:])
    return c1, c21/12 + c22/24


@pytest.mark.parametrize('words,coefficients,expected', [
    (['X0', 'Z0'], [1, 1], (1, 1/2)),
    (['X0', 'Z0'], [1, 2], (2, 5/3)),
    (['Z0', 'X0'], [2, 1], (2, 4/3)),
    (['X0', 'Z0', 'X0 X1'], [1, 1, 1], (2, 4/3)),
])
def test_ordered_nested_coefficient_against_dense_oracle(words, coefficients, expected):
    terms = [QubitOperator(word, coefficient) for word, coefficient in zip(words, coefficients)]
    oracle = dense_coefficients(terms)
    np.testing.assert_allclose(oracle, expected, atol=1e-14)
    actual = trotter_error_estimator_fast(terms, 10, mode='exact')
    np.testing.assert_allclose(actual, oracle, atol=1e-14)


@pytest.mark.parametrize('n', [2, 3, 5])
def test_uniform_rejection_sampler_keeps_each_ordered_triple(n):
    from qhat.analysis.trotter_coefficients_fast import sample_ordered_triples
    population = np.array(list(product(range(n), repeat=3)))
    class EnumeratedDraws:
        offset = 0
        def randint(self, low, high, size):
            assert (low, high) == (0, n)
            count = size[0]
            result = population[np.arange(self.offset, self.offset+count) % len(population)]
            self.offset += count
            return result
    expected = {(i,j,k) for k in range(n) for i in range(k+1,n) for j in range(k+1,n)}
    actual = sample_ordered_triples(n, len(expected), EnumeratedDraws())
    assert len(actual) == len(expected)
    assert set(map(tuple, actual)) == expected


def test_uniform_pair_sampler_preserves_equal_pair_weights():
    from collections import Counter
    from qhat.analysis.trotter_coefficients_fast import sample_ordered_pairs
    class EnumeratedDraws:
        calls = 0
        def randint(self, low, high, size):
            self.calls += 1
            return np.repeat(np.arange(4), 3) if self.calls == 1 else np.tile(np.arange(3), 4)
    counts = Counter(map(tuple, sample_ordered_pairs(4, 12, EnumeratedDraws())))
    assert counts == {(i,j):2 for i in range(4) for j in range(i+1,4)}


def test_reference_monte_carlo_includes_two_term_diagonal():
    from qhat.analysis.trotter_coefficients import trotter_error_estimator
    terms = [QubitOperator('X0'), QubitOperator('Z0', 2)]
    np.testing.assert_allclose(trotter_error_estimator(terms, .03, batch_size=10), (2, 5/3))


def test_sampled_tracking_preserves_swapped_nested_terms():
    import qhat.analysis.trotter_coefficients_fast as fast
    terms = [QubitOperator('X0'), QubitOperator('Z0'), QubitOperator('X0 X1')]
    np.random.seed(276)
    actual = fast.trotter_error_estimator_fast(terms, 1, batch_size=100, mode='monte_carlo')
    np.testing.assert_allclose(actual, dense_coefficients(terms), atol=1e-14)


def test_matrix_free_real_vectors_and_matrix_columns_keep_y_phases():
    dense = pauli_dict_to_matrix({'YI':1,'IX':.3},2)
    operator = PauliStringOperator({'YI':1,'IX':.3},2)
    vector = np.array([1,0,0,0])
    np.testing.assert_allclose(operator @ vector,dense @ vector)
    np.testing.assert_allclose(operator @ np.eye(4),dense)
    np.testing.assert_allclose(operator.H @ vector,dense.conj().T @ vector)
