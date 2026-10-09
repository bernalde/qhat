"""Independent small-matrix checks for Albert's imported molecular schedules.

The dense fermion oracle uses occupation-basis ladder action directly. It does
not share normal-ordering or Jordan--Wigner coefficient accumulation with the
benchmark. Tests marked xfail document unmodified historical precision bugs.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from openfermion import FermionOperator, InteractionOperator, QubitOperator
from openfermion import get_sparse_operator, hermitian_conjugated, jordan_wigner

from qhat.analysis import benchmark_L_sweep_trotter as dense
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as matrix_free
from qhat.analysis import benchmark_b2_coloring_robustness as coloring
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis import benchmark_fermionic_structure_ablation as ablation
from qhat.analysis.fermionic_coefficients import (
    fermion_from_interaction, jordan_wigner_molecular, normal_ordered_molecular,
)


TOL = 1e-12


@lru_cache(None)
def ladder(mode, creation, n_qubits):
    result = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    bit = 1 << (n_qubits - 1 - mode)
    for source in range(2**n_qubits):
        occupied = bool(source & bit)
        if occupied != bool(creation):
            parity = (source >> (n_qubits - mode)).bit_count() % 2
            result[source ^ bit, source] = (-1) ** parity
    return result


def fermion_matrix(operator, n_qubits):
    result = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for word, coefficient in operator.terms.items():
        term = np.eye(2**n_qubits, dtype=complex)
        for mode, creation in word:
            term = term @ ladder(mode, creation, n_qubits)
        result += coefficient * term
    return result


def pauli_matrix(key, n_qubits):
    # Unit coefficient avoids OpenFermion's hidden pruning in this oracle.
    return get_sparse_operator(QubitOperator(key, 1.0), n_qubits=n_qubits).toarray()


def qubit_matrix(operator, n_qubits):
    return sum((coefficient * pauli_matrix(key, n_qubits)
                for key, coefficient in operator.terms.items()),
               np.zeros((2**n_qubits, 2**n_qubits), dtype=complex))


def formula(matrices, order, dt):
    if order == 4:
        p = 1 / (4 - 4 ** (1 / 3))
        outer = formula(matrices, 2, p * dt)
        return outer @ outer @ formula(matrices, 2, (1 - 4 * p) * dt) @ outer @ outer
    factors = [expm(-1j * dt * h / (2 if order == 2 else 1)) for h in matrices]
    if order == 2:
        factors += factors[::-1]
    result = np.eye(len(matrices[0]), dtype=complex)
    for factor in factors:
        result = factor @ result
    return result


def raw_molecular_fixture():
    n0, n1, n2 = [FermionOperator(((i, 1), (i, 0)), 1) for i in range(3)]
    hop = FermionOperator('0^ 1', 0.7) + FermionOperator('1^ 0', 0.7)
    double = FermionOperator('3^ 2^ 1 0', 0.23)
    # a0 a0^ has an exact scalar contraction; n2*hop has a repeated spectator.
    return (FermionOperator((), 0.17) + 1.1 * n0 + 0.3 * n1
            + 1.3 * FermionOperator('0 0^') * n1 + hop + 0.4 * n2 * hop
            + double + hermitian_conjugated(double))


def owner_data(fermion_hamiltonian):
    mapped = jordan_wigner(fermion_hamiltonian)
    coefficients = {key: value for key, value in mapped.terms.items() if key and abs(value) > TOL}
    keys = list(coefficients)
    parents = dense.build_hermitian_fermion_terms(fermion_hamiltonian, TOL)
    descendants = coloring.precompute_fermion_to_pauli_indices(
        parents, coefficients, {key: i for i, key in enumerate(keys)}, TOL,
    )
    return mapped, coefficients, keys, parents, descendants


@pytest.fixture(scope='module')
def molecular_data():
    raw = raw_molecular_fixture()
    cleaned = dense.clean_fermion_operator(raw, TOL)
    return raw, cleaned, owner_data(cleaned)


@pytest.mark.parametrize('order', [1, 2, 4])
@pytest.mark.parametrize('steps', [1, 3])
@pytest.mark.parametrize('parallel_threshold', [1, 100000])
def test_numba_chronological_evolution_matches_independent_dense_formula(order, steps, parallel_threshold):
    keys = [((0, 'Y'), (2, 'X')), ((1, 'Z'),), (), ((0, 'X'), (1, 'Y'))]
    coefficients = dict(zip(keys, [0.4, -0.7, 0.19, 0.3]))
    matrices = [coefficients[key] * pauli_matrix(key, 3) for key in keys]
    rng = np.random.default_rng(197)
    initial = rng.normal(size=8) + 1j * rng.normal(size=8)
    initial /= np.linalg.norm(initial)
    compiled = matrix_free.compile_ordered_terms(keys, coefficients, 3, TOL)
    actual, count = matrix_free.evolve_trotter_state(initial, compiled, order, steps, 0.43, parallel_threshold)
    expected = np.linalg.matrix_power(formula(matrices, order, 0.43 / steps), steps) @ initial
    np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=0)
    assert count == dense.nominal_exponential_count(len(keys), order, steps)


def test_normal_ordering_retains_contractions_and_final_identity(molecular_data):
    raw, cleaned, (mapped, coefficients, _, parents, _) = molecular_data
    expected = fermion_matrix(raw, 4)
    np.testing.assert_allclose(fermion_matrix(cleaned, 4), expected, atol=1e-14, rtol=0)
    np.testing.assert_allclose(qubit_matrix(mapped, 4), expected, atol=1e-14, rtol=0)
    scalar = mapped.terms.get((), 0)
    assert abs(scalar - raw.terms.get((), 0)) > 0.1  # Scalar is not just the raw input constant.
    consolidated = sum((value * pauli_matrix(key, 4) for key, value in coefficients.items()), np.zeros_like(expected))
    np.testing.assert_allclose(consolidated + scalar * np.eye(16), expected, atol=1e-14, rtol=0)
    parent_sum = sum((fermion_matrix(parent.operator, 4) for parent in parents), np.zeros_like(expected))
    np.testing.assert_allclose(parent_sum + cleaned.terms.get((), 0) * np.eye(16), expected, atol=1e-14, rtol=0)


@pytest.mark.parametrize('parent_method', ['signed_ascending', 'magnitude_descending'])
def test_consolidated_owned_blocks_conserve_number_and_reconstruct(molecular_data, parent_method):
    _, _, (_, coefficients, keys, parents, descendants) = molecular_data
    order = baseline.fermionic_term_order_indices(parents, parent_method, TOL)
    buckets, owners, fallback = ablation.build_reference_parent_buckets(order, descendants, len(keys))
    assert not fallback
    assert np.all(owners >= 0)
    number = np.diag([index.bit_count() for index in range(16)])
    for parent, indices in buckets:
        matrices = [coefficients[keys[i]] * pauli_matrix(keys[i], 4) for i in indices]
        block = sum(matrices)
        np.testing.assert_allclose(block @ number, number @ block, atol=1e-14, rtol=0)
        for a in matrices:
            for b in matrices:
                np.testing.assert_allclose(a @ b, b @ a, atol=1e-14, rtol=0)
    flattened = ablation.flatten_reference_buckets(buckets, fallback)
    assert sorted(flattened) == list(range(len(keys)))
    np.testing.assert_allclose(sum(coefficients[keys[i]] * pauli_matrix(keys[i], 4) for i in flattened),
                               sum(value * pauli_matrix(key, 4) for key, value in coefficients.items()),
                               atol=1e-14, rtol=0)


@pytest.mark.parametrize('order', [1, 2, 4])
def test_within_owner_permutation_is_exact_for_every_formula(molecular_data, order):
    _, _, (_, coefficients, keys, parents, descendants) = molecular_data
    parent_order = baseline.fermionic_term_order_indices(parents, 'signed_ascending', TOL)
    buckets, _, fallback = ablation.build_reference_parent_buckets(parent_order, descendants, len(keys))
    original = ablation.flatten_reference_buckets(buckets, fallback)
    changed = ablation.flatten_reference_buckets([(parent, indices[::-1]) for parent, indices in buckets], fallback)
    assert changed != original
    blocks = [sum(coefficients[keys[i]] * pauli_matrix(keys[i], 4) for i in indices) for _, indices in buckets]
    expected = formula(blocks, order, 0.37)
    for indices in (original, changed):
        factors = [(dense.real_coefficient(coefficients[keys[i]], TOL), pauli_matrix(keys[i], 4)) for i in indices]
        actual = {1: dense.first_order_slice, 2: dense.second_order_slice, 4: dense.fourth_order_slice}[order](factors, 0.37, np.eye(16))
        np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=0)


def test_original_parent_product_is_not_the_owned_block_product():
    # Manuscript witness: n0, hopping01, n0*n1. Both sums are the same H.
    n0 = FermionOperator('0^ 0')
    hop = FermionOperator('0^ 1') + FermionOperator('1^ 0')
    density = dense.clean_fermion_operator(n0 * FermionOperator('1^ 1'), TOL)
    clean = dense.clean_fermion_operator(n0 + hop + density, TOL)
    _, coefficients, keys, parents, descendants = owner_data(clean)
    buckets, _, fallback = ablation.build_reference_parent_buckets(range(len(parents)), descendants, len(keys))
    assert not fallback
    owned = [sum(coefficients[keys[i]] * pauli_matrix(keys[i], 2) for i in indices) for _, indices in buckets]
    original = [fermion_matrix(parent.operator, 2) for parent in parents]
    original = [h - np.trace(h) / 4 * np.eye(4) for h in original]
    np.testing.assert_allclose(sum(owned), sum(original), atol=1e-14, rtol=0)
    assert np.linalg.norm(formula(owned, 1, 0.4) - formula(original, 1, 0.4)) > 1e-2


def test_hf_bch_diagnostic_matches_dense_commutators_for_all_five_orders(molecular_data):
    _, cleaned, (_, coefficients, keys, _, _) = molecular_data
    schedules = baseline.build_deterministic_orderings(cleaned, keys, coefficients, 4, TOL)
    graph, _ = dense.build_pauli_noncommutation_graph(keys)
    evaluator = coloring.HFCommutatorEvaluator.build(keys, coefficients, graph, 4, 2, TOL)
    hf = np.zeros(16, dtype=complex)
    hf[12] = 1  # |1100>, independent of the benchmark's basis-index helper.
    matrices = {key: coefficients[key] * pauli_matrix(key, 4) for key in keys}
    for schedule in schedules.values():
        bch = sum((matrices[a] @ matrices[b] - matrices[b] @ matrices[a]
                   for i, a in enumerate(schedule) for b in schedule[i + 1:]), np.zeros((16, 16), complex))
        expected = np.linalg.norm(bch @ hf)
        assert expected > 1e-3
        actual = evaluator.evaluate([keys.index(key) for key in schedule])
        assert actual == pytest.approx(expected, abs=2e-14)


@pytest.mark.xfail(strict=True, reason='Historical normal_ordered uses OpenFermion implicit 1e-8 pruning before requested 1e-12 compression.')
def test_legacy_requested_cutoff_preserves_a_larger_coefficient():
    operator = FermionOperator('0^ 0', 5e-10)
    cleaned = dense.clean_fermion_operator(operator, TOL)
    np.testing.assert_allclose(fermion_matrix(cleaned, 1), fermion_matrix(operator, 1), atol=TOL, rtol=0)


@pytest.mark.xfail(strict=True, reason='Historical fermionic commutator arithmetic drops genuine 1e-10 edges at implicit 1e-8 cutoff.')
def test_legacy_graph_detects_commutator_above_requested_cutoff():
    left = FermionOperator('0^ 0', 1e-5)
    right = FermionOperator('0^ 1', 1e-5) + FermionOperator('1^ 0', 1e-5)
    a, b = fermion_matrix(left, 2), fermion_matrix(right, 2)
    assert np.linalg.norm(a @ b - b @ a, 2) > 10 * TOL
    assert dense.fermionic_terms_noncommute(left, right, TOL)


@pytest.mark.parametrize('coefficient', [5e-10, -7e-11, 2e-12])
def test_explicit_transforms_preserve_sub_openfermion_cutoff_weights(coefficient):
    operator = FermionOperator('0^ 0', coefficient)
    normalized = normal_ordered_molecular(operator, atol=0)
    mapped = jordan_wigner_molecular(operator, atol=0)
    expected = fermion_matrix(operator, 1)
    np.testing.assert_allclose(fermion_matrix(normalized, 1), expected, atol=1e-25, rtol=0)
    np.testing.assert_allclose(qubit_matrix(mapped, 1), expected, atol=1e-25, rtol=0)


def test_explicit_normal_order_contractions_cancellation_and_final_cutoff():
    raw = FermionOperator()
    raw.terms = {((0, 0), (0, 1)): 1.0, ((0, 1), (0, 0)): 1.0, (): -(1 - 5e-10)}
    normalized = normal_ordered_molecular(raw, atol=TOL)
    assert set(normalized.terms) == {()}
    np.testing.assert_allclose(fermion_matrix(normalized, 1), fermion_matrix(raw, 1), atol=1e-16, rtol=0)
    raw.terms = {((0, 1), (1, 0)): 0.7 * TOL, ((1, 0), (0, 1)): -0.7 * TOL}
    assert normal_ordered_molecular(raw, atol=TOL).terms == {((0, 1), (1, 0)): 1.4 * TOL}
    assert normal_ordered_molecular(raw, atol=2 * TOL).terms == {}


def test_explicit_jw_combines_before_cutoff_and_keeps_scalar():
    raw = FermionOperator()
    raw.terms = {((0, 1), (0, 0)): 1.4 * TOL, ((1, 1), (1, 0)): 1.4 * TOL}
    mapped = jordan_wigner_molecular(raw, atol=TOL)
    assert mapped.terms == {(): 1.4 * TOL}  # Two individually sub-cutoff scalar contributions combine.


def test_explicit_transform_preserves_tiny_tensor_entries_and_no_extra_half():
    one = np.array([[5e-10, 0], [0, -3e-10]])
    two = np.zeros((2,) * 4)
    two[0, 1, 1, 0] = 7e-10
    interaction = InteractionOperator(2e-10, one, two)
    raw = fermion_from_interaction(interaction)
    assert raw.terms[((0, 1), (1, 1), (1, 0), (0, 0))] == 7e-10
    mapped = jordan_wigner_molecular(raw, atol=0)
    np.testing.assert_allclose(qubit_matrix(mapped, 2), fermion_matrix(raw, 2), atol=1e-24, rtol=0)


def test_explicit_transform_reconstructs_frozen_molecular_integrals():
    path = next((Path(__file__).parent / 'fixtures' / 'H2S_s1.00').glob('*.tensors.npz'))
    interaction, _ = dense.load_interaction_operator(path)
    # Four-mode restriction of actual supplied integrals keeps this oracle tiny.
    restricted = InteractionOperator(interaction.constant, interaction.one_body_tensor[:4, :4],
                                     interaction.two_body_tensor[:4, :4, :4, :4])
    raw = fermion_from_interaction(restricted)
    normalized = normal_ordered_molecular(raw, atol=0)
    mapped = jordan_wigner_molecular(raw, atol=0)
    expected = fermion_matrix(raw, 4)
    actual = qubit_matrix(mapped, 4)
    np.testing.assert_allclose(fermion_matrix(normalized, 4), expected, atol=2e-13, rtol=0)
    # This fixture includes an approximately -400 Hartree scalar. Allow a
    # small machine-relative summation error while tiny-coefficient tests above
    # independently use absolute tolerances down to 1e-25.
    np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=2e-15)
    np.testing.assert_allclose(actual, actual.conj().T, atol=2e-13, rtol=0)
    number = np.diag([index.bit_count() for index in range(16)])
    np.testing.assert_allclose(actual @ number, number @ actual, atol=2e-13, rtol=0)


@pytest.mark.parametrize('transform', [normal_ordered_molecular, jordan_wigner_molecular])
def test_explicit_transforms_reject_unsupported_rank_and_invalid_cutoff(transform):
    with pytest.raises(ValueError, match='at most four'):
        transform(FermionOperator('0^ 1^ 2^ 2 1 0'), atol=TOL)
    with pytest.raises(ValueError, match='nonnegative'):
        transform(FermionOperator('0^ 0'), atol=-1)
    with pytest.raises(ValueError, match='finite'):
        transform(FermionOperator('0^ 0', np.inf), atol=TOL)


def test_explicit_ownership_preserves_tiny_molecular_hamiltonian_and_blocks():
    from qhat.analysis.fermionic_ownership import build_molecular_ownership

    raw = raw_molecular_fixture()
    raw.terms = {key: value * 1e-9 for key, value in raw.terms.items()}
    schedule = build_molecular_ownership(raw, atol=1e-15)
    expected = fermion_matrix(raw, 4)
    actual = sum((value * pauli_matrix(key, 4) for key, value in schedule.final_coefficients.items()), np.zeros((16, 16), complex))
    np.testing.assert_allclose(actual + schedule.scalar * np.eye(16), expected, atol=3e-24, rtol=0)
    np.testing.assert_allclose(fermion_matrix(schedule.hamiltonian, 4), expected, atol=3e-24, rtol=0)
    assert not schedule.fallback_pauli_keys
    assert schedule.diagnostics['noncommuting_owned_block_count'] == 0
    assert schedule.diagnostics['max_owned_block_number_commutator_l1'] < 1e-24
    number = np.diag([index.bit_count() for index in range(16)])
    for block in schedule.owned_blocks:
        h = sum(schedule.final_coefficients[key] * pauli_matrix(key, 4) for key in block.pauli_keys)
        np.testing.assert_allclose(h @ number, number @ h, atol=3e-24, rtol=0)
        for key in block.pauli_keys:
            p = pauli_matrix(key, 4)
            # Independently recover each parent contribution by Hilbert--Schmidt
            # projection of the direct occupation-basis fermion matrix.
            first = next(index for index in schedule.parent_order
                         if abs(np.trace(p @ fermion_matrix(schedule.parents[index].operator, 4)) / 16) > 1e-24)
            assert first == block.parent_index


def test_explicit_ownership_reports_small_hermitian_mismatch_without_averaging():
    from qhat.analysis.fermionic_ownership import build_molecular_ownership

    raw = FermionOperator('0^ 1', 0.2) + FermionOperator('1^ 0', 0.2 + 5e-14)
    schedule = build_molecular_ownership(raw, atol=TOL)
    assert schedule.hamiltonian.terms == raw.terms
    assert schedule.parents[0].operator.terms == raw.terms
    assert schedule.diagnostics['parent_hermiticity_residual_l1'] > 0
    assert schedule.diagnostics['discarded_pauli_coefficient_l1'] > 0


@pytest.mark.parametrize('bad, message', [
    (FermionOperator('0^ 1', 0.2) + FermionOperator('1^ 0', 0.21), 'Hermitian-adjoint'),
    (FermionOperator('0^ 1', 0.2j) + FermionOperator('1^ 0', -0.2j), 'real canonical'),
    (FermionOperator('0^') + FermionOperator('0'), 'number-conserving'),
])
def test_explicit_ownership_rejects_inputs_outside_stated_scope(bad, message):
    from qhat.analysis.fermionic_ownership import build_molecular_ownership

    with pytest.raises(ValueError, match=message):
        build_molecular_ownership(bad, atol=TOL)


def test_particle_number_commutator_diagnostic_has_correct_sign_and_cancellation():
    from qhat.analysis.fermionic_ownership import particle_number_commutator_l1

    # X alone changes number; XX+YY is a conserving hopping bundle.
    assert particle_number_commutator_l1({((0, 'X'),): 0.3}) == pytest.approx(0.3)
    assert particle_number_commutator_l1({((0, 'X'), (1, 'X')): 0.3,
                                         ((0, 'Y'), (1, 'Y')): 0.3}) == 0
