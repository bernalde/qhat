"""Non-tautological controls for low-energy selection and full-state evolution."""
import json

import numpy as np
import pytest
from scipy.linalg import expm

from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis.state_metrics import state_error_metrics
from qhat.numerical_audit import low_energy_state_sensitivity as experiment


def test_degenerate_basis_is_invariant_to_solver_unitary_rotation():
    rng = np.random.default_rng(2903)
    basis, _ = np.linalg.qr(rng.normal(size=(8,3)) + 1j*rng.normal(size=(8,3)))
    rotation, _ = np.linalg.qr(rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)))
    a = experiment.canonical_subspace_basis(basis)
    b = experiment.canonical_subspace_basis(basis @ rotation)
    np.testing.assert_allclose(a,b,atol=2e-14)
    np.testing.assert_allclose(a.conj().T@a,np.eye(3),atol=2e-14)
    np.testing.assert_allclose(a@a.conj().T,basis@basis.conj().T,atol=2e-14)


def test_ground_next_degeneracy_is_retained_and_coherent_phase_is_fixed():
    matrix = np.diag([3.,-2.,-2.,7.])
    states,info,_,_ = experiment.low_energy_states(matrix)
    np.testing.assert_allclose(states['sector_ground'],[0,1,0,0],atol=1e-14)
    np.testing.assert_allclose(states['sector_next'],[0,0,1,0],atol=1e-14)
    np.testing.assert_allclose(states['equal_coherent'],np.array([0,1,1,0])/np.sqrt(2),atol=1e-14)
    assert info['selected_vectors_same_numerical_eigenspace']
    assert info['first_gap_hartree'] == 0


def test_cluster_uses_total_span_not_transitive_adjacent_gaps():
    matrix = np.diag([0.,.75e-10,1.5e-10,1.])
    _,info,_,_ = experiment.low_energy_states(matrix)
    assert [(c['start'],c['stop']) for c in info['numerical_clusters']] == [(0,2),(2,3),(3,4)]


def test_phase_is_invariant_to_global_complex_phase():
    vector = np.array([.5j,-.5,.5,-.5j])
    canonical = experiment.canonical_phase(vector)
    np.testing.assert_allclose(canonical,experiment.canonical_phase(np.exp(.74j)*vector),atol=1e-14)
    assert canonical[0].real > 0
    assert abs(canonical[0].imag) < 1e-14


def test_corrected_design_hashes_reject_coefficient_or_schedule_mutation():
    path = experiment.SOURCE_RUN/'CH2O_s1.00.design.json'
    design = json.loads(path.read_text())
    coeffs,schedules,_ = experiment.decode_frozen_design(design)
    assert len(coeffs) == design['pauli_terms']
    assert all(set(order) == set(coeffs) for order in schedules.values())
    mutated = json.loads(path.read_text())
    mutated['identity_free_coefficients'][0][1] = float(1.).hex()
    with pytest.raises(RuntimeError,match='Hamiltonian hash'):
        experiment.decode_frozen_design(mutated)
    mutated = json.loads(path.read_text())
    order = mutated['pauli_orders']['block_signed']
    order[0],order[1] = order[1],order[0]
    with pytest.raises(RuntimeError,match='schedule hash'):
        experiment.decode_frozen_design(mutated)


def test_full_kernel_matches_literal_product_and_retains_sector_leakage():
    # XX+YY conserves number; splitting its pieces around Z leaks at finite dt.
    x = np.array([[0,1],[1,0]],complex)
    y = np.array([[0,-1j],[1j,0]],complex)
    z = np.diag([1.,-1.])
    matrices = [np.kron(x,x),.37*np.kron(z,np.eye(2)),np.kron(y,y)]
    coeffs = {((0,'X'),(1,'X')):1.,((0,'Z'),):.37,((0,'Y'),(1,'Y')):1.}
    state = np.array([1,0,0,0],complex)
    terms = evolution.compile_ordered_terms(list(coeffs),coeffs,2,1e-12)
    actual,count = evolution.evolve_trotter_state(state,terms,1,3,.71,2**30)
    expected = state.copy()
    for _ in range(3):
        for matrix in matrices:
            expected = expm(-1j*.71/3*matrix)@expected
    np.testing.assert_allclose(actual,expected,atol=2e-14)
    assert count == 9
    assert abs(actual[3])**2 > 1e-5
    exact = expm(-1j*.71*sum(matrices))@state
    metrics = state_error_metrics(exact,actual)
    normalized = actual/np.linalg.norm(actual)
    assert metrics['one_minus_overlap'] > 1e-5
    assert abs(metrics['one_minus_overlap']-(1-abs(np.vdot(exact,normalized)))) < 1e-14
