"""Independent CAR density matrices, sector states and charged constructor paths."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
import state_preparation as sp
import transfer_common as k


def annihilation(q, n):
    factors = []
    for p in reversed(range(n)):
        factors.append(np.array([[0, 1], [0, 0]]) if p == q else
                       np.diag([1, -1]) if p < q else np.eye(2))
    result = np.array([[1.]])
    for factor in factors:
        result = np.kron(result, factor)
    return result


def literal(word, n):
    letters = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
               'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    labels = dict(word)
    result = np.array([[1.]])
    for q in reversed(range(n)):
        result = np.kron(result, letters[labels.get(q, 'I')])
    return result


def toy():
    terms = [([[(0, 'Z')], '3']), ([[(1, 'Z')], '-2']),
             ([[(2, 'Z')], '1']), ([[(3, 'Z')], '-1']),
             ([[(0, 'Z'), (1, 'Z')], '5']),
             ([[(0, 'X'), (2, 'X')], '2']), ([[(0, 'Y'), (2, 'Y')], '2']),
             ([[(1, 'X'), (3, 'X')], '1']), ([[(1, 'Y'), (3, 'Y')], '1'])]
    return dict(case='reference_toy', n=4, electrons=2, denominator='8', coefficients=terms,
                blocks=[[t] for t in terms[:5]] + [terms[5:7], terms[7:9]])


def test_RDM_matches_independent_full_JW_matrices():
    psi = np.array([1 + .3j, -.2, .7j, .4, -.5j, .1, .3j, .6] * 2)
    psi /= np.linalg.norm(psi)
    gamma, natural, defect, label = sp.occupations(psi, 4)
    operators = [annihilation(q, 4) for q in range(4)]
    expected = np.array([[np.vdot(psi, a.conj().T @ b @ psi) for b in operators] for a in operators])
    assert np.max(abs(gamma - expected)) < 1e-14
    assert natural == pytest.approx(np.linalg.eigvalsh(expected), abs=1e-14)
    assert defect == pytest.approx(np.trace(expected - expected @ expected).real, abs=1e-14)


def test_participation_is_not_a_correlation_label():
    slater = np.zeros(16, complex)
    slater[3] = slater[6] = 1 / np.sqrt(2)
    _, _, defect, label = sp.occupations(slater, 4)
    assert defect < 1e-12 and label == 'Slater-like'
    correlated = np.zeros(16, complex)
    correlated[3] = correlated[12] = 1 / np.sqrt(2)
    _, _, defect, label = sp.occupations(correlated, 4)
    assert defect == pytest.approx(1.) and label == 'correlated'


def test_all_five_state_constructors_match_literal_eigenspaces_and_phases(tmp_path):
    target = toy()
    path = tmp_path / 'model.json'
    k.write(path, target)
    H = sum(int(c) / 8 * literal(w, 4) for w, c in target['coefficients'])
    sector = [3, 6, 9, 12]
    values, vectors = np.linalg.eigh(H[np.ix_(sector, sector)])
    produced = {}
    for label in k.STATES:
        psi, record = sp.construct(path, k.sha(path), label)
        produced[label] = psi
        assert record['metadata']['sector'] == sector
        assert abs(np.linalg.norm(psi) - 1) < 1e-14
        assert record['CPU_seconds'] >= sum(record['phases'].values()) - 1e-12
        if label == 'HF':
            assert psi[3] == 1
            assert record['phases']['eigensystem_CPU_seconds'] == 0
        elif label in ('ground', 'excited'):
            at = 0 if label == 'ground' else 1
            assert abs(np.vdot(vectors[:, at], psi[sector])) == pytest.approx(1., abs=1e-12)
            assert np.linalg.norm(H @ psi - values[at] * psi) < 2e-12
    g, e = produced['ground'], produced['excited']
    assert produced['coherent_real'] == pytest.approx((g + e) / np.sqrt(2), abs=1e-14)
    assert produced['coherent_imag'] == pytest.approx((g + 1j * e) / np.sqrt(2), abs=1e-14)
    assert np.linalg.norm(produced['coherent_imag'].imag) > .1


def test_spin_leaking_unit_is_refused_before_sector_projection(tmp_path):
    target = toy()
    target['blocks'][5][1][1] = '1'
    target['coefficients'][6][1] = '1'
    path = tmp_path / 'model.json'
    k.write(path, target)
    with pytest.raises(ValueError, match='Exact complete-unit spin-sector conservation'):
        sp.construct(path, k.sha(path), 'ground')


def test_oversized_register_refused_before_dense_allocation(tmp_path):
    target = toy()
    target['n'] = 30
    path = tmp_path / 'model.json'
    k.write(path, target)
    with pytest.raises(ValueError, match='before allocation'):
        sp.construct(path, k.sha(path), 'ground')


def test_independent_tensor_CAR_index_convention():
    from verify_references import tensor_hamiltonian
    one = np.diag([2., 5.])
    two = np.zeros((2, 2, 2, 2))
    two[0, 1, 1, 0] = 3.
    two[0, 1, 0, 1] = 2.
    H, _ = tensor_hamiltonian(one, two, .7)
    assert H == pytest.approx(np.diag([.7, 2.7, 5.7, 8.7]))


def test_geometry_identity_ignores_rigid_motion_and_atom_order():
    from geometry_identity import signature, same
    first = [['N', [0., 0., 0.]], ['H', [1., 0., 0.]],
             ['H', [0., 2., 0.]], ['H', [0., 0., 3.]]]
    rotated = [[atom, [2. - y, 3. + x, 4. + z]] for atom, (x, y, z) in reversed(first)]
    assert same(signature(first), signature(rotated))
    rotated[0][1][0] += .001
    assert not same(signature(first), signature(rotated))


@pytest.mark.parametrize('partition', ['conserving', 'pauli'])
def test_first_order_adapter_reduces_to_verified_HF_objective(partition):
    import training
    target = toy()
    psi = np.eye(16, dtype=complex)[:, 3]
    current = training.FirstOrderContext(target, psi, partition)
    old = training.prior.FirstOrderContext(target, k.TRAIN_TIME, partition)
    for order in [list(range(current.m)), list(reversed(range(current.m)))]:
        actual = current.score(order, k.TRAIN_TIME, 270)
        expected = old.score(order, k.TRAIN_TIME, 270)
        assert actual['steps'] == expected['steps']
        assert actual['JD1'] == pytest.approx(expected['JD1'], rel=1e-12, abs=1e-18)


def test_first_order_coherent_adapter_matches_literal_integral():
    import training
    from scipy.integrate import quad_vec
    target = toy()
    psi = np.zeros(16, complex)
    psi[[3, 6, 9, 12]] = [1, 1j, -.2, .3j]
    psi /= np.linalg.norm(psi)
    context = training.FirstOrderContext(target, psi, 'conserving')
    matrices = [sum(int(c) / 8 * literal(w, 4) for w, c in block) for block in target['blocks']]
    prefix = np.zeros((16, 16), complex)
    C = prefix.copy()
    for G in matrices:
        C += -.5j * (prefix @ G - G @ prefix)
        prefix += G
    diagonal = np.diag(prefix).real
    def action(t):
        return np.exp(1j * t * diagonal) * (C @ (np.exp(-1j * t * diagonal) * psi))
    Qpsi, _ = quad_vec(action, 0, k.TRAIN_TIME, epsabs=1e-13)
    Qpsi -= psi * np.vdot(psi, Qpsi)
    assert context.score(list(range(context.m)), k.TRAIN_TIME, 270)['JD1'] == pytest.approx(
        np.vdot(Qpsi, Qpsi).real, rel=1e-12, abs=1e-18)


def test_cubic_training_charges_warm_constructor_and_preparation_once(tmp_path):
    import training
    target = toy()
    path = tmp_path / 'model.json'
    k.write(path, target)
    psi = np.eye(16, dtype=complex)[:, 3]
    result = training.cubic(path, k.sha(path), psi, reference_CPU_charge=.125)
    assert result['warm_start']['reference_preparation_CPU_charge'] == 0
    assert result['standalone_reference_and_ordering_CPU_charge'] == pytest.approx(
        result['actual_local_CPU_seconds'] + .125)
    assert result['actual_local_CPU_seconds'] >= result['warm_start']['actual_local_CPU_seconds']
    assert result['record']['order'] == result['result']['record']['order']
    assert result['actual_local_CPU_seconds'] >= sum(result['phases'].values()) - 1e-12
