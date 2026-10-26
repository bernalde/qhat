"""Independent full signed operator and coherent integral checks."""
import itertools
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
from scipy.integrate import quad_vec
from cubic_response import Context, Interrupted, prior


def model():
    terms = [([[0, 'Z']], '3'), ([[1, 'Z']], '-2'), ([[0, 'X']], '2'),
             ([[0, 'X'], [1, 'Z']], '1'), ([[1, 'X']], '1')]
    return dict(n=2, electrons=1, denominator='7', coefficients=terms,
                blocks=[[term] for term in terms])


def literal(word):
    pool = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]), 'Z': np.diag([1, -1])}
    labels = dict(word)
    return np.kron(pool[labels.get(1, 'I')], pool[labels.get(0, 'I')])


def reference_operator(target, order):
    matrices = [sum(int(c) / 7 * literal(w) for w, c in block) for block in target['blocks']]
    rest = np.zeros((4, 4))
    defect = np.zeros_like(rest)
    comm = lambda a, b: a @ b - b @ a
    for index in reversed(order):
        outer = matrices[index]
        defect += -comm(outer, comm(outer, rest)) / 24 + comm(rest, comm(rest, outer)) / 12
        rest += outer
    return defect


def state():
    psi = np.array([1, 1j, -.2, .3j])
    return psi / np.linalg.norm(psi)


def test_all_signed_columns_match_independent_nested_commutators():
    target = model()
    context = Context(target, state())
    for order in itertools.permutations(range(5)):
        assert np.max(abs(context.cubic_columns(order) - reference_operator(target, order))) < 2e-14


@pytest.mark.parametrize('phase', [1, 1j, -1, -1j])
def test_coherent_transport_matches_actual_quadrature(phase):
    target = model()
    psi = np.array([1, phase, -.2, .3j], complex)
    psi /= np.linalg.norm(psi)
    context = Context(target, psi)
    order = [2, 0, 4, 1, 3]
    defect = reference_operator(target, order)
    energies = sum(int(c) / 7 * literal(w) for w, c in target['coefficients'] if all(p == 'Z' for q, p in w)).diagonal()
    def action(t):
        return np.exp(1j * t * energies) * (defect @ (np.exp(-1j * t * energies) * psi))
    transported, _ = quad_vec(action, 0, .83, epsabs=1e-13)
    transported -= psi * np.vdot(psi, transported)
    expected = float(np.vdot(transported, transported).real)
    assert context.score(order, .83, 150)['JD2'] == pytest.approx(expected, rel=1e-12, abs=1e-18)


def test_coherent_variance_is_not_incoherent_mixture():
    target = model()
    psi = state()
    order = [2, 0, 4, 1, 3]
    coherent = Context(target, psi).score(order, .83, 150)['JD2']
    mixture = sum(abs(psi[i])**2 * Context(target, np.eye(4, dtype=complex)[:, i]).score(order, .83, 150)['JD2']
                  for i in range(4))
    assert abs(coherent - mixture) > 1e-7


def test_HF_reduction_matches_legacy_full_signed_column_and_score():
    target = model()
    psi = np.eye(4, dtype=complex)[:, 1]
    general = Context(target, psi)
    old = prior.Context(target)
    for order in itertools.permutations(range(5)):
        current = general.score(order, .83, 150)
        previous = old.score(order, .83, 150)
        assert current['steps'] == previous['steps']
        assert np.max(abs(np.array(current['cubic_columns'])[:, 0] - previous['cubic_coefficients'])) < 2e-14
        assert current['JD2'] == pytest.approx(previous['JD2'], rel=1e-12, abs=1e-18)
        assert current['score'] == pytest.approx(previous['score'], rel=1e-12, abs=1e-22)


def test_interruption_preserves_completed_stage_count_and_refuses_oversize():
    context = Context(model(), state())
    with pytest.raises(Interrupted) as error:
        context.cubic_columns([0, 1, 2, 3, 4], deadline=0., clock=lambda: 1.)
    assert error.value.completed_stages == 0
    target = model()
    target['n'] = 30
    with pytest.raises(ValueError, match='before allocation'):
        Context(target, state())
    target = model()
    target['n'] = 6
    psi = np.ones(64) / 8
    with pytest.raises(ValueError, match='support bound'):
        Context(target, psi)
