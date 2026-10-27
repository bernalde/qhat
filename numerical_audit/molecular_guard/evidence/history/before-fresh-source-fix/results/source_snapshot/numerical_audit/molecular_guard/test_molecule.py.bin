"""Independent finite polynomial, Pauli-action, BCH and propagation checks."""
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import pytest
from intervals import Real, Complex
import exact as ex
import propagation as p


def toy_model():
    blocks = [[[[[0,'Z']],'3']], [[[[2,'Z']],'1']],
              [[[[0,'X'],[1,'Z'],[2,'X']],'2'], [[[0,'Y'],[1,'Z'],[2,'Y']],'2']]]
    return dict(n=4,electrons=2,denominator='8',blocks=blocks,
                coefficients=[term for block in blocks for term in block])


def test_complete_actions_match_literal_local_pauli_multiplication():
    # Independent local matrices with integer Gaussian phases, no mask formula.
    local = {'X':(((0,0),(1,0)),((1,0),(0,0))),
             'Y':(((0,0),(0,-1)),((0,1),(0,0))),
             'Z':(((1,0),(0,0)),((0,0),(-1,0)))}
    model = toy_model()
    data = ex.from_model(model)
    for block,matrix in zip(model['blocks'],data['matrices']):
        for i,w in enumerate(data['sector']):
            for j,v in enumerate(data['sector']):
                real = imag = 0
                for word,c in block:
                    a,b = 1,0
                    occupied = {q:letter for q,letter in word}
                    for q in range(4):
                        r,s = (local[occupied[q]][(w>>q)&1][(v>>q)&1] if q in occupied
                               else (int(((w>>q)&1)==((v>>q)&1)),0))
                        a,b = a*r-b*s,a*s+b*r
                    real += int(c)*a
                    imag += int(c)*b
                assert imag == 0
                assert matrix[i][j] == F(real,8)
    assert data['sector'] == [3,6,9,12]
    assert data['matrices'][2][1][0] == F(-1,2)


def test_incomplete_flip_class_is_rejected_before_sector_projection():
    model = toy_model()
    model['blocks'][2] = model['blocks'][2][:1]
    model['coefficients'] = [term for block in model['blocks'] for term in block]
    with pytest.raises(ValueError,match='full-space spin-sector'):
        ex.from_model(model)


def test_dimensions_are_refused_before_large_allocation():
    model = toy_model()
    model['n'] = 10**9
    with pytest.raises(ValueError,match='before allocation'):
        ex.from_model(model)


def test_physical_bch_sign_and_exact_epsilon_decomposition():
    d = [[F(1),F(0)],[F(0),F(-1)]]
    v = [[F(0),F(2)],[F(2),F(0)]]
    data = dict(matrices=[d,v],diag=[0],sector=[0,1])
    a,b = ex.bch(data,[0,1])
    assert a == [[0,-2],[2,0]] and b == ex.zero(2)
    assert ex.bch(data,[1,0])[0] == ex.times(a,-1)
    # Add a noncommuting second off-diagonal unit on a third determinant.
    ds = [[F(1),F(0),F(0)],[F(0),F(-2),F(0)],[F(0),F(0),F(3)]]
    va = [[F(0),F(2),F(0)],[F(2),F(0),F(0)],[F(0),F(0),F(0)]]
    vb = [[F(0),F(0),F(0)],[F(0),F(0),F(3)],[F(0),F(3),F(0)]]
    k1,k2 = ex.bch(dict(matrices=[ds,va,vb],diag=[0],sector=[0,1,2]),[0,1,2])
    assert k2[0][2] == -3 and k2[2][0] == 3
    epsilon = F(2,7)
    full = ex.times(ex.add(ex.add(ex.commutator(ds,ex.times(va,epsilon)),
                                 ex.commutator(ds,ex.times(vb,epsilon))),
                          ex.commutator(ex.times(va,epsilon),ex.times(vb,epsilon))),F(-1,2))
    assert full == ex.add(ex.times(k1,epsilon),ex.times(k2,epsilon**2))


def test_matrix_polynomial_has_hand_computable_coefficients():
    theta = F(1,8)
    h = [[F(0),F(1,4)],[F(1,4),F(0)]]
    a,tail,details = p.exp_matrix(h,F(1,2),terms=4,bits=64)
    assert details['scaling'] == 1 and tail >= F(1,8)
    cosine,sine = 1-theta**2/2,-theta+theta**3/6
    assert a[0][0].real.contains(cosine)
    assert a[1][1].real.contains(cosine)
    assert a[0][1].imag.contains(sine)
    assert a[1][0].imag.contains(sine)
    assert a[0][1].real.contains(0) and a[0][0].imag.contains(0)
    structured,structured_tail = p.exp_unit(h,F(1,2),terms=4,bits=64)
    assert structured[0][0].real.contains(cosine)
    assert structured[0][1].imag.contains(sine)
    assert structured_tail == tail


def test_matrix_and_structured_routes_agree_with_independent_scipy():
    import numpy as np
    from scipy.linalg import expm
    data = ex.from_model(toy_model())
    for h in data['matrices']:
        generic,gtail,_ = p.exp_matrix(h,F(3,5))
        structured,stail = p.exp_unit(h,F(3,5))
        truth = expm(-.6j*np.array(h,dtype=float))
        for result,tail in ((generic,gtail),(structured,stail)):
            for i in range(len(h)):
                for j in range(len(h)):
                    re = sum(result[i][j].real.bounds())/2
                    im = sum(result[i][j].imag.bounds())/2
                    assert abs(complex(float(re),float(im))-truth[i,j]) < 2e-15


def test_written_product_is_applied_in_reverse_order():
    import numpy as np
    from scipy.linalg import expm
    d = [[F(1),F(0)],[F(0),F(-1)]]
    v = [[F(0),F(1)],[F(1),F(0)]]
    units = [p.exp_unit(h,F(3,5))[0] for h in (d,v)]
    initial = [Complex.exact(1),Complex.exact()]
    state = initial
    for unit in reversed(units):
        state = p.matvec(unit,state)
    actual = np.array([complex(float(sum(z.real.bounds())/2),float(sum(z.imag.bounds())/2)) for z in state])
    ud,uv = [expm(-.6j*np.array(h,dtype=float)) for h in (d,v)]
    assert np.max(abs(actual-(ud@uv)[:,0])) < 1e-15
    assert np.max(abs(actual-(uv@ud)[:,0])) > .1


def test_zero_gap_transport_and_zero_denominator_guard():
    data = dict(sector=[0,1],hf_index=0,matrices=[ex.zero(2),[[F(0),F(1)],[F(1),F(0)]]],diag=[0])
    k1 = [[F(0),F(-2)],[F(2),F(0)]]
    a,c,details = p.leading_response(data,k1,ex.zero(2),F(3,5))
    assert a.contains(F(6,5))
    assert details['components'][0]['gap'] == '0'
    # No division by a gap in response norms is needed, including a tied slope.
    guard = p.response_guard(Real.exact(0),F(1),F(1,4))
    assert guard.lo == 0 and guard.contains(F(1,16))


def test_analytic_infidelity_error_and_unitarity_intersection():
    # True states both equal |0>; polynomial state carries a deliberate norm error.
    reference = [Complex.exact(1),Complex.exact()]
    state = [Complex.exact(F(9,8)),Complex.exact()]
    error,details = p.infidelity(reference,F(0),state,F(1,8))
    assert error.contains(0) and error.lo == 0
    assert F(details['analytic_infidelity_radius']) >= F(1,8)*(2+F(1,8))
    assert p.combine(F(1,8),F(1,16),64) >= (1+F(1,8))*(1+F(1,16))-1


def test_floating_point_coefficients_are_refused():
    model = toy_model()
    model['blocks'][0][0][1] = 3.1
    with pytest.raises(ValueError,match='Exact integer'):
        ex.from_model(model)


def test_interval_rankings_require_strict_separation():
    from study import strict_winner
    assert strict_winner(dict(current=Real(1,2,16),diagonal_prefix=Real(3,4,16))) == 'current'
    assert strict_winner(dict(current=Real(3,4,16),diagonal_prefix=Real(1,2,16))) == 'diagonal_prefix'
    assert strict_winner(dict(current=Real(1,2,16),diagonal_prefix=Real(2,3,16))) == 'unresolved'


def test_writes_preserve_completed_and_interrupted_evidence(tmp_path):
    from study import write
    target = tmp_path/'evidence.json'
    write(target,dict(value=1))
    with pytest.raises(ValueError,match='existing attempt'):
        write(target,dict(value=2))
    assert target.read_text().find('1') >= 0
    other = tmp_path/'interrupted.json'
    other.with_suffix('.json.partial').write_text('retained')
    with pytest.raises(ValueError,match='interrupted write'):
        write(other,dict(value=2))
    assert other.with_suffix('.json.partial').read_text() == 'retained'
