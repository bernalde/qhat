"""Hand-derived checks of the verifier, including deliberate narrow certificates."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import pytest
import independent as ind


def test_fraction_interval_extrema_and_sqrt():
    a = ind.Q(F(-3,7),F(5,11))
    b = ind.Q(F(-13,17),F(19,23))
    out = a*b
    for x in (a.lo,a.hi):
        for y in (b.lo,b.hi):
            assert out.lo <= x*y <= out.hi
    root = ind.Q.point(F(2,3)).sqrt()
    assert root.lo**2 <= F(2,3) <= root.hi**2
    assert a.square().lo == 0


def test_generic_matrix_exponential_with_independent_decimal_series():
    value,tail = ind.generic_exponential([[F(0),F(1,4)],[F(1,4),F(0)]],F(1,2))
    with localcontext() as context:
        context.prec = 110
        x = Decimal(1)/8
        cosine = ct = Decimal(1)
        sine = st = x
        for k in range(1,150):
            ct *= -x*x/((2*k-1)*(2*k))
            st *= -x*x/((2*k)*(2*k+1))
            cosine += ct
            sine += st
    real = value[0][0][0].expanded(tail)
    imag = value[0][1][1].expanded(tail)
    assert real.lo <= F(cosine) <= real.hi
    assert imag.lo <= -F(sine) <= imag.hi
    assert tail < F(1,1 << 150)


def test_literal_full_space_action_and_conservation_refusal():
    blocks = [[[[[0,'Z']],'3']], [[[[2,'Z']],'1']],
              [[[[0,'X'],[1,'Z'],[2,'X']],'2'], [[[0,'Y'],[1,'Z'],[2,'Y']],'2']]]
    model = dict(n=6,electrons=2,denominator='8',blocks=blocks,
                 coefficients=[term for block in blocks for term in block])
    matrices,sector,hf,diagonal = ind.literal_matrices(model,F(1))
    assert len(sector) == 9 and sector[hf] == 3 and diagonal == [0,1]
    assert matrices[2][sector.index(6)][hf] == F(-1,2)
    model['blocks'][2] = model['blocks'][2][:1]
    model['coefficients'] = [term for block in model['blocks'] for term in block]
    with pytest.raises(ValueError,match='Full-space sector'):
        ind.literal_matrices(model,F(1))


def test_independent_bch_has_physical_sign():
    a = [[F(1),F(0)],[F(0),F(-1)]]
    b = [[F(0),F(2)],[F(2),F(0)]]
    k1,k2 = ind.independent_bch([a,b],[0,1],[0])
    assert k1 == [[0,-2],[2,0]] and k2 == [[0,0],[0,0]]


def test_independent_enclosure_must_fit_inside_certificate():
    raw = dict(lo=str(1 << 191),hi=str(1 << 192),bits=192)
    ind.contained(ind.Q(F(2,3),F(3,4)),raw)
    with pytest.raises(ValueError,match='outside producer certificate'):
        ind.contained(ind.Q(F(1,3),F(3,4)),raw)
    collapsed = dict(lo=str(1 << 191),hi=str(1 << 191),bits=192)
    with pytest.raises(ValueError,match='outside producer certificate'):
        ind.contained(ind.Q(F(1,2),F(1,2)+F(1,1 << 160)),collapsed)


def test_independent_overlap_error_includes_analytic_allowance():
    reference = [ind.zpoint(1),ind.zpoint()]
    state = [ind.zpoint(F(9,8)),ind.zpoint()]
    error = ind.squared_error(reference,F(0),state,F(1,8))
    assert error.lo == 0 and error.hi >= 0
