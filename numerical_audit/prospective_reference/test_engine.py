"""Independent small-matrix checks before prospective molecular evaluation."""
import itertools
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.linalg import eigh,expm
from scipy.integrate import quad_vec
import engine as e
from campaign import current_placement,target_matrix
from verify import literal


def units():
    return [{((0,'Z'),):2},{((2,'Z'),):3},{((0,'Z'),(1,'Z')):1},
            {((0,'X'),(2,'X')):1,((0,'Y'),(2,'Y')):1},
            {((1,'X'),(3,'X')):2,((1,'Y'),(3,'Y')):2}]


def dense(block):
    return sum(c*literal(w,4)[0].toarray() for w,c in block.items())


def full_response(C,A,v):
    w,U=eigh(A);gaps=w[:,None]-w[None,:]
    Q=U@((np.exp(.5j*gaps)*np.sinc(gaps/(2*np.pi)))*(U.conj().T@C@U))@U.conj().T
    col=Q[:,v].copy();col[v]=0
    return col


def test_two_level_response_against_direct_quadrature():
    rng=np.random.default_rng(82716);v=3;u=12
    diagonal={((0,'Z'),):2,((2,'Z'),):3,((0,'Z'),(1,'Z')):1}
    raw=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16));C=(raw+raw.conj().T)/2
    ref=e.Reference(diagonal,v,1.,1,u,2)
    actual=ref.response([{w:C[w,col] for w in range(16)} for col in (v,u)])
    A=dense(diagonal).astype(complex);A[v,u]=A[u,v]=2
    direct,err=quad_vec(lambda t:expm(1j*t*A)@C@expm(-1j*t*A)[:,v],0,1,epsabs=1e-11)
    direct[v]=0
    observed=np.array([actual.get(w,0) for w in range(16)])
    assert np.linalg.norm(observed-direct)<2e-11
    zero=e.Reference(diagonal,v,1.,1,u,0).response([{w:C[w,v] for w in range(16)}])
    simple=e.Reference(diagonal,v,1.,1).response([{w:C[w,v] for w in range(16)}])
    assert zero==simple


@pytest.mark.parametrize('corrected',[False,True])
def test_every_small_permutation_against_literal_commutators(corrected):
    blocks=units();v=3
    diagonal={w:c for b in blocks[:3] for w,c in b.items()}
    partner,coupling=e.choose_partner(blocks,v)
    ref=e.Reference(diagonal,v,1.,1,partner if corrected else None,coupling)
    features=e.Features(blocks,ref);G=[dense(b) for b in blocks]
    A=dense(diagonal).astype(complex)
    if corrected:A[v,partner]=A[partner,v]=coupling
    for order in itertools.permutations(range(len(blocks))):
        C=sum((-.5j*(G[i]@G[j]-G[j]@G[i]) for at,i in enumerate(order) for j in order[at+1:]),start=np.zeros_like(A))
        response=full_response(C,A,v);expected=float(np.vdot(response,response).real)
        assert features.score(features.amplitude(order))==pytest.approx(expected,rel=2e-12,abs=2e-12)


def test_exact_current_rule_and_corrected_placement_keep_target_and_off_order():
    blocks=units();diag=[1,0,2];off=[4,3];v=3
    D={w:c for i in diag for w,c in blocks[i].items()}
    features=e.Features(blocks,e.Reference(D,v,1.,1))
    result=current_placement(blocks,diag,off,v,1)
    assert [i for i in result['order'] if i in off]==off
    assert result['score']==pytest.approx(features.score(features.amplitude(result['order'])))
    assert result['evaluations']==(len(diag)+1)*(len(off)+1)
    u,h=e.choose_partner(blocks,v);corrected=e.Features(blocks,e.Reference(D,v,1.,1,u,h))
    placed=e.placement(blocks,diag,off,corrected)
    assert [i for i in placed['order'] if i in off]==off
    assert sorted(placed['order'])==list(range(len(blocks)))
    assert placed['score']<=placed['initial_score']+1e-12
    assert placed['evaluations']==result['evaluations']
    assert all(x['after']<=x['before'] for x in placed['history'])


def test_insertion_budget_crossing_and_nonincrease():
    blocks=units();features=e.Features(blocks,e.Reference({((0,'Z'),):1},3,1.,1))
    starts=[list(range(5)),list(reversed(range(5)))]
    for budget in (0,1,7,60):
        result=e.insertion(starts,features,budget)
        assert result['evaluations']<=budget
        assert result['score']<=min(features.score(features.amplitude(s)) for s in starts)+1e-12
        assert result['score']==pytest.approx(features.score(features.amplitude(result['order'])))
    order=starts[0];a=features.amplitude(order);features.crossing(a,order[0],order[1])
    swapped=[order[1],order[0]]+order[2:]
    assert np.linalg.norm(a-features.amplitude(swapped))<1e-12


@pytest.mark.parametrize('lengths',[[3],[2,4],[1,2,3]])
@pytest.mark.parametrize('steps',[1,3])
def test_symmetric_accounting_and_total_duration(lengths,steps):
    ids,fractions=e.pulses(lengths,steps,True);U=sum(lengths)
    count=U if len(lengths)==1 else steps*(2*U-lengths[0]-lengths[-1])+lengths[0]
    assert len(ids)==count
    assert np.max(abs(np.bincount(ids,weights=fractions,minlength=U)-1))<1e-14
    assert e.feasible_steps(lengths,count,True)==(1 if len(lengths)==1 else steps)
    assert e.feasible_steps(lengths,min(lengths)-1,True)==0


def test_limits_fail_before_large_allocation_or_search():
    with pytest.raises(ValueError):target_matrix({'n':9})
    features=e.Features(units(),e.Reference({((0,'Z'),):1},3,1.,1))
    with pytest.raises(ValueError):e.insertion([list(range(5))],features,2_000_001)
    with pytest.raises(ValueError):e.feasible_steps([1],-1)
    with pytest.raises(ValueError):e.pulses([1],0)
