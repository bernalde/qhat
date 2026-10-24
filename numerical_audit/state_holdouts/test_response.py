"""Independent mathematical checks for the generalized coherent objective."""
import itertools
import numpy as np
import pytest
from scipy.integrate import quad_vec
from response import Features,canonical_eigensystem,determinant as e
from verify import literal,independent_score,replay,independent_occupations
from study import propagate,encode_schedule,natural_occupations


def fixture():
    blocks=[{((0,'Z'),):3},{((1,'Z'),):-2},{((0,'X'),):2},
            {((0,'X'),(1,'Z')):1},{((1,'X'),):1}]
    diagonal={**blocks[0],**blocks[1]}
    psi=np.array([1,1j,-.2,.3j],complex);psi/=np.linalg.norm(psi)
    return blocks,diagonal,psi


def test_complex_response_quadrature_and_updates():
    blocks,D,psi=fixture();T=.83;den=7
    f=Features(blocks,D,psi,T,den,2)
    matrices=[sum((c/den)*literal(w,2) for w,c in b.items()).toarray() for b in blocks]
    d=sum((c/den)*literal(w,2) for w,c in D.items()).diagonal()
    for order in itertools.permutations(range(len(blocks))):
        C=sum((-.5j)*(matrices[i]@matrices[j]-matrices[j]@matrices[i])
              for a,i in enumerate(order) for j in order[a+1:])
        def action(t):
            return np.exp(1j*t*d)*(C@(np.exp(-1j*t*d)*psi))
        integrated,_=quad_vec(action,0,T,epsabs=1e-12)
        integrated-=psi*np.vdot(psi,integrated)
        assert np.allclose(f.amplitude(order),integrated,atol=1e-12)
        for at in range(len(order)-1):
            swap=list(order);swap[at],swap[at+1]=swap[at+1],swap[at]
            updated=f.amplitude(order);f.crossing(updated,order[at],order[at+1])
            assert np.allclose(updated,f.amplitude(swap),atol=1e-12)


def test_interference_is_not_mixture():
    blocks,D,psi=fixture();order=list(range(len(blocks)))
    coherent=Features(blocks,D,psi,.83,7,2).score(Features(blocks,D,psi,.83,7,2).amplitude(order))
    mixture=0.
    for i in range(4):
        ket=np.eye(4,dtype=complex)[:,i];f=Features(blocks,D,ket,.83,7,2)
        mixture+=abs(psi[i])**2*f.score(f.amplitude(order))
    assert abs(coherent-mixture)>1e-5


def test_determinant_reduction():
    blocks,D,_=fixture();psi=np.eye(4,dtype=complex)[:,1]
    general=Features(blocks,D,psi,.83,7,2)
    old=e.Features(blocks,e.Reference(D,1,.83,7))
    for order in itertools.permutations(range(len(blocks))):
        assert np.isclose(general.score(general.amplitude(order)),old.score(old.amplitude(order)),atol=1e-15)


def test_degenerate_phase_and_basis():
    H=np.diag([1.,0.,0.,2.]);values,U=canonical_eigensystem(H)
    assert np.array_equal(values,[0,0,1,2])
    assert np.allclose(U,np.eye(4)[:,[1,2,0,3]])


def test_full_space_and_fusion():
    blocks,D,psi=fixture();total={w:c for b in blocks for w,c in b.items()}
    model=dict(n=2,denominator='7',coefficients=[[w,str(c)] for w,c in total.items()],
               blocks=[[[w,str(c)] for w,c in b.items()] for b in blocks])
    for formula in (1,2):
        method=dict(order=list(range(5)),partition='conserving',formula=formula,
                    unit_lengths=[1]*5,schedule=[[w,str(c)] for w,c in total.items()])
        actual,count=propagate(psi,method,model,.83,3)
        assert np.allclose(actual,replay(model,method,psi,.83,3),atol=1e-13)
        assert count==(15 if formula==1 else 25)
        f=Features(blocks,D,psi,.83,7,2)
        assert np.isclose(independent_score(model,method,psi,.83),f.score(f.amplitude(method['order'])),atol=1e-14)


@pytest.mark.parametrize('n,m,T,den',[(9,2,1,1),(2,129,1,1),(2,2,0,1),(2,2,1,0)])
def test_guard_before_allocation(n,m,T,den):
    with pytest.raises(ValueError):Features([{((0,'Z'),):1}]*m,{((0,'Z'),):1},[],T,den,n)


def test_normalization_and_order_guards():
    blocks,D,psi=fixture()
    with pytest.raises(ValueError):Features(blocks,D,2*psi,1,7,2)
    f=Features(blocks,D,psi,1,7,2)
    with pytest.raises(ValueError):f.amplitude([0,1,2,3,3])
    with pytest.raises(ValueError):e.insertion([[0,1,2,3,4]],f,2_000_001)
    with pytest.raises(ValueError):Features(blocks,D,psi,float('nan'),7,2)
    with pytest.raises(ValueError):Features(blocks,D,psi*np.nan,1,7,2)
    with pytest.raises(ValueError):Features([{((2,'X'),):1}],D,psi,1,7,2)


def test_schedule_serialization_matches_frozen_json():
    assert encode_schedule([(((0,'X'),(1,'Z')),3)])==[[[[0,'X'],[1,'Z']],'3']]


def test_correlation_distinguishes_rotated_determinant():
    for indices,expected in (([3],0.),([1,2],0.),([3,12],1.)):
        psi=np.zeros(16,complex);psi[indices]=1/np.sqrt(len(indices))
        occupations,defect=natural_occupations(psi,4)
        other,other_defect=independent_occupations(psi,4)
        assert abs(defect-expected)<1e-12
        assert np.allclose(occupations,other,atol=1e-12) and abs(defect-other_defect)<1e-12
