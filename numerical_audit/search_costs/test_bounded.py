"""Budget failure paths, factorial initialization and exact physical model checks."""
from fractions import Fraction
import itertools
import numpy as np
import pytest
from scipy import sparse
from bounded import Features,Exhausted,insertion,run_arm,engine as e
from models import synthetic


def car(key,n):
    rows=[];columns=[];values=[]
    for v in range(1<<n):
        out=v;amplitude=1
        for p,creation in reversed(key):
            occupied=(out>>p)&1
            if occupied==creation:break
            amplitude*=(-1)**((out&((1<<p)-1)).bit_count());out^=1<<p
        else:rows.append(out);columns.append(v);values.append(amplitude)
    return sparse.csr_matrix((values,(rows,columns)),shape=(1<<n,1<<n),dtype=complex)


def pauli(word,n):
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    labels=dict(word);out=sparse.csr_matrix([[1.+0j]])
    for q in reversed(range(n)):out=sparse.kron(out,matrices[labels.get(q,'I')],format='csr')
    return out


def decode(raw):return {tuple((q,p) for q,p in w):int(c) for w,c in raw}


def test_exact_scaling_model_against_car():
    n=4;model=synthetic(n);den=int(model['denominator']);dimension=1<<n
    total=decode(model['coefficients']);target=sum((c/den)*pauli(w,n) for w,c in total.items())
    direct=sparse.csr_matrix((dimension,dimension),dtype=complex)
    for p in range(n):
        epsilon=Fraction(2*(p//2+1)+(p%2),20)
        direct+=float(epsilon)*car(((p,1),(p,0)),n)
        for q in range(p+1,n):
            if p%2==q%2:
                t=-Fraction(3,10*(1+(q-p)//2))
                direct+=float(t)*(car(((p,1),(q,0)),n)+car(((q,1),(p,0)),n))
            U=Fraction(3,5) if p//2==q//2 else Fraction(1,10*(1+q-p))
            direct+=float(U)*(car(((p,1),(p,0)),n)@car(((q,1),(q,0)),n))
    assert np.max(abs((direct-target-float(Fraction(model['removed_scalar_exact']))*sparse.eye(dimension)).toarray()))<1e-14
    blocks=[decode(b) for b in model['blocks']]
    assert len({w for b in blocks for w in b})==len(total)
    for b in blocks:
        G=sum((c/den)*pauli(w,n) for w,c in b.items())
        for parity in (0,1):
            N=sparse.diags([sum((v>>q)&1 for q in range(parity,n,2)) for v in range(dimension)],dtype=float)
            assert np.max(abs((G@N-N@G).toarray()))<1e-14


def fixture():
    blocks=[{((0,'Z'),):2},{((1,'Z'),):1},{((0,'X'),):3}]
    return blocks,{**blocks[0],**blocks[1]}


def test_bounded_features_match_prior_definition():
    blocks,D=fixture();reference=e.Reference(D,0,1.,7)
    features=Features(blocks,reference,float('inf'))
    previous=e.Features(blocks,reference)
    for order in itertools.permutations(range(3)):
        assert np.isclose(features.score(features.amplitude(order)),previous.score(previous.amplitude(order)),atol=1e-15)


def test_refuse_setup_before_any_allocation():
    blocks,D=fixture()
    with pytest.raises(Exhausted) as result:Features(blocks,e.Reference(D,0,1.,7),1,clock=lambda:2)
    assert result.value.phase=='before feature setup' and result.value.pairs==0
    with pytest.raises(ValueError):Features([blocks[0]]*257,e.Reference(D,0,1.,7),float('inf'))


def test_partial_feature_setup_is_retained():
    blocks,D=fixture();ticks=iter([0.,0.,0.,2.])
    with pytest.raises(Exhausted) as result:Features(blocks,e.Reference(D,0,1.,7),1.,clock=lambda:next(ticks))
    assert result.value.pairs==2


def test_incomplete_sweep_counts_actual_positions():
    blocks,D=fixture();features=Features(blocks,e.Reference(D,0,1.,7),float('inf'))
    result=insertion([[0,1,2],[2,1,0]],features,float('inf'),evaluation_cap=2)
    assert result['evaluations']==2
    assert all(r['evaluations']==1 and len(r['sweeps'])==1 and not r['sweeps'][0]['complete'] for r in result['starts'])
    assert sorted(result['order'])==[0,1,2]


def test_initializer_charge_cannot_disappear():
    blocks,D=fixture();record=run_arm(blocks,D,0,7,'JD','placement',[2,1,0],.5,.1,.6,100)
    assert record['status'].startswith('CPU overrun retained: initializer charge')
    assert record['inherited_placement_CPU_charge']==.6 and record['common_setup_CPU_charge']==.1
    assert record['score'] is None and record['evaluations']==0 and record['order']==[2,1,0]
    assert not record['within_CPU_budget']


def test_objective_and_initializer_are_separate_factors():
    blocks,D=fixture()
    for objective in ('JD','J0'):
        result=run_arm(blocks,D,0,7,objective,'placement',[2,1,0],2.,0.,.01,100)
        assert result['objective']==objective and result['initialization']=='placement'
        assert result['inherited_placement_CPU_charge']==.01
        assert result['within_CPU_budget'] and result['starts'][0]['order']


@pytest.mark.parametrize('n',(0,2,14,100))
def test_refuse_undeclared_exponential_dimensions(n):
    with pytest.raises(ValueError):synthetic(n)


def test_independent_objectives_and_unfused_circuits():
    from study import common,CircuitTarget
    from verify_costs import Reference
    model=synthetic(4);context=common(model);target=CircuitTarget(model);reference=Reference(model)
    for partition in ('conserving','pauli'):
        blocks=context['paulis'] if partition=='pauli' else context['blocks']
        orders=[list(range(len(blocks))),list(reversed(range(len(blocks))))]
        for order in orders:
            record=dict(partition=partition,order=order,objective='JD',
                unit_lengths=[len(blocks[i]) for i in order],
                schedule=[[[[q,p] for q,p in w],str(c)] for i in order for w,c in sorted(blocks[i].items())])
            scores=reference.scores(record)
            for objective in ('JD','J0'):
                diagonal=context['D'] if objective=='JD' else {():0}
                features=Features(blocks,e.Reference(diagonal,context['v'],1.,context['den']),float('inf'))
                assert np.isclose(features.score(features.amplitude(order)),scores[objective],rtol=1e-12,atol=1e-16)
            for formula in (1,2):
                steps=e.feasible_steps(record['unit_lengths'],10*len(model['coefficients']),formula==2)
                actual,count=target.propagate(record,steps,formula)
                independent=reference.replay(record,steps,formula)
                assert count<=10*len(model['coefficients'])
                assert np.max(abs(actual-independent))<2e-13
                assert abs(np.linalg.norm(actual)-1)<2e-13
