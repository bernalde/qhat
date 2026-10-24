"""Independent small-matrix checks before prospective molecular evaluation."""
import itertools
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.linalg import eigh,expm
from scipy.integrate import quad_vec
import engine as e
import campaign
from campaign import current_placement,target_matrix
from verify import (literal,tensor_hamiltonian,tensor_bridge,tensor_product_oracle,
                    check_current_ledger,check_reference_edge,independent_decision)


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
    check_current_ledger(result,result['evaluations'],len(diag),len(off))
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


def test_tensor_bridge_detects_off_diagonal_loss_with_unchanged_HF_energy():
    one=np.array([[1.,.2],[.2,2.]]);two=np.zeros((2,2,2,2));constant=.4
    raw=tensor_hamiltonian(one,two,constant)
    policy=dict(input_constant=constant,removed_parent_pauli_scalar=1.5,
                canonical_hermitian_adjustment_fermionic_l1=0.,discarded_parent_fermionic_l1=0.,parent_cutoff=1e-12)
    target=raw-(constant+1.5)*np.eye(4)
    tensor_bridge(one,two,constant,target,policy,1,1.4)
    changed=one.copy();changed[0,1]+=.02;changed[1,0]+=.02
    assert tensor_hamiltonian(changed,two,constant)[1,1]==raw[1,1]
    with pytest.raises(ValueError,match='Raw tensors'):
        tensor_bridge(changed,two,constant,target,policy,1,1.4)


@pytest.mark.parametrize('n',[2,4])
def test_conventional_tensor_product_to_occupation_indexing(n):
    words=([tuple((q,p) for q,p in enumerate(labels) if p!='I')
            for labels in itertools.product('IXYZ',repeat=2)] if n==2 else
           [((0,'Y'),(2,'Y'),(3,'Z')),((0,'Z'),(1,'X'),(3,'Y'))])
    for word in words:
        assert np.array_equal(literal(word,n)[0].toarray(),tensor_product_oracle(word,n))


@pytest.mark.parametrize('phase,marker',[
    ('prepare','ALL_SCHEDULES_FROZEN.json'),('prepare','designs'),('prepare','rows'),
    ('run','references'),('run','states'),('run','rows'),('run','run_summary.json')])
def test_frozen_phase_outputs_cannot_be_overwritten(tmp_path,phase,marker):
    target=tmp_path/marker
    if '.' in marker:target.write_text('preserve original evidence')
    else:
        target.mkdir();(target/'evidence.txt').write_text('preserve original evidence')
    before={str(p.relative_to(tmp_path)):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    with pytest.raises(FileExistsError):getattr(campaign,phase)(SimpleNamespace(output=tmp_path))
    after={str(p.relative_to(tmp_path)):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    assert before==after


def test_reference_guards_and_missing_partner_are_rejected():
    with pytest.raises(ValueError):e.Reference({((0,'X'),):1},3,1.,1)
    with pytest.raises(ValueError):e.Reference({((0,'Z'),):1},3,1.,1,3,1)
    with pytest.raises(ValueError):e.choose_partner([units()[3],units()[3]],3)
    check_reference_edge(dict(partner=None,coupling_numerator='0'),[])
    with pytest.raises(ValueError):
        check_reference_edge(dict(partner=None,coupling_numerator='0'),[(-2,5,6,2)])


def test_decision_predicates_and_independent_summary_with_floor_losses():
    from decisions import compute,ARMS
    current=[.1,.1,1e-13,1e-13,.1,.1]
    corrected=[.05,.1+5e-13,2e-13,6e-12,.3,.02]
    contexts={};rows=[]
    for point,cost in enumerate((1.,1.5,2.,2.5,1.,1.)):
        case='case'+str(point)
        contexts[case,'HF']={arm:dict(setup_cpu_seconds=cost if arm=='pair_reference' else 1.,
                                    search_cpu_seconds=0.) for arm in ARMS}
        for budget in (5,10,20,50):
            for arm in ARMS:
                rows.append(dict(case=case,state='HF',arm=arm,budget_units=budget,
                    one_minus_overlap=corrected[point] if arm=='pair_reference' else current[point]))
    result=compute(rows,contexts)
    assert result==independent_decision(rows,contexts)
    assert (result['primary']['wins'],result['primary']['ties'],result['primary']['losses'])==(2,2,2)
    assert result['primary']['floor_exclusions']==2
    assert result['maximum_strict_loss_ratio']==6
    assert result['median_construction_CPU_ratio']==1.25
    assert result['conditions']==dict(at_least_four_wins=False,no_loss_exceeds_twofold=False,
                                     median_construction_CPU_at_most_twofold=True)
    assert not result['correction_passes']
