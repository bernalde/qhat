"""Hand-derived Strang coefficients and an independent commutator oracle."""
import itertools
from pathlib import Path
import cubic as cubic_module
import independent as independent_module
assert Path(cubic_module.__file__).resolve().parent==Path(__file__).resolve().parent
assert Path(independent_module.__file__).resolve().parent==Path(__file__).resolve().parent
import numpy as np
import pytest
from cubic import Context,FirstOrderContext,Interrupted,insertion,engine
from independent import Reference


def algebra_model(three=False):
    # Deliberately nonconserving small algebra fixtures, not study targets.
    terms=[([[(0,'Z')]],2),([[(0,'X')]],3)]
    if three:terms.append(([[(0,'Z'),(1,'Z')]],5))
    coefficients=[[[[q,p] for q,p in words[0]],str(c)] for words,c in terms]
    return dict(n=2,electrons=1,denominator='7',coefficients=coefficients,blocks=[[r] for r in coefficients])


def test_hand_derived_two_unit_cubic_and_new_orientation():
    model=algebra_model();context=Context(model);independent=Reference(model)
    # A=(2/7)Z, B=(3/7)X, input|01>. For A/2 B A/2,
    # A3=-[A,[A,B]]/24+[B,[B,A]]/12=-(2/343)X+(6/343)Z.
    expected=np.array([-2/343,-6/343,0.,0.])
    reverse=np.array([4/343,3/343,0.,0.])
    assert np.allclose(context.cubic([0,1]),expected,rtol=1e-12,atol=1e-16)
    assert np.allclose(context.cubic([1,0]),reverse,rtol=1e-12,atol=1e-16)
    assert np.allclose(independent.cubic([0,1]),expected,rtol=1e-12,atol=1e-16)
    assert np.allclose(independent.cubic([1,0]),reverse,rtol=1e-12,atol=1e-16)
    better=context.score([0,1],1.,20);worse=context.score([1,0],1.,20)
    assert np.isclose(worse['score'],4*better['score'],rtol=1e-11)
    # First-order response reversal has equal norm; this is a distinct objective.
    blocks=context.blocks;features=engine.Features(blocks,engine.Reference(blocks[0],context.v,1.,7))
    assert np.isclose(features.score(features.amplitude([0,1])),features.score(features.amplitude([1,0])))


def test_three_unit_polynomial_against_nested_commutators():
    model=algebra_model(True);context=Context(model);independent=Reference(model)
    for order in itertools.permutations(range(3)):
        actual=context.score(order,.83,90);other=independent.score(order,.83,90)
        assert np.allclose(actual['cubic_coefficients'],other['cubic_coefficients'],rtol=1e-11,atol=1e-16)
        assert actual['steps']==other['steps']
        assert np.isclose(actual['score'],other['score'],rtol=1e-11,atol=1e-22)


def test_first_order_control_is_recomputed_at_the_same_time():
    model=algebra_model(True);reference=Reference(model)
    for T in (.83,2.11):
        for partition in ('conserving','pauli'):
            context=FirstOrderContext(model,T,partition)
            for order in itertools.permutations(range(context.m)):
                actual=context.score(order,T,90);other=reference.first_score(order,T,90,partition)
                assert actual['steps']==other['steps']
                assert np.isclose(actual['score'],other['score'],rtol=1e-11,atol=1e-20)


def test_complete_hopping_unit_and_boundary_affordability():
    diagonal=[[[[0,'Z']],'2']]
    hopping=[[[[0,'X'],[1,'Z'],[2,'X']],'1'],[[[0,'Y'],[1,'Z'],[2,'Y']],'1']]
    model=dict(n=4,electrons=2,denominator='7',coefficients=diagonal+hopping,blocks=[diagonal,hopping])
    context=Context(model);independent=Reference(model);independent.conservation()
    records=[context.score(order,2.11,31) for order in ([0,1],[1,0])]
    assert [r['steps'] for r in records]==[10,9]
    for order,actual in zip(([0,1],[1,0]),records):
        other=independent.score(order,2.11,31)
        assert actual['steps']==other['steps']
        assert np.allclose(actual['cubic_coefficients'],other['cubic_coefficients'],rtol=1e-11,atol=1e-16)


def test_improving_partial_sweep_is_accepted_and_rechecked():
    model=algebra_model();context=Context(model);independent=Reference(model)
    seed=context.score([1,0],1.,20)
    result=insertion(context,[seed,seed,seed],1.,20,evaluation_cap=3)
    assert result['attempted']==result['completed']==3
    assert result['attempted_rescores']==result['accepted_rescores']==3
    for start in result['starts']:
        assert start['attempted']==start['completed']==1 and len(start['accepted'])==1
        assert len(start['sweeps'])==1 and not start['sweeps'][0]['complete']
        assert start['record']['score']<start['initial_score'] and start['record']['order']==[0,1]
        assert np.isclose(start['record']['score'],independent.score([0,1],1.,20)['score'],rtol=1e-11)


def test_interrupted_candidate_has_actual_stage_count():
    context=Context(algebra_model())
    with pytest.raises(Interrupted) as error:context.cubic([0,1],deadline=0.,clock=lambda:1.)
    assert error.value.completed_stages==0


@pytest.mark.parametrize('constructor',[Context,Reference])
def test_register_limit_before_other_model_access(constructor):
    with pytest.raises(ValueError):constructor({'n':1000000})


@pytest.mark.parametrize('order',[[0,0],[False,1],[0.,1.]])
def test_nonintegral_or_repeated_order_refused(order):
    context=Context(algebra_model())
    with pytest.raises(ValueError):context.cubic(order)


def test_no_quadratic_feature_pool_or_oversized_allocation(monkeypatch):
    import cubic
    calls=[]
    def forbidden(*args,**kwargs):calls.append(1);raise AssertionError('array allocation reached')
    monkeypatch.setattr(cubic.np,'arange',forbidden)
    model=dict(n=2,electrons=1,denominator='1',blocks=[None]*129,coefficients=[])
    with pytest.raises(ValueError):Context(model)
    assert not calls


def test_interrupted_accepted_rescore_keeps_incumbent_and_stage_detail(monkeypatch):
    context=Context(algebra_model());seed=context.score([1,0],1.,20)
    class Clock:
        now=0.;inside=False;stages=0
        def __call__(self):
            if self.inside:
                self.stages+=1
                if self.stages>2:self.now=10.
            return self.now
    clock=Clock();original=context.score;calls=0
    def score(order,T,cap,deadline=float('inf'),clock=clock):
        nonlocal calls
        calls+=1;clock.inside=calls==3
        try:return original(order,T,cap,deadline,clock)
        finally:clock.inside=False
    monkeypatch.setattr(context,'score',score)
    result=insertion(context,[seed,seed,seed],1.,20,clock=clock)
    first=result['starts'][0]
    assert first['termination']=='CPU during accepted rescore'
    assert first['attempted']==first['completed']==2
    assert result['attempted_rescores']==1 and result['accepted_rescores']==0
    assert first['record']==result['record']==seed and not first['accepted']
    assert first['interrupted_rescore_completed_stages']==2
    assert all(start['attempted']==0 for start in result['starts'][1:])


def test_single_unit_cap_is_enforced_by_both_independent_and_producer_scores():
    hopping=[[[[0,'X'],[1,'Z'],[2,'X']],'1'],[[[0,'Y'],[1,'Z'],[2,'Y']],'1']]
    model=dict(n=4,electrons=2,denominator='7',coefficients=hopping,blocks=[hopping])
    for constructor in (Context,Reference):
        context=constructor(model)
        with pytest.raises(ValueError):context.score([0],.83,1)
        assert context.score([0],.83,2)['steps']==1
