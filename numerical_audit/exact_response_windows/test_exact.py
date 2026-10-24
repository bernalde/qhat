"""Regressions for exact coefficient representation and restricted domains."""
import copy
import itertools
import math
import pytest
from exact import amplitudes,objective,solve,freeze_domain,validate_model,pair_data
from verify import integer_action,full_amplitudes,independent_minimum,independent_pairs,Reference,diagonal_response_weight


def fixture():
    # Three responses to one destination can cancel. Fill the domain with
    # six inert labels so the complete 9! enumeration has degenerate optima.
    pairs=[[0,1,2,'7'],[0,2,2,'11'],[1,2,2,'13']]
    order=list(range(9))
    return dict(width=9,window=order,warm_order=order,outside_prefix=[],outside_suffix=[],
        pairs=pairs,integer_weights=[[2,'3']],objective_denominator='5',offset=0)


def test_exact_global_minimum_and_independent_coverage():
    d=fixture();c=solve(d,120)
    expected=min((objective(list(p)+list(range(3,9)),d),p) for p in itertools.permutations(range(3)))[0]
    assert c['complete'] and c['evaluated']==math.factorial(9)
    assert int(c['incumbent_integer'])==expected<int(c['warm_start_integer'])
    assert c['lower_bound_integer']==c['incumbent_integer'] and c['absolute_gap_integer']=='0'
    other=independent_minimum(d['warm_order'],0,9,d['pairs'],{2:3},120)
    for k in ('incumbent_integer','selected_window','minimizers','evaluated'):assert c[k]==other[k]
    assert sum(3*b*b for b in full_amplitudes(c['order'],d['pairs']).values())==expected


def test_outside_order_and_cross_window_pairs_remain_fixed():
    d=fixture();d.update(window=list(range(1,10)),warm_order=list(range(11)),outside_prefix=[0],outside_suffix=[10],offset=1)
    d['pairs']=[[0,1,2,'7'],[0,10,2,'-3'],[1,2,2,'11'],[2,10,2,'13'],[1,10,2,'-15']]
    c=solve(d,120);assert c['order'][0]==0 and c['order'][-1]==10
    assert int(c['incumbent_integer'])<int(c['warm_start_integer'])
    other=independent_minimum(d['warm_order'],1,9,d['pairs'],{2:3},120)
    assert c['incumbent_integer']==other['incumbent_integer'] and c['selected_window']==other['selected_window']


def test_interruption_retains_incumbent_and_valid_lower_bound():
    c=solve(fixture(),0)
    assert not c['complete'] and c['evaluated']==0
    assert c['lower_bound_integer']=='0' and c['incumbent_integer']==c['warm_start_integer']
    assert c['status']=='incomplete' and 'CPU limit' in c['termination']


@pytest.mark.parametrize('order',[[0]*9,[False]+list(range(1,9)),list(range(8))])
def test_bad_domain_rejected(order):
    d=fixture();d['warm_order']=order
    with pytest.raises(ValueError):solve(d)


def test_arbitrary_precision_pauli_actions_and_supplied_dyadic_weights():
    coefficient=(1<<100)+1
    assert integer_action({((0,'Y'),(1,'Y')):coefficient},0)==(3,-coefficient)
    assert integer_action({((0,'Y'),(1,'Y')):coefficient},1)==(2,coefficient)
    d=fixture();d['integer_weights']=[[2,str((1<<70)+3)]]
    a=amplitudes(d['warm_order'],d['pairs']);assert objective(d['warm_order'],d)==((1<<70)+3)*a[2]**2


def test_hand_computed_commutator_and_normalization():
    # This algebra fixture deliberately uses X0, not a physical conserving panel.
    # Z0 and Z0Z1 each contribute -2 to the |0001> commutator column.
    supports=[[(0,'Z')],[(0,'Z'),(1,'Z')],[(0,'X')],[(1,'Z')],[(2,'Z')],[(3,'Z')],
        [(1,'Z'),(2,'Z')],[(1,'Z'),(3,'Z')],[(2,'Z'),(3,'Z')]]
    raw=[[[q,p] for q,p in word] for word in supports]
    model=dict(n=4,electrons=0,denominator='7',coefficients=[[w,'1'] for w in raw],blocks=[[[w,'1']] for w in raw])
    data=pair_data(model);pairs,_=independent_pairs(model)
    assert pairs==data['pairs']==[[0,2,1,'-2'],[1,2,1,'-2']]
    assert data['energy_gap_numerators']==[[1,'-4']]
    _,num,divisor=data['weight_ratios'][0];rho=int(num)/int(divisor)
    # sin(x)^2=(1-cos(2x))/2, x=gap/(2den)=-2/7.
    expected=(1-math.cos(4/7))/(2*(2/7)**2)
    assert abs(rho-expected)<5e-15 and abs(diagonal_response_weight(-4,7)-rho)<5e-15
    assert diagonal_response_weight(0,7)==1. and diagonal_response_weight(4,7)==diagonal_response_weight(-4,7)
    weight=int(data['integer_weights'][0][1]);assert objective(list(range(9)),data)==16*weight
    assert objective([0,2,1]+list(range(3,9)),data)==0
    assert int(data['objective_denominator'])==4*7**4*int(data['weight_scale'])


@pytest.mark.parametrize('call',[validate_model,independent_pairs,Reference])
def test_oversized_register_refused_before_model_access_or_allocation(call):
    # Other fields deliberately absent: the bound must stop before reading them.
    with pytest.raises(ValueError):call({'n':1000000})
