"""Analytic small-register oracles and nonmonotone target-selection regressions."""
from fractions import Fraction
import cmath,json
from pathlib import Path
import numpy as np
import pytest
from compiler import Compiler,Gate,VERSION,write_trace
from reference import Reference,metrics
import study
from verify import independent_choices


def model(n=2,electrons=0,word=((0,'Z'),),coefficient=1,denominator=1):
    term=[[[q,p] for q,p in word],str(coefficient)]
    return dict(n=n,electrons=electrons,denominator=str(denominator),coefficients=[term],blocks=[[term]])


def test_literal_primitive_replay_and_independent_parallel_depth(tmp_path):
    compiler=Compiler(2)
    for gate in [Gate('h',0),Gate('cx',0,1),Gate('rz',1,angle=Fraction(2,3)),Gate('h',0)]:compiler.add(gate)
    path=tmp_path/'bell.ndjson.gz';identity={'fixture':'Bell phase and disjoint final gates'}
    write_trace(path,2,compiler.gates,identity)
    state,costs=Reference(model()).audit_compiled_trace(path,identity,VERSION)
    a,b=cmath.exp(-1j/3),cmath.exp(1j/3)
    np.testing.assert_allclose(state,np.array([a,a,b,-b])/2,atol=5e-16,rtol=0)
    assert costs==dict(gates=4,single_qubit_gates=3,two_qubit_gates=1,unit_duration_ASAP_depth=3,CX_dependency_depth=1)


def test_scaled_taylor_target_and_unfused_formula_against_closed_form():
    reference=Reference(model(n=1,word=((0,'X'),),coefficient=3,denominator=7))
    T=2.11;expected=np.array([np.cos(3*T/7),-1j*np.sin(3*T/7)])
    np.testing.assert_allclose(reference.exact_action(T),expected,atol=2e-14,rtol=0)
    for formula in (1,2):
        record=dict(formula=formula,schedule=[[[[0,'X']],'3']],unit_lengths=[1])
        for steps in (1,4,128):np.testing.assert_allclose(reference.ideal_formula(record,steps,T),expected,atol=2e-13,rtol=0)


def test_joint_spin_sector_is_distinct_from_particle_number_and_global_phase():
    target=model(n=2,electrons=1);exact=np.array([0,1,0,0],complex)
    phase=metrics(exact,1j*exact,target)
    assert phase['error']==phase['infidelity']==phase['number_leakage']==phase['joint_sector_leakage']==0
    wrong_spin=metrics(exact,np.array([0,0,1,0],complex),target)
    assert wrong_spin['number_leakage']==0 and wrong_spin['joint_sector_leakage']==1 and wrong_spin['error']==1


def selection_fixture(tmp_path,all_unmet=False):
    for name in ('curves','candidates','rows'):(tmp_path/name).mkdir()
    curve=dict(curve='probe',case='probe',T=.83,arm='cubic_insertion',formula=2)
    (tmp_path/'curves/probe.json').write_text(json.dumps(curve))
    errors={1:.02,2:.009,4:.004,8:.001,16:.0005,64:.01,128:.002}
    resources={1:(100,10,150),2:(20,10,40),4:(10,20,30),8:(15,5,25),16:(50,40,100),64:(1,1,2),128:(10,10,35)}
    candidates={};rows={};independent={}
    for steps in study.GRID:
        label='probe.r'+format(steps,'03d');candidate=dict(candidate=label,status='censored' if steps==32 else 'complete')
        candidates[label]=candidate;(tmp_path/'candidates'/(label+'.json')).write_text(json.dumps(candidate))
        if steps==32:continue
        cx,depth,gates=resources[steps];row=dict(candidate=label,steps=steps,error=.02 if all_unmet else errors[steps],
            compiled_resources=dict(two_qubit_gates=cx,unit_duration_ASAP_depth=depth,gates=gates))
        rows[label]=row;independent[label]={'error':row['error']};(tmp_path/'rows'/(label+'.json')).write_text(json.dumps(row))
    return {'probe':curve},candidates,rows,independent


def test_target_selection_retains_ambiguity_censoring_and_nonmonotone_cost(tmp_path):
    curves,candidates,rows,errors=selection_fixture(tmp_path)
    actual=independent_choices(curves,candidates,rows,errors)
    assert actual==study.selections(tmp_path)
    loose,tight=actual
    assert loose['selected_candidate']=='probe.r128' and loose['qualifying_candidates']==5
    assert loose['ambiguous_candidates']==['probe.r064'] and loose['censored_candidates']==1
    assert tight['selected_candidate']=='probe.r016' and tight['qualifying_candidates']==1
    assert tight['ambiguous_candidates']==['probe.r008'] and tight['censored_candidates']==1


def test_unmet_target_keeps_all_candidates(tmp_path):
    curves,candidates,rows,errors=selection_fixture(tmp_path,all_unmet=True)
    actual=independent_choices(curves,candidates,rows,errors)
    assert actual==study.selections(tmp_path)
    assert len(actual)==2 and all(x['status']=='unmet' and x['selected_candidate'] is None and x['feasible_candidates']==7 and x['censored_candidates']==1 for x in actual)


def test_independent_trace_rejects_wrong_identity_and_noncanonical_angle(tmp_path):
    import gzip
    path=tmp_path/'invalid.ndjson.gz';header=dict(compiler=VERSION,n=2,identity={'expected':True})
    with gzip.open(path,'wt') as stream:
        stream.write(json.dumps(header)+'\n'+json.dumps(['rz',0,'2','6'])+'\n')
    reference=Reference(model())
    with pytest.raises(ValueError,match='header'):reference.audit_compiled_trace(path,{'expected':False},VERSION)
    with pytest.raises(ValueError,match='Noncanonical'):reference.audit_compiled_trace(path,header['identity'],VERSION)


def test_noncommuting_unit_is_refused_before_boundary_fusion():
    bad=dict(n=1,electrons=0,denominator='1',coefficients=[[[[0,'X']],'1'],[[[0,'Z']],'1']],
        blocks=[[[[[0,'X']],'1'],[[[0,'Z']],'1']]])
    with pytest.raises(ValueError,match='commuting'):study.validate_model(bad)
    with pytest.raises(ValueError,match='commutation'):Reference(bad)


@pytest.mark.parametrize('implementation',['producer','literal'])
def test_deep_hadamard_trace_preserves_raw_norm_without_renormalizing(tmp_path,implementation):
    # H^2=I, so this deep circuit has a hand-known state and norm. Binary64
    # repeated basis changes lose over2e-12 here; no intermediate normalization
    # is permitted to conceal the drift measured by the independent metric.
    path=tmp_path/'deep.ndjson.gz';identity={'fixture':'50000literalHadamards'}
    write_trace(path,1,[Gate('h',0)]*50_000,identity)
    if implementation=='producer':state=study.PrimitiveReplay(model(n=1)).replay(path,identity)
    else:state,_=Reference(model(n=1)).audit_compiled_trace(path,identity,VERSION)
    assert abs(float(np.linalg.norm(state))-1)<=2e-12
    np.testing.assert_allclose(state,np.array([1,0]),atol=2e-12,rtol=0)
