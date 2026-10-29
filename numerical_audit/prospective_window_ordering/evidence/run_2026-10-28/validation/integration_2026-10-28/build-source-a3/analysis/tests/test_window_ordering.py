"""Normal API/CLI tests: real native witnesses, cheap owned fake workers.

The fakes exercise process/transport refusals, never substitute for the separate
real-Julia coupled tests. They use a native-scored tiny fixture as a witness.
"""
from dataclasses import asdict
from fractions import Fraction
from itertools import permutations
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap
import time

import pytest

from qhat.analysis import sparse_response_ordering as native
from qhat.analysis import window_ordering as window
from qhat.analysis.train_ordering import main


def test_installed_asset_manifest_binds_the_product_target():
    tomllib=pytest.importorskip('tomllib')
    assets=window.asset_path().parent
    manifest=tomllib.loads((assets/'SOURCE.toml').read_text())
    assert manifest['version']==1
    assert set(manifest['files'])=={'ResponseCore.jl','WindowSearch.jl','worker.jl','tests.jl','Project.toml','README.md'}
    for name,digest in manifest['files'].items():
        assert hashlib.sha256((assets/name).read_bytes()).hexdigest()==digest


def tiny_blocks():
    # Original lengths 1,2,1,2,1,2; identity is still a charged rotation.
    return ({(): Fraction(7,100)},
            {((0,'X'),(1,'X')): Fraction(1,2), ((0,'Y'),(1,'Y')): Fraction(1,2)},
            {((0,'Z'),): Fraction(3,10)},
            {((1,'X'),(2,'X')): Fraction(1,5), ((1,'Y'),(2,'Y')): Fraction(1,5)},
            {((1,'Z'),(2,'Z')): Fraction(13,100)},
            {((0,'Z'),): Fraction(11,100), ((2,'Z'),): Fraction(-19,100)})


@pytest.fixture
def witness():
    blocks=tiny_blocks(); context=native.SparseResponse(blocks,3,2,time=Fraction(37,100))
    scored=[]
    for interior in permutations((1,2,3,4)):
        order=(5,)+interior+(0,)
        scored.append((context.score(order,80).value,order))
    best=min(scored); worst=max(scored)
    assert best[0]<worst[0]-1e-12
    return blocks,worst[1],best[1],context.score(best[1],80)


def fake_worker(tmp_path, mode='valid', candidate=None, score=None):
    # argv is the actual Julia argv contract. The executable records one PID,
    # flags and the raw authenticated request for assertions across the boundary.
    executable=tmp_path/('worker_'+mode)
    record=tmp_path/'record.json'
    program='''\
#!/usr/bin/env python3
import hashlib,json,os,pathlib,signal,sys,time,tomllib
request_path, checkpoint_path = map(pathlib.Path,sys.argv[-4:-2])
raw=request_path.read_bytes();request=tomllib.loads(raw.decode())
record=pathlib.Path(RECORD)
record.write_text(json.dumps(dict(pid=os.getpid(),argv=sys.argv,request=request)))
mode=MODE
if mode=='crash': sys.exit(9)
if mode=='missing_reply': sys.exit(0)
if mode=='malformed': checkpoint_path.write_text('this = [ broken');sys.exit(0)
if mode=='deep_reply': checkpoint_path.write_text('x='+ '['*500+'0'+']'*500);sys.exit(0)
if mode=='deep_literal_reply': checkpoint_path.write_text("prefix='"+']'*500+"'\\nx="+'['*500+'0'+']'*500);sys.exit(0)
if mode=='startup_burn':
    while True: pass
if mode=='burn':
    started=time.process_time()
    while time.process_time()-started<.06: pass
order=CANDIDATE if CANDIDATE is not None else request['order']
score=SCORE if SCORE is not None else request['parent_score']
if mode=='budget_changed':
    assert request['max_evaluations']>=4, 'child lacks a parent/baseline/trial/changed-checkpoint prefix'
    order=list(request['order']);order[1],order[2]=order[2],order[1]
    score=request['parent_score'] # Native rescore, not this claim, remains authoritative.
lengths=request['word_lengths'];steps=(request['budget']-lengths[order[0]])//(2*sum(lengths)-lengths[order[0]]-lengths[order[-1]])
changed=order!=request['order']
windows=[[1,len(order)-2]] if changed else []
if mode=='outside': order=list(reversed(request['order']));windows=[[1,len(order)-2]]
digest=hashlib.sha256(raw).hexdigest() if mode!='hash' else '0'*64
values=dict(version=1,request_sha256=digest,order=order,word_lengths=lengths,
            score=float(score),cubic_response=float(score)*steps**4,steps=steps,status='complete',
            checked_updates=len(windows),windows=windows)
text=''.join(k+'='+json.dumps(v)+'\\n' for k,v in values.items())
text+='history=['+('{positions='+json.dumps(windows[0])+',baseline_labels='+json.dumps(request['order'][1:-1])+',proposal_labels='+json.dumps(order[1:-1])+',accepted=true}' if changed else '')+']\\n'
evaluation=1;full_scores=1;flat_candidates=0
if mode in ('budget_changed','budget_saturated','budget_excess'):
    assert request['max_evaluations']>=3, 'serialized child allocation must admit its first trial'
    evaluation=request['max_evaluations']+(1 if mode=='budget_excess' else 0)
    full_scores=2 if changed else 1
    flat_candidates=evaluation-full_scores
text+='[counters]\\nevaluation='+str(evaluation)+'\\naction=0\\natom=0\\nfull_score='+str(full_scores)+'\\nflat_candidate='+str(flat_candidates)+'\\nexact_witness=0\\ninsertion_move=0\\n[phases]\\nsetup_CPU_seconds=0.0\\n'
checkpoint_path.write_text(text)
if mode=='checkpoint_crash': sys.exit(7)
if mode=='checkpoint_burn':
    while True: pass
'''
    program=program.replace('RECORD',repr(str(record))).replace('MODE',repr(mode))
    program=program.replace('#!/usr/bin/env python3','#!'+sys.executable)
    program=program.replace('CANDIDATE',repr(list(candidate)) if candidate else 'None')
    program=program.replace('SCORE',repr(score.value) if score else 'None')
    executable.write_text(textwrap.dedent(program)); executable.chmod(0o700)
    return executable,record


def invoke(blocks,parent,executable,**kwargs):
    options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=str(executable))
    return native.train_ordering(blocks,3,2,time=Fraction(37,100),rotation_budget=80,
        cpu_seconds=2,max_evaluations=100,max_moves=0,starts=[parent],method='window',
        window_options=options,**kwargs)


def test_opt_in_roundtrip_keeps_labels_lengths_cost_and_native_authority(tmp_path,witness):
    blocks,parent,candidate,score=witness
    executable,record=fake_worker(tmp_path,candidate=candidate,score=score)
    result=invoke(blocks,parent,executable)
    assert result.order==candidate
    assert result.score==score
    assert result.stats['window_status']=='verified_improvement'
    assert result.stats['window_history'][0]['baseline_labels']==list(parent[1:-1])
    assert result.stats['window_history'][0]['proposal_labels']==list(candidate[1:-1])
    assert result.order[0]==parent[0] and result.order[-1]==parent[-1]
    request=json.loads(record.read_text())
    assert request['request']['order']==list(parent)
    assert request['request']['reference']==2 # arbitrary basis state, not HF reconstruction
    assert request['request']['word_lengths']==[1,2,1,2,1,2]
    assert Fraction(request['request']['time'])==Fraction(float(Fraction(37,100)))
    assert request['request']['budget']==80
    assert '--startup-file=no' in request['argv'] and '--threads=1' in request['argv']
    assert '--heap-size-hint=250M' in request['argv']
    lengths=request['request']['word_lengths'];r=result.score.steps
    assert r*(2*sum(lengths)-lengths[parent[0]]-lengths[parent[-1]])+lengths[parent[0]]<=80
    assert result.evaluations<=100
    assert result.CPU_seconds==pytest.approx(result.stats['parent_process_CPU_seconds']+result.stats['child_process_CPU_seconds'])
    with pytest.raises(ProcessLookupError): os.kill(request['pid'],0)
    assert set(asdict(result))=={'order','score','CPU_seconds','evaluations','accepted_moves','status','stats'}


@pytest.mark.parametrize('mode,expected',[('hash','invalid_worker_reply'),('malformed','invalid_worker_reply'),
    ('deep_reply','invalid_worker_reply'),('outside','invalid_worker_reply'),
    ('deep_literal_reply','invalid_worker_reply'),
    ('crash','child_failure'),('missing_reply','missing_worker_reply')])
def test_bad_or_crashed_worker_retains_checked_parent(tmp_path,witness,mode,expected):
    blocks,parent,_,_=witness
    executable,_=fake_worker(tmp_path,mode)
    result=invoke(blocks,parent,executable)
    assert result.order==parent and result.score is not None
    assert result.stats['window_status']==expected
    assert result.stats['child_process_CPU_seconds']>0


def test_native_verification_rejects_worker_claim(tmp_path,witness):
    blocks,parent,candidate,score=witness
    # Worker claims a different score for the correct order: full native score is
    # still run and the Julia/native disagreement must retain the checked parent.
    executable,_=fake_worker(tmp_path,candidate=candidate,score=native.Score(score.value/2,score.cubic_response/2,score.steps))
    result=invoke(blocks,parent,executable)
    assert result.order==parent
    assert result.stats['window_status']=='native_julia_score_mismatch'
    assert result.stats['verification_CPU_seconds']>0


@pytest.mark.parametrize('mode',['checkpoint_crash','checkpoint_burn'])
def test_fully_checked_checkpoint_survives_later_child_failure(tmp_path,witness,mode):
    blocks,parent,candidate,score=witness
    executable,_=fake_worker(tmp_path,mode,candidate,score)
    result=invoke(blocks,parent,executable)
    assert result.order==candidate and result.score==score
    assert result.stats['window_status'].startswith('verified_checkpoint_after_')
    assert result.stats['child_process_CPU_seconds']>0


def test_owned_cpu_receipt_excludes_an_unrelated_child(tmp_path,witness):
    blocks,parent,_,_=witness
    executable,_=fake_worker(tmp_path,'burn')
    # A separately reaped child makes global RUSAGE_CHILDREN deltas unsuitable.
    other=subprocess.Popen([sys.executable,'-c','import time; a=time.process_time();\nwhile time.process_time()-a<.12: pass'])
    other.wait()
    result=invoke(blocks,parent,executable)
    charged=result.stats['child_process_CPU_seconds']
    assert .05<=charged<.16
    assert result.CPU_seconds>=charged


def test_startup_watchdog_kills_and_reaps_owned_child(tmp_path,witness):
    blocks,parent,_,_=witness
    executable,record=fake_worker(tmp_path,'startup_burn')
    result=invoke(blocks,parent,executable)
    assert result.order==parent and result.score is not None
    assert result.stats['window_status'] in ('child_cpu_limit','whole_call_limit','child_wall_limit')
    assert result.stats['child_returncode']==-9
    assert result.CPU_seconds<2.20 # sampling/kill/cleanup tolerance, no cap extension
    with pytest.raises(ProcessLookupError): os.kill(json.loads(record.read_text())['pid'],0)


def test_final_verification_deadline_preserves_native_parent(tmp_path,witness,monkeypatch):
    blocks,parent,candidate,score=witness
    executable,_=fake_worker(tmp_path,candidate=candidate,score=score)
    original=native.SparseResponse.score
    def score_with_failure(self,order,budget,deadline=float('inf')):
        if tuple(order)==candidate:
            raise native.SearchDeadline('verification interrupted')
        return original(self,order,budget,deadline)
    monkeypatch.setattr(native.SparseResponse,'score',score_with_failure)
    result=invoke(blocks,parent,executable)
    assert result.order==parent and result.score is not None
    assert result.stats['window_status']=='whole_call_limit'


def test_context_setup_deadline_preserves_completed_parent(tmp_path,witness,monkeypatch):
    blocks,parent,_,_=witness
    executable,_=fake_worker(tmp_path)
    original=native.SparseResponse.__init__;calls=0
    def setup(self,*args,**kwargs):
        nonlocal calls
        calls+=1
        if calls==2: raise native.SearchDeadline('second context setup interrupted')
        original(self,*args,**kwargs)
    monkeypatch.setattr(native.SparseResponse,'__init__',setup)
    result=invoke(blocks,parent,executable)
    assert result.order==parent and result.score is not None
    assert result.stats['window_status']=='whole_call_limit'
    assert not native._EXTERNAL_BUDGET_CHECK.get()


def test_no_child_is_launched_without_a_measured_verification_reserve(tmp_path,witness,monkeypatch):
    blocks,parent,_,_=witness
    original=native.SparseResponse.score;calls=0
    def measured_slow_score(self,*args,**kwargs):
        nonlocal calls
        calls+=1
        result=original(self,*args,**kwargs)
        if calls==2:
            began=time.process_time()
            while time.process_time()-began<.05: pass
        return result
    def forbidden(*args,**kwargs): raise AssertionError('worker launched without reserve')
    monkeypatch.setattr(native.SparseResponse,'score',measured_slow_score)
    monkeypatch.setattr(window.subprocess,'Popen',forbidden)
    result=native.train_ordering(blocks,3,2,time=.37,rotation_budget=80,cpu_seconds=.10,
        max_evaluations=100,max_moves=0,starts=[parent],method='window',
        window_options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=sys.executable))
    assert result.order==parent and result.score is not None
    assert result.stats['window_status']=='verification_reserve_unavailable'


def test_no_completed_parent_is_never_zero(tmp_path,witness):
    blocks,parent,_,_=witness
    result=native.train_ordering(blocks,3,2,time=.37,rotation_budget=80,
        max_evaluations=0,starts=[parent],method='window')
    assert result.order==parent and result.score is None
    assert result.stats['window_status']=='parent_unscored'


def test_evaluation_allocation_formula_keeps_existing_fraction_and_absolute_cap():
    for total,expected in ((0,0),(1,1),(5,5),(6,1),(7,2),(8,3),(12,6),(20,10),(3000,1798)):
        assert window._parent_evaluation_cap(total,.60)==expected
    assert window._parent_evaluation_cap(6,.01)==1
    assert window._parent_evaluation_cap(12,.20)==2
    for total in range(6,80):
        parent=window._parent_evaluation_cap(total,.60)
        assert 1<=parent<=total-2-3
        assert parent+2+3<=total


@pytest.mark.parametrize('total',range(6))
def test_insufficient_evaluation_pool_keeps_real_parent_and_never_launches(tmp_path,total,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('insufficient evaluation pool launched a child')
    monkeypatch.setattr(window.subprocess,'Popen',forbidden)
    parent=(5,3,1,4,2,0)
    result=native.train_ordering(tiny_blocks(),3,2,time=Fraction(37,100),rotation_budget=80,
        cpu_seconds=2,max_evaluations=total,max_moves=0,starts=[parent],method='window',
        window_options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=sys.executable))
    assert result.order==parent
    assert result.evaluations==(1 if total else 0)
    assert (result.score is None)==(total==0)
    stats=result.stats
    assert stats['requested_evaluations']==total
    assert stats['parent_allocated_evaluations']==total
    assert stats['parent_evaluations']==result.evaluations<=total
    assert stats['native_reserved_evaluations']==0
    assert stats['native_parent_revalidation_evaluations']==0
    assert stats['native_final_reserved_evaluations']==0
    assert stats['child_allocated_evaluations']==0
    assert stats['child_completed_evaluations'] is None
    assert stats['window_status']==('parent_unscored' if total==0 else 'insufficient_evaluation_budget')


def test_real_sparse_parent_saturates_reserved_cap_then_real_native_checks_changed_reply(tmp_path,monkeypatch):
    blocks=tiny_blocks();parent=(5,3,1,4,2,0)
    executable,record=fake_worker(tmp_path,'budget_changed')
    original=native.SparseResponse.score
    native_calls=[]
    def recorded_score(self,order,budget,deadline=float('inf')):
        native_calls.append((id(self),tuple(order)))
        return original(self,order,budget,deadline)
    monkeypatch.setattr(native.SparseResponse,'score',recorded_score)
    # For this six-unit fixture the real first native insertion sweep has30
    # trials; its six-evaluation cap is reached before that sweep completes.
    # No TrainingResult/evaluation count is manufactured by the test.
    result=native.train_ordering(blocks,3,2,time=Fraction(37,100),rotation_budget=80,
        cpu_seconds=3,max_evaluations=12,max_moves=10,starts=[parent],method='window',
        window_options=window.WindowOptions(width=4,parent_method='sparse',julia_executable=str(executable)))
    request=json.loads(record.read_text())['request']
    stats=result.stats
    assert stats['parent_allocated_evaluations']==6
    assert stats['parent_evaluations']==6 and stats['parent_evaluation_cap_saturated']
    assert stats['parent_status']=='evaluation_limit'
    assert stats['native_reserved_evaluations']==2
    assert stats['native_parent_revalidation_evaluations']==1
    assert stats['native_final_reserved_evaluations']==1
    assert request['max_evaluations']==stats['child_allocated_evaluations']==4
    assert stats['child_completed_evaluations']==4
    assert stats['parent_evaluations']+stats['native_parent_revalidation_evaluations']+request['max_evaluations']+stats['native_final_reserved_evaluations']==12
    assert result.evaluations==12
    assert stats['native_verification_attempts']==stats['native_verifications']==1
    assert len(native_calls)==8 # Six real parent scores, revalidation, final.
    expected=list(request['order']);expected[1],expected[2]=expected[2],expected[1]
    assert native_calls[-1][1]==tuple(expected)
    assert native_calls[-2][1]==tuple(request['order'])
    assert expected[0]==request['order'][0] and expected[-1]==request['order'][-1]
    assert stats['window_status'] in ('native_julia_score_mismatch','native_no_resolved_improvement','verified_improvement')
    assert stats['child_evaluation_scope']=='checkpoint_prefix'
    assert 'uncheckpointed' in stats['evaluation_accounting_scope']
    # Recheck the retained result independently; this extra fixture check is
    # outside the construction call and is not attributed to its evaluation cap.
    check=native.SparseResponse(blocks,3,2,time=Fraction(37,100)).score(result.order,80)
    assert result.score==check


@pytest.mark.parametrize('total,parent_cap',[(6,1),(20,10),(3000,1798)])
def test_early_parent_releases_unused_evaluations_to_saturated_child_prefix(tmp_path,total,parent_cap):
    parent=(5,3,1,4,2,0)
    executable,record=fake_worker(tmp_path,'budget_saturated')
    result=native.train_ordering(tiny_blocks(),3,2,time=Fraction(37,100),rotation_budget=80,
        cpu_seconds=2,max_evaluations=total,max_moves=0,starts=[parent],method='window',
        window_options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=str(executable)))
    stats=result.stats;request=json.loads(record.read_text())['request']
    assert stats['parent_allocated_evaluations']==parent_cap
    assert stats['parent_evaluations']==1
    assert stats['parent_allocation_fraction']==.60
    assert 0<stats['parent_allocated_CPU_seconds']<=2*.60
    assert stats['native_parent_revalidation_evaluations']==1
    assert request['max_evaluations']==stats['child_allocated_evaluations']==total-3
    assert request['max_evaluations']>=3 # Parent full score, flat baseline, first trial.
    assert stats['child_completed_evaluations']==request['max_evaluations']
    assert stats['parent_evaluations']+stats['native_parent_revalidation_evaluations']+request['max_evaluations']+stats['native_final_reserved_evaluations']==total
    assert result.evaluations==total-1 # Unchanged order leaves final reserve unused.
    assert stats['native_verification_attempts']==0
    assert result.order==parent and result.score is not None
    assert stats['child_evaluation_scope']=='checkpoint_prefix'


def test_child_reply_cannot_extend_the_serialized_evaluation_allocation(tmp_path):
    executable,record=fake_worker(tmp_path,'budget_excess')
    result=native.train_ordering(tiny_blocks(),3,2,time=Fraction(37,100),rotation_budget=80,
        cpu_seconds=2,max_evaluations=8,max_moves=0,starts=[(5,3,1,4,2,0)],method='window',
        window_options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=str(executable)))
    request=json.loads(record.read_text())['request'];stats=result.stats
    assert request['max_evaluations']==stats['child_allocated_evaluations']==5
    assert stats['window_status']=='invalid_worker_reply'
    assert 'allocated evaluation cap' in stats['refusal_reason']
    assert stats['child_completed_evaluations'] is None
    assert stats['native_verification_attempts']==0
    assert result.evaluations==2 and result.score is not None
    assert result.evaluations<=8


def test_normalized_underflow_is_unscored_instead_of_completed_zero():
    blocks=({((0,'X'),):.5},{((0,'Z'),):.3})
    result=native.train_ordering(blocks,1,0,time=1e-150,rotation_budget=10**18,
        method='window',max_moves=0,window_options=window.WindowOptions(width=2,parent_method='magnitude'))
    assert result.score is None and result.status=='numerical_limit'
    assert result.stats['window_status']=='numerical_limit'


def test_missing_julia_and_unsupported_domains_are_explicit(tmp_path,witness):
    blocks,parent,_,_=witness
    result=invoke(blocks,parent,tmp_path/'does-not-exist')
    assert result.order==parent and result.score is not None
    assert result.stats['window_status']=='julia_unavailable'
    coherent={1:math.sqrt(.5),2:math.sqrt(.5)}
    result=native.train_ordering(blocks,3,coherent,rotation_budget=80,method='window',max_moves=0,
        window_options=window.WindowOptions(width=4,parent_method='magnitude'))
    assert result.stats['window_status']=='unsupported_coherent_reference'
    mixed=list(blocks);mixed[1]={((0,'X'),):.5,((1,'X'),):.5}
    result=native.train_ordering(mixed,3,2,rotation_budget=80,method='window',max_moves=0,
        window_options=window.WindowOptions(width=4,parent_method='magnitude'))
    assert result.stats['window_status']=='unsupported_complex_or_mixed_flip_unit'
    with pytest.raises(ValueError,match='sparse or magnitude'):
        window.WindowOptions(parent_method='current')
    with pytest.raises(ValueError,match='afford'):
        native.train_ordering(blocks,3,2,rotation_budget=0,method='window')


def test_width_is_not_silently_shrunk_and_default_never_launches(tmp_path,witness,monkeypatch):
    blocks,parent,_,_=witness
    def forbidden(*args,**kwargs): raise AssertionError('unexpected Julia process')
    monkeypatch.setattr(window.subprocess,'Popen',forbidden)
    result=native.train_ordering(blocks,3,2,rotation_budget=80,max_moves=0,max_evaluations=1)
    assert result.score is not None
    result=native.train_ordering(blocks,3,2,rotation_budget=80,max_moves=0,method='window',
        window_options=window.WindowOptions(width=8,parent_method='magnitude'))
    assert result.stats['window_status']=='window_unavailable'


def test_normal_cli_window_selection_and_schema(tmp_path,witness,capsys):
    blocks,_,_,_=witness
    # Common denominator only belongs to the existing model JSON boundary.
    model=dict(n=3,electrons=1,denominator='100',blocks=[
        [[list(map(list,word)),str(c.numerator*(100//c.denominator))] for word,c in block.items()]
        for block in blocks])
    path=tmp_path/'model.json';path.write_text(json.dumps(model))
    assert main([str(path),'--method','window','--window-width','4','--window-parent','magnitude',
        '--time','1.6','--rotation-budget','41','--julia-executable',str(tmp_path/'missing')])==0
    data=json.loads(capsys.readouterr().out)
    assert data['training']['stats']['window_status']=='julia_unavailable'
    assert data['training']['score'] is not None
    assert data['time']=='8/5' and data['rotation_budget']==41
    assert 'not a finite-error bound' in data['score_scope']


@pytest.mark.parametrize('duration,budget',[(Fraction(37,100),17),(Fraction(8,5),32)])
def test_real_julia_normal_api_and_native_parity(duration,budget):
    julia=os.environ.get('QHAT_JULIA_EXECUTABLE') or shutil.which('julia')
    if not julia: pytest.skip('pre-existing Julia binary unavailable')
    blocks=tiny_blocks();parent=(5,3,1,4,2,0)
    result=native.train_ordering(blocks,3,2,time=duration,rotation_budget=budget,
        starts=[parent],max_moves=0,max_evaluations=300,cpu_seconds=12,method='window',
        window_options=window.WindowOptions(width=4,parent_method='magnitude',julia_executable=julia,wall_seconds=48))
    assert result.stats['window_status'] in ('complete','verified_improvement','native_no_resolved_improvement')
    assert result.stats['child_returncode']==0
    context=native.SparseResponse(blocks,3,2,time=duration)
    assert result.score.value==pytest.approx(context.score(result.order,budget).value,rel=1e-11,abs=1e-28)
    assert result.order[0]==parent[0] and result.order[-1]==parent[-1]
    assert result.CPU_seconds>result.stats['child_process_CPU_seconds']>0


def test_installed_wheel_normal_cli_and_asset_discovery(tmp_path,witness):
    wheel_path=os.environ.get('QHAT_WINDOW_WHEEL')
    if not wheel_path: pytest.skip('set QHAT_WINDOW_WHEEL to an exact-source --no-deps wheel')
    installed=tmp_path/'installed';installed.mkdir()
    subprocess.run([sys.executable,'-m','pip','install','--no-deps','--no-index','--target',str(installed),wheel_path],check=True,capture_output=True,text=True)
    blocks,_,_,_=witness
    model=dict(n=3,electrons=1,denominator='100',blocks=[
        [[list(map(list,word)),str(c.numerator*(100//c.denominator))] for word,c in block.items()]
        for block in blocks])
    model_path=tmp_path/'model.json';model_path.write_text(json.dumps(model))
    executable,record=fake_worker(tmp_path)
    env=os.environ.copy();env['PYTHONPATH']=str(installed);env['PYTHONDONTWRITEBYTECODE']='1'
    completed=subprocess.run([sys.executable,'-S','-m','qhat.analysis.train_ordering',str(model_path),
        '--method','window','--window-width','4','--window-parent','magnitude',
        '--time','.37','--rotation-budget','80','--julia-executable',str(executable)],
        cwd=tmp_path,env=env,check=True,capture_output=True,text=True)
    result=json.loads(completed.stdout)
    assert result['training']['stats']['window_status']=='complete'
    argv=json.loads(record.read_text())['argv']
    worker=Path(argv[-5])
    assert worker.is_relative_to(installed)
    for asset in ('worker.jl','ResponseCore.jl','WindowSearch.jl','Project.toml','SOURCE.toml'):
        assert (worker.parent/asset).is_file()
    assert 'not a finite-error bound' in result['score_scope']
    julia=os.environ.get('QHAT_JULIA_EXECUTABLE') or shutil.which('julia')
    if julia:
        # Exercise actual Julia from the installed package and the real CLI,
        # outside the source checkout and without its import/site machinery.
        actual=subprocess.run([sys.executable,'-S','-m','qhat.analysis.train_ordering',str(model_path),
            '--method','window','--window-width','4','--window-parent','magnitude',
            '--time','.37','--rotation-budget','80','--julia-executable',str(julia),
            '--cpu-seconds','12','--window-wall-seconds','48','--max-evaluations','300'],
            cwd=tmp_path,env=env,check=True,capture_output=True,text=True)
        trained=json.loads(actual.stdout)['training']
        assert trained['stats']['window_status'] in ('complete','verified_improvement','native_no_resolved_improvement')
        assert trained['stats']['child_returncode']==0
        expected=native.SparseResponse(blocks,3,1,time=.37).score(trained['order'],80)
        assert trained['score']['value']==pytest.approx(expected.value,rel=1e-11,abs=1e-28)
