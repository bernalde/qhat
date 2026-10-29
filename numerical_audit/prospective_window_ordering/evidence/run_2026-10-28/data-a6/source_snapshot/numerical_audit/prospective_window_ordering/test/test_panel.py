"""Panel authentication/serialization only: synthetic records, no chemistry/kernels.

The source-map fixture hashes real source files and mocks only the map inventory;
no scientific source, shared freeze, environment or existing file is modified.
Guard receipts are explicit fake metadata, never measured resource evidence.
"""
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import tomllib

import pytest

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import panel


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    return dict(path=str(path),sha256=sha(path))


@pytest.fixture
def bundle(tmp_path,monkeypatch):
    root=tmp_path/'artifacts';root.mkdir()
    registry=json.loads((HERE/'REGISTRY.json').read_text())
    source_files=[HERE/'REGISTRY.json',HERE/'panel.py',HERE/'study_common.py']
    source_map={str(path.relative_to(panel.k.ROOT)):sha(path) for path in source_files}
    monkeypatch.setattr(panel.k,'source_map',lambda:dict(source_map))
    def ref(path): return dict(path=str(path.relative_to(root)),sha256=sha(path))
    def put(path,data): write_json(path,data);return ref(path)
    freeze_path=root/'SOURCE_FROZEN.json'
    freeze=put(freeze_path,dict(source_commit='a'*40,sources=source_map,registry=registry))
    def model_files(directory,data):
        directory.mkdir(parents=True,exist_ok=True)
        json_path=directory/'model.json';put(json_path,data)
        # Literal original block order is 10Z,2Z. Global Pauli order is 2Z,10Z.
        text='\n'.join([
            'case='+json.dumps(data['case']),'n=14','electrons=4',
            'partition="conserving"','denominator='+json.dumps(data['denominator']),
            'source_json_sha256='+json.dumps(sha(json_path)),'baseline_order=[1,2]',
            '[[units]]','words=["10Z"]','numerators='+json.dumps([data['blocks'][0][0][1]]),
            '[[units]]','words=["2Z"]','numerators='+json.dumps([data['blocks'][1][0][1]])])+'\n'
        toml_path=directory/'model.toml';toml_path.write_text(text)
        return ref(json_path),ref(toml_path)
    original={};successes={}
    for spec in registry['cases']:
        directory=root/'inputs'/spec['case'];directory.mkdir(parents=True)
        tensors=directory/'tensors.npz';tensors.write_bytes(b'synthetic-opaque-tensor-identity')
        orbitals=directory/'orbitals.npz';orbitals.write_bytes(b'synthetic-opaque-orbital-identity')
        word2=[[[2,'Z']],'3'];word10=[[[10,'Z']],'5']
        data=dict(spec,basis=registry['basis'],unit=registry['unit'],charge=0,spin=0,
            denominator='8',coefficients=[word2,word10],blocks=[[word10],[word2]],
            tensors_sha256=sha(tensors),orbitals_sha256=sha(orbitals))
        model_ref,toml_ref=model_files(directory,data)
        original[spec['case']]=(data,model_ref,toml_ref)
        successes[spec['case']]={'model.json':model_ref['sha256'],'model.toml':toml_ref['sha256'],
            'tensors.npz':sha(tensors),'orbitals.npz':sha(orbitals)}
    parent_path=root/'INPUTS_FROZEN.json'
    parent=put(parent_path,dict(status='complete',expected_cases=4,successes=successes,
                               refusals={},protocol_sha256=freeze['sha256']))
    # Equal word counts make the prescribed lexical tie decisive.
    chosen=min(original)
    before,parent_model,_=original[chosen]
    altered=json.loads((root/parent_model['path']).read_text())
    altered.update(case=chosen+'_preflight_7over5',denominator='40',
        preflight_parent_case=chosen,preflight_parent_model_sha256=parent_model['sha256'],preflight_scale='7/5')
    for block in [altered['coefficients'],*altered['blocks']]:
        for entry in block: entry[1]=str(7*int(entry[1]))
    model_ref,toml_ref=model_files(root/'preflight'/'inputs'/altered['case'],altered)
    input_freeze=put(root/'preflight'/'INPUTS_FROZEN.json',dict(status='complete',expected_cases=1,
        preflight_only=True,parent_inventory_sha256=parent['sha256'],
        successes={altered['case']:{'model.json':model_ref['sha256'],'model.toml':toml_ref['sha256']}}))
    entries=[]
    def receipt(arm,output,*,suffix='a1',cpu=.2,wall=1.,exit_code=0):
        path=root/'logs'/f'{arm}.{suffix}.resources.json'
        command=[sys.executable,str(HERE/'constructors.py'),str(root/model_ref['path']),model_ref['sha256'],
                 arm,str(output),'--repetition','1','--source-freeze',str(freeze_path)]
        body=dict(unit=f'synthetic-{arm}-{suffix}',settings=dict(CPUQuota='25%',CPUWeight='5',IOWeight='10',
            Nice='19',MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0',TasksMax='32',
            MemoryAccounting='yes',CPUAccounting='yes',RuntimeMaxSec='150'),
            environment={key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                            'NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS')},
            command=command,budget=dict(cpu_seconds=30.,wall_seconds=150,disk_mib=16),
            elapsed=wall,final_accounting=dict(CPU_seconds=cpu,elapsed_seconds=wall,
                CPU_accounting_available=True,process_CPU_timing_available=True),
            exit_code=exit_code,supervisor_exit_code=exit_code,termination=None,guard_error=None,
            cleanup_errors=[],frozen_seconds=0.)
        return put(path,body)
    for arm in panel.k.ARMS:
        partition='pauli' if arm=='pauli' else 'conserving'
        # Independent layout fixture: global Pauli labels use numeric tuple
        # sorting, whereas conserving labels keep ORIGINAL 10Z,2Z block order.
        by_label=[(2,21),(10,35)] if partition=='pauli' else [(10,35),(2,21)]
        order=[1,0]
        layout=[[[[by_label[i][0],'Z']],str(Fraction(by_label[i][1],40))] for i in order]
        kind=('retained_first_order_response' if arm in ('first','pauli') else
              'diagonal_reference_cubic_response' if arm in ('cubic','window_sparse','window_magnitude') else None)
        result=dict(order=order,partition=partition,status='complete',score=0. if kind else None,
            evaluations=1 if kind else 0,accepted_moves=0,stats={})
        if kind=='diagonal_reference_cubic_response':
            result['score_components']=dict(value=0.,cubic_response=0.,steps=29)
        if arm.startswith('window'):
            result['stats']=dict(window_status='complete',child_process_CPU_seconds=.01)
        record=dict(schema='prospective-window-constructor-v1',artifact_complete=True,arm=arm,repetition=1,
            case=altered['case'],partition=partition,model_sha256=model_ref['sha256'],
            source_freeze_sha256=freeze['sha256'],training_time=1.,training_rotation_budget=60,
            total_CPU_ceiling=64.,constructor=result,score_kind=kind,original_word_count=2,unit_lengths=[1,1],
            pulse_layout_sha256=hashlib.sha256(json.dumps(layout,separators=(',',':')).encode()).hexdigest(),
            process_CPU_before_output=.1,wall_seconds_before_output=.5,
            physical_error_read=False,response_is_error_bound=False)
        output=root/'constructors'/f'{arm}.r1.a1.json';record_ref=put(output,record)
        receipt_ref=receipt(arm,output)
        entries.append(dict(arm=arm,repetition=1,selected_attempt=1,
                            attempts=[dict(receipt=receipt_ref,record=record_ref)]))
    index=dict(schema='prospective-window-panel-inventory-v1',mode='preflight',case=altered['case'],
        source_freeze=freeze,input_freeze=input_freeze,parent_input_freeze=parent,
        model_json=model_ref,model_toml=toml_ref,execution_root=str(root),execution_source_root=str(panel.k.ROOT),
        constructors=entries)
    index_path=root/'PANEL_INVENTORY.json'
    whole_path=root/'ALL_ORDERS_FROZEN.json'
    whole_control=dict(auto=True)
    def refresh_whole():
        artifacts=[]
        for entry in index['constructors']:
            final=entry['attempts'][-1]
            status=json.loads((root/final['receipt']['path']).read_text())['exit_code']
            artifacts.append(dict(case=index['case'],arm=entry['arm'],repetition=entry['repetition'],
                receipt=final['receipt'],record=final.get('record'),exit_code=status))
        value=dict(schema='prospective-window-all-orders-freeze-v1',mode=index['mode'],
            source_freeze=index['source_freeze'],plan_sha256='b'*64,inventory_complete=True,
            expected_calls=len(artifacts),artifacts=artifacts,
            all_calls_have_artifacts=all(row['record'] is not None for row in artifacts),
            physical_error_read=False,no_selected_repetition=True)
        index['all_orders_freeze']=put(whole_path,value)
    refresh_whole()
    def publish():
        if whole_control['auto']: refresh_whole()
        put(index_path,index)
        return index_path
    return dict(root=root,index=index,publish=publish,put=put,ref=ref,receipt=receipt,
                registry=registry,original=original,parent_path=parent_path,
                sources=source_map,output=tmp_path/'panel',whole_control=whole_control,
                whole_path=whole_path,refresh_whole=refresh_whole)


def assemble(bundle):
    return panel.assemble(bundle['publish'](),bundle['output'])


def test_complete_roster_original_labels_numeric_pauli_and_scored_zero(bundle):
    report=assemble(bundle)
    assert report['status']=='complete' and report['candidate_panel_complete']
    rows=tomllib.loads((bundle['output']/'CANDIDATES.toml').read_text())['candidates']
    assert len(rows)==44 and len({row['name'] for row in rows})==44
    assert {tuple(row['order']) for row in rows}=={(2,1)}
    pauli=[row for row in rows if row['arm']=='pauli']
    assert len(pauli)==4 and all(row['partition']=='pauli' for row in pauli)
    retained=next(row for row in report['constructors'] if row['arm']=='pauli')
    record=json.loads((bundle['root']/retained['selected_record']['path']).read_text())
    expected=[[[[10,'Z']],'7/8'],[[[2,'Z']],'21/40']]
    assert record['pulse_layout_sha256']==hashlib.sha256(json.dumps(expected,separators=(',',':')).encode()).hexdigest()
    assert all(row['constructor_scored'] and row['constructor_score']==0. for row in pauli)
    assert all('constructor_score' not in row and not row['constructor_scored']
               for row in rows if row['arm']=='canonical')
    assert report['guard_CPU_seconds']==pytest.approx(8*.2)
    assert report['candidates']['sha256']==sha(bundle['output']/'CANDIDATES.toml')
    assert report['all_orders_freeze']['inventory_complete']
    assert report['all_orders_freeze']['censored_call_count']==0
    assert json.loads((bundle['output']/'PANEL.json').read_text())['status']=='complete'
    with pytest.raises(ValueError,match='Fresh'):
        assemble(bundle)


def test_failed_attempt_costs_remain_charged_once_not_per_physical_row(bundle):
    entry=bundle['index']['constructors'][0]
    missing=bundle['root']/'constructors'/'canonical.failed.json'
    failed=bundle['receipt']('canonical',missing,suffix='failed',cpu=.15,wall=.8,exit_code=2)
    entry['attempts'].insert(0,dict(receipt=failed,record=None,censor_reason='process exited before final record'))
    entry['selected_attempt']=2
    report=assemble(bundle)
    assert report['status']=='complete'
    arm=report['arm_totals']['canonical']
    assert arm['CPU_seconds']==pytest.approx(.35)
    assert arm['wall_seconds']==pytest.approx(1.8)
    assert report['guard_CPU_seconds']==pytest.approx(8*.2+.15)
    retained=next(row for row in report['constructors'] if row['arm']=='canonical')
    assert retained['attempts'][0]['exit_code']==2
    assert retained['attempts'][0]['status']=='censored_no_constructor_record'
    assert retained['selected_attempt']==2


def test_censor_never_becomes_reduced_candidate_panel(bundle):
    entry=next(row for row in bundle['index']['constructors'] if row['arm']=='cubic')
    (bundle['root']/entry['attempts'][0]['record']['path']).unlink()
    entry['attempts'][0]['record']=None;entry['selected_attempt']=None
    entry['censor_reason']='constructor never produced a complete original-label order'
    report=assemble(bundle)
    assert report['status']=='incomplete' and not report['candidate_panel_complete']
    assert len(report['constructors'])==8 and report['available_logical_rows']==40
    assert report['arm_totals']['cubic']['CPU_seconds']==pytest.approx(.2)
    assert not (bundle['output']/'CANDIDATES.toml').exists()
    assert (bundle['output']/'PANEL.json').exists()


@pytest.mark.parametrize('tamper',['order','hash','source'])
def test_tampered_order_artifact_hash_or_frozen_source_refuses(bundle,monkeypatch,tamper):
    if tamper=='source':
        changed=dict(bundle['sources']);changed['numerical_audit/prospective_window_ordering/panel.py']='f'*64
        monkeypatch.setattr(panel.k,'source_map',lambda:dict(changed))
    else:
        entry=next(row for row in bundle['index']['constructors'] if row['arm']=='cubic')
        descriptor=entry['attempts'][0]['record']
        if tamper=='hash':
            descriptor['sha256']='0'*64
        else:
            path=bundle['root']/descriptor['path'];record=json.loads(path.read_text())
            record['constructor']['order']=[0,0]
            entry['attempts'][0]['record']=bundle['put'](path,record)
    report=assemble(bundle)
    assert report['status']=='incomplete' and not report['candidate_panel_complete']
    assert not (bundle['output']/'CANDIDATES.toml').exists()
    if tamper=='order':
        assert report['arm_totals']['cubic']['CPU_seconds']==pytest.approx(.2)
        retained=next(row for row in report['constructors'] if row['arm']=='cubic')
        assert retained['attempts'][0]['constructor']==entry['attempts'][0]['record']


def test_unknown_receipt_cpu_is_unknown_not_completed_zero(bundle):
    entry=bundle['index']['constructors'][0];descriptor=entry['attempts'][0]['receipt']
    path=bundle['root']/descriptor['path'];receipt=json.loads(path.read_text())
    receipt['final_accounting'].update(CPU_seconds=None,CPU_accounting_available=False,process_CPU_timing_available=False)
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    report=assemble(bundle)
    assert report['status']=='incomplete' and report['guard_CPU_seconds'] is None
    assert report['arm_totals']['canonical']['CPU_seconds'] is None
    assert not report['all_guard_accounting_complete']


def test_wrong_pauli_string_sort_is_not_authenticated(bundle):
    entry=next(row for row in bundle['index']['constructors'] if row['arm']=='pauli')
    descriptor=entry['attempts'][0]['record'];path=bundle['root']/descriptor['path']
    record=json.loads(path.read_text())
    # This is the erroneous unpadded lexical 10Z,2Z Pauli labeling under order1,0.
    wrong=[[[[2,'Z']],'21/40'],[[[10,'Z']],'7/8']]
    record['pulse_layout_sha256']=hashlib.sha256(json.dumps(wrong,separators=(',',':')).encode()).hexdigest()
    entry['attempts'][0]['record']=bundle['put'](path,record)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']


def test_reused_guard_receipt_is_not_charged_as_another_arm(bundle):
    entries=bundle['index']['constructors']
    entries[1]['attempts'][0]['receipt']=copy.deepcopy(entries[0]['attempts'][0]['receipt'])
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert report['arm_totals']['magnitude']['CPU_seconds'] is None


def test_literal_env_prefix_and_campaign_tmpdir_are_supported(bundle):
    for entry in bundle['index']['constructors']:
        descriptor=entry['attempts'][0]['receipt'];path=bundle['root']/descriptor['path']
        receipt=json.loads(path.read_text())
        receipt['command']=['/usr/bin/env','TMPDIR='+str(bundle['root']/'tmp'),
            'JULIA_NUM_THREADS=1','JULIA_NUM_GC_THREADS=1','JULIA_PKG_OFFLINE=true',*receipt['command']]
        entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    assert assemble(bundle)['candidate_panel_complete']


@pytest.mark.parametrize('failure',['missing','wrong_source','bad_artifact','incomplete_roster','selected_descriptor'])
def test_whole_order_freeze_is_required_and_hash_bound_without_losing_costs(bundle,failure):
    bundle['whole_control']['auto']=False
    index=bundle['index'];path=bundle['whole_path']
    whole=json.loads(path.read_text())
    if failure=='missing':
        index.pop('all_orders_freeze')
    else:
        if failure=='wrong_source':
            whole['source_freeze']=dict(index['source_freeze'],sha256='0'*64)
        elif failure=='bad_artifact':
            whole['artifacts'][0]['receipt']=dict(whole['artifacts'][0]['receipt'],sha256='0'*64)
        elif failure=='incomplete_roster':
            whole['artifacts'].pop()
        else:
            # Identical valid bytes under a different path still do not match
            # the selected descriptor frozen for this current case.
            row=whole['artifacts'][0]
            original=json.loads((bundle['root']/row['record']['path']).read_text())
            row['record']=bundle['put'](bundle['root']/'copied-record.json',original)
        index['all_orders_freeze']=bundle['put'](path,whole)
    report=assemble(bundle)
    assert not report['candidate_panel_complete'] and report['status']=='incomplete'
    assert report['all_guard_accounting_complete']
    assert report['guard_CPU_seconds']==pytest.approx(8*.2)
    assert not (bundle['output']/'CANDIDATES.toml').exists()
    assert any(reason.get('phase')=='all_orders_freeze' for reason in report['failures'])


def main_bundle(bundle):
    root=bundle['root'];index=bundle['index'];templates=copy.deepcopy(index['constructors'])
    parent_ref=index['parent_input_freeze']
    parent=json.loads((root/parent_ref['path']).read_text())
    current=min(parent['successes'])
    current_entries=[];artifacts=[]
    for case in sorted(parent['successes']):
        model_path=root/'inputs'/case/'model.json'
        for repetition in (1,2):
            for template in templates:
                arm=template['arm'];before=template['attempts'][0]
                original=json.loads((root/before['record']['path']).read_text())
                original.update(case=case,repetition=repetition,model_sha256=sha(model_path))
                labels=[(2,3),(10,5)] if arm=='pauli' else [(10,5),(2,3)]
                layout=[[[[labels[i][0],'Z']],str(Fraction(labels[i][1],8))]
                        for i in original['constructor']['order']]
                original['pulse_layout_sha256']=hashlib.sha256(json.dumps(layout,separators=(',',':')).encode()).hexdigest()
                record_path=root/'main-records'/case/f'{arm}.r{repetition}.json'
                record_ref=bundle['put'](record_path,original)
                receipt=json.loads((root/before['receipt']['path']).read_text())
                receipt['unit']=f'synthetic-main-{case}-{arm}-r{repetition}'
                command=receipt['command'];script=command.index(str(HERE/'constructors.py'))
                command[script+1]=str(model_path);command[script+2]=sha(model_path)
                command[script+4]=str(record_path);command[command.index('--repetition')+1]=str(repetition)
                receipt_ref=bundle['put'](root/'main-logs'/case/f'{arm}.r{repetition}.json',receipt)
                row=dict(case=case,arm=arm,repetition=repetition,receipt=receipt_ref,record=record_ref,exit_code=0)
                artifacts.append(row)
                if case==current:
                    current_entries.append(dict(arm=arm,repetition=repetition,selected_attempt=1,
                        attempts=[dict(receipt=receipt_ref,record=record_ref)]))
    index.update(mode='main',case=current,input_freeze=parent_ref,constructors=current_entries,
        model_json=bundle['ref'](root/'inputs'/current/'model.json'),
        model_toml=bundle['ref'](root/'inputs'/current/'model.toml'))
    index.pop('parent_input_freeze')
    whole=dict(schema='prospective-window-all-orders-freeze-v1',mode='main',
        source_freeze=index['source_freeze'],plan_sha256='b'*64,inventory_complete=True,
        expected_calls=64,artifacts=artifacts,all_calls_have_artifacts=True,
        physical_error_read=False,no_selected_repetition=True)
    index['all_orders_freeze']=bundle['put'](bundle['whole_path'],whole)
    bundle['whole_control']['auto']=False
    return whole


def test_main_whole_freeze_preserves_other_case_censor_without_campaign_completion_claim(bundle):
    whole=main_bundle(bundle);root=bundle['root'];index=bundle['index'];current=index['case']
    artifacts=whole['artifacts']
    censored=next(row for row in artifacts if row['case']!=current)
    (root/censored['record']['path']).unlink()
    censored['record']=None;censored['censor_reason']='guarded constructor produced no final record'
    whole['all_calls_have_artifacts']=False
    index['all_orders_freeze']=bundle['put'](bundle['whole_path'],whole)
    bundle['whole_control']['auto']=False
    report=assemble(bundle)
    assert report['status']=='complete' and report['candidate_panel_complete']
    assert report['available_logical_rows']==88
    global_check=report['all_orders_freeze']
    assert global_check['expected_calls']==64 and global_check['censored_call_count']==1
    assert global_check['censored_calls'][0]['case']==censored['case']
    assert not global_check['all_calls_have_artifacts']
    candidate=tomllib.loads((bundle['output']/'CANDIDATES.toml').read_text())
    assert not candidate['whole_constructor_artifacts_present']
    assert candidate['all_orders_freeze_sha256']==index['all_orders_freeze']['sha256']


def test_all_four_parent_tomls_are_authenticated_before_preflight_selection(bundle):
    index=bundle['index'];root=bundle['root'];path=bundle['parent_path']
    frozen=json.loads(path.read_text());other=max(frozen['successes'])
    toml_path=path.parent/'inputs'/other/'model.toml'
    toml_path.write_text(toml_path.read_text().replace('10Z','9Z'))
    frozen['successes'][other]['model.toml']=sha(toml_path)
    index['parent_input_freeze']=bundle['put'](path,frozen)
    altered_path=root/index['input_freeze']['path'];altered=json.loads(altered_path.read_text())
    altered['parent_inventory_sha256']=index['parent_input_freeze']['sha256']
    index['input_freeze']=bundle['put'](altered_path,altered)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert any('TOML original words' in f['message'] for f in report['failures'])
    assert (bundle['output']/'PANEL.json').is_file()


def test_frozen_case_inventory_key_cannot_be_another_registry_case(bundle):
    root=bundle['root'];path=bundle['parent_path'];frozen=json.loads(path.read_text())
    first,last=sorted(frozen['successes'])[::3]
    source=path.parent/'inputs'/last/'model.json';target=path.parent/'inputs'/first/'model.json'
    target.write_bytes(source.read_bytes());frozen['successes'][first]['model.json']=sha(target)
    with pytest.raises(ValueError,match='inventory key/model case'):
        panel.frozen_case(root,path,frozen,first,bundle['registry'])


@pytest.mark.parametrize('filename',['model.json','model.toml','tensors.npz','orbitals.npz'])
def test_every_frozen_case_payload_refuses_symlink_escape(bundle,tmp_path,filename):
    root=bundle['root'];path=bundle['parent_path'];frozen=json.loads(path.read_text())
    case=max(frozen['successes']);target=path.parent/'inputs'/case/filename
    outside=tmp_path/('external-'+filename);outside.write_bytes(target.read_bytes())
    target.unlink();target.symlink_to(outside)
    with pytest.raises(ValueError,match='symlink escape'):
        panel.frozen_case(root,path,frozen,case,bundle['registry'])


@pytest.mark.parametrize('tamper',['record','profile','command'])
def test_main_other_case_receipt_and_full_constructor_semantics_are_required(bundle,tamper):
    whole=main_bundle(bundle);root=bundle['root'];index=bundle['index']
    row=next(r for r in whole['artifacts'] if r['case']!=index['case'])
    descriptor=row['record'] if tamper=='record' else row['receipt']
    path=root/descriptor['path'];data=json.loads(path.read_text())
    if tamper=='record': data.pop('constructor')
    elif tamper=='profile': data['settings']['CPUQuota']='50%'
    else: data['command'][data['command'].index(row['arm'])]='cubic' if row['arm']!='cubic' else 'canonical'
    row['record' if tamper=='record' else 'receipt']=bundle['put'](path,data)
    index['all_orders_freeze']=bundle['put'](bundle['whole_path'],whole)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert any(f.get('phase')=='all_orders_freeze' for f in report['failures'])
    assert report['all_guard_accounting_complete']


def test_whole_roster_cannot_reuse_a_current_nonfinal_service(bundle):
    whole=main_bundle(bundle);root=bundle['root'];index=bundle['index']
    other=next(r for r in whole['artifacts'] if r['case']!=index['case'])
    other_receipt=json.loads((root/other['receipt']['path']).read_text())
    entry=index['constructors'][0];source=entry['attempts'][0]['receipt']
    failed=json.loads((root/source['path']).read_text())
    failed.update(unit=other_receipt['unit'],exit_code=2,supervisor_exit_code=2)
    failed['command'][failed['command'].index(str(HERE/'constructors.py'))+4]=str(root/'failed-distinct.json')
    failed_ref=bundle['put'](root/'earlier.resources.json',failed)
    entry['attempts'].insert(0,dict(receipt=failed_ref,record=None,censor_reason='failed earlier call'))
    entry['selected_attempt']=2
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert any(f.get('phase')=='all_orders_freeze' and 'unit' in f['message'].lower()
               for f in report['failures'])


def test_not_launched_is_an_unperformed_audit_not_a_free_completed_method(bundle):
    entry=bundle['index']['constructors'][0]
    entry.update(attempts=[],selected_attempt=None,not_launched=True,censor_reason='never dispatched')
    bundle['whole_control']['auto']=False
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    row=report['constructors'][0]
    assert row['status']=='unperformed' and not row['cost_complete']
    assert row['guard_CPU_seconds'] is None and row['guard_wall_seconds'] is None
    assert any('actual receipt required' in f['message'] for f in report['failures'])
    assert (bundle['output']/'PANEL.json').is_file()


def test_stale_earlier_output_cannot_be_selected_from_a_later_failed_attempt(bundle):
    entry=bundle['index']['constructors'][0];first=entry['attempts'][0]
    root=bundle['root'];destination=root/first['record']['path']
    second=bundle['receipt']('canonical',destination,suffix='a2',cpu=.1,wall=1.,exit_code=2)
    entry['attempts'].append(dict(receipt=second,record=copy.deepcopy(first['record'])))
    entry['selected_attempt']=2
    report=assemble(bundle)
    retained=report['constructors'][0]
    assert not report['candidate_panel_complete']
    assert retained['known_guard_CPU_seconds']==pytest.approx(.3)
    assert 'destination' in retained['attempts'][1]['failure']['message']
    assert destination.is_file() and retained['attempts'][1]['constructor']==first['record']


@pytest.mark.parametrize('field,value',[('exit_code',2),('supervisor_exit_code',2),('termination','cpu_limit')])
def test_outer_failure_retains_record_and_cost_but_not_an_admitted_order(bundle,field,value):
    entry=bundle['index']['constructors'][0];path=bundle['root']/entry['attempts'][0]['receipt']['path']
    receipt=json.loads(path.read_text());receipt[field]=value
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    report=assemble(bundle);retained=report['constructors'][0]
    assert not report['candidate_panel_complete']
    assert retained['attempts'][0]['status']=='failed_outer_record_retained'
    assert retained['guard_CPU_seconds']==pytest.approx(.2)
    assert (bundle['root']/entry['attempts'][0]['record']['path']).is_file()


@pytest.mark.parametrize('key',['PYTHONPATH','LD_PRELOAD','JULIA_LOAD_PATH'])
def test_command_environment_cannot_redirect_imports(bundle,key):
    entry=bundle['index']['constructors'][0];path=bundle['root']/entry['attempts'][0]['receipt']['path']
    receipt=json.loads(path.read_text());receipt['command']=['/usr/bin/env',key+'=/outside',*receipt['command']]
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert 'environment override' in report['constructors'][0]['attempts'][0]['failure']['message']


def test_effective_tmpdir_and_guard_cache_paths_are_checked(bundle):
    entry=bundle['index']['constructors'][0];path=bundle['root']/entry['attempts'][0]['receipt']['path']
    receipt=json.loads(path.read_text());receipt['environment']['TMPDIR']='/outside'
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert 'outside declared execution_root' in report['constructors'][0]['attempts'][0]['failure']['message']


@pytest.mark.parametrize('field',['settings','final_accounting'])
def test_malformed_nested_receipt_writes_refusal_audit(bundle,field):
    entry=bundle['index']['constructors'][0];path=bundle['root']/entry['attempts'][0]['receipt']['path']
    receipt=json.loads(path.read_text());receipt[field]=None
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    assert not assemble(bundle)['candidate_panel_complete']
    assert (bundle['output']/'PANEL.json').is_file()


def test_nonobject_attempt_is_caught_inside_the_audit_boundary(bundle):
    bundle['index']['constructors'][0]['attempts']=[None]
    bundle['whole_control']['auto']=False
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert 'Attempt object required' in report['constructors'][0]['attempts'][0]['failure']['message']
    assert (bundle['output']/'PANEL.json').is_file()


def test_nonobject_selected_input_hashes_write_a_refusal_audit(bundle):
    root=bundle['root'];index=bundle['index'];path=root/index['input_freeze']['path']
    frozen=json.loads(path.read_text());frozen['successes'][index['case']]=None
    index['input_freeze']=bundle['put'](path,frozen)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert any('hash object' in f['message'] for f in report['failures'])
    assert (bundle['output']/'PANEL.json').is_file()


def test_toml_unit_objects_are_checked_before_attribute_access(bundle):
    root=bundle['root'];index=bundle['index'];model=json.loads((root/index['model_json']['path']).read_text())
    raw=('case='+json.dumps(model['case'])+'\nn=14\nelectrons=4\npartition="conserving"\n'+
         'denominator="40"\nsource_json_sha256='+json.dumps(index['model_json']['sha256'])+
         '\nbaseline_order=[1,2]\nunits=[0,0]\n').encode()
    with pytest.raises(ValueError,match='TOML unit object'):
        panel.authenticate_toml(raw,model,index['model_json']['sha256'],panel.k.blocks(model))


@pytest.mark.parametrize('tamper',['missing_score','missing_components','wall','evaluations','moves','window_status'])
def test_selected_constructor_metadata_must_match_declared_semantics(bundle,tamper):
    entry=next(e for e in bundle['index']['constructors'] if e['arm']=='window_sparse')
    path=bundle['root']/entry['attempts'][0]['record']['path'];record=json.loads(path.read_text())
    if tamper=='missing_score': record['constructor']['score']=None;record['constructor']['score_components']=None
    elif tamper=='missing_components': record['constructor']['score_components']=None
    elif tamper=='wall': record['wall_seconds_before_output']=2.
    elif tamper=='evaluations': record['constructor']['evaluations']=True
    elif tamper=='moves': record['constructor']['accepted_moves']=-1
    else: record['constructor']['stats']['window_status']=[]
    entry['attempts'][0]['record']=bundle['put'](path,record)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert (bundle['output']/'PANEL.json').is_file()


def test_allocation_rounding_tolerance_does_not_extend_actual_totals():
    attempts=[dict(CPU_ceiling=64.,wall_ceiling=320.,CPU_seconds=.1,wall_seconds=.1),
              dict(CPU_ceiling=63.9,wall_ceiling=319.9,CPU_seconds=.2,wall_seconds=.2),
              dict(CPU_ceiling=63.7,wall_ceiling=319.7,CPU_seconds=0.,wall_seconds=0.)]
    assert panel.allocations_preserve_budget(attempts)
    attempts[-1]['CPU_ceiling']=63.70000001
    assert not panel.allocations_preserve_budget(attempts)


def test_preflight_tensor_identity_cannot_be_replaced_under_new_hashes(bundle):
    root=bundle['root'];index=bundle['index'];path=root/index['model_json']['path']
    model=json.loads(path.read_text());old=index['model_json']['sha256']
    model['tensors_sha256']='f'*64
    index['model_json']=bundle['put'](path,model)
    toml_path=root/index['model_toml']['path']
    toml_path.write_text(toml_path.read_text().replace(old,index['model_json']['sha256']))
    index['model_toml']=bundle['ref'](toml_path)
    freeze_path=root/index['input_freeze']['path'];frozen=json.loads(freeze_path.read_text())
    frozen['successes'][index['case']].update({'model.json':index['model_json']['sha256'],
                                             'model.toml':index['model_toml']['sha256']})
    index['input_freeze']=bundle['put'](freeze_path,frozen)
    protocol=json.loads((root/index['source_freeze']['path']).read_text())
    with pytest.raises(ValueError,match='tensor/orbital identities'):
        panel.input_binding(root,index,protocol,bundle['registry'])


def test_preflight_intrablock_order_is_part_of_the_supplied_target(bundle):
    root=bundle['root'];index=bundle['index'];parent_path=bundle['parent_path']
    frozen=json.loads(parent_path.read_text());chosen=min(frozen['successes'])
    parent_model_path=parent_path.parent/'inputs'/chosen/'model.json'
    parent=json.loads(parent_model_path.read_text())
    # One original block has literal10Z,2Z word order; the altered block swaps
    # only that order. Exact coefficients, unit membership and counts coincide.
    parent['blocks']=[[parent['coefficients'][1],parent['coefficients'][0]]]
    parent_ref=bundle['put'](parent_model_path,parent)
    def one_block_toml(path,model,json_hash,numerators):
        path.write_text('\n'.join(['case='+json.dumps(model['case']),'n=14','electrons=4',
            'partition="conserving"','denominator='+json.dumps(model['denominator']),
            'source_json_sha256='+json.dumps(json_hash),'baseline_order=[1]','[[units]]',
            'words=["2Z","10Z"]','numerators='+json.dumps(numerators)])+'\n')
        return bundle['ref'](path)
    parent_toml=one_block_toml(parent_model_path.with_suffix('.toml'),parent,parent_ref['sha256'],['3','5'])
    frozen['successes'][chosen].update({'model.json':parent_ref['sha256'],'model.toml':parent_toml['sha256']})
    index['parent_input_freeze']=bundle['put'](parent_path,frozen)
    altered_path=root/index['model_json']['path'];altered=json.loads(altered_path.read_text())
    altered['blocks']=[[altered['coefficients'][0],altered['coefficients'][1]]]
    altered['preflight_parent_model_sha256']=parent_ref['sha256']
    index['model_json']=bundle['put'](altered_path,altered)
    index['model_toml']=one_block_toml(altered_path.with_suffix('.toml'),altered,index['model_json']['sha256'],['21','35'])
    input_path=root/index['input_freeze']['path'];preflight=json.loads(input_path.read_text())
    preflight['parent_inventory_sha256']=index['parent_input_freeze']['sha256']
    preflight['successes'][index['case']].update({'model.json':index['model_json']['sha256'],
                                               'model.toml':index['model_toml']['sha256']})
    index['input_freeze']=bundle['put'](input_path,preflight)
    protocol=json.loads((root/index['source_freeze']['path']).read_text())
    with pytest.raises(ValueError,match='within-block word ordering'):
        panel.input_binding(root,index,protocol,bundle['registry'])


def test_guard_cache_environment_is_preserved_without_whitelisting_all_inherited_keys(bundle):
    for entry in bundle['index']['constructors']:
        path=bundle['root']/entry['attempts'][0]['receipt']['path'];receipt=json.loads(path.read_text())
        receipt['environment'].update(NUMBA_CACHE_DIR=str(bundle['root']/'logs'/'numba-cache'),
            MPLCONFIGDIR=str(bundle['root']/'logs'/'matplotlib-cache'),TMPDIR=str(bundle['root']/'tmp'))
        entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    assert assemble(bundle)['candidate_panel_complete']


@pytest.mark.parametrize('boundary',['CPU','wall'])
def test_actual_call_ceiling_has_no_allocation_rounding_extension(bundle,boundary):
    entry=bundle['index']['constructors'][0];path=bundle['root']/entry['attempts'][0]['receipt']['path']
    receipt=json.loads(path.read_text())
    if boundary=='CPU':
        receipt['budget']['cpu_seconds']=64.;receipt['final_accounting']['CPU_seconds']=64.+5e-10
    else:
        receipt['budget']['wall_seconds']=320.;receipt['settings']['RuntimeMaxSec']='320'
        receipt['elapsed']=receipt['final_accounting']['elapsed_seconds']=320.+5e-10
    entry['attempts'][0]['receipt']=bundle['put'](path,receipt)
    report=assemble(bundle)
    assert not report['candidate_panel_complete']
    assert report['constructors'][0]['within_total_construction_ceiling'] is False


def test_generated_preflight_preserves_parent_payload_identity_for_panel(bundle):
    import preflight
    root=bundle['root'];index=copy.deepcopy(bundle['index'])
    # Use the real altered-input producer, rather than the fixture's separately
    # authored altered model. Opaque parent payloads remain synthetic metadata.
    (root/'INPUT_PROTOCOL_FROZEN.json').write_bytes((root/index['source_freeze']['path']).read_bytes())
    output=root/'generated-preflight'
    preflight.alter(root,output)
    inventory=json.loads((output/'INPUTS_FROZEN.json').read_text())
    case=next(iter(inventory['successes']))
    model_path=output/'inputs'/case/'model.json'
    index.update(input_freeze=bundle['ref'](output/'INPUTS_FROZEN.json'),
                 model_json=bundle['ref'](model_path),model_toml=bundle['ref'](model_path.with_suffix('.toml')))
    protocol=json.loads((root/index['source_freeze']['path']).read_text())
    _,_,model,_,_=panel.input_binding(root,index,protocol,bundle['registry'])
    parent=bundle['original'][inventory['selected_case']][0]
    assert model['tensors_sha256']==parent['tensors_sha256']
    assert model['orbitals_sha256']==parent['orbitals_sha256']
