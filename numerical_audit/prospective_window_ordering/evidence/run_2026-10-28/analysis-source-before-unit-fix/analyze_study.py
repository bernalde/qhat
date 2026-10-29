"""Authenticated descriptive analysis of the fixed prospective window panel.

Usage: analyze_study.py ROOT NEW_SUMMARY_JSON [--mode main|preflight]
Reads ROOT/ANALYSIS_INDEX.<mode>.json; every referenced artifact is hash bound.
No optimizer, chemistry, propagation, package installation or subprocess runs.
This separate evidence tool does not alter the frozen scientific source map.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import statistics
import tomllib

SCHEMA='prospective-window-analysis-index-v1'
MAIN_CASES=('LiH_r1p487_631g_14m4e','LiH_r2p113_631g_14m4e',
            'BeH2_r1p369_631g_14m4e','BeH2_r1p907_631g_14m4e')
ARMS=('canonical','magnitude','current','first','pauli','cubic','window_sparse','window_magnitude')
FORMULAS={arm:(1,2) if arm in ('canonical','current','first') else (2,) for arm in ARMS}
TIMES=(.4,1.6)
CAPS=(8,24)
WINDOWS=('window_sparse','window_magnitude')
BASELINES=('cubic','current','magnitude','canonical')
METRICS=('projective_error','infidelity','norm_drift','number_leakage','spin_leakage')
EPS=1e-12
STATE_TOL=1e-10
HEX=re.compile(r'^[0-9a-f]{64}$')
MAX_FILE=32*1024**2
MAX_HASH_BYTES=512*1024**2
MAX_ARTIFACTS=8192
MAX_OUTPUT=32*1024**2
PROFILE=dict(CPUQuota='25%',CPUWeight='5',IOWeight='10',Nice='19',MemoryHigh='448M',
             MemoryMax='512M',MemorySwapMax='0',MemoryAccounting='yes',CPUAccounting='yes')


def require(condition,message):
    if not condition:
        raise ValueError(message)


def integer(value,name,low=0,high=10**9):
    require(type(value) is int and low<=value<=high,name+': integer domain')
    return value


def finite(value,name,nonnegative=True):
    require(type(value) in (int,float),name+': numeric scalar required')
    try: result=float(value)
    except OverflowError as error: raise ValueError(name+': overflow') from error
    require(math.isfinite(result) and (not nonnegative or result>=0),name+': finite numeric domain')
    return result


def close(a,b,tolerance=1e-9):
    return abs(a-b)<=tolerance+1e-12*max(abs(a),abs(b))


def issue(error):
    return dict(type=type(error).__name__,message=str(error)[:1000])


def decode_json(raw):
    def pairs(items):
        out={}
        for key,value in items:
            require(key not in out,'Duplicate JSON key: '+key);out[key]=value
        return out
    def floating(value):
        result=float(value);require(math.isfinite(result),'Nonfinite JSON float');return result
    def constant(value): raise ValueError('Nonfinite JSON constant: '+value)
    result=json.loads(raw,object_pairs_hook=pairs,parse_float=floating,parse_constant=constant)
    require(isinstance(result,dict),'Top-level JSON object required')
    return result


class Artifacts:
    def __init__(self,root):
        self.root=Path(root).resolve();self.hashes={};self.bytes_hashed=0

    def path(self,spelling,base=None):
        require(isinstance(spelling,str) and 0<len(spelling)<=2048,'Bounded artifact path required')
        rel=PurePosixPath(spelling)
        require(not rel.is_absolute() and '..' not in rel.parts and '\\' not in spelling,'Artifact path escape')
        target=(self.root if base is None else Path(base)).joinpath(*rel.parts).resolve()
        require(target.is_relative_to(self.root),'Artifact symlink escape')
        return target

    def relative(self,path): return str(Path(path).resolve().relative_to(self.root))

    def descriptor(self,path,sha256):
        return dict(path=self.relative(path),sha256=sha256)

    def read_path(self,path,expected=None,maximum=MAX_FILE,keep=True):
        path=Path(path).resolve();require(path.is_relative_to(self.root),'Artifact outside campaign root')
        require(path.is_file() and path.stat().st_size<=maximum,'Missing/oversize artifact: '+self.relative(path))
        name=self.relative(path)
        require(expected is None or isinstance(expected,str) and HEX.fullmatch(expected),'SHA-256 spelling')
        if name in self.hashes:
            require(expected is None or self.hashes[name]['sha256']==expected,'Conflicting artifact hash: '+name)
        else:
            require(len(self.hashes)<MAX_ARTIFACTS,'Artifact-count ceiling')
            h=hashlib.sha256();count=0
            with path.open('rb') as stream:
                while chunk:=stream.read(1024**2):
                    count+=len(chunk);self.bytes_hashed+=len(chunk)
                    require(count<=maximum and self.bytes_hashed<=MAX_HASH_BYTES,'Hashing byte ceiling')
                    h.update(chunk)
            observed=h.hexdigest()
            require(expected is None or observed==expected,'Artifact hash differs: '+name)
            self.hashes[name]=dict(sha256=observed,bytes=count)
        if not keep: return path,None
        with path.open('rb') as stream: raw=stream.read(maximum+1)
        require(len(raw)<=maximum and hashlib.sha256(raw).hexdigest()==self.hashes[name]['sha256'],
                'Artifact changed after authentication: '+name)
        return path,raw

    def read(self,descriptor,kind='json',maximum=MAX_FILE):
        require(isinstance(descriptor,dict) and {'path','sha256'}<=set(descriptor),'Hash-bound descriptor required')
        path,raw=self.read_path(self.path(descriptor['path']),descriptor['sha256'],maximum,keep=kind!='binary')
        if 'bytes' in descriptor: require(descriptor['bytes']==self.hashes[self.relative(path)]['bytes'],'Artifact byte-size binding')
        if kind=='binary': return path,None
        value=decode_json(raw) if kind=='json' else tomllib.loads(raw.decode('utf-8'))
        require(isinstance(value,dict),'Artifact table required')
        return path,value

    def child(self,parent,spelling,sha256,kind='binary',maximum=MAX_FILE):
        return self.read(self.descriptor(self.path(spelling,Path(parent).parent),sha256),kind,maximum)


def row_key(row,repetitions):
    arm=row['arm'];require(arm in ARMS,'Unknown arm')
    rep=integer(row['repetition'],'repetition',1,repetitions)
    formula=integer(row['formula'],'formula',1,2);require(formula in FORMULAS[arm],'Frozen formula policy differs')
    t=finite(row['time'],'time');cap=integer(row['budget_multiplier'],'budget multiplier',8,24)
    require(t in TIMES and cap in CAPS,'Frozen physical condition differs')
    return arm,rep,formula,t,cap


def expected_keys(repetitions):
    return [(arm,rep,formula,t,cap) for rep in range(1,repetitions+1) for arm in ARMS
            for formula in FORMULAS[arm] for t in TIMES for cap in CAPS]


def keyed(rows,repetitions,maximum=88):
    require(isinstance(rows,list) and len(rows)<=maximum,'Logical row ceiling')
    out={}
    for row in rows:
        require(isinstance(row,dict),'Logical row table required')
        key=row_key(row,repetitions);require(key not in out,'Duplicate logical row identity');out[key]=row
    return out


def call_key(row,cases,repetitions):
    require(row['case'] in cases and row['arm'] in ARMS,'Unknown constructor identity')
    return row['case'],row['arm'],integer(row['repetition'],'constructor repetition',1,repetitions)


def coeff(value):
    require(type(value) in (str,int),'Integer coefficient spelling')
    if isinstance(value,str): require(len(value)<=1280 and re.fullmatch(r'[+-]?[0-9]+',value),'Coefficient spelling bound')
    result=int(value);require(result.bit_length()<=4096,'Coefficient bit bound');return result


def words(raw):
    require(isinstance(raw,list) and len(raw)<=8192,'Word list bound')
    out={}
    for word,value in raw:
        require(isinstance(word,list) and len(word)<=14,'Word support bound')
        support=[]
        for q,letter in word:
            require(type(q) is int and 0<=q<14 and letter in ('X','Y','Z'),'Pauli site/letter')
            support.append((q,letter))
        key=tuple(support)
        require(key==tuple(sorted(key)) and len({q for q,_ in key})==len(key) and key not in out,'Canonical unique word')
        numerator=coeff(value);require(numerator!=0,'Zero supplied coefficient');out[key]=numerator
    return out


def model_definition(model,raw_toml,json_hash):
    require(type(model.get('n')) is int and model['n']==14 and
            type(model.get('electrons')) is int and model['electrons']==4,'Frozen 14-mode/4-electron target')
    denominator=coeff(model['denominator']);require(denominator>0,'Positive common denominator')
    total=words(model['coefficients']);blocks=[words(block) for block in model['blocks']]
    require(0<len(total)<=8192 and 0<len(blocks)<=1024 and all(blocks) and
            sum(map(len,blocks))==len(total) and {w:v for block in blocks for w,v in block.items()}==total,
            'Original block partition differs')
    for field in ('case','n','electrons'):
        require(raw_toml.get(field)==model[field],'JSON/TOML identity differs: '+field)
    require(raw_toml.get('source_json_sha256')==json_hash and raw_toml.get('partition')=='conserving' and
            coeff(raw_toml['denominator'])==denominator and
            raw_toml.get('baseline_order')==list(range(1,len(blocks)+1)),'Original JSON/TOML target binding')
    require(len(raw_toml['units'])==len(blocks),'Original TOML block cardinality')
    for compiled,original in zip(raw_toml['units'],blocks):
        ordered=sorted(original.items())
        require(compiled['words']==[','.join(str(q)+letter for q,letter in w) for w,_ in ordered] and
                [coeff(v) for v in compiled['numerators']]==[v for _,v in ordered],
                'Original numeric-word/TOML coefficient mapping differs')
    return dict(model=model,total=total,blocks=blocks,U=len(total),lengths=[len(block) for block in blocks])


def source_and_inputs(reader,index,mode,summary):
    protocol_path,protocol=reader.read(index['input_protocol'])
    require(protocol_path.name=='INPUT_PROTOCOL_FROZEN.json','Explicit generation protocol artifact')
    require(isinstance(protocol.get('source_commit'),str) and re.fullmatch(r'[0-9a-f]{40}',protocol['source_commit']),
            'Full frozen scientific source commit required')
    if 'expected_source_commit' in index:
        require(index['expected_source_commit']==protocol['source_commit'],'Declared source commit differs')
    source_map=protocol['sources'];require(isinstance(source_map,dict) and 0<len(source_map)<=4096,'Frozen source map bound')
    for name,sha256 in sorted(source_map.items()):
        reader.read_path(reader.path(name,protocol_path.parent/'source_snapshot'),sha256,keep=False)
    registry=protocol['registry'];specs=registry['cases']
    require(len(specs)==4 and {spec['case'] for spec in specs}==set(MAIN_CASES) and
            all(spec['n']==14 and spec['electrons']==4 for spec in specs),'Original four-case registry')
    registry_name='numerical_audit/prospective_window_ordering/REGISTRY.json'
    _,snapshot_registry=reader.read(reader.descriptor(reader.path(registry_name,protocol_path.parent/'source_snapshot'),source_map[registry_name]))
    require(snapshot_registry==registry,'Source/registry snapshot differs')
    config=protocol['configuration'];physical=config['physical']
    require(config['training_time']==1 and config['training_cap_times_original_words']==30 and
            config['fresh_repetitions']==2 and config['arms']==list(ARMS) and
            physical['times']==list(TIMES) and physical['cap_times_original_words']==list(CAPS) and
            physical['formulas']=={arm:list(FORMULAS[arm]) for arm in ARMS} and
            physical['win_tie_ambiguity']==EPS,'Frozen comparison configuration differs')
    summary['scientific_source']=dict(commit=protocol['source_commit'],input_protocol=index['input_protocol'],
        sources=source_map,registry=registry,configuration=config,
        PySCF_protocol_memory_fields={name:protocol[name] for name in
            ('internal_PySCF_max_memory_MB','internal_PySCF_memory_MiB') if name in protocol})
    models={};input_failures={};input_freeze=None
    if index.get('input_freeze') is not None:
        input_path,input_freeze=reader.read(index['input_freeze'])
        require(input_path.parent==protocol_path.parent and input_path.name=='INPUTS_FROZEN.json' and
                input_freeze.get('protocol_sha256')==index['input_protocol']['sha256'] and
                input_freeze.get('expected_cases')==4,'Original input freeze linkage')
        require(set(input_freeze.get('successes',{}))|set(input_freeze.get('refusals',{}))=={spec['case'] for spec in specs},
                'Input freeze must retain every declared case')
        require(not (set(input_freeze.get('successes',{}))&set(input_freeze.get('refusals',{}))),'Input success/refusal overlap')
        for spec in specs:
            case=spec['case'];hashes=input_freeze.get('successes',{}).get(case)
            if hashes is None:
                refusal=input_freeze.get('refusals',{}).get(case,{})
                if isinstance(refusal,dict) and {'path','sha256'}<=set(refusal):
                    _,record=reader.child(input_path,refusal['path'],refusal['sha256'],'json')
                    input_failures[case]=dict(descriptor=refusal,record=record)
                else: input_failures[case]=refusal
                continue
            try:
                base=input_path.parent/'inputs'/case
                jp,j=reader.read(reader.descriptor(base/'model.json',hashes['model.json']),maximum=8*1024**2)
                tp,t=reader.read(reader.descriptor(base/'model.toml',hashes['model.toml']),'toml',8*1024**2)
                require(all(j.get(key)==value for key,value in spec.items()),'Input/registry case identity')
                for field in ('basis','unit','charge','spin'): require(j.get(field)==registry[field],'Input/global registry identity')
                for filename,key in (('tensors.npz','tensors_sha256'),('orbitals.npz','orbitals_sha256')):
                    require(hashes[filename]==j[key],'SCF tensor/orbital hash link')
                    reader.read(reader.descriptor(base/filename,hashes[filename]),'binary')
                definition=model_definition(j,t,hashes['model.json'])
                definition.update(json_hash=hashes['model.json'],toml_hash=hashes['model.toml'],json_path=reader.relative(jp),toml_path=reader.relative(tp))
                models[case]=definition
            except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
                input_failures[case]=issue(error)
    else:
        input_failures={spec['case']:'No final INPUTS_FROZEN artifact supplied' for spec in specs}
    original_complete=(input_freeze is not None and input_freeze.get('status')=='complete' and len(models)==4 and not input_failures)
    summary['inputs']=dict(original_complete=original_complete,refusals_or_invalid=input_failures,
        frozen_status=None if input_freeze is None else input_freeze.get('status'))
    if mode=='main': return protocol,registry,models,[spec['case'] for spec in specs],original_complete
    altered_ref=index.get('preflight_input_freeze')
    if altered_ref is None:
        summary['inputs']['preflight_failure']='No altered-input freeze supplied'
        return protocol,registry,{},['unperformed_preflight'],False
    ap,altered=reader.read(altered_ref)
    require(altered.get('expected_cases')==1 and altered.get('preflight_only') is True and
            altered.get('parent_inventory_sha256')==index['input_freeze']['sha256'] and
            altered.get('protocol_sha256')==index['input_protocol']['sha256'] and
            Path(altered.get('source_protocol_path','')).resolve()==protocol_path and
            not altered.get('refusals'),'Altered-input inventory binding')
    names=list(altered.get('successes',{}));require(len(names)==1,'One altered preflight target required')
    case=names[0]
    if not original_complete:
        return protocol,registry,{},[case],False
    chosen=min(models,key=lambda name:(-models[name]['U'],name));parent=models[chosen]
    hashes=altered['successes'][case];base=ap.parent/'inputs'/case
    jp,j=reader.read(reader.descriptor(base/'model.json',hashes['model.json']),maximum=8*1024**2)
    tp,t=reader.read(reader.descriptor(base/'model.toml',hashes['model.toml']),'toml',8*1024**2)
    definition=model_definition(j,t,hashes['model.json'])
    require(all(j.get(field)==parent['model'].get(field) for field in
        ('n','electrons','family','geometry','basis','unit','charge','spin','active_spatial','frozen_spatial')),
        'Altered-input geometry/basis/active-space substitution')
    require(case==chosen+'_preflight_7over5' and j.get('preflight_parent_case')==chosen and
            j.get('preflight_parent_model_sha256')==parent['json_hash'] and j.get('preflight_scale')=='7/5' and
            coeff(j['denominator'])==5*coeff(parent['model']['denominator']) and
            definition['total']=={w:7*v for w,v in parent['total'].items()} and
            definition['blocks']==[{w:7*v for w,v in block.items()} for block in parent['blocks']],
            'Exact largest-word altered-input target differs')
    definition.update(json_hash=hashes['model.json'],toml_hash=hashes['model.toml'],json_path=reader.relative(jp),toml_path=reader.relative(tp))
    return protocol,registry,{case:definition},[case],altered.get('status')=='complete'


def receipt(reader,descriptor):
    _,data=reader.read(descriptor)
    final=data.get('final_accounting',{})
    require(all(data.get('settings',{}).get(name)==value for name,value in PROFILE.items()),'Mandatory guard resource profile differs')
    cpu=None if final.get('CPU_seconds') is None else finite(final['CPU_seconds'],'guard CPU')
    wall=None if final.get('elapsed_seconds') is None else finite(final['elapsed_seconds'],'guard wall')
    complete=(final.get('CPU_accounting_available') is True and final.get('process_CPU_timing_available') is True and
              data.get('guard_error') is None and not data.get('cleanup_errors') and cpu is not None and wall is not None)
    return dict(descriptor=descriptor,unit=data.get('unit'),CPU_seconds=cpu,wall_seconds=wall,
        accounting_complete=complete,exit_code=data.get('exit_code'),termination=data.get('termination'),
        guard_error=data.get('guard_error'),cleanup_errors=data.get('cleanup_errors',[]),
        frozen_seconds=data.get('frozen_seconds'),budget=data.get('budget'),command=data.get('command'))


def global_inventory(reader,index,protocol,cases,models,mode,summary):
    repetitions=2 if mode=='main' else 1
    expected={(case,arm,rep) for case in cases for arm in ARMS for rep in range(1,repetitions+1)}
    plan_map={};global_map={}
    if index.get('plan') is not None:
        _,plan=reader.read(index['plan'])
        require(plan.get('schema')=='prospective-window-constructor-plan-v1' and plan.get('mode')==mode and
                plan.get('expected_calls')==len(expected) and plan.get('source_snapshot')==protocol['sources'] and
                plan.get('source_freeze')==index['input_protocol'] and plan.get('physical_error_read') is False and
                plan.get('input_freeze')==index.get('preflight_input_freeze' if mode=='preflight' else 'input_freeze'),
                'Planned call/source identity differs')
        entries=plan['entries'];require(isinstance(entries,list) and len(entries)==len(expected),'Exact planned call roster')
        for entry in entries:
            key=call_key(entry,cases,repetitions);require(key not in plan_map,'Duplicate planned call');plan_map[key]=entry
        require(set(plan_map)==expected,'Missing planned call')
        summary['plan']=dict(status='authenticated',descriptor=index['plan'],expected_calls=len(expected))
    else: summary['plan']=dict(status='unavailable',expected_calls=len(expected))
    if index.get('all_orders_freeze') is not None:
        _,whole=reader.read(index['all_orders_freeze'])
        require(plan_map and whole.get('schema')=='prospective-window-all-orders-freeze-v1' and
                whole.get('mode')==mode and whole.get('inventory_complete') is True and
                whole.get('source_freeze')==index['input_protocol'] and
                whole.get('plan_sha256')==index['plan']['sha256'] and whole.get('expected_calls')==len(expected) and
                whole.get('physical_error_read') is False and whole.get('no_selected_repetition') is True,
                'Whole-order freeze/plan/source linkage differs')
        artifacts=whole['artifacts'];require(isinstance(artifacts,list) and len(artifacts)==len(expected),'Whole-order call roster bound')
        units=set();hashes=set();missing=[]
        for row in artifacts:
            key=call_key(row,cases,repetitions);require(key not in global_map,'Duplicate whole-order identity')
            rec=receipt(reader,row['receipt']);require(rec['unit'] not in units and row['receipt']['sha256'] not in hashes,'Guard receipt reused globally')
            units.add(rec['unit']);hashes.add(row['receipt']['sha256'])
            require(rec['command']==plan_map[key]['command'] and rec['exit_code']==row.get('exit_code'), 'Executed/planned command or exit differs')
            item=dict(artifact=row,receipt=rec)
            if row.get('record') is not None:
                _,record=reader.read(row['record'])
                require(record.get('schema')=='prospective-window-constructor-v1' and record.get('artifact_complete') is True and
                        (record.get('case'),record.get('arm'),record.get('repetition'))==key and
                        record.get('source_freeze_sha256')==index['input_protocol']['sha256'] and
                        key[0] in models and record.get('model_sha256')==models[key[0]]['json_hash'] and
                        record.get('physical_error_read') is False,'Whole-order constructor target/source identity')
                item['record']=record
            else: missing.append(dict(case=key[0],arm=key[1],repetition=key[2],receipt=row['receipt']))
            global_map[key]=item
        require(set(global_map)==expected and type(whole.get('all_calls_have_artifacts')) is bool and
                whole['all_calls_have_artifacts']==(not missing),'Whole-order missing/presence identity')
        checkpoints=whole.get('retained_temporary_and_checkpoint_artifacts',[])
        require(isinstance(checkpoints,list) and len(checkpoints)<=512,'Checkpoint roster bound')
        for descriptor in checkpoints: reader.read(descriptor,'binary',256*1024**2)
        summary['whole_orders']=dict(status='authenticated',descriptor=index['all_orders_freeze'],expected_calls=len(expected),
            artifacts_present=len(expected)-len(missing),censored=missing,all_calls_have_artifacts=not missing,
            scope='Complete attempted-call roster, not proof of completed optimizer stages or physical comparison')
    else: summary['whole_orders']=dict(status='unavailable',expected_calls=len(expected),censored=[])
    return expected,plan_map,global_map


def objective_rows(check,candidates,repetitions):
    logical=check.get('logical_rows',[]);require(isinstance(logical,list) and len(logical)<=88,'Objective logical row bound')
    checks=check.get('checks',[]);require(isinstance(checks,list) and len(checks)<=256,'Objective unique task bound')
    output={};seen_indices=set()
    for row in logical:
        key=row_key(row,repetitions);require(key in candidates and key not in output,'Objective identity differs')
        candidate=candidates[key]
        require(row.get('name')==candidate['name'] and
                ('partition' not in row or row['partition']==candidate['partition']),'Objective candidate identity differs')
        position=integer(row['candidate_index'],'objective candidate index',1,88)
        require(position<=len(candidates) and list(candidates)[position-1]==key and
                position not in seen_indices,'Objective candidate position/identity differs');seen_indices.add(position)
        if not candidate.get('constructor_scored',False):
            require('constructor_score' not in candidate and row.get('status')=='unscored_censored' and
                    'check_index' not in row,'Unscored constructor/objective claim differs')
            output[key]=dict(status='unscored_no_numeric_claim',accepted=True,driver_status=row['status'])
            continue
        if 'check_index' not in row:
            require(row.get('status') not in ('verified','unscored_censored') and isinstance(row.get('status'),str),
                    'Scored objective must retain its refusal status or numeric check identity')
            output[key]=dict(status=row['status'],accepted=False,reason=row.get('reason'),
                claimed_score=candidate['constructor_score'],check_index=None)
            continue
        index=integer(row['check_index'],'objective check index',1,len(checks));record=checks[index-1]
        require(record.get('check_index')==index and record.get('partition')==candidate['partition'] and
                record.get('order')==candidate['order'] and record.get('score_kind')==candidate['score_kind'] and
                finite(record['constructor_score'],'claimed objective')==finite(candidate['constructor_score'],'candidate objective'),
                'Objective check/order/score identity differs')
        require(row.get('status')==record.get('status'),'Objective logical/check status differs')
        accepted=(record.get('status')=='verified' and record.get('verified') is True and record.get('within_tolerance') is True)
        output[key]=dict(status=record.get('status'),accepted=accepted,claimed_score=record['constructor_score'],
            replayed_score=record.get('replayed_score'),absolute_difference=record.get('absolute_difference'),
            numerical_floor=record.get('numerical_floor'),check_index=index,reason=record.get('reason'))
    return output


def metric_values(row):
    values={name:finite(row[name],name) for name in METRICS}
    require(all(values[name]<=1+1e-10 for name in ('projective_error','infidelity','number_leakage','spin_leakage')) and
            values['norm_drift']<=1e-9,'Physical metric domain/norm check')
    return values


def affordability(model,candidate,physical):
    lengths=model['lengths'] if candidate['partition']=='conserving' else [1]*model['U']
    order=candidate['order']
    require(isinstance(order,list) and all(type(i) is int for i in order) and sorted(order)==list(range(1,len(lengths)+1)),
            'Complete original one-based order')
    ordered=[lengths[i-1] for i in order];budget=physical['budget_multiplier']*model['U'];formula=physical['formula']
    r=(1 if budget>=sum(lengths) else 0) if len(lengths)==1 else (
        budget//sum(lengths) if formula==1 else (budget-ordered[0])//(2*sum(lengths)-ordered[0]-ordered[-1]))
    rotations=sum(lengths) if len(lengths)==1 else (r*sum(lengths) if formula==1 else
        r*(2*sum(lengths)-ordered[0]-ordered[-1])+ordered[0])
    require(r>=1 and physical.get('rotation_budget')==budget and physical.get('steps')==r and
            physical.get('ideal_rotations')==rotations<=budget,'Original-word affordability/count differs')


def panel_costs(reader,panel,inventory,case,repetitions,global_map):
    expected={(arm,rep) for arm in ARMS for rep in range(1,repetitions+1)}
    entries={(entry['arm'],entry['repetition']):entry for entry in inventory['constructors']}
    retained={(entry['arm'],entry['repetition']):entry for entry in panel['constructors']}
    require(len(entries)==len(inventory['constructors']) and len(retained)==len(panel['constructors']) and
            set(entries)==expected and set(retained)==expected,'Panel constructor roster differs')
    costs={};all_receipts=set();all_units=set()
    for arm,rep in sorted(expected):
        key=(arm,rep);item=retained[key];ledger=entries[key];attempts=item.get('attempts',[])
        require(len(attempts)==len(ledger['attempts']) and len(attempts)<=32,'Complete attempt-ledger identity')
        receipts=[];complete=True;known_cpu=known_wall=0.;hashes=set();within=True
        for recorded,original in zip(attempts,ledger['attempts']):
            require(recorded.get('receipt')==original.get('receipt') and recorded.get('constructor')==original.get('record'),
                    'PANEL/constructor attempt descriptor differs')
            require(original['receipt']['sha256'] not in hashes,'Repeated receipt within a constructor ledger')
            hashes.add(original['receipt']['sha256'])
            measured=receipt(reader,original['receipt']);receipts.append(measured)
            require(original['receipt']['sha256'] not in all_receipts and measured['unit'] not in all_units,
                    'One guard receipt/unit charged to multiple constructor attempts')
            all_receipts.add(original['receipt']['sha256']);all_units.add(measured['unit'])
            budget=measured['budget'];cpu_cap=finite(budget['cpu_seconds'],'constructor attempt CPU ceiling')
            wall_cap=finite(budget['wall_seconds'],'constructor attempt wall ceiling')
            require(cpu_cap>0 and wall_cap>0 and cpu_cap<=64-known_cpu and wall_cap<=320-known_wall,
                    'Constructor retry allocation silently extends the whole-call ceiling')
            within &= (measured['CPU_seconds'] is not None and measured['wall_seconds'] is not None and
                       measured['CPU_seconds']<=cpu_cap and measured['wall_seconds']<=wall_cap)
            if measured['CPU_seconds'] is not None: known_cpu+=measured['CPU_seconds']
            if measured['wall_seconds'] is not None: known_wall+=measured['wall_seconds']
            complete &= measured['accounting_complete']
        if not attempts:
            complete=False
        if complete:
            require(item.get('cost_complete') is True and close(finite(item['guard_CPU_seconds'],'panel CPU'),known_cpu) and
                    close(finite(item['guard_wall_seconds'],'panel wall'),known_wall),'PANEL/whole guard costs differ')
            require(item.get('within_total_construction_ceiling')==bool(within and known_cpu<=64 and known_wall<=320),
                    'PANEL cumulative construction-ceiling claim differs from raw receipts')
        else: require(item.get('guard_CPU_seconds') is None and item.get('guard_wall_seconds') is None,
                      'Unknown guard cost must remain unknown')
        final=ledger['attempts'][-1] if ledger['attempts'] else None
        if (case,arm,rep) in global_map:
            frozen=global_map[(case,arm,rep)]['artifact']
            require(final is not None and final['receipt']==frozen['receipt'] and final.get('record')==frozen.get('record'),
                    'Case/global final call differs')
        selected=item.get('selected_record')
        if item.get('status')=='retained_order':
            require(selected is not None and selected==final.get('record') and
                    item.get('selected_attempt')==ledger.get('selected_attempt')==len(ledger['attempts']),
                    'Final selected constructor descriptor differs')
            _,record=reader.read(selected)
            require((record.get('case'),record.get('arm'),record.get('repetition'))==(case,arm,rep),'Selected constructor identity differs')
            require(item['constructor']['order']==record['constructor']['order'] and
                    item['constructor']['score']==record['constructor'].get('score'),'PANEL/selected constructor witness differs')
        costs[key]=dict(case=case,arm=arm,repetition=rep,status=item.get('status'),
            constructor=item.get('constructor'),reason=item.get('reason'),selected_record=selected,
            CPU_seconds=known_cpu if complete else None,wall_seconds=known_wall if complete else None,
            known_CPU_seconds=known_cpu,known_wall_seconds=known_wall,accounting_complete=complete,
            attempts=receipts,source='PANEL complete per-call attempt ledger, checked against raw guard receipts',
            within_limits=item.get('within_total_construction_ceiling'))
    return costs


def case_artifacts(reader,entry,index,case,model,repetitions,global_map,protocol,summary):
    result=dict(case=case,status='partial',graph_valid=False,failures=[],objective_status=None,independent_status=None,physical_status=None,
                primary_panel_complete=False,independent_panel_complete=False,objective_claims_complete=False)
    costs={};candidates={};physical={};independent={};objectives={};references=[]
    if entry is None:
        result['failures'].append(dict(message='No analysis-index entry for the declared case'))
        return result,costs,candidates,physical,independent,objectives,references
    try:
        inventory_path,inventory=reader.read(entry['inventory'])
        require(inventory_path.parent==reader.root and inventory.get('case')==case and inventory.get('mode')==index['mode'] and
                inventory.get('source_freeze')==index['input_protocol'] and
                inventory.get('schema')=='prospective-window-panel-inventory-v1' and
                inventory.get('input_freeze')==index.get('preflight_input_freeze' if index['mode']=='preflight' else 'input_freeze') and
                inventory.get('model_json')==dict(path=model['json_path'],sha256=model['json_hash']) and
                inventory.get('model_toml')==dict(path=model['toml_path'],sha256=model['toml_hash']) and
                inventory.get('all_orders_freeze')==index.get('all_orders_freeze'),'Case/global inventory/source linkage')
        panel_path,panel=reader.read(entry['panel'])
        require(panel.get('schema')=='prospective-window-panel-audit-v1' and
                panel.get('inventory_sha256')==entry['inventory']['sha256'] and panel.get('case')==case and
                panel.get('mode')==index['mode'] and panel.get('source_freeze')==index['input_protocol'] and
                panel.get('source_commit')==protocol['source_commit'] and
                panel.get('model_json')==inventory['model_json'] and panel.get('model_toml')==inventory['model_toml'],
                'PANEL/inventory/source linkage')
        costs=panel_costs(reader,panel,inventory,case,repetitions,global_map)
        result['panel_status']=panel.get('status');result['censor_reason']=entry.get('censor_reason')
        if entry.get('candidates') is None:
            require(panel.get('candidate_panel_complete') is not True,'Complete PANEL must retain its candidates')
            return result,costs,candidates,physical,independent,objectives,references
        cp,c=reader.read(entry['candidates'],'toml',8*1024**2)
        declared=panel['candidates']
        require(reader.path(declared['path'],panel_path.parent)==cp and declared['sha256']==entry['candidates']['sha256'] and
                panel.get('status')=='complete' and panel.get('candidate_panel_complete') is True and
                c.get('schema')=='prospective-window-candidates-v1' and c.get('status')=='complete' and
                c.get('case')==case and c.get('mode')==index['mode'] and c.get('input_sha256')==model['toml_hash'] and
                c.get('source_json_sha256')==model['json_hash'] and
                c.get('source_commit')==protocol['source_commit'] and c.get('n')==14 and c.get('hf')=='15' and
                c.get('expected_logical_rows')==44*repetitions and
                c.get('input_freeze_sha256')==inventory['input_freeze']['sha256'] and
                c.get('inventory_sha256')==entry['inventory']['sha256'] and
                c.get('source_freeze_sha256')==index['input_protocol']['sha256'] and
                index.get('all_orders_freeze') is not None and c.get('all_orders_freeze_sha256')==index['all_orders_freeze']['sha256'],
                'Frozen candidate/input/source/PANEL linkage')
        candidates=keyed(c['candidates'],repetitions)
        require(set(candidates)==set(expected_keys(repetitions)),'Complete candidate matrix required')
        for key,row in candidates.items():
            arm,rep=key[:2];cost=costs[(arm,rep)];ctor=cost['constructor']
            require(cost['status']=='retained_order' and cost['accounting_complete'] and cost['within_limits'] is True and
                    row.get('constructor_record_sha256')==cost['selected_record']['sha256'] and
                    row.get('order')==[i+1 for i in ctor['order']] and
                    row.get('partition')==('pauli' if arm=='pauli' else 'conserving') and
                    row.get('constructor_scored')==(ctor['score'] is not None), 'Candidate/constructor order/cost identity')
            require(close(finite(row['constructor_guard_CPU_seconds'],'candidate CPU'),cost['CPU_seconds']) and
                    close(finite(row['constructor_guard_wall_seconds'],'candidate wall'),cost['wall_seconds']),
                    'Candidate construction cost is not the full attempt ledger')
            if ctor['score'] is not None:
                finite(ctor['score'],'constructor score')
                require(row.get('constructor_score')==ctor['score'] and row.get('score_kind')==ctor['score_kind'],
                        'Candidate constructor score claim differs')
                require(ctor['score_kind'] in ('retained_first_order_response','diagonal_reference_cubic_response'),
                        'Unknown scored constructor objective')
        if entry.get('objective_check') is not None:
            _,objective=reader.read(entry['objective_check'],'toml')
            require(objective.get('schema')=='prospective-window-objective-check-v1' and
                    objective.get('input_sha256')==model['toml_hash'] and
                    objective.get('candidates_sha256')==entry['candidates']['sha256'] and
                    objective.get('case')==case and objective.get('n')==14 and objective.get('hf')=='15' and
                    objective.get('training_time')==1 and objective.get('training_budget_multiplier')==30 and
                    finite(objective['relative_tolerance'],'objective relative tolerance')<=1e-8 and
                    finite(objective['absolute_tolerance'],'objective absolute tolerance')<=1e-24,
                    'Objective source/target/training/tolerance linkage')
            require(set(objective.get('source_hashes',{}))=={'objective_check.jl','src/WindowStudyIO.jl',
                '../julia_cubic_feasibility/src/ResponseCore.jl','../scalable_ordering/training.py',
                '../../analysis/sparse_response_ordering.py'},'Complete objective replay source-hash inventory required')
            for name,sha256 in objective.get('source_hashes',{}).items():
                normalized=PurePosixPath(os.path.normpath('numerical_audit/prospective_window_ordering/'+name))
                require(protocol['sources'].get(str(normalized))==sha256,'Objective replay scientific source differs')
            objectives=objective_rows(objective,candidates,repetitions)
            result['objective_status']=objective.get('status')
            scored=any(row.get('constructor_scored') is True for row in candidates.values())
            result['objective_claims_complete']=(set(objectives)==set(candidates) and all(v['accepted'] for v in objectives.values()) and
                (objective.get('complete') is True if scored else objective.get('status')=='no_scored_orders_supplied'))
        if entry.get('physical') is None: return result,costs,candidates,physical,independent,objectives,references
        pp,p=reader.read(entry['physical'],'toml')
        require(p.get('case')==case and p.get('mode')==index['mode'] and p.get('input_sha256')==model['toml_hash'] and
                p.get('source_json_sha256')==model['json_hash'] and p.get('orders_sha256')==entry['candidates']['sha256'] and
                p.get('n')==14 and p.get('hf')=='15' and p.get('expected_logical_rows')==44*repetitions,
                'Physical manifest/frozen candidate target differs')
        physical=keyed(p['rows'],repetitions);result['physical_status']=p.get('status')
        result['primary_panel_complete']=(p.get('status')=='complete' and p.get('physical_panel_complete') is True and
                                           set(physical)==set(candidates))
        for key,row in physical.items():
            require(key in candidates and row.get('name')==candidates[key]['name'] and row.get('order')==candidates[key]['order'] and
                    row.get('partition')==candidates[key]['partition'],'Physical/frozen candidate order differs')
            affordability(model,candidates[key],row);metric_values(row)
            for prefix in ('actual','reference'):
                state_path,_=reader.child(pp,row[prefix+'_file'],row[prefix+'_sha256'],'binary',16*16384)
                require(reader.hashes[reader.relative(state_path)]['bytes']==16*16384,'Full physical state exact file length')
        if entry.get('independent_check') is None: return result,costs,candidates,physical,independent,objectives,references
        _,check=reader.read(entry['independent_check'],'toml')
        require(check.get('input_toml_sha256')==model['toml_hash'] and check.get('source_json_sha256')==model['json_hash'] and
                check.get('manifest_sha256')==entry['physical']['sha256'] and check.get('case')==case and
                check.get('n')==14 and check.get('hf')=='15' and
                check.get('checker_source_sha256')==protocol['sources']['numerical_audit/prospective_window_ordering/independent.jl'],
                'Independent checker/source/physical binding differs')
        require(finite(check['state_norm_tolerance'],'state check tolerance')<=STATE_TOL and
                finite(check['absolute_metric_tolerance'],'metric check tolerance')<=EPS,'Independent tolerances exceed frozen limits')
        independent=keyed(check.get('rows',[]),repetitions);references=check.get('reference_checks',[])
        require(isinstance(references,list) and len(references)<=88,'Independent reference-check bound')
        valid_refs=set()
        for reference in references:
            if reference.get('status')=='complete':
                require(finite(reference['state_norm_residual'],'reference residual')<=STATE_TOL and
                        finite(reference['phase_aligned_residual'],'phase-aligned reference residual')<=STATE_TOL,
                        'Independent reference residual exceeds declared tolerance')
                valid_refs.add((reference['time'],reference['reference_file'],reference['reference_sha256']))
        for key,row in independent.items():
            require(key in physical and row.get('name')==physical[key]['name'],'Independent physical row identity differs')
            if row.get('status')!='complete': continue
            primary=physical[key]
            require(row.get('saved_actual_sha256')==primary['actual_sha256'] and
                    (key[3],primary['reference_file'],primary['reference_sha256']) in valid_refs,
                    'Independent state/reference hash coverage differs')
            for metric,value in metric_values(primary).items():
                require(abs(finite(row['metrics'][metric],'independent metric')-value)<=EPS,'Independent metric differs')
            require(finite(primary['reference_norm_drift'],'primary reference norm drift')<=1e-9 and
                    abs(finite(row['metrics']['reference_norm_drift'],'independent reference norm drift')-
                        primary['reference_norm_drift'])<=EPS,'Independent reference norm drift differs')
            hard=key[3]==1.6 and key[4]==24
            require(row.get('independent_unfused_replay') is hard,'Independent propagation scope differs')
            if hard: require(finite(row['unfused_state_norm_residual'],'unfused state residual')<=STATE_TOL,'Unfused replay residual')
        result['independent_status']=check.get('status')
        result['independent_panel_complete']=(check.get('status')=='complete' and not check.get('failures') and
            check.get('checked_logical_rows')==len(candidates) and set(independent)==set(candidates))
        result['status']='authenticated';result['graph_valid']=True
    except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
        result['failures'].append(issue(error))
        # Retain quarantined raw rows, but never derive from an invalid graph.
    return result,costs,candidates,physical,independent,objectives,references


def partial_costs(reader,index,cases,repetitions,global_map,existing):
    records=index.get('partial_attempts',[])
    require(isinstance(records,list) and len(records)<=2048,'Partial attempt bound')
    grouped=defaultdict(list);seen=set()
    for row in records:
        require(row['receipt']['sha256'] not in seen,'Partial attempt receipt repeated')
        seen.add(row['receipt']['sha256'])
        key=call_key(row,cases,repetitions);measured=receipt(reader,row['receipt'])
        grouped[key].append(measured)
        if row.get('record') is not None: reader.read(row['record'])
    for key,attempts in grouped.items():
        if key in existing:
            require({a['descriptor']['sha256'] for a in attempts}<=
                    {a['descriptor']['sha256'] for a in existing[key]['attempts']},'Partial ledger omits a PANEL attempt linkage')
            continue
        complete=all(a['accounting_complete'] for a in attempts)
        existing[key]=dict(case=key[0],arm=key[1],repetition=key[2],status='partial_attempts_only',constructor=None,
            attempts=attempts,accounting_complete=False,provided_receipts_accounting_complete=complete,
            CPU_seconds=None,wall_seconds=None,
            known_CPU_seconds=sum(a['CPU_seconds'] for a in attempts if a['CPU_seconds'] is not None),
            known_wall_seconds=sum(a['wall_seconds'] for a in attempts if a['wall_seconds'] is not None),
            source='Hash-bound partial attempt ledger; total whole-call cost remains unknown without a complete PANEL ledger')
    for key,item in global_map.items():
        if key not in existing:
            a=item['receipt'];record=item.get('record');body=None if record is None else record.get('constructor')
            ctor=None if not isinstance(body,dict) else dict(status=body.get('status'),score=body.get('score'),
                score_kind=record.get('score_kind'),order=body.get('order'),reason=body.get('reason'),
                window_status=body.get('stats',{}).get('window_status'))
            existing[key]=dict(case=key[0],arm=key[1],repetition=key[2],status='whole_freeze_without_PANEL',
                constructor=ctor,attempts=[a],accounting_complete=False,CPU_seconds=None,wall_seconds=None,
                known_CPU_seconds=a['CPU_seconds'] if a['CPU_seconds'] is not None else 0.,
                known_wall_seconds=a['wall_seconds'] if a['wall_seconds'] is not None else 0.,
                source='Whole-freeze final receipt only; earlier attempts unknown without PANEL')


def row_record(case,key,case_result,cost,candidate,primary,independent,objective):
    arm,rep,formula,t,cap=key
    out=dict(case=case,arm=arm,repetition=rep,formula=formula,time=t,budget_multiplier=cap,
        path=arm+'_f'+str(formula),status='censored',physical_checked=False,objective_accepted=False,
        constructor_CPU_seconds=None if cost is None else cost['CPU_seconds'],
        constructor_wall_seconds=None if cost is None else cost['wall_seconds'],
        constructor_status=None if cost is None or cost.get('constructor') is None else cost['constructor'].get('status'),
        constructor_window_status=None if cost is None or cost.get('constructor') is None else cost['constructor'].get('window_status'),
        constructor_scored=None if cost is None or cost.get('constructor') is None else cost['constructor'].get('score') is not None,
        candidate=candidate,physical=primary,independent=independent,objective=objective,
        case_graph_status=case_result['status'],case_failures=case_result['failures'])
    if case_result.get('graph_valid') is True and primary is not None and independent is not None and independent.get('status')=='complete':
        out.update(status='physical_checked',physical_checked=True,metrics=metric_values(primary),
            physical_floor=primary['projective_error']<=EPS,objective_accepted=bool(objective and objective.get('accepted')))
        out['fully_validated']=out['objective_accepted']
    else:
        out['censor_reason']=case_result.get('censor_reason') or 'Missing or incomplete authenticated physical/independent row'
        out['fully_validated']=False
    return out


def classify(delta): return 'win' if delta < -EPS else 'loss' if delta > EPS else 'tie'


def paired_derivations(raw,cases,repetitions):
    lookup={(r['case'],r['arm'],r['repetition'],r['formula'],r['time'],r['budget_multiplier']):r for r in raw}
    pairs=[];medians=[];summaries=[]
    for window in WINDOWS:
        for baseline in BASELINES:
            observations=[]
            for case in cases:
                for rep in range(1,repetitions+1):
                    for t in TIMES:
                        for cap in CAPS:
                            w=lookup[(case,window,rep,2,t,cap)];b=lookup[(case,baseline,rep,2,t,cap)]
                            record=dict(window=window,baseline=baseline,case=case,repetition=rep,formula=2,time=t,budget_multiplier=cap,
                                window_row_status=w['status'],baseline_row_status=b['status'],classification='censored',
                                error_delta=None,infidelity_delta=None,CPU_delta=None,wall_delta=None,
                                baseline_over_window_error_ratio=None,ratio_suppression_reason='Unknown/incomplete physical metric',
                                fully_validated=False)
                            if w['physical_checked'] and b['physical_checked']:
                                we=w['metrics']['projective_error'];be=b['metrics']['projective_error'];delta=we-be
                                record.update(error_delta=delta,infidelity_delta=w['metrics']['infidelity']-b['metrics']['infidelity'],
                                    classification=classify(delta),window_error=we,baseline_error=be,
                                    window_floor=we<=EPS,baseline_floor=be<=EPS,
                                    fully_validated=w['fully_validated'] and b['fully_validated'])
                                if min(we,be)<=EPS: record['ratio_suppression_reason']='At or below absolute physical floor/ambiguity'
                                elif not record['fully_validated']: record['ratio_suppression_reason']='Objective check censored or unverified'
                                else:
                                    ratio=be/we
                                    if math.isfinite(ratio): record.update(baseline_over_window_error_ratio=ratio,ratio_suppression_reason=None)
                                    else: record['ratio_suppression_reason']='Ratio overflows'
                            for field,target in (('constructor_CPU_seconds','CPU_delta'),('constructor_wall_seconds','wall_delta')):
                                if w[field] is not None and b[field] is not None: record[target]=w[field]-b[field]
                            pairs.append(record);observations.append(record)
            complete_groups=[]
            for case in cases:
                for t in TIMES:
                    for cap in CAPS:
                        group=[p for p in observations if p['case']==case and p['time']==t and p['budget_multiplier']==cap]
                        available=[p for p in group if p['error_delta'] is not None]
                        validated=[p for p in group if p['fully_validated']]
                        complete=len(available)==repetitions
                        value=statistics.median(p['error_delta'] for p in available) if available else None
                        median=dict(window=window,baseline=baseline,case=case,formula=2,time=t,budget_multiplier=cap,
                            expected_repetitions=repetitions,available_repetitions=[p['repetition'] for p in available],
                            complete_repetitions=complete,median_paired_error_delta=value if complete else None,
                            observed_only_median_paired_error_delta=value,
                            fully_validated_median_paired_error_delta=(statistics.median(p['error_delta'] for p in validated)
                                if len(validated)==repetitions else None),
                            definition='Subtract matching repetition errors first; then median of those paired deltas. No repetition selected.')
                        for field in ('CPU_delta','wall_delta'):
                            values=[p[field] for p in group if p[field] is not None]
                            median['median_paired_'+field]=(statistics.median(values) if len(values)==repetitions else None)
                        medians.append(median)
                        if complete: complete_groups.append(value)
            known=[p for p in observations if p['error_delta'] is not None]
            counts=Counter(p['classification'] for p in observations)
            summaries.append(dict(window=window,baseline=baseline,formula=2,expected_pairs=len(observations),
                physically_checked_pairs=len(known),fully_validated_pairs=sum(p['fully_validated'] for p in observations),
                win=counts['win'],tie=counts['tie'],loss=counts['loss'],censored=counts['censored'],
                complete_repetition_groups=len(complete_groups),expected_repetition_groups=len(cases)*len(TIMES)*len(CAPS),
                complete_panel_median_of_paired_deltas=statistics.median(p['error_delta'] for p in known) if len(known)==len(observations) else None,
                ratios_retained=sum(p['baseline_over_window_error_ratio'] is not None for p in observations),
                scope='Descriptive fixed-instance outcomes; two fresh processes are not chemistry samples or a confidence interval.'))
    return pairs,medians,summaries


def pareto_sets(raw,cases,repetitions):
    output=[]
    for case in cases:
        for rep in range(1,repetitions+1):
            for t in TIMES:
                for cap in CAPS:
                    group=[r for r in raw if (r['case'],r['repetition'],r['time'],r['budget_multiplier'])==(case,rep,t,cap)]
                    record=dict(case=case,repetition=rep,time=t,budget_multiplier=cap,expected_paths=11,sets={})
                    for axis,field in (('CPU_seconds','constructor_CPU_seconds'),('wall_seconds','constructor_wall_seconds')):
                        eligible=[r for r in group if r['physical_checked'] and r[field] is not None]
                        dominated={};points=[]
                        for row in eligible:
                            error=row['metrics']['projective_error'];cost=row[field]
                            dominators=[other['path'] for other in eligible if other['path']!=row['path'] and
                                other[field]<=cost and other['metrics']['projective_error']<=error+EPS and
                                (other[field]<cost or other['metrics']['projective_error']<error-EPS)]
                            dominated[row['path']]=dominators
                            points.append(dict(path=row['path'],cost=cost,projective_error=error,physical_floor=error<=EPS,
                                fully_validated=row['fully_validated'],objective_status=None if row['objective'] is None else row['objective']['status']))
                        record['sets'][axis]=dict(points=points,pareto_paths=[p['path'] for p in points if not dominated[p['path']]],
                            dominated_by=dominated,excluded=[dict(path=r['path'],reason='Unverified physical row or unknown whole-call cost')
                                for r in group if r not in eligible],
                            dominance='Recorded cost <= and error <= +1e-12; strictly lower recorded cost or error lower by >1e-12. No measurement uncertainty/optimality guarantee.')
                    output.append(record)
    return output


def analyze(root,mode):
    reader=Artifacts(root)
    index_path=reader.path('ANALYSIS_INDEX.'+mode+'.json')
    _,raw=reader.read_path(index_path);index=decode_json(raw)
    require(index.get('schema')==SCHEMA and index.get('mode')==mode,'Analysis index schema/mode')
    repetitions=2 if mode=='main' else 1
    summary=dict(schema='prospective-window-report-data-v1',mode=mode,status='partial_or_censored',
        analysis_index_sha256=reader.hashes[reader.relative(index_path)]['sha256'],failures=[],cases=[],raw_rows=[],constructor_costs=[],
        ambiguity=EPS,ratio_floor=EPS,source_artifact_SHA256_map={},
        analyzer=dict(sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),frozen_with_scientific_sources=False,
            scope='Separate evidence/report-data tool; source-only numerical aggregation, independently recomputed by study owner'),
        interpretation=dict(independent_replay_scope='Every distinct hard T1.6/24U identity only; all other conditions independently check saved-state metrics',
            repetitions='Two fresh-process repetitions, not independent chemistry observations or population confidence intervals',
            ratio_policy='Suppress at/below1e-12 physical error or unverified/unknown operands; retain floor and loss records',
            PySCF_memory='PySCF max_memory=256 is its documented MB soft guidance, recorded by internal_PySCF_max_memory_MB. Any obsolete historical MiB-named label is not an exact unit conversion or hard limit. Owning cgroup512MiB/no-swap is authoritative.',
            limits='No time/accuracy universality, finite-error certificate, compiled-gate claim, or confidence interval. Pareto membership is descriptive for recorded valid costs.'))
    cases=[];models={};input_complete=False;global_map={};cost_map={};partial_loaded=False
    try:
        protocol,registry,models,cases,input_complete=source_and_inputs(reader,index,mode,summary)
        expected_calls,plan_map,global_map=global_inventory(reader,index,protocol,cases,models,mode,summary)
        supplied=index.get('cases',[]);require(isinstance(supplied,list) and len(supplied)<=len(cases),'Case analysis-index bound')
        bycase={entry['case']:entry for entry in supplied}
        require(len(bycase)==len(supplied) and set(bycase)<=set(cases),'Duplicate/foreign analysis case')
        for case in cases:
            if case not in models:
                result=dict(case=case,status='input_censored',failures=[dict(message='Authenticated model/input unavailable')])
                costs={};candidates={};physical={};independent={};objectives={};references=[]
            else:
                result,costs,candidates,physical,independent,objectives,references=case_artifacts(reader,bycase.get(case),index,
                    case,models[case],repetitions,global_map,protocol,summary)
            summary['cases'].append(result)
            for key,value in costs.items(): cost_map[(case,*key)]=value
            for key in expected_keys(repetitions):
                summary['raw_rows'].append(row_record(case,key,result,costs.get(key[:2]),candidates.get(key),physical.get(key),independent.get(key),objectives.get(key)))
        partial_costs(reader,index,cases,repetitions,global_map,cost_map)
        partial_loaded=True
        # Fill constructor provenance for censored raw cells after partial ledgers.
        for row in summary['raw_rows']:
            cost=cost_map.get((row['case'],row['arm'],row['repetition']))
            if cost is not None:
                row['constructor_CPU_seconds']=cost['CPU_seconds'];row['constructor_wall_seconds']=cost['wall_seconds']
                ctor=cost.get('constructor')
                if ctor is not None:
                    row['constructor_status']=ctor.get('status');row['constructor_window_status']=ctor.get('window_status')
                    row['constructor_scored']=ctor.get('score') is not None
        pairs,medians,comparisons=paired_derivations(summary['raw_rows'],cases,repetitions)
        summary.update(paired_observations=pairs,paired_repetition_medians=medians,comparisons=comparisons,
            pareto_memberships=pareto_sets(summary['raw_rows'],cases,repetitions))
        expected_count=len(cases)*44*repetitions
        verified=sum(row['fully_validated'] for row in summary['raw_rows'])
        physical_count=sum(row['physical_checked'] for row in summary['raw_rows'])
        constructor_complete=(set(cost_map)==expected_calls and all(c['accounting_complete'] for c in cost_map.values()))
        all_artifacts=summary.get('whole_orders',{}).get('all_calls_have_artifacts') is True
        complete=(input_complete and bool(global_map) and all_artifacts and physical_count==expected_count and
                  verified==expected_count and constructor_complete and len(bycase)==len(cases) and
                  all(case['status']=='authenticated' and not case['failures'] and case['primary_panel_complete'] and
                      case['independent_panel_complete'] and case['objective_claims_complete'] for case in summary['cases']))
        summary['completion']=dict(expected_logical_rows=expected_count,authenticated_independent_physical_rows=physical_count,
            fully_validated_rows=verified,missing_or_censored_rows=expected_count-physical_count,
            objective_unverified_rows=physical_count-verified,expected_constructor_calls=len(expected_calls),
            known_constructor_calls=sum(c['accounting_complete'] for c in cost_map.values()),
            source_and_inputs_complete=input_complete,whole_call_artifacts_present=all_artifacts,
            complete_352_main_claim_allowed=complete and mode=='main' and expected_count==352)
        if complete: summary['status']='complete_declared_panel'
    except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
        summary['failures'].append(issue(error))
        summary['completion']=dict(expected_logical_rows=352 if mode=='main' else 44,
            complete_352_main_claim_allowed=False,source_or_artifact_validation_failed=True)
    if not cases:
        cases=list(MAIN_CASES) if mode=='main' else [index.get('cases',[{}])[0].get('case','unperformed_preflight')
            if index.get('cases') else 'unperformed_preflight']
    observed={(r['case'],r['arm'],r['repetition'],r['formula'],r['time'],r['budget_multiplier']) for r in summary['raw_rows']}
    for case in cases:
        for key in expected_keys(repetitions):
            if (case,*key) not in observed:
                invalid=dict(status='artifact_validation_failed',failures=summary['failures'])
                summary['raw_rows'].append(row_record(case,key,invalid,None,None,None,None,None))
    if not partial_loaded:
        try: partial_costs(reader,index,cases,repetitions,global_map,cost_map)
        except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
            summary['failures'].append(dict(phase='partial_costs',**issue(error)))
    for row in summary['raw_rows']:
        cost=cost_map.get((row['case'],row['arm'],row['repetition']))
        if cost is not None:
            row['constructor_CPU_seconds']=cost['CPU_seconds'];row['constructor_wall_seconds']=cost['wall_seconds']
            if cost.get('constructor') is not None:
                row['constructor_status']=cost['constructor'].get('status')
                row['constructor_window_status']=cost['constructor'].get('window_status')
                row['constructor_scored']=cost['constructor'].get('score') is not None
    # Rebuild after partial-cost recovery so every derived cost uses the final ledger.
    pairs,medians,comparisons=paired_derivations(summary['raw_rows'],cases,repetitions)
    summary.update(paired_observations=pairs,paired_repetition_medians=medians,comparisons=comparisons,
        pareto_memberships=pareto_sets(summary['raw_rows'],cases,repetitions))
    summary['constructor_costs']=[cost_map[key] for key in sorted(cost_map)]
    expected_total=(64 if mode=='main' else 8)
    known=[cost for cost in cost_map.values() if cost['accounting_complete']]
    summary['construction_totals']=dict(expected_calls=expected_total,known_calls=len(known),
        CPU_seconds=sum(c['CPU_seconds'] for c in known) if len(known)==expected_total else None,
        wall_seconds=sum(c['wall_seconds'] for c in known) if len(known)==expected_total else None,
        known_CPU_seconds=sum(c.get('known_CPU_seconds',0.) for c in cost_map.values()),
        known_wall_seconds=sum(c.get('known_wall_seconds',0.) for c in cost_map.values()),
        cost_scope='Whole-call guard receipts per arm/repetition including every PANEL attempt; never summed over repeated physical rows')
    summary['construction_by_arm']={}
    for arm in ARMS:
        calls=[cost for cost in cost_map.values() if cost['arm']==arm]
        known_calls=[cost for cost in calls if cost['accounting_complete']]
        expected_arm=len(cases)*repetitions
        summary['construction_by_arm'][arm]=dict(expected_calls=expected_arm,indexed_calls=len(calls),
            known_calls=len(known_calls),CPU_seconds=sum(c['CPU_seconds'] for c in known_calls) if len(known_calls)==expected_arm else None,
            wall_seconds=sum(c['wall_seconds'] for c in known_calls) if len(known_calls)==expected_arm else None,
            known_CPU_seconds=sum(c.get('known_CPU_seconds',0.) for c in calls),
            known_wall_seconds=sum(c.get('known_wall_seconds',0.) for c in calls),
            call_status_counts=dict(Counter(c['status'] for c in calls)))
    phase=[];seen=set();construction_hashes={attempt['descriptor']['sha256'] for call in cost_map.values() for attempt in call['attempts']}
    phase_records=index.get('phase_receipts',[])
    if not isinstance(phase_records,list) or len(phase_records)>2048:
        summary['failures'].append(dict(phase='phase_receipts',message='Phase receipt bound'));phase_records=[]
    for item in phase_records:
        try:
            require(isinstance(item,dict) and isinstance(item.get('phase'),str),'Phase receipt metadata')
            measurement=receipt(reader,item['receipt']);sha256=item['receipt']['sha256']
            require(sha256 not in seen,'Phase receipt duplicated');seen.add(sha256)
            phase.append(dict(phase=item['phase'],case=item.get('case'),measurement=measurement,
                already_counted_in_construction=sha256 in construction_hashes))
        except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
            summary['failures'].append(dict(phase='phase_receipts',**issue(error)))
    summary['other_phase_costs']=phase
    summary['other_phase_totals']={}
    for name in sorted({entry['phase'] for entry in phase}):
        group=[entry['measurement'] for entry in phase if entry['phase']==name and not entry['already_counted_in_construction']]
        complete=all(entry['accounting_complete'] for entry in group)
        summary['other_phase_totals'][name]=dict(indexed_attempts=len(group),accounting_complete=complete,
            CPU_seconds=sum(entry['CPU_seconds'] for entry in group) if complete else None,
            wall_seconds=sum(entry['wall_seconds'] for entry in group) if complete else None,
            known_CPU_seconds=sum(entry['CPU_seconds'] for entry in group if entry['CPU_seconds'] is not None),
            known_wall_seconds=sum(entry['wall_seconds'] for entry in group if entry['wall_seconds'] is not None),
            scope='Provided unique phase-receipt inventory only; no unprovided phase cost is assumed zero')
    historical=index.get('historical_artifacts',[])
    if not isinstance(historical,list) or len(historical)>1024:
        summary['failures'].append(dict(phase='historical_artifacts',message='Historical artifact bound'));historical=[]
    summary['historical_artifacts']=[]
    for item in historical:
        try:
            kind=item.get('format','binary');require(kind in ('json','toml','binary'),'Historical artifact format')
            _,value=reader.read(item['artifact'],kind)
            summary['historical_artifacts'].append(dict(phase=item.get('phase'),case=item.get('case'),
                descriptor=item['artifact'],record=value,scope='Historical/refused preparation artifact; not mixed into current physical targets'))
        except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
            summary['failures'].append(dict(phase='historical_artifacts',**issue(error)))
    if summary['failures'] or any(not entry['measurement']['accounting_complete'] for entry in phase):
        summary['status']='partial_or_censored';summary['completion']['complete_352_main_claim_allowed']=False
    summary['source_artifact_SHA256_map']=reader.hashes
    summary['bytes_hashed']=reader.bytes_hashed
    return summary


def write_new(path,value):
    path=Path(path);partial=path.with_name(path.name+'.partial')
    require(not path.exists() and not partial.exists(),'Append-only new summary/partial required')
    def serializable(value):
        if isinstance(value,float) and not math.isfinite(value):
            return dict(nonfinite_quarantined_value=repr(value))
        if isinstance(value,dict): return {key:serializable(item) for key,item in value.items()}
        if isinstance(value,(tuple,list)): return [serializable(item) for item in value]
        return value
    raw=(json.dumps(serializable(value),indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    require(len(raw)<=MAX_OUTPUT,'Summary output byte ceiling')
    path.parent.mkdir(parents=True,exist_ok=True)
    with partial.open('xb') as stream: stream.write(raw)
    os.link(partial,path);partial.unlink()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--mode',choices=('main','preflight'),default='main')
    args=parser.parse_args(argv)
    result=analyze(args.root,args.mode);write_new(args.output,result)
    return 0 if result['status']=='complete_declared_panel' else 2


if __name__=='__main__':
    raise SystemExit(main())
