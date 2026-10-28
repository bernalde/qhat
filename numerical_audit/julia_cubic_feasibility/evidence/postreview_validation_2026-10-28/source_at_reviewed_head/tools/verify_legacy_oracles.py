"""Independent checks using retained legacy exact/matrix/replay oracles.

This invokes existing verified Python implementations unchanged. New solvers
and physical propagation in the delivered prototype are Julia.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import argparse
import importlib.util
import json
import sys
import time
import tomllib

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('output',type=Path)
parser.add_argument('--study-root',type=Path,default=Path(__file__).resolve().parents[1])
parser.add_argument('--repo-root',type=Path)
args=parser.parse_args()
ROOT=args.study_root.resolve()
WORK=args.repo_root.resolve() if args.repo_root else ROOT.parents[1]
RESULTS=ROOT/'evidence/results'
LEGACY_MODELS=ROOT/'evidence/inputs/legacy_models'
OLD=WORK/'numerical_audit/exact_response_windows/evidence/results'

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    assert Path(module.__file__).resolve()==path.resolve()
    return module

def toml(path):return tomllib.loads(Path(path).read_text())

def record(model,order,formula=2):
    zero=[i-1 for i in order]
    schedule=[];lengths=[]
    for i in zero:
        block=model['blocks'][i];lengths.append(len(block))
        schedule.extend([[word,c] for word,c in sorted(block)])
    return dict(order=zero,formula=formula,partition='conserving',schedule=schedule,unit_lengths=lengths)

def close(a,b,label):
    assert abs(a-b)<=2e-12,(label,a,b,abs(a-b))

began=time.process_time()
out=args.output;assert not out.exists()
exact=load_module('archived_exact_oracle',WORK/'numerical_audit/exact_response_windows/exact.py')
old_verify=load_module('archived_physical_oracle',WORK/'numerical_audit/exact_response_windows/verify.py')
archive=toml(RESULTS/'ARCHIVE_MAIN.a1.toml')
physical_archive=toml(RESULTS/'PHYSICAL_ARCHIVE.a1.toml')
checked_methods=checked_old_physical=0
for case in archive['cases']:
    name=case['case'];domain=json.loads((OLD/'domains'/f'{name}.json').read_text())
    for sample in case['first_calls']+case['warm_samples']:
        solver=sample['solver'];assert solver['complete']
        full=solver['validation']['full_order']
        value=exact.objective(full,domain)
        assert value==int(case['expected_integer_minimum'])
        assert Fraction(solver['cost'].replace('//','/'))==value
        assert full[:domain['offset']]==domain['outside_prefix']
        assert full[domain['offset']+domain['width']:]==domain['outside_suffix']
        checked_methods+=1
    rows=next(c['rows'] for c in physical_archive['cases'] if c['case']==name)
    for row in rows:
        if row['name']=='dp':continue
        index=(0 if row['name']=='parent' else 1)+(4 if row['formula']==2 else 0)
        old=json.loads((OLD/'rows'/f'{name}.{index}.b{row["budget_multiplier"]}.json').read_text())
        close(row['projective_error'],old['error'],name+' archived physical')
        assert row['steps']==old['steps'] and row['ideal_rotations']==old['rotations']
        checked_old_physical+=1

common=load_module('common',WORK/'numerical_audit/scalable_ordering/common.py')
reference_module=load_module('large_literal_reference_oracle',WORK/'numerical_audit/scalable_ordering/reference.py')
physical=toml(RESULTS/'PHYSICAL_CUBIC.a1.toml')
objective_checks=[];physical_checks=[]
for short,name in (('CH4','CH4_r1p09'),('CO','CO_r1p13')):
    model=json.loads((LEGACY_MODELS/f'{name}.json').read_text())
    window=toml(RESULTS/f'{short}_WINDOW.a2.toml')
    control=toml(RESULTS/f'{name}.CONTINUATION.a1.toml')
    assert hashlib.sha256((LEGACY_MODELS/f'{name}.json').read_bytes()).hexdigest()==control['model_sha256']
    ref=reference_module.Reference(model)
    exact_state=ref.exact_action(1.5)
    rows=next(c['rows'] for c in physical['cases'] if c['case']==name)
    for arm,order,wanted in (('window_exact',window['new_full_order'],window['new_J']),
                              ('legacy_continuation',control['order'],control['independently_checked_score'])):
        rec=record(model,order)
        observed=ref.objective(rec,kind='cubic')['value']
        assert abs(observed-wanted)<=max(2e-18,2e-7*abs(wanted)),(name,arm,'objective',observed,wanted)
        objective_checks.append(dict(case=name,arm=arm,value=observed,expected=wanted))
        row=next(r for r in rows if r['name']==arm and r['time']==1.5 and r['budget_multiplier']==10)
        actual=ref.ideal_formula(rec,row['steps'],1.5)
        metric=old_verify.physical_metrics(exact_state,actual,model)
        close(metric['error'],row['projective_error'],name+' '+arm+' physical')
        close(metric['number_leakage'],row['number_leakage'],name+' leakage')
        physical_checks.append(dict(case=name,arm=arm,error=metric['error'],julia_error=row['projective_error'],absolute_difference=abs(metric['error']-row['projective_error'])))
    checkpoint=out.with_name(out.name+'.'+name+'.json')
    checkpoint.write_text(json.dumps(dict(objective_checks=objective_checks,physical_checks=physical_checks),indent=2)+'\n')
    del ref
summary=dict(status='passed',archived_method_outputs_checked=checked_methods,
    archived_physical_rows_checked_against_retained_independent_results=checked_old_physical,
    fresh_literal_nested_commutator_objectives=objective_checks,
    fresh_unfused_full_space_physical_replays=physical_checks,CPU_seconds=time.process_time()-began,
    scope='Exact archived certificates allcases/methods; fresh matrix/replay checks of window and whole-order continuation atT1.5cap10U forboth12modecases. Othernewphysicalsettings useJulia validatedkernel; no claim of independent replay of all68rows.')
out.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ('status','archived_method_outputs_checked','archived_physical_rows_checked_against_retained_independent_results','CPU_seconds')}))
