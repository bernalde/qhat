#!/usr/bin/env python3
"""Audit all saved rows and replay a predeclared sample independently.

All small closures (dimension <=256) are replayed fully. For larger cases,
replay every deterministic b1/b2 circuit, the first random seed of each family,
and optimized first-order b10 and Pauli-JD b100. Selection uses no errors.
Incomplete campaigns fail rather than producing a successful partial summary.
"""
import argparse,json,resource,time
from pathlib import Path
from core import *

def select_replay(row,dim):
    if dim<=256:return True
    name=row['ordering'];b=row['budget_units']
    if name.startswith('random_'):return name.endswith('_00')
    if b in (1,2):return True
    if name in ('class_refined_response','class_refined_J0','class_pauli_JD','class_pauli_JD_budget') and b==10:return True
    return name=='class_pauli_JD' and b==100

def get_design(root,case):
    for suffix in ('.json','.design.json'):
        p=root/'designs'/(case+suffix)
        if p.exists():return p,json.loads(p.read_text())
    raise FileNotFoundError('Missing design '+case)

def get_rows(root,case):
    paths=list((root/'results').glob(case+'.*.b*.json'))+list((root/'rows'/case).glob('*.json'))
    require(paths,'No rows '+case)
    return [(p,json.loads(p.read_text())) for p in sorted(paths)]

def state_path(root,rowpath,row):
    p=root/row['state_file']
    if not p.exists():p=rowpath.parent/row['state_file']
    require(p.exists(),'Missing state '+str(p));return p

def check_case(root,case,reference_dir,output):
    dp,d=get_design(root,case);pairs=get_rows(root,case);dh=sha(dp)
    with np.load(reference_dir/(case+'.npz')) as f:initial=f['initial'];exact=f['exact']
    action=OccupationAction(d);require(initial.shape==(action.dim,),'State dimension mismatch')
    methods=d['methods'];U=len(d['coefficients']);records=[];started=time.monotonic();seen=set()
    for name,m in methods.items():
        require(object_sha(m['schedule'])==m['schedule_sha256'],'Schedule digest mismatch '+name)
        require(schedule_sum(m['schedule'])==action.total,'Schedule changes exact target '+name)
    for rowpath,row in pairs:
        name=row['ordering'];b=int(row['budget_units']);token=(name,b)
        require(token not in seen,'Duplicate result');seen.add(token)
        require(name in methods,'Unknown method');m=methods[name];r=int(row['steps']);formula=int(row['formula_order']);schedule=m['schedule']
        require(row['extension_design_sha256']==dh,'Row design hash mismatch')
        require(row['schedule_sha256']==m['schedule_sha256'],'Row schedule hash mismatch')
        require(float(row['T'])==1. and formula==m.get('formula_order',1),'Target time/formula mismatch')
        budget=U*b;require(row['rotation_budget']==budget,'Wrong budget')
        if formula==1:expected_steps=budget//len(schedule);rotations=r*len(schedule)
        else:
            lengths=m.get('unit_lengths',[1]*len(schedule));first,last=lengths[0],lengths[-1]
            expected_steps=(budget-first)//(2*len(schedule)-first-last);expected_steps=max(0,expected_steps)
            rotations=r*(2*len(schedule)-first-last)+first if r else 0
        require(r==expected_steps and rotations==row['pauli_rotations'] and rotations<=budget,'Incorrect matched cost/step count '+str(token))
        if row['status']=='infeasible':
            require(r==0,'Feasible row marked infeasible');records.append({'ordering':name,'budget_units':b,'status':'infeasible','row_sha256':sha(rowpath)});continue
        require(row['status']=='complete' and r>0,'Bad row status')
        sp=state_path(root,rowpath,row);require(sha(sp)==row['state_sha256'],'State hash mismatch')
        with np.load(sp) as f:
            actual=f[row.get('state_key','actual')]
            require(actual.shape==initial.shape and np.all(np.isfinite(actual)),'Malformed output state')
            if formula==2 and 'pulse_term_indices' in f:
                ids=f['pulse_term_indices'];nums=f['pulse_duration_numerators'];den=int(f['pulse_duration_denominator'])
                require(len(ids)==rotations and len(nums)==rotations,'Wrong pulse array lengths')
                require(np.all((0<=ids)&(ids<len(schedule))) and np.all(nums>0) and den==2*r,'Bad pulse indices/durations')
                require(np.all(np.bincount(ids,weights=nums,minlength=len(schedule))==den),'Pulse duration does not sum to T')
        metrics=action.metrics(exact,actual,int(d['electrons']))
        differences={k:abs(metrics[k]-float(row[k])) for k in ('one_minus_overlap','infidelity','particle_number_leakage','spin_sector_leakage')}
        require(max(differences.values())<2e-12,'Saved physical metric mismatch '+str(token))
        require(abs(metrics['norm']-float(row['approx_norm']))<1e-12 and abs(metrics['norm']-1)<1e-8,'Unacceptable or misreported norm drift')
        delta=None
        if select_replay(row,action.dim):
            replay=action.replay(initial,schedule,r,formula)
            delta=float(np.linalg.norm(replay-actual));require(delta<2e-11,'Independent circuit mismatch '+str(token)+': '+str(delta))
        records.append({'ordering':name,'budget_units':b,'status':'passed','metric_difference':differences,'state_distance':delta,'row_sha256':sha(rowpath),'state_sha256':sha(sp)})
    expected={(name,b) for name,m in methods.items() for b in ((100,) if name.startswith('random_') else ((2,5,10,20,50,100) if name.startswith('strang_') else (1,2,5,10,20,50,100)))}
    # Published V2 Strang includes a deliberate b1 infeasible record.
    extra={(name,1) for name,m in methods.items() if m.get('formula_order',1)==2 and not name.startswith('strang_')}
    expected|=extra
    require(seen==expected,f'Incomplete/extra rows for {case}: missing {expected-seen}, extra {seen-expected}')
    result={'case':case,'status':'passed','rows':len(records),'replayed':sum(x.get('state_distance') is not None for x in records),'design_sha256':dh,'source_sha256':sha(__file__),'core_sha256':sha(Path(__file__).with_name('core.py')),'seconds':time.monotonic()-started,'checks':records}
    p=output/(case+'.json');p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.partial');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(p)
    print('VERIFIED',root.name,case,len(records),'replayed',result['replayed'],'seconds',result['seconds'],flush=True)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--references',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='*');a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    freeze=next((f for f in (a.campaign/'freeze.json',a.campaign/'ALL_SCHEDULES_FROZEN.json') if f.exists()),None);require(freeze is not None,'No frozen campaign manifest')
    frozen=json.loads(freeze.read_text());cases=a.cases or frozen['cases'];records=[]
    extension_root=next(parent for parent in (a.campaign,*a.campaign.parents) if (parent/'common.py').exists())
    for path,h in frozen.get('source_sha256',{}).items():
        full=extension_root/path
        require(full.exists() and sha(full)==h,'Frozen source changed '+path)
    for case in cases:
        dp,d=get_design(a.campaign,case);require(sha(dp)==frozen['design_sha256'][case],'Frozen design changed '+case)
        saved=a.output/(case+'.json')
        if saved.exists():
            record=json.loads(saved.read_text());require(record['source_sha256']==sha(__file__) and record['core_sha256']==sha(Path(__file__).with_name('core.py')) and record['design_sha256']==sha(dp),'Audit checkpoint no longer valid')
            current={(r['ordering'],r['budget_units']):(p,r) for p,r in get_rows(a.campaign,case)}
            require(len(current)==len(record['checks']),'Rows changed after audit')
            for checked in record['checks']:
                p,r=current[(checked['ordering'],checked['budget_units'])];require(sha(p)==checked['row_sha256'],'Audited row changed')
                if 'state_sha256' in checked:require(sha(state_path(a.campaign,p,r))==checked['state_sha256'],'Audited state changed')
            records.append(record)
        else:records.append(check_case(a.campaign,case,a.references,a.output))
    checks=[x for rec in records for x in rec['checks']];result={'status':'passed','campaign':a.campaign.name,'case_count':len(records),'rows':len(checks),'replayed':sum(x.get('state_distance') is not None for x in checks),'max_state_distance':max((x['state_distance'] for x in checks if x.get('state_distance') is not None),default=0.),'max_metric_difference':max((max(x['metric_difference'].values()) for x in checks if 'metric_difference' in x),default=0.),'freeze_sha256':sha(freeze),'replay_selection':__doc__,'cases':cases}
    a.output.mkdir(parents=True,exist_ok=True);(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
