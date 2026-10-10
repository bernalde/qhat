#!/usr/bin/env python3
"""Resumable published-comparator campaign; all plans freeze before propagation."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,math,resource,sys,time
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import *
import algorithms,extras

def fingerprint(d):
    co=[int(c) for w,c in d['coefficients']];g=math.gcd(*map(abs,co));return object_sha([[w,str(int(c)//g)] for w,c in d['coefficients']]),g

def prepare(case,campaign):
    base=HERE/campaign;p=base/'designs'/f'{case}.json'
    if p.exists():return json.loads(p.read_text())
    started=time.monotonic();d,states=load_case(case);total=block(d['coefficients'])
    if campaign=='extra':methods=extras.build_extra_schedules(d)
    else:
        fp,g=fingerprint(d);cache=base/'scale_cache'/f'{fp}.json'
        if cache.exists():
            prior=json.loads(cache.read_text());oldg=int(prior['gcd']);require(g%oldg==0,'Nonintegral uniform coefficient scaling');factor=g//oldg;m=prior['method'];m['schedule']=[(word(w),int(c)*factor) for w,c in m['schedule']];meta=m['metadata'];meta['history']=[[i,p,str(int(v)*factor**3)] for i,p,v in meta['history']];meta['exact_uniform_scaling_reuse']={'source_case':prior['case'],'integer_multiplier':factor,'normalized_coefficients_sha256':fp};meta['latest_native_seconds']=0.
        else:
            m=extras.native_v2(d,base/'native'/case);write_atomic(cache,{'case':case,'gcd':str(g),'method':{'schedule':encode_schedule(m['schedule']),'metadata':m['metadata']}})
        methods={'tranter_v2_l1_first_order':m,'tranter_v2_l1_strang':{'schedule':m['schedule'],'metadata':{**m['metadata'],'formula':'chronological Strang; adjacent identical boundary and central Pauli rotations merged'},'formula_order':2}}
    frozen={}
    for name,m in methods.items():
        accum={}
        for w,c in m['schedule']:accum[w]=accum.get(w,0)+c
        require({w:c for w,c in accum.items() if c}==total,'Schedule changes target '+name)
        enc=encode_schedule(m['schedule']);frozen[name]={'schedule':enc,'metadata':m['metadata'],'formula_order':m.get('formula_order',1),'terms':len(enc),'schedule_sha256':object_sha(enc)}
    d['methods']=frozen;d['extension_panel']='development' if case in DEV_CASES else 'prospective_followup';d['extension_preparation']={'seconds':time.monotonic()-started};write_atomic(p,d);print('PREPARED',campaign,case,round(time.monotonic()-started,3),flush=True);return d

def freeze(campaign):
    base=HERE/campaign
    record={'cases':CASES,'design_sha256':{case:sha(base/'designs'/f'{case}.json') for case in CASES},'source_sha256':{str(p.relative_to(HERE.parent)):sha(p) for p in (Path(__file__),HERE/'algorithms.py',HERE/'extras.py',HERE/'v2_exact.cpp',HERE/'PROTOCOL.md',HERE.parent/'common.py')},'expected_rows':len(CASES)*(4 if campaign=='extra' else 2)*len(BUDGETS),'budget_units':BUDGETS,'status':'all schedules frozen before first campaign propagation'}
    p=base/'freeze.json'
    if p.exists():require(json.loads(p.read_text())==record,'Freeze changed')
    else:write_atomic(p,record)

def run_case(case,campaign):
    import numpy as np
    from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics
    base=HERE/campaign;frozen=json.loads((base/'freeze.json').read_text());p=base/'designs'/f'{case}.json';require(sha(p)==frozen['design_sha256'][case],'Design mutated')
    for rel,h in frozen['source_sha256'].items():require(sha(HERE.parent/rel)==h,'Source mutated '+rel)
    d=json.loads(p.read_text());raw,states=load_case(case);require(d['extension_input_provenance']==raw['extension_input_provenance'],'Input changed');comp=compile_case(d,states);U=len(d['coefficients']);out=base/'results';out.mkdir(exist_ok=True);rows=[]
    for name,m in d['methods'].items():
        seq=decode_schedule(m['schedule']);terms=[comp['term'](w,c) for w,c in seq];R=len(terms);order=m['formula_order']
        for b in BUDGETS:
            budget=U*b;r=budget//R if order==1 else (budget-1)//(2*(R-1));rot=r*R if order==1 else (2*r*(R-1)+1 if r>0 else 0);rp=out/f'{case}.{name}.b{b}.json';sp=rp.with_suffix('.npz')
            if rp.exists():
                row=json.loads(rp.read_text());require(row['extension_design_sha256']==frozen['design_sha256'][case] and row['schedule_sha256']==m['schedule_sha256'],'Checkpoint mismatch')
                if row['status']=='complete':require(sha(sp)==row['state_sha256'],'State checkpoint corrupt')
                rows.append(row);continue
            row={'case':case,'panel':d['extension_panel'],'ordering':name,'T':1.,'formula_order':order,'budget_units':b,'rotation_budget':budget,'steps':r,'role':'matched_rotation_budget','extension_design_sha256':frozen['design_sha256'][case],'schedule_sha256':m['schedule_sha256'],'pauli_rotations':rot}
            if not r:row['status']='infeasible';write_atomic(rp,row);rows.append(row);continue
            started=time.monotonic()
            if order==1:actual=evolve(comp['initial'],terms,1.,r)
            else:
                # Merge center and inter-step boundary rotations physically.
                pulses=[]
                def emit(i,factor):
                    t=terms[i];pulses.append((t[0]*factor/(2*r),t[1],t[2],t[3]))
                emit(0,1)
                for step in range(r):
                    for i in range(1,R-1):emit(i,1)
                    emit(R-1,2)
                    for i in range(R-2,0,-1):emit(i,1)
                    emit(0,1 if step==r-1 else 2)
                require(len(pulses)==rot,'Strang cost mismatch');actual=evolve(comp['initial'],pulses,1.,1)
            row.update(status='complete',seconds=time.monotonic()-started,state_file=sp.name,state_key='actual',**physical_metrics(comp['exact'],actual,comp['number_mask'],comp['spin_mask']));tmp=sp.with_suffix('.partial.npz');np.savez_compressed(tmp,actual=actual);os.replace(tmp,sp);row['state_sha256']=sha(sp);write_atomic(rp,row);rows.append(row);print('EVOLVE',campaign,case,name,b,row['one_minus_overlap'],round(row['seconds'],3),flush=True)
    write_atomic(out/f'{case}.completion.json',{'rows':len(rows),'status':'complete'})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare','freeze','run','collect'));ap.add_argument('campaign',choices=('extra','v2'));ap.add_argument('--cases',nargs='*',default=CASES);args=ap.parse_args();resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base=HERE/args.campaign
    if args.action=='prepare':
        for case in args.cases:prepare(case,args.campaign)
    elif args.action=='freeze':freeze(args.campaign)
    elif args.action=='run':
        for case in args.cases:run_case(case,args.campaign)
    else:
        rows=[]
        for case in CASES:
            require((base/'results'/f'{case}.completion.json').exists(),'Incomplete '+case)
            rows.extend(json.loads(p.read_text()) for p in sorted((base/'results').glob(case+'.*.b*.json')))
        require(len(rows)==len(CASES)*(4 if args.campaign=='extra' else 2)*7,'Wrong row count');write_atomic(base/'results.json',rows);write_atomic(base/'completion.json',{'rows':len(rows),'complete':sum(r['status']=='complete' for r in rows),'infeasible':sum(r['status']=='infeasible' for r in rows)})
if __name__=='__main__':main()
