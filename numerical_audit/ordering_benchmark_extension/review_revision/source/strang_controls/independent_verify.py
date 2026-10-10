#!/usr/bin/env python3
"""Independent all-row replay with physical q0-MSB HF convention.

This supersedes only the HF guard in the original frozen run.py verifier.
Production sources, schedules and propagated states are unchanged. See
CHECKER_CORRECTION.md for the exact correction and provenance.
"""
import os,resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,gc,json,sys,time
from pathlib import Path
import run as r
ROOT=r.ROOT
CASES=r.CASES
BUDGETS=r.BUDGETS
require=r.require
sha=r.sha
write=r.write
sources=r.sources
verified_design=r.verified_design
def verify(cases):
    import numpy as np
    import independent_core as ind
    require(not any(n=='qhat' or n.startswith('qhat.') for n in sys.modules),'Independent checker imported QHAT')
    started=time.monotonic();allreports=[]
    for case in cases:
        d=verified_design(case);model=ind.OccupationAction(d)
        with np.load(ROOT/'inputs'/f'{case}.npz') as z:initial=z['initial'].copy();exact=z['exact'].copy()
        fresh,nnz=model.reference(initial);reference_distance=float(np.linalg.norm(fresh-exact));require(reference_distance<3e-11,'Independent exact reference mismatch')
        require(np.count_nonzero(initial)==1 and int(model.basis[np.argmax(abs(initial))])==((1<<int(d['electrons']))-1)<<(int(d['n_qubits'])-int(d['electrons'])),'Not frozen HF')
        records=[]
        for name,method in d['methods'].items():
            require(ind.schedule_sum(method['schedule'])==model.total,'Target sum mismatch')
            schedule=method['schedule'];offset=0
            for length in method['unit_lengths']:
                unit=schedule[offset:offset+length];offset+=length
                if 'refined' in name:require(not ind.number_commutator(unit),'Conserving unit violates N')
            for b in BUDGETS:
                rp=ROOT/'rows'/case/f'{name}.b{b}.json';row=json.loads(rp.read_text());sp=ROOT/row['state_file'];require(sha(sp)==row['state_sha256'],'State hash mismatch');out=ROOT/'verification'/case/f'{name}.b{b}.json'
                if out.exists():
                    old=json.loads(out.read_text());require(old['row_sha256']==sha(rp) and old.get('verification_source_sha256')==sha(__file__),'Verification row/source changed');records.append(old);continue
                with np.load(sp) as z:actual=z['actual'].copy()
                before=time.monotonic();replay=model.replay(initial,schedule,row['steps'],2);delta=float(np.linalg.norm(replay-actual));seconds=time.monotonic()-before;require(delta<3e-11,'Unmerged independent replay mismatch')
                metrics=model.metrics(fresh,actual,d['electrons']);differences={k:abs(float(row[k])-v) for k,v in metrics.items() if k in row};require(max(differences.values(),default=0)<3e-11,'Independent metric mismatch');require(abs(metrics['norm']-1)<2e-9,'State normalization drift')
                record={'case':case,'ordering':name,'budget_units':b,'state_distance':delta,'reference_distance':reference_distance,'metric_differences':differences,'independent_metrics':metrics,'seconds':seconds,'reused':row['reused'],'row_sha256':sha(rp),'state_sha256':sha(sp),'source_sha256':sources(),'verification_source_sha256':sha(__file__)};write(out,record);records.append(record)
                print('VERIFY',case,name,b,delta,'seconds',seconds,flush=True)
        report={'case':case,'rows':len(records),'replays':len(records),'new_state_replays':sum(not r['reused'] for r in records),'reference_distance':reference_distance,'reference_nnz':nnz,'max_state_distance':max(r['state_distance'] for r in records),'max_metric_difference':max(max(r['metric_differences'].values(),default=0) for r in records),'records':records};write(ROOT/'verification'/f'{case}.json',report);allreports.append(report);gc.collect()
    reports=[json.loads((ROOT/'verification'/f'{case}.json').read_text()) for case in CASES if (ROOT/'verification'/f'{case}.json').exists()]
    if len(reports)==14:write(ROOT/'verification.json',{'status':'passed','cases':14,'replays':sum(r['replays'] for r in reports),'new_state_replays':sum(r['new_state_replays'] for r in reports),'max_state_distance':max(r['max_state_distance'] for r in reports),'max_metric_difference':max(r['max_metric_difference'] for r in reports),'max_reference_distance':max(r['reference_distance'] for r in reports),'seconds':time.monotonic()-started,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'verification_source_sha256':sha(__file__),'scope':'All 420 new states and 84 reused baseline states replayed with independently derived occupation-basis Pauli action and unmerged Strang palindromes; all 14 reference states recomputed independently.'})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cases',nargs='*',default=CASES);a=p.parse_args();require(set(a.cases)<=set(CASES),'Unknown cases');verify(a.cases)
