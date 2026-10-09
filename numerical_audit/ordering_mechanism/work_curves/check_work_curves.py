#!/usr/bin/env python3
"""Read-only provenance, fixed-order, metric and threshold audit; no QHAT."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json,hashlib
from pathlib import Path
import numpy as np
def require(b,msg):
    if not b:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def block(raw):return {tuple((int(q),s) for q,s in word):int(v) for word,v in raw}
def schedule(d,name):
    total=block(d['coefficients'])
    if name.startswith('jw_'):return [(tuple((int(q),s) for q,s in word),total[tuple((int(q),s) for q,s in word)]) for word in d['jw_orders'][name]]
    bs=[block(b) for b in d['blocks'][name.split('_')[0]]]
    return [(k,v) for i in d['orders'][name] for k,v in bs[i].items()]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    m=json.loads((a.run/'manifest.json').read_text());cases=m['cases'];names=('class_refined_response','class_lex','owned_signed','jw_signed','jw_magnitude');grid=(1,2,5,10,20,50,100)
    for rel,h in m['source_sha256'].items():require(sha(a.run/rel)==h,'Changed source '+rel)
    require(sha(a.run/'frozen_candidate.py')==sha(a.candidate/'run.py'),'Changed frozen candidate helper')
    rows=json.loads((a.run/'combined_results.json').read_text());fresh=json.loads((a.run/'results.json').read_text())
    require(len(rows)==280 and len(fresh)==240 and rows[:240]==fresh,'Wrong row cardinalities or prefix')
    require(len({(r['case'],r['ordering'],r['steps']) for r in rows})==280,'Duplicate combined rows')
    source_rows={a.candidate:json.loads((a.candidate/'results.json').read_text()),a.baseline:json.loads((a.baseline/'results.json').read_text())}
    source_hash={p:sha(p/'results.json') for p in source_rows};maximum_metric_difference=0.;casechecks=[];arraychecks=0
    for case in cases:
        d=json.loads((a.run/f'{case}.design.json').read_text());cd=json.loads((a.candidate/f'{case}.design.json').read_text());bd=json.loads((a.baseline/f'{case}.design.json').read_text())
        require(d['parent_candidate_design_sha256']==sha(a.candidate/f'{case}.design.json'),'Candidate design provenance mismatch')
        for k in ('coefficients','denominator','electrons','n_qubits','affine_basis_indices','target_sha256','initial_sha256','tensor_sha256'):
            require(d[k]==cd[k]==bd[k],'Target or initial-state identity mismatch '+k)
        require(set(d['orders'])==set(names[:3]) and set(d['jw_orders'])==set(names[3:]),'Inactive orders retained')
        planned={(r['ordering'],r['steps']) for r in d['planned_runs']}
        require(planned=={(n,r) for n in names for r in grid[:-1]},'Changed finite-step grid')
        total=block(d['coefficients']);U=len(total);seq={name:schedule(d,name) for name in names}
        for name in names:
            baseline_design=cd if name=='class_refined_response' else bd
            require(seq[name]==schedule(baseline_design,name),'Frozen word/coefficient order changed: '+name)
            require(len(seq[name])==U and dict(seq[name])==total,'Unequal Pauli multiset')
        archives={p:np.load(p/f'{case}.states.npz',allow_pickle=False) for p in (a.run,a.candidate,a.baseline)}
        try:
            for k in ('initial','exact','exact_dense'):
                require(np.array_equal(archives[a.run][k],archives[a.candidate][k]) and np.array_equal(archives[a.run][k],archives[a.baseline][k]),'Reference arrays changed')
            exact=archives[a.run]['exact'];e=exact/np.linalg.norm(exact)
            for row in (r for r in rows if r['case']==case):
                name,r=row['ordering'],row['steps'];require(r in grid and name in names and row['T']==1.,'Unexpected row')
                require(row['pauli_rotations']==U*r,'Wrong word-rotation cost')
                if r==100:
                    src=a.candidate if name=='class_refined_response' else a.baseline
                    require(row['reused_source_results_sha256']==source_hash[src],'Reused results provenance mismatch')
                    original=next(x for x in source_rows[src] if x['case']==case and x['ordering']==name and x['steps']==100)
                    require({k:v for k,v in row.items() if k!='reused_source_results_sha256'}==original,'Reused row was modified')
                else:src=a.run
                state=archives[src][f"{name}__r{r}__{row['role']}"];v=state/np.linalg.norm(state);z=np.vdot(e,v);aligned=v if abs(z)==0 else v*z.conjugate()/abs(z)
                E=.5*float(np.vdot(e-aligned,e-aligned).real);delta=abs(E-row['one_minus_overlap']);maximum_metric_difference=max(maximum_metric_difference,delta)
                require(delta<1e-16,'Saved-state metric disagreement');arraychecks+=1
        finally:
            for f in archives.values():f.close()
        casechecks.append({'case':case,'paulis_per_step':U,'exact_same_H_state_references':True,'all_word_orders_exactly_match_sources':True})
    table=json.loads((a.run/'thresholds.json').read_text());require(len(table)==24,'Wrong threshold rows')
    threshold_claim=[]
    for row in table:
        require(row['case'] in cases and row['threshold'] in (1e-4,1e-6,1e-8),'Unexpected threshold')
        for name in names:
            admissible=[r for r in rows if r['case']==row['case'] and r['ordering']==name and r['one_minus_overlap']<=row['threshold']]
            winner=min(admissible,key=lambda r:r['steps']) if admissible else None
            expected=None if winner is None else {'steps':winner['steps'],'rotations':winner['pauli_rotations'],'E':winner['one_minus_overlap']}
            require(row[name]==expected,'Incorrect least-tested threshold aggregation')
        if row['threshold']==1e-6:
            c=row['class_refined_response'];b=row['jw_signed'];factor=b['rotations']/c['rotations']
            require(2<=factor<=5,'Claimed work factor fails')
            threshold_claim.append({'case':row['case'],'candidate_r':c['steps'],'jw_signed_r':b['steps'],'candidate_rotations':c['rotations'],'jw_signed_rotations':b['rotations'],'jw_over_candidate_work':factor})
    outcomes={}
    for r in grid:
        outcomes[str(r)]={}
        for other in names[1:]:
            counts={'wins':0,'ties':0,'losses':0}
            for case in cases:
                candidate=next(x for x in rows if x['case']==case and x['ordering']==names[0] and x['steps']==r)['one_minus_overlap']
                baseline=next(x for x in rows if x['case']==case and x['ordering']==other and x['steps']==r)['one_minus_overlap']
                diff=candidate-baseline;counts['ties' if abs(diff)<=1e-12 else ('wins' if diff<0 else 'losses')]+=1
            outcomes[str(r)][other]=counts
            require(counts=={'wins':8,'ties':0,'losses':0},'All-eight candidate-win claim fails')
    result={'status':'passed','qhat_imports':False,'run_manifest_sha256':sha(a.run/'manifest.json'),'script_sha256':sha(Path(__file__)),
            'new_rows':240,'reused_r100_rows':40,'saved_state_metrics_checked':arraychecks,'max_metric_absolute_difference':maximum_metric_difference,
            'scope':'Independent read-only exact target, schedule and reference provenance, saved-state metric and aggregation audit; no new propagation',
            'cases':casechecks,'threshold_1e_minus6':threshold_claim,'fixed_step_outcomes':outcomes}
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','threshold_1e_minus6','fixed_step_outcomes')}))
if __name__=='__main__':main()
