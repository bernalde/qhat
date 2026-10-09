#!/usr/bin/env python3
"""Exploratory exact signed BCH-l1 adjacent swaps; separate ten-row supplement."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-class-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-class-numba')
import argparse
from collections import Counter
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics,sha,write,object_sha
from exact_blocks import commutator,add_into

HERE=Path(__file__).resolve().parent


def block(raw):return {tuple((int(q),str(p)) for q,p in key):int(value) for key,value in raw}
def records(b):return [[key,str(value)] for key,value in sorted(b.items())]


def improve(blocks):
    order=list(range(len(blocks)));cache={};initial={}
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            c=commutator(blocks[i],blocks[j])
            if c:cache[i,j]=c;add_into(initial,c)
    current=dict(initial);objective=sum(map(abs,current.values()));accepted=[]
    for iteration in range(2):
        for at in range(len(order)-1):
            i,j=order[at:at+2];c=cache.get((min(i,j),max(i,j)),{})
            factor=-2 if i<j else 2
            candidate=objective+sum(abs(current.get(k,0)+factor*v)-abs(current.get(k,0)) for k,v in c.items())
            if candidate<objective:
                add_into(current,c,factor);order[at:at+2]=[j,i]
                accepted.append({'pass':iteration,'position':at,'old_pair':[i,j],
                    'old_objective':str(objective),'new_objective':str(candidate)})
                objective=candidate
    # Independently reconstruct the full signed sum from the resulting order.
    final={}
    for at,i in enumerate(order):
        for j in order[at+1:]:
            add_into(final,cache.get((min(i,j),max(i,j)),{}),1 if i<j else -1)
    if final!=current or objective!=sum(map(abs,final.values())):raise RuntimeError('Adjacent BCH update does not match full recomputation.')
    return order,initial,final,accepted


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));a.output.mkdir(parents=True,exist_ok=False)
    parent_manifest=json.loads((a.parent_run/'manifest.json').read_text());cases=parent_manifest['cases']
    if len(cases)!=8:raise ValueError('Supplement requires original eight-case campaign.')
    sources=[Path(__file__),HERE/'exact_blocks.py',HERE/'bch_supplement_protocol.md',HERE/'recheck_saved.py']
    source_hashes={str(s):sha(s) for s in sources}
    for source in sources:
        target=a.output/'source_snapshot'/source.name;target.parent.mkdir(exist_ok=True);shutil.copy2(source,target)
    write(a.output/'manifest.json',{'version':'exploratory-class-bch-l1-v1','cases':cases,'expected_rows':10,
        'parent_manifest_sha256':sha(a.parent_run/'manifest.json'),'source_sha256':source_hashes,
        'scope':'Exploratory supplement after suffix-l1 diagnosis; no tuning to errors; separate from original106rows',
        'initial_order':'class_lex','adjacent_passes':2,'acceptance':'strict exact integer globalBCH-l1 decrease',
        'methods':['class_bch_l1'],'T':1.,'steps':100,'anchor_case':'HCN_s1.00','anchor_steps':[50,200]})
    rows=[];analyses=[];began=time.monotonic()
    oldrows=json.loads((a.parent_run/'results.json').read_text())
    for case in cases:
        design_path=a.parent_run/f'{case}.design.json';d=json.loads(design_path.read_text())
        den=int(d['denominator']);blocks=[block(b) for b in d['blocks']['class']]
        order,initialB,finalB,accepted=improve(blocks)
        old_norm=sum(map(abs,initialB.values()));new_norm=sum(map(abs,finalB.values()))
        planned=[{'ordering':'class_bch_l1','T':1.,'steps':r,'role':'exploratory_primary' if r==100 else 'exploratory_convergence'}
                 for r in ([100,50,200] if case=='HCN_s1.00' else [100])]
        supplement={'case':case,'parent_design_sha256':sha(design_path),'parent_states_sha256':sha(a.parent_run/f'{case}.states.npz'),
            'initial_B_coefficients':records(initialB),'final_B_coefficients':records(finalB),
            'initial_l1_integer':str(old_norm),'final_l1_integer':str(new_norm),
            'initial_l1':old_norm/(den*den),'final_l1':new_norm/(den*den),
            'final_over_initial':new_norm/old_norm if old_norm else None,
            'accepted_swaps':accepted,'order':order,'planned_runs':planned}
        # Complete coefficients, schedule and objective are saved before state evaluation.
        write(a.output/f'{case}.bch_design.json',supplement)
        d['orders']['class_bch_l1']=order;d['planned_runs']=planned
        d['source_parent_design_sha256']=sha(design_path)
        d['schedule_sha256']['class_bch_l1']=object_sha([[k,str(v)] for i in order for k,v in blocks[i].items()])
        write(a.output/f'{case}.design.json',d)
        n=d['n_qubits'];keys=sorted(set().union(*(set(block(b)) for b in d['blocks']['parent'])))
        fullunit=evolution.compile_ordered_terms(keys,{k:1. for k in keys},n,0.)
        states,reduced,_=ref.affine_space(fullunit,d['affine_basis_indices'][0])
        if list(map(int,states))!=d['affine_basis_indices']:raise RuntimeError('Closure changed.')
        unit=dict(zip(keys,reduced))
        terms=[(v/den*unit[k][0],unit[k][1],unit[k][2],unit[k][3]) for i in order for k,v in blocks[i].items()]
        number=np.array([int(x).bit_count()==d['electrons'] for x in states])
        physical=set(map(int,evolution.spin_sector_basis_indices(n,d['electrons'])))
        spin=np.array([int(x) in physical for x in states])
        with np.load(a.parent_run/f'{case}.states.npz') as saved:
            initial=saved['initial'];exact=saved['exact'];dense=saved['exact_dense']
        state_record={'initial':initial,'exact':exact,'exact_dense':dense}
        for run in planned:
            r=run['steps'];role=run['role'];start=time.monotonic();approx=evolve(initial,terms,1.,r)
            row={'case':case,'ordering':'class_bch_l1','T':1.,'steps':r,'role':role,
                 'target_sha256':d['target_sha256'],'schedule_sha256':d['schedule_sha256']['class_bch_l1'],
                 'pauli_rotations':r*len(terms),'seconds':time.monotonic()-start,
                 'bch_l1_initial':supplement['initial_l1'],'bch_l1_final':supplement['final_l1'],
                 **physical_metrics(exact,approx,number,spin)}
            old=next(x for x in oldrows if x['case']==case and x['ordering']=='class_lex' and x['steps']==r)
            row['class_lex_E']=old['one_minus_overlap'];row['error_ratio_over_class_lex']=row['one_minus_overlap']/old['one_minus_overlap']
            rows.append(row);write(a.output/f'{case}.class_bch_l1.r{r}.{role}.json',row)
            state_record[f'class_bch_l1__r{r}__{role}']=approx
            print('BCH',len(rows),case,r,'swaps',len(accepted),'normratio',supplement['final_over_initial'],'Eratio',row['error_ratio_over_class_lex'],flush=True)
        np.savez_compressed(a.output/f'{case}.states.npz',**state_record);analyses.append(supplement)
    if len(rows)!=10:raise RuntimeError('Incomplete supplement.')
    primary=[r for r in rows if r['role']=='exploratory_primary']
    outcomes=Counter('tie' if abs(r['one_minus_overlap']-r['class_lex_E'])<=1e-12 else ('win' if r['one_minus_overlap']<r['class_lex_E'] else 'loss') for r in primary)
    summary={'rows':len(rows),'primary_rows':len(primary),'win_tie_loss_vs_class_lex':dict(outcomes),
        'median_error_ratio_over_class_lex':float(np.median([r['error_ratio_over_class_lex'] for r in primary])),
        'bch_l1_ratio_range':[min(x['final_over_initial'] for x in analyses),max(x['final_over_initial'] for x in analyses)],
        'max_number_leakage':max(r['particle_number_leakage'] for r in rows),
        'interpretation':'Exploratory fixed-rule surrogate comparison; decreasingBCH-l1doesnotguaranteedecreasingfinite-timeoverlaperror.'}
    write(a.output/'results.json',rows);write(a.output/'summary.json',summary)
    lines=['# Exploratory signed BCH-l1 adjacent-swap supplement','',
        'This ten-row supplement followed the suffix-bound diagnosis. It is separate from the original 106-row campaign. All schedules were frozen before their errors were evaluated; two strict-improvement adjacent passes started from class_lex.','',
        '| Case | Accepted swaps | BCH-l1 new/old | E new/lex at T=1,r=100 |','|---|---:|---:|---:|']
    for d in analyses:
        row=next(r for r in primary if r['case']==d['case'])
        lines.append(f'| {d["case"]} | {len(d["accepted_swaps"])} | {d["final_over_initial"]:.6g} | {row["error_ratio_over_class_lex"]:.6g} |')
    lines+=['','The exact coefficient-l1 BCH surrogate decreases in every accepted swap. This is not a finite-time accuracy guarantee or a transported-variance calculation. The outcome is descriptive for the fixed HF state on four molecular families.','']
    (a.output/'REPORT.md').write_text('\n'.join(lines))
    for s,digest in source_hashes.items():
        if sha(s)!=digest:raise RuntimeError('Source changed during supplement.')
    write(a.output/'completion.json',{'status':'complete','rows':len(rows),'seconds':time.monotonic()-began})

if __name__=='__main__':main()
