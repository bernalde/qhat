#!/usr/bin/env python3
"""Freeze optimized schedules, then perform equal-Pauli-rotation-budget trials."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,gc,json,sys,time,resource
HERE=Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR',str(HERE/'cache'/'numba'))
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'cache'/'matplotlib'))
sys.path.insert(0,str(HERE.parent))
import numpy as np
import common as c
import optimize as opt
METHODS=('class_refined_response','class_refined_J0','class_pauli_JD','class_pauli_JD_budget')
STRANG=('strang_class_refined_response','strang_jw_magnitude','strang_jw_signed','strang_class_lex')

def encoded_block(b):return c.encode_schedule(list(b.items()))
def source_hashes():return {str(p.relative_to(HERE.parent)):c.sha(p) for p in (HERE/'run.py',HERE/'optimize.py',HERE/'frozen_candidate.py',HERE/'PROTOCOL.md',HERE.parent/'common.py')}

def optimize_two(case,name,bs,xs,pairs,rho,den,budget=None):
    lex=list(range(len(bs)));mag=sorted(lex,key=lambda i:(-sum(abs(v) for v in bs[i].values()),xs[i]))
    cap=None if budget is None else budget//2
    starts=[]
    for index,start in enumerate((lex,mag)):
        limit=cap if index==0 or budget is None else budget-cap
        starts.append(opt.search(start,pairs,rho,den,max_evals=limit,checkpoint=HERE/'checkpoints'/f'{case}.{name}.start{index}.json'))
    chosen=opt.choose(starts)
    return {'order':chosen['order'],'J':chosen['J'],'accepted_moves':chosen['history'],'starts':starts,'evaluations':sum(s['evaluations'] for s in starts),'rescored':sum(s['rescored'] for s in starts),'search_seconds':sum(s['seconds'] for s in starts),'evaluation_budget':budget,'nonzero_pair_count':len(pairs),'units':len(bs),'rho_sha256':c.object_sha(rho),'unit_partition': 'individual Pauli words' if 'pauli' in name else 'individual diagonal Paulis and complete nonzero flip classes'}

def prepare(case):
    output=HERE/'designs'/f'{case}.design.json'
    if output.exists():
        d=json.loads(output.read_text());c.require(d['source_sha256']==source_hashes(),'Preparation source changed after freeze');print('SKIP PREPARED',case,flush=True);return
    begin=time.monotonic();n=c.normalized_case(case);raw=n['design'];den=n['denominator'];hf=(1<<raw['electrons'])-1;bs=n['blocks']['refined'];xs=n['refined_flip_masks']
    before=time.monotonic();pairs,rho=opt.pair_setup(bs,xs,hf,den);setup_refined=time.monotonic()-before
    reports={}
    for name,weights in (('class_refined_response',rho),('class_refined_J0',{k:1. for k in rho})):
        reports[name]=optimize_two(case,name,bs,xs,pairs,weights,den)
        reports[name]['pair_setup_seconds']=setup_refined
    c.require(reports['class_refined_response']['order']==n['candidate_order'],'Stable refined search differs from archived frozen candidate')
    reference_budget=reports['class_refined_response']['evaluations']
    del pairs,rho;gc.collect()
    pbs=[{w:n['total'][w]} for w in sorted(n['total'],key=lambda w:(c.masks(w)[0],w))];px=[c.masks(next(iter(b)))[0] for b in pbs]
    before=time.monotonic();pairs,rho=opt.pair_setup(pbs,px,hf,den,pauli=True);setup_pauli=time.monotonic()-before
    for name,budget in (('class_pauli_JD',None),('class_pauli_JD_budget',reference_budget)):
        reports[name]=optimize_two(case,name,pbs,px,pairs,rho,den,budget)
        reports[name]['pair_setup_seconds']=setup_pauli
    if case in c.DEV_CASES:
        prev=c.MECH/'objective_ablation/results/class_pauli_JD'/f'{case}.design.json';old=json.loads(prev.read_text())
        c.require(reports['class_pauli_JD']['order']==old['orders']['class_pauli_JD'],'Stable Pauli search differs from archived development order')
        prev=c.MECH/'objective_ablation/results/class_refined_J0'/f'{case}.design.json';old=json.loads(prev.read_text())
        c.require(reports['class_refined_J0']['order']==old['orders']['class_refined_J0'],'Stable J0 search differs from archived development order')
    del pairs,rho;gc.collect()
    methods={}
    def add(name,units,order,formula,objective=None):
        scheduled=[(w,v) for i in order for w,v in units[i].items()];encoded=c.encode_schedule(scheduled)
        c.require(dict(scheduled)==n['total'] and len(scheduled)==len(n['total']),'Method changes exact Pauli target')
        methods[name]={'schedule':encoded,'schedule_sha256':c.object_sha(encoded),'unit_order':order,'unit_lengths':[len(units[i]) for i in order],'formula_order':formula,'objective':objective,'terms':len(scheduled)}
    for name in METHODS:add(name,pbs if 'pauli' in name else bs,reports[name]['order'],1,reports[name]['J'])
    add('strang_class_refined_response',bs,reports['class_refined_response']['order'],2,reports['class_refined_response']['J'])
    for name in ('jw_magnitude','jw_signed'):
        units=[{w:n['total'][w]} for w in n['jw_orders'][name]];add('strang_'+name,units,list(range(len(units))),2)
    cb=n['blocks']['class'];add('strang_class_lex',cb,list(range(len(cb))),2)
    d={k:raw[k] for k in ('case','n_qubits','electrons','denominator','coefficients','affine_basis_indices','tensor_sha256','state')}
    d.update(blocks={'parent':raw['blocks']['parent'],'refined':[encoded_block(b) for b in bs],'pauli':[encoded_block(b) for b in pbs],'class':raw['blocks']['class']},extension_input_provenance=raw['extension_input_provenance'],panel='development' if case in c.DEV_CASES else 'prospective',methods=methods,optimizer_reports=reports,source_sha256=source_hashes(),preparation_seconds=time.monotonic()-begin,product_convention='Schedules list chronological factors; later factors multiply on the left.',T=1.)
    c.write_atomic(output,d);print('FROZEN',case,'seconds',d['preparation_seconds'],'eval_budget',reference_budget,flush=True)

def freeze():
    paths=[HERE/'designs'/f'{case}.design.json' for case in c.CASES]
    c.require(all(p.exists() for p in paths),'Some schedules have not been frozen')
    value={'source_sha256':source_hashes(),'design_sha256':{case:c.sha(p) for case,p in zip(c.CASES,paths)},'cases':c.CASES,'expected_rows':728,'T':1.,'budget_units':list(c.BUDGETS),'strang_budget_units':list(c.BUDGETS[1:]),'cost':'Pauli rotations; commuting boundary blocks fused exactly; steps are largest feasible within each budget','formula_order_method_names':{'1':list(METHODS),'2':list(STRANG)}}
    out=HERE/'ALL_SCHEDULES_FROZEN.json'
    if out.exists():c.require(json.loads(out.read_text())==value,'Freeze manifest changed')
    else:c.write_atomic(out,value)

def strang_pulses(lengths,r):
    U=sum(lengths)
    if len(lengths)==1:return np.arange(U,dtype=np.int32),np.full(U,2*r,dtype=np.int32)
    offsets=np.cumsum([0]+lengths);blocks=[list(range(int(offsets[i]),int(offsets[i+1]))) for i in range(len(lengths))];ids=[];dur=[]
    def append(unit,mult):ids.extend(unit);dur.extend([mult]*len(unit))
    append(blocks[0],1)
    for step in range(r):
        for b in blocks[1:-1]:append(b,1)
        append(blocks[-1],2)
        for b in reversed(blocks[1:-1]):append(b,1)
        append(blocks[0],1 if step==r-1 else 2)
    ii=np.asarray(ids,dtype=np.int32);dd=np.asarray(dur,dtype=np.int32)
    c.require(len(ii)==r*(2*U-lengths[0]-lengths[-1])+lengths[0],'Strang pulse count mismatch')
    totals=np.bincount(ii,weights=dd,minlength=U);c.require(np.all(totals==2*r),'Strang per-word duration mismatch')
    return ii,dd

def run_case(case):
    manifest=json.loads((HERE/'ALL_SCHEDULES_FROZEN.json').read_text());path=HERE/'designs'/f'{case}.design.json';d=json.loads(path.read_text())
    c.require(manifest['source_sha256']==source_hashes() and manifest['design_sha256'][case]==c.sha(path),'Frozen sources/design changed')
    raw,states=c.load_case(case);compiled=c.compile_case(raw,states)
    from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    for name,method in d['methods'].items():
        schedule=c.decode_schedule(method['schedule']);U=len(schedule);terms=[compiled['term'](w,v) for w,v in schedule]
        for b in (c.BUDGETS if method['formula_order']==1 else c.BUDGETS[1:]):
            rowpath=HERE/'rows'/case/f'{name}.b{b}.json';statepath=HERE/'states'/case/f'{name}.b{b}.npz'
            if rowpath.exists():
                row=json.loads(rowpath.read_text());c.require(statepath.exists() and row['state_sha256']==c.sha(statepath),'Incomplete existing row');continue
            B=b*U;lengths=method['unit_lengths'];r=b if method['formula_order']==1 else (B-lengths[0])//(2*U-lengths[0]-lengths[-1])
            c.require(r>=1,'Unexpected infeasible scheduled Strang budget')
            before=time.monotonic();pulses={}
            if method['formula_order']==1:actual=evolve(compiled['initial'],terms,1.,r);rotations=r*U
            else:
                ii,dd=strang_pulses(lengths,r);actual=compiled['initial'].copy();work=np.empty_like(actual)
                for i,mult in zip(ii,dd):actual,work=ev.apply_one_exponential(actual,work,terms[int(i)],int(mult)/(2*r),2**30)
                rotations=len(ii);pulses={'pulse_term_indices':ii,'pulse_duration_numerators':dd,'pulse_duration_denominator':np.array(2*r)}
            c.require(rotations<=B,'Exceeded Pauli-rotation budget')
            statepath.parent.mkdir(parents=True,exist_ok=True);temp=statepath.with_suffix('.partial.npz');np.savez_compressed(temp,actual=actual,**pulses);os.replace(temp,statepath)
            row={'case':case,'panel':d['panel'],'ordering':name,'T':1.,'formula_order':method['formula_order'],'budget_units':b,'rotation_budget':B,'steps':r,'role':'matched_rotation_budget','extension_design_sha256':c.sha(path),'schedule_sha256':method['schedule_sha256'],'status':'complete','pauli_rotations':rotations,'seconds':time.monotonic()-before,'state_file':str(statepath.relative_to(HERE)),'state_key':'actual','state_sha256':c.sha(statepath),'surrogate_J':method['objective'],**physical_metrics(compiled['exact'],actual,compiled['number_mask'],compiled['spin_mask'])}
            c.write_atomic(rowpath,row);print('EVOLVE',case,name,b,r,rotations,row['one_minus_overlap'],'seconds',row['seconds'],flush=True)

def main():
    a=argparse.ArgumentParser();a.add_argument('command',choices=['prepare','freeze','run','finish']);a.add_argument('--cases',nargs='*',default=c.CASES);args=a.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    if args.command=='prepare':
        for case in args.cases:prepare(case)
    elif args.command=='freeze':freeze()
    elif args.command=='run':
        for case in args.cases:run_case(case)
    else:
        rows=[json.loads(p.read_text()) for p in sorted((HERE/'rows').rglob('*.json'))];c.require(len(rows)==728,'Wrong final row count');c.write_atomic(HERE/'results.json',rows);c.write_atomic(HERE/'completion.json',{'rows':len(rows),'designs':14,'source_sha256':source_hashes(),'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
if __name__=='__main__':main()
