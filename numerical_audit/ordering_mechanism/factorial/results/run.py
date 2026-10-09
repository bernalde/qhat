#!/usr/bin/env python3
"""Complete the frozen 2x2 development comparison; see PROTOCOL.md."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import argparse, importlib.util, json, shutil, time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def load(name):
 spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
p=load('refined_source');old=load('complete_source')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 sources={name:p.sha(HERE/name) for name in ('run.py','PROTOCOL.md','refined_source.py','complete_source.py')}
 for name in sources:shutil.copyfile(HERE/name,a.output/name)
 origcases=json.loads((a.parent_run/'manifest.json').read_text())['cases'];prepared=[];cases=[]
 for orig in origcases:
  for refined in (False,True):
   arm='refined_adjacent' if refined else 'complete_insertion';case=orig+'__'+arm;cases.append(case)
   d=json.loads((a.parent_run/f'{orig}.design.json').read_text());den=int(d['denominator']);bs=[p.block(b) for b in d['blocks']['class']];xs=d['class_flip_masks'];hf=(1<<d['electrons'])-1
   if refined:
    newbs=[];newxs=[]
    for b,x in zip(bs,xs):
     for unit in ([{w:v} for w,v in sorted(b.items())] if x==0 else [b]):newbs.append(unit);newxs.append(x)
    bs,xs=newbs,newxs
   pairs,rho=p.pair_setup(bs,xs,hf,den);lex=list(range(len(bs)));mag=sorted(lex,key=lambda i:(-sum(abs(v) for v in bs[i].values()),xs[i]))
   optimizer=old.optimize if refined else p.optimize;start=time.monotonic();candidates=[optimizer(s,pairs,rho,den) for s in (lex,mag)];elapsed=time.monotonic()-start;winner=candidates[0]
   for c in candidates[1:]:
    tol=max(1e-30,1e-12*min(winner[1],c[1]))
    if c[1]<winner[1]-tol or (abs(c[1]-winner[1])<=tol and tuple(c[0])<tuple(winner[0])):winner=c
   order,J,history=winner;name='class_'+arm
   d['case']=case;d['source_case']=orig;d['source_design_sha256']=p.sha(a.parent_run/f'{orig}.design.json')
   d['blocks']['class']=[[[word,str(value)] for word,value in sorted(b.items())] for b in bs];d['class_flip_masks']=xs;d['orders']={name:order};d['jw_orders']={}
   for field in ('graphs','class_suffix_l1_upper_bounds','suffix_candidate_construction','schedule_sha256'):d.pop(field,None)
   d['unit_partition']='split diagonal words' if refined else 'complete classes'
   d['planned_runs']=[{'ordering':name,'steps':100,'role':'factorial','T':1.}]
   p.write(a.output/f'{case}.design.json',d);p.write(a.output/f'{case}.candidate.json',{'source_case':orig,'arm':arm,'J_D':J,'order':order,'moves':history,'optimizer_seconds':elapsed,'source_sha256':sources})
   prepared.append((case,orig,d,bs,J))
 p.write(a.output/'manifest.json',{'cases':cases,'source_sha256':sources,'expected_rows':16,'scope':'post-development factorial control'})
 rows=[]
 for case,orig,d,bs,J in prepared:
  den=int(d['denominator']);n=d['n_qubits'];total=p.block(d['coefficients']);keys=sorted(set(total).union(*(p.block(b) for b in d['blocks']['parent'])))
  raw=p.ev.compile_ordered_terms(keys,{k:1. for k in keys},n,0.);basis,reduced,_=p.ref.affine_space(raw,d['affine_basis_indices'][0]);p.require(list(map(int,basis))==d['affine_basis_indices'],'Closure changed');unit=dict(zip(keys,reduced))
  def term(k,v):return (v/den*unit[k][0],unit[k][1],unit[k][2],unit[k][3])
  with np.load(a.parent_run/f'{orig}.states.npz') as saved:initial=saved['initial'];exact=saved['exact'];dense=saved['exact_dense']
  nr=np.array([int(v).bit_count()==d['electrons'] for v in basis]);ss=set(map(int,p.ev.spin_sector_basis_indices(n,d['electrons'])));sr=np.array([int(v) in ss for v in basis])
  run=d['planned_runs'][0];name=run['ordering'];terms=[term(k,v) for i in d['orders'][name] for k,v in bs[i].items()]
  p.require(len(terms)==len(total),'Rotation count changed');start=time.monotonic();actual=p.evolve(initial,terms,1.,100)
  row={'case':case,'source_case':orig,'ordering':name,'T':1.,'steps':100,'role':'factorial','pauli_rotations':len(terms)*100,'seconds':time.monotonic()-start,'surrogate_J':J,**p.physical_metrics(exact,actual,nr,sr)}
  p.write(a.output/f'{case}.{name}.r100.factorial.json',row);rows.append(row);np.savez_compressed(a.output/f'{case}.states.npz',initial=initial,exact=exact,exact_dense=dense,**{f'{name}__r100__factorial':actual})
  print(case,row['one_minus_overlap'],flush=True)
 p.write(a.output/'results.json',rows);p.require(sources=={name:p.sha(HERE/name) for name in sources},'Source changed');p.write(a.output/'completion.json',{'rows':len(rows),'status':'complete'})
if __name__=='__main__':main()
