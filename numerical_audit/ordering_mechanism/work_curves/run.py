#!/usr/bin/env python3
"""Fixed-order finite-step controls; not further algorithm development."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import argparse,importlib.util,json,shutil,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('candidate',HERE/'frozen_candidate.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--candidate-run',type=Path,required=True);ap.add_argument('--baseline-run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 sources={s:p.sha(HERE/s) for s in ('run.py','PROTOCOL.md','frozen_candidate.py')}
 for s in sources:shutil.copyfile(HERE/s,a.output/s)
 cases=json.loads((a.candidate_run/'manifest.json').read_text())['cases'];prepared=[];names=('class_refined_response','class_lex','owned_signed','jw_signed','jw_magnitude');rs=(1,2,5,10,20,50)
 for case in cases:
  d=json.loads((a.candidate_run/f'{case}.design.json').read_text());d['orders']={n:d['orders'][n] for n in names if not n.startswith('jw_')};d['jw_orders']={n:d['jw_orders'][n] for n in names if n.startswith('jw_')}
  for field in ('graphs','class_suffix_l1_upper_bounds','suffix_candidate_construction','schedule_sha256','float_schedule_sum_target_l1'):d.pop(field,None)
  d['planned_runs']=[{'ordering':name,'steps':r,'role':'work_curve','T':1.} for name in names for r in rs];d['parent_candidate_design_sha256']=p.sha(a.candidate_run/f'{case}.design.json')
  p.write(a.output/f'{case}.design.json',d);prepared.append((case,d))
 p.write(a.output/'manifest.json',{'cases':cases,'source_sha256':sources,'expected_rows':240,'scope':'fixed-order development-case work curves'})
 rows=[]
 for case,d in prepared:
  den=int(d['denominator']);n=d['n_qubits'];total=p.block(d['coefficients']);blocks={k:[p.block(b) for b in v] for k,v in d['blocks'].items()};keys=sorted(set(total).union(*(blocks['parent'])))
  raw=p.ev.compile_ordered_terms(keys,{k:1. for k in keys},n,0.);basis,reduced,_=p.ref.affine_space(raw,d['affine_basis_indices'][0]);p.require(list(map(int,basis))==d['affine_basis_indices'],'Closure changed');unit=dict(zip(keys,reduced))
  def term(k,v):return(v/den*unit[k][0],unit[k][1],unit[k][2],unit[k][3])
  with np.load(a.candidate_run/f'{case}.states.npz') as saved:initial=saved['initial'];exact=saved['exact'];dense=saved['exact_dense']
  nr=np.array([int(v).bit_count()==d['electrons'] for v in basis]);ss=set(map(int,p.ev.spin_sector_basis_indices(n,d['electrons'])));sr=np.array([int(v) in ss for v in basis]);states={'initial':initial,'exact':exact,'exact_dense':dense}
  for run in d['planned_runs']:
   name,r=run['ordering'],run['steps']
   if name.startswith('jw_'):
    words=[tuple((int(q),s) for q,s in w) for w in d['jw_orders'][name]];schedule=[(w,total[w]) for w in words]
   else:schedule=[(w,v) for i in d['orders'][name] for w,v in blocks[name.split('_')[0]][i].items()]
   p.require(len(schedule)==len(total) and dict(schedule)==total,'Changed rotation multiset');terms=[term(w,v) for w,v in schedule];start=time.monotonic();actual=p.evolve(initial,terms,1.,r)
   row={'case':case,'ordering':name,'T':1.,'steps':r,'role':'work_curve','pauli_rotations':len(terms)*r,'seconds':time.monotonic()-start,**p.physical_metrics(exact,actual,nr,sr)}
   p.write(a.output/f'{case}.{name}.r{r}.work_curve.json',row);rows.append(row);states[f'{name}__r{r}__work_curve']=actual
  np.savez_compressed(a.output/f'{case}.states.npz',**states);print('DONE',case,len(rows),flush=True)
 p.write(a.output/'results.json',rows)
 for source,allowed in ((a.candidate_run,{'class_refined_response'}),(a.baseline_run,set(names)-{'class_refined_response'})):
  for row in json.loads((source/'results.json').read_text()):
   if row['steps']==100 and row['ordering'] in allowed:rows.append({**row,'reused_source_results_sha256':p.sha(source/'results.json')})
 p.write(a.output/'combined_results.json',rows);thresholds=[]
 for case in cases:
  for threshold in (1e-4,1e-6,1e-8):
   row={'case':case,'threshold':threshold,'scope':'least tested grid point; no interpolation'}
   for name in names:
    acceptable=[r for r in rows if r['case']==case and r['ordering']==name and r['one_minus_overlap']<=threshold];best=min(acceptable,key=lambda r:r['steps']) if acceptable else None;row[name]=None if best is None else {'steps':best['steps'],'rotations':best['pauli_rotations'],'E':best['one_minus_overlap']}
   thresholds.append(row)
 p.write(a.output/'thresholds.json',thresholds);p.require(sources=={s:p.sha(HERE/s) for s in sources},'Source changed');p.write(a.output/'completion.json',{'new_rows':240,'combined_rows':len(rows),'status':'complete'})
if __name__=='__main__':main()
