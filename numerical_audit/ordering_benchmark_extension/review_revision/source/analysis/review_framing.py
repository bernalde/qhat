#!/usr/bin/env python3
"""Recompute review corrections from frozen raw rows; no propagator imports."""
import argparse, collections, csv, hashlib, json, math, statistics
from pathlib import Path
FLOOR=1e-12
CAND='class_refined_response'; STRANG='strang_class_refined_response'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(v,m):
 if not v: raise ValueError(m)
def compare(a,b):
 ds=[x-y for x,y in zip(a,b)];rr=[x/y for x,y in zip(a,b) if y>FLOOR]
 return dict(wins=sum(d<-FLOOR for d in ds),ties=sum(abs(d)<=FLOOR for d in ds),losses=sum(d>FLOOR for d in ds),median_ratio=statistics.median(rr) if rr else None,ratio_cases=len(rr))
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 inputs={};rows=[]
 for name in ('fixed','optimized','published/extra','published/v2'):
  folder=a.root/name;paths=sorted((folder/'rows').rglob('*.json')) if (folder/'rows').exists() else sorted((folder/'results').glob('*.b*.json'))
  for path in paths:
   r=json.loads(path.read_text());r.update(campaign=name,row_path=str(path.relative_to(a.root)));rows.append(r);inputs[str(path.relative_to(a.root))]=sha(path)
 require(len(rows)==4312,'Frozen campaign row count')
 ix={(r['case'],r['ordering'],int(r['budget_units'])):r for r in rows};require(len(ix)==len(rows),'Duplicate row')
 cases=list(json.loads((a.root/'optimized/ALL_SCHEDULES_FROZEN.json').read_text())['cases']);require(len(cases)==14,'Case count')
 def errs(name,b):return [ix[(c,name,b)]['one_minus_overlap'] for c in cases]
 names=sorted({r['ordering'] for r in rows if r['formula_order']==2});bands={}
 for b in (2,5,10,20,50,100):
  bands[str(b)]={n:compare(errs(CAND,b),errs(n,b)) for n in names}
 work=json.loads((a.root/'optimized/boundary_sensitivity.json').read_text());inputs['optimized/boundary_sensitivity.json']=sha(a.root/'optimized/boundary_sensitivity.json')
 wrows=work['cost_rows'];threshold=[]
 for eps in (1e-4,1e-6,1e-8):
  for c in cases:
   def best(n):
    valid=[r for r in wrows if r['case']==c and r['ordering']==n and r['one_minus_overlap']<=eps]
    return min(valid,key=lambda r:r['uniform_endpoint_rotations']) if valid else None
   s=best(STRANG)
   for n in names:
    if n==STRANG:continue
    z=best(n)
    require(s is not None and z is not None,'Threshold not reached')
    delta=s['original_budget_units']-z['original_budget_units'];status='lower' if delta<0 else 'higher' if delta>0 else 'same'
    threshold.append(dict(case=c,epsilon=eps,comparator=n,conserving_budget=s['original_budget_units'],comparator_budget=z['original_budget_units'],level=status,cost_ratio_comparator_over_conserving=z['uniform_endpoint_rotations']/s['uniform_endpoint_rotations'],conserving_cost=s['uniform_endpoint_rotations'],comparator_cost=z['uniform_endpoint_rotations']))
 totals=dict(collections.Counter(t['level'] for t in threshold));require(len(threshold)==168,'Threshold pairs')
 thrtable=[]
 for eps in (1e-4,1e-6,1e-8):
  for n in names:
   if n!=STRANG:thrtable.append(dict(epsilon=eps,comparator=n,**{s:sum(t['level']==s for t in threshold if t['epsilon']==eps and t['comparator']==n) for s in ('lower','same','higher')}))
 search=[]
 for c in cases:
  path=a.root/'optimized/designs'/f'{c}.design.json';d=json.loads(path.read_text());inputs[str(path.relative_to(a.root))]=sha(path);reports=d['optimizer_reports'];cp=reports['class_pauli_JD_budget'];cu=reports[CAND];u=cp['units']
  # Both starts receive half the total cap; a full insertion sweep costs u(u-1).
  search.append(dict(case=c,conserving_units=cu['units'],pauli_units=u,unit_ratio=u/cu['units'],conserving_moves_per_start=[len(s['history']) for s in cu['starts']],capped_moves_per_start=[len(s['history']) for s in cp['starts']],cap_per_start_as_full_sweep=cp['evaluation_budget']/(2*u*(u-1))))
 j0={str(b):compare(errs(CAND,b),errs('class_refined_J0',b)) for b in (2,5,10,20,50,100)}
 j0['10_per_case']=[dict(case=c,ratio=x/y,candidate=x,unweighted=y) for c,x,y in zip(cases,errs(CAND,10),errs('class_refined_J0',10))]
 rnd=[]
 for family in ('owned','refined','pauli'):
  for c in cases:
   vals=sorted(ix[(c,f'random_{family}_{k:02d}',100)]['one_minus_overlap'] for k in range(20));ce=ix[(c,CAND,100)]['one_minus_overlap']
   def q(t):
    pos=19*t;l=math.floor(pos);h=math.ceil(pos);return vals[l]+(vals[h]-vals[l])*(pos-l)
   rnd.append(dict(case=c,family=family,count=20,min=vals[0],q25=q(.25),median=q(.5),q75=q(.75),max=vals[-1],candidate=ce,strict_wins=sum(ce-x<-FLOOR for x in vals),candidate_over_best=ce/vals[0]))
 result=dict(cases=cases,rows=len(rows),completed=sum(r['status']=='complete' for r in rows),infeasible=dict(collections.Counter(r['campaign'] for r in rows if r['status']=='infeasible')),first_order_vs_strang=bands,threshold_counts=totals,threshold_comparisons=thrtable,threshold_pairs=threshold,search=search,weighted_vs_unweighted=j0,random=rnd,inputs=inputs,source_sha256=sha(Path(__file__)))
 (a.output/'framing_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 for filename,items in [('threshold_levels.csv',threshold),('search_depth.csv',search),('random_distributions.csv',rnd)]:
  with (a.output/filename).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=items[0].keys());w.writeheader();w.writerows(items)
 # Human-facing derivative labels change, never the frozen raw rows.
 with (a.root/'analysis/all_rows.csv').open() as f:reader=csv.DictReader(f);fields=reader.fieldnames;derived=list(reader)
 for r in derived:
  if r['case'] in cases[8:]:r['panel']='followup'
 with (a.output/'all_rows.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(derived)
 text=['# Review framing checks','','All statistics below are recomputed from the frozen production row files; no propagation or new order selection occurs. Original manifests and rows are unchanged.','','Threshold levels (conserving lower / same / higher): '+str(totals),'','| Target | Comparator | Lower / same / higher |','|---|---|---|']
 for t in thrtable:text.append(f"| {t['epsilon']:g} | {t['comparator']} | {t['lower']}/{t['same']}/{t['higher']} |")
 text+=['','The same-level fused rotation counts remain in threshold_levels.csv. They are exact costs of tested circuits, not resolved differences between minimum costs required to reach the target.','',f"First-order candidate beats {sum(t['strict_wins'] for t in rnd)} of 840 fixed random outcomes. Random quartiles are in random_distributions.csv.",'','Per-case search depth, floor-aware formula-order comparisons and unweighted controls are in framing_checks.json.']
 (a.output/'FRAMING_CHECKS.md').write_text('\n'.join(text)+'\n')
 print(json.dumps({k:result[k] for k in ('threshold_counts','infeasible','search')},indent=2))
if __name__=='__main__':main()
