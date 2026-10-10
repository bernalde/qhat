#!/usr/bin/env python3
"""Supplementary provenance and reporting audit, with active gates under -O."""
import os,resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import json
from pathlib import Path
import run as r
from analyze import masks
ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results.json').read_text());r.require(len(rows)==504,'Incomplete rows');freeze=json.loads((ROOT/'ALL_PLANS_FROZEN.json').read_text());seen=set();maximum_norm_report_difference=0;commuting_units=0
for case in r.CASES:
 d=r.verified_design(case);design_sha=r.sha(ROOT/'designs'/f'{case}.json')
 for name,method in d['methods'].items():
  r.require(r.objsha(method['schedule'])==method['schedule_sha256'],'Bad method digest');offset=0
  for length in method['unit_lengths']:
   unit=method['schedule'][offset:offset+length];offset+=length;mm=[masks(word) for word,v in unit]
   for i,(x,z) in enumerate(mm):
    for xx,zz in mm[:i]:r.require(((x&zz)^(z&xx)).bit_count()%2==0,'Noncommuting declared unit')
   commuting_units+=1
 for row in [v for v in rows if v['case']==case]:
  name=row['ordering'];b=row['budget_units'];key=(case,name,b);r.require(key not in seen,'Duplicate row');seen.add(key);method=d['methods'][name];rp=ROOT/'rows'/case/f'{name}.b{b}.json'
  r.require(row==json.loads(rp.read_text()),'Aggregate row drift');r.require(row['extension_design_sha256']==design_sha,'Wrong row design');r.require(row['schedule_sha256']==method['schedule_sha256'],'Wrong row schedule');r.require(row['T']==1 and row['formula_order']==2 and row['status']=='complete','Unexpected row scope');r.require(row['state_sha256']==r.sha(ROOT/row['state_file']),'State digest mismatch')
  check=json.loads((ROOT/'verification'/case/f'{name}.b{b}.json').read_text());r.require(check['row_sha256']==r.sha(rp) and check['state_sha256']==row['state_sha256'] and check['source_sha256']==r.sources() and check['verification_source_sha256']==r.sha(ROOT/'independent_verify.py'),'Stale independent verification')
  diff=abs(check['independent_metrics']['norm']-row['approx_norm']);maximum_norm_report_difference=max(maximum_norm_report_difference,diff);r.require(diff<1e-12,'Reported norm differs from independent norm')
  r.require(row['reused']==method['reused'],'Reuse status mismatch')
  if row['reused']:r.require(row['source_original_state_sha256']==row['state_sha256'],'Reused baseline changed')
r.require(len(seen)==504,'Missing unique rows')
r.write(ROOT/'record_audit.json',{'status':'passed','rows':504,'declared_units_checked_pairwise_commuting':commuting_units,'maximum_norm_report_difference':maximum_norm_report_difference,'source_sha256':r.sha(__file__),'freeze_sha256':r.sha(ROOT/'ALL_PLANS_FROZEN.json'),'results_sha256':r.sha(ROOT/'results.json'),'verification_sha256':r.sha(ROOT/'verification.json'),'scope':'Every unique row matches frozen design and schedule, saved state and independent verification hashes; exact pairwise commutation checked inside every declared unit; reported norms checked independently; reused states unchanged.'})
print('PASS 504 record chains and declared-unit commutation')
