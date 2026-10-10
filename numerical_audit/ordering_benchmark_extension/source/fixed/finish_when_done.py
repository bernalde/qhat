#!/usr/bin/env python3
"""Collect durable results after all capped worker job records finish."""
import json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import *
while True:
 jobs=[json.loads(p.read_text()) for p in sorted(HERE.glob('job_*.json'))]
 require(len(jobs)==3,'Expected three worker records')
 require(all(j['status']!='failed' for j in jobs),'A benchmark worker failed; inspect its durable log')
 if all(j['status']=='complete' for j in jobs):break
 time.sleep(10)
subprocess.run([sys.executable,str(HERE/'run.py'),'collect'],check=True)
frozen=json.loads((HERE/'freeze.json').read_text())
for name,h in frozen['source_sha256'].items():require(sha(HERE.parent/name)==h,'Frozen source changed '+name)
for case,h in frozen['design_sha256'].items():require(sha(HERE/'designs'/f'{case}.json')==h,'Frozen design changed '+case)
provenance=json.loads((HERE/'execution_provenance.json').read_text())
for name,item in provenance['qhat_sources'].items():require(sha(QHAT/name)==item['sha256'],'Propagation dependency changed '+name)
write_atomic(HERE/'final_execution_check.json',{'status':'pass','freeze_sha256':sha(HERE/'freeze.json'),'source_files_unchanged':len(frozen['source_sha256']),'designs_unchanged':len(frozen['design_sha256']),'qhat_dependencies_unchanged':len(provenance['qhat_sources']),'completed_workers':3,'results_sha256':sha(HERE/'results.json')})
print((HERE/'completion.json').read_text(),flush=True)
