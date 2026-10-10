#!/usr/bin/env python3
"""Fresh-directory production and verification smoke test of one full fixture."""
import os,resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import argparse,json,shutil,subprocess,sys,time
from pathlib import Path
import run as r
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--archive-root',type=Path,required=True);p.add_argument('--fresh-directory',type=Path,required=True);a=p.parse_args()
while not (root/'READY.json').exists():time.sleep(2)
r.require(not a.fresh_directory.exists(),'Fresh-directory destination already exists');a.fresh_directory.mkdir(parents=True)
for name in ('run.py','independent_core.py','PROTOCOL.md','independent_verify.py'):shutil.copyfile(root/name,a.fresh_directory/name)
commands=[['run.py','freeze','--archive-root',str(a.archive_root.resolve())],['run.py','run','--archive-root',str(a.archive_root.resolve()),'--cases','CH2O_s1.00'],['independent_verify.py','--cases','CH2O_s1.00']]
records=[]
for ix,command in enumerate(commands):
 log=a.fresh_directory/f'command{ix}.log';before=time.monotonic()
 with log.open('w') as out:subprocess.run([sys.executable,'-O',str(a.fresh_directory/command[0]),*command[1:]],cwd=a.fresh_directory,stdout=out,stderr=subprocess.STDOUT,check=True)
 records.append({'command':command,'seconds':time.monotonic()-before,'log_sha256':r.sha(log)})
import numpy as np
case='CH2O_s1.00';original=json.loads((root/'designs'/f'{case}.json').read_text());fresh=json.loads((a.fresh_directory/'designs'/f'{case}.json').read_text());r.require(original==fresh,'Fresh frozen plan changed');checks=[]
for name in original['methods']:
 for b in r.BUDGETS:
  rel=Path('states')/case/f'{name}.b{b}.npz'
  with np.load(root/rel) as z:x=z['actual'].copy()
  with np.load(a.fresh_directory/rel) as z:y=z['actual'].copy()
  delta=float(np.linalg.norm(x-y));r.require(delta<3e-11,'Portable production mismatch');checks.append({'ordering':name,'budget_units':b,'state_distance':delta,'bitwise_equal':bool(np.array_equal(x,y))})
newchecks=json.loads((a.fresh_directory/'verification'/f'{case}.json').read_text());r.require(newchecks['replays']==36,'Portable independent coverage incomplete')
r.write(root/'portability.json',{'status':'passed','case':case,'fresh_production_states':30,'reused_states':6,'independent_replays':36,'max_state_distance':max(x['state_distance'] for x in checks),'bitwise_equal_states':sum(x['bitwise_equal'] for x in checks),'independent_max_state_distance':newchecks['max_state_distance'],'commands':records,'checks':checks,'source_sha256':r.sha(__file__),'scope':'Fresh-directory execution using copied identical sources and the original immutable archive with exclusively bound QHAT imports; exact pinned existing environment, not a new package installation.'})
print('PORTABILITY PASS: 30 new states + 6 reused; 36 independent replays',flush=True)
