#!/usr/bin/env python3
"""Run two NH3 cases and the smaller cases in three single-thread workers."""
import datetime,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import CASES,sha,write_atomic
jobs={'small':[c for c in CASES if not c.startswith('nh3_')],'nh3_coeff1p00':['nh3_coeff1p00'],'nh3_coeff1p25':['nh3_coeff1p25']}
processes=[]
for label,cases in jobs.items():
 env=dict(os.environ)
 for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
 for key,folder in (('NUMBA_CACHE_DIR','numba'),('MPLCONFIGDIR','mpl')):
  path=HERE/'runtime_cache'/label/folder;path.mkdir(parents=True,exist_ok=True);env[key]=str(path)
 cmd=[sys.executable,str(HERE/'run.py'),'run','--cases',*cases];log=HERE/f'run_{label}.log';handle=log.open('a');proc=subprocess.Popen(cmd,stdout=handle,stderr=subprocess.STDOUT,env=env)
 record={'label':label,'cases':cases,'pid':proc.pid,'command':cmd,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'freeze_sha256':sha(HERE/'freeze.json'),'status':'running','log':log.name,'threads':1,'address_space_limit_GiB':6}
 write_atomic(HERE/f'job_{label}.json',record);processes.append((proc,handle,record))
failed=False
for proc,handle,record in processes:
 code=proc.wait();handle.close();record.update(status='complete' if code==0 else 'failed',returncode=code,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());write_atomic(HERE/f'job_{record["label"]}.json',record);failed|=code!=0
sys.exit(1 if failed else 0)
