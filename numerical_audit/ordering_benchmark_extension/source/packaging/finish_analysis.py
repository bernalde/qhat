#!/usr/bin/env python3
"""Continue the independent audit and analysis after the durable V2 campaign ends."""
import json,os,subprocess,sys,time
from pathlib import Path
root=Path(__file__).resolve().parents[1]
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['NUMBA_CACHE_DIR']=str(root/'independent/numba-cache')
completion=root/'published/v2/completion.json'
while not completion.exists():time.sleep(10)
record=json.loads(completion.read_text())
if record.get('rows')!=196:raise RuntimeError('V2 completion has wrong count')
commands=[['-O','independent/check_campaign.py','--campaign','published/v2','--references','independent/references','--output','independent/v2'],['analysis/summarize.py','--root','.','--output','analysis'],['analysis/timings.py','--root','.','--output','analysis'],['analysis/render.py','--analysis','analysis','--output','analysis/figures']]
steps=[]
for cmd in commands:
 start=time.monotonic();r=subprocess.run([sys.executable,*cmd],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);name=Path(cmd[1] if cmd[0]=='-O' else cmd[0]).stem;log=root/'packaging'/('final_'+name+'.log');log.write_text(r.stdout);steps.append({'command':cmd,'returncode':r.returncode,'seconds':time.monotonic()-start});print(name,r.returncode,r.stdout[-2500:],flush=True)
 if r.returncode:raise RuntimeError('Final stage failed: '+str(log))
(root/'packaging/analysis_ready.json').write_text(json.dumps({'status':'complete','steps':steps},indent=2)+'\n')
