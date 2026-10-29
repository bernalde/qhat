"""Root-only serial execution of a hash-stable bounded stage plan."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('plan',type=Path);a=p.parse_args();raw=a.plan.read_bytes();digest=hashlib.sha256(raw).hexdigest();plan=json.loads(raw)
root=Path(plan['root']);dispatcher=Path(__file__).with_name('run_guarded_job.py')
for job in plan['jobs']:
 assert hashlib.sha256(a.plan.read_bytes()).hexdigest()==digest
 final=root/'jobs'/(job['label']+'.FINAL.json')
 if final.exists():
  assert json.loads(final.read_text())['exit_code']==0,'Retained failure prevents dependent work'
  continue
 print('STAGE',job['label'],flush=True)
 cmd=[sys.executable,str(dispatcher),'--root',str(root),'--phase',job['phase'],'--label',job['label'],'--cpu',str(job['cpu']),'--wall',str(job['wall']),'--cwd',plan['cwd'],'--',*job['command']]
 status=subprocess.run(cmd).returncode
 if status or not final.is_file():raise SystemExit(status or 2)
 print('STAGE_COMPLETE',job['label'],flush=True)
