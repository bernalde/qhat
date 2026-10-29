"""Root-owned serialized dispatch of an immutable complete constructor plan."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('plan',type=Path);p.add_argument('--count',type=int,default=1);a=p.parse_args()
raw=a.plan.read_bytes();frozen_hash=hashlib.sha256(raw).hexdigest();plan=json.loads(raw)
assert plan['mode'] in ('preflight','main') and len(plan['entries'])==plan['expected_calls']
dispatcher=Path(__file__).with_name('run_guarded_job.py');remaining=a.count
for job in plan['entries']:
    if (a.root/'jobs'/(job['label']+'.FINAL.json')).is_file():continue
    if remaining<=0:break
    assert hashlib.sha256(a.plan.read_bytes()).hexdigest()==frozen_hash
    cmd=[sys.executable,str(dispatcher),'--root',str(a.root),'--phase','preflight' if plan['mode']=='preflight' else 'constructors','--label',job['label'],'--cpu',str(job['attempt_CPU_seconds']),'--wall',str(int(job['attempt_wall_seconds'])),'--cwd',job['cwd'],'--',*job['command']]
    print('DISPATCH',job['case'],job['arm'],job['repetition'],flush=True)
    status=subprocess.run(cmd).returncode
    final=a.root/'jobs'/(job['label']+'.FINAL.json')
    if not final.is_file():raise RuntimeError('Unreconciled guard failure; no next job')
    print('RETAINED',job['label'],'exit',status,flush=True)
    remaining-=1
