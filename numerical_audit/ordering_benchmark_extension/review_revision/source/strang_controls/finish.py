#!/usr/bin/env python3
"""Continue serial verification/analysis after the production checkpoint exists."""
import resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
import json,subprocess,sys,time
from pathlib import Path
root=Path(__file__).resolve().parent
while not (root/'completion.json').exists():time.sleep(2)
for script,args,log in [('analyze.py',[],'analysis.log'),('independent_verify.py',[],'independent_replay.log'),('audit_records.py',[],'record_audit.log')]:
 with (root/log).open('w') as out:subprocess.run([sys.executable,'-O',str(root/script),*args],cwd=root,stdout=out,stderr=subprocess.STDOUT,check=True)
 print('DONE',script,args,flush=True)
(root/'READY.json').write_text(json.dumps({'status':'complete','production':json.loads((root/'completion.json').read_text()),'verification':json.loads((root/'verification.json').read_text()),'record_audit':json.loads((root/'record_audit.json').read_text())},indent=2)+'\n')
