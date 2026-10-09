#!/usr/bin/env python3
"""Behavior checks: wrong checkout, changed resume protocol, optimized tamper."""
import argparse,json,os,shutil,subprocess,tempfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--qhat-root',type=Path,default=Path('/home/bernalde/repos/qhat'));ap.add_argument('--output',type=Path,default=ROOT/'guard_checks.json');args=ap.parse_args()
PY=sys.executable
ENV=dict(os.environ,PYTHONPATH=str(args.qhat_root.resolve().parent))
DATA=args.qhat_root/'numerical_audit/h2o_recovery/data'
PILOT=args.qhat_root/'numerical_audit/h2o_recovery/validation'
checks=[]
def rejection(name,args,phrase):
 p=subprocess.run([PY,*args],env=ENV,capture_output=True,text=True,timeout=45)
 if p.returncode==0 or phrase not in p.stderr:raise RuntimeError(f'Guard failed {name}: {p.returncode}, {p.stderr}')
 checks.append({'name':name,'return_code':p.returncode,'rejection':phrase})
with tempfile.TemporaryDirectory(prefix='h2o-guards-') as tmp:
 p=Path(tmp)
 rejection('low_energy_rejects_wrong_import_root',[str(ROOT/'run_h2o_low_energy.py'),'--qhat-root','/tmp/not-qhat','--sweep',str(ROOT/'results'),'--output',str(p/'wrong-root')],'Imported unintended QHAT checkout')
 if (p/'wrong-root').exists():raise RuntimeError('Wrong-root check wrote outputs')
 target=p/'bad-protocol';target.mkdir();protocol=json.loads((ROOT/'results/protocol.json').read_text());protocol['T']=999;(target/'protocol.json').write_text(json.dumps(protocol))
 rejection('resume_rejects_changed_protocol',[str(ROOT/'run_h2o_sweep.py'),'--qhat-root',str(args.qhat_root),'--data',str(DATA),'--pilot',str(PILOT),'--output',str(target),'--max-qubits','4'],'Resume protocol/source/input mismatch')
 target=p/'tamper';target.mkdir();shutil.copyfile(ROOT/'results/protocol.json',target/'protocol.json');shutil.copytree(ROOT/'results/source',target/'source')
 first=sorted((ROOT/'results').glob('H2O*/results.json'))[0].parent;shutil.copytree(first,target/first.name)
 states=target/first.name/'states.npz';raw=bytearray(states.read_bytes());raw[-1]^=1;states.write_bytes(raw)
 rejection('python_O_rejects_tampered_states',['-O',str(ROOT/'recheck_h2o_states.py'),'--sweep',str(target),'--low-energy',str(ROOT/'low_energy'),'--data',str(DATA),'--output',str(p/'bad-check.json')],'State archive hash mismatch')
with args.output.open('x') as f:json.dump({'passed':True,'checks':checks},f,indent=2);f.write('\n')
print(json.dumps(checks,indent=2))
