#!/usr/bin/env python3
"""Exercise archive integrity rejection with assertions disabled."""
import argparse,hashlib,json,subprocess,sys,tempfile,zipfile
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();verifier=Path(__file__).with_name('verify_archive.py');evidence=[]
 with tempfile.TemporaryDirectory(prefix='ordering-archive-tests-') as tmp:
  base=Path(tmp);data=b'original scientific evidence\n';manifest={'files':{'extension/file.txt':{'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}}}
  for kind in ['valid','changed_member','missing_member','unlisted_member','duplicate_member','unsafe_path','wrong_external_hash']:
   path=base/(kind+'.zip');m=json.loads(json.dumps(manifest));name='extension/file.txt'
   if kind=='unsafe_path':m['files']['../escape']=m['files'].pop(name);name='../escape'
   with zipfile.ZipFile(path,'w') as z:
    if kind!='missing_member':z.writestr(name,b'corrupted scientific evidence\n' if kind=='changed_member' else data)
    if kind=='duplicate_member':z.writestr(name,data)
    if kind=='unlisted_member':z.writestr('extension/unrecorded',b'new')
    z.writestr('artifact_manifest.json',json.dumps(m))
   cmd=[sys.executable,'-O',str(verifier),str(path)]
   if kind=='wrong_external_hash':cmd+=['--expected-sha256','0'*64]
   r=subprocess.run(cmd,text=True,capture_output=True);expected=kind=='valid'
   if (r.returncode==0)!=expected:raise RuntimeError('Verifier regression '+kind+': '+r.stdout+r.stderr)
   evidence.append({'probe':kind,'expected_acceptance':expected,'returncode':r.returncode,'message':r.stdout.strip() if expected else r.stderr.strip().splitlines()[-1]})
 result={'status':'passed','python':sys.version,'verifier_sha256':hashlib.sha256(verifier.read_bytes()).hexdigest(),'probes':evidence};a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','probes':len(evidence)}))
if __name__=='__main__':main()
