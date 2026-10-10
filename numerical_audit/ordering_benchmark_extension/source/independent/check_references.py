#!/usr/bin/env python3
"""Independent target reference and replay smoke checks on every frozen input."""
import argparse,io,json,os,time,zipfile,resource
from pathlib import Path
from core import OccupationAction,require,sha,block,key
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--qhat-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));a.output.mkdir(parents=True,exist_ok=True)
    dev=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]
    pros=[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
    base=a.qhat_root/'numerical_audit/followup_classes/results/comparison';archive=a.qhat_root/'numerical_audit/ordering_mechanism/heldout/heldout_validation.zip'
    records=[]
    for case in dev+pros:
        output=a.output/f'{case}.json'
        if output.exists():
            old=json.loads(output.read_text());require(old['source_sha256']==sha(__file__) and old['core_sha256']==sha(Path(__file__).with_name('core.py')),'Verifier source changed; choose new output');records.append(old);continue
        start=time.monotonic()
        if case in dev:
            dp=base/f'{case}.design.json';sp=base/f'{case}.states.npz';db=dp.read_bytes();sb=sp.read_bytes()
        else:
            with zipfile.ZipFile(archive) as z:db=z.read(f'results/{case}.design.json');sb=z.read(f'results/{case}.states.npz')
        d=json.loads(db)
        with np.load(io.BytesIO(sb)) as z:arrays={k:z[k].copy() for k in z.files}
        action=OccupationAction(d);reference,nnz=action.reference(arrays['initial']);distance=float(np.linalg.norm(reference-arrays['exact']))
        require(distance<2e-11,f'{case} reference mismatch: {distance}')
        replay=None
        if case in dev:
            total=block(d['coefficients']);schedule=[(w,str(total[key(w)])) for w in d['jw_orders']['jw_magnitude']]
            actual=action.replay(arrays['initial'],schedule,100);replay=float(np.linalg.norm(actual-arrays['jw_magnitude__r100__equal_steps']))
            require(replay<2e-11,f'{case} replay smoke mismatch: {replay}')
        result={'case':case,'status':'passed','dimension':action.dim,'hamiltonian_nnz':int(nnz),'reference_distance':distance,'jw_magnitude_r100_state_distance':replay,'seconds':time.monotonic()-start,'source_sha256':sha(__file__),'core_sha256':sha(Path(__file__).with_name('core.py')),'input_design_sha256':__import__('hashlib').sha256(db).hexdigest(),'input_states_sha256':__import__('hashlib').sha256(sb).hexdigest()}
        np.savez_compressed(a.output/f'{case}.npz',exact=reference,initial=arrays['initial'])
        output.write_text(json.dumps(result,indent=2)+'\n');records.append(result);print(case,action.dim,distance,replay,flush=True)
    (a.output/'summary.json').write_text(json.dumps({'status':'passed','cases':len(records),'max_reference_distance':max(x['reference_distance'] for x in records),'checks':records},indent=2)+'\n')
if __name__=='__main__':main()
