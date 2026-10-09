#!/usr/bin/env python3
"""Independent NumPy-only hash, state-error and leakage replay; no QHAT imports."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np

def require(x,msg):
    if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def objsha(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def arrsha(arr):
    return objsha({'shape':arr.shape,'dtype':arr.dtype.str,'raw_sha256':hashlib.sha256(arr.tobytes()).hexdigest()})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    summaries=[];maxmetric=maxleak=maxnorm=0.;count=0
    for directory in sorted(args.root.glob('*q')):
        if not (directory/'results.json').exists():continue
        result=json.loads((directory/'results.json').read_text());design=json.loads((directory/'frozen_state_design.json').read_text())
        require(sha(args.root/'PROTOCOL.json')==design['protocol_sha256'],'Protocol hash mismatch')
        require(sha(args.root/'run_state_control.py')==design['script_sha256'],'Execution source hash mismatch')
        for file,digest in (('states.npz',result['states_file_sha256']),('frozen_state_design.json',result['frozen_state_design_sha256']),
                            ('reference_checks.json',result['reference_checks_sha256']),('initial_states.npz',design['initial_states_file_sha256'])):
            require(sha(directory/file)==digest,f'File hash mismatch {directory/file}')
        require(sha(directory/'source_design.json')==design['source_design_sha256'],'Source design hash mismatch')
        n=design['n_qubits'];alpha=sum(1<<(n-1-q) for q in range(0,n,2));beta=sum(1<<(n-1-q) for q in range(1,n,2))
        indices=range(1<<n);number=np.array([i.bit_count()==6 for i in indices]);spin=np.array([(i&alpha).bit_count()==3 and (i&beta).bit_count()==3 for i in indices])
        with np.load(directory/'states.npz',allow_pickle=False) as saved,np.load(directory/'initial_states.npz',allow_pickle=False) as initial:
            for key,digest in result['state_arrays_sha256'].items():require(arrsha(saved[key])==digest,f'Array hash mismatch {key}')
            for name in design['state_info']:
                require(np.array_equal(saved['initial_'+name],initial[name]),'Frozen initial state changed')
                require(arrsha(initial[name])==design['state_info'][name]['state_sha256'],'Initial state array hash mismatch')
            groups={}
            for row in result['rows']:
                name=row['state'];method=row['ordering'];u=saved['exact_'+name];v=saved[method+'_'+name]
                require(arrsha(v)==row['state_sha256'],'Row state hash mismatch')
                nu=np.linalg.norm(u);nv=np.linalg.norm(v);u=u/nu;v=v/nv;overlap=np.vdot(u,v)
                require(abs(overlap)>0,'Unexpected orthogonal states')
                aligned=v*np.conjugate(overlap)/abs(overlap)
                E=.5*float(np.vdot(u-aligned,u-aligned).real)
                infidelity=E*(2-E)
                p=float(np.vdot(v[~number],v[~number]).real);s=float(np.vdot(v[~spin],v[~spin]).real)
                md=max(abs(E-row['one_minus_overlap']),abs(infidelity-row['infidelity']))
                ld=max(abs(p-row['particle_number_leakage']),abs(s-row['spin_sector_leakage']))
                require(md<1e-16,'Independent metric mismatch');require(ld<1e-20,'Independent leakage mismatch')
                maxmetric=max(maxmetric,md);maxleak=max(maxleak,ld);maxnorm=max(maxnorm,abs(nu-1),abs(nv-1));count+=1
                groups.setdefault(name,{})[method]=E
            for name,errors in groups.items():
                own=errors['fermionic_signed']
                summaries.append({'case':design['case'],'state':name,'errors':errors,
                    'owned_over_comparator':{method:own/errors[method] for method in ('jw_signed','jw_magnitude')},
                    'outcome':{method:'tie' if abs(own-errors[method])<=1e-12 else ('win' if own<errors[method] else 'loss') for method in ('jw_signed','jw_magnitude')}})
    require(count>0,'No completed runs')
    result={'completed_new_evolutions':count,'replayed':'Saved vectors, hashes, normalized phase-insensitive errors and leakage. No independent full product-formula propagation.',
            'max_independent_metric_difference':maxmetric,'max_independent_leakage_difference':maxleak,'max_norm_defect':maxnorm,'comparisons':summaries}
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
