#!/usr/bin/env python3
"""Independent NumPy/SciPy replay; no QHAT, OpenFermion or Numba imports.

Read only frozen designs and saved arrays. Independently reconstruct Pauli
actions in the explicitly saved closure, reference evolution and all circuits.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time
import numpy as np
from scipy.linalg import eigh


def key(raw):return tuple((int(q),str(p)) for q,p in raw)

def block(raw):return {key(k):int(v) for k,v in raw}

def pauli_data(word,basis,n):
    flip=0;z=0;ys=0
    for q,p in word:
        if p in ('X','Y'):flip ^= 1<<(n-1-q)
        if p in ('Z','Y'):z ^= 1<<(n-1-q)
        if p=='Y':ys+=1
    coordinates={int(v):i for i,v in enumerate(basis)}
    source=np.array([coordinates[int(v)^flip] for v in basis])
    phase=(1j**ys)*np.array([(-1 if ((int(v)^flip)&z).bit_count()%2 else 1) for v in basis])
    return source,phase

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    if args.output.exists():raise ValueError('Output exists.')
    manifest=json.loads((args.run/'manifest.json').read_text());checks=[];start=time.monotonic()
    for case in manifest['cases']:
        design=json.loads((args.run/f'{case}.design.json').read_text());den=int(design['denominator']);total=block(design['coefficients'])
        n=design['n_qubits'];basis=design['affine_basis_indices'];dim=len(basis)
        blocks={kind:[block(b) for b in bs] for kind,bs in design['blocks'].items()}
        words=set(total).union(*(set(b) for b in blocks['parent']))
        action={k:pauli_data(k,basis,n) for k in words}
        H=np.zeros((dim,dim),complex)
        for k,v in total.items():
            src,phase=action[k];H[np.arange(dim),src] += (v/den)*phase
        if np.linalg.norm(H-H.conj().T)>1e-11:raise ValueError('Non-Hermitian saved target.')
        with np.load(args.run/f'{case}.states.npz') as saved:
            initial=saved['initial'];values,vectors=eigh(H)
            exact=vectors@(np.exp(-1j*values)*(vectors.conj().T@initial))
            exact_delta=float(np.linalg.norm(exact-saved['exact']))
            if exact_delta>2e-11:raise RuntimeError('Independent exact reference fails.')
            for run in design['planned_runs']:
                name,r,role=run['ordering'],run['steps'],run['role']
                if name.startswith('jw_'):
                    schedule=[(key(k),total[key(k)]) for k in design['jw_orders'][name]]
                else:
                    kind=name.split('_')[0]
                    schedule=[(k,v) for i in design['orders'][name] for k,v in blocks[kind][i].items()]
                prepared=[]
                for k,v in schedule:
                    src,phase=action[k];angle=(v/den)/r
                    prepared.append((src,np.cos(angle),-1j*np.sin(angle)*phase))
                state=initial.copy()
                for _ in range(r):
                    for src,c,s in prepared:state=c*state+s*state[src]
                archived=saved[f'{name}__r{r}__{role}'];distance=float(np.linalg.norm(state-archived))
                if distance>2e-11:raise RuntimeError(f'Independent replay failed {case}/{name}/{r}: {distance}')
                a=exact/np.linalg.norm(exact);b=state/np.linalg.norm(state);z=np.vdot(a,b)
                phase=z/abs(z) if abs(z) else 1.
                error=.5*float(np.vdot(b-phase*a,b-phase*a).real)
                row=json.loads((args.run/f'{case}.{name}.r{r}.{role}.json').read_text())
                error_delta=abs(error-row['one_minus_overlap'])
                if error_delta>2e-12:raise RuntimeError('Independent stable metric differs.')
                checks.append({'case':case,'ordering':name,'r':r,'role':role,'state_distance':distance,
                    'error_absolute_difference':error_delta,'reference_distance':exact_delta})
        print('RECHECKED',case,len(checks),flush=True)
    result={'status':'passed','rows':len(checks),'checks':checks,'seconds':time.monotonic()-start,
        'max_state_distance':max(c['state_distance'] for c in checks),
        'max_metric_absolute_difference':max(c['error_absolute_difference'] for c in checks),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Independent NumPy/SciPy full-closure coefficient/reference/circuit replay; no QHAT imports. It checks these saved designs, not the physical molecular tensor-generation assumptions.'}
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':main()
