#!/usr/bin/env python3
"""Independently remeasure saved sweep/low-energy states, no QHAT imports."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np

def require(c,msg):
    if not c:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def object_sha(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def array_sha(a):
    a=np.asarray(a);return object_sha({'shape':a.shape,'dtype':a.dtype.str,'raw_sha256':hashlib.sha256(a.tobytes()).hexdigest()})
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sweep',type=Path,required=True);ap.add_argument('--low-energy',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--data',type=Path,default=Path('/home/bernalde/repos/qhat/numerical_audit/h2o_recovery/data'));a=ap.parse_args()
    protocol=json.loads((a.sweep/'protocol.json').read_text())
    require(sha(a.data/'manifest.json')==protocol['reconstruction_manifest_sha256'],'Input manifest changed')
    for rel,digest in protocol['source_sha256'].items():require(sha(a.sweep/'source'/rel)==digest,'Frozen source changed')
    lp=json.loads((a.low_energy/'protocol.json').read_text())
    require(sha(a.low_energy/'run_h2o_low_energy.py')==lp['script_sha256'],'Low-energy source changed')
    require(sha(a.low_energy/'LOW_ENERGY_PROTOCOL.md')==lp['protocol_sha256'],'Low-energy protocol changed')
    require(sha(a.low_energy/'run_h2o_sweep.py')==lp['shared_helper_sha256'],'Low-energy helper changed')
    require(sha(a.low_energy/'low_energy_state_sensitivity.py')==lp['state_preparation_source_sha256'],'Low-energy preparation changed')
    checked=0;arrays=0;max_E_delta=0.;max_I_delta=0.;max_leak_delta=0.;max_norm_defect=0.;records=[]
    for is_low,root in ((False,a.sweep),(True,a.low_energy)):
        for path in sorted(root.glob('H2O*/results.json')):
            out=path.parent;r=json.loads(path.read_text());design=json.loads((out/'design.json').read_text())
            require(sha(out/'states.npz')==r['states_sha256'],'State archive hash mismatch')
            require(sha(out/'design.json')==r['design_sha256'],'Design hash mismatch')
            if not is_low:require(sha(out/'checks.json')==r['checks_sha256'],'Checks hash mismatch')
            with np.load(out/'states.npz',allow_pickle=False) as saved:
                require(set(saved.files)==set(r['state_array_sha256']),'Array membership changed')
                for key,digest in r['state_array_sha256'].items():require(array_sha(saved[key])==digest,'State array hash mismatch');arrays+=1
                n=8 if is_low else design['n_qubits'];ne=4 if is_low else design['n_electrons']
                basis=np.arange(2**n) if is_low else saved['affine_basis_indices']
                require(len(set(map(int,basis)))==len(basis),'Duplicate embedding index')
                if not is_low:
                    tensor=a.data/design['basis']/design['tensor_policy']/(design['case']+'.tensors.npz')
                    require(sha(tensor)==design['tensor_sha256'],'Tensor hash changed')
                    require(list(map(int,basis))==design['affine_basis_indices'],'Embedding differs from design')
                    co={tuple((int(q),p) for q,p in key):complex(float.fromhex(re),float.fromhex(im)) for key,re,im in design['identity_free_coefficients']}
                    require(object_sha(design['identity_free_coefficients'])==design['hamiltonian_sha256'],'H hash mismatch')
                    keys=[tuple((int(q),p) for q,p in key) for key in design['pauli_order_key_reference']]
                    for method,order in design['pauli_orders'].items():
                        require(sorted(order)==list(range(len(keys))),'Nonpermutation schedule')
                        require(object_sha([[keys[i],co[keys[i]].real.hex()] for i in order])==design['schedule_sha256'][method],'Schedule hash mismatch')
                    points=set(map(int,basis));flips=set(sum(1<<(n-1-q) for q,p in k if p in 'XY') for k in keys)
                    require(all(all((int(s)^f) in points for s in basis) for f in flips),'Embedding is not Pauli invariant')
                am=sum(1<<(n-1-p) for p in range(0,n,2));bm=sum(1<<(n-1-p) for p in range(1,n,2))
                nm=np.array([int(s).bit_count()==ne for s in basis]);sm=np.array([(int(s)&am).bit_count()==ne//2 and (int(s)&bm).bit_count()==ne//2 for s in basis])
                for row in r['rows']:
                    exact=saved[f'exact_{row["state"]}'] if is_low else saved['exact']
                    approx=saved[f'{row["ordering"]}_{row["state"]}'] if is_low else saved[f'{row["ordering"]}_r{row["steps"]}']
                    require(array_sha(approx)==row['state_array_sha256'],'Row state hash mismatch')
                    en=float(np.sqrt(np.vdot(exact,exact).real));an=float(np.sqrt(np.vdot(approx,approx).real));x=exact/en;y=approx/an
                    # Projector residual gives infidelity directly without phase alignment
                    # or subtraction from one; convert to E by rationalization.
                    residual=y-x*np.vdot(x,y);I=float(np.vdot(residual,residual).real)
                    require(-1e-14<=I<=1+1e-12,'Invalid projector residual')
                    E=I/(1+math.sqrt(max(0.,1-I)))
                    de=abs(E-row['one_minus_overlap']);di=abs(I-row['infidelity'])
                    require(de<2e-12*max(E,1e-16)+1e-22,'Independent E disagrees')
                    require(di<4e-12*max(I,1e-16)+2e-22,'Independent infidelity disagrees')
                    max_E_delta=max(max_E_delta,de);max_I_delta=max(max_I_delta,di)
                    for mask,key in ((nm,'particle_number_leakage'),(sm,'spin_sector_leakage')):
                        leak=float(sum(abs(v)**2 for v in y[~mask]));delta=abs(leak-row[key]);max_leak_delta=max(max_leak_delta,delta)
                        require(delta<1e-18+1e-12*leak,'Leakage mismatch')
                    max_norm_defect=max(max_norm_defect,abs(en-1),abs(an-1));checked+=1
            records.append({'path':str(path),'sha256':sha(path)})
    require(checked==162,f'Expected 162 states, found {checked}')
    report={'passed':True,'checked_approximate_states':checked,'checked_saved_arrays':arrays,
            'independent_E_max_absolute_difference':max_E_delta,'independent_infidelity_max_absolute_difference':max_I_delta,
            'leakage_max_absolute_difference':max_leak_delta,'maximum_norm_defect':max_norm_defect,
            'checked_files':records,'checker_sha256':sha(__file__),'qhat_imported':False,
            'scope':'State/hash/metric/closure verification, not independent re-evolution of all large states.'}
    with a.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='checked_files'}))
if __name__=='__main__':main()
