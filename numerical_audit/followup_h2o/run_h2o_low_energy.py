#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='2'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-h2o-numba');os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-h2o-mpl')
import argparse,json,hashlib,shutil,resource,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from qhat.numerical_audit.low_energy_state_sensitivity import low_energy_states
from qhat.numerical_audit import low_energy_state_sensitivity as prep
from run_h2o_sweep import require,sha,array_sha,write,metrics,ev,METHODS,progress

CASES=('H2O_s-1.00_sto-6g_as-004-004','H2O_s-1.00_hgbs-5_as-004-004')

def pauli(key,n=8):
    d={'I':np.eye(2,dtype=complex),'X':np.array([[0,1],[1,0]],complex),
       'Y':np.array([[0,-1j],[1j,0]],complex),'Z':np.diag([1.,-1.]).astype(complex)}
    ops=dict(key);result=np.ones((1,1),complex)
    for p in range(n):result=np.kron(result,d[ops.get(p,'I')])
    return result

def spin_squared(n):
    sp=np.zeros((2**n,2**n),complex)
    for b in range(2**n):
        for p in range(0,n,2):
            x=b;sgn=1
            for mode,creation in ((p+1,False),(p,True)):
                bit=1<<(n-1-mode)
                if bool(x&bit)==creation:sgn=0;break
                if (x>>(n-mode)).bit_count()%2:sgn=-sgn
                x^=bit
            if sgn:sp[x,b]+=sgn
    am=sum(1<<(n-1-p) for p in range(0,n,2));bm=sum(1<<(n-1-p) for p in range(1,n,2))
    sz=np.array([((b&am).bit_count()-(b&bm).bit_count())/2 for b in range(2**n)])
    return np.diag(sz**2)+.5*(sp@sp.conj().T+sp.conj().T@sp)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sweep',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--qhat-root',type=Path,default=Path('/home/bernalde/repos/qhat'));a=ap.parse_args()
    require(all(Path(m.__file__).resolve().is_relative_to(a.qhat_root.resolve()) for m in (prep,ev)), 'Imported unintended QHAT checkout')
    a.output.mkdir(parents=True,exist_ok=False)
    for p in (Path(__file__),Path(__file__).with_name('LOW_ENERGY_PROTOCOL.md'),Path(__file__).with_name('run_h2o_sweep.py'),Path(prep.__file__),Path(prep.__file__).with_name('initial_state_sensitivity.py')):
        shutil.copyfile(p,a.output/p.name)
    write(a.output/'protocol.json',{'script_sha256':sha(__file__),'protocol_sha256':sha(Path(__file__).with_name('LOW_ENERGY_PROTOCOL.md')),
                                 'state_preparation_source_sha256':sha(prep.__file__),'shared_helper_sha256':sha(Path(__file__).with_name('run_h2o_sweep.py')),'qhat_root':str(a.qhat_root.resolve()),'cases':CASES,'rows':36,'T':1.,'r':100})
    n=8;dim=256;S2=spin_squared(n)
    am=sum(1<<(n-1-p) for p in range(0,n,2));bm=sum(1<<(n-1-p) for p in range(1,n,2))
    sm=np.array([(b&am).bit_count()==2 and (b&bm).bit_count()==2 for b in range(dim)])
    nm=np.array([b.bit_count()==4 for b in range(dim)]);sector=np.flatnonzero(sm)
    allrows=[]
    for case in CASES:
      for policy in ('historical_policy','unpruned'):
        out=a.output/f'{case}__{policy}';out.mkdir();source=a.sweep/out.name/'design.json';design=json.loads(source.read_text())
        co={tuple((int(q),p) for q,p in key):complex(float.fromhex(r),float.fromhex(i)) for key,r,i in design['identity_free_coefficients']}
        keys=[tuple((int(q),p) for q,p in key) for key in design['pauli_order_key_reference']]
        orders={name:[keys[i] for i in order] for name,order in design['pauli_orders'].items()}
        P={k:pauli(k,n) for k in keys};H=sum((co[k]*P[k] for k in keys),np.zeros((dim,dim),complex))
        require(np.linalg.norm(H-H.conj().T)<1e-12,'Nonhermitian H')
        small,spectrum,values,vectors=low_energy_states(H[np.ix_(sector,sector)])
        full={name:np.zeros(dim,complex) for name in small}
        info={}
        for name,vec in small.items():
            v=full[name];v[sector]=vec;s=float(np.vdot(v,S2@v).real);e=float(np.vdot(v,H@v).real)
            info[name]={'state_sha256':array_sha(v),'S2_expectation':s,'S2_variance':float(np.linalg.norm(S2@v-s*v)**2),
                        'energy_identity_free_hartree':e,'energy_std_hartree':float(np.linalg.norm(H@v-e*v))}
        state_design={'source_design_sha256':sha(source),'case':case,'tensor_policy':policy,'spectrum':spectrum,'states':info,
                      'sector_dimension':len(sector),'state_preparation':'Full fixed Nalpha=Nbeta=2 sector, total spin measured',
                      'S2_commutator_frobenius':float(np.linalg.norm(H@S2-S2@H)),
                      'off_sector_frobenius':float(np.linalg.norm(H[~sm][:,sm]))}
        write(out/'design.json',state_design);np.savez_compressed(out/'initial_states.npz',**full)
        energies,vectors=eigh(H);saved={f'initial_{k}':v for k,v in full.items()};rows=[];checks={}
        exact={}
        for name,v in full.items():
            exact[name]=expm_multiply(-1j*sparse.csr_matrix(H),v,traceA=-1j*np.trace(H))
            alternate=vectors@(np.exp(-1j*energies)*(vectors.conj().T@v))
            residual=float(np.linalg.norm(exact[name]-alternate));require(residual<1e-10,'Reference disagreement')
            checks[f'{name}_reference_dense_eigh_distance']=residual;saved[f'exact_{name}']=exact[name]
        for method,order in orders.items():
            U=np.eye(dim,dtype=complex)
            for key in order:
                angle=co[key].real/100;U=np.cos(angle)*U-1j*np.sin(angle)*(P[key]@U)
            Un=np.linalg.matrix_power(U,100)
            compiled=ev.compile_ordered_terms(order,co,n,1e-12)
            for name,v in full.items():
                actual,count=ev.evolve_trotter_state(v,compiled,1,100,1.,2**30)
                alternate=Un@v;residual=float(np.linalg.norm(actual-alternate));require(residual<1e-10,'Dense product disagreement')
                row={'case':case,'tensor_policy':policy,'state':name,'ordering':method,'T':1.,'steps':100,'n_qubits':8,'n_electrons':4,
                     'nominal_exponentials':count,'dense_product_state_distance':residual,'state_array_sha256':array_sha(actual),
                     **metrics(exact[name],actual,nm,sm)}
                saved[f'{method}_{name}']=actual;rows.append(row)
        np.savez_compressed(out/'states.npz',**saved)
        result={'rows':rows,'checks':checks,'states_sha256':sha(out/'states.npz'),'state_array_sha256':{k:array_sha(v) for k,v in saved.items()},
                'design_sha256':sha(out/'design.json'),'initial_states_sha256':sha(out/'initial_states.npz')}
        write(out/'results.json',result);allrows+=rows;progress(event='low_energy_complete',case=case,policy=policy,spin=info)
    write(a.output/'rows.json',allrows)
if __name__=='__main__':main()
