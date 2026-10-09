#!/usr/bin/env python3
"""Frozen H2O reconstruction sweep. See PROTOCOL.md; no historical-input claim."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[k]='2'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-h2o-numba')
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-h2o-mpl')
import argparse,hashlib,json,resource,shutil,sys,time,importlib.metadata
from pathlib import Path
cap=12*1024**3
soft,hard=resource.getrlimit(resource.RLIMIT_AS)
resource.setrlimit(resource.RLIMIT_AS,(min(cap,hard) if hard>=0 else cap,hard))
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as bl
from qhat.analysis.fermionic_coefficients import fermion_from_interaction,PRECISION_POLICY_VERSION
from qhat.analysis.fermionic_ownership import build_molecular_ownership
from qhat.analysis.state_metrics import state_error_metrics
from qhat.numerical_audit import reference_tools as ref

METHODS=('fermionic_signed','jw_signed','jw_magnitude')
ANCHORS=('H2O_s-1.00_sto-6g_as-004-004','H2O_s-1.00_hgbs-5_as-004-008','H2O_s-1.00_hgbs-5_as-006-010')
ATOL=1e-12

def require(x,msg):
    if not x: raise RuntimeError(msg)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def object_sha(o): return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def array_sha(a):
    a=np.asarray(a)
    return object_sha({'shape':a.shape,'dtype':a.dtype.str,'raw_sha256':hashlib.sha256(a.tobytes()).hexdigest()})
def write(p,o):
    with Path(p).open('x') as f: json.dump(o,f,indent=2,allow_nan=False);f.write('\n')
def progress(**x): print(json.dumps(x),flush=True)
def canonical(co): return [[k,complex(v).real.hex(),complex(v).imag.hex()] for k,v in sorted(co.items())]
def metrics(exact,actual,nm,sm):
    m=state_error_metrics(exact,actual)
    require(max(abs(m['exact_norm']-1),abs(m['approx_norm']-1))<1e-8,'Excess norm drift')
    z=actual/m['approx_norm']
    m.update(particle_number_leakage=float(np.vdot(z[~nm],z[~nm]).real),
             spin_sector_leakage=float(np.vdot(z[~sm],z[~sm]).real),
             naive_unnormalized_one_minus_overlap=max(0.,1.-min(1.,float(abs(np.vdot(exact,actual))))),
             below_analysis_floor=m['one_minus_overlap']<=ATOL)
    return m

def sources(qhat):
    result={'run_h2o_sweep.py':Path(__file__),'PROTOCOL.md':Path(__file__).with_name('PROTOCOL.md')}
    for rel in ('analysis/benchmark_b2_active_spaces_matrix_free.py','analysis/benchmark_b2_signed_coefficient_baseline.py',
                'analysis/fermionic_coefficients.py','analysis/fermionic_ownership.py','analysis/state_metrics.py',
                'analysis/benchmark_L_sweep_trotter.py','numerical_audit/reference_tools.py'):
        result[rel]=qhat/rel
    return result

def prepare(case,policy,path,out,pilot_root):
    started=time.monotonic(); interaction,n=bl.load_interaction_operator(path)
    electrons=int(case.split('_as-')[1].split('-')[0]); hf=((1<<electrons)-1)<<(n-electrons)
    owned=build_molecular_ownership(fermion_from_interaction(interaction),atol=ATOL,ordering='signed_ascending')
    require(not owned.fallback_pauli_keys and not owned.diagnostics['noncommuting_owned_block_count'],'Invalid ownership')
    co=owned.final_coefficients; keys=list(co)
    require(all(complex(c).imag==0 for c in co.values()),'Nonreal Pauli coefficient')
    schedules={'fermionic_signed':list(owned.pauli_order),
               'jw_signed':bl.signed_coefficient_lexicographic_order(keys,co,n,ATOL),
               'jw_magnitude':bl.magnitude_descending_lexicographic_order(keys,co,n,ATOL)}
    require(all(len(order)==len(keys) and set(order)==set(keys) for order in schedules.values()),'Changed common H')
    raw=ev.compile_ordered_terms(keys,co,n,ATOL)
    require(len(raw)==len(keys),'Compilation changed coefficients')
    states,reduced,basis=ref.affine_space(raw,hf)
    index={key:i for i,key in enumerate(keys)}
    compiled={name:[reduced[index[k]] for k in order] for name,order in schedules.items()}
    nm=np.array([int(b).bit_count()==electrons for b in states])
    am=sum(1<<(n-1-p) for p in range(0,n,2)); bm=sum(1<<(n-1-p) for p in range(1,n,2))
    sm=np.array([(int(b)&am).bit_count()==electrons//2 and (int(b)&bm).bit_count()==electrons//2 for b in states])
    initial=np.zeros(len(states),complex);initial[0]=1
    design={'case':case,'tensor_policy':policy,'n_qubits':n,'n_electrons':electrons,
            'basis':case.split('_')[2],'tensor_sha256':sha(path),'tensor_path':str(path),
            'pauli_terms':len(keys),'owned_blocks':len(owned.owned_blocks),'coefficient_cutoff':ATOL,
            'coefficient_policy':PRECISION_POLICY_VERSION,'identity_free_coefficients':canonical(co),
            'hamiltonian_sha256':object_sha(canonical(co)),
            'discarded_global_scalar':[complex(owned.scalar).real,complex(owned.scalar).imag],
            'pauli_order_key_reference':keys,'pauli_orders':{name:[index[k] for k in order] for name,order in schedules.items()},
            'schedule_sha256':{name:object_sha([[k,complex(co[k]).real.hex()] for k in order]) for name,order in schedules.items()},
            'closed_dimension':len(states),'full_hilbert_dimension':2**n,'affine_flip_basis':basis,
            'affine_basis_indices':list(map(int,states)),'physical_spin_sector_in_closure_dimension':int(sum(sm)),
            'initial_state':'HF first n_electrons modes occupied, q0 most significant; no postselection',
            'initial_array_sha256':array_sha(initial),'order_convention':'chronological; first listed factor acts first',
            'planned_steps':[50,100,200] if case in ANCHORS else [100],
            'ownership_diagnostics':owned.diagnostics}
    write(out/'design.json',design)
    progress(event='prepared',case=case,policy=policy,n=n,paulis=len(keys),closed_dimension=len(states),seconds=time.monotonic()-started)
    H,_=ref.assemble_hamiltonian(reduced,len(states))
    checks={'hermiticity_frobenius':float(sparse.linalg.norm(H-H.getH())),
            'off_spin_sector_frobenius':float(sparse.linalg.norm(H[~sm][:,sm])),
            'reference_nnz':H.nnz}
    # Check an arbitrary sample of Pauli rotations against the full Hilbert space,
    # retaining all leakage; each large case checks the affine sign convention.
    chosen=np.unique(np.linspace(0,len(keys)-1,min(12,len(keys)),dtype=int))
    small,_=ev.evolve_trotter_state(initial,[reduced[i] for i in chosen],1,1,.031,2**30)
    full0=np.zeros(2**n,complex);full0[hf]=1
    full,_=ev.evolve_trotter_state(full0,[raw[i] for i in chosen],1,1,.031,2**30)
    embed=np.zeros_like(full);embed[states]=small
    checks['sample_pauli_full_vs_affine_distance']=float(np.linalg.norm(full-embed))
    require(checks['sample_pauli_full_vs_affine_distance']<1e-12,'Affine sign/embedding failure')
    del full,full0,embed
    exact_start=time.monotonic();trace=-1j*H.diagonal().sum()
    exact=expm_multiply(-1j*H,initial,traceA=trace)
    half=expm_multiply(-.5j*H,initial,traceA=trace/2)
    split=expm_multiply(-.5j*H,half,traceA=trace/2)
    checks['reference_split_time_state_distance']=float(np.linalg.norm(split-exact))
    checks['reference_norm']=float(np.linalg.norm(exact))
    checks['reference_elapsed_seconds']=time.monotonic()-exact_start
    require(checks['reference_split_time_state_distance']<1e-10,'Reference split-time disagreement')
    if n<=6:
        pilot=pilot_root/f'{case}__{policy}.states.npz'
        with np.load(pilot,allow_pickle=False) as f:
            checks['reference_vs_independent_dense_pilot_distance']=float(np.linalg.norm(exact-f['exact_retained'][states]))
            require(checks['reference_vs_independent_dense_pilot_distance']<1e-10,'Reference disagrees with dense pilot')
    write(out/'checks.json',checks)
    progress(event='reference_ready',case=case,policy=policy,checks=checks)
    return initial,exact,compiled,nm,sm,design,states

def run_one(case,policy,path,out,pilot_root):
    started=time.monotonic()
    initial,exact,compiled,nm,sm,design,states=prepare(case,policy,path,out,pilot_root)
    saved={'initial':initial,'exact':exact,'affine_basis_indices':states}
    rows=[]
    for steps in design['planned_steps']:
        for name in METHODS:
            t=time.monotonic();actual,count=ev.evolve_trotter_state(initial,compiled[name],1,steps,1.,2**30)
            require(count==steps*design['pauli_terms'],'Incorrect exponential count')
            row={'case':case,'basis':design['basis'],'tensor_policy':policy,'n_qubits':design['n_qubits'],
                 'n_electrons':design['n_electrons'],'T':1.,'steps':steps,'role':'primary' if steps==100 else 'convergence_anchor',
                 'ordering':name,'nominal_exponentials':count,'elapsed_seconds':time.monotonic()-t,
                 'state_array_sha256':array_sha(actual),**metrics(exact,actual,nm,sm)}
            if design['n_qubits']<=6:
                with np.load(pilot_root/f'{case}__{policy}.states.npz',allow_pickle=False) as f:
                    row['vs_independent_dense_pilot_state_distance']=float(np.linalg.norm(actual-f[f'{name}_r{steps}_dense'][states]))
                    require(row['vs_independent_dense_pilot_state_distance']<1e-10,'Compiled evolution disagrees with dense pilot')
            saved[f'{name}_r{steps}']=actual;rows.append(row)
            progress(event='evolved',case=case,policy=policy,ordering=name,r=steps,E=row['one_minus_overlap'],seconds=row['elapsed_seconds'])
    np.savez_compressed(out/'states.npz',**saved)
    result={'case':case,'tensor_policy':policy,'rows':rows,'elapsed_seconds':time.monotonic()-started,
            'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'states_sha256':sha(out/'states.npz'),'state_array_sha256':{k:array_sha(v) for k,v in saved.items()},
            'design_sha256':sha(out/'design.json'),'checks_sha256':sha(out/'checks.json')}
    write(out/'results.json',result)
    require(result['peak_rss_KiB']<10*1024**2,'Peak RSS exceeded target cap')
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--qhat-root',type=Path,default=Path('/home/bernalde/repos/qhat'))
    ap.add_argument('--data',type=Path,default=Path('/home/bernalde/repos/qhat/numerical_audit/h2o_recovery/data'))
    ap.add_argument('--pilot',type=Path,default=Path('/home/bernalde/repos/qhat/numerical_audit/h2o_recovery/validation'))
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--max-qubits',type=int,default=20)
    args=ap.parse_args();qhat=args.qhat_root.resolve()
    require(all(Path(m.__file__).resolve().is_relative_to(qhat) for m in (ev,bl,ref)),'Imported unintended QHAT checkout')
    src=sources(qhat);source_hashes={k:sha(p) for k,p in src.items()}
    manifest_path=args.data/'manifest.json';manifest=json.loads(manifest_path.read_text())
    protocol={'version':'h2o-reconstruction-sweep-v1','protocol_sha256':source_hashes['PROTOCOL.md'],
              'source_sha256':source_hashes,'reconstruction_manifest_sha256':sha(manifest_path),
              'methods':METHODS,'anchors':ANCHORS,'primary_steps':100,'T':1.,'planned_rows':126,
              'coefficient_cutoff':ATOL,'tie_absolute':1e-12,'analysis_floor':1e-12,
              'python':sys.version,'packages':{n:importlib.metadata.version(n) for n in ('numpy','scipy','numba','openfermion')},
              'address_space_limit_bytes':resource.getrlimit(resource.RLIMIT_AS)[0]}
    args.output.mkdir(parents=True,exist_ok=True)
    if (args.output/'protocol.json').exists():
        require(json.loads((args.output/'protocol.json').read_text())==json.loads(json.dumps(protocol)),'Resume protocol/source/input mismatch')
    else:
        write(args.output/'protocol.json',protocol)
        for rel,p in src.items():
            target=args.output/'source'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    planned=[]
    for basis,bm in manifest['bases'].items():
        for rel,digest in bm['files_sha256'].items():
            if rel.endswith('.tensors.npz'):
                path=args.data/basis/rel;require(sha(path)==digest,f'Changed input {path}')
                case=path.name.removesuffix('.tensors.npz');n=sum(map(int,case.split('_as-')[1].split('-')))
                planned.append((n,basis,case,rel.split('/')[0],path))
    require(len(planned)==30,'Expected 30 tensor files')
    for n,basis,case,policy,path in sorted(planned):
        if n>args.max_qubits:continue
        out=args.output/f'{case}__{policy}'
        if (out/'results.json').exists():
            old=json.loads((out/'results.json').read_text())
            for file,key in (('states.npz','states_sha256'),('design.json','design_sha256'),('checks.json','checks_sha256')):
                require(sha(out/file)==old[key],f'Changed completed {out/file}')
            progress(event='verified_previous_case',case=case,policy=policy);continue
        require(not out.exists(),f'Incomplete result directory, inspect before restarting {out}')
        require({k:sha(p) for k,p in src.items()}==source_hashes,'Source changed midrun')
        out.mkdir();run_one(case,policy,path,out,args.pilot)
    progress(event='batch_complete',max_qubits=args.max_qubits,output=str(args.output))
if __name__=='__main__':main()
