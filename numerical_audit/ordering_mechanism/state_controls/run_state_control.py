#!/usr/bin/env python3
"""Frozen larger-state control. See PROTOCOL.json before interpreting results."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='2'
os.environ['NUMBA_CACHE_DIR']='/tmp/ordering-large-state-control/numba-cache'
os.environ['MPLCONFIGDIR']='/tmp/ordering-large-state-control/mpl-cache'
import argparse, hashlib, json, resource, shutil, subprocess, time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, expm_multiply
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
from qhat.analysis.state_metrics import state_error_metrics

ROOT=Path('/home/bernalde/repos/qhat')
METHODS=('fermionic_signed','jw_signed','jw_magnitude')
CLUSTER=1e-10

def require(condition,message):
    if not condition: raise RuntimeError(message)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def objsha(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def arrsha(arr):
    arr=np.asarray(arr)
    return objsha({'shape':arr.shape,'dtype':arr.dtype.str,'raw_sha256':hashlib.sha256(arr.tobytes()).hexdigest()})
def write(path,obj):
    with Path(path).open('x') as out: json.dump(obj,out,indent=2,allow_nan=False);out.write('\n')
def event(**obj): print(json.dumps(obj),flush=True)

def phase(vector):
    vector=np.asarray(vector,dtype=complex).copy();mag=abs(vector)
    pivot=np.flatnonzero(mag>=max(mag)-1e-12)[0]
    vector*=np.conjugate(vector[pivot])/abs(vector[pivot])
    return vector

def first_canonical(vectors):
    # First coordinate projection of the complete numerical eigenspace.
    # No materialized projector is needed for the first canonical vector.
    for index in range(vectors.shape[0]):
        candidate=vectors@vectors[index,:].conj()
        norm=np.linalg.norm(candidate)
        if norm>1e-10:return phase(candidate/norm)
    raise RuntimeError('Cannot canonicalize empty eigenspace')

def select_states(H):
    count=min(16,len(H));elapsed=time.monotonic()
    while True:
        vals,vecs=eigh(H,subset_by_index=(0,count-1),driver='evr')
        clusters=[];start=0
        while start<count:
            end=start+1
            while end<count and vals[end]-vals[start]<=CLUSTER:end+=1
            clusters.append((start,end));start=end
        if len(clusters)>=3 or (len(clusters)>=2 and count==len(H)):break
        require(count<len(H),'No next distinct energy exists')
        count=min(2*count,len(H))
    selected=[first_canonical(vecs[:,s:e]) for s,e in clusters[:2]]
    require(abs(np.vdot(*selected))<1e-10,'Nonorthogonal selected vectors')
    independent=eigsh(sparse.csr_matrix(H),k=min(4,len(H)-1),which='SA',
                      v0=np.sin(np.arange(len(H))+0.271),tol=1e-13,maxiter=100000,
                      return_eigenvectors=False)
    independent=np.sort(independent)
    require(max(abs(independent-vals[:len(independent)]))<1e-10,'Dense/ARPACK eigenvalue mismatch')
    return {'sector_ground':selected[0], 'sector_next_distinct_energy':selected[1],
            'equal_coherent':(selected[0]+selected[1])/np.sqrt(2)}, {
       'computed_low_eigenvalues':vals.tolist(),
       'selected_clusters':[{'start':s,'stop':e,'minimum':float(vals[s]),'maximum':float(vals[e-1])} for s,e in clusters[:2]],
       'next_means':'First canonical vector of next distinct energy cluster, not second listed eigenvector.',
       'independent_arpack_low_eigenvalues':independent.tolist(),
       'max_dense_arpack_eigenvalue_difference':float(max(abs(independent-vals[:len(independent)]))),
       'elapsed_seconds':time.monotonic()-elapsed}

def independent_hamiltonian(coeff,n):
    # Written Pauli letters act on each computational-basis column directly.
    cols=np.arange(1<<n,dtype=np.int64);groups={}
    for key,c in coeff.items():
        flip=0;amplitude=np.full(len(cols),c,dtype=complex)
        for q,letter in key:
            bit=1<<(n-1-q);occupied=(cols&bit)!=0
            if letter in ('X','Y'):flip^=bit
            if letter=='Y':amplitude*=np.where(occupied,-1j,1j)
            if letter=='Z':amplitude*=np.where(occupied,-1,1)
        if flip not in groups:groups[flip]=amplitude
        else:groups[flip]+=amplitude
    rows=[];columns=[];data=[]
    for flip,amplitude in groups.items():
        nonzero=np.flatnonzero(amplitude!=0)
        rows.append(nonzero^flip);columns.append(nonzero);data.append(amplitude[nonzero])
    return sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(columns))),
                             shape=(1<<n,1<<n)).tocsr()

def tensor_axis_pauli(v,key,n):
    matrices={'X':np.array([[0,1],[1,0]],complex),'Y':np.array([[0,-1j],[1j,0]],complex),
              'Z':np.diag([1,-1]).astype(complex)}
    state=v.reshape((2,)*n)
    for q,letter in key:
        state=np.moveaxis(np.tensordot(matrices[letter],state,axes=(1,q)),0,q)
    return state.reshape(-1)

def spin_tools(n):
    dim=1<<n;indices=np.arange(dim,dtype=np.int64)
    alpha=sum(1<<(n-1-q) for q in range(0,n,2));beta=sum(1<<(n-1-q) for q in range(1,n,2))
    na=np.array([int(b&alpha).bit_count() for b in indices]);nb=np.array([int(b&beta).bit_count() for b in indices])
    rows=[];cols=[]
    for q in range(0,n,2):
        a=1<<(n-1-q);b=1<<(n-2-q)
        source=indices[((indices&a)==0)&((indices&b)!=0)]
        rows.extend(source^(a|b));cols.extend(source)
    raising=sparse.coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(dim,dim)).tocsr()
    lowering=raising.T
    sz2=((na-nb)/2)**2
    def S2(v):return sz2*v+.5*(raising@(lowering@v)+lowering@(raising@v))
    return na,nb,S2

def diagnostics(v,H,S2,number_mask,spin_mask):
    z=v/np.linalg.norm(v);hv=H@z;s2v=S2(z)
    energy=float(np.vdot(z,hv).real);spin=float(np.vdot(z,s2v).real)
    return {'norm':float(np.linalg.norm(v)), 'energy_identity_free_hartree':energy,
            'energy_std_hartree':float(np.linalg.norm(hv-energy*z)),
            'S2_expectation':spin,'S2_variance':float(np.linalg.norm(s2v-spin*z)**2),
            'particle_number_leakage':float(np.vdot(z[~number_mask],z[~number_mask]).real),
            'spin_sector_leakage':float(np.vdot(z[~spin_mask],z[~spin_mask]).real)}

def run(case,out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'numerical_audit/followup_h2o/results'/f'{case}__historical_policy'
    design=json.loads((source/'design.json').read_text());n=design['n_qubits'];dim=1<<n
    require(n in (14,16) and design['n_electrons']==6,'Unplanned case')
    coeff={tuple((int(q),p) for q,p in key):complex(float.fromhex(r),float.fromhex(i))
           for key,r,i in design['identity_free_coefficients']}
    canonical=[[key,c.real.hex(),c.imag.hex()] for key,c in sorted(coeff.items())]
    require(objsha(canonical)==design['hamiltonian_sha256'],'Hamiltonian hash mismatch')
    require(all(c.imag==0 for c in coeff.values()),'Unexpected complex coefficients')
    keys=[tuple((int(q),p) for q,p in key) for key in design['pauli_order_key_reference']]
    schedules={name:[keys[i] for i in design['pauli_orders'][name]] for name in METHODS}
    for name,order in schedules.items():
        require(objsha([[key,coeff[key].real.hex()] for key in order])==design['schedule_sha256'][name],'Schedule mismatch')
        require(len(order)==len(coeff) and set(order)==set(coeff),'Changed common H')
    shutil.copyfile(source/'design.json',out/'source_design.json')
    compiled={name:ev.compile_ordered_terms(order,coeff,n,1e-12) for name,order in schedules.items()}
    H=independent_hamiltonian(coeff,n);na,nb,S2=spin_tools(n)
    sm=(na==3)&(nb==3);nm=(na+nb)==6;sector=np.flatnonzero(sm)
    checks={'hermiticity_frobenius':float(sparse.linalg.norm(H-H.getH())),
            'off_spin_sector_frobenius':float(sparse.linalg.norm(H[~sm][:,sm])),
            'full_dimension':dim,'sector_dimension':len(sector),'H_nnz':H.nnz}
    require(checks['hermiticity_frobenius']<1e-10,'Nonhermitian target')
    sub=H[sector][:,sector].toarray()
    small,spectrum=select_states(sub);del sub
    full={}
    for name,v in small.items():
        w=np.zeros(dim,complex);w[sector]=v;full[name]=w
    hf=np.zeros(dim,complex);hf[((1<<6)-1)<<(n-6)]=1
    state_info={name:{**diagnostics(v,H,S2,nm,sm),'state_sha256':arrsha(v),
                      'HF_fidelity':float(abs(np.vdot(hf,v))**2)} for name,v in full.items()}
    for name in ('sector_ground','sector_next_distinct_energy'):
        require(state_info[name]['energy_std_hartree']<2e-10,'Full H eigenvector residual too large')
    np.savez_compressed(out/'initial_states.npz',sector_indices=sector,**full)
    # This file and every initial vector are frozen before approximate evolutions.
    frozen={'protocol_sha256':sha(Path(__file__).with_name('PROTOCOL.json')),'source_design_sha256':sha(source/'design.json'),
            'script_sha256':sha(__file__),'case':case,'n_qubits':n,'n_electrons':6,'time':1.,'steps':100,
            'state_info':state_info,'spectrum':spectrum,'initial_states_file_sha256':sha(out/'initial_states.npz'),
            'hamiltonian_sha256':design['hamiltonian_sha256'],'ordered_schedule_sha256':design['schedule_sha256'],
            'propagation_space':'Full Hilbert space, no projection/postselection','checks_before_evolution':checks}
    write(out/'frozen_state_design.json',frozen)
    event(event='states_frozen',case=case,states=state_info,seconds=time.monotonic()-started)
    # Independent direct tensor-axis rotation checks on actual nontrivial input.
    chosen=np.unique(np.linspace(0,len(keys)-1,12,dtype=int));v=full['equal_coherent'].copy()
    direct=v.copy();terms=[];duration=.037
    for i in chosen:
        key=keys[i];angle=coeff[key].real*duration
        direct=np.cos(angle)*direct-1j*np.sin(angle)*tensor_axis_pauli(direct,key,n)
        terms+=ev.compile_ordered_terms([key],coeff,n,1e-12)
    native,_=ev.evolve_trotter_state(v,terms,1,1,duration,2**30)
    checks['12_rotation_tensor_axis_check_distance']=float(np.linalg.norm(direct-native))
    require(checks['12_rotation_tensor_axis_check_distance']<1e-11,'Pauli tensor-axis mismatch')
    # Tensor-axis Hamiltonian action independent of grouped sparse assembly.
    independent=np.zeros(dim,complex)
    for key,c in coeff.items():independent+=c*tensor_axis_pauli(v,key,n)
    checks['tensor_axis_H_action_distance']=float(np.linalg.norm(independent-H@v))
    require(checks['tensor_axis_H_action_distance']<1e-10,'Independent H action mismatch')
    checks['H_S2_commutator_on_coherent_norm']=float(np.linalg.norm(H@S2(v)-S2(H@v)))
    exact={};saved={f'initial_{name}':v for name,v in full.items()};reference_info={}
    trace=-1j*H.diagonal().sum()
    for name,v in full.items():
        t=time.monotonic();exact[name]=expm_multiply(-1j*H,v,traceA=trace)
        half=expm_multiply(-.5j*H,v,traceA=trace/2)
        split=expm_multiply(-.5j*H,half,traceA=trace/2)
        if name!='equal_coherent':analytic=np.exp(-1j*state_info[name]['energy_identity_free_hartree'])*v
        else:analytic=sum(np.exp(-1j*state_info[k]['energy_identity_free_hartree'])*full[k] for k in ('sector_ground','sector_next_distinct_energy'))/np.sqrt(2)
        info={**diagnostics(exact[name],H,S2,nm,sm),'split_time_reference_distance':float(np.linalg.norm(exact[name]-split)),
              'eigenphase_prediction_distance':float(np.linalg.norm(exact[name]-analytic)),
              'elapsed_seconds':time.monotonic()-t,'state_sha256':arrsha(exact[name])}
        require(max(info['split_time_reference_distance'],info['eigenphase_prediction_distance'])<5e-10,'Reference validation failed')
        reference_info[name]=info;saved[f'exact_{name}']=exact[name]
    write(out/'reference_checks.json',{'checks':checks,'references':reference_info})
    rows=[]
    for method in METHODS:
        for name,v in full.items():
            t=time.monotonic();actual,count=ev.evolve_trotter_state(v,compiled[method],1,100,1.,2**30)
            m=state_error_metrics(exact[name],actual);d=diagnostics(actual,H,S2,nm,sm)
            require(max(abs(m['exact_norm']-1),abs(m['approx_norm']-1))<1e-8,'Norm drift')
            require(count==100*len(keys),'Work mismatch')
            row={'case':case,'state':name,'ordering':method,'time':1.,'steps':100,'n_qubits':n,
                 'nominal_exponentials':count,'elapsed_seconds':time.monotonic()-t,'state_sha256':arrsha(actual),
                 'energy_drift_vs_exact_hartree':d['energy_identity_free_hartree']-reference_info[name]['energy_identity_free_hartree'],
                 **m,**d}
            rows.append(row);saved[f'{method}_{name}']=actual
            event(event='evolved',case=case,state=name,ordering=method,E=row['one_minus_overlap'],seconds=row['elapsed_seconds'])
    # Historical HF states are a comparator, explicitly not new propagations.
    archive=json.loads((source/'results.json').read_text());require(sha(source/'states.npz')==archive['states_sha256'],'Archived HF states hash mismatch')
    hfrows=[]
    with np.load(source/'states.npz',allow_pickle=False) as old:
        indices=old['affine_basis_indices'];hfexact=np.zeros(dim,complex);hfexact[indices]=old['exact']
        hfref=diagnostics(hfexact,H,S2,nm,sm)
        for method in METHODS:
            actual=np.zeros(dim,complex);actual[indices]=old[f'{method}_r100']
            matching=next(r for r in archive['rows'] if r['ordering']==method and r['steps']==100)
            m=state_error_metrics(hfexact,actual);d=diagnostics(actual,H,S2,nm,sm)
            require(abs(m['one_minus_overlap']-matching['one_minus_overlap'])<1e-20,'HF replay metric differs')
            hfrows.append({'case':case,'state':'HF_archived','ordering':method,'time':1.,'steps':100,
                           'source_results_sha256':sha(source/'results.json'), 'source_states_sha256':sha(source/'states.npz'),
                           'energy_drift_vs_exact_hartree':d['energy_identity_free_hartree']-hfref['energy_identity_free_hartree'],**m,**d})
    np.savez_compressed(out/'states.npz',**saved)
    result={'rows':rows,'archived_HF_comparator_rows':hfrows,'states_file_sha256':sha(out/'states.npz'),
            'state_arrays_sha256':{k:arrsha(v) for k,v in saved.items()},'frozen_state_design_sha256':sha(out/'frozen_state_design.json'),
            'reference_checks_sha256':sha(out/'reference_checks.json'),'elapsed_seconds':time.monotonic()-started,
            'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    require(sha(out/'initial_states.npz')==frozen['initial_states_file_sha256'],'Initial states changed')
    write(out/'results.json',result)
    event(event='complete',case=case,elapsed_seconds=result['elapsed_seconds'],peak_rss_KiB=result['peak_rss_KiB'])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    require(Path(ev.__file__).resolve().is_relative_to(ROOT.resolve()),'Unintended QHAT import')
    require(subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()=='0be7d4857537668eec957dabaa8438ea5c0d33f2','Wrong QHAT revision')
    protocol=json.loads(Path(__file__).with_name('PROTOCOL.json').read_text());require(args.case in protocol['case_plan'],'Case not predeclared')
    run(args.case,args.output)
