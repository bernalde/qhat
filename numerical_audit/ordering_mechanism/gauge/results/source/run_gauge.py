#!/usr/bin/env python3
"""Bounded fixed-pattern orbital-sign diagnostic; see PROTOCOL.md."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[k]='2'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/ordering-gauge-numba')
os.environ.setdefault('MPLCONFIGDIR','/tmp/ordering-gauge-mpl')
import argparse,hashlib,json,resource,shutil,subprocess,sys,time,importlib.metadata
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from openfermion import InteractionOperator

CASES=('H2O_s-1.00_sto-6g_as-006-004','H2O_s-1.00_hgbs-5_as-004-008',
       'H2O_s-1.00_hgbs-5_as-006-008','H2O_s-1.00_hgbs-5_as-006-010')
METHODS=('fermionic_signed','jw_signed','jw_magnitude')
PIN='0be7d4857537668eec957dabaa8438ea5c0d33f2'
ATOL=1e-12
def require(x,msg):
    if not x: raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def objsha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def arrsha(a):
    a=np.asarray(a)
    return objsha({'shape':a.shape,'dtype':a.dtype.str,'raw_sha256':hashlib.sha256(a.tobytes()).hexdigest()})
def dump(p,o):
    with Path(p).open('x') as f:json.dump(o,f,indent=2,allow_nan=False);f.write('\n')
def keyjson(k):return [[int(p),a] for p,a in k]
def canonical(co):return [[keyjson(k),complex(v).real.hex(),complex(v).imag.hex()] for k,v in sorted(co.items())]
def patterns(m):
    predicates=(lambda p:False,lambda p:p==0,lambda p:p==m-1,lambda p:p%2==1,
                lambda p:p<m//2,lambda p:p%3==1)
    ans=[[(-1 if f(p) else 1) for p in range(m)] for f in predicates]
    require(len({tuple(a) for a in ans})==6,'Non-distinct frozen patterns')
    return ans
def phases(states,spin):
    n=len(spin); mask=sum(1<<(n-1-p) for p,s in enumerate(spin) if s<0)
    return np.array([(-1 if (int(s)&mask).bit_count()%2 else 1) for s in states],dtype=np.int8)
def psign(k,spin):
    s=1
    for p,a in k:
        if a in ('X','Y'):s*=int(spin[p])
    return s
def block_record(owned,index):
    return [{'parent':int(b.parent_index),'paulis':[index[k] for k in b.pauli_keys]} for b in owned.owned_blocks]
def independent_action(k,n):
    """Direct local actions on basis columns, no QHAT descriptors or affine map."""
    col=np.arange(1<<n,dtype=np.int64);row=col.copy();v=np.ones(col.size,complex)
    for p,a in k:
        bit=(col>>(n-1-p))&1
        if a=='X':row^=1<<(n-1-p)
        elif a=='Y':row^=1<<(n-1-p);v*=1j*(1-2*bit)
        elif a=='Z':v*=1-2*bit
        else:raise RuntimeError(a)
    return row,v
def independent_product(initial,order,co,n,actions):
    state=initial.copy()
    for _ in range(100):
        for k in order:
            row,v=actions[k];a=float(complex(co[k]).real)/100
            applied=np.empty_like(state);applied[row]=v*state
            state=np.cos(a)*state-1j*np.sin(a)*applied
    return state
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--qhat-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    qhat=args.qhat_root.resolve();sys.path.insert(0,str(qhat.parent))
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as bl
    from qhat.analysis import fermionic_coefficients as fc
    from qhat.analysis import fermionic_ownership as fo
    from qhat.analysis import state_metrics as met
    from qhat.numerical_audit import reference_tools as ref
    mods=(ev,bl,fc,fo,met,ref)
    require(all(Path(m.__file__).resolve().is_relative_to(qhat) for m in mods),'Wrong QHAT import root')
    require(subprocess.check_output(['git','rev-parse','HEAD'],cwd=qhat,text=True).strip()==PIN,'Wrong QHAT revision')
    paths={str(Path(m.__file__).resolve().relative_to(qhat)):Path(m.__file__).resolve() for m in mods}
    paths['analysis/benchmark_L_sweep_trotter.py']=qhat/'analysis/benchmark_L_sweep_trotter.py'
    paths['run_gauge.py']=Path(__file__).resolve();paths['PROTOCOL.md']=Path(__file__).with_name('PROTOCOL.md')
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    protocol={'source_sha256':{k:sha(p) for k,p in paths.items()},'qhat_commit':PIN,'cases':CASES,
              'methods':METHODS,'T':1.,'r':100,'tensor_policy':'historical_policy','cutoff':ATOL,
              'pattern_description':'PROTOCOL.md','python':sys.version,
              'packages':{k:importlib.metadata.version(k) for k in ('numpy','scipy','numba','openfermion')},
              'planned_resorted_rows':72,'planned_fixed_control_rows':60,'max_qubits':16,'max_threads':2}
    dump(out/'protocol.json',protocol)
    for k,p in paths.items():
        dest=out/'source'/k;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    rows=[];checks=[]
    manifest=json.loads((qhat/'numerical_audit/h2o_recovery/data/manifest.json').read_text())
    for case in CASES:
        basis=case.split('_')[2];electrons=int(case.split('_as-')[1].split('-')[0])
        inp=qhat/'numerical_audit/h2o_recovery/data'/basis/'historical_policy'/f'{case}.tensors.npz'
        require(sha(inp)==manifest['bases'][basis]['files_sha256'][f'historical_policy/{inp.name}'],'Changed tensor input')
        old=qhat/'numerical_audit/followup_h2o/results'/f'{case}__historical_policy'
        oldresult=json.loads((old/'results.json').read_text());olddesign=json.loads((old/'design.json').read_text())
        require(sha(old/'states.npz')==oldresult['states_sha256'],'Old state archive hash mismatch')
        require(sha(old/'design.json')==oldresult['design_sha256'],'Old design hash mismatch')
        with np.load(old/'states.npz',allow_pickle=False) as a:oldarrays={k:a[k] for k in a.files}
        for k,a in oldarrays.items():require(arrsha(a)==oldresult['state_array_sha256'][k],'Old array hash mismatch')
        interaction,n=bl.load_interaction_operator(inp);hf=((1<<electrons)-1)<<(n-electrons)
        base=fo.build_molecular_ownership(fc.fermion_from_interaction(interaction),atol=ATOL,ordering='signed_ascending')
        baseco=base.final_coefficients;keys=list(baseco);index={k:i for i,k in enumerate(keys)}
        require(canonical(baseco)==olddesign['identity_free_coefficients'],'Original retained H differs from saved design')
        baseorders={'fermionic_signed':list(base.pauli_order),'jw_signed':bl.signed_coefficient_lexicographic_order(keys,baseco,n,ATOL),
                    'jw_magnitude':bl.magnitude_descending_lexicographic_order(keys,baseco,n,ATOL)}
        baseblocks=block_record(base,index)
        raw=ev.compile_ordered_terms(keys,baseco,n,ATOL);states,reduced,flipbasis=ref.affine_space(raw,hf)
        require(np.array_equal(states,oldarrays['affine_basis_indices']),'Original affine closure changed')
        H,_=ref.assemble_hamiltonian(reduced,len(states))
        initial=oldarrays['initial'];exact=oldarrays['exact'];nm=np.array([int(s).bit_count()==electrons for s in states])
        am=sum(1<<(n-1-p) for p in range(0,n,2));bm=sum(1<<(n-1-p) for p in range(1,n,2))
        sm=np.array([(int(s)&am).bit_count()==electrons//2 and (int(s)&bm).bit_count()==electrons//2 for s in states])
        casedir=out/case;casedir.mkdir();shutil.copyfile(inp,casedir/'original.tensors.npz')
        casecheck={'case':case,'n_qubits':n,'closed_dimension':len(states),'input_sha256':sha(inp),
                   'original_state_archive_sha256':sha(old/'states.npz'),'original_design_sha256':sha(old/'design.json')}
        if n==10:
            actions={k:independent_action(k,n) for k in keys};dense=np.zeros((1<<n,1<<n),complex);cols=np.arange(1<<n)
            for k,c in baseco.items():
                dst,v=actions[k];dense[dst,cols]+=c*v
            fullinitial=np.zeros(1<<n,complex);fullinitial[states]=initial
            vals,vec=eigh(dense,check_finite=True);denseexact=vec@(np.exp(-1j*vals)*(vec.conj().T@fullinitial))
            embed=np.zeros(1<<n,complex);embed[states]=exact
            casecheck['independent_dense_reference_distance']=float(np.linalg.norm(denseexact-embed))
            require(casecheck['independent_dense_reference_distance']<1e-10,'Independent dense reference mismatch')
            casecheck['independent_dense_H_distance']=float(np.linalg.norm(dense[np.ix_(states,states)]-H.toarray()))
            require(casecheck['independent_dense_H_distance']<1e-12,'Independent dense H mismatch')
        for gauge,spatial in enumerate(patterns(n//2)):
            gdir=casedir/f'gauge{gauge}';gdir.mkdir();spin=np.repeat(np.array(spatial,dtype=np.int8),2)
            one=interaction.one_body_tensor*spin[:,None]*spin[None,:]
            two=interaction.two_body_tensor*spin[:,None,None,None]*spin[None,:,None,None]*spin[None,None,:,None]*spin[None,None,None,:]
            shifted=InteractionOperator(interaction.constant,one,two)
            np.savez_compressed(gdir/'transformed.tensors.npz',constant=interaction.constant,one_body=one,two_body=two)
            owned=fo.build_molecular_ownership(fc.fermion_from_interaction(shifted),atol=ATOL,ordering='signed_ascending')
            co=owned.final_coefficients
            require(set(co)==set(baseco),'Gauge changed retained Pauli keys')
            require(all(co[k]==psign(k,spin)*baseco[k] for k in keys),'Gauge coefficients not exact signed baseline')
            require(owned.scalar==base.scalar,'Gauge scalar changed')
            require(not owned.fallback_pauli_keys and not owned.diagnostics['noncommuting_owned_block_count'],'Invalid ownership')
            orders={'fermionic_signed':list(owned.pauli_order),'jw_signed':bl.signed_coefficient_lexicographic_order(list(co),co,n,ATOL),
                    'jw_magnitude':bl.magnitude_descending_lexicographic_order(list(co),co,n,ATOL)}
            require(orders['jw_magnitude']==baseorders['jw_magnitude'],'Magnitude schedule not gauge invariant')
            rawg=ev.compile_ordered_terms(keys,co,n,ATOL);stateg,reducedg,basisg=ref.affine_space(rawg,hf)
            require(np.array_equal(stateg,states) and basisg==flipbasis,'Gauge changed affine closure')
            Hg,_=ref.assemble_hamiltonian(reducedg,len(states));d=phases(states,spin)
            conjugated=H.multiply(d[:,None]).multiply(d[None,:]).tocsr();difference=Hg-conjugated
            require(difference.nnz==0 or np.max(np.abs(difference.data))==0,'Sparse H conjugation not exact')
            # A deliberate wrong coefficient must be rejected by the exact gauge certificate.
            flipped=[k for k in keys if psign(k,spin)<0]
            negative_control=None
            if flipped:
                k=flipped[0];wrong=dict(co);wrong[k]=-wrong[k]
                negative_control=not all(wrong[k]==psign(k,spin)*baseco[k] for k in keys)
                require(negative_control,'Wrong-sign negative control was accepted')
            initg=d*initial;exactg=d*exact;saved={'affine_basis_indices':states,'gauge_diagonal':d,'initial':initg,'exact':exactg}
            blocks=block_record(owned,index);baseowner={p:b['parent'] for b in baseblocks for p in b['paulis']}
            newowner={p:b['parent'] for b in blocks for p in b['paulis']}
            design={'case':case,'gauge':gauge,'spatial_signs':spatial,'spin_signs':spin.tolist(),
                    'n_qubits':n,'n_electrons':electrons,'closed_dimension':len(states),'pauli_terms':len(keys),
                    'identity_free_coefficients':canonical(co),'pauli_key_reference':[keyjson(k) for k in keys],
                    'original_orders':{m:[index[k] for k in baseorders[m]] for m in METHODS},
                    'resorted_orders':{m:[index[k] for k in orders[m]] for m in METHODS},
                    'original_blocks':baseblocks,'resorted_blocks':blocks,
                    'parent_order':list(owned.parent_order),'original_parent_order':list(base.parent_order),
                    'owner_changed_paulis':sum(newowner[k]!=v for k,v in baseowner.items()),
                    'block_partition_changed':sorted(sorted(b['paulis']) for b in blocks)!=sorted(sorted(b['paulis']) for b in baseblocks),
                    'exact_retained_coefficient_covariance':True,'exact_sparse_H_covariance':True,
                    'wrong_sign_negative_control_rejected':negative_control,
                    'coefficient_sign_changes':len(flipped),'original_tensor_sha256':sha(inp),
                    'transformed_tensor_sha256':sha(gdir/'transformed.tensors.npz'),'ownership_diagnostics':owned.diagnostics}
            grows=[]
            for name in METHODS:
                modes=('resorted',) if gauge==0 else ('resorted','fixed_control')
                for mode in modes:
                    schedule=orders[name] if mode=='resorted' else baseorders[name]
                    compiled=[reducedg[index[k]] for k in schedule];start=time.monotonic()
                    actual,count=ev.evolve_trotter_state(initg,compiled,1,100,1.,2**30)
                    require(count==100*len(keys),'Wrong work count')
                    metric=met.state_error_metrics(exactg,actual);normalized=actual/metric['approx_norm']
                    row={'case':case,'gauge':gauge,'spatial_signs':spatial,'ordering':name,'mode':mode,
                         'T':1.,'r':100,'n_qubits':n,'elapsed_seconds':time.monotonic()-start,
                         'nominal_exponentials':count,'changed_order_positions':sum(a!=b for a,b in zip(schedule,baseorders[name])),
                         'schedule_sha256':objsha([index[k] for k in schedule]),'owner_changed_paulis':design['owner_changed_paulis'],
                         'particle_number_leakage':float(np.vdot(normalized[~nm],normalized[~nm]).real),
                         'spin_sector_leakage':float(np.vdot(normalized[~sm],normalized[~sm]).real),**metric}
                    expected=d*oldarrays[f'{name}_r100']
                    if gauge==0 or mode=='fixed_control' or name=='jw_magnitude':
                        row['matched_conjugated_state_distance']=float(np.linalg.norm(actual-expected))
                        require(row['matched_conjugated_state_distance']<1e-10,'Fixed/magnitude schedule covariance failed')
                    if n==10 and mode=='resorted':
                        fullg=np.zeros(1<<n,complex);fullg[states]=initg
                        independent=independent_product(fullg,schedule,co,n,actions)
                        embed=np.zeros(1<<n,complex);embed[states]=actual
                        row['independent_full_product_distance']=float(np.linalg.norm(independent-embed))
                        require(row['independent_full_product_distance']<1e-10,'Independent full product disagrees')
                    saved[f'{name}_{mode}']=actual;grows.append(row);rows.append(row)
                    print(json.dumps({'case':case,'gauge':gauge,'ordering':name,'mode':mode,'E':row['one_minus_overlap']}),flush=True)
            np.savez_compressed(gdir/'states.npz',**saved);dump(gdir/'design.json',design)
            dump(gdir/'results.json',{'rows':grows,'state_file_sha256':sha(gdir/'states.npz'),
                 'array_sha256':{k:arrsha(a) for k,a in saved.items()},'design_sha256':sha(gdir/'design.json')})
        checks.append(casecheck);dump(casedir/'checks.json',casecheck)
    require(len(rows)==132,'Wrong campaign row count')
    require(all(sha(p)==protocol['source_sha256'][k] for k,p in paths.items()),'Source changed mid-run')
    summary={'rows':rows,'checks':checks,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    require(summary['peak_rss_KiB']<8*1024**2,'Exceeded8GiB RSS')
    dump(out/'results.json',summary)
    filehashes={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    dump(out/'artifact_manifest.json',filehashes)
    print(json.dumps({'completed_rows':len(rows),'peak_rss_KiB':summary['peak_rss_KiB']}),flush=True)
if __name__=='__main__':main()
