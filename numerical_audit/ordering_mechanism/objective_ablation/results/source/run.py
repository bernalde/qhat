#!/usr/bin/env python3
"""Fixed objective/partition ablations; see PROTOCOL.md."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/ordering-ablation-numba')
os.environ.setdefault('MPLCONFIGDIR','/tmp/ordering-ablation-mpl')
import argparse,hashlib,importlib.util,json,resource,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
def require(b,msg):
    if not b:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def encoded(b):return [[[[int(q),p] for q,p in word],str(value)] for word,value in sorted(b.items())]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent-run',type=Path,required=True);ap.add_argument('--phase2-run',type=Path,required=True)
    ap.add_argument('--qhat-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    qhat=a.qhat_root.resolve();sys.path.insert(0,str(qhat.parent));here=Path(__file__).resolve().parent
    spec=importlib.util.spec_from_file_location('phase2_reference',here/'phase2_reference.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    require(all(Path(m.__file__).resolve().is_relative_to(qhat) for m in (helper.ev,helper.ref)),'Unintended QHAT imports')
    require(subprocess.check_output(['git','rev-parse','HEAD'],cwd=qhat,text=True).strip()=='0be7d4857537668eec957dabaa8438ea5c0d33f2','Wrong QHAT commit')
    require(sha(here/'phase2_reference.py')==sha(a.phase2_run/'run.py'),'Helper is not frozen phase2 source')
    sources={p.name:p for p in (Path(__file__).resolve(),here/'phase2_reference.py',here/'PROTOCOL.md')}
    for rel in ('analysis/benchmark_b2_active_spaces_matrix_free.py','analysis/benchmark_b2_signed_coefficient_baseline.py','analysis/benchmark_L_sweep_trotter.py','analysis/state_metrics.py','numerical_audit/reference_tools.py','numerical_audit/run_frozen_audit.py'):
        sources[rel]=qhat/rel
    hashes={k:sha(p) for k,p in sources.items()};cases=json.loads((a.parent_run/'manifest.json').read_text())['cases'];a.output.mkdir(parents=True,exist_ok=False)
    for rel,p in sources.items():
        dst=a.output/'source'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
    manifest={'cases':cases,'source_sha256':hashes,'parent_manifest_sha256':sha(a.parent_run/'manifest.json'),
              'phase2_manifest_sha256':sha(a.phase2_run/'manifest.json'),'qhat_commit':'0be7d4857537668eec957dabaa8438ea5c0d33f2',
              'methods':['class_refined_J0','class_pauli_JD'],'expected_rows':16,'T':1.,'r':100,'threads':1,
              'status':'exploratory objective/partition ablation; not held-out'}
    write(a.output/'manifest.json',manifest);prepared=[];started=time.monotonic()
    for mode in manifest['methods']:
        arm=a.output/mode;arm.mkdir();write(arm/'manifest.json',{'cases':cases,'source_manifest_sha256':sha(a.output/'manifest.json')})
        for case in cases:
            old=json.loads((a.parent_run/f'{case}.design.json').read_text());phase2=json.loads((a.phase2_run/f'{case}.design.json').read_text())
            require(old['coefficients']==phase2['coefficients'] and old['denominator']==phase2['denominator'],'Changed phase2 target')
            total=helper.block(old['coefficients']);den=int(old['denominator']);hf=(1<<old['electrons'])-1
            if mode=='class_refined_J0':bs=[helper.block(b) for b in phase2['blocks']['class']]
            else:bs=[{k:total[k]} for k in sorted(total,key=lambda k:(helper.masks(k)[0],k))]
            xs=[helper.masks(next(iter(b)))[0] for b in bs];seen={}
            for b,x in zip(bs,xs):
                for k,v in b.items():
                    require(k not in seen,'Repeated emitted Pauli');seen[k]=v
                    require(helper.masks(k)[0]==x and helper.masks(k)[2]%2==0,'Non-real or mixed-mask unit')
            require(seen==total,'Partition changes exact target')
            pairs,rho=helper.pair_setup(bs,xs,hf,den)
            require(all(dest!=0 for dest,value in pairs.values()),'Nonzero HF expectation requires variance subtraction')
            if mode=='class_refined_J0':rho={k:1. for k in rho}
            lex=list(range(len(bs)));mag=sorted(lex,key=lambda i:(-sum(abs(v) for v in bs[i].values()),xs[i]))
            candidates=[helper.optimize(o,pairs,rho,den) for o in (lex,mag)];winner=candidates[0]
            for c in candidates[1:]:
                tol=max(1e-30,1e-12*min(winner[1],c[1]))
                if c[1]<winner[1]-tol or (abs(c[1]-winner[1])<=tol and tuple(c[0])<tuple(winner[0])):winner=c
            order,J,history=winner;require(sorted(order)==lex,'Candidate not permutation')
            d={k:old[k] for k in ('case','n_qubits','electrons','denominator','coefficients','affine_basis_indices','tensor_sha256','target_sha256','initial_sha256','coefficient_policy','product_convention')}
            d.update(blocks={'parent':old['blocks']['parent'],'class':[encoded(b) for b in bs]},class_flip_masks=xs,orders={mode:order},
                     planned_runs=[{'ordering':mode,'steps':100,'role':'ablation','T':1.}],
                     auxiliary_parent_blocks_purpose='Only to reconstruct the original common Pauli-invariant closure; not emitted',
                     unit_partition='individual diagonal Paulis and complete nonzero masks' if mode=='class_refined_J0' else 'individual Pauli words',
                     unit_number_conservation_expected=mode=='class_refined_J0',objective='J0' if mode=='class_refined_J0' else 'JD',
                     source_parent_design_sha256=sha(a.parent_run/f'{case}.design.json'),source_phase2_design_sha256=sha(a.phase2_run/f'{case}.design.json'))
            candidate={'case':case,'method':mode,'J':J,'order':order,'accepted_moves':history,
                       'starts':[{'order':c[0],'J':c[1],'accepted_moves':len(c[2])} for c in candidates],
                       'units':len(bs),'nonzero_pair_count':len(pairs),'duplicate_nonzero_masks':len([x for x in xs if x])-len(set(x for x in xs if x)),
                       'zero_destination_pair_amplitudes_exactly_zero':True,'preparation_seconds':time.monotonic()-started}
            write(arm/f'{case}.candidate.json',candidate);write(arm/f'{case}.design.json',d)
            prepared.append((mode,case,d,bs,J))
            print('FROZEN',mode,case,'J',J,'units',len(bs),'seconds',time.monotonic()-started,flush=True)
    # Selection is finished for every method/case before loading comparator errors.
    write(a.output/'all_candidates_frozen.json',{'count':len(prepared),'elapsed_seconds':time.monotonic()-started})
    baseline=json.loads((a.parent_run/'results.json').read_text());prior=json.loads((a.phase2_run/'results.json').read_text());rows=[]
    for mode,case,d,bs,J in prepared:
        arm=a.output/mode;den=int(d['denominator']);n=d['n_qubits'];total=helper.block(d['coefficients'])
        keys=sorted(set(total).union(*(helper.block(b) for b in d['blocks']['parent'])))
        raw=helper.ev.compile_ordered_terms(keys,{k:1. for k in keys},n,0.);basis,reduced,_=helper.ref.affine_space(raw,d['affine_basis_indices'][0])
        require(list(map(int,basis))==d['affine_basis_indices'],'Changed affine closure');units=dict(zip(keys,reduced))
        terms=[(v/den*units[k][0],units[k][1],units[k][2],units[k][3]) for i in d['orders'][mode] for k,v in bs[i].items()]
        require(len(terms)==len(total),'Changed exponential count')
        with np.load(a.parent_run/f'{case}.states.npz',allow_pickle=False) as f:initial=f['initial'];exact=f['exact'];dense=f['exact_dense']
        nm=np.array([int(v).bit_count()==d['electrons'] for v in basis]);sector=set(map(int,helper.ev.spin_sector_basis_indices(n,d['electrons'])));sm=np.array([int(v) in sector for v in basis])
        t=time.monotonic();actual=helper.evolve(initial,terms,1.,100)
        row={'case':case,'ordering':mode,'T':1.,'steps':100,'role':'ablation','pauli_rotations':100*len(terms),
             'surrogate_J':J,'seconds':time.monotonic()-t,**helper.physical_metrics(exact,actual,nm,sm)}
        for name in ('class_lex','owned_signed','jw_signed','jw_magnitude'):
            old=next(r for r in baseline if r['case']==case and r['ordering']==name and r['steps']==100)
            row[name+'_E']=old['one_minus_overlap'];row['ratio_over_'+name]=row['one_minus_overlap']/old['one_minus_overlap']
        fixed=next(r for r in prior if r['case']==case and r['ordering']=='class_refined_response' and r['steps']==100)
        row['fixed_refined_JD_E']=fixed['one_minus_overlap'];row['ratio_over_fixed_refined_JD']=row['one_minus_overlap']/fixed['one_minus_overlap']
        write(arm/f'{case}.{mode}.r100.ablation.json',row);rows.append(row)
        np.savez_compressed(arm/f'{case}.states.npz',initial=initial,exact=exact,exact_dense=dense,**{mode+'__r100__ablation':actual})
        print('EVOLVE',mode,case,row['one_minus_overlap'],'overJD',row['ratio_over_fixed_refined_JD'],flush=True)
    write(a.output/'results.json',rows)
    summary={}
    for mode in manifest['methods']:
        part=[r for r in rows if r['ordering']==mode];report={}
        for b in ('jw_signed','jw_magnitude','fixed_refined_JD'):
            outcomes={'win':0,'tie':0,'loss':0}
            for r in part:
                delta=r['one_minus_overlap']-r[b+'_E'];outcomes['tie' if abs(delta)<=1e-12 else ('win' if delta<0 else 'loss')]+=1
            report[b]={'outcomes':outcomes,'median_ratio':float(np.median([r['ratio_over_'+b] for r in part]))}
        summary[mode]=report
    write(a.output/'summary.json',summary)
    require(hashes=={k:sha(p) for k,p in sources.items()},'Sources changed during run')
    write(a.output/'completion.json',{'rows':len(rows),'seconds':time.monotonic()-started,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    write(a.output/'artifact_manifest.json',{str(p.relative_to(a.output)):sha(p) for p in sorted(a.output.rglob('*')) if p.is_file()})
    print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
