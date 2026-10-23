"""Freeze, execute, and summarize the prospective three-geometry comparison."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
            'NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
CASES={
    'Li2_r2p95':dict(geometry=[['Li',[0.,0.,0.]],['Li',[0.,0.,2.95]]],electrons=2,n=6),
    'F2_r1p55':dict(geometry=[['F',[0.,0.,0.]],['F',[0.,0.,1.55]]],electrons=6,n=8),
    'H2O_scale1p10':dict(geometry=[['O',[0.,0.,0.]],['H',[.8332948775591,0.,.6457876002479]],
                                ['H',[-.8332948775591,0.,.6457876002479]]],electrons=4,n=8),
}
BUDGETS=(5,10,20,50)
ORIGINAL_ONE_PASS_SHA='3cfaced88f188d418ad5fac1be585478ad7894d6d4da53142489f82c8e82fcb6'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_sha(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.partial')
    temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(path)


def require(ok,message):
    if not ok:raise ValueError(message)


def encoded(block):
    return [[w,str(c)] for w,c in sorted(block.items())]


def decoded(raw):
    return {tuple((int(q),p) for q,p in w):int(c) for w,c in raw}


def bootstrap(root):
    sys.path.insert(0,str(root.parent))
    from qhat.numerical_audit.followup_classes import exact_blocks as eb
    from qhat.analysis import fermionic_coefficients as fc
    require(Path(eb.__file__).resolve().is_relative_to(root.resolve()),'Wrong QHAT import')
    return eb,fc


def sources(root):
    return dict(local={p.name:sha(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md','.txt')},
                qhat={name:sha(root/name) for name in (
                    'analysis/fermionic_coefficients.py','numerical_audit/followup_classes/exact_blocks.py')})


def freeze_protocol(a):
    require(not a.output.exists(),'Use a fresh output directory')
    require(sha(HERE/'one_pass.py')==ORIGINAL_ONE_PASS_SHA,'Original one-pass implementation changed')
    a.output.mkdir(parents=True)
    for case,spec in CASES.items():
        require(spec['n']<=8 and spec['n']%2==0 and spec['electrons']%2==0,
                'Unsupported dimensions before allocation')
    env={k:importlib.metadata.version(k) for k in
         ('numpy','scipy','numba','llvmlite','openfermion','openfermionpyscf','pyscf','h5py')}
    commit=subprocess.check_output(['git','-C',str(a.qhat_root),'rev-parse','HEAD'],text=True).strip()
    for name in sources(a.qhat_root)['local']:
        relative=str((HERE/name).relative_to(a.qhat_root))
        require(subprocess.check_output(['git','-C',str(a.qhat_root),'show',commit+':'+relative])
                ==(HERE/name).read_bytes(),'Commit scientific source before generating inputs')
    write(a.output/'PROTOCOL_FROZEN.json',dict(cases=CASES,budgets=BUDGETS,T=1.,
        source_sha256=sources(a.qhat_root),environment=env,
        source_commit=commit,python=sys.version,platform=platform.platform(),
        one_pass_original_sha256=sha(HERE/'one_pass.py')))


def check_protocol(a):
    f=json.loads((a.output/'PROTOCOL_FROZEN.json').read_text())
    require(f['source_sha256']==sources(a.qhat_root),'Scientific source changed after freeze')
    require(f['cases']==CASES and f['budgets']==list(BUDGETS),'Protocol changed')
    require(f['environment']=={k:importlib.metadata.version(k) for k in f['environment']},
            'Numerical environment changed after freeze')
    return f


def generate(a):
    import numpy as np
    from pyscf import gto,scf,lib
    from openfermionpyscf._run_pyscf import compute_integrals
    from openfermion.ops.representations import get_active_space_integrals
    from openfermion import InteractionOperator
    eb,fc=bootstrap(a.qhat_root)
    freeze_protocol(a);lib.num_threads(1)
    for case,spec in CASES.items():
        output=a.output/'inputs'/case;output.mkdir(parents=True)
        started=time.monotonic();cpu_started=time.process_time()
        mol=gto.M(atom=spec['geometry'],basis='sto-6g',unit='Angstrom',
                  charge=0,spin=0,symmetry=False,verbose=0,max_memory=1024)
        mean=scf.RHF(mol);mean.conv_tol=1e-11;mean.max_memory=1024
        energy=mean.kernel()
        require(mean.converged,'RHF did not converge')
        require(np.array_equal(mean.mo_occ,np.r_[np.full(mol.nelectron//2,2.),
                    np.zeros(len(mean.mo_occ)-mol.nelectron//2)]),'Unexpected RHF occupation')
        one,two=compute_integrals(mol,mean)
        frozen=(mol.nelectron-spec['electrons'])//2
        active=list(range(frozen,frozen+spec['n']//2))
        require(active[-1]<len(one),'Insufficient orbitals for active window')
        core,h,g=get_active_space_integrals(one,two,list(range(frozen)),active)
        n=spec['n'];occupied=spec['electrons'];h1=np.zeros((n,n));h2=np.zeros((n,n,n,n))
        for sigma in (0,1):
            h1[sigma::2,sigma::2]=h
            for tau in (0,1):h2[sigma::2,tau::2,tau::2,sigma::2]=g/2
        constant=float(mol.energy_nuc()+core)
        contraction=constant+sum(h1[p,p] for p in range(occupied))+sum(
            h2[p,q,q,p]-h2[p,q,p,q] for p in range(occupied) for q in range(occupied))
        require(abs(contraction-energy)<1e-8,'Active-space HF contraction mismatch')
        ortho=float(np.max(abs(mean.mo_coeff.T@mean.get_ovlp()@mean.mo_coeff-np.eye(len(one)))))
        require(ortho<1e-8,'Orbitals not S-orthonormal')
        np.savez_compressed(output/'tensors.npz',constant=constant,one_body=h1,two_body=h2)
        np.savez_compressed(output/'orbitals.npz',mo_coeff=mean.mo_coeff,mo_energy=mean.mo_energy,
                            mo_occ=mean.mo_occ,ao_overlap=mean.get_ovlp(),ao_hcore=mean.get_hcore())
        raw=fc.fermion_from_interaction(InteractionOperator(constant,h1,h2))
        parents,total,den,policy=eb.exact_parents(raw,1e-12)
        _,_,classes,xs=eb.decompose(parents,total)
        blocks=[]
        for b,x in zip(classes,xs):blocks.extend([{w:c} for w,c in sorted(b.items())] if x==0 else [b])
        require(all(eb.internally_commuting(b) and not eb.number_commutator(b) for b in blocks),
                'Completed unit does not conserve number/commute internally')
        retained=[]
        for parent in parents:
            keys=parent['components'];first=keys[0];weight=parent['weight']
            comp=[[first,str(weight)]]
            if len(keys)>1:
                partner,phase=eb.adjoint(first)
                require(partner==keys[1],'Parent reconstruction mismatch')
                comp.append([partner,str(phase*weight)])
            retained.append(dict(components=comp,removed_scalar=str(parent['scalar'])))
        specout=dict(case=case,**spec,denominator=str(den),coefficients=encoded(total),
            blocks=[encoded(b) for b in blocks],fermionic_parents=retained,policy=policy,
            frozen_spatial=list(range(frozen)),active_spatial=active,
            hf_energy=float(energy),hf_contraction=float(contraction),orthogonality_error=ortho,
            scf_converged=True,scf_cycles=int(mean.cycles),generation_seconds=time.monotonic()-started,
            generation_cpu_seconds=time.process_time()-cpu_started,total_electrons=int(mol.nelectron),
            tensors_sha256=sha(output/'tensors.npz'),orbitals_sha256=sha(output/'orbitals.npz'))
        write(output/'model.json',specout)
        print('GENERATED',case,'qubits',n,'words',len(total),'units',len(blocks),flush=True)
    write(a.output/'INPUTS_FROZEN.json',{case:sha(a.output/'inputs'/case/'model.json') for case in CASES})


def current_placement(blocks,diag,off,v,den):
    import engine as e
    import one_pass
    import numpy as np
    cpu=time.process_time();began=time.monotonic()
    desc=[e.descriptor(b) for b in blocks]
    masks=[desc[i][0] for i in off]
    off_a={}
    for at,i in enumerate(off):
        for j in off[at+1:]:
            x,y=desc[i][0],desc[j][0]
            b=e.value(desc[j],v)*e.value(desc[i],v^y)-e.value(desc[i],v)*e.value(desc[j],v^x)
            d=x^y
            if d and b:off_a[d]=off_a.get(d,0)+b
    rows=[[e.value(desc[j],v)*(e.value(desc[i],v^desc[j][0])-e.value(desc[i],v))
           for j in off] for i in diag]
    D={w:c for i in diag for w,c in blocks[i].items()};dd=e.descriptor(D)
    destinations=set(off_a)|set(masks)
    rho={d:float(np.sinc((e.value(dd,v^d)-e.value(dd,v))/(2*np.pi*den))**2)
         for d in destinations}
    ratios={d:w.as_integer_ratio() for d,w in rho.items()}
    common=max((b for a,b in ratios.values()),default=1)
    weights={d:a*(common//b) for d,(a,b) in ratios.items()}
    setup_cpu=time.process_time()-cpu;setup_wall=time.monotonic()-began
    cpu=time.process_time();began=time.monotonic()
    scores=[]
    for slot in range(len(off)+1):
        a=dict(off_a)
        for row in rows:
            for j,(d,b) in enumerate(zip(masks,row)):a[d]=a.get(d,0)+(b if j>=slot else -b)
        scores.append(one_pass.score(a,weights))
    start=min(range(len(scores)),key=lambda s:(scores[s],s))
    slots,amplitude,history=one_pass.place(off_a,rows,masks,start,weights)
    order=[]
    for slot in range(len(off)+1):
        order.extend(i for i,s in zip(diag,slots) if s==slot)
        if slot<len(off):order.append(off[slot])
    score=.25*sum(rho[d]*(b/(den*den))**2 for d,b in amplitude.items())
    return dict(order=order,score=score,evaluations=(len(diag)+1)*(len(off)+1),
        start_slot=start,slots=slots,exact_score_numerator=str(one_pass.score(amplitude,weights)),
        history=history,initial_slot_scores=[str(s) for s in scores],
        exact_weight_denominator=str(common),setup_cpu_seconds=setup_cpu,setup_seconds=setup_wall,
        search_cpu_seconds=time.process_time()-cpu,seconds=time.monotonic()-began)


def fresh_phase(output,phase):
    """Preserve every frozen/partial attempt instead of overwriting evidence."""
    names=('designs','ALL_SCHEDULES_FROZEN.json','references','rows','states','run_summary.json','reference_timings.json',
           'verification','verification_summary.json') if phase=='prepare' else (
           'references','rows','states','run_summary.json','reference_timings.json','verification','verification_summary.json')
    if any((output/name).exists() for name in names):
        raise FileExistsError('Phase artifacts already exist; preserve this attempt and use fresh output')


def prepare(a):
    fresh_phase(a.output,'prepare')
    import numpy as np
    import engine as e
    check_protocol(a);inputs=json.loads((a.output/'INPUTS_FROZEN.json').read_text());designs={}
    for case in CASES:
        modelpath=a.output/'inputs'/case/'model.json';require(sha(modelpath)==inputs[case],'Input changed')
        d=json.loads(modelpath.read_text());n=d['n'];require(n<=8,'Dimension guard')
        blocks=[decoded(b) for b in d['blocks']];den=int(d['denominator']);hf=(1<<d['electrons'])-1
        diag=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]==0),
                    key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
        off=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]),
                   key=lambda i:(-sum(abs(c) for c in blocks[i].values()),e.descriptor(blocks[i])[0]))
        D={w:c for i in diag for w,c in blocks[i].items()}
        count=lambda v,parity:sum((v>>q)&1 for q in range(parity,n,2))
        sector=[v for v in range(1<<n) if all(count(v,p)==count(hf,p) for p in (0,1))]
        near=min((v for v in sector if v!=hf),key=lambda v:((v^hf).bit_count(),v))
        total=decoded(d['coefficients']);words=sorted(total);paulis=[{w:total[w]} for w in words]
        points=[]
        for label,v in (('HF',hf),('near',near)):
            methods={};current=current_placement(blocks,diag,off,v,den);budget=current['evaluations']
            methods['current']={**current,'partition':'conserving','formula':1}
            cpu=time.process_time();wall=time.monotonic()
            partner,coupling=e.choose_partner(blocks,v)
            corrected_features=e.Features(blocks,e.Reference(D,v,1.,den,partner,coupling))
            setup_cpu=time.process_time()-cpu;setup_seconds=time.monotonic()-wall
            cpu=time.process_time();corrected=e.placement(blocks,diag,off,corrected_features)
            methods['pair_reference']={**corrected,'partition':'conserving','formula':1,
                'partner':partner,'coupling_numerator':str(coupling),'setup_cpu_seconds':setup_cpu,
                'setup_seconds':setup_seconds,'search_cpu_seconds':time.process_time()-cpu}
            for name,bs,unweighted,big in (
                ('optimized_JD',blocks,False,False),('optimized_J0',blocks,True,False),
                ('pauli_JD',paulis,False,False),('pauli_JD_large',paulis,False,True)):
                cpu=time.process_time();wall=time.monotonic()
                reference=e.Reference({():0} if unweighted else D,v,1.,den)
                features=e.Features(bs,reference)
                starts=[list(range(len(bs))),sorted(range(len(bs)),key=lambda i:(-sum(abs(c) for c in bs[i].values()),i))]
                setup_cpu=time.process_time()-cpu;setup_seconds=time.monotonic()-wall
                cap=20*len(bs)*(len(bs)-1) if big else budget
                cpu=time.process_time();result=e.insertion(starts,features,cap)
                methods[name]={**result,'partition':'pauli' if bs is paulis else 'conserving','formula':1,
                    'setup_cpu_seconds':setup_cpu,'setup_seconds':setup_seconds,
                    'search_cpu_seconds':time.process_time()-cpu}
            for name,key in (('jw_signed',lambda w:(total[w],w)),('jw_magnitude',lambda w:(-abs(total[w]),w))):
                cpu=time.process_time();wall=time.monotonic()
                order=[words.index(w) for w in sorted(words,key=key)]
                methods[name]=dict(partition='pauli',formula=1,order=order,
                    evaluations=0,setup_cpu_seconds=time.process_time()-cpu,
                    setup_seconds=time.monotonic()-wall,search_cpu_seconds=0.)
            methods['symmetric_JD']={**methods['optimized_JD'],'formula':2,'reuses_order':'optimized_JD'}
            for name,method in methods.items():
                bs=paulis if method['partition']=='pauli' else blocks
                schedule=[(w,c) for i in method['order'] for w,c in sorted(bs[i].items())]
                require(dict(schedule)==total and len(schedule)==len(total),'Schedule changed target')
                method['schedule']=[[w,str(c)] for w,c in schedule]
                method['schedule_sha256']=object_sha(method['schedule'])
                method['unit_lengths']=[len(bs[i]) for i in method['order']]
            points.append(dict(state=label,occupation=v,methods=methods,position_budget=budget))
            print('PREPARED',case,label,'matched positions',budget,'larger Pauli',methods['pauli_JD_large']['evaluations'],flush=True)
        design=dict(case=case,input_sha256=sha(modelpath),points=points,n=n,electrons=d['electrons'],
                    denominator=d['denominator'],coefficients=d['coefficients'],blocks=d['blocks'])
        path=a.output/'designs'/(case+'.json');write(path,design);designs[case]=sha(path)
    write(a.output/'ALL_SCHEDULES_FROZEN.json',dict(source_sha256=sources(a.qhat_root),
        protocol_sha256=sha(a.output/'PROTOCOL_FROZEN.json'),input_sha256=inputs,design_sha256=designs,
        expected_rows=216,arms=9,states_per_geometry=2,budgets=BUDGETS))


def frozen(a):
    check_protocol(a);f=json.loads((a.output/'ALL_SCHEDULES_FROZEN.json').read_text())
    require(f['source_sha256']==sources(a.qhat_root),'Frozen source changed')
    require(f['protocol_sha256']==sha(a.output/'PROTOCOL_FROZEN.json'),'Protocol changed')
    for case,h in f['design_sha256'].items():
        require(sha(a.output/'designs'/(case+'.json'))==h,'Design changed')
        require(sha(a.output/'inputs'/case/'model.json')==f['input_sha256'][case],'Input changed')
    return f


def target_matrix(d):
    import numpy as np
    import engine as e
    n=d['n'];require(n<=8,'Dimension guard before matrix allocation')
    matrix=np.zeros((1<<n,1<<n),dtype=complex);den=int(d['denominator'])
    for word,c in decoded(d['coefficients']).items():
        desc=e.descriptor({word:c})
        for v in range(1<<n):matrix[v^desc[0],v]+=e.value(desc,v)/den
    return matrix


def run(a):
    fresh_phase(a.output,'run')
    import numpy as np
    from scipy.linalg import eigh
    from numba import njit
    import engine as e
    f=frozen(a)
    @njit(cache=True)
    def propagate(initial,ids,fractions,coefficients,xs,phases):
        actual=initial.copy();work=np.empty_like(actual)
        for term,fraction in zip(ids,fractions):
            theta=coefficients[term]*fraction
            cs=np.cos(theta);sn=-1j*np.sin(theta);x=xs[term]
            for v in range(len(actual)):work[v^x]=cs*actual[v^x]+sn*phases[term,v]*actual[v]
            actual,work=work,actual
        return actual
    # Compile before physical timing. This warmup is not a molecular trial.
    began=time.monotonic();cpu=time.process_time()
    propagate(np.array([1.+0j,0j]),np.array([0],np.int32),np.array([.1]),np.array([1.]),
              np.array([1],np.int32),np.ones((1,2)))
    warmup=dict(seconds=time.monotonic()-began,cpu_seconds=time.process_time()-cpu)
    rows=[];reference_timings={}
    for case in CASES:
        d=json.loads((a.output/'designs'/(case+'.json')).read_text());n=d['n'];require(n<=8,'Dimension guard')
        began=time.monotonic();cpu=time.process_time()
        H=target_matrix(d);eigenvalues,U=eigh(H);basis=np.arange(1<<n)
        reference_timings[case]=dict(target_and_eigensystem_seconds=time.monotonic()-began,
                                    target_and_eigensystem_cpu_seconds=time.process_time()-cpu,points={})
        for point in d['points']:
            initial=np.zeros(1<<n,dtype=complex);initial[point['occupation']]=1
            began=time.monotonic();cpu=time.process_time()
            exact=U@(np.exp(-1j*eigenvalues)*(U.conj().T@initial))
            reference_timings[case]['points'][point['state']]=dict(
                exact_action_seconds=time.monotonic()-began,exact_action_cpu_seconds=time.process_time()-cpu)
            prefix=case+'.'+point['state'];reference=a.output/'references'/(prefix+'.npz')
            reference.parent.mkdir(exist_ok=True);np.savez_compressed(reference,initial=initial,exact=exact)
            number=np.array([int(v).bit_count()==d['electrons'] for v in basis])
            counts=lambda v,p:sum((v>>q)&1 for q in range(p,n,2))
            spin=np.array([all(counts(int(v),p)==counts(point['occupation'],p) for p in (0,1)) for v in basis])
            for name,method in point['methods'].items():
                schedule=[(tuple((q,p) for q,p in w),int(c)) for w,c in method['schedule']]
                descriptors=[e.descriptor({w:c}) for w,c in schedule]
                coefficients=np.array([c/int(d['denominator']) for w,c in schedule])
                xs=np.array([desc[0] for desc in descriptors],dtype=np.int32)
                phases=np.array([[e.value(e.descriptor({w:1}),v) for v in basis] for w,c in schedule],dtype=float)
                for budget_units in BUDGETS:
                    budget=budget_units*len(schedule);symmetric=method['formula']==2
                    steps=e.feasible_steps(method['unit_lengths'],budget,symmetric)
                    require(steps>=1,'Unexpected infeasible budget')
                    ids,fractions=e.pulses(method['unit_lengths'],steps,symmetric)
                    require(len(ids)<=budget,'Exceeded physical rotation budget')
                    began=time.monotonic();actual=propagate(initial,ids,fractions,coefficients,xs,phases)
                    elapsed=time.monotonic()-began
                    require(abs(np.linalg.norm(actual)-1)<1e-8,'Norm drift')
                    phi=exact/np.linalg.norm(exact);psi=actual/np.linalg.norm(actual)
                    overlap=np.vdot(phi,psi);residual=psi-phi*overlap
                    infidelity=float(np.vdot(residual,residual).real)
                    error=infidelity/(1+abs(overlap))
                    label=prefix+'.'+name+'.b'+str(budget_units)
                    statepath=a.output/'states'/(label+'.npz');statepath.parent.mkdir(exist_ok=True)
                    np.savez_compressed(statepath,actual=actual)
                    row=dict(case=case,state=point['state'],arm=name,budget_units=budget_units,
                        steps=steps,formula=method['formula'],rotations=len(ids),rotation_budget=budget,
                        one_minus_overlap=error,infidelity=infidelity,norm=float(np.linalg.norm(actual)),
                        number_leakage=float(np.vdot(psi[~number],psi[~number]).real),
                        spin_projection_leakage=float(np.vdot(psi[~spin],psi[~spin]).real),
                        propagation_seconds=elapsed,schedule_sha256=method['schedule_sha256'],
                        state_file=str(statepath.relative_to(a.output)),state_sha256=sha(statepath),
                        reference_sha256=sha(reference),design_sha256=f['design_sha256'][case])
                    write(a.output/'rows'/(label+'.json'),row);rows.append(row)
                print('EVOLVED',prefix,name,flush=True)
    require(len(rows)==216,'Incomplete physical grid')
    write(a.output/'reference_timings.json',reference_timings)
    write(a.output/'run_summary.json',dict(rows=len(rows),synthetic_JIT_warmup=warmup))


def summarize(a):
    from decisions import compute
    f=frozen(a);verification=json.loads((a.output/'verification_summary.json').read_text())
    require(verification['status']=='passed' and verification['new_replays']==216,'Incomplete independent verification')
    rows=[];contexts={}
    for case in CASES:
        audit=json.loads((a.output/'verification'/(case+'.json')).read_text())
        require(audit['status']=='passed' and audit['design_sha256']==f['design_sha256'][case],'Stale case audit')
        for item in audit['rows']:
            path=a.output/'rows'/(item['row']+'.json')
            require(sha(path)==item['row_sha256'],'Verified row changed')
            row=json.loads(path.read_text())
            require(sha(a.output/row['state_file'])==row['state_sha256']==item['state_sha256'],'Verified state changed')
            rows.append(row)
        design=json.loads((a.output/'designs'/(case+'.json')).read_text())
        for point in design['points']:contexts[case,point['state']]=point['methods']
    result=compute(rows,contexts)
    independent=json.loads((a.output/'independent_decision.json').read_text())
    require(result==independent,'Independent decision-rule/contrast verification mismatch')
    write(a.output/'summary.json',dict(status='passed',rows=len(rows),decision=result,
        verification=verification,scientific_source_commit=json.loads((a.output/'PROTOCOL_FROZEN.json').read_text())['source_commit']))
    print(json.dumps({k:result[k] for k in ('primary','conditions','correction_passes',
        'median_construction_CPU_ratio','maximum_strict_loss_ratio')},indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=('generate','prepare','run','summarize'))
    p.add_argument('--qhat-root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();bootstrap(a.qhat_root)
    {'generate':generate,'prepare':prepare,'run':run,'summarize':summarize}[a.phase](a)


if __name__=='__main__':main()
