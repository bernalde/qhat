"""Frozen correlated/coherent-state panel; phase separation forbids error tuning."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from response import Features,canonical_eigensystem,determinant as e

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'prospective_reference'/'evidence'/'results'
CASES=('Li2_r2p95','F2_r1p55','H2O_scale1p10')
TIMES=(.37,1.73)
BUDGETS=(10,30)
STATES=('ground','excited','coherent_real','coherent_imag')
ARMS=('hf_current','state_current','state_JD','state_pauli_JD','hf_pauli_large',
      'hf_symmetric','state_symmetric','jw_signed','jw_magnitude')
FLOOR=1e-12


def require(condition,message):
    if not condition:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    part=path.with_suffix(path.suffix+'.partial')
    part.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');part.replace(path)


def decoded(raw):return {tuple((int(q),p) for q,p in w):int(c) for w,c in raw}


def encode(block):return [[w,str(c)] for w,c in sorted(block.items())]


def sources():
    files=[HERE/p for p in ('PROTOCOL.md','response.py','study.py','verify.py')]
    files+=[HERE.parent/'prospective_reference'/'engine.py']
    return {str(p.relative_to(HERE.parent)):sha(p) for p in files}


def source_commit():
    return subprocess.check_output(['git','-C',str(HERE),'rev-parse','HEAD'],text=True).strip()


def matrix(model):
    n=model['n'];require(0<n<=8 and len(model['coefficients'])<=128,'Target dimension guard')
    dim=1<<n;basis=np.arange(dim);H=np.zeros((dim,dim),complex)
    for word,c in decoded(model['coefficients']).items():
        d=e.descriptor({word:c});require(d[0]<dim,'Invalid target mask')
        H[basis^d[0],basis]+=[e.value(d,int(v))/int(model['denominator']) for v in basis]
    require(np.max(abs(H-H.conj().T))<1e-12,'Target not Hermitian')
    return H


def sector(model):
    n=model['n'];require(0<n<=8,'State dimension guard')
    hf=(1<<model['electrons'])-1
    counts=lambda v,p:sum((v>>q)&1 for q in range(p,n,2))
    return [v for v in range(1<<n) if all(counts(v,p)==counts(hf,p) for p in (0,1))]


def natural_occupations(psi,n):
    require(0<n<=8 and len(psi)==1<<n,'RDM size guard')
    gamma=np.zeros((n,n),complex)
    for v,amplitude in enumerate(psi):
        for q in range(n):
            if not (v>>q)&1:continue
            w=v^(1<<q);sign=(-1)**((v&((1<<q)-1)).bit_count())
            for p in range(n):
                if not (w>>p)&1:
                    factor=sign*(-1)**((w&((1<<p)-1)).bit_count())
                    gamma[p,q]+=np.conj(psi[w^(1<<p)])*amplitude*factor
    occupations=np.linalg.eigvalsh(gamma)
    return occupations,float(np.sum(occupations*(1-occupations)))


def states(model):
    cpu=time.process_time();wall=time.monotonic()
    H=matrix(model);ids=sector(model)
    values,U=canonical_eigensystem(H[np.ix_(ids,ids)])
    g=np.zeros(len(H),complex);h=g.copy()
    g[ids]=U[:,0];h[ids]=U[:,1]
    vectors=np.column_stack([g,h,(g+h)/np.sqrt(2),(g+1j*h)/np.sqrt(2)])
    metadata=dict(sector=ids,energies=values.tolist(),gap=float(values[1]-values[0]),
        labels={},cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
    for at,label in enumerate(STATES):
        psi=vectors[:,at];weights=abs(psi)**2
        occupations,defect=natural_occupations(psi,model['n'])
        metadata['labels'][label]=dict(participation_ratio=float(1/np.sum(weights**2)),
            maximum_determinant_weight=float(max(weights)),norm=float(np.linalg.norm(psi)),
            determinant_like=bool(1/np.sum(weights**2)<=1+1e-8),
            natural_occupations=occupations.tolist(),one_RDM_idempotency_defect=defect,
            correlated_beyond_single_Slater=bool(defect>1e-8))
        if at<2:metadata['labels'][label]['eigen_residual']=float(np.linalg.norm(H@psi-values[at]*psi))
    metadata.update(cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
    return vectors,metadata


def timed(call):
    cpu=time.process_time();wall=time.monotonic();result=call()
    return result,dict(cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)


def construct(model,old,psi,T):
    blocks=[decoded(b) for b in model['blocks']];total=decoded(model['coefficients'])
    require(len(blocks)<=128 and len(total)<=128,'Ordering size guard')
    paulis=[{w:c} for w,c in sorted(total.items())];den=int(model['denominator']);n=model['n']
    diag=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]==0),
                key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
    off=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]),
               key=lambda i:(-sum(abs(c) for c in blocks[i].values()),e.descriptor(blocks[i])[0]))
    D={w:c for i in diag for w,c in blocks[i].items()}
    cons,cons_cost=timed(lambda:Features(blocks,D,psi,T,den,n))
    paul,paul_cost=timed(lambda:Features(paulis,D,psi,T,den,n))
    placement,placement_cost=timed(lambda:e.placement(blocks,diag,off,cons))
    starts=lambda bs:[list(range(len(bs))),sorted(range(len(bs)),key=lambda i:(-sum(abs(c) for c in bs[i].values()),i))]
    insertion,insert_cost=timed(lambda:e.insertion(starts(blocks),cons,4*len(blocks)*(len(blocks)-1)))
    pinsert,pinsert_cost=timed(lambda:e.insertion(starts(paulis),paul,4*len(paulis)*(len(paulis)-1)))
    result={}
    for name,method,cost in (('state_current',placement,placement_cost),('state_JD',insertion,insert_cost),
                            ('state_pauli_JD',pinsert,pinsert_cost)):
        result[name]={**method,'partition':'pauli' if name=='state_pauli_JD' else 'conserving',
                      'formula':1,'search_cost':cost,'training':'current state and T'}
    frozen=next(p for p in old['points'] if p['state']=='HF')['methods']
    for name,previous in (('hf_current','current'),('hf_pauli_large','pauli_JD_large'),
                          ('hf_symmetric','symmetric_JD'),('jw_signed','jw_signed'),('jw_magnitude','jw_magnitude')):
        prior=frozen[previous]
        result[name]=dict(order=prior['order'],partition=prior['partition'],formula=prior['formula'],
            training='frozen previous HF/T=1; no new search',evaluations=0,
            previous_method=previous,previous_schedule_sha256=prior['schedule_sha256'],
            historical_setup_cpu_seconds=prior['setup_cpu_seconds'],
            historical_search_cpu_seconds=prior['search_cpu_seconds'])
    result['state_symmetric']={**result['state_JD'],'formula':2,'reuses_order':'state_JD'}
    for name in ARMS:
        method=result[name];bs=paulis if method['partition']=='pauli' else blocks
        feature=paul if method['partition']=='pauli' else cons
        require(sorted(method['order'])==list(range(len(bs))),'Order is not a permutation')
        schedule=[(w,c) for i in method['order'] for w,c in sorted(bs[i].items())]
        require(dict(schedule)==total and len(schedule)==len(total),'Schedule changed target')
        method['score']=feature.score(feature.amplitude(method['order']))
        method['schedule']=encode_schedule(schedule)
        method['unit_lengths']=[len(bs[i]) for i in method['order']]
        if name.startswith('hf_') or name.startswith('jw_'):
            require(method['schedule']==frozen[method['previous_method']]['schedule'],'Frozen transfer order changed')
    return dict(methods=result,setup_costs={'conserving_features':cons_cost,'pauli_features':paul_cost},
        feature_storage_bytes={'conserving':cons.storage_bytes,'pauli':paul.storage_bytes})


def encode_schedule(schedule):return [[[[q,p] for q,p in w],str(c)] for w,c in schedule]


def bind(output):
    freeze=json.loads((output/'SCHEDULES_FROZEN.json').read_text())
    require(freeze['sources']==sources(),'Active source differs from frozen design')
    for path,h in freeze['artifacts'].items():require(sha(output/path)==h,'Frozen artifact changed: '+path)
    for path,h in freeze['inputs'].items():require(sha(OLD/path)==h,'Previous input/design changed')
    protocol=json.loads((output/'PROTOCOL_FROZEN.json').read_text())
    require(protocol['environment']=={p:importlib.metadata.version(p) for p in protocol['environment']},
            'Numerical environment changed')
    return freeze


def prepare(output):
    require(not output.exists(),'Use a fresh output directory; preserve partial attempts')
    for relative in sources():
        p=HERE.parent/relative
        source=subprocess.check_output(['git','-C',str(HERE),'show',source_commit()+':'+str(p.relative_to(HERE.parents[1]))])
        require(source==p.read_bytes(),'Commit active scientific source before freezing')
    output.mkdir(parents=True)
    protocol=dict(sources=sources(),source_commit=source_commit(),cases=CASES,times=TIMES,
        budgets=BUDGETS,states=STATES,arms=ARMS,python=sys.version,
        environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')},
        expected_schedules=216,expected_rows=432,protocol_sha256=sha(HERE/'PROTOCOL.md'))
    write(output/'PROTOCOL_FROZEN.json',protocol)
    prior_inputs=json.loads((OLD/'INPUTS_FROZEN.json').read_text())
    prior_designs=json.loads((OLD/'ALL_SCHEDULES_FROZEN.json').read_text())['design_sha256']
    inputs={};artifacts={'PROTOCOL_FROZEN.json':sha(output/'PROTOCOL_FROZEN.json')}
    for case in CASES:
        mp=OLD/'inputs'/case/'model.json';dp=OLD/'designs'/(case+'.json')
        require(sha(mp)==prior_inputs[case] and sha(dp)==prior_designs[case],'Previous evidence authentication failed')
        model=json.loads(mp.read_text());old=json.loads(dp.read_text())
        inputs[str(mp.relative_to(OLD))]=sha(mp);inputs[str(dp.relative_to(OLD))]=sha(dp)
        vectors,meta=states(model)
        statepath=output/'states'/(case+'.npz');statepath.parent.mkdir(exist_ok=True)
        np.savez_compressed(statepath,vectors=vectors)
        design=dict(case=case,state_file=str(statepath.relative_to(output)),state_sha256=sha(statepath),
            state_preparation=meta,points=[],input_sha256=sha(mp),previous_design_sha256=sha(dp))
        for at,label in enumerate(STATES):
            for T in TIMES:
                point=construct(model,old,vectors[:,at],T)
                point.update(state=label,T=T,state_index=at)
                design['points'].append(point)
                print('PREPARED',case,label,T,flush=True)
        path=output/'designs'/(case+'.json');write(path,design)
        artifacts[str(path.relative_to(output))]=sha(path)
        artifacts[str(statepath.relative_to(output))]=sha(statepath)
    write(output/'SCHEDULES_FROZEN.json',dict(sources=sources(),inputs=inputs,artifacts=artifacts,
        schedules=216,expected_rows=432))


def propagate(psi,method,model,T,steps):
    psi=psi.copy();basis=np.arange(len(psi));den=int(model['denominator'])
    descriptors=[e.descriptor({w:c}) for w,c in decoded_schedule(method['schedule'])]
    phases=[np.array([e.value(d,int(v))/den for v in basis]) for d in descriptors]
    ids,fractions=e.pulses(method['unit_lengths'],steps,method['formula']==2)
    for index,fraction in zip(ids,fractions):
        d=descriptors[index];theta=T*fraction*phases[index]
        psi=np.cos(theta)*psi-1j*np.sin(theta[basis^d[0]])*psi[basis^d[0]]
    return psi,len(ids)


def decoded_schedule(raw):return [(tuple((int(q),p) for q,p in w),int(c)) for w,c in raw]


def metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));phi=exact/np.linalg.norm(exact);psi=actual/norm
    overlap=np.vdot(phi,psi);residual=psi-phi*overlap
    infidelity=float(np.vdot(residual,residual).real)
    correct=np.array([i in set(sector(model)) for i in range(len(psi))])
    correct_number=np.array([i.bit_count()==model['electrons'] for i in range(len(psi))])
    return dict(error=infidelity/(1+abs(overlap)),infidelity=infidelity,norm=norm,
        number_leakage=float(np.vdot(psi[~correct_number],psi[~correct_number]).real),
        spin_projection_leakage=float(np.vdot(psi[~correct],psi[~correct]).real))


def run(output):
    bind(output);require(not (output/'rows').exists(),'Preserve existing physical attempts')
    (output/'rows').mkdir();(output/'propagated').mkdir()
    count=0
    for case in CASES:
        model=json.loads((OLD/'inputs'/case/'model.json').read_text())
        design=json.loads((output/'designs'/(case+'.json')).read_text())
        with np.load(output/design['state_file']) as data:vectors=data['vectors']
        (values,U),exact_setup=timed(lambda:np.linalg.eigh(matrix(model)))
        for at,point in enumerate(design['points']):
            initial=vectors[:,point['state_index']];T=point['T']
            exact=U@(np.exp(-1j*T*values)*(U.conj().T@initial))
            for name in ARMS:
                method=point['methods'][name]
                for multiplier in BUDGETS:
                    budget=multiplier*len(model['coefficients'])
                    steps=e.feasible_steps(method['unit_lengths'],budget,method['formula']==2)
                    (actual,rotations),cost=timed(lambda:propagate(initial,method,model,T,steps))
                    require(rotations<=budget and abs(np.linalg.norm(actual)-1)<1e-8,'Invalid physical circuit')
                    label=case+'.'+str(at)+'.'+name+'.b'+str(multiplier)
                    sp=output/'propagated'/(label+'.npz');np.savez_compressed(sp,actual=actual)
                    row=dict(case=case,state=point['state'],T=T,point=at,arm=name,budget_units=multiplier,
                        steps=steps,rotations=rotations,rotation_cap=budget,formula=method['formula'],
                        propagation_cost=cost,reference_setup_cost=exact_setup,**metrics(exact,actual,model),
                        state_file=str(sp.relative_to(output)),state_sha256=sha(sp),
                        design_sha256=sha(output/'designs'/(case+'.json')))
                    write(output/'rows'/(label+'.json'),row);count+=1
            print('EVOLVED',case,point['state'],T,flush=True)
    require(count==432,'Incomplete panel');write(output/'RUN_COMPLETE.json',dict(rows=count))


def summary(output):
    bind(output)
    verified=json.loads((output/'VERIFIED.json').read_text())
    require(verified['rows']==432 and verified['status']=='passed','Independent verification incomplete')
    require(verified['schedule_freeze_sha256']==sha(output/'SCHEDULES_FROZEN.json') and
            verified['sources']==sources(),'Stale independent verification')
    for path,h in verified['artifacts'].items():require(sha(output/path)==h,'Verified result changed')
    rows=[json.loads(p.read_text()) for p in sorted((output/'rows').glob('*.json'))]
    require(len(rows)==432,'Incomplete physical rows')
    lookup={(r['case'],r['state'],r['T'],r['arm'],r['budget_units']):r for r in rows}
    require(len(lookup)==432,'Duplicate rows')
    coverage=[];contrasts=[]
    for state in STATES:
        for T in TIMES:
            informative=[c for c in CASES if lookup[c,state,T,'hf_current',10]['error']>FLOOR]
            coverage.append(dict(state=state,T=T,informative_families=informative,
                sufficient=len(informative)>=2,status='covered' if len(informative)>=2 else 'inconclusive'))
    for case_group in [(c,) for c in CASES]+[CASES]:
        for state in STATES:
            for T in TIMES:
                for budget in BUDGETS:
                    for arm in ARMS[1:]:
                        pairs=[(lookup[c,state,T,arm,budget]['error'],lookup[c,state,T,'hf_current',budget]['error']) for c in case_group]
                        differences=[a-b for a,b in pairs];ratios=[a/b for a,b in pairs if b>FLOOR]
                        contrasts.append(dict(families=list(case_group),state=state,T=T,budget_units=budget,arm=arm,
                            wins=sum(d < -FLOOR for d in differences),ties=sum(abs(d)<=FLOOR for d in differences),
                            losses=sum(d>FLOOR for d in differences),floor_exclusions=len(pairs)-len(ratios),
                            median_ratio=float(np.median(ratios)) if ratios else None))
    write(output/'SUMMARY.json',dict(rows=432,coverage=coverage,contrasts=contrasts,
        all_strata_covered=all(c['sufficient'] for c in coverage),floor=FLOOR))


def preflight(output):
    require(not output.exists(),'Fresh preflight output required');output.mkdir(parents=True)
    start=time.monotonic();cpu=time.process_time()
    model=json.loads((OLD/'inputs'/'H2O_scale1p10'/'model.json').read_text())
    blocks=[decoded(b) for b in model['blocks']]
    # Shape-preserving synthetic coefficients, not any holdout ordering errors.
    altered=[];total={}
    for i,b in enumerate(blocks):
        factor=(i%7+1)/8
        b={w:int(c)*int(8*factor) for w,c in b.items()}
        altered.append(b);total.update(b)
    model['blocks']=[encode(b) for b in altered];model['coefficients']=encode(total)
    model['denominator']=str(int(model['denominator'])*8)
    vectors,meta=states(model)
    old=json.loads((OLD/'designs'/'H2O_scale1p10.json').read_text())
    # Frozen transfer orders retain word identity; coefficients now synthetic.
    old=json.loads(json.dumps(old))
    for p in old['points']:
        for m in p['methods'].values():m['schedule']=[[w,str(total[tuple((q,s) for q,s in w)])] for w,c in m['schedule']]
    point=construct(model,old,vectors[:,3],TIMES[-1]);method=point['methods']['state_pauli_JD']
    actual,rotations=propagate(vectors[:,3],method,model,TIMES[-1],30)
    from verify import independent_score,replay
    score=independent_score(model,method,vectors[:,3],TIMES[-1])
    require(abs(score-method['score'])<1e-10,'Preflight response mismatch')
    other=replay(model,method,vectors[:,3],TIMES[-1],30)
    require(np.max(abs(actual-other))<2e-10,'Preflight circuit mismatch')
    write(output/'PREFLIGHT.json',dict(status='passed',n=model['n'],words=len(total),
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-start,
        largest_pair_storage_bytes=point['feature_storage_bytes'],rotations=rotations,
        state_preparation=meta,source_sha256=sources(),scope='synthetic altered coefficients; no panel errors'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=('prepare','run','summary','preflight'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    globals()[a.phase](a.output.resolve())


if __name__=='__main__':main()
