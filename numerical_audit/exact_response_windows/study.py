"""Freeze exact-window domains, certify them, then measure separate circuits."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from fractions import Fraction
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from exact import OLD,e,decode,require,freeze_domain,solve,validate_model

HERE=Path(__file__).resolve().parent
CASES=('Li2_r2p95','F2_r1p55','H2O_scale1p10')
BUDGETS=(10,30)
ARMS=(('heuristic',1),('certified_window',1),('current',1),('pauli_JD_large',1),('heuristic',2),('certified_window',2))
SOURCE_FILES=('exact.py','study.py','verify.py','PROTOCOL.md')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def sources():
    return {str(p.relative_to(HERE.parent)):sha(p) for p in [HERE/x for x in SOURCE_FILES]+[OLD/'engine.py']}


def old_case(case):
    prior=OLD/'evidence/results';freeze=json.loads((prior/'ALL_SCHEDULES_FROZEN.json').read_text())
    mp=prior/'inputs'/case/'model.json';dp=prior/'designs'/(case+'.json')
    require(sha(mp)==freeze['input_sha256'][case] and sha(dp)==freeze['design_sha256'][case],'Prior input/design authentication')
    model=json.loads(mp.read_text());design=json.loads(dp.read_text());validate_model(model)
    require(design['input_sha256']==sha(mp),'Prior target/design binding')
    points=[p for p in design['points'] if p['state']=='HF']
    require(len(points)==1,'Missing unique frozen HF point')
    point=points[0];require(point['occupation']==(1<<model['electrons'])-1,'HF label differs')
    return model,point,dict(input_sha256=sha(mp),design_sha256=sha(dp),prior_freeze_sha256=sha(prior/'ALL_SCHEDULES_FROZEN.json'))


def prepare(output,preflight=False):
    require(not output.exists(),'Fresh output required; preserve partial attempts')
    output.mkdir(parents=True);commit=subprocess.check_output(['git','-C',str(HERE),'rev-parse','HEAD'],text=True).strip()
    for name in SOURCE_FILES:
        rel=str((HERE/name).relative_to(HERE.parents[1]))
        require(subprocess.check_output(['git','-C',str(HERE),'show',commit+':'+rel])==(HERE/name).read_bytes(),'Commit scientific source before freezing')
    write(output/'PROTOCOL_FROZEN.json',dict(source_commit=commit,sources=sources(),cases=CASES if not preflight else [CASES[-1]],
        width=9,T=1.,budgets=BUDGETS,arms=ARMS,expected_circuits=36 if not preflight else 12,
        preflight=preflight,python=sys.version,environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')}))
    artifacts={}
    for case in (CASES[-1:] if preflight else CASES):
        model,point,provenance=old_case(case)
        if preflight:
            # Alter coefficients without changing the predeclared register/units.
            model['coefficients']=[[w,str(7*int(c))] for w,c in model['coefficients']]
            model['blocks']=[[[w,str(7*int(c))] for w,c in b] for b in model['blocks']]
            model['denominator']=str(5*int(model['denominator']))
            model['removed_scalar_exact']=str(Fraction(model['removed_scalar_exact'])*Fraction(7,5))
            provenance['preflight_coefficient_scale']='7/5'
        domain=freeze_domain(model,point['methods']['optimized_JD']['order'])
        domain.update(case=case,prior=provenance,
            prior_construction={arm:{k:v for k,v in point['methods'][arm].items() if k.endswith('seconds') or k=='evaluations'}
                for arm in ('optimized_JD','current','pauli_JD_large')},
            control_orders={arm:point['methods'][arm]['order'] for arm in ('current','pauli_JD_large')})
        for folder,record in (('inputs',model),('domains',domain)):
            path=output/folder/(case+'.json');write(path,record);artifacts[str(path.relative_to(output))]=sha(path)
    artifacts['PROTOCOL_FROZEN.json']=sha(output/'PROTOCOL_FROZEN.json')
    write(output/'DOMAINS_FROZEN.json',dict(sources=sources(),artifacts=artifacts,
        policy='All supplied coefficients, weights, domains, warm starts and outside orders fixed before enumeration and physical errors'))


def bind(output):
    freeze=json.loads((output/'DOMAINS_FROZEN.json').read_text())
    require(freeze['sources']==sources(),'Scientific source changed since domain freeze')
    for p,h in freeze['artifacts'].items():require(sha(output/p)==h,'Frozen domain/input changed')
    protocol=json.loads((output/'PROTOCOL_FROZEN.json').read_text())
    require(protocol['environment']=={p:importlib.metadata.version(p) for p in protocol['environment']},'Pinned environment changed')
    return protocol


def record(model,order,partition):
    blocks=validate_model(model)
    if partition=='pauli':blocks=[{w:c} for w,c in sorted(decode(model['coefficients']).items())]
    require(sorted(order)==list(range(len(blocks))),'Control/solution order is not a permutation')
    return dict(order=order,partition=partition,unit_lengths=[len(blocks[i]) for i in order],
        schedule=[[[[q,p] for q,p in w],str(c)] for i in order for w,c in sorted(blocks[i].items())])


def certify(output):
    bind(output);require(not (output/'certificates').exists(),'Retain previous enumeration attempt')
    artifacts={}
    for path in sorted((output/'domains').glob('*.json')):
        domain=json.loads(path.read_text());model=json.loads((output/'inputs'/path.name).read_text())
        result=solve(domain,cpu_limit=120);result['domain_sha256']=sha(path)
        cp=output/'certificates'/path.name;write(cp,result)
        require(result['complete'],'Incomplete window: retain incumbent/lower bound; no physical run')
        orders={'heuristic':domain['warm_order'],'certified_window':result['order'],**domain['control_orders']}
        design=dict(case=domain['case'],records=[dict(arm=arm,formula=formula,
            **record(model,orders[arm],'pauli' if arm=='pauli_JD_large' else 'conserving')) for arm,formula in ARMS],
            domain_sha256=sha(path),certificate_sha256=sha(cp))
        dp=output/'designs'/path.name;write(dp,design)
        artifacts[str(cp.relative_to(output))]=sha(cp);artifacts[str(dp.relative_to(output))]=sha(dp)
        print('CERTIFIED',domain['case'],result['evaluated'],result['initial_JD'],result['incumbent_JD'],flush=True)
    artifacts['DOMAINS_FROZEN.json']=sha(output/'DOMAINS_FROZEN.json')
    write(output/'SCHEDULES_FROZEN.json',dict(artifacts=artifacts,sources=sources(),
        policy='All certified schedules and controls frozen before any physical propagation'))


def bind_schedules(output):
    protocol=bind(output);freeze=json.loads((output/'SCHEDULES_FROZEN.json').read_text())
    require(freeze['sources']==sources(),'Scientific source changed')
    for p,h in freeze['artifacts'].items():require(sha(output/p)==h,'Certified schedule changed')
    return protocol


class Target:
    def __init__(self,model):
        blocks=validate_model(model);self.model=model;self.dim=1<<model['n'];self.v=(1<<model['electrons'])-1
        basis=np.arange(self.dim);self.actions={};den=int(model['denominator'])
        for word,c in decode(model['coefficients']).items():
            desc=e.descriptor({word:1});x=desc[0]
            phase=np.array([e.value(desc,int(v)) for v in basis],float)
            self.actions[word]=(c/den,basis^x,phase)

    def hamiltonian(self):
        from scipy import sparse
        basis=np.arange(self.dim);H=sparse.csr_matrix((self.dim,self.dim),dtype=complex)
        for c,indices,phase in self.actions.values():H+=sparse.csr_matrix((c*phase,(indices,basis)),shape=H.shape)
        return H

    def propagate(self,rec,steps):
        psi=np.zeros(self.dim,complex);psi[self.v]=1
        ids,dt=e.pulses(rec['unit_lengths'],steps,rec['formula']==2)
        words=[tuple((q,p) for q,p in w) for w,c in rec['schedule']]
        for index,fraction in zip(ids,dt):
            c,indices,phase=self.actions[words[index]];theta=c*fraction
            psi=np.cos(theta)*psi-1j*np.sin(theta)*phase[indices]*psi[indices]
        return psi,len(ids)


def metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));a=actual/norm;z=exact/np.linalg.norm(exact);overlap=np.vdot(z,a)
    residual=a-overlap*z;loss=float(np.vdot(residual,residual).real);hf=(1<<model['electrons'])-1;n=model['n']
    target=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)]
    number=np.array([v.bit_count()!=model['electrons'] for v in range(len(a))])
    joint=np.array([[sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=target for v in range(len(a))])
    return dict(error=loss/(1+abs(overlap)),infidelity=loss,norm=norm,
        number_leakage=float(np.vdot(a[number],a[number]).real),joint_sector_leakage=float(np.vdot(a[joint],a[joint]).real))


def run(output):
    protocol=bind_schedules(output);require(not (output/'rows').exists(),'Retain prior physical attempt')
    (output/'rows').mkdir();(output/'states').mkdir();count=0
    from scipy.sparse.linalg import expm_multiply
    for path in sorted((output/'designs').glob('*.json')):
        design=json.loads(path.read_text());model=json.loads((output/'inputs'/path.name).read_text());target=Target(model)
        initial=np.zeros(target.dim,complex);initial[target.v]=1;H=target.hamiltonian()
        cpu=time.process_time();wall=time.monotonic();exact=expm_multiply(-1j*H,initial)
        reference_cost=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
        for index,rec in enumerate(design['records']):
            for budget in BUDGETS:
                cap=budget*len(model['coefficients']);steps=e.feasible_steps(rec['unit_lengths'],cap,rec['formula']==2)
                require(steps>=1,'No affordable physical circuit')
                cpu=time.process_time();wall=time.monotonic();actual,rotations=target.propagate(rec,steps)
                cost=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
                label=f"{design['case']}.{index}.b{budget}";sp=output/'states'/(label+'.npz')
                np.savez_compressed(sp,actual=actual)
                row=dict(case=design['case'],record_index=index,arm=rec['arm'],formula=rec['formula'],budget_units=budget,
                    steps=steps,rotations=rotations,rotation_cap=cap,**metrics(exact,actual,model),
                    propagation_cost=cost,reference_cost=reference_cost,design_sha256=sha(path),state_file=str(sp.relative_to(output)),state_sha256=sha(sp))
                require(rotations<=cap and abs(row['norm']-1)<1e-8,'Physical cap/norm violation')
                write(output/'rows'/(label+'.json'),row);count+=1
        print('EVOLVED',design['case'],flush=True)
    require(count==protocol['expected_circuits'],'Incomplete physical grid');write(output/'RUN_COMPLETE.json',dict(rows=count))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','certify','run'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--preflight',action='store_true');a=p.parse_args()
    if a.phase=='prepare':prepare(a.output.resolve(),a.preflight)
    else:require(not a.preflight,'Preflight flag applies only to prepare');globals()[a.phase](a.output.resolve())


if __name__=='__main__':main()
