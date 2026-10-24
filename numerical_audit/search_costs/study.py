"""Freeze controlled-cost searches, then separately evaluate ideal circuits."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import hashlib
import importlib.util
import importlib.metadata
import json
import resource
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from bounded import run_arm,engine as e,OLD
from models import synthetic

HERE=Path(__file__).resolve().parent
MOLECULAR=('Li2_r2p95','F2_r1p55','H2O_scale1p10')
SIZES=(4,6,8,10,12)
CPU_BUDGETS=(.5,2.)
REPETITIONS=3
ARMS=tuple((p,j,i) for p in ('conserving','pauli') for j in ('JD','J0') for i in ('cold','placement'))


def require(condition,message):
    if not condition:raise ValueError(message)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    part=path.with_suffix(path.suffix+'.partial')
    part.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');part.replace(path)


def decode(raw):return {tuple((int(q),p) for q,p in w):int(c) for w,c in raw}


def encoded(schedule):return [[[[q,p] for q,p in w],str(c)] for w,c in schedule]


def sources():
    files=[HERE/p for p in ('PROTOCOL.md','bounded.py','models.py','study.py','verify_costs.py')]
    files+=[OLD/p for p in ('engine.py','campaign.py','one_pass.py')]
    return {str(p.relative_to(HERE.parent)):sha(p) for p in files}


def load_previous_placement():
    # The prior campaign imports engine/one_pass lazily inside its unchanged
    # current_placement function. Resolve and authenticate those exact modules.
    sys.path.append(str(OLD))
    for name in ('engine','one_pass'):
        module=__import__(name)
        require(Path(module.__file__).resolve()==(OLD/(name+'.py')).resolve(),'Wrong prior helper import')
    spec=importlib.util.spec_from_file_location('search_cost_prior_campaign',OLD/'campaign.py')
    campaign=importlib.util.module_from_spec(spec);spec.loader.exec_module(campaign)
    return campaign.current_placement


def _common_body(model):
    cpu=time.process_time();wall=time.monotonic()
    if isinstance(model,Path):model=json.loads(model.read_text())
    require(model['n']<=12 and len(model['coefficients'])<=256,'Declared model bounds')
    blocks=[decode(b) for b in model['blocks']];total=decode(model['coefficients'])
    require(len(total)==sum(map(len,blocks)) and {w:c for b in blocks for w,c in b.items()}==total,'Units do not partition target')
    paulis=[{w:c} for w,c in sorted(total.items())]
    diag=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]==0),
                key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
    off=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]),
               key=lambda i:(-sum(abs(c) for c in blocks[i].values()),e.descriptor(blocks[i])[0]))
    D={w:c for i in diag for w,c in blocks[i].items()}
    return dict(blocks=blocks,paulis=paulis,diag=diag,off=off,D=D,den=int(model['denominator']),
        v=(1<<model['electrons'])-1,cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)


def common(model):
    # The constructor's complete return, including temporary teardown, is timed.
    cpu=time.process_time();wall=time.monotonic();context=_common_body(model)
    elapsed=time.process_time()-cpu
    context.update(internal_snapshot_CPU_seconds=context['cpu_seconds'],
        cpu_seconds=elapsed,wall_seconds=time.monotonic()-wall,cost_timing='through constructor return')
    return context


def _initializer_body(context,placement):
    cpu=time.process_time();wall=time.monotonic()
    result=placement(context['blocks'],context['diag'],context['off'],context['v'],context['den'])
    placement_CPU=time.process_time()-cpu;conversion=time.process_time()
    index={next(iter(b)):i for i,b in enumerate(context['paulis'])}
    flattened=[index[w] for i in result['order'] for w in sorted(context['blocks'][i])]
    ended=time.process_time();total_CPU=ended-cpu;conversion_CPU=ended-conversion
    return dict(conserving=result['order'],pauli=flattened,result=result,
                cpu_seconds=total_CPU,wall_seconds=time.monotonic()-wall,
                placement_CPU_seconds=placement_CPU,pauli_conversion_CPU_seconds=conversion_CPU,
                charges={'conserving':placement_CPU,'pauli':total_CPU})


def initializer(context,placement):
    cpu=time.process_time();wall=time.monotonic();warm=_initializer_body(context,placement)
    elapsed=time.process_time()-cpu
    warm.update(internal_snapshot_CPU_seconds=warm['cpu_seconds'],cpu_seconds=elapsed,
        wall_seconds=time.monotonic()-wall,cost_timing='through constructor return',
        pauli_conversion_CPU_seconds=elapsed-warm['placement_CPU_seconds'])
    warm['charges']['pauli']=elapsed
    return warm


def measured_arm(*args,**kwargs):
    """Time run_arm through return, including its feature/reference teardown.

    Results collection after this boundary is benchmark overhead and is retained
    separately in the whole-process ledger. Every warmup and measured arm uses
    this same boundary; an overrun is retained rather than relabelled feasible.
    """
    clock=kwargs.get('clock',time.process_time)
    began=clock();record=run_arm(*args,**kwargs);elapsed=clock()-began
    internal=record['actual_local_CPU_seconds']
    require(elapsed>=internal,'Through-return CPU cannot precede internal snapshot')
    record.update(cost_timing='through run_arm return',internal_snapshot_CPU_seconds=internal,
        return_time_CPU_seconds=elapsed-internal,actual_local_CPU_seconds=elapsed,
        accounting_overhead_CPU_seconds=elapsed-sum(record[k] for k in
            ('initializer_CPU_seconds','setup_CPU_seconds','search_CPU_seconds')))
    total=record['common_setup_CPU_charge']+record['inherited_placement_CPU_charge']+elapsed
    record.update(total_charged_CPU_seconds=total,CPU_overrun_seconds=max(0.,total-record['CPU_budget']),
        within_CPU_budget=total<=record['CPU_budget'])
    if not record['within_CPU_budget'] and not record['status'].startswith('CPU overrun retained: '):
        record['status']='CPU overrun retained: '+record['status']
    return record


def decorate(result,context,partition):
    blocks=context['paulis'] if partition=='pauli' else context['blocks']
    order=result['order'];require(sorted(order)==list(range(len(blocks))),'Search returned non-permutation')
    schedule=[(w,c) for i in order for w,c in sorted(blocks[i].items())]
    result.update(partition=partition,schedule=encoded(schedule),unit_lengths=[len(blocks[i]) for i in order])
    return result


def benchmark_model(model_path,placement):
    model=json.loads(model_path.read_text())
    warmup=[];context=common(model_path);warm=initializer(context,placement)
    # Same warmup CPU/evaluation policy for every factorial arm. No JIT is used.
    for partition,objective,initial in ARMS:
        bs=context['paulis'] if partition=='pauli' else context['blocks']
        result=measured_arm(bs,context['D'],context['v'],context['den'],objective,initial,
            warm[partition],CPU_BUDGETS[0],context['cpu_seconds'],warm['charges'][partition],20*len(bs)*(len(bs)-1))
        warmup.append(dict(partition=partition,objective=objective,initialization=initial,record=result))
    records=[];rep_setup=[]
    for rep in range(REPETITIONS):
        context=common(model_path);warm=initializer(context,placement)
        rep_setup.append(dict(repetition=rep,common_setup_cpu_seconds=context['cpu_seconds'],
            common_setup_wall_seconds=context['wall_seconds'],placement=warm))
        rotated=ARMS[rep:]+ARMS[:rep]
        for budget in CPU_BUDGETS:
            for ordinal,(partition,objective,initial) in enumerate(rotated):
                bs=context['paulis'] if partition=='pauli' else context['blocks']
                result=measured_arm(bs,context['D'],context['v'],context['den'],objective,initial,
                    warm[partition],budget,context['cpu_seconds'],warm['charges'][partition],20*len(bs)*(len(bs)-1))
                result.update(arm=partition+'_'+objective+'_'+initial,repetition=rep,CPU_budget=budget,
                              run_ordinal=ordinal,objective=objective,initialization=initial)
                records.append(decorate(result,context,partition))
        bs=context['paulis']
        larger=measured_arm(bs,context['D'],context['v'],context['den'],'JD','cold',warm['pauli'],
            8.,context['cpu_seconds'],warm['charges']['pauli'],20*len(bs)*(len(bs)-1))
        larger.update(arm='pauli_JD_larger_CPU',repetition=rep,CPU_budget=8.,run_ordinal=8,
                      label='unequal classical allowance; outside matched-CPU factorial')
        records.append(decorate(larger,context,'pauli'))
        baseline=dict(arm='current_placement',repetition=rep,CPU_budget=None,order=warm['conserving'],
            score=warm['result']['score'],evaluations=warm['result']['evaluations'],objective='JD',
            initialization='current placement',status='baseline separate timing',
            total_charged_CPU_seconds=context['cpu_seconds']+warm['charges']['conserving'],
            setup_CPU_seconds=warm['result']['setup_cpu_seconds'],
            search_CPU_seconds=warm['result']['search_cpu_seconds'])
        records.append(decorate(baseline,context,'conserving'))
    require(len(records)==54,'Incomplete declared search panel for input')
    return dict(case=model['case'],warmups=warmup,repetitions=rep_setup,records=records,
        measured_record_count=54,formulas=(1,2),T=1.,rotation_cap_units=10,
        process_peak_RSS_KiB_at_input_end=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def prepare(output):
    require(not output.exists(),'Fresh output required; retain partial attempts')
    output.mkdir(parents=True)
    commit=subprocess.check_output(['git','-C',str(HERE),'rev-parse','HEAD'],text=True).strip()
    for relative in sources():
        path=HERE.parent/relative
        blob=subprocess.check_output(['git','-C',str(HERE),'show',commit+':'+str(path.relative_to(HERE.parents[1]))])
        require(blob==path.read_bytes(),'Commit active scientific source before freeze')
    write(output/'PROTOCOL_FROZEN.json',dict(source_commit=commit,sources=sources(),CPU_budgets=CPU_BUDGETS,
        repetitions=REPETITIONS,arms=ARMS,synthetic_sizes=SIZES,molecular_cases=MOLECULAR,
        expected_orders=432,expected_circuits=864,python=sys.version,
        environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')}))
    models=[];prior_hashes={}
    old_inputs=json.loads((OLD/'evidence/results/INPUTS_FROZEN.json').read_text())
    for case in MOLECULAR:
        path=OLD/'evidence/results/inputs'/case/'model.json'
        require(sha(path)==old_inputs[case],'Molecular target authentication failed')
        model=json.loads(path.read_text());models.append(model);prior_hashes[case]=sha(path)
    models.extend(synthetic(n) for n in SIZES)
    placement=load_previous_placement();artifacts={'PROTOCOL_FROZEN.json':sha(output/'PROTOCOL_FROZEN.json')}
    for model in models:
        mp=output/'inputs'/(model['case']+'.json');write(mp,model)
        design=benchmark_model(mp,placement);design['input_sha256']=sha(mp)
        path=output/'designs'/(model['case']+'.json');write(path,design)
        for p in (mp,path):artifacts[str(p.relative_to(output))]=sha(p)
        print('FROZE_SEARCH',model['case'],'words',len(model['coefficients']),'orders',len(design['records']),flush=True)
    write(output/'SCHEDULES_FROZEN.json',dict(sources=sources(),artifacts=artifacts,prior_models=prior_hashes,
        orders=432,expected_circuits=864))


def bind(output):
    freeze=json.loads((output/'SCHEDULES_FROZEN.json').read_text())
    require(freeze['sources']==sources(),'Active scientific source changed')
    protocol=json.loads((output/'PROTOCOL_FROZEN.json').read_text())
    require(protocol['environment']=={p:importlib.metadata.version(p) for p in protocol['environment']},'Pinned numerical environment changed')
    for p,h in freeze['artifacts'].items():require(sha(output/p)==h,'Frozen model/design changed')
    for case,h in freeze['prior_models'].items():require(sha(OLD/'evidence/results/inputs'/case/'model.json')==h,'Prior model changed')
    return freeze


class CircuitTarget:
    def __init__(self,model):
        n=model['n'];require(0<n<=12 and len(model['coefficients'])<=256,'Circuit allocation guard')
        self.dim=1<<n;basis=np.arange(self.dim);self.actions={};self.H=None
        self.v=(1<<model['electrons'])-1
        for word,c in decode(model['coefficients']).items():
            desc=e.descriptor({word:c});x=desc[0]
            require(x<self.dim,'Word outside declared register')
            phase=np.array([e.value(e.descriptor({word:1}),int(v)) for v in basis],float)
            self.actions[word]=(c/int(model['denominator']),basis^x,phase)
        self.payload_bytes=sum(ids.nbytes+phase.nbytes for c,ids,phase in self.actions.values())
        require(self.payload_bytes<=64*1024**2,'Circuit action cache exceeds 64 MiB')

    def hamiltonian(self):
        from scipy import sparse
        basis=np.arange(self.dim);H=sparse.csr_matrix((self.dim,self.dim),dtype=complex)
        for c,permutation,phase in self.actions.values():
            H+=sparse.csr_matrix((c*phase,(permutation,basis)),shape=H.shape)
        return H

    def propagate(self,record,steps,formula):
        psi=np.zeros(self.dim,complex);psi[self.v]=1
        ids,fractions=e.pulses(record['unit_lengths'],steps,formula==2)
        words=[tuple((q,p) for q,p in w) for w,c in record['schedule']]
        for index,fraction in zip(ids,fractions):
            c,permutation,phase=self.actions[words[index]];theta=c*fraction
            psi=np.cos(theta)*psi-1j*np.sin(theta)*phase[permutation]*psi[permutation]
        return psi,len(ids)


def physical_metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));psi=actual/norm;phi=exact/np.linalg.norm(exact)
    overlap=np.vdot(phi,psi);residual=psi-phi*overlap;infidelity=float(np.vdot(residual,residual).real)
    n=model['n'];hf=(1<<model['electrons'])-1
    counts=lambda v,p:sum((v>>q)&1 for q in range(p,n,2))
    correct_number=np.array([v.bit_count()==model['electrons'] for v in range(len(psi))])
    correct_joint=np.array([all(counts(v,p)==counts(hf,p) for p in (0,1)) for v in range(len(psi))])
    return dict(error=infidelity/(1+abs(overlap)),infidelity=infidelity,norm=norm,
        number_leakage=float(np.vdot(psi[~correct_number],psi[~correct_number]).real),
        joint_sector_leakage=float(np.vdot(psi[~correct_joint],psi[~correct_joint]).real))


def run(output):
    bind(output);require(not (output/'rows').exists(),'Retain existing physical attempts')
    (output/'rows').mkdir();(output/'states').mkdir();count=0
    from scipy.sparse.linalg import expm_multiply
    for path in sorted((output/'designs').glob('*.json')):
        design=json.loads(path.read_text());model=json.loads((output/'inputs'/(design['case']+'.json')).read_text())
        target=CircuitTarget(model);initial=np.zeros(target.dim,complex);initial[target.v]=1
        H=target.hamiltonian();cpu=time.process_time();wall=time.monotonic();exact=expm_multiply(-1j*H,initial)
        reference_cost=dict(cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
        for index,record in enumerate(design['records']):
            for formula in (1,2):
                budget=10*len(model['coefficients']);steps=e.feasible_steps(record['unit_lengths'],budget,formula==2)
                cpu=time.process_time();wall=time.monotonic();actual,rotations=target.propagate(record,steps,formula)
                cost=dict(cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
                label=design['case']+'.'+str(index)+'.f'+str(formula)
                sp=output/'states'/(label+'.npz');np.savez_compressed(sp,actual=actual)
                row=dict(case=design['case'],record_index=index,arm=record['arm'],repetition=record['repetition'],
                    CPU_budget=record['CPU_budget'],formula=formula,steps=steps,rotations=rotations,rotation_cap=budget,
                    **physical_metrics(exact,actual,model),propagation_cost=cost,reference_cost=reference_cost,
                    state_file=str(sp.relative_to(output)),state_sha256=sha(sp),design_sha256=sha(path),
                    circuit_action_cache_bytes=target.payload_bytes)
                require(rotations<=budget and abs(row['norm']-1)<1e-8,'Physical budget/norm failure')
                write(output/'rows'/(label+'.json'),row);count+=1
        print('EVOLVED',design['case'],flush=True)
    require(count==864,'Incomplete fixed circuit grid');write(output/'RUN_COMPLETE.json',dict(rows=count))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','run'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();globals()[a.phase](a.output.resolve())


if __name__=='__main__':main()
