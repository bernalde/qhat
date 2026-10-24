"""Freeze symmetric-specific schedules, then separately propagate full states."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from fractions import Fraction
import argparse,hashlib,importlib.metadata,importlib.util,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from cubic import Context,FirstOrderContext,engine as e,OLD,decode,validate_model,require,insertion

HERE=Path(__file__).resolve().parent
CASES=('Li2_r2p95','F2_r1p55','H2O_scale1p10')
TIMES=(.83,2.11)
BUDGETS=(10,30)
ARMS=(('current',1),('first_JD_sameT',1),('pauli_JD_large_sameT',1),
      ('current',2),('first_JD_sameT',2),('pauli_JD_large_sameT',2),
      ('cubic_orientation',2),('cubic_insertion',2))
CORE=('cubic.py','independent.py','study.py','verify.py','PROTOCOL.md','DERIVATION.md')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temporary=path.with_suffix('.partial')
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def sources():
    return {str(p.relative_to(HERE.parent)):sha(p) for p in [HERE/n for n in CORE]+[OLD/n for n in ('engine.py','campaign.py','one_pass.py')]}


def old_model(case):
    prior=OLD/'evidence/results';f=json.loads((prior/'INPUTS_FROZEN.json').read_text());path=prior/'inputs'/case/'model.json'
    require(sha(path)==f[case],'Retained target authentication differs');return json.loads(path.read_text()),sha(path)


def placement_function():
    sys.path.insert(0,str(OLD))
    for name in ('engine','one_pass'):
        module=__import__(name);require(Path(module.__file__).resolve()==(OLD/(name+'.py')).resolve(),'Wrong inherited helper')
    spec=importlib.util.spec_from_file_location('symmetric_response_prior_campaign',OLD/'campaign.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module.current_placement


def _base(model):
    blocks=validate_model(model)
    diag=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]==0),
        key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
    off=sorted((i for i,b in enumerate(blocks) if e.descriptor(b)[0]),
        key=lambda i:(-sum(abs(c) for c in blocks[i].values()),e.descriptor(blocks[i])[0]))
    return dict(blocks=blocks,diag=diag,off=off,v=(1<<model['electrons'])-1,den=int(model['denominator']))


def _first(model,T,partition,current):
    phases={};cpu=time.process_time();context=FirstOrderContext(model,T,partition)
    phases['method_setup_CPU_seconds']=time.process_time()-cpu
    cpu=time.process_time()
    if partition=='conserving':warm=list(current)
    else:
        original=[decode(b) for b in model['blocks']];index={next(iter(b)):i for i,b in enumerate(context.blocks)}
        warm=[index[w] for i in current for w in sorted(original[i])]
    label=list(range(context.m));l1=sorted(label,key=lambda i:(-sum(abs(c) for c in context.blocks[i].values()),i))
    orders=[label,l1,warm];phases['initializer_CPU_seconds']=time.process_time()-cpu
    cpu=time.process_time();cap=10*len(model['coefficients'])
    seeds=[context.score(order,T,cap) for order in orders]
    phases['seed_scores_CPU_seconds']=time.process_time()-cpu
    allowance=2. if partition=='conserving' else 8.
    evaluations=1000 if partition=='conserving' else 20*context.m*(context.m-1)
    moves=5 if partition=='conserving' else 10
    cpu=time.process_time();result=insertion(context,seeds,T,cap,allowance,evaluations,moves)
    phases['search_CPU_seconds']=time.process_time()-cpu
    cpu=time.process_time();del context;phases['context_cleanup_CPU_seconds']=time.process_time()-cpu
    return dict(result=result,seeds=seeds,partition=partition,T=T,phases=phases)


def first_method(model,T,partition,current):
    cpu=time.process_time();wall=time.monotonic();record=_first(model,T,partition,current)
    elapsed=time.process_time()-cpu
    record.update(actual_local_CPU_seconds=elapsed,wall_seconds=time.monotonic()-wall,
        bookkeeping_CPU_seconds=elapsed-sum(record['phases'].values()),cost_timing='through constructor return')
    return record


def _orientation(context,first,current,T,cap):
    cpu=time.process_time();orders=[list(first),list(reversed(first)),list(current)]
    initializer_CPU=time.process_time()-cpu
    cpu=time.process_time();seeds=[context.score(order,T,cap) for order in orders]
    seed_CPU=time.process_time()-cpu;orientation=min(seeds,key=lambda r:(r['score'],tuple(r['order'])))
    return dict(seeds=seeds,orientation=orientation,phases=dict(initializer_CPU_seconds=initializer_CPU,seed_scores_CPU_seconds=seed_CPU))


def orientation_method(context,first,current,T,cap):
    cpu=time.process_time();record=_orientation(context,first,current,T,cap);elapsed=time.process_time()-cpu
    record.update(actual_local_CPU_seconds=elapsed,bookkeeping_CPU_seconds=elapsed-sum(record['phases'].values()),cost_timing='through constructor return')
    return record


def _cubic_method(context,first,current,T,cap):
    orientation=orientation_method(context,first,current,T,cap);seeds=orientation['seeds']
    cpu=time.process_time();result=insertion(context,seeds,T,cap,2.,1000,5)
    search_CPU=time.process_time()-cpu
    return dict(seeds=seeds,orientation=orientation['orientation'],result=result,
        orientation_constructor_CPU_seconds=orientation['actual_local_CPU_seconds'],
        phases={**orientation['phases'],'orientation_bookkeeping_CPU_seconds':orientation['bookkeeping_CPU_seconds'],'search_CPU_seconds':search_CPU})


def cubic_method(context,first,current,T,cap):
    cpu=time.process_time();wall=time.monotonic();record=_cubic_method(context,first,current,T,cap)
    elapsed=time.process_time()-cpu
    record.update(actual_local_CPU_seconds=elapsed,wall_seconds=time.monotonic()-wall,
        bookkeeping_CPU_seconds=elapsed-sum(record['phases'].values()),cost_timing='through constructor return')
    return record


def circuit(model,order,partition,arm,formula,T,budget):
    blocks=[decode(b) for b in model['blocks']]
    if partition=='pauli':blocks=[{w:c} for w,c in sorted(decode(model['coefficients']).items())]
    require(all(type(i) is int for i in order) and sorted(order)==list(range(len(blocks))),'Invalid circuit order')
    return dict(arm=arm,formula=formula,T=T,budget_units=budget,partition=partition,order=order,
        unit_lengths=[len(blocks[i]) for i in order],
        schedule=[[[[q,p] for q,p in w],str(c)] for i in order for w,c in sorted(blocks[i].items())])


def design_model(model,placement,times,model_load_CPU=0.):
    cpu=time.process_time();base=_base(model);base_setup=time.process_time()-cpu
    cpu=time.process_time();current=placement(base['blocks'],base['diag'],base['off'],base['v'],base['den'])
    current_CPU=time.process_time()-cpu
    cpu=time.process_time();context=Context(model);cache_setup=time.process_time()-cpu
    cache_bytes=context.cache_bytes;records=[];time_records=[]
    for T in times:
        first=first_method(model,T,'conserving',current['order'])
        pauli=first_method(model,T,'pauli',current['order'])
        cells=[]
        for budget in BUDGETS:
            cap=budget*len(model['coefficients'])
            new=cubic_method(context,first['result']['record']['order'],current['order'],T,cap)
            orders={'current':current['order'],'first_JD_sameT':first['result']['record']['order'],
                'pauli_JD_large_sameT':pauli['result']['record']['order'],
                'cubic_orientation':new['orientation']['order'],'cubic_insertion':new['result']['record']['order']}
            start=len(records)
            for arm,formula in ARMS:
                records.append(circuit(model,orders[arm],'pauli' if arm=='pauli_JD_large_sameT' else 'conserving',arm,formula,T,budget))
            cells.append(dict(budget_units=budget,first_record_index=start,cubic=new))
        time_records.append(dict(T=T,first=first,pauli=pauli,cells=cells))
    cpu=time.process_time();del context;cache_cleanup=time.process_time()-cpu
    cpu=time.process_time();del base;base_cleanup=time.process_time()-cpu
    base_charge=model_load_CPU+base_setup+base_cleanup;cache_charge=cache_setup+cache_cleanup
    for tr in time_records:
        first_charge=base_charge+current_CPU+tr['first']['actual_local_CPU_seconds']
        pauli_charge=base_charge+current_CPU+tr['pauli']['actual_local_CPU_seconds']
        for cell in tr['cells']:
            new=cell['cubic'];orientation=first_charge+cache_charge+new['orientation_constructor_CPU_seconds']
            full=first_charge+cache_charge+new['actual_local_CPU_seconds']
            for record in records[cell['first_record_index']:cell['first_record_index']+8]:
                charge={'current':base_charge+current_CPU,'first_JD_sameT':first_charge,
                    'pauli_JD_large_sameT':pauli_charge,'cubic_orientation':orientation,'cubic_insertion':full}[record['arm']]
                record['standalone_construction_CPU_charge']=charge
    return dict(case=model['case'],records=records,time_records=time_records,current=current,
        construction=dict(model_load_CPU_seconds=model_load_CPU,base_setup_CPU_seconds=base_setup,base_cleanup_CPU_seconds=base_cleanup,
            legacy_current_CPU_seconds=current_CPU,legacy_current_reference_time=1.,
            cubic_cache_setup_CPU_seconds=cache_setup,cubic_cache_cleanup_CPU_seconds=cache_cleanup,cubic_cache_bytes=cache_bytes),
        scope='Fresh same-time first-order controls,originallegacycurrentT1,newsymmetriccubicorientation/insertion;allconstructorcosts charged,notmatchedtotalCPU')


def prepare(output,preflight=False):
    require(not output.exists(),'Fresh output required; preserve partial attempts');output.mkdir(parents=True)
    commit=subprocess.check_output(['git','-C',str(HERE),'rev-parse','HEAD'],text=True).strip()
    if not preflight:
        for name in CORE:
            relative=str((HERE/name).relative_to(HERE.parents[1]))
            require(subprocess.check_output(['git','-C',str(HERE),'show',commit+':'+relative])==(HERE/name).read_bytes(),'Commit scientific source before freeze')
    cases=CASES[-1:] if preflight else CASES;times=TIMES[-1:] if preflight else TIMES
    write(output/'PROTOCOL_FROZEN.json',dict(source_commit=commit,sources=sources(),cases=cases,times=times,budgets=BUDGETS,
        arms=ARMS,preflight=preflight,expected_circuits=16 if preflight else 96,python=sys.version,
        environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')}))
    placement=placement_function();artifacts={'PROTOCOL_FROZEN.json':sha(output/'PROTOCOL_FROZEN.json')}
    prior={}
    for case in cases:
        cpu=time.process_time();model,digest=old_model(case);model_load_CPU=time.process_time()-cpu;prior[case]=digest
        if preflight:
            model['case']='preflight_altered_'+case
            model['coefficients']=[[w,str(7*int(c))] for w,c in model['coefficients']]
            model['blocks']=[[[w,str(7*int(c))] for w,c in b] for b in model['blocks']]
            model['denominator']=str(5*int(model['denominator']))
            if 'removed_scalar_exact' in model:
                model['removed_scalar_exact']=str(Fraction(model['removed_scalar_exact'])*Fraction(7,5))
        path=output/'inputs'/(model['case']+'.json');write(path,model)
        design=design_model(model,placement,times,model_load_CPU);design['input_sha256']=sha(path)
        dp=output/'designs'/path.name;write(dp,design)
        for p in (path,dp):artifacts[str(p.relative_to(output))]=sha(p)
        print('FROZE_SYMMETRIC',model['case'],len(design['records']),'schedules',flush=True)
    write(output/'SCHEDULES_FROZEN.json',dict(sources=sources(),artifacts=artifacts,prior_models=prior,
        policy='Allinputs,times,controls,searchpolicies,costs andcircuits fixedbefore anynewphysical errors'))


def bind(output, *, allow_archived_sources=False):
    freeze=json.loads((output/'SCHEDULES_FROZEN.json').read_text())
    protocol=json.loads((output/'PROTOCOL_FROZEN.json').read_text())
    current=sources()
    require(set(freeze['sources'])==set(current),'Frozen source domain changed')
    if freeze['sources']!=current:
        require(allow_archived_sources,'Scientific source changed since freeze')
        archive=output/'source_snapshot'
        retained=json.loads((archive/'MANIFEST.json').read_text())
        require(retained['source_commit']==protocol['source_commit'] and retained['sources']==freeze['sources'],
            'Archived producing source identity differs')
        inventory={str(p.relative_to(archive)) for p in archive.rglob('*') if p.is_file()}
        require(inventory==set(freeze['sources'])|{'MANIFEST.json'},'Archived source inventory differs')
        for name,digest in freeze['sources'].items():
            require(sha(archive/name)==digest,'Archived producing source changed: '+name)
        # Only independent verification opts into this path. Physical propagation
        # still requires its exact live producing sources. The verifier records
        # its own current versions separately from the archived producer.

    for p,h in freeze['artifacts'].items():require(sha(output/p)==h,'Frozen input/design changed')
    for case,h in freeze['prior_models'].items():require(old_model(case)[1]==h,'Prior target authentication changed')
    protocol=json.loads((output/'PROTOCOL_FROZEN.json').read_text())
    require(protocol['environment']=={p:importlib.metadata.version(p) for p in protocol['environment']},'Pinned environment changed')
    return protocol


class Target:
    def __init__(self,model):
        validate_model(model);self.model=model;self.dimension=1<<model['n'];self.v=(1<<model['electrons'])-1
        basis=np.arange(self.dimension);self.actions={}
        for word,c in decode(model['coefficients']).items():
            d=e.descriptor({word:1});ids=basis^d[0]
            values=np.array([e.value(d,int(v)) for v in basis],float)
            self.actions[word]=(c/int(model['denominator']),ids,values)

    def hamiltonian(self):
        from scipy import sparse
        basis=np.arange(self.dimension);H=sparse.csr_matrix((self.dimension,self.dimension),dtype=complex)
        for c,ids,values in self.actions.values():H+=sparse.csr_matrix((c*values,(ids,basis)),shape=H.shape)
        return H

    def propagate(self,record,steps,T):
        state=np.zeros(self.dimension,complex);state[self.v]=1
        ids,fractions=e.pulses(record['unit_lengths'],steps,record['formula']==2)
        words=[tuple((q,p) for q,p in w) for w,c in record['schedule']]
        for index,fraction in zip(ids,fractions):
            c,permutation,phase=self.actions[words[index]];theta=c*T*fraction
            state=np.cos(theta)*state-1j*np.sin(theta)*phase[permutation]*state[permutation]
        return state,len(ids)


def metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));a=actual/norm;z=exact/np.linalg.norm(exact);overlap=np.vdot(z,a)
    residual=a-overlap*z;loss=float(np.vdot(residual,residual).real);n=model['n'];hf=(1<<model['electrons'])-1
    sector=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)]
    number=np.array([v.bit_count()!=model['electrons'] for v in range(len(a))])
    joint=np.array([[sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=sector for v in range(len(a))])
    return dict(error=loss/(1+abs(overlap)),infidelity=loss,norm=norm,
        number_leakage=float(np.vdot(a[number],a[number]).real),joint_sector_leakage=float(np.vdot(a[joint],a[joint]).real))


def run(output):
    protocol=bind(output);require(not (output/'rows').exists(),'Fresh physical phase required; preserve partial attempts')
    (output/'rows').mkdir();(output/'states').mkdir();count=0
    from scipy.sparse.linalg import expm_multiply
    for path in sorted((output/'designs').glob('*.json')):
        design=json.loads(path.read_text());model=json.loads((output/'inputs'/path.name).read_text());target=Target(model)
        H=target.hamiltonian();initial=np.zeros(target.dimension,complex);initial[target.v]=1;references={}
        for record in design['records']:
            T=record['T']
            if T not in references:
                cpu=time.process_time();wall=time.monotonic();state=expm_multiply(-1j*T*H,initial)
                references[T]=(state,dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall))
        for index,record in enumerate(design['records']):
            cap=record['budget_units']*len(model['coefficients']);steps=e.feasible_steps(record['unit_lengths'],cap,record['formula']==2)
            cpu=time.process_time();wall=time.monotonic();actual,rotations=target.propagate(record,steps,record['T'])
            cost=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
            label=model['case']+'.'+str(index);sp=output/'states'/(label+'.npz');np.savez_compressed(sp,actual=actual)
            exact,refcost=references[record['T']]
            row={k:record[k] for k in ('arm','formula','T','budget_units','standalone_construction_CPU_charge')}
            row.update(case=model['case'],record_index=index,steps=steps,rotations=rotations,rotation_cap=cap,
                **metrics(exact,actual,model),propagation_cost=cost,reference_cost=refcost,
                design_sha256=sha(path),state_file=str(sp.relative_to(output)),state_sha256=sha(sp))
            require(rotations<=cap and abs(row['norm']-1)<1e-8,'Physical cap/norm violation')
            write(output/'rows'/(label+'.json'),row);count+=1
        print('EVOLVED_SYMMETRIC',model['case'],flush=True)
    require(count==protocol['expected_circuits'],'Incomplete fixed physical grid');write(output/'RUN_COMPLETE.json',dict(rows=count))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','run'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--preflight',action='store_true');a=p.parse_args()
    if a.phase=='prepare':prepare(a.output.resolve(),a.preflight)
    else:require(not a.preflight,'Preflight flag applies to prepare only');run(a.output.resolve())
