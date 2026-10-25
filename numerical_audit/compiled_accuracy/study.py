"""Freeze a declared exact compiler grid before any compiled-state errors."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from fractions import Fraction
from functools import lru_cache
import argparse,hashlib,importlib.metadata,importlib.util,json,math,re,shutil,subprocess,sys,time
from pathlib import Path,PurePosixPath
import numpy as np
from compiler import Compiler,CompilationLimit,VERSION,resources,write_trace,read_trace,require

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
OLD=HERE.parent/'prospective_reference'
spec=importlib.util.spec_from_file_location('compiled_accuracy_prior_engine',OLD/'engine.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
CASES=('Li2_r2p95','F2_r1p55','H2O_scale1p10')
TIMES=(.83,2.11)
GRID=(1,2,4,8,16,32,64,128)
TARGETS=(.01,.001)
AMBIGUITY=2e-12
ARMS=(('current',1),('first_JD_sameT',1),('pauli_JD_large_sameT',1),
      ('current',2),('first_JD_sameT',2),('pauli_JD_large_sameT',2),
      ('cubic_orientation',2),('cubic_insertion',2))
PRODUCER=('compiler.py','study.py','PROTOCOL.md')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temporary=path.with_suffix('.partial')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)
def relative(name):
    p=PurePosixPath(name);require(not p.is_absolute() and '..' not in p.parts and str(p)==name,'Unsafe relative artifact path')
    return name
def source_files():
    return sorted([p for p in HERE.iterdir() if p.suffix in ('.py','.md','.txt')]+[OLD/'engine.py'])
def sources():return {str(p.relative_to(REPO)):sha(p) for p in source_files()}
def producer_sources():
    return {str(p.relative_to(REPO)):sha(p) for p in [HERE/n for n in PRODUCER]+[OLD/'engine.py']}
def verifier_sources():return {n:sha(HERE/n) for n in ('reference.py','verify.py')}

def authenticate_upstream(evidence,commit,preflight):
    """Authenticate source and payload; preflight data are explicitly untracked."""
    require(re.fullmatch('[0-9a-f]{40}',commit) is not None,'Full ordering source commit required')
    manifest=read(evidence/'MANIFEST.json');expected=set(manifest['files'])|{'MANIFEST.json'}
    actual={str(p.relative_to(evidence)) for p in evidence.rglob('*') if p.is_file()}
    require(actual==expected and len(actual)<=2000,'Ordering delivery inventory differs')
    require(sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file())<64*1024**2,'Ordering payload storage bound')
    tree={}
    if not preflight:
        prefix='numerical_audit/symmetric_response/evidence/'
        listing=subprocess.check_output(['git','-C',str(REPO),'ls-tree','-rz',commit,'--',prefix])
        for item in listing.split(b'\0'):
            if not item:continue
            metadata,name=item.split(b'\t',1);tree[name.decode()[len(prefix):]]=metadata.split()[2].decode()
        require(set(tree)==actual,'Merged ordering commit has a different payload inventory')
        parents=subprocess.check_output(['git','-C',str(REPO),'show','-s','--format=%P',commit],text=True).split()
        require(len(parents)>=2,'Main requires the retained actual ordering merge commit')
    for name in actual:
        relative(name);path=evidence/name;require(not path.is_symlink(),'Ordering symlink refused')
        data=path.read_bytes()
        if name!='MANIFEST.json':require(hashlib.sha256(data).hexdigest()==manifest['files'][name],'Ordering artifact changed: '+name)
        if tree:require(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==tree[name],'Ordering Git blob differs: '+name)
    original_sources={}
    for name,digest in manifest['sources'].items():
        relative(name);data=subprocess.check_output(['git','-C',str(REPO),'show',commit+':numerical_audit/'+name])
        require(hashlib.sha256(data).hexdigest()==digest,'Ordering source does not match its pinned commit')
        original_sources[name]=data
    protocol=read(evidence/'results/PROTOCOL_FROZEN.json');verified=read(evidence/'results/VERIFIED.json')
    require(protocol['preflight']==preflight and protocol['times']==([2.11] if preflight else list(TIMES)),'Ordering study scope differs')
    require(verified['status']=='passed' and verified['decision_agreement'],'Ordering prerequisite is unverified')
    require(verified['rows']==(16 if preflight else 96) and verified['comparisons']==(56 if preflight else 336),'Ordering panel incomplete')
    require(verified['schedule_freeze_sha256']==sha(evidence/'results/SCHEDULES_FROZEN.json'),'Ordering freeze differs')
    return manifest,original_sources,dict(ordering_repository='bernalde/qhat',ordering_commit=commit,
        dependency_kind='development_preflight' if preflight else 'merged',payload_Git_verified=not preflight,
        manifest_sha256=sha(evidence/'MANIFEST.json'),producing_source_commit=protocol['source_commit'],
        limitation='Preflight outputs were retained outside Git; their hashes, source revision and altered model are authenticated separately' if preflight else
        'Git content is authenticated locally; coordinator verifies the actual remote merge before publication')

def validate_model(model):
    n=model['n'];require(type(n) is int and 1<=n<=8,'Register bound before allocation')
    require(type(model['electrons']) is int and 0<=model['electrons']<=n,'Invalid occupation')
    require(0<len(model['blocks'])<=128 and 0<len(model['coefficients'])<=128,'Term/unit bound')
    require(type(model['denominator']) in (str,int) and int(model['denominator'])>0,'Positive integer denominator')
    terms={}
    for raw,c in model['coefficients']:
        word=tuple((q,p) for q,p in raw)
        require(word and len(dict(word))==len(word) and all(type(q) is int and 0<=q<n and p in 'XYZ' for q,p in word),'Scalar-removed valid support required')
        require(sum(p=='Y' for q,p in word)%2==0,'The inherited producing engine requires real even-Y target terms')
        require(word not in terms and type(c) in (str,int),'Unique integer Pauli numerators required');terms[word]=int(c)
    supplied=[(tuple((q,p) for q,p in w),int(c)) for block in model['blocks'] for w,c in block]
    require(all(model['blocks']) and len(supplied)==len(terms) and dict(supplied)==terms,'Complete units do not partition the target')
    for block in model['blocks']:
        words=[dict(w) for w,c in block]
        for i,P in enumerate(words):
            for Q in words[i+1:]:
                require(sum(P[q]!=Q[q] for q in P.keys()&Q.keys())%2==0,'Boundary fusion requires internally commuting unit Pauli words')
    return terms

def original_model(case,preflight):
    path=OLD/'evidence/results/inputs'/case/'model.json';manifest=read(OLD/'evidence/results/INPUTS_FROZEN.json')
    require(sha(path)==manifest[case],'Original target authentication differs');model=read(path)
    if preflight:
        model['case']='preflight_altered_'+case
        model['coefficients']=[[w,str(7*int(c))] for w,c in model['coefficients']]
        model['blocks']=[[[w,str(7*int(c))] for w,c in block] for block in model['blocks']]
        model['denominator']=str(5*int(model['denominator']))
        if 'removed_scalar_exact' in model:model['removed_scalar_exact']=str(Fraction(model['removed_scalar_exact'])*Fraction(7,5))
    return model

def validate_record(model,record):
    terms=validate_model(model);partition='pauli' if record['arm']=='pauli_JD_large_sameT' else 'conserving'
    blocks=[{tuple((q,p) for q,p in w):int(c) for w,c in block} for block in model['blocks']]
    if partition=='pauli':blocks=[{w:c} for w,c in sorted(terms.items())]
    order=record['order'];require(all(type(i) is int for i in order) and sorted(order)==list(range(len(blocks))),'Integral complete order required')
    lengths=[len(blocks[i]) for i in order]
    schedule=[[[[q,p] for q,p in word],str(c)] for i in order for word,c in sorted(blocks[i].items())]
    require(record['partition']==partition and record['unit_lengths']==lengths and record['schedule']==schedule,'Frozen unit circuit identity differs')
    require((record['arm'],record['formula']) in ARMS and record['budget_units']==30,'Only all eight fixed30U training arms are allowed')

def fused_units(lengths,steps,formula):
    require(type(steps) is int and steps in GRID and formula in (1,2),'Declared step grid required')
    m=len(lengths);last=None;duration=Fraction(0)
    for _ in range(steps):
        sequence=([(i,Fraction(1,steps)) for i in range(m)] if formula==1 else
            [(i,Fraction(1,2*steps)) for i in range(m-1)]+[(m-1,Fraction(1,steps))]+[(i,Fraction(1,2*steps)) for i in reversed(range(m-1))])
        for unit,tau in sequence:
            if unit==last:duration+=tau
            else:
                if last is not None:yield last,duration
                last,duration=unit,tau
    if last is not None:yield last,duration

def compile_candidate(model,record,steps):
    validate_record(model,record);compiler=Compiler(model['n']);cpu=time.process_time();wall=time.monotonic();rotations=0
    try:
        compiler.prepare_hf(model['electrons']);T=Fraction(*float(record['T']).as_integer_ratio());den=int(model['denominator'])
        offsets=[0]
        for size in record['unit_lengths']:offsets.append(offsets[-1]+size)
        for unit,tau in fused_units(record['unit_lengths'],steps,record['formula']):
            for at in range(offsets[unit],offsets[unit+1]):
                word,c=record['schedule'][at];compiler.rotation([(q,p) for q,p in word],Fraction(int(c),den)*T*tau);rotations+=1
    except CompilationLimit as error:
        return compiler,dict(status='censored',reason=str(error),generated=compiler.generated,stored=len(compiler.gates),
            completed_ideal_rotations=rotations,lowering_CPU_seconds=time.process_time()-cpu,lowering_wall_seconds=time.monotonic()-wall)
    return compiler,dict(status='complete',reason='complete fixed-policy compilation',generated=compiler.generated,stored=len(compiler.gates),
        ideal_Pauli_rotations=rotations,lowering_CPU_seconds=time.process_time()-cpu,lowering_wall_seconds=time.monotonic()-wall)

def _construct(model,record,steps,label,curve,output):
    compiler,result=compile_candidate(model,record,steps)
    identity=dict(candidate=label,curve=curve,case=model['case'],T=record['T'],T_ratio=list(float(record['T']).as_integer_ratio()),steps=steps,
        formula=record['formula'],arm=record['arm'],compilation_status=result['status'],
        model_sha256=sha(output/'models'/(model['case']+'.json')),
        curve_sha256=sha(output/'curves'/(curve+'.json')),protocol_sha256=sha(output/'PROTOCOL_FROZEN.json'))
    cpu=time.process_time();wall=time.monotonic();cost=resources(model['n'],compiler.gates)
    result['resource_counting']=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
    complete=result['status']=='complete';suffix='.ndjson.gz' if complete else '.partial.ndjson.gz'
    trace=output/'traces'/(label+suffix);trace.parent.mkdir(exist_ok=True)
    serialization={};digest=write_trace(trace,model['n'],compiler.gates,identity,measurements=serialization)
    if complete:result.update(compiled_resources=cost,trace_file=str(trace.relative_to(output)),trace_sha256=digest,trace_bytes=trace.stat().st_size)
    else:result.update(partial_compiled_resources=cost,partial_trace_file=str(trace.relative_to(output)),partial_trace_sha256=digest,partial_trace_bytes=trace.stat().st_size)
    result.update(identity=identity,candidate=label,curve=curve,steps=steps,case=model['case'],T=record['T'],arm=record['arm'],formula=record['formula'],
        inherited_ordering_CPU_charge=record['standalone_construction_CPU_charge'],**serialization)
    del compiler
    return result


def construct_candidate(model,record,steps,label,curve,output):
    cpu=time.process_time();wall=time.monotonic();result=_construct(model,record,steps,label,curve,output)
    elapsed=time.process_time()-cpu
    accounted=result['lowering_CPU_seconds']+sum(result[k]['CPU_seconds'] for k in ('resource_counting','serialization','hashing'))
    result.update(actual_local_constructor_CPU_seconds=elapsed,constructor_wall_seconds=time.monotonic()-wall,
        bookkeeping_CPU_seconds=elapsed-accounted,cost_timing='through complete constructor return')
    return result


def prepare(output,evidence,ordering_commit,preflight=False):
    began_prepare=time.process_time();wall_prepare=time.monotonic();constructor_CPU=0.
    require(not output.exists(),'Fresh output required; retain interrupted attempts');output.mkdir(parents=True)
    source_commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()
    source_hashes=sources()
    for path in source_files():
        relative_path=str(path.relative_to(REPO));data=subprocess.check_output(['git','-C',str(REPO),'show',source_commit+':'+relative_path])
        require(data==path.read_bytes(),'Commit every study source before execution')
        target=output/'source_snapshot'/relative_path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    write(output/'source_snapshot/MANIFEST.json',dict(source_commit=source_commit,sources=source_hashes))
    _,upstream_sources,provenance=authenticate_upstream(evidence,ordering_commit,preflight)
    shutil.copytree(evidence,output/'upstream/evidence')
    for name,data in upstream_sources.items():
        target=output/'upstream/source_snapshot'/relative(name);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    write(output/'upstream/PROVENANCE.json',provenance)
    cases=CASES[-1:] if preflight else CASES;times=TIMES[-1:] if preflight else TIMES
    write(output/'PROTOCOL_FROZEN.json',dict(source_commit=source_commit,sources=source_hashes,producer_sources=producer_sources(),
        prepared_verification_sources=verifier_sources(),upstream=provenance,preflight=preflight,cases=cases,times=times,
        arms=ARMS,training_rotation_cap_units=30,step_grid=GRID,accuracy_targets=TARGETS,threshold_ambiguity=AMBIGUITY,
        expected_curves=8 if preflight else 48,expected_candidates=64 if preflight else 384,compiler=VERSION,
        compiler_sha256=sha(HERE/'compiler.py'),python=sys.version,environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')}))
    count=feasible=censored=0
    for original_case in cases:
        case=('preflight_altered_' if preflight else '')+original_case
        model=read(evidence/'results/inputs'/(case+'.json'));require(model==original_model(original_case,preflight),'Supplied target changed')
        validate_model(model);write(output/'models'/(case+'.json'),model);design=read(evidence/'results/designs'/(case+'.json'))
        for ti,T in enumerate(times):
            records=[r for r in design['records'] if r['T']==T and r['budget_units']==30]
            require([(r['arm'],r['formula']) for r in records]==list(ARMS),'Fixed8armtraining cell differs')
            for ai,record in enumerate(records):
                validate_record(model,record);curve=case+'.T'+format(T,'.2f')+'.'+record['arm']+'.f'+str(record['formula'])
                write(output/'curves'/(curve+'.json'),dict(curve=curve,case=case,T=T,arm=record['arm'],formula=record['formula'],record=record,
                    upstream_design_sha256=sha(evidence/'results/designs'/(case+'.json')),model_sha256=sha(output/'models'/(case+'.json'))))
                for steps in GRID:
                    label=curve+'.r'+format(steps,'03d')
                    result=construct_candidate(model,record,steps,label,curve,output)
                    constructor_CPU+=result['actual_local_constructor_CPU_seconds']
                    if result['status']=='complete':feasible+=1
                    else:censored+=1
                    write(output/'candidates'/(label+'.json'),result);count+=1
                print('FROZE_COMPILED_CURVE',curve,flush=True)
    protocol=read(output/'PROTOCOL_FROZEN.json');require(count==protocol['expected_candidates'],'Incomplete predeclared compilation grid')
    total=time.process_time()-began_prepare
    write(output/'PREPARE_COSTS.json',dict(CPU_seconds=total,wall_seconds=time.monotonic()-wall_prepare,constructor_CPU_seconds=constructor_CPU,
        common_authentication_model_curve_and_metadata_CPU_seconds=total-constructor_CPU,
        scope='Through all source/upstream authentication,curve records and candidate JSON writes; final cost/freeze writes and interpreter startup remain in whole-process resources'))
    artifacts={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    write(output/'SCHEDULES_FROZEN.json',dict(artifacts=artifacts,candidates=count,feasible=feasible,censored=censored,
        policy='Everyfixedgridcandidateandexactfeasibletracefrozenbeforeanycompiledstateerror'))

def bind(output,*,independent=False):
    protocol=read(output/'PROTOCOL_FROZEN.json');frozen=read(output/'SCHEDULES_FROZEN.json')
    if protocol['producer_sources']!=producer_sources():
        require(independent,'Producing source changed since compilation freeze')
    snapshot=output/'source_snapshot';manifest=read(snapshot/'MANIFEST.json')
    require(re.fullmatch('[0-9a-f]{40}',protocol['source_commit']) is not None,'Full producing revision required')
    require(manifest['source_commit']==protocol['source_commit'] and manifest['sources']==protocol['sources'],'Original compiler source identity differs')
    expected_producer={'numerical_audit/compiled_accuracy/'+n for n in PRODUCER}|{'numerical_audit/prospective_reference/engine.py'}
    require(set(protocol['producer_sources'])==expected_producer and
        protocol['producer_sources']=={n:manifest['sources'].get(n) for n in expected_producer},'Producing source map is not a complete archived subset')
    require(manifest['sources'] and all(n=='numerical_audit/prospective_reference/engine.py' or
        (PurePosixPath(n).parent==PurePosixPath('numerical_audit/compiled_accuracy') and PurePosixPath(n).suffix in ('.py','.md','.txt'))
        for n in manifest['sources']),'Source archive domain differs')
    actual_snapshot={str(p.relative_to(snapshot)) for p in snapshot.rglob('*') if p.is_file()}
    require(actual_snapshot==set(manifest['sources'])|{'MANIFEST.json'},'Source archive inventory differs')
    for name,digest in manifest['sources'].items():
        path=snapshot/relative(name)
        require(not path.is_symlink() and sha(path)==digest,'Original source archive changed')
    for name,digest in frozen['artifacts'].items():require(sha(output/relative(name))==digest,'Frozen compiled artifact changed: '+name)
    require(set(protocol['environment'])=={'numpy','scipy'} and protocol['environment']=={p:importlib.metadata.version(p) for p in protocol['environment']},'Pinned numerical environment differs')
    require(protocol['step_grid']==list(GRID) and protocol['accuracy_targets']==list(TARGETS) and protocol['threshold_ambiguity']==AMBIGUITY,'Fixed accuracy/grid policy differs')
    return protocol,frozen

class PrimitiveReplay:
    def __init__(self,model):
        self.model=model;self.terms=validate_model(model);self.n=model['n'];self.dimension=1<<self.n
        basis=np.arange(self.dimension);self.low={};self.high={};self.cx={}
        for q in range(self.n):
            self.low[q]=basis[(basis&(1<<q))==0];self.high[q]=self.low[q]^(1<<q)
        for q in range(self.n):
            for t in range(self.n):
                if abs(q-t)==1:self.cx[q,t]=basis^(((basis>>q)&1)<<t)
    @lru_cache(maxsize=128)
    def phase(self,numerator,denominator):
        angle=float(Fraction(numerator,denominator))/2;require(math.isfinite(angle),'Nonfinite exact angle')
        return complex(math.cos(angle),-math.sin(angle))
    def replay(self,path,identity):
        iterator=read_trace(path);header=next(iterator)
        require(header['n']==self.n and header['identity']==identity,'Compiled trace identity differs')
        state=np.zeros(self.dimension,complex);state[0]=1
        for g in iterator:
            if g.kind=='cx':state=state[self.cx[g.a,g.b]]
            elif g.kind=='h':
                low=state[self.low[g.a]].copy();high=state[self.high[g.a]].copy()
                state[self.low[g.a]]=(low+high)/math.sqrt(2);state[self.high[g.a]]=(low-high)/math.sqrt(2)
            elif g.kind=='s':state[self.high[g.a]]*=1j
            elif g.kind=='sdg':state[self.high[g.a]]*=-1j
            else:
                phase=self.phase(g.angle.numerator,g.angle.denominator)
                state[self.low[g.a]]*=phase;state[self.high[g.a]]*=phase.conjugate()
        return state
    def hamiltonian(self):
        from scipy import sparse
        basis=np.arange(self.dimension);H=sparse.csr_matrix((self.dimension,self.dimension),dtype=complex)
        for word,c in self.terms.items():
            descriptor=engine.descriptor({word:1});ids=basis^descriptor[0]
            values=np.array([engine.value(descriptor,int(v)) for v in basis],float)
            H+=sparse.csr_matrix((c/int(self.model['denominator'])*values,(ids,basis)),shape=H.shape)
        return H

def metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));require(math.isfinite(norm) and norm>0,'Invalid compiled state norm')
    a=actual/norm;z=exact/np.linalg.norm(exact);inner=np.vdot(z,a);residual=a-inner*z;loss=float(np.vdot(residual,residual).real)
    n=model['n'];v=(1<<model['electrons'])-1;sector=[sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]
    number=np.array([w.bit_count()!=model['electrons'] for w in range(len(a))]);joint=np.array([[sum((w>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=sector for w in range(len(a))])
    return dict(error=loss/(1+abs(inner)),infidelity=loss,norm=norm,number_leakage=float(np.vdot(a[number],a[number]).real),joint_sector_leakage=float(np.vdot(a[joint],a[joint]).real))

def threshold(error,target):
    return 'ambiguous' if abs(error-target)<=AMBIGUITY else 'met' if error<target else 'not_met'

def selections(output):
    selected=[]
    for path in sorted((output/'curves').glob('*.json')):
        curve=read(path);candidates=[read(output/'candidates'/(curve['curve']+'.r'+format(r,'03d')+'.json')) for r in GRID]
        for target in TARGETS:
            qualifying=[];ambiguous=[]
            for candidate in candidates:
                if candidate['status']!='complete':continue
                row=read(output/'rows'/(candidate['candidate']+'.json'));classification=threshold(row['error'],target)
                if classification=='met':qualifying.append(row)
                elif classification=='ambiguous':ambiguous.append(row['candidate'])
            winner=min(qualifying,key=lambda x:(x['compiled_resources']['two_qubit_gates'],x['compiled_resources']['unit_duration_ASAP_depth'],x['compiled_resources']['gates'],x['steps'])) if qualifying else None
            selected.append(dict(curve=curve['curve'],case=curve['case'],T=curve['T'],arm=curve['arm'],formula=curve['formula'],target=target,
                status='met' if winner else 'unmet',selected_candidate=winner['candidate'] if winner else None,
                ambiguous_candidates=ambiguous,censored_candidates=sum(c['status']=='censored' for c in candidates),
                feasible_candidates=sum(c['status']=='complete' for c in candidates),qualifying_candidates=len(qualifying)))
    return selected

def run(output):
    began_run=time.process_time();wall_run=time.monotonic()
    protocol,frozen=bind(output);require(not (output/'rows').exists(),'Fresh physical phase required; retain partial attempts')
    (output/'rows').mkdir();(output/'states').mkdir();count=0
    from scipy.sparse.linalg import expm_multiply
    for model_path in sorted((output/'models').glob('*.json')):
        model=read(model_path);context=PrimitiveReplay(model);H=context.hamiltonian();initial=np.zeros(context.dimension,complex);initial[(1<<model['electrons'])-1]=1
        references={}
        for T in protocol['times']:
            cpu=time.process_time();wall=time.monotonic();exact=expm_multiply(-1j*T*H,initial)
            references[T]=(exact,dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall))
        for path in sorted((output/'candidates').glob(model['case']+'.*.json')):
            candidate=read(path)
            if candidate['status']=='censored':continue
            cpu=time.process_time();wall=time.monotonic();state=context.replay(output/candidate['trace_file'],candidate['identity'])
            replay=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall);exact,reference_cost=references[candidate['T']]
            label=candidate['candidate'];state_path=output/'states'/(label+'.npz');np.savez_compressed(state_path,actual=state)
            row={k:candidate[k] for k in ('candidate','curve','case','T','arm','formula','steps','ideal_Pauli_rotations','compiled_resources','inherited_ordering_CPU_charge')}
            row.update(**metrics(exact,state,model),compiled_replay_cost=replay,reference_cost=reference_cost,
                candidate_sha256=sha(path),trace_sha256=candidate['trace_sha256'],state_file=str(state_path.relative_to(output)),state_sha256=sha(state_path))
            require(abs(row['norm']-1)<1e-8,'Compiled circuit norm drift exceeds guard');write(output/'rows'/(label+'.json'),row);count+=1
        print('EVOLVED_COMPILED',model['case'],flush=True)
    require(count==frozen['feasible'],'Incomplete feasible compiled-state grid')
    write(output/'SELECTIONS.json',selections(output))
    write(output/'RUN_COMPLETE.json',dict(candidates=frozen['candidates'],feasible=count,censored=frozen['censored'],target_rows=2*protocol['expected_curves']))
    write(output/'RUN_COSTS.json',dict(CPU_seconds=time.process_time()-began_run,wall_seconds=time.monotonic()-wall_run,scope='Allbinding,reference,compiledreplay,metrics,state/rowwritesandselection;interpreter/finalrecordwritesinprocessresources'))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('phase',choices=('prepare','run'))
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--ordering-evidence',type=Path)
    parser.add_argument('--ordering-commit');parser.add_argument('--preflight',action='store_true');args=parser.parse_args()
    if args.phase=='prepare':
        require(args.ordering_evidence is not None and args.ordering_commit is not None,'Pinned ordering evidence/commit required')
        prepare(args.output.resolve(),args.ordering_evidence.resolve(),args.ordering_commit,args.preflight)
    else:require(not args.preflight,'Preflight flag only applies to preparation');run(args.output.resolve())
