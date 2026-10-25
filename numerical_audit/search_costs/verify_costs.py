"""Independent literal-Pauli response, cost-ledger audit and unfused replay."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse


def decode(raw):return {tuple((q,p) for q,p in w):int(c) for w,c in raw}


def literal(word,n):
    if not 0<n<=12 or len(dict(word))!=len(word) or any(not 0<=q<n or p not in 'XYZ' for q,p in word):
        raise ValueError('Independent Pauli allocation/word guard')
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    labels=dict(word);out=sparse.csr_matrix([[1.+0j]])
    for q in reversed(range(n)):out=sparse.kron(out,matrices[labels.get(q,'I')],format='csr')
    return out


class Reference:
    def __init__(self,model):
        n=model['n'];P=len(model['coefficients'])
        if not 0<n<=12 or not 0<P<=256:raise ValueError('Independent reference size guard')
        dim=1<<n
        if P*dim*24>64*1024**2:raise ValueError('Literal sparse pool storage guard')
        self.model=model;self.dim=dim;self.v=(1<<model['electrons'])-1
        self.pool={w:literal(w,n) for w in decode(model['coefficients'])}
        self.H=sparse.csr_matrix((dim,dim),dtype=complex)
        for w,c in decode(model['coefficients']).items():self.H+=(c/int(model['denominator']))*self.pool[w]
        self.initial=np.zeros(dim,complex);self.initial[self.v]=1
        charges=[sparse.diags([sum((v>>q)&1 for q in range(p,n,2)) for v in range(dim)],format='csr',dtype=float) for p in (0,1)]
        for block in model['blocks']:
            terms=[(w,c) for w,c in decode(block).items()]
            G=sum((c/int(model['denominator']))*self.pool[w] for w,c in terms)
            for N in charges:
                comm=G@N-N@G
                if comm.nnz and max(abs(comm.data))>2e-12:raise ValueError('Completed unit does not conserve declared spin occupations')
            for at,(w,c) in enumerate(terms):
                for z,d in terms[at+1:]:
                    comm=self.pool[w]@self.pool[z]-self.pool[z]@self.pool[w]
                    if comm.nnz and max(abs(comm.data))>1e-13:raise ValueError('Unit flattening is not internally commuting')

    def exact_action(self):
        # Independent scaled Taylor reference. ||H||_2 <= ||H||_infinity for
        # Hermitian H. Each substep has norm <=.25; 24 Taylor terms make the
        # analytic truncation negligible here. Binary64 roundoff is checked
        # against the separately produced SciPy exponential, not certified.
        bound=float(max(np.asarray(abs(self.H).sum(axis=1)).ravel()))
        steps=max(1,int(np.ceil(bound/.25)));state=self.initial.copy()
        for _ in range(steps):
            term=state.copy();result=state.copy()
            for k in range(1,25):term=(-1j/steps/k)*(self.H@term);result+=term
            state=result
        return state

    def scores(self,record):
        model=self.model
        bs=([{w:c} for w,c in sorted(decode(model['coefficients']).items())]
             if record['partition']=='pauli' else [decode(b) for b in model['blocks']])
        prefix=sparse.csr_matrix(self.H.shape,dtype=complex);column=np.zeros(self.dim,complex)
        for index in record['order']:
            G=sparse.csr_matrix(self.H.shape,dtype=complex)
            for w,c in bs[index].items():G+=(c/int(model['denominator']))*self.pool[w]
            column+=(-.5j)*(prefix@(G@self.initial)-G@(prefix@self.initial));prefix+=G
        gaps=self.H.diagonal().real-float(self.H[self.v,self.v].real)
        column[self.v]=0
        unweighted=float(sum(abs(column)**2))
        column*=np.exp(.5j*gaps)*np.sinc(gaps/(2*np.pi))
        return dict(JD=float(sum(abs(column)**2)),J0=unweighted)

    def score(self,record):return self.scores(record)[record['objective']]

    def replay(self,record,steps,formula):
        schedule=[(tuple((q,p) for q,p in w),int(c)/int(self.model['denominator'])) for w,c in record['schedule']]
        offsets=np.cumsum([0]+record['unit_lengths']);actual=self.initial.copy()
        def apply(unit,dt):
            nonlocal actual
            for k in range(offsets[unit],offsets[unit+1]):
                word,c=schedule[k];theta=c*dt
                actual=np.cos(theta)*actual-1j*np.sin(theta)*(self.pool[word]@actual)
        for _ in range(steps):
            if formula==1:
                for unit in range(len(offsets)-1):apply(unit,1/steps)
            elif formula==2:
                for unit in range(len(offsets)-2):apply(unit,1/(2*steps))
                apply(len(offsets)-2,1/steps)
                for unit in reversed(range(len(offsets)-2)):apply(unit,1/(2*steps))
            else:raise ValueError('Unsupported formula')
        return actual


def metrics(exact,actual,model):
    norm=float(np.linalg.norm(actual));a=actual/norm;z=exact/np.linalg.norm(exact)
    overlap=np.vdot(z,a);residual=a-overlap*z;loss=float(sum(abs(residual)**2))
    n=model['n'];hf=(1<<model['electrons'])-1
    target=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)]
    number=joint=0.
    for v in range(len(a)):
        if v.bit_count()!=model['electrons']:number+=abs(a[v])**2
        if [sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=target:joint+=abs(a[v])**2
    return dict(error=loss/(1+abs(overlap)),infidelity=loss,norm=norm,number_leakage=float(number),joint_sector_leakage=float(joint))


def audit_ledger(record,unit_count):
    if record['arm']=='current_placement':return
    if record['cost_timing']!='through run_arm return':raise ValueError('Incomplete construction timing boundary')
    if record['return_time_CPU_seconds']<0:raise ValueError('Negative return/teardown time')
    if abs(record['internal_snapshot_CPU_seconds']+record['return_time_CPU_seconds']-record['actual_local_CPU_seconds'])>1e-12:
        raise ValueError('Through-return timing does not close')
    total=record['common_setup_CPU_charge']+record['inherited_placement_CPU_charge']+record['actual_local_CPU_seconds']
    if abs(total-record['total_charged_CPU_seconds'])>1e-12:raise ValueError('Construction CPU charge does not close')
    local=record['initializer_CPU_seconds']+record['setup_CPU_seconds']+record['search_CPU_seconds']+record['accounting_overhead_CPU_seconds']
    if abs(local-record['actual_local_CPU_seconds'])>1e-12:raise ValueError('Local CPU phase accounting does not close')
    if record['within_CPU_budget']!=(total<=record['CPU_budget']):raise ValueError('CPU feasibility mislabelled')
    if record['initialization']=='cold' and record['inherited_placement_CPU_charge']!=0:raise ValueError('Cold arm has warm charge')
    if record['initialization']=='placement' and record['inherited_placement_CPU_charge']<=0:raise ValueError('Uncharged warm initializer')
    if record['evaluations']>record['evaluation_cap']:raise ValueError('Evaluation cap exceeded')
    if 'starts' in record:
        if len(record['starts'])!=2 or sum(s['evaluations'] for s in record['starts'])!=record['evaluations']:
            raise ValueError('Start work accounting differs')
        for start in record['starts']:
            if sum(s['evaluations'] for s in start['sweeps'])!=start['evaluations']:raise ValueError('Sweep work omitted')
            if len(start['accepted'])>10:raise ValueError('Move cap exceeded')
            for sweep in start['sweeps']:
                if sweep['complete'] and sweep['evaluations']!=unit_count*(unit_count-1):raise ValueError('Incomplete sweep called complete')
                if not 0<=sweep['evaluations']<=unit_count*(unit_count-1):raise ValueError('Invalid examined positions')


def verify(output):
    from study import bind,sha,write,require,ARMS,CPU_BUDGETS,REPETITIONS
    freeze=bind(output);artifacts={};count=0;objective_errors=[];physical_errors=[];max_state=max_metric=0.
    for path in sorted((output/'designs').glob('*.json')):
        design=json.loads(path.read_text());model=json.loads((output/'inputs'/(design['case']+'.json')).read_text())
        reference=Reference(model);exact=reference.exact_action();records=design['records']
        require(len(records)==54,'Incomplete input search panel')
        require(len(design['repetitions'])==3,'Incomplete setup repetitions')
        conserving=[decode(b) for b in model['blocks']]
        pauli_words=sorted(decode(model['coefficients']))
        word_index={w:i for i,w in enumerate(pauli_words)}
        for setup in design['repetitions']:
            warm=setup['placement']
            require(sorted(warm['conserving'])==list(range(len(conserving))),'Placement is not a conserving permutation')
            flattened=[word_index[w] for i in warm['conserving'] for w in sorted(conserving[i])]
            require(warm['pauli']==flattened,'Pauli warm start is not the compatible complete-unit flattening')
            require(warm['charges']['conserving']==warm['placement_CPU_seconds'] and warm['charges']['pauli']==warm['cpu_seconds'],'Initializer charges do not bind measured phases')
        expected={(p+'_'+j+'_'+i,b,rep) for p,j,i in ARMS for b in CPU_BUDGETS for rep in range(REPETITIONS)}
        expected|={(name,b,rep) for name,b in (('pauli_JD_larger_CPU',8.),('current_placement',None)) for rep in range(REPETITIONS)}
        require({(r['arm'],r['CPU_budget'],r['repetition']) for r in records}==expected,'Missing factorial/control records')
        for index,record in enumerate(records):
            if record['arm']=='current_placement':expected_identity=('conserving','JD','current placement')
            elif record['arm']=='pauli_JD_larger_CPU':expected_identity=('pauli','JD','cold')
            else:expected_identity=tuple(record['arm'].split('_'))
            require(tuple(record[k] for k in ('partition','objective','initialization'))==expected_identity,'Named factorial/control identity changed')
            bs=([{w:c} for w,c in sorted(decode(model['coefficients']).items())]
                 if record['partition']=='pauli' else [decode(b) for b in model['blocks']])
            require(sorted(record['order'])==list(range(len(bs))),'Nonintegral/nonpermutation order')
            expanded=[(w,c) for i in record['order'] for w,c in sorted(bs[i].items())]
            supplied=[(tuple((q,p) for q,p in w),int(c)) for w,c in record['schedule']]
            require(expanded==supplied and dict(supplied)==decode(model['coefficients']),'Schedule changes target')
            require(record['unit_lengths']==[len(bs[i]) for i in record['order']],'Unit boundaries changed')
            setup=design['repetitions'][record['repetition']]
            require(setup['repetition']==record['repetition'],'Repetition/setup binding differs')
            warm=setup['placement'];common_charge=setup['common_setup_cpu_seconds']
            if record['arm']=='current_placement':
                require(record['order']==warm['conserving'],'Current baseline initializer changed')
                require(record['total_charged_CPU_seconds']==common_charge+warm['charges']['conserving'],'Current baseline CPU charge changed')
            else:
                expected_warm=warm['charges'][record['partition']] if record['initialization']=='placement' else 0.
                require(record['common_setup_CPU_charge']==common_charge and record['inherited_placement_CPU_charge']==expected_warm,'Unbound setup/initializer charge')
                l1=sorted(range(len(bs)),key=lambda i:(-sum(abs(c) for c in bs[i].values()),i))
                initial=warm[record['partition']] if record['initialization']=='placement' else list(range(len(bs)))
                if 'starts' in record:
                    require([x['initial_order'] for x in record['starts']]==[initial,l1],'Warm/cold starts changed')
                else:require(record['order']==initial,'Incomplete setup changed its fallback initializer')
            audit_ledger(record,len(bs));scores=reference.scores(record);score=scores[record['objective']]
            if record['score'] is not None:
                require(np.isclose(score,record['score'],rtol=2e-9,atol=5e-18),'Independent surrogate score differs')
            objective_errors.append(dict(case=design['case'],record_index=index,validation_score=score,
                validation_JD=scores['JD'],validation_J0=scores['J0'],
                search_score=record['score'],search_status=record['status']))
            for formula in (1,2):
                label=design['case']+'.'+str(index)+'.f'+str(formula);rp=output/'rows'/(label+'.json')
                row=json.loads(rp.read_text());lengths=record['unit_lengths'];cap=10*len(model['coefficients'])
                require(row['design_sha256']==sha(path),'Physical row uses a different search design')
                if formula==1:steps=cap//sum(lengths);rotations=steps*sum(lengths)
                elif len(lengths)==1:steps=1;rotations=sum(lengths)
                else:
                    steps=(cap-lengths[0])//(2*sum(lengths)-lengths[0]-lengths[-1])
                    rotations=steps*(2*sum(lengths)-lengths[0]-lengths[-1])+lengths[0]
                expected_row=(design['case'],index,record['arm'],record['repetition'],record['CPU_budget'],formula,steps,rotations,cap)
                require(tuple(row[k] for k in ('case','record_index','arm','repetition','CPU_budget','formula','steps','rotations','rotation_cap'))==expected_row,'Physical comparator identity differs')
                require(sha(output/row['state_file'])==row['state_sha256'],'Retained state changed')
                with np.load(output/row['state_file']) as data:actual=data['actual']
                other=reference.replay(record,steps,formula);delta=float(np.max(abs(actual-other)));max_state=max(max_state,delta)
                require(delta<2e-10,'Independent unfused propagation differs')
                values=metrics(exact,other,model)
                physical_errors.append(dict(case=design['case'],repetition=record['repetition'],formula=formula,
                    record_index=index,producer_error=row['error'],independent_error=values['error']))
                for key,value in values.items():
                    delta=abs(row[key]-value);max_metric=max(max_metric,delta)
                    require(delta<2e-12,'Independent physical metric differs')
                artifacts[str(rp.relative_to(output))]=sha(rp);artifacts[row['state_file']]=sha(output/row['state_file']);count+=1
        print('VERIFIED',design['case'],flush=True)
    require(count==864,'Incomplete physical audit')
    ambiguities=[];classification=lambda d:'loss' if d>1e-12 else 'win' if d < -1e-12 else 'tie'
    for row in physical_errors:
        if (row['producer_error']>1e-12)!=(row['independent_error']>1e-12):
            ambiguities.append(dict(kind='floor',row=row))
    groups={}
    for row in physical_errors:groups.setdefault((row['case'],row['repetition']),[]).append(row)
    for key,items in groups.items():
        for at,a in enumerate(items):
            for b in items[at+1:]:
                if classification(a['producer_error']-b['producer_error'])!=classification(a['independent_error']-b['independent_error']):
                    ambiguities.append(dict(kind='win_tie_loss',group=list(key),records=[[a['record_index'],a['formula']],[b['record_index'],b['formula']]]))
    write(output/'INDEPENDENT_ERRORS.json',dict(rows=physical_errors,decision_agreement=not ambiguities,ambiguities=ambiguities))
    artifacts['INDEPENDENT_ERRORS.json']=sha(output/'INDEPENDENT_ERRORS.json')
    write(output/'VALIDATION_SCORES.json',objective_errors);artifacts['VALIDATION_SCORES.json']=sha(output/'VALIDATION_SCORES.json')
    write(output/'VERIFIED.json',dict(status='passed' if not ambiguities else 'numerically_ambiguous',
        decision_agreement=not ambiguities,rows=count,artifacts=artifacts,max_state_delta=max_state,
        max_metric_delta=max_metric,schedule_freeze_sha256=sha(output/'SCHEDULES_FROZEN.json'),sources=freeze['sources']))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    verify(p.parse_args().output.resolve())
