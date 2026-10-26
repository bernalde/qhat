"""Independent objective, search-ledger and full-space physical verification."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
from independent import Reference,require


def metric(exact,actual,model):
    norm=float(np.sqrt(sum(abs(actual)**2)));a=actual/norm;z=exact/np.sqrt(sum(abs(exact)**2));inner=np.vdot(z,a)
    residual=a-z*inner;loss=float(sum(abs(residual)**2));n=model['n'];hf=(1<<model['electrons'])-1
    sector=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)];number=joint=0.
    for v in range(len(a)):
        if v.bit_count()!=model['electrons']:number+=abs(a[v])**2
        if [sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=sector:joint+=abs(a[v])**2
    return dict(error=loss/(1+abs(inner)),infidelity=loss,norm=norm,number_leakage=float(number),joint_sector_leakage=float(joint))


def time_ledger(record):
    require(record['cost_timing']=='through constructor return','Incomplete constructor timing')
    pieces=list(record['phases'].values())+[record['bookkeeping_CPU_seconds']]
    require(all(np.isfinite(x) and x>=-1e-12 for x in pieces),'Invalid construction phase cost')
    require(abs(sum(pieces)-record['actual_local_CPU_seconds'])<1e-12,'Constructor phase accounting does not close')


def affordability(lengths,cap,formula):
    U=sum(lengths)
    if formula==1:return cap//U,(cap//U)*U
    if len(lengths)==1:return 1,U
    denominator=2*U-lengths[0]-lengths[-1]
    steps=(cap-lengths[0])//denominator
    return steps,steps*denominator+lengths[0]


def audit_search(result,seeds,score,m,CPU,evaluations,moves):
    require(result['nominal_search_CPU_seconds']==CPU and result['evaluation_cap']==evaluations,'Search allowance changed')
    require(result['max_moves_per_start']==moves and len(result['starts'])==3,'Search starts/moves changed')
    attempted=completed=accepted_count=0
    for index,(start,seed) in enumerate(zip(result['starts'],seeds)):
        score(seed)
        require(start['initial_order']==seed['order'] and start['initial_score']==seed['score'],'Initial seed changed')
        limit=evaluations//3+(index<evaluations%3)
        require(0<=start['completed']<=start['attempted']<=limit,'Per-start work allowance exceeded')
        require(len(start['accepted'])<=moves and len(start['sweeps'])<=moves,'Move/sweep cap exceeded')
        require(sum(s['attempted'] for s in start['sweeps'])==start['attempted'],'Attempted work does not close')
        require(sum(s['completed'] for s in start['sweeps'])==start['completed'],'Completed work does not close')
        for sweep in start['sweeps']:
            require(0<=sweep['completed']<=sweep['attempted']<=m*(m-1),'Invalid sweep work')
            require(sweep['attempted']-sweep['completed']<=1,'Multiple uncompleted candidates in one sweep')
            require(not sweep['complete'] or sweep['attempted']==sweep['completed']==m*(m-1),'False full-sweep label')
        interrupted_rescore=start.get('interrupted_rescore_completed_stages')
        if start['termination']=='CPU during accepted rescore':
            require(type(interrupted_rescore) is int and 0<=interrupted_rescore<2*m-1,
                'Interrupted rescore stage detail missing/invalid')
        else:require(interrupted_rescore is None,'False interrupted-rescore stage label')
        previous=seed;last_evaluations=0
        for accepted in start['accepted']:
            current=accepted['record'];score(current);move=accepted['move']
            require(len(move)==2 and all(type(x) is int and 0<=x<m for x in move) and move[0]!=move[1],'Invalid insertion move')
            order=list(previous['order']);unit=order.pop(move[0]);order.insert(move[1],unit)
            require(order==current['order'],'Accepted order is inconsistent with insertion')
            require(current['score']<previous['score']-max(1e-30,1e-12*abs(previous['score'])),'Nonimproving accepted move')
            require(last_evaluations<=accepted['evaluations']<=start['completed'],'Invalid accepted work count')
            last_evaluations=accepted['evaluations'];previous=current;accepted_count+=1
        require(start['record']==previous,'Final start differs from retained accepted history')
        score(start['record']);attempted+=start['attempted'];completed+=start['completed']
    require((result['attempted'],result['completed'])==(attempted,completed),'Total candidate work differs')
    require(accepted_count<=result['accepted_rescores']<=result['attempted_rescores']<=accepted_count+3,'Rescore accounting differs')
    best=min((s['record'] for s in result['starts']),key=lambda x:(x['score'],tuple(x['order'])))
    require(result['record']==best,'Final search winner differs')


def verify(output):
    from study import bind,old_model,ARMS,sha,write,sources
    protocol=bind(output,allow_archived_sources=True);freeze=json.loads((output/'SCHEDULES_FROZEN.json').read_text())
    read=lambda p:json.loads(p.read_text());all_rows=[];objective_records=[];max_state=max_metric=max_cubic=0.
    expected_cases=set(protocol['cases']);observed_cases=set()
    for path in sorted((output/'designs').glob('*.json')):
        design=read(path);model=read(output/'inputs'/path.name);case=model['case']
        original_case=case.removeprefix('preflight_altered_');observed_cases.add(original_case)
        original,digest=old_model(original_case)
        if protocol['preflight']:
            original['case']='preflight_altered_'+original_case
            original['coefficients']=[[w,str(7*int(c))] for w,c in original['coefficients']]
            original['blocks']=[[[w,str(7*int(c))] for w,c in b] for b in original['blocks']]
            original['denominator']=str(5*int(original['denominator']))
            if 'removed_scalar_exact' in original:
                original['removed_scalar_exact']=str(Fraction(original['removed_scalar_exact'])*Fraction(7,5))
        require(model==original and design['input_sha256']==sha(output/'inputs'/path.name),'Authenticated input changed')
        require(digest==freeze['prior_models'][original_case],'Prior input binding changed')
        reference=Reference(model);reference.conservation();cache={}
        def checked(record,T,cap,partition='conserving',kind='cubic'):
            nonlocal max_cubic
            order=record['order'];key=(tuple(order),T,cap,partition,kind)
            if key not in cache:
                other=(reference.score(order,T,cap) if kind=='cubic' else reference.first_score(order,T,cap,partition))
                cache[key]=other
            other=cache[key]
            for field in ('JD2','score') if kind=='cubic' else ('JD1','score'):
                require(np.isclose(record[field],other[field],rtol=2e-9,atol=5e-18),'Independent objective differs: '+field)
            require(record['steps']==other['steps'],'Objective affordable-step count differs')
            if kind=='cubic':
                delta=float(np.max(abs(np.array(record['cubic_coefficients'])-other['cubic_coefficients'])))
                max_cubic=max(max_cubic,delta);require(delta<2e-11,'Independent full cubic coefficient differs')
        construction=design['construction'];require(construction['cubic_cache_bytes']<=8*1024**2,'Cubic cache exceeds bound')
        require(all(np.isfinite(v) and v>=0 for v in construction.values()),'Invalid shared cost')
        require(construction['legacy_current_reference_time']==1.,'Legacy baseline relabelled as same-time optimization')
        base=sum(construction[k] for k in ('model_load_CPU_seconds','base_setup_CPU_seconds','base_cleanup_CPU_seconds'))
        current_charge=base+construction['legacy_current_CPU_seconds']
        cache_charge=construction['cubic_cache_setup_CPU_seconds']+construction['cubic_cache_cleanup_CPU_seconds']
        require([r['T'] for r in design['time_records']]==protocol['times'],'Time coverage differs')
        require(len(design['records'])==len(protocol['times'])*2*8,'Missing circuit comparators')
        for tindex,tr in enumerate(design['time_records']):
            T=tr['T'];current=design['current']['order'];first_order=tr['first']['result']['record']['order']
            for partition,key in (('conserving','first'),('pauli','pauli')):
                record=tr[key];time_ledger(record);m=len(reference.blocks) if partition=='conserving' else len(reference.paulis)
                require(record['partition']==partition and record['T']==T,'First-order control identity changed')
                bs=([{tuple((q,p) for q,p in w):int(c) for w,c in raw} for raw in model['blocks']] if partition=='conserving' else
                    [{w:int(c)} for w,c in sorted((tuple((q,p) for q,p in w),c) for w,c in model['coefficients'])])
                l1=sorted(range(m),key=lambda i:(-sum(abs(c) for c in bs[i].values()),i))
                warm=list(current)
                if partition=='pauli':
                    index={next(iter(b)):i for i,b in enumerate(bs)}
                    warm=[index[w] for i in current for w in sorted({tuple((q,p) for q,p in w):int(c) for w,c in model['blocks'][i]})]
                require([r['order'] for r in record['seeds']]==[list(range(m)),l1,warm],'Strong-control initializers differ')
                score=lambda r:checked(r,T,10*len(model['coefficients']),partition,'first')
                audit_search(record['result'],record['seeds'],score,m,2. if partition=='conserving' else 8.,1000 if partition=='conserving' else 20*m*(m-1),5 if partition=='conserving' else 10)
            require([c['budget_units'] for c in tr['cells']]==protocol['budgets'],'Budget coverage differs')
            for bindex,cell in enumerate(tr['cells']):
                cap=cell['budget_units']*len(model['coefficients']);new=cell['cubic'];time_ledger(new)
                require(cell['first_record_index']==8*(2*tindex+bindex),'Circuit cell indexing differs')
                require([s['order'] for s in new['seeds']]==[first_order,list(reversed(first_order)),current],'Cubic seed policy changed')
                score=lambda r:checked(r,T,cap)
                for seed in new['seeds']:score(seed)
                require(new['orientation']==min(new['seeds'],key=lambda x:(x['score'],tuple(x['order']))),'Orientation rule changed')
                require(abs(new['orientation_constructor_CPU_seconds']-sum(v for k,v in new['phases'].items() if k!='search_CPU_seconds'))<1e-12,'Orientation cost does not close')
                audit_search(new['result'],new['seeds'],score,len(reference.blocks),2.,1000,5)
                orders={'current':current,'first_JD_sameT':first_order,'pauli_JD_large_sameT':tr['pauli']['result']['record']['order'],
                    'cubic_orientation':new['orientation']['order'],'cubic_insertion':new['result']['record']['order']}
                first_charge=current_charge+tr['first']['actual_local_CPU_seconds']
                charges={'current':current_charge,'first_JD_sameT':first_charge,'pauli_JD_large_sameT':current_charge+tr['pauli']['actual_local_CPU_seconds'],
                    'cubic_orientation':first_charge+cache_charge+new['orientation_constructor_CPU_seconds'],
                    'cubic_insertion':first_charge+cache_charge+new['actual_local_CPU_seconds']}
                records=design['records'][cell['first_record_index']:cell['first_record_index']+8]
                require([(r['arm'],r['formula']) for r in records]==list(ARMS),'Comparator identities changed')
                for rec in records:
                    require(rec['order']==orders[rec['arm']] and rec['T']==T and rec['budget_units']==cell['budget_units'],'Comparator order/time/budget changed')
                    require(abs(rec['standalone_construction_CPU_charge']-charges[rec['arm']])<1e-12,'Standalone cost attribution differs')
        references={T:reference.exact_action(T) for T in protocol['times']}
        for index,rec in enumerate(design['records']):
            partition='pauli' if rec['arm']=='pauli_JD_large_sameT' else 'conserving';require(rec['partition']==partition,'Partition changed')
            bs=([{tuple((q,p) for q,p in w):int(c) for w,c in raw} for raw in model['blocks']] if partition=='conserving' else
                [{w:int(c)} for w,c in sorted((tuple((q,p) for q,p in w),c) for w,c in model['coefficients'])])
            require(all(type(i) is int for i in rec['order']) and sorted(rec['order'])==list(range(len(bs))),'Invalid physical permutation')
            lengths=[len(bs[i]) for i in rec['order']];schedule=[[[[q,p] for q,p in w],str(c)] for i in rec['order'] for w,c in sorted(bs[i].items())]
            require(rec['unit_lengths']==lengths and rec['schedule']==schedule,'Flattened complete-unit circuit differs')
            label=case+'.'+str(index);row=read(output/'rows'/(label+'.json'));sp=output/row['state_file']
            require(row['design_sha256']==sha(path) and row['state_sha256']==sha(sp),'Physical artifact binding differs')
            for key in ('arm','formula','T','budget_units','standalone_construction_CPU_charge'):require(row[key]==rec[key],'Physical row identity changed')
            require(row['case']==case and row['record_index']==index,'Physical row index differs')
            cap=rec['budget_units']*len(model['coefficients']);steps,rotations=affordability(lengths,cap,rec['formula'])
            require(steps>=1 and (row['rotation_cap'],row['steps'],row['rotations'])==(cap,steps,rotations),'Actual circuit affordability differs')
            actual=reference.replay(rec,steps,rec['T'])
            with np.load(sp) as data:producer=data['actual']
            delta=float(np.max(abs(actual-producer)));max_state=max(max_state,delta);require(delta<2e-10,'Unfused full-space replay differs')
            measured=metric(references[rec['T']],actual,model)
            for key,value in measured.items():
                delta=abs(value-row[key]);max_metric=max(max_metric,delta);require(delta<2e-12,'Independent physical metric differs')
            require((measured['error']<=1e-12)==(row['error']<=1e-12),'Error floor decision differs')
            all_rows.append(dict(case=case,record_index=index,T=rec['T'],budget_units=rec['budget_units'],arm=rec['arm'],formula=rec['formula'],
                error=row['error'],independent_error=measured['error']))
        objective_records.append(dict(case=case,distinct_objectives=len(cache),status='passed'))
        print('VERIFIED_SYMMETRIC',case,len(design['records']),'circuits',flush=True)
    require(observed_cases==expected_cases and len(all_rows)==protocol['expected_circuits'],'Incomplete independent grid')
    require(len(list((output/'rows').glob('*.json')))==len(all_rows)==len(list((output/'states').glob('*.npz'))),'Extra/missing physical artifact')
    classify=lambda x:'win' if x < -1e-12 else 'loss' if x>1e-12 else 'tie';decisions=[]
    for at,a in enumerate(all_rows):
        for c in all_rows[at+1:]:
            if (a['case'],a['T'],a['budget_units'])!=(c['case'],c['T'],c['budget_units']):continue
            require(classify(a['error']-c['error'])==classify(a['independent_error']-c['independent_error']),'Pairwise accuracy decision differs')
            require((min(a['error'],c['error'])<=1e-12)==(min(a['independent_error'],c['independent_error'])<=1e-12),'Paired ratio-floor decision differs')
            decisions.append(dict(case=a['case'],T=a['T'],budget_units=a['budget_units'],left=a['record_index'],right=c['record_index'],classification=classify(a['error']-c['error'])))
    require(len(decisions)==28*len(expected_cases)*len(protocol['times'])*2,'Incomplete pairwise decisions')
    write(output/'INDEPENDENT_ERRORS.json',dict(rows=all_rows,decisions=decisions,decision_agreement=True))
    write(output/'INDEPENDENT_OBJECTIVES.json',objective_records)
    artifacts={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='VERIFIED.json'}
    write(output/'VERIFIED.json',dict(status='passed',rows=len(all_rows),decision_agreement=True,comparisons=len(decisions),
        max_state_delta=max_state,max_metric_delta=max_metric,max_cubic_column_delta=max_cubic,
        schedule_freeze_sha256=sha(output/'SCHEDULES_FROZEN.json'),sources=freeze['sources'],
        verification_sources=sources(),artifacts=artifacts))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    verify(p.parse_args().output.resolve())
