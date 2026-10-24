"""CPU/evaluation-capped response insertion with explicit initialization charges."""
import importlib.util
from pathlib import Path
import time
import numpy as np

OLD=Path(__file__).resolve().parents[1]/'prospective_reference'
spec=importlib.util.spec_from_file_location('cost_study_determinant_engine',OLD/'engine.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
RESCORE_RESERVATION=.1
CANDIDATE_MARGIN=.002


class Exhausted(RuntimeError):
    def __init__(self,phase,pairs=0):self.phase=phase;self.pairs=pairs


class Features:
    """Sparse one-determinant response; no 2^n array or dense Hamiltonian."""
    amplitude=engine.Features.amplitude
    crossing=engine.Features.crossing
    score=staticmethod(engine.Features.score)

    def __init__(self,blocks,reference,deadline,clock=time.process_time):
        m=len(blocks)
        if not 0<m<=256 or m*(m-1)//2*1024>64*1024**2:
            raise ValueError('Conservative pair/container storage guard exceeded')
        self.storage_upper_bound_bytes=m*(m-1)//2*1024
        if clock()>=deadline:raise Exhausted('before feature setup')
        desc=[engine.descriptor(b) for b in blocks];maps={};destinations=set();visited=0
        for i in range(m):
            for j in range(i+1,m):
                if clock()>=deadline:raise Exhausted('pair feature setup',visited)
                visited+=1;x,y=desc[i][0],desc[j][0];v=reference.v
                b=(engine.value(desc[j],v)*engine.value(desc[i],v^y)
                   -engine.value(desc[i],v)*engine.value(desc[j],v^x))
                if not b:continue
                response=reference.response([{v^x^y:-.5j*b/(reference.den**2)}])
                if response:maps[i,j]=response;destinations.update(response)
        self.destinations=sorted(destinations);columns={w:k for k,w in enumerate(self.destinations)}
        self.pairs={}
        for pair,response in maps.items():
            if clock()>=deadline:raise Exhausted('sparse feature indexing',visited)
            self.pairs[pair]=(np.array([columns[w] for w in response],np.int32),np.array(list(response.values()),complex))
        self.m=m;self.visited_pairs=visited
        self.array_payload_bytes=sum(ids.nbytes+a.nbytes for ids,a in self.pairs.values())


def insertion(starts,features,deadline,evaluation_cap,max_moves=10,clock=time.process_time):
    if len(starts)!=2 or not 0<=evaluation_cap<=2_000_000 or not 0<=max_moves<=10:
        raise ValueError('Bounded two-start search required')
    remaining=max(0.,deadline-clock());per_start=remaining/2
    results=[];used=0;rescores=[]
    for start_index,start in enumerate(starts):
        order=list(start);cap=evaluation_cap//2+(start_index<evaluation_cap%2)
        local_deadline=min(deadline,clock()+per_start)
        if clock()+RESCORE_RESERVATION>=local_deadline:
            results.append(dict(order=order,initial_order=list(start),score=None,evaluations=0,accepted=[],sweeps=[],termination='insufficient initial-score reservation'))
            continue
        before=clock();a=features.amplitude(order);score=features.score(a);rescores.append(clock()-before)
        initial_score=score;count=0;accepted=[];sweeps=[];termination='move cap'
        for sweep in range(max_moves):
            best=score;move=None;visited=0;complete=True;reason='local optimum'
            for source in range(len(order)):
                for direction in (-1,1):
                    trial=a.copy();at=source
                    while 0<=at+direction<len(order):
                        if count>=cap or clock()+RESCORE_RESERVATION+CANDIDATE_MARGIN>=local_deadline:
                            complete=False;reason='evaluation cap' if count>=cap else 'CPU/rescore reservation'
                            break
                        crossed=order[at+direction]
                        left,right=((crossed,order[source]) if direction<0 else (order[source],crossed))
                        features.crossing(trial,left,right);at+=direction;count+=1;visited+=1
                        candidate=features.score(trial)
                        if candidate<best-max(1e-30,1e-12*best):best=candidate;move=(source,at)
                    if not complete:break
                if not complete:break
            sweeps.append(dict(sweep=sweep,evaluations=visited,complete=complete,reason=reason))
            if move is None:termination=reason;break
            if clock()+RESCORE_RESERVATION>=local_deadline:
                termination='CPU/rescore reservation before accepted move';break
            candidate=list(order);item=candidate.pop(move[0]);candidate.insert(move[1],item)
            before=clock();actual=features.amplitude(candidate);actual_score=features.score(actual);rescores.append(clock()-before)
            if actual_score>=score-max(1e-30,1e-12*score):termination='roundoff rejected move';break
            order,a,score=candidate,actual,actual_score
            accepted.append(dict(move=list(move),score=score,evaluations=count))
            if not complete:termination=reason;break
        used+=count
        results.append(dict(order=order,initial_order=list(start),score=score,initial_score=initial_score,evaluations=count,
            accepted=accepted,sweeps=sweeps,termination=termination))
    valid=[r for r in results if r['score'] is not None]
    chosen=min(valid,key=lambda r:(r['score'],tuple(r['order']))) if valid else results[0]
    return dict(order=chosen['order'],score=chosen['score'],evaluations=used,starts=results,
        full_rescore_CPU_seconds=rescores,rescore_reservation_seconds=RESCORE_RESERVATION,
        rescore_reservation_violated=any(x>RESCORE_RESERVATION for x in rescores))


def run_arm(blocks,D,v,den,objective,initialization,warm_order,CPU_budget,common_charge,
            warm_charge,evaluation_cap,clock=time.process_time):
    if objective not in ('JD','J0') or initialization not in ('cold','placement') or not 0<CPU_budget<=8:
        raise ValueError('Unsupported benchmark arm/budget')
    if not 0<len(blocks)<=256 or not np.all(np.isfinite([common_charge,warm_charge])) or min(common_charge,warm_charge)<0:
        raise ValueError('Invalid size or setup/initializer charge')
    began=clock();charged_warm=warm_charge if initialization=='placement' else 0.
    deadline=began+CPU_budget-common_charge-charged_warm
    init_began=clock()
    l1=sorted(range(len(blocks)),key=lambda i:(-sum(abs(c) for c in blocks[i].values()),i))
    starts=[list(range(len(blocks))) if initialization=='cold' else list(warm_order),l1]
    if any(sorted(s)!=list(range(len(blocks))) for s in starts):raise ValueError('Initializer is not a unit permutation')
    initialization_CPU=clock()-init_began
    ledger=dict(CPU_budget=CPU_budget,common_setup_CPU_charge=common_charge,
        inherited_placement_CPU_charge=charged_warm,initializer_CPU_seconds=initialization_CPU,
        objective=objective,initialization=initialization,evaluation_cap=evaluation_cap,
        status='completed',setup_CPU_seconds=0.,search_CPU_seconds=0.,evaluations=0)
    if clock()>=deadline:
        ledger.update(status='initializer charge exhausts CPU budget',order=starts[0],score=None)
    else:
        setup=clock()
        try:
            reference=engine.Reference(D if objective=='JD' else {():0},v,1.,den)
            features=Features(blocks,reference,deadline,clock)
        except Exhausted as exhausted:
            ledger.update(status='incomplete '+exhausted.phase,setup_pairs_visited=exhausted.pairs,
                          order=starts[0],score=None,setup_CPU_seconds=clock()-setup)
        else:
            ledger.update(setup_CPU_seconds=clock()-setup,setup_pairs_visited=features.visited_pairs,
                nonzero_pair_count=len(features.pairs),destinations=len(features.destinations),
                array_payload_bytes=features.array_payload_bytes,storage_upper_bound_bytes=features.storage_upper_bound_bytes)
            search=clock();result=insertion(starts,features,deadline,evaluation_cap,clock=clock)
            ledger.update(result,search_CPU_seconds=clock()-search)
            if result['score'] is None:ledger['status']='insufficient search-score reservation'
    ledger['actual_local_CPU_seconds']=clock()-began
    ledger['accounting_overhead_CPU_seconds']=ledger['actual_local_CPU_seconds']-(
        ledger['initializer_CPU_seconds']+ledger['setup_CPU_seconds']+ledger['search_CPU_seconds'])
    ledger['total_charged_CPU_seconds']=common_charge+charged_warm+ledger['actual_local_CPU_seconds']
    ledger['CPU_overrun_seconds']=max(0.,ledger['total_charged_CPU_seconds']-CPU_budget)
    ledger['within_CPU_budget']=ledger['CPU_overrun_seconds']==0
    if not ledger['within_CPU_budget']:ledger['status']='CPU overrun retained: '+ledger['status']
    return ledger
