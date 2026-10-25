"""Symmetric-step cubic-defect response and bounded insertion search.

No physical error or exact-Hamiltonian time propagation is used by the search.
The actual symmetric-step Taylor coefficient is evaluated on the HF column.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import importlib.util
import math
from pathlib import Path
import time
import numpy as np

OLD=Path(__file__).resolve().parents[1]/'prospective_reference'
spec=importlib.util.spec_from_file_location('symmetric_response_prior_engine',OLD/'engine.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
RESCORE_RESERVATION=.02
CANDIDATE_MARGIN=.002


def require(ok,message):
    if not ok:raise ValueError(message)


def decode(raw):
    out={}
    for word,c in raw:
        require(type(c) in (int,str),'Integer coefficient numerators required')
        w=tuple((q,p) for q,p in word);require(w not in out,'Repeated Pauli word')
        out[w]=int(c)
    return out


def validate_model(model):
    n=model['n'];require(type(n) is int and 0<n<=8,'Register bound before allocation')
    require(type(model['electrons']) is int and 0<=model['electrons']<=n,'Invalid HF occupation')
    require(0<len(model['blocks'])<=128 and 0<len(model['coefficients'])<=128,'Unit/term bound before allocation')
    require(type(model['denominator']) in (int,str) and int(model['denominator'])>0,'Positive integer denominator required')
    blocks=[decode(raw) for raw in model['blocks']];total=decode(model['coefficients'])
    require(all(blocks) and sum(map(len,blocks))==len(total) and {w:c for b in blocks for w,c in b.items()}==total,'Units do not partition target')
    for block in blocks:
        for word in block:
            require(len(dict(word))==len(word) and all(type(q) is int and 0<=q<n and p in 'XYZ' for q,p in word),'Invalid Pauli support')
        engine.descriptor(block)
    return blocks


def permutation(order,m):
    require(all(type(i) is int for i in order) and sorted(order)==list(range(m)),'Integral unit permutation required')


def symmetric_units(order):
    require(bool(order),'Nonempty symmetric step required')
    return [(i,.5) for i in order[:-1]]+[(order[-1],1.)]+[(i,.5) for i in reversed(order[:-1])]


class Interrupted(RuntimeError):
    def __init__(self,completed_stages):self.completed_stages=completed_stages


class Context:
    """Real one-flip unit actions and at most eight-mode Taylor arrays."""
    def __init__(self,model):
        blocks=validate_model(model);n=model['n'];dimension=1<<n
        require(len(blocks)*dimension*24<=8*1024**2,'Action-cache storage guard')
        self.model=model;self.blocks=blocks;self.m=len(blocks);self.dimension=dimension
        self.den=int(model['denominator']);self.v=(1<<model['electrons'])-1
        self.lengths=[len(b) for b in blocks];self.words=len(model['coefficients'])
        basis=np.arange(dimension);self.actions=[]
        for b in blocks:
            desc=engine.descriptor(b);permutation=basis^desc[0]
            values=np.array([engine.value(desc,int(v))/self.den for v in basis],float)
            require(np.all(np.isfinite(values)),'Nonfinite action coefficient')
            self.actions.append((permutation,values))
        D={w:c for b in blocks if engine.descriptor(b)[0]==0 for w,c in b.items()}
        diagonal=engine.descriptor(D or {():0})
        gap_numerators=[engine.value(diagonal,v)-engine.value(diagonal,self.v) for v in range(dimension)]
        self.gaps=np.array([g/self.den for g in gap_numerators],float)
        self.initial=np.zeros(dimension);self.initial[self.v]=1
        self.H1=self.hamiltonian_action(self.initial)
        self.H2=self.hamiltonian_action(self.H1)/2
        self.H3=self.hamiltonian_action(self.hamiltonian_action(self.H1))/6
        require(all(np.all(np.isfinite(a)) for a in (self.gaps,self.H1,self.H2,self.H3)),'Nonfinite Taylor reference')
        self.cache_bytes=sum(p.nbytes+v.nbytes for p,v in self.actions)+sum(a.nbytes for a in (self.gaps,self.initial,self.H1,self.H2,self.H3))
        require(self.cache_bytes<=8*1024**2,'Actual action cache exceeds bound')

    def apply(self,unit,state):
        ids,values=self.actions[unit]
        return values[ids]*state[ids]

    def hamiltonian_action(self,state):
        out=np.zeros_like(state)
        for unit in range(self.m):out+=self.apply(unit,state)
        return out

    def cubic(self,order,deadline=float('inf'),clock=time.process_time):
        permutation(order,self.m)
        q=[self.initial.copy()]+[np.zeros(self.dimension) for _ in range(3)]
        completed=0
        for unit,tau in symmetric_units(order):
            if clock()>=deadline:raise Interrupted(completed)
            out=[a.copy() for a in q];powers=list(q);factor=1.
            for degree in range(1,4):
                factor*=tau/degree
                for k in range(4-degree):
                    powers[k]=self.apply(unit,powers[k]);out[k+degree]+=factor*powers[k]
            q=out;completed+=1
        require(np.max(abs(q[1]-self.H1))<2e-11 and np.max(abs(q[2]-self.H2))<2e-11,'Symmetric lower Taylor coefficients differ')
        return q[3]-self.H3

    def score(self,order,T,rotation_cap,deadline=float('inf'),clock=time.process_time):
        require(math.isfinite(T) and T>0 and type(rotation_cap) is int and rotation_cap>0,'Invalid time/rotation cap')
        coefficient=self.cubic(order,deadline,clock)
        weights=T*T*np.sinc(T*self.gaps/(2*np.pi))**2
        weights[self.v]=0
        JD2=float(np.dot(weights,coefficient*coefficient))
        require(math.isfinite(JD2) and JD2>=0,'Nonfinite/negative cubic response score')
        lengths=[self.lengths[i] for i in order];steps=engine.feasible_steps(lengths,rotation_cap,True)
        require(steps>=1,'No affordable symmetric step')
        return dict(order=list(order),JD2=JD2,score=JD2/(steps**4),steps=steps,
            cubic_coefficients=coefficient.tolist())


class FirstOrderContext:
    """The established first-order response at the same declared physical T.

    This adapter uses the same bounded insertion policy as the new objective.
    Its setup and initializers are recorded separately by the study controller.
    """
    def __init__(self,model,T,partition):
        blocks=validate_model(model);total=decode(model['coefficients'])
        require(partition in ('conserving','pauli') and math.isfinite(T) and T>0,'Invalid first-order control')
        self.blocks=blocks if partition=='conserving' else [{w:c} for w,c in sorted(total.items())]
        self.m=len(self.blocks);self.v=(1<<model['electrons'])-1;self.T=T
        diagonal={w:c for b in blocks if engine.descriptor(b)[0]==0 for w,c in b.items()}
        self.features=engine.Features(self.blocks,engine.Reference(diagonal or {():0},self.v,T,int(model['denominator'])))
        self.lengths=[len(b) for b in self.blocks]

    def score(self,order,T,rotation_cap,deadline=float('inf'),clock=time.process_time):
        require(T==self.T,'Control reference time changed')
        if clock()>=deadline:raise Interrupted(0)
        permutation(order,self.m);value=self.features.score(self.features.amplitude(order))
        steps=rotation_cap//sum(self.lengths)
        require(steps>=1,'No affordable first-order step')
        # The first-order count is order independent; this positive common
        # normalization cannot change its within-cell ordering preference.
        return dict(order=list(order),JD1=value,score=value/(steps**2),steps=steps)


def insertion(context,seeds,T,rotation_cap,CPU_seconds=2.,evaluation_cap=1000,max_moves=5,clock=time.process_time):
    """Three fixed starts; complete evaluation of every visited candidate.

    Initial seed scores belong to method setup and are passed explicitly.
    Equal per-start CPU slices do not carry unused time to the next start.
    A best improving interrupted-sweep candidate may be independently rescored
    and accepted. Actual incomplete/overrun work remains in the returned ledger.
    """
    require(len(seeds)==3 and 0<CPU_seconds<=8 and 0<=evaluation_cap<=325120 and 0<=max_moves<=10,'Bounded declared search policy required')
    for seed in seeds:permutation(seed['order'],context.m)
    began=clock();global_deadline=began+CPU_seconds;starts=[];total_attempted=total_completed=0;rescores=0;attempted_rescores=0
    for index,seed in enumerate(seeds):
        started=clock();deadline=min(global_deadline,started+CPU_seconds/3)
        limit=evaluation_cap//3+(index<evaluation_cap%3)
        current=dict(seed);initial_score=current['score'];attempted=completed=0;accepted=[];sweeps=[];termination='move cap'
        for sweep in range(max_moves):
            best=current;move=None;visited=done=0;complete=True;reason='local optimum';partial_stages=None
            for source in range(context.m):
                for destination in range(context.m):
                    if source==destination:continue
                    if attempted>=limit or clock()+RESCORE_RESERVATION+CANDIDATE_MARGIN>=deadline:
                        complete=False;reason='evaluation cap' if attempted>=limit else 'CPU/rescore reservation';break
                    order=list(current['order']);unit=order.pop(source);order.insert(destination,unit)
                    attempted+=1;visited+=1
                    try:trial=context.score(order,T,rotation_cap,deadline-RESCORE_RESERVATION,clock)
                    except Interrupted as interrupted:
                        complete=False;reason='CPU during candidate';partial_stages=interrupted.completed_stages;break
                    completed+=1;done+=1
                    if trial['score']<best['score']-max(1e-30,1e-12*abs(best['score'])):
                        best=trial;move=[source,destination]
                if not complete:break
            sweeps.append(dict(sweep=sweep,attempted=visited,completed=done,complete=complete,reason=reason,partial_candidate_stages=partial_stages))
            if move is None:termination=reason;break
            if clock()+CANDIDATE_MARGIN>=deadline:termination='CPU before accepted rescore';break
            attempted_rescores+=1
            try:actual=context.score(best['order'],T,rotation_cap,deadline,clock)
            except Interrupted:termination='CPU during accepted rescore';break
            rescores+=1
            if actual['score']>=current['score']-max(1e-30,1e-12*abs(current['score'])):
                termination='roundoff rejected move';break
            current=actual;accepted.append(dict(move=move,evaluations=completed,record=dict(current)))
            if not complete:termination=reason;break
        starts.append(dict(initial_order=seed['order'],initial_score=initial_score,record=current,
            attempted=attempted,completed=completed,accepted=accepted,sweeps=sweeps,termination=termination))
        total_attempted+=attempted;total_completed+=completed
    best=min([s['record'] for s in starts],key=lambda r:(r['score'],tuple(r['order'])))
    return dict(record=best,starts=starts,attempted=total_attempted,completed=total_completed,
        evaluation_cap=evaluation_cap,attempted_rescores=attempted_rescores,accepted_rescores=rescores,nominal_search_CPU_seconds=CPU_seconds,
        max_moves_per_start=max_moves,stopping_policy='Three equal CPU slices,no unused carry;explicit evaluation/move limits;actualthrough-returntiming suppliedbycaller')
