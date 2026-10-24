"""Exact rational optimization of one fixed contiguous response-order window.

No Hamiltonian propagation is imported here. The finite domain and supplied
binary64 sinc-squared weights are frozen by study.py before enumeration.
"""
import itertools
import math
import time
from pathlib import Path
import importlib.util
import numpy as np

OLD=Path(__file__).resolve().parents[1]/'prospective_reference'
spec=importlib.util.spec_from_file_location('exact_window_prior_engine',OLD/'engine.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)


def require(condition,message):
    if not condition:raise ValueError(message)


def decode(raw):return {tuple((int(q),p) for q,p in w):int(c) for w,c in raw}


def validate_model(model):
    require(type(model['n']) is int and 0<model['n']<=8,'Register bound')
    require(type(model['electrons']) is int and 0<=model['electrons']<=model['n'],'Occupation bound')
    require(9<=len(model['blocks'])<=128 and 0<len(model['coefficients'])<=128,'Model storage bound')
    blocks=[decode(raw) for raw in model['blocks']];total=decode(model['coefficients'])
    require(int(model['denominator'])>0,'Positive coefficient denominator required')
    require(sum(map(len,blocks))==len(total) and {w:c for b in blocks for w,c in b.items()}==total,'Partition/target mismatch')
    for b in blocks:
        for w in b:
            require(len(dict(w))==len(w) and all(type(q) is int and 0<=q<model['n'] and p in 'XYZ' for q,p in w),'Invalid word')
        e.descriptor(b)
    return blocks


def pair_data(model):
    """Integer [Gi,Gj]|HF> amplitudes, with the parallel component removed."""
    blocks=validate_model(model);descs=[e.descriptor(b) for b in blocks]
    v=(1<<model['electrons'])-1;den=int(model['denominator']);pairs=[]
    D={w:c for b,d in zip(blocks,descs) if d[0]==0 for w,c in b.items()}
    dd=e.descriptor(D);destinations=set()
    for i,di in enumerate(descs):
        for j in range(i+1,len(descs)):
            dj=descs[j];w=v^di[0]^dj[0]
            b=e.value(dj,v)*e.value(di,v^dj[0])-e.value(di,v)*e.value(dj,v^di[0])
            if w!=v and b:pairs.append([i,j,w,str(b)]);destinations.add(w)
    ratios=[];energies=[]
    for w in sorted(destinations):
        gap=e.value(dd,w)-e.value(dd,v)
        rho=float(np.sinc(gap/(2*np.pi*den))**2)
        require(math.isfinite(rho) and rho>=0,'Invalid supplied weight')
        numerator,divisor=rho.as_integer_ratio()
        ratios.append([w,str(numerator),str(divisor)]);energies.append([w,str(gap)])
    scale=max((int(r[2]) for r in ratios),default=1)
    weights=[[w,str(int(n)*(scale//int(d)))] for w,n,d in ratios]
    return dict(pairs=pairs,weight_ratios=ratios,integer_weights=weights,
        weight_scale=str(scale),objective_denominator=str(4*den**4*scale),energy_gap_numerators=energies)


def amplitudes(order,pairs):
    pos={unit:at for at,unit in enumerate(order)};out={}
    for i,j,w,b in pairs:
        if i in pos and j in pos:out[w]=out.get(w,0)+(int(b) if pos[i]<pos[j] else -int(b))
    return out


def objective(order,data):
    a=amplitudes(order,data['pairs'])
    return sum(int(weight)*a.get(w,0)**2 for w,weight in data['integer_weights'])


def freeze_domain(model,warm_order,width=9):
    cpu=time.process_time();wall=time.monotonic();blocks=validate_model(model)
    require(width==9 and all(type(i) is int for i in warm_order) and sorted(warm_order)==list(range(len(blocks))),'Declared nine-unit permutation required')
    data=pair_data(model);counts=[]
    for offset in range(len(warm_order)-width+1):
        selected=set(warm_order[offset:offset+width])
        counts.append(sum(i in selected and j in selected for i,j,w,b in data['pairs']))
    offset=max(range(len(counts)),key=lambda k:(counts[k],-k))
    data.update(warm_order=list(warm_order),width=width,offset=offset,
        window=warm_order[offset:offset+width],outside_prefix=warm_order[:offset],outside_suffix=warm_order[offset+width:],
        selection='Maximum nonzero internal response pair count; earliest offset on ties; no physical error used',
        all_window_pair_counts=counts,domain_size=math.factorial(width),
        setup_CPU_seconds=time.process_time()-cpu,setup_wall_seconds=time.monotonic()-wall)
    return data


def solve(domain,cpu_limit=120):
    """Complete itertools enumeration, exact integer objective and tie breaking.

    Incomplete work returns lower bound zero and its actual incumbent, never a
    certificate. study.py refuses to proceed to physical errors in that case.
    """
    cpu=time.process_time();wall=time.monotonic();window=domain['window'];m=domain['width']
    require(m==9 and len(window)==m and len(set(window))==m,'Invalid frozen window')
    order=domain['warm_order'];require(order==domain['outside_prefix']+window+domain['outside_suffix'],'Outside/window inconsistency')
    require(all(type(i) is int for i in order) and sorted(order)==list(range(len(order))),'Invalid full permutation')
    selected=set(window);pairs=domain['pairs'];inside=[p for p in pairs if p[0] in selected and p[1] in selected]
    full=amplitudes(order,pairs);initial=amplitudes(window,inside);base={w:b-initial.get(w,0) for w,b in full.items()}
    weights={w:int(b) for w,b in domain['integer_weights']};affected=sorted({w for i,j,w,b in inside})
    constant=sum(weights[w]*b*b for w,b in base.items() if w not in affected)
    internals=[(i,j,w,int(b)) for i,j,w,b in inside]
    first=objective(order,domain);best=first;winner=list(window);ties=0;count=0
    termination='complete enumeration'
    for permutation in itertools.permutations(sorted(window)):
        if count%256==0 and time.process_time()-cpu>=cpu_limit:
            termination='CPU limit reached; incomplete domain';break
        pos={unit:at for at,unit in enumerate(permutation)}
        a={w:base.get(w,0) for w in affected}
        for i,j,w,b in internals:a[w]+=b if pos[i]<pos[j] else -b
        score=constant+sum(weights[w]*b*b for w,b in a.items());count+=1
        if score<best:best=score;winner=list(permutation);ties=1
        elif score==best:
            ties+=1
            if tuple(permutation)<tuple(winner):winner=list(permutation)
    optimal=domain['outside_prefix']+winner+domain['outside_suffix']
    require(objective(optimal,domain)==best,'Full reconstructed incumbent differs')
    complete=count==math.factorial(m);lower=best if complete else 0
    divisor=int(domain['objective_denominator'])
    return dict(status='certified restricted optimum' if complete else 'incomplete',termination=termination,
        domain_size=math.factorial(m),evaluated=count,complete=complete,minimizers=ties if complete else None,
        warm_start_integer=str(first),incumbent_integer=str(best),lower_bound_integer=str(lower),
        absolute_gap_integer=str(best-lower),relative_gap=0. if complete or best==0 else 1.,
        objective_denominator=str(divisor),initial_JD=first/divisor,incumbent_JD=best/divisor,
        absolute_optimality_tolerance='0',relative_optimality_tolerance='0',integrality_tolerance='0',
        selected_window=winner,order=optimal,CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        CPU_limit_seconds=cpu_limit,solver='Python itertools complete enumeration with exact integers',
        scope='Fixed contiguous nine-unit window, outside order fixed, exact supplied binary64 weight ratios only')
