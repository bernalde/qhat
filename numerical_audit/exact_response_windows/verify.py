"""Independent exact-integer Pauli actions, DFS certificate and unfused replay."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy import sparse


def check(ok,message):
    if not ok:raise ValueError(message)


def decode(raw):return {tuple((q,p) for q,p in w):int(c) for w,c in raw}


def validate(model):
    n=model['n'];check(type(n) is int and 0<n<=8,'Independent register bound')
    check(type(model['electrons']) is int and 0<=model['electrons']<=n,'Independent occupation bound')
    check(int(model['denominator'])>0 and 9<=len(model['blocks'])<=128 and 0<len(model['coefficients'])<=128,'Independent storage bound')
    for raw in [model['coefficients']]+model['blocks']:
        check(bool(raw),'Empty supplied Pauli collection')
        for word,c in raw:
            check(len(word)<=n and len({q for q,p in word})==len(word),'Invalid Pauli support')
            check(all(type(q) is int and 0<=q<n and p in 'XYZ' for q,p in word),'Word outside independent register')
            int(c)
    return n


def integer_action(block,occupation):
    """Multiply individual Pauli factors using an integer i-exponent."""
    destination=None;amplitude=0
    for word,coefficient in block.items():
        target=occupation;exponent=0
        for q,pauli in word:
            bit=(occupation>>q)&1
            if pauli=='Z':exponent+=2*bit
            elif pauli=='X':target^=1<<q
            elif pauli=='Y':target^=1<<q;exponent+=1+2*bit
            else:raise ValueError('Unexpected Pauli factor')
        check(exponent%2==0,'Declared real units required')
        if destination is None:destination=target
        check(destination==target,'Unequal flip destinations in unit')
        amplitude+=coefficient*(-1 if exponent%4 else 1)
    check(destination is not None,'Empty unit')
    return destination,amplitude


def independent_pairs(model):
    validate(model)
    blocks=[decode(b) for b in model['blocks']];v=(1<<model['electrons'])-1
    pairs=[]
    for i,Gi in enumerate(blocks):
        for j in range(i+1,len(blocks)):
            Gj=blocks[j];vj,cj=integer_action(Gj,v);vi,ci=integer_action(Gi,v)
            wij,cij=integer_action(Gi,vj);wji,cji=integer_action(Gj,vi)
            check(wij==wji,'Commutator destinations differ')
            value=cj*cij-ci*cji
            if wij!=v and value:pairs.append([i,j,wij,str(value)])
    diagonal={w:c for b in blocks for w,c in b.items() if all(p=='Z' for q,p in w)}
    # Every diagonal Pauli appears in a singleton diagonal unit in this panel.
    for b in blocks:
        kinds=[all(p=='Z' for q,p in w) for w in b]
        check(all(kinds) or not any(kinds),'Mixed diagonal/off-diagonal unit')
    return pairs,diagonal


def full_amplitudes(order,pairs):
    check(all(type(i) is int for i in order) and sorted(order)==list(range(len(order))),'Nonintegral or invalid order')
    positions=[0]*len(order)
    for at,i in enumerate(order):positions[i]=at
    result={}
    for i,j,w,b in pairs:result[w]=result.get(w,0)+(int(b) if positions[i]<positions[j] else -int(b))
    return result


def independent_minimum(order,offset,width,pairs,weights,cpu_limit=120):
    """DFS adds a unit on the right; pairs acquire their sign once, on entry.

    Lehmer ranks cover the complete domain and a bitset rejects duplicate leaves.
    This shares neither permutations nor producer amplitude decomposition.
    """
    cpu=time.process_time();wall=time.monotonic();window=order[offset:offset+width];chosen=set(window)
    check(1<=width<=9 and len(chosen)==width,'Independent domain bound')
    destinations=sorted(weights);ids={w:k for k,w in enumerate(destinations)}
    a=[0]*len(destinations);positions={unit:at for at,unit in enumerate(order)};couplings={}
    for i,j,w,b in pairs:
        value=int(b);k=ids[w]
        if i in chosen and j in chosen:
            couplings[i,j]=(k,value);couplings[j,i]=(k,-value)
        else:a[k]+=value if positions[i]<positions[j] else -value
    total=math.factorial(width);seen=bytearray((total+7)//8);leaves=0;best=None;winner=None;ties=0
    weight_vector=[weights[w] for w in destinations]
    def visit(prefix,remaining,amplitude,rank):
        nonlocal leaves,best,winner,ties
        if not remaining:
            if leaves%256==0:check(time.process_time()-cpu<cpu_limit,'Independent enumeration CPU allowance exhausted')
            byte,bit=divmod(rank,8);check(not seen[byte]&(1<<bit),'Duplicate feasible order')
            seen[byte]|=1<<bit;leaves+=1
            score=sum(w*b*b for w,b in zip(weight_vector,amplitude))
            if best is None or score<best:best=score;winner=prefix;ties=1
            elif score==best:
                ties+=1
                if prefix<winner:winner=prefix
            return
        stride=math.factorial(len(remaining)-1)
        for index,new in enumerate(remaining):
            next_a=list(amplitude)
            for earlier in prefix:
                response=couplings.get((earlier,new))
                if response is not None:k,b=response;next_a[k]+=b
            visit(prefix+(new,),remaining[:index]+remaining[index+1:],next_a,rank+index*stride)
    visit((),tuple(sorted(window)),a,0)
    check(leaves==total and all(seen[k//8]&(1<<(k%8)) for k in range(total)),'Incomplete rank coverage')
    return dict(incumbent_integer=str(best),selected_window=list(winner),minimizers=ties,evaluated=leaves,
        rank_bitset_sha256=hashlib.sha256(seen).hexdigest(),CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)


def verify_certificate(model,domain,certificate):
    pairs,diagonal=independent_pairs(model);check(pairs==domain['pairs'],'Integer commutator responses differ')
    weights={};v=(1<<model['electrons'])-1;den=int(model['denominator']);gaps=[]
    check(len(domain['weight_ratios'])==len({p[2] for p in pairs}),'Missing supplied destination weight')
    scale=max((int(d) for w,n,d in domain['weight_ratios']),default=1)
    for w,n,d in domain['weight_ratios']:
        numerator,divisor=int(n),int(d);rho=numerator/divisor
        check(rho.as_integer_ratio()==(numerator,divisor) and 0<=rho<=1,'Noncanonical binary64 weight')
        gap=integer_action(diagonal,w)[1]-integer_action(diagonal,v)[1];x=gap/(2*math.pi*den)
        other=(math.sin(x)/x)**2 if x else 1.
        check(abs(rho-other)<5e-15,'Supplied weight differs from independent sinc check')
        weights[w]=numerator*(scale//divisor);gaps.append([w,str(gap)])
    check(sorted(weights)==sorted({p[2] for p in pairs}) and len(weights)==len(domain['weight_ratios']),'Duplicate/missing weight')
    check(domain['energy_gap_numerators']==gaps and domain['weight_scale']==str(scale),'Energy/scale mismatch')
    check(domain['integer_weights']==[[w,str(b)] for w,b in sorted(weights.items())],'Integer weights differ')
    check(domain['objective_denominator']==str(4*den**4*scale),'Objective normalization differs')
    order=domain['warm_order'];full_amplitudes(order,pairs);width=domain['width'];offset=domain['offset']
    check(width==9 and domain['window']==order[offset:offset+width],'Window identity changed')
    check(order==domain['outside_prefix']+domain['window']+domain['outside_suffix'],'Outside order changed')
    counts=[]
    for k in range(len(order)-width+1):
        ids=set(order[k:k+width]);counts.append(sum(i in ids and j in ids for i,j,w,b in pairs))
    check(counts==domain['all_window_pair_counts'] and offset==max(range(len(counts)),key=lambda k:(counts[k],-k)),'Window selection differs')
    other=independent_minimum(order,offset,width,pairs,weights)
    check(certificate['complete'] is True and certificate['termination']=='complete enumeration','Unfinished certificate')
    for k in ('incumbent_integer','selected_window','minimizers','evaluated'):check(other[k]==certificate[k],'Independent optimum differs: '+k)
    check(certificate['evaluated']==domain['domain_size']==math.factorial(width),'Wrong certified domain')
    optimum=domain['outside_prefix']+other['selected_window']+domain['outside_suffix']
    check(certificate['order']==optimum,'Full optimum reconstruction differs');full_amplitudes(optimum,pairs)
    for key,full in (('warm_start_integer',order),('incumbent_integer',optimum)):
        amplitude=full_amplitudes(full,pairs);score=sum(weights[w]*b*b for w,b in amplitude.items())
        check(str(score)==certificate[key],'Independent full-order objective differs')
    check(certificate['lower_bound_integer']==certificate['incumbent_integer'] and certificate['absolute_gap_integer']=='0' and certificate['relative_gap']==0.,'Invalid bound/gap')
    check(all(certificate[k]=='0' for k in ('absolute_optimality_tolerance','relative_optimality_tolerance','integrality_tolerance')),'Nonzero tolerances')
    check(certificate['objective_denominator']==domain['objective_denominator'],'Certificate normalization changed')
    divisor=int(domain['objective_denominator'])
    check(certificate['initial_JD']==int(certificate['warm_start_integer'])/divisor and certificate['incumbent_JD']==int(certificate['incumbent_integer'])/divisor,'Displayed surrogate differs')
    return other


def literal(word,n):
    check(0<n<=8,'Independent register bound');labels=dict(word);out=sparse.csr_matrix([[1.+0j]])
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    for q in reversed(range(n)):out=sparse.kron(out,matrices[labels.get(q,'I')],format='csr')
    return out


class Reference:
    def __init__(self,model):
        self.model=model;n=validate(model)
        self.dim=1<<n;self.v=(1<<model['electrons'])-1;den=int(model['denominator'])
        self.pool={w:literal(w,n) for w in decode(model['coefficients'])}
        self.H=sum((c/den)*self.pool[w] for w,c in decode(model['coefficients']).items())
        self.initial=np.zeros(self.dim,complex);self.initial[self.v]=1
        charges=[sparse.diags([sum((v>>q)&1 for q in range(p,n,2)) for v in range(self.dim)],dtype=float) for p in (0,1)]
        for raw in model['blocks']:
            block=decode(raw);G=sum((c/den)*self.pool[w] for w,c in block.items())
            for N in charges:
                comm=G@N-N@G
                check(not comm.nnz or max(abs(comm.data))<2e-12,'Unit spin-sector violation')
            words=list(block)
            for i,w in enumerate(words):
                for z in words[i+1:]:check((self.pool[w]@self.pool[z]-self.pool[z]@self.pool[w]).nnz==0,'Unit is not internally commuting')

    def exact_action(self):
        bound=float(max(np.asarray(abs(self.H).sum(axis=1)).ravel()));steps=max(1,math.ceil(bound/.25))
        state=self.initial.copy()
        for _ in range(steps):
            term=state.copy();total=state.copy()
            for k in range(1,25):term=(-1j/steps/k)*(self.H@term);total+=term
            state=total
        return state

    def score(self,record,domain):
        bs=([decode(b) for b in self.model['blocks']] if record['partition']=='conserving' else
            [{w:c} for w,c in sorted(decode(self.model['coefficients']).items())])
        den=int(self.model['denominator']);prefix=sparse.csr_matrix(self.H.shape,dtype=complex);column=np.zeros(self.dim,complex)
        for i in record['order']:
            G=sum((c/den)*self.pool[w] for w,c in bs[i].items())
            column+=(-.5j)*(prefix@(G@self.initial)-G@(prefix@self.initial));prefix+=G
        column[self.v]=0;gaps=self.H.diagonal().real-self.H[self.v,self.v].real
        column*=np.sinc(gaps/(2*np.pi));return float(np.vdot(column,column).real)

    def replay(self,rec,steps):
        schedule=[(tuple((q,p) for q,p in w),int(c)/int(self.model['denominator'])) for w,c in rec['schedule']]
        offsets=np.cumsum([0]+rec['unit_lengths']);a=self.initial.copy()
        def apply(unit,dt):
            nonlocal a
            for k in range(offsets[unit],offsets[unit+1]):
                word,c=schedule[k];theta=c*dt;a=np.cos(theta)*a-1j*np.sin(theta)*(self.pool[word]@a)
        m=len(rec['unit_lengths'])
        for _ in range(steps):
            if rec['formula']==1:
                for unit in range(m):apply(unit,1/steps)
            else:
                for unit in range(m-1):apply(unit,1/(2*steps))
                apply(m-1,1/steps)
                for unit in reversed(range(m-1)):apply(unit,1/(2*steps))
        return a


def physical_metrics(exact,actual,model):
    norm=float(np.sqrt(sum(abs(actual)**2)));a=actual/norm;z=exact/np.sqrt(sum(abs(exact)**2));inner=np.vdot(z,a)
    residual=a-z*inner;loss=float(sum(abs(residual)**2));n=model['n'];hf=(1<<model['electrons'])-1
    target=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)];number=joint=0.
    for v in range(len(a)):
        if v.bit_count()!=model['electrons']:number+=abs(a[v])**2
        if [sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=target:joint+=abs(a[v])**2
    return dict(error=loss/(1+abs(inner)),infidelity=loss,norm=norm,number_leakage=float(number),joint_sector_leakage=float(joint))


def verify(output):
    from study import bind_schedules,sha,write,old_case,BUDGETS,ARMS
    protocol=bind_schedules(output);certificates=[];physical=[];maxstate=maxmetric=maxscore=0.;count=0;artifacts={}
    for path in sorted((output/'designs').glob('*.json')):
        case=path.stem;model=json.loads((output/'inputs'/path.name).read_text());domain=json.loads((output/'domains'/path.name).read_text())
        certificate=json.loads((output/'certificates'/path.name).read_text());design=json.loads(path.read_text())
        check(certificate['domain_sha256']==sha(output/'domains'/path.name),'Certificate domain binding differs')
        check(design['domain_sha256']==certificate['domain_sha256'] and design['certificate_sha256']==sha(output/'certificates'/path.name),'Design certificate binding differs')
        original,point,prior=old_case(case)
        check(domain['warm_order']==point['methods']['optimized_JD']['order'],'Warm start changed')
        check(domain['control_orders']=={arm:point['methods'][arm]['order'] for arm in ('current','pauli_JD_large')},'Strong controls changed')
        for k,v in prior.items():check(domain['prior'][k]==v,'Prior provenance differs')
        if protocol['preflight']:
            original['coefficients']=[[w,str(7*int(c))] for w,c in original['coefficients']]
            original['blocks']=[[[w,str(7*int(c))] for w,c in b] for b in original['blocks']]
            original['denominator']=str(5*int(original['denominator']))
        check(model==original,'Target differs from authenticated input')
        independently=verify_certificate(model,domain,certificate);certificates.append(dict(case=case,**independently))
        reference=Reference(model);exact=reference.exact_action();check(len(design['records'])==6,'Missing controls')
        check([(r['arm'],r['formula']) for r in design['records']]==list(ARMS),'Comparator identities changed')
        for index,rec in enumerate(design['records']):
            expected=(certificate['order'] if rec['arm']=='certified_window' else domain['warm_order'] if rec['arm']=='heuristic' else domain['control_orders'][rec['arm']])
            check(rec['order']==expected,'Comparator order changed')
            partition='pauli' if rec['arm']=='pauli_JD_large' else 'conserving';check(rec['partition']==partition,'Partition changed')
            bs=([decode(b) for b in model['blocks']] if partition=='conserving' else [{w:c} for w,c in sorted(decode(model['coefficients']).items())])
            check(all(type(i) is int for i in rec['order']) and sorted(rec['order'])==list(range(len(bs))),'Invalid circuit order')
            expected_schedule=[[[[q,p] for q,p in w],str(c)] for i in rec['order'] for w,c in sorted(bs[i].items())]
            check(rec['schedule']==expected_schedule and rec['unit_lengths']==[len(bs[i]) for i in rec['order']],'Circuit target/unit boundaries changed')
            if rec['arm'] in ('heuristic','certified_window'):
                value=reference.score(rec,domain);expected_score=certificate['initial_JD'] if rec['arm']=='heuristic' else certificate['incumbent_JD']
                maxscore=max(maxscore,abs(value-expected_score));check(np.isclose(value,expected_score,rtol=3e-10,atol=5e-18),'Literal surrogate differs')
            for budget in BUDGETS:
                label=f'{case}.{index}.b{budget}';rp=output/'rows'/(label+'.json');row=json.loads(rp.read_text())
                lengths=rec['unit_lengths'];cap=budget*len(model['coefficients'])
                if rec['formula']==1:steps=cap//sum(lengths);rotations=steps*sum(lengths)
                else:
                    steps=(cap-lengths[0])//(2*sum(lengths)-lengths[0]-lengths[-1])
                    rotations=steps*(2*sum(lengths)-lengths[0]-lengths[-1])+lengths[0]
                expected_identity=(case,index,rec['arm'],rec['formula'],budget,steps,rotations,cap)
                check(tuple(row[k] for k in ('case','record_index','arm','formula','budget_units','steps','rotations','rotation_cap'))==expected_identity,'Physical row identity/accounting differs')
                check(row['design_sha256']==sha(path) and row['state_file']=='states/'+label+'.npz' and row['state_sha256']==sha(output/row['state_file']),'State/design binding differs')
                with np.load(output/row['state_file']) as states:actual=states['actual']
                other=reference.replay(rec,steps);delta=float(np.max(abs(actual-other)));maxstate=max(maxstate,delta)
                check(delta<2e-10,'Independent unfused propagation differs');values=physical_metrics(exact,other,model)
                for key,value in values.items():maxmetric=max(maxmetric,abs(row[key]-value));check(abs(row[key]-value)<2e-12,'Physical metric differs')
                physical.append(dict(case=case,record_index=index,budget_units=budget,producer_error=row['error'],independent_error=values['error']))
                artifacts[str(rp.relative_to(output))]=sha(rp);artifacts[row['state_file']]=sha(output/row['state_file']);count+=1
        print('VERIFIED',case,independently['evaluated'],flush=True)
    check(count==protocol['expected_circuits'],'Incomplete independent grid')
    classification=lambda d:'win' if d < -1e-12 else 'loss' if d>1e-12 else 'tie';ambiguities=[]
    for row in physical:
        if (row['producer_error']>1e-12)!=(row['independent_error']>1e-12):ambiguities.append(dict(kind='floor',row=row))
    for at,a in enumerate(physical):
        for b in physical[at+1:]:
            if (a['case'],a['budget_units'])!=(b['case'],b['budget_units']):continue
            if classification(a['producer_error']-b['producer_error'])!=classification(a['independent_error']-b['independent_error']):ambiguities.append(dict(kind='WLT',a=a,b=b))
    write(output/'INDEPENDENT_ERRORS.json',dict(rows=physical,decision_agreement=not ambiguities,ambiguities=ambiguities))
    write(output/'INDEPENDENT_CERTIFICATES.json',certificates)
    for name in ('INDEPENDENT_ERRORS.json','INDEPENDENT_CERTIFICATES.json'):artifacts[name]=sha(output/name)
    write(output/'VERIFIED.json',dict(status='passed' if not ambiguities else 'numerically ambiguous',rows=count,
        certificates=len(certificates),decision_agreement=not ambiguities,max_state_delta=maxstate,max_metric_delta=maxmetric,
        max_surrogate_delta=maxscore,schedule_freeze_sha256=sha(output/'SCHEDULES_FROZEN.json'),artifacts=artifacts))
    check(not ambiguities,'Independent reporting decisions differ')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);verify(p.parse_args().output.resolve())
