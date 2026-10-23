"""Independent literal-Pauli/CAR, response, reference and circuit verification.

Does not import QHAT or the search/propagation engine. Uses saved frozen inputs.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
            'NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from campaign import sha,object_sha,write,require,decoded,frozen,CASES,BUDGETS


def literal(word,n,coefficient=1):
    require(0<n<=8,'Refuse large literal-Pauli allocation')
    dest=[];phase=[]
    for col in range(1<<n):
        row=col;factor=1+0j
        for q,p in word:
            bit=(col>>q)&1
            if p=='X':row^=1<<q
            elif p=='Y':row^=1<<q;factor*=(-1j if bit else 1j)
            elif p=='Z':factor*=(-1 if bit else 1)
            else:raise ValueError('Unknown Pauli label')
        dest.append(row);phase.append(factor)
    matrix=sparse.csr_matrix((coefficient*np.array(phase),(dest,range(1<<n))),shape=(1<<n,1<<n))
    return matrix,np.array(dest),np.array(phase)


def car(key,n):
    require(n<=8,'Refuse large CAR allocation')
    rows=[];cols=[];values=[]
    for col in range(1<<n):
        row=col;factor=1
        for q,creation in reversed(key):
            occupied=(row>>q)&1
            if occupied==creation:break
            factor*=(-1 if (row&((1<<q)-1)).bit_count()%2 else 1)
            row^=1<<q
        else:rows.append(row);cols.append(col);values.append(factor)
    return sparse.csr_matrix((values,(rows,cols)),shape=(1<<n,1<<n),dtype=complex)


def tensor_hamiltonian(one,two,constant):
    """Map saved InteractionOperator tensors directly through CAR, without QHAT."""
    n=len(one)
    require(0<n<=8 and one.shape==(n,n) and two.shape==(n,n,n,n),
            'Unsupported tensor dimensions before CAR allocation')
    require(np.all(np.isfinite(one)) and np.all(np.isfinite(two)) and np.isfinite(constant),
            'Nonfinite molecular tensors')
    H=constant*sparse.eye(1<<n,dtype=complex,format='csr')
    for p,q in np.argwhere(one!=0):
        H+=one[p,q]*car(((int(p),1),(int(q),0)),n)
    for p,q,r,s in np.argwhere(two!=0):
        H+=two[p,q,r,s]*car(((int(p),1),(int(q),1),(int(r),0),(int(s),0)),n)
    return H.toarray()


def tensor_bridge(one,two,constant,target,policy,hf,hf_energy):
    raw=tensor_hamiltonian(one,two,constant)
    require(constant==policy['input_constant'],'Input scalar convention mismatch')
    scalar=constant+policy['removed_parent_pauli_scalar']
    difference=float(np.max(abs(raw-target-scalar*np.eye(len(raw)))))
    allowance=policy['canonical_hermitian_adjustment_fermionic_l1']+policy['discarded_parent_fermionic_l1']
    require(policy['parent_cutoff']==1e-12 and 0<=allowance<=2e-8,
            'Unexpected target-policy adjustment on the bounded panel')
    require(difference<=allowance+2e-10,'Raw tensors do not match retained target and declared policy')
    hf_difference=abs(raw[hf,hf]-hf_energy)
    require(hf_difference<1e-8,'Independent CAR HF matrix element mismatch')
    return dict(raw_tensor_target_difference=difference,target_policy_l1_allowance=allowance,
                independent_CAR_HF_difference=float(hf_difference))


def tensor_product_oracle(word,n):
    """Conventional q=0 leftmost tensor factors, permuted to little-bit indexing."""
    require(0<n<=8,'Refuse large tensor-product oracle')
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),
              'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    labels=dict(word);matrix=np.ones((1,1),complex)
    for q in range(n):matrix=np.kron(matrix,matrices[labels.get(q,'I')])
    permutation=[int(f'{v:0{n}b}'[::-1],2) for v in range(1<<n)]
    return matrix[np.ix_(permutation,permutation)]


def check_current_ledger(method,position_budget,diagonal_count,off_count):
    initial=list(map(int,method['initial_slot_scores']));history=method['history']
    require(len(initial)==off_count+1 and len(history)==diagonal_count,'Incomplete exact search ledger')
    start=min(range(len(initial)),key=lambda i:(initial[i],i))
    require(method['start_slot']==start,'Initial slot is not the earliest exact minimum')
    previous=initial[start];slots=[]
    for index,item in enumerate(history):
        scores=list(map(int,item['slot_scores']));selected=min(range(len(scores)),key=lambda i:(scores[i],i))
        require(len(scores)==off_count+1 and item['index']==index and item['old_slot']==start,
                'Incorrect coordinate scan identity/count')
        require(int(item['old_score'])==previous==scores[start] and item['new_slot']==selected,
                'Exact coordinate history is inconsistent')
        require(int(item['new_score'])==scores[selected]<=previous,'Exact coordinate score increased')
        previous=scores[selected];slots.append(selected)
    require(previous==int(method['exact_score_numerator']) and slots==method['slots'],
            'Exact search terminal score/slots mismatch')
    require(position_budget==(diagonal_count+1)*(off_count+1)==method['evaluations'],
            'Exact position count does not match recorded scans')


def check_reference_edge(method,choices):
    expected=min(choices)[2:] if choices else (None,0)
    require((method['partner'],int(method['coupling_numerator']))==expected,
            'Reference partner was not selected by the frozen rule')


def independent_decision(rows,contexts):
    """Independent predicates and median implementation, without decisions.py."""
    arms=('current','pair_reference','optimized_JD','optimized_J0','pauli_JD',
          'pauli_JD_large','jw_signed','jw_magnitude','symmetric_JD')
    errors={(r['case'],r['state'],r['arm'],r['budget_units']):r['one_minus_overlap'] for r in rows}
    require(len(errors)==len(rows)==216,'Duplicate/incomplete independent decision grid')
    median=lambda values: (sorted(values)[len(values)//2] if len(values)%2 else
                           (sorted(values)[len(values)//2-1]+sorted(values)[len(values)//2])/2)
    contrasts={}
    for budget in (5,10,20,50):
        block={}
        for numerator in ('current','pair_reference'):
            for denominator in arms:
                if numerator==denominator:continue
                pairs=[(errors[c,s,numerator,budget],errors[c,s,denominator,budget]) for c,s in contexts]
                differences=[a-b for a,b in pairs]
                wins=sum(d < -1e-12 for d in differences);losses=sum(d > 1e-12 for d in differences)
                ratios=[a/b for a,b in pairs if b>1e-12]
                block[numerator+'_over_'+denominator]=dict(n=len(pairs),wins=wins,ties=len(pairs)-wins-losses,
                    losses=losses,ratios_retained=len(ratios),floor_exclusions=len(pairs)-len(ratios),
                    median_ratio=median(ratios) if ratios else None,max_ratio=max(ratios) if ratios else None)
        contrasts[str(budget)]=block
    costs=[];losses=[]
    for (case,state),methods in sorted(contexts.items()):
        current=methods['current']['setup_cpu_seconds']+methods['current']['search_cpu_seconds']
        corrected=methods['pair_reference']['setup_cpu_seconds']+methods['pair_reference']['search_cpu_seconds']
        require(current>0 and corrected>=0,'Invalid independent construction costs')
        costs.append(dict(case=case,state=state,current_cpu=current,corrected_cpu=corrected,ratio=corrected/current))
        a=errors[case,state,'pair_reference',10];b=errors[case,state,'current',10]
        if a-b>1e-12:losses.append(dict(case=case,state=state,ratio=a/(b if b>1e-12 else 1e-12)))
    primary=contrasts['10']['pair_reference_over_current'];median_cpu=median([x['ratio'] for x in costs])
    max_loss=max([x['ratio'] for x in losses],default=1.)
    conditions=dict(at_least_four_wins=primary['wins']>=4,no_loss_exceeds_twofold=max_loss<=2,
                    median_construction_CPU_at_most_twofold=median_cpu<=2)
    return dict(contrasts=contrasts,costs=costs,primary=primary,strict_loss_ratios=losses,
        median_construction_CPU_ratio=median_cpu,maximum_strict_loss_ratio=max_loss,conditions=conditions,
        correction_passes=all(conditions.values()),tie_and_denominator_floor=1e-12,
        loss_guard='strict losses use corrected/max(current,1e-12); descriptive ratios exclude denominator<=1e-12')


def bch_taylor(units,order):
    n=units[0].shape[0]
    first=np.zeros((n,n),complex);second=np.zeros_like(first)
    for i in order:
        G=units[i]
        second-=1j*(G@first)+.5*(G@G).toarray()
        first-=1j*G.toarray()
    H=1j*first
    return -1j*(second+.5*H@H)


def response(C,A,v):
    eigenvalues,U=eigh(A)
    gaps=eigenvalues[:,None]-eigenvalues[None,:]
    factors=np.exp(.5j*gaps)*np.sinc(gaps/(2*np.pi))
    Q=U@(factors*(U.conj().T@C@U))@U.conj().T
    column=Q[:,v].copy();column[v]=0
    return float(np.vdot(column,column).real)


def independent_pulses(lengths,steps,symmetric):
    offsets=np.cumsum([0]+lengths)
    if not symmetric:
        return [(i,1/steps) for _ in range(steps) for i in range(offsets[-1])]
    if len(lengths)==1:return [(i,1.) for i in range(offsets[-1])]
    # Unfused symmetric execution is mathematically identical, and its pulse
    # accounting is checked separately against fused boundary counts.
    out=[]
    for _ in range(steps):
        for unit in range(len(lengths)-1):out.extend((i,1/(2*steps)) for i in range(offsets[unit],offsets[unit+1]))
        out.extend((i,1/steps) for i in range(offsets[-2],offsets[-1]))
        for unit in reversed(range(len(lengths)-1)):out.extend((i,1/(2*steps)) for i in range(offsets[unit],offsets[unit+1]))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--qhat-root',type=Path,required=True);a=p.parse_args();f=frozen(a)
    require(not any(k=='qhat' or k.startswith('qhat.') for k in sys.modules),'Verifier imported QHAT')
    verified=[];raw_rows=[];contexts={};max_state=max_metric=max_reference=max_target=max_score=0.;references=0
    for case in CASES:
        d=json.loads((a.output/'designs'/(case+'.json')).read_text());n=d['n'];require(n<=8,'Size guard')
        model=json.loads((a.output/'inputs'/case/'model.json').read_text());den=int(d['denominator']);dim=1<<n
        require(sha(a.output/'inputs'/case/'tensors.npz')==model['tensors_sha256'],'Tensor hash changed')
        require(sha(a.output/'inputs'/case/'orbitals.npz')==model['orbitals_sha256'],'Orbital hash changed')
        total=decoded(d['coefficients']);P={w:literal(w,n) for w in total}
        oracle_words=[((0,'Y'),),((0,'Y'),(2,'Y'),(3,'Z'))] if n>=4 else [((0,'Y'),)]
        for word in oracle_words:
            require(np.array_equal(literal(word,n)[0].toarray(),tensor_product_oracle(word,n)),
                    'Tensor-product/occupation indexing bridge failed')
        H=sum((c/den)*P[w][0] for w,c in total.items()).toarray()
        fermionic=sparse.csr_matrix((dim,dim),dtype=complex)
        for parent in model['fermionic_parents']:
            for key,c in parent['components']:fermionic+=(int(c)/den)*car(key,n)
            fermionic-=(int(parent['removed_scalar'])/den)*sparse.eye(dim)
        target_delta=float(np.max(abs(fermionic.toarray()-H)));max_target=max(max_target,target_delta)
        require(target_delta<2e-11,'Retained CAR/Pauli target mismatch')
        with np.load(a.output/'inputs'/case/'tensors.npz') as z:
            occupied=d['electrons'];h=z['one_body'];g=z['two_body'];constant=float(z['constant'])
            contraction=constant+sum(h[q,q] for q in range(occupied))+sum(
                g[q,r,r,q]-g[q,r,q,r] for q in range(occupied) for r in range(occupied))
            require(abs(contraction-model['hf_energy'])<1e-8,'Independent tensor/HF contraction mismatch')
            bridge=tensor_bridge(h,g,constant,H,model['policy'],(1<<occupied)-1,model['hf_energy'])
        with np.load(a.output/'inputs'/case/'orbitals.npz') as z:
            C=z['mo_coeff'];S=z['ao_overlap'];occupancies=z['mo_occ']
            require(np.max(abs(C.T@S@C-np.eye(C.shape[1])))<1e-8,'Saved orbital orthogonality mismatch')
            require(np.sum(occupancies)==model['total_electrons'],'Saved total electron count mismatch')
            occupied_spatial=model['total_electrons']//2
            require(np.array_equal(occupancies,np.r_[np.full(occupied_spatial,2.),
                np.zeros(len(occupancies)-occupied_spatial)]),'Saved RHF occupations mismatch')
        blocks=[decoded(b) for b in d['blocks']]
        units=[sum((c/den)*P[w][0] for w,c in b.items()) for b in blocks]
        paulis=[(total[w]/den)*P[w][0] for w in sorted(total)]
        occupations=np.array([v.bit_count() for v in range(dim)])
        counts=[np.array([sum((v>>q)&1 for q in range(parity,n,2)) for v in range(dim)]) for parity in (0,1)]
        for block,G in zip(blocks,units):
            require(np.max(abs(G.toarray()-G.toarray().conj().T))<1e-12,'Non-Hermitian unit')
            for col in range(dim):
                amplitudes={}
                for w,c in block.items():
                    _,dest,phase=P[w]
                    require(phase[col].imag==0,'Expected real molecular unit')
                    row=int(dest[col]);amplitudes[row]=amplitudes.get(row,0)+c*int(phase[col].real)
                require(all(occupations[row]==occupations[col] for row,c in amplitudes.items() if c),
                        'Exact completed unit changes number')
                require(all(all(k[row]==k[col] for k in counts) for row,c in amplitudes.items() if c),
                        'Exact completed unit changes spin projection')
        D=np.diag(H.diagonal());case_rows=[];score_rows=[]
        flips=[sum(1<<q for q,p in next(iter(b)) if p in 'XY') for b in blocks]
        diagonal_ids=sorted((i for i,x in enumerate(flips) if x==0),
            key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
        off_ids=sorted((i for i,x in enumerate(flips) if x),
            key=lambda i:(-sum(abs(c) for c in blocks[i].values()),flips[i]))
        hf=(1<<d['electrons'])-1
        near=min((u for u in range(dim) if u!=hf and all(k[u]==k[hf] for k in counts)),
                 key=lambda u:((u^hf).bit_count(),u))
        require(len(d['points'])==2 and {x['state'] for x in d['points']}=={'HF','near'},'Wrong state selection grid')
        for point in d['points']:
            v=point['occupation'];prefix=case+'.'+point['state'];refpath=a.output/'references'/(prefix+'.npz')
            require(v==({'HF':hf,'near':near}[point['state']]),'Wrong declared HF/near determinant')
            require(point['position_budget']==(len(diagonal_ids)+1)*(len(off_ids)+1),'Wrong position budget')
            require(set(point['methods'])=={'current','pair_reference','optimized_JD','optimized_J0','pauli_JD',
                'pauli_JD_large','jw_signed','jw_magnitude','symmetric_JD'},'Wrong arm set')
            contexts[case,point['state']]=point['methods']
            with np.load(refpath) as z:initial=z['initial'];exact=z['exact']
            require(np.count_nonzero(initial)==1 and initial[v]==1,'Wrong input determinant')
            check=expm_multiply(-1j*sparse.csr_matrix(H),initial,traceA=-1j*np.trace(H))
            reference_delta=float(np.linalg.norm(check-exact));max_reference=max(max_reference,reference_delta);references+=1
            require(reference_delta<2e-10,'Exact reference mismatch')
            number=occupations==d['electrons'];spin=(counts[0]==counts[0][v])&(counts[1]==counts[1][v])
            integer_column={}
            for w,c in total.items():
                _,dest,phase=P[w];row=int(dest[v])
                integer_column[row]=integer_column.get(row,0)+c*int(phase[v].real)
            choices=[(-abs(c),row^v,row,c) for row,c in integer_column.items() if row!=v and c]
            for name,method in point['methods'].items():
                bs=paulis if method['partition']=='pauli' else units
                require(sorted(method['order'])==list(range(len(bs))),'Invalid permutation')
                schedule=[(tuple((q,p) for q,p in w),int(c)) for w,c in method['schedule']]
                require(dict(schedule)==total and len(schedule)==len(total),'Schedule target mismatch')
                require(object_sha(method['schedule'])==method['schedule_sha256'],'Schedule hash mismatch')
                expected=[(w,c) for i in method['order'] for w,c in sorted(
                    ({list(sorted(total))[i]:total[list(sorted(total))[i]]} if method['partition']=='pauli' else blocks[i]).items())]
                require(schedule==expected,'Unit order/schedule mismatch')
                require(method['unit_lengths']==[1 if method['partition']=='pauli' else len(blocks[i])
                    for i in method['order']],'Incorrect completed-unit boundaries')
                require(method['formula']==(2 if name=='symmetric_JD' else 1),'Wrong product formula')
                if name in ('current','pair_reference'):
                    require([i for i in method['order'] if flips[i]]==off_ids,'Off-unit order changed')
                    if name=='pair_reference':
                        require([x['unit'] for x in method['history']]==diagonal_ids,'Corrected diagonal priority changed')
                if name=='symmetric_JD':
                    base=point['methods']['optimized_JD']
                    require(method['order']==base['order'] and method['schedule_sha256']==base['schedule_sha256'],
                            'Symmetric control does not reuse declared order')
                if name in ('jw_signed','jw_magnitude'):
                    key=(lambda w:(total[w],w)) if name=='jw_signed' else (lambda w:(-abs(total[w]),w))
                    require([w for w,c in schedule]==sorted(total,key=key),'Wrong fixed Pauli sort')
                if name not in ('jw_signed','jw_magnitude','symmetric_JD'):
                    C=bch_taylor(bs,method['order']);A=D.copy()
                    if name=='optimized_J0':A[:]=0
                    if name=='pair_reference':
                        check_reference_edge(method,choices)
                        if method['partner'] is not None:
                            u=method['partner'];h=int(method['coupling_numerator'])/den
                            require(abs(H[u,v]-h)<1e-12,'Incorrect reference edge')
                            A[u,v]=A[v,u]=h
                    score=response(C,A,v);score_delta=abs(score-method['score']);max_score=max(max_score,score_delta)
                    require(score_delta<max(2e-11,2e-8*score),'Independent response score mismatch')
                    score_rows.append(dict(state=point['state'],arm=name,score=score,difference=score_delta))
                if name in ('optimized_JD','optimized_J0','pauli_JD'):
                    require(method['evaluations']<=point['position_budget'],'Matched search budget exceeded')
                if name in ('current','pair_reference'):
                    require(method['evaluations']==point['position_budget'],'Placement scan count mismatch')
                if name=='current':
                    diagonal_count=len(diagonal_ids)
                    check_current_ledger(method,point['position_budget'],diagonal_count,len(blocks)-diagonal_count)
                    expected_order=[]
                    for slot in range(len(off_ids)+1):
                        expected_order.extend(i for i,s in zip(diagonal_ids,method['slots']) if s==slot)
                        if slot<len(off_ids):expected_order.append(off_ids[slot])
                    require(method['order']==expected_order,'Current exact diagonal priority/slot assignment mismatch')
                    exact_score=int(method['exact_score_numerator'])/(4*int(method['exact_weight_denominator'])*den**4)
                    require(abs(exact_score-method['score'])<max(1e-14,1e-12*exact_score),
                            'Reported exact current score has the wrong normalization')
                if name=='pauli_JD_large':require(method['evaluations']<=20*len(paulis)*(len(paulis)-1),'Large search budget exceeded')
                coefficients=[c/den for w,c in schedule]
                for budget_units in BUDGETS:
                    label=prefix+'.'+name+'.b'+str(budget_units);path=a.output/'rows'/(label+'.json')
                    row=json.loads(path.read_text());statepath=a.output/row['state_file']
                    require(sha(statepath)==row['state_sha256'] and sha(refpath)==row['reference_sha256'],'State/reference hash changed')
                    require(row['design_sha256']==f['design_sha256'][case] and row['schedule_sha256']==method['schedule_sha256'],'Row identity mismatch')
                    lengths=method['unit_lengths'];Uwords=len(total);symmetric=method['formula']==2
                    steps=row['steps'];budget=budget_units*Uwords
                    expected_count=(Uwords if symmetric and len(lengths)==1 else
                        steps*(2*Uwords-lengths[0]-lengths[-1])+lengths[0] if symmetric else steps*Uwords)
                    require(row['rotation_budget']==budget and row['rotations']==expected_count<=budget,'Incorrect rotation budget')
                    if not symmetric:require(steps==budget//Uwords,'Wrong first-order steps')
                    elif len(lengths)>1:require(expected_count+2*Uwords-lengths[0]-lengths[-1]>budget,'Symmetric steps not maximal')
                    actual=initial.copy()
                    for term,fraction in independent_pulses(lengths,steps,symmetric):
                        word,c=schedule[term];_,dest,phase=P[word]
                        transformed=np.empty_like(actual);transformed[dest]=phase*actual
                        theta=coefficients[term]*fraction
                        actual=np.cos(theta)*actual-1j*np.sin(theta)*transformed
                    with np.load(statepath) as z:saved=z['actual']
                    delta=float(np.linalg.norm(actual-saved));max_state=max(max_state,delta)
                    require(delta<2e-10,'Independent circuit replay mismatch')
                    phi=exact/np.linalg.norm(exact);psi=saved/np.linalg.norm(saved);overlap=np.vdot(phi,psi)
                    residual=psi-phi*overlap;infidelity=float(np.vdot(residual,residual).real)
                    phase=overlap/abs(overlap) if overlap else 1.
                    aligned=psi-phase*phi
                    metrics=dict(one_minus_overlap=.5*float(np.vdot(aligned,aligned).real),infidelity=infidelity,
                        norm=float(np.linalg.norm(saved)),number_leakage=float(np.sum(abs(psi[~number])**2)),
                        spin_projection_leakage=float(np.sum(abs(psi[~spin])**2)))
                    metric_delta=max(abs(metrics[k]-row[k]) for k in metrics);max_metric=max(max_metric,metric_delta)
                    require(metric_delta<2e-12 and abs(metrics['norm']-1)<1e-8,'Metric/norm mismatch')
                    item=dict(row=label,row_sha256=sha(path),state_sha256=row['state_sha256'],
                              state_distance=delta,metric_difference=metric_delta)
                    case_rows.append(item);verified.append(item)
                    raw_rows.append(row)
            print('VERIFIED',case,point['state'],'rows',len(verified),flush=True)
        write(a.output/'verification'/(case+'.json'),dict(status='passed',rows=case_rows,scores=score_rows,
            target_difference=target_delta,tensor_bridge=bridge,design_sha256=f['design_sha256'][case]))
    require(len(verified)==216 and references==6,'Incomplete declared grid')
    require({x.name for x in (a.output/'rows').glob('*.json')}=={x['row']+'.json' for x in verified},
            'Unaccounted physical row files')
    write(a.output/'independent_decision.json',independent_decision(raw_rows,contexts))
    write(a.output/'verification_summary.json',dict(status='passed',new_replays=len(verified),references=references,
        max_state_distance=max_state,max_metric_difference=max_metric,max_reference_distance=max_reference,
        max_target_difference=max_target,max_score_difference=max_score))


if __name__=='__main__':main()
