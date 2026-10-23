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
    verified=[];max_state=max_metric=max_reference=max_target=max_score=0.;references=0
    for case in CASES:
        d=json.loads((a.output/'designs'/(case+'.json')).read_text());n=d['n'];require(n<=8,'Size guard')
        model=json.loads((a.output/'inputs'/case/'model.json').read_text());den=int(d['denominator']);dim=1<<n
        require(sha(a.output/'inputs'/case/'tensors.npz')==model['tensors_sha256'],'Tensor hash changed')
        require(sha(a.output/'inputs'/case/'orbitals.npz')==model['orbitals_sha256'],'Orbital hash changed')
        total=decoded(d['coefficients']);P={w:literal(w,n) for w in total}
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
        for point in d['points']:
            v=point['occupation'];prefix=case+'.'+point['state'];refpath=a.output/'references'/(prefix+'.npz')
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
            expected_edge=min(choices) if choices else None
            for name,method in point['methods'].items():
                bs=paulis if method['partition']=='pauli' else units
                require(sorted(method['order'])==list(range(len(bs))),'Invalid permutation')
                schedule=[(tuple((q,p) for q,p in w),int(c)) for w,c in method['schedule']]
                require(dict(schedule)==total and len(schedule)==len(total),'Schedule target mismatch')
                require(object_sha(method['schedule'])==method['schedule_sha256'],'Schedule hash mismatch')
                expected=[(w,c) for i in method['order'] for w,c in sorted(
                    ({list(sorted(total))[i]:total[list(sorted(total))[i]]} if method['partition']=='pauli' else blocks[i]).items())]
                require(schedule==expected,'Unit order/schedule mismatch')
                if name not in ('jw_signed','jw_magnitude','symmetric_JD'):
                    C=bch_taylor(bs,method['order']);A=D.copy()
                    if name=='optimized_J0':A[:]=0
                    if name=='pair_reference' and method['partner'] is not None:
                        u=method['partner'];h=int(method['coupling_numerator'])/den
                        require(expected_edge is not None and (u,int(method['coupling_numerator']))==expected_edge[2:],
                                'Reference partner was not selected by the frozen rule')
                        require(abs(H[u,v]-h)<1e-12,'Incorrect reference edge')
                        A[u,v]=A[v,u]=h
                    score=response(C,A,v);score_delta=abs(score-method['score']);max_score=max(max_score,score_delta)
                    require(score_delta<max(2e-11,2e-8*score),'Independent response score mismatch')
                    score_rows.append(dict(state=point['state'],arm=name,score=score,difference=score_delta))
                if name in ('optimized_JD','optimized_J0','pauli_JD'):
                    require(method['evaluations']<=point['position_budget'],'Matched search budget exceeded')
                if name in ('current','pair_reference'):
                    require(method['evaluations']==point['position_budget'],'Placement scan count mismatch')
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
            print('VERIFIED',case,point['state'],'rows',len(verified),flush=True)
        write(a.output/'verification'/(case+'.json'),dict(status='passed',rows=case_rows,scores=score_rows,
            target_difference=target_delta,design_sha256=f['design_sha256'][case]))
    require(len(verified)==216 and references==6,'Incomplete declared grid')
    write(a.output/'verification_summary.json',dict(status='passed',new_replays=len(verified),references=references,
        max_state_distance=max_state,max_metric_difference=max_metric,max_reference_distance=max_reference,
        max_target_difference=max_target,max_score_difference=max_score))


if __name__=='__main__':main()
