"""Literal Pauli/BCH response and unfused full-space replay, independent of search.

Uses study only for frozen IO/schema, never its response or pulse implementation.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def literal(word,n):
    if not 0<n<=8:raise ValueError('Independent literal size limit')
    local={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),
           'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    labels=dict(word);matrix=sparse.csr_matrix([[1.+0j]])
    for q in reversed(range(n)):matrix=sparse.kron(matrix,local[labels.get(q,'I')],format='csr')
    return matrix


def decoded(raw):return {tuple((int(q),p) for q,p in w):int(c) for w,c in raw}


def units(model,method):
    n=model['n'];den=int(model['denominator'])
    blocks=([{w:c} for w,c in sorted(decoded(model['coefficients']).items())]
            if method['partition']=='pauli' else [decoded(b) for b in model['blocks']])
    if not 0<n<=8 or not 0<len(blocks)<=128:raise ValueError('Independent unit size limit')
    out=[]
    for index in method['order']:
        G=sparse.csr_matrix((1<<n,1<<n),dtype=complex)
        for w,c in blocks[index].items():G+=(c/den)*literal(w,n)
        out.append(G)
    return out


def independent_score(model,method,psi,T):
    gs=units(model,method);dim=len(psi)
    prefix=sparse.csr_matrix((dim,dim),dtype=complex);C=np.zeros((dim,dim),complex)
    for G in gs:
        C+=(-.5j)*(prefix@G-G@prefix).toarray();prefix+=G
    d=prefix.diagonal().real
    gaps=d[:,None]-d[None,:]
    Q=C*(T*np.exp(.5j*T*gaps)*np.sinc(T*gaps/(2*np.pi)))
    z=Q@psi;z-=psi*np.vdot(psi,z)
    return float(np.vdot(z,z).real)


def replay(model,method,psi,T,steps):
    """Unfused symmetric formula; direct sparse Pauli action, no phase masks."""
    if steps<1 or model['n']>8:raise ValueError('Invalid independent circuit')
    flat=[(tuple((q,p) for q,p in w),int(c)/int(model['denominator'])) for w,c in method['schedule']]
    lengths=method['unit_lengths'];offsets=np.cumsum([0]+lengths)
    matrices=[literal(w,model['n']) for w,c in flat]
    actual=psi.copy()
    def apply(unit,dt):
        nonlocal actual
        for index in range(offsets[unit],offsets[unit+1]):
            theta=flat[index][1]*dt
            actual=np.cos(theta)*actual-1j*np.sin(theta)*(matrices[index]@actual)
    for _ in range(steps):
        if method['formula']==1:
            for u in range(len(lengths)):apply(u,T/steps)
        else:
            for u in range(len(lengths)-1):apply(u,T/(2*steps))
            apply(len(lengths)-1,T/steps)
            for u in reversed(range(len(lengths)-1)):apply(u,T/(2*steps))
    return actual


def independent_metrics(exact,actual,n,electrons):
    actual=actual/np.linalg.norm(actual);exact=exact/np.linalg.norm(exact)
    overlap=np.vdot(exact,actual);residual=actual-exact*overlap
    infidelity=float(sum(abs(residual)**2))
    number=sum(abs(actual[i])**2 for i in range(1<<n) if i.bit_count()!=electrons)
    hf=(1<<electrons)-1
    target=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)]
    spin=sum(abs(actual[i])**2 for i in range(1<<n)
             if [sum((i>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=target)
    return dict(error=infidelity/(1+abs(overlap)),infidelity=infidelity,
                number_leakage=float(number),spin_projection_leakage=float(spin))


def verify(output):
    from study import bind,sha,write,require,OLD,CASES,STATES,TIMES,BUDGETS,ARMS
    freeze=bind(output)
    require(json.loads((output/'RUN_COMPLETE.json').read_text())['rows']==432,'Run incomplete')
    artifacts={};count=0;max_state=max_score=max_metric=0.
    for case in CASES:
        model=json.loads((OLD/'inputs'/case/'model.json').read_text())
        design=json.loads((output/'designs'/(case+'.json')).read_text())
        n=model['n'];require(n<=8,'Size guard');dim=1<<n
        H=sparse.csr_matrix((dim,dim),dtype=complex)
        for w,c in decoded(model['coefficients']).items():H+=(c/int(model['denominator']))*literal(w,n)
        with np.load(output/design['state_file']) as data:vectors=data['vectors']
        require(vectors.shape==(dim,4),'Missing full-space states')
        hf=(1<<model['electrons'])-1
        counts=lambda v,p:sum((v>>q)&1 for q in range(p,n,2))
        ids=[i for i in range(dim) if all(counts(i,p)==counts(hf,p) for p in (0,1))]
        require(ids==design['state_preparation']['sector'],'Wrong state sector')
        eigvals=np.linalg.eigvalsh(H.toarray()[np.ix_(ids,ids)])
        g,h=vectors[:,:2].T
        require(abs(np.vdot(g,h))<1e-12 and max(abs(np.linalg.norm(vectors,axis=0)-1))<1e-12,'State orthogonality/norm')
        require(np.max(abs(vectors[:,2]-(g+h)/np.sqrt(2)))<1e-13 and
                np.max(abs(vectors[:,3]-(g+1j*h)/np.sqrt(2)))<1e-13,'Coherent phases changed')
        for at,label in enumerate(STATES):
            psi=vectors[:,at]
            require(sum(abs(psi[i])**2 for i in range(dim) if i not in ids)<1e-24,'State leaks outside preparation sector')
            if at<2:
                residual=float(np.linalg.norm(H@psi-eigvals[at]*psi))
                require(residual<2e-10,'State is not claimed sector eigenstate')
                require(abs(psi[np.argmax(abs(psi))].imag)<1e-12 and psi[np.argmax(abs(psi))].real>0,'State phase convention')
                require(abs(residual-design['state_preparation']['labels'][label]['eigen_residual'])<2e-10,'State residual metadata')
            meta=design['state_preparation']['labels'][label]
            ipr=float(1/sum(abs(psi)**4));require(abs(ipr-meta['participation_ratio'])<1e-10,'Correlation metadata')
            require(meta['determinant_like']==(ipr<=1+1e-8),'Correlation label mismatch')
        require([(p['state'],p['T']) for p in design['points']]==[(s,T) for s in STATES for T in TIMES],'Missing or selected state/time inputs')
        old=json.loads((OLD/'designs'/(case+'.json')).read_text())
        previous=next(p for p in old['points'] if p['state']=='HF')['methods']
        for at,point in enumerate(design['points']):
            psi=vectors[:,point['state_index']];T=point['T']
            exact=expm_multiply(-1j*T*H,psi)
            require(set(point['methods'])==set(ARMS),'Comparator missing')
            for name in ARMS:
                method=point['methods'][name]
                bs=([{w:c} for w,c in sorted(decoded(model['coefficients']).items())]
                    if method['partition']=='pauli' else [decoded(b) for b in model['blocks']])
                require(sorted(method['order'])==list(range(len(bs))),'Order not integral permutation')
                expanded=[(w,c) for i in method['order'] for w,c in sorted(bs[i].items())]
                supplied=[(tuple((q,p) for q,p in w),int(c)) for w,c in method['schedule']]
                require(expanded==supplied and dict(supplied)==decoded(model['coefficients']),'Schedule/target mismatch')
                require(method['unit_lengths']==[len(bs[i]) for i in method['order']],'Unit boundaries changed')
                if 'previous_method' in method:
                    require(method['schedule']==previous[method['previous_method']]['schedule'],'Transfer changed HF-trained order')
                if name=='state_symmetric':require(method['order']==point['methods']['state_JD']['order'],'Symmetric reuse order changed')
                score=independent_score(model,method,psi,T);delta=abs(score-method['score']);max_score=max(max_score,delta)
                require(np.isclose(score,method['score'],rtol=2e-9,atol=5e-18),'Complex-state response score mismatch')
                for budget in BUDGETS:
                    label=case+'.'+str(at)+'.'+name+'.b'+str(budget)
                    rp=output/'rows'/(label+'.json');row=json.loads(rp.read_text())
                    cap=budget*len(model['coefficients']);lengths=method['unit_lengths']
                    if method['formula']==1:steps=cap//sum(lengths);rotations=steps*sum(lengths)
                    elif len(lengths)==1:steps=1;rotations=sum(lengths)
                    else:
                        steps=(cap-lengths[0])//(2*sum(lengths)-lengths[0]-lengths[-1])
                        rotations=steps*(2*sum(lengths)-lengths[0]-lengths[-1])+lengths[0]
                    expected=(case,point['state'],T,at,name,budget,method['formula'],steps,rotations,cap)
                    got=tuple(row[k] for k in ('case','state','T','point','arm','budget_units','formula','steps','rotations','rotation_cap'))
                    require(got==expected and rotations<=cap and row['design_sha256']==sha(output/'designs'/(case+'.json')),'Circuit ledger mismatch')
                    require(sha(output/row['state_file'])==row['state_sha256'],'Propagated state changed')
                    with np.load(output/row['state_file']) as data:actual=data['actual']
                    reference=replay(model,method,psi,T,steps)
                    difference=float(np.max(abs(reference-actual)));max_state=max(max_state,difference)
                    require(difference<2e-10,'Independent unfused circuit mismatch')
                    metrics=independent_metrics(exact,reference,n,model['electrons'])
                    for key,value in metrics.items():
                        difference=abs(value-row[key]);max_metric=max(max_metric,difference)
                        require(difference<2e-12,'Independent physical metric mismatch')
                    require(abs(np.linalg.norm(actual)-row['norm'])<1e-12 and abs(row['norm']-1)<1e-8,'Norm ledger mismatch')
                    artifacts[str(rp.relative_to(output))]=sha(rp)
                    artifacts[row['state_file']]=sha(output/row['state_file']);count+=1
            print('VERIFIED',case,point['state'],T,flush=True)
    expected_rows={p for p in artifacts if p.startswith('rows/')}
    require(expected_rows=={str(p.relative_to(output)) for p in (output/'rows').glob('*.json')},'Extra or missing physical rows')
    require(count==432,'Incomplete independent replay')
    write(output/'VERIFIED.json',dict(status='passed',rows=count,max_state_delta=max_state,
        max_score_delta=max_score,max_metric_delta=max_metric,artifacts=artifacts,
        sources=freeze['sources'],schedule_freeze_sha256=sha(output/'SCHEDULES_FROZEN.json')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    verify(p.parse_args().output.resolve())
