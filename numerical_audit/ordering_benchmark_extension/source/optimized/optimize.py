"""Exact integer commutator amplitudes with stable floating objective selection.

The mathematical search matches frozen_candidate.optimize. Incremental floats only
screen candidates; every possible improvement is rescored from its integer
amplitudes and a sum of nonnegative terms before acceptance. Saved checkpoints
include each accepted sweep, counts, and search configuration.
"""
import json, math, time
from pathlib import Path
import numpy as np

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def atomic(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');t.replace(p)
def masks(w):
    x=z=y=0
    for q,p in w:
        if p in 'XY':x|=1<<q
        if p in 'ZY':z|=1<<q
        y+=p=='Y'
    return x,z,y
def fval(b,v):
    ans=0
    for w,c in b.items():
        _,z,y=masks(w);require(y%2==0,'Non-real term')
        ans+=c*(-1 if ((z&v).bit_count()+y//2)%2 else 1)
    return ans

def pair_setup(bs,xs,hf,den,pauli=False):
    """Preserves the frozen dictionary ordering and diagonal-gap arithmetic."""
    f0=[fval(b,hf) for b in bs];pairs={};m=len(bs)
    if pauli:
        zs=[masks(next(iter(b)))[1] for b in bs]
        for i in range(m):
            for j in range(i+1,m):
                a=(zs[i]&xs[j]).bit_count()%2;b=(zs[j]&xs[i]).bit_count()%2
                if a!=b:pairs[i,j]=(xs[i]^xs[j],2*f0[i]*f0[j]*(1 if a==0 else -1))
    else:
        desc=[[(masks(w)[1],masks(w)[2]//2,c) for w,c in b.items()] for b in bs]
        fs=[[sum(c*(-1 if ((z&(hf^x)).bit_count()+y)%2 else 1) for z,y,c in desc_i) for x in xs] for desc_i in desc]
        for i in range(m):
            for j in range(i+1,m):
                v=f0[j]*fs[i][j]-f0[i]*fs[j][i]
                if v:pairs[i,j]=(xs[i]^xs[j],v)
    diag={}
    for b,x in zip(bs,xs):
        if x==0:
            for w,v in b.items():diag[w]=diag.get(w,0)+v
    e0=fval(diag,hf)/den
    rho={d:float(np.sinc(((fval(diag,hf^d)/den-e0)/2)/np.pi)**2) for d in dict.fromkeys(d for d,_ in pairs.values())}
    require(all(d!=0 for d,v in pairs.values()),'Expected zero HF mean')
    return pairs,rho

def amplitudes(order,pairs):
    pos={i:p for p,i in enumerate(order)};a={}
    for (i,j),(d,b) in pairs.items():a[d]=a.get(d,0)+(b if pos[i]<pos[j] else -b)
    return a

def cost(a,rho,den):return .25*sum(rho[d]*(v/(den*den))**2 for d,v in a.items())

def search(start,pairs,rho,den,max_evals=None,checkpoint=None):
    order=list(start);history=[];evals=0;rescored=0;start_pass=0;elapsed_old=0.;stopped=False
    if checkpoint and Path(checkpoint).exists():
        saved=json.loads(Path(checkpoint).read_text());require(saved['start']==list(start) and saved['max_evals']==max_evals,'Stale checkpoint')
        if saved['complete']:return saved
        order=saved['order'];history=saved['history'];evals=saved['evaluations'];rescored=saved['rescored'];start_pass=saved['next_pass'];elapsed_old=saved['seconds']
    began=time.monotonic();a=amplitudes(order,pairs);value=cost(a,rho,den);den2=den*den;m=len(order)
    def save(complete,next_pass):
        out={'start':list(start),'order':order,'J':value,'history':history,'evaluations':evals,'rescored':rescored,'next_pass':next_pass,'max_evals':max_evals,'complete':complete,'seconds':elapsed_old+time.monotonic()-began}
        if checkpoint:atomic(checkpoint,out)
        return out
    for sweep in range(start_pass,10):
        best=value;best_move=None;best_amp=None
        for pos in range(m):
            if stopped:break
            moving=order[pos]
            for direction in (-1,1):
                trial=dict(a);tv=value;absolute_updates=0.;steps=0
                at=pos
                while 0<=at+direction<m:
                    if max_evals is not None and evals>=max_evals:stopped=True;break
                    crossed=order[at+direction];i,j=(crossed,moving) if direction<0 else (moving,crossed)
                    item=pairs.get((min(i,j),max(i,j)));evals+=1;steps+=1;at+=direction
                    if not item:continue
                    d,b=item;b=b if i<j else -b;old=trial.get(d,0);new=old-2*b
                    change=.25*rho[d]*((new/den2)**2-(old/den2)**2)
                    tv+=change;absolute_updates+=abs(change);trial[d]=new
                    # A standard sum-roundoff bound, deliberately conservative.
                    guard=32*np.finfo(float).eps*(steps+2)*(abs(value)+absolute_updates+abs(tv))
                    tol=max(1e-30,1e-12*best)
                    if tv<best-tol+guard:
                        exact=cost(trial,rho,den);rescored+=1
                        if exact<best-max(1e-30,1e-12*best):
                            best=exact;best_move=(pos,at);best_amp=dict(trial)
                if stopped:break
        if best_move is None:
            return save(True,sweep+1)
        old=value;pos,dst=best_move;item=order.pop(pos);order.insert(dst,item)
        a=amplitudes(order,pairs);value=cost(a,rho,den)
        require(a==best_amp and value==best,'Exact selected insertion mismatch')
        require(value<old-max(1e-30,1e-12*old),'Accepted non-improvement')
        history.append({'pass':sweep,'move':list(best_move),'old_J':old,'new_J':value,'evaluations_cumulative':evals})
        out=save(stopped or sweep==9,sweep+1)
        print('ACCEPT',Path(checkpoint).name if checkpoint else '',sweep,best_move,value,'evaluations',evals,'seconds',out['seconds'],flush=True)
        if stopped:return out
    return save(True,10)

def choose(candidates):
    winner=candidates[0]
    for c in candidates[1:]:
        tol=max(1e-30,1e-12*min(winner['J'],c['J']))
        if c['J']<winner['J']-tol or (abs(c['J']-winner['J'])<=tol and tuple(c['order'])<tuple(winner['order'])):winner=c
    return winner
