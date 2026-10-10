#!/usr/bin/env python3
"""Independent full-rescoring and unfused-product checks for optimized code."""
import sys,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
import numpy as np
import optimize as opt
import run
import frozen_candidate as frozen

def main():
    n=c.normalized_case('HCN_s1.00');den=n['denominator'];hf=(1<<n['design']['electrons'])-1
    pbs=[{w:n['total'][w]} for w in sorted(n['total'],key=lambda w:(c.masks(w)[0],w))];px=[c.masks(next(iter(b)))[0] for b in pbs]
    checks={};t=time.monotonic()
    pairs,rho=opt.pair_setup(pbs,px,hf,den,True);fp,fr=frozen.pair_setup(pbs,px,hf,den)
    c.require(pairs==fp and rho==fr,'Specialized pair setup differs from frozen class algebra')
    checks['specialized_pair_setup_exact_equality']={'units':len(pbs),'nonzero_pairs':len(pairs)}
    # Select 64 terms spread throughout the molecular input, so the test includes
    # diagonal and off-diagonal words rather than an all-diagonal prefix.
    inds=np.linspace(0,len(pbs)-1,64,dtype=int).tolist();bs=[pbs[i] for i in inds];xs=[px[i] for i in inds]
    pairs,rho=opt.pair_setup(bs,xs,hf,den,True);evaluated=0
    for start in (list(range(64)),sorted(range(64),key=lambda i:(-sum(abs(v) for v in bs[i].values()),xs[i]))):
        a=opt.amplitudes(start,pairs);best=opt.cost(a,rho,den);bestmove=None
        for pos in range(64):
            for direction in (-1,1):
                trial=dict(a);working=list(start);at=pos
                while 0<=at+direction<64:
                    left=min(at,at+direction);i,j=working[left:left+2];item=pairs.get((min(i,j),max(i,j)))
                    if item:
                        d,b=item;b=b if i<j else -b;trial[d]=trial.get(d,0)-2*b
                    working[left:left+2]=[j,i];at+=direction
                    exact=opt.amplitudes(working,pairs);c.require(trial==exact,'Incremental integer insertion action differs from complete pair sum')
                    score=opt.cost(exact,rho,den);evaluated+=1
                    if score<best-max(1e-30,1e-12*best):best=score;bestmove=[pos,at]
        actual=opt.search(start,pairs,rho,den)
        c.require((actual['history'][0]['move'] if actual['history'] else None)==bestmove,'Stable optimizer first move differs from full oracle')
        if bestmove is not None:c.require(actual['history'][0]['new_J']==best,'Stable optimizer objective differs from full oracle')
    checks['exact_insertion_oracle']={'candidate_positions':evaluated,'starts':2,'molecular_terms':64}
    d=json.loads((HERE/'designs/CH2O_s1.00.design.json').read_text());raw,states=c.load_case('CH2O_s1.00');compiled=c.compile_case(raw,states)
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    residuals={}
    for name,method in d['methods'].items():
        if method['formula_order']!=2:continue
        terms=[compiled['term'](w,v) for w,v in c.decode_schedule(method['schedule'])];r=3
        ii,dd=run.strang_pulses(method['unit_lengths'],r);actual=compiled['initial'].copy();work=np.empty_like(actual)
        for i,mult in zip(ii,dd):actual,work=ev.apply_one_exponential(actual,work,terms[int(i)],int(mult)/(2*r),2**30)
        expected=compiled['initial'].copy();work=np.empty_like(expected)
        for step in range(r):
            for term in terms:expected,work=ev.apply_one_exponential(expected,work,term,1/(2*r),2**30)
            for term in reversed(terms):expected,work=ev.apply_one_exponential(expected,work,term,1/(2*r),2**30)
        error=float(np.linalg.norm(actual-expected));c.require(error<1e-12,'Fused block Strang differs from unfused palindrome');residuals[name]=error
    checks['strang_unmerged_product_norm_difference']=residuals;checks['seconds']=time.monotonic()-t
    checks['source_sha256']={p.name:c.sha(p) for p in (HERE/'check_oracles.py',HERE/'optimize.py',HERE/'run.py',HERE/'frozen_candidate.py')}
    c.write_atomic(HERE/'oracle_checks.json',checks);print(json.dumps(checks,indent=2),flush=True)
if __name__=='__main__':main()
