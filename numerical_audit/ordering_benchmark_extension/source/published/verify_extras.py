"""Independent small-system checks of extra baselines and Strang pulse counts."""
import itertools,json,sys
from pathlib import Path
import numpy as np
import scipy.linalg as la
import algorithms as a
import extras
from verify_algorithms import matrix

def main():
    rng=np.random.default_rng(93145);words=[tuple((q,p) for q,p in enumerate(s) if p!='I') for s in itertools.product('IXYZ',repeat=2)];ms=list(map(a.masks,words));mats=[matrix(w,2) for w in words];checks=0
    for _ in range(20):
        ids=list(map(int,rng.choice(16,12,replace=False)));co=[int(v) or 1 for v in rng.integers(-9,10,12)];adj=a.adjacency([ms[i] for i in ids]);c,r=a.appendix_orders(co,adj)
        for expected,reverse in [(c,False),(r,True)]:
            remaining=list(range(12));order=[]
            while remaining:
                scores={i:sum(np.allclose(mats[ids[i]]@mats[ids[j]],mats[ids[j]]@mats[ids[i]],atol=0,rtol=0) for j in (remaining if reverse else order)) for i in remaining}
                nxt=min(remaining,key=lambda i:((-scores[i] if reverse else scores[i]),-abs(co[i]),i));order.append(nxt);remaining.remove(nxt)
            if order!=expected:raise RuntimeError('AppendixA dense mismatch')
            checks+=1
    def W(s):return tuple((i,p) for i,p in enumerate(s) if p!='I')
    parent=[{W('ZII'):-2},{W('IZI'):-2},{W('ZII'):-1,W('IZI'):-1,W('ZZI'):1},{W('XXI'):2,W('YYI'):2},{W('XXI'):1,W('YYI'):1,W('XXZ'):-1,W('YYZ'):-1},{W('IXX'):2,W('IYY'):2},{W('IXX'):1,W('IYY'):1,W('ZXX'):-1,W('ZYY'):-1}]
    components=[[[0,1],[0,0]],[[1,1],[1,0]],[[0,1],[1,1],[1,0],[0,0]],[[0,1],[1,0]],[[0,1],[2,1],[2,0],[1,0]],[[1,1],[2,0]],[[0,1],[1,1],[2,0],[0,0]]]
    total={}
    for b in parent:
        for w,c in b.items():total[w]=total.get(w,0)+c
    d={'blocks':{'parent':[[[w,str(c)] for w,c in b.items()] for b in parent]},'parent_metadata':[{'components':[c]} for c in components],'coefficients':[[w,str(c)] for w,c in total.items() if c]}
    lit,cons,meta=a.hastings(d);res=[]
    for t in (.01,.1,.4,1.):
        U=np.eye(8,dtype=complex);V=U.copy()
        for w,c in lit:U=la.expm(-1j*t*c/4*matrix(w,3))@U
        for w,c in cons:V=la.expm(-1j*t*c/4*matrix(w,3))@V
        res.append(float(np.linalg.norm(U-V)))
    if max(res)>1e-12:raise RuntimeError('Hastings consolidation changes unitary')
    # The emitted merged pulse sequence equals the literal repeated palindrome.
    strang=[]
    for r in (1,2,5):
        terms=[.7*mats[5],-.2*mats[11],.4*mats[14]];U=np.eye(4,dtype=complex);V=U.copy()
        for step in range(r):
            for m in terms+terms[::-1]:U=la.expm(-1j*m/(2*r))@U
        pulses=[(0,1)]
        for step in range(r):pulses +=[(1,1),(2,2),(1,1),(0,1 if step==r-1 else 2)]
        for i,f in pulses:V=la.expm(-1j*terms[i]*f/(2*r))@V
        if len(pulses)!=2*r*(len(terms)-1)+1:raise RuntimeError('Strang rotation count')
        strang.append(float(np.linalg.norm(U-V)))
    if max(strang)>1e-12:raise RuntimeError('Merged Strang differs')
    report={'appendix_A_dense_cases':checks,'hastings_consolidation_dense_residuals':res,'strang_merged_unmerged_residuals':strang,'third_party_code_redistributed':False}
    Path(__file__).with_name('extra_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
