"""Independent dense-matrix and NetworkX checks for recovered ordering rules."""
import itertools,json,tempfile,subprocess,hashlib
from pathlib import Path
import numpy as np
import scipy.linalg as la
import networkx as nx
import algorithms as a
P={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
def matrix(w,n):
    d=dict(w);m=np.array([[1]])
    for q in range(n):m=np.kron(m,P[d.get(q,'I')])
    return m
def comm(A,B):return A@B-B@A

def fullW(order, mats):
    out=np.zeros_like(mats[0],dtype=complex)
    for b in range(len(order)):
        for c in range(b):
            bc=comm(mats[order[b]],mats[order[c]])
            for k in range(b+1):out+=(1 if k==b else 2)*comm(mats[order[k]],bc)
    return out

def coefficients(M,base):
    values=[np.trace(b@M)/len(M) for b in base]
    if max(abs(v.imag) for v in values)>1e-7:raise RuntimeError('Nonreal dense coefficients')
    return [round(v.real) for v in values]
def main():
    rng=np.random.default_rng(491730);n=2
    words=[tuple((q,p) for q,p in enumerate(s) if p!='I') for s in itertools.product('IXYZ',repeat=n)]
    matrices=[matrix(w,n) for w in words];ms=list(map(a.masks,words));product_checks=0
    for i,j in itertools.product(range(len(words)),repeat=2):
        x,z=ms[i][0]^ms[j][0],ms[i][1]^ms[j][1];k=ms.index((x,z));expected=(1j**a.phase(ms[i],ms[j]))*matrices[k]
        if not np.allclose(matrices[i]@matrices[j],expected,atol=0,rtol=0):raise RuntimeError('Pauli phase mismatch')
        product_checks+=1
    graph_checks=0
    for size in range(3,17):
        ids=list(rng.choice(len(words),size,replace=False));am=[ms[i] for i in ids];adj=a.adjacency(am);G=nx.Graph();G.add_nodes_from(range(size));G.add_edges_from((i,j) for i in range(size) for j in range(i) if a.anti(am[i],am[j]))
        groups=[];remaining=set(G)
        while remaining:
            eligible=set(remaining);group=[]
            while eligible:
                induced=G.subgraph(eligible);i=min(eligible,key=lambda j:(induced.degree(j),j));group.append(i);eligible-=set(induced[i])|{i}
            groups.append(sorted(group));remaining-=set(group)
        if sorted(map(sorted,a.independent_groups(adj)))!=sorted(groups):raise RuntimeError('NetworkX independent_set mismatch')
        graph_checks+=1
    v2_checks=0;bch_checks=0;strang_errors=[]
    with tempfile.TemporaryDirectory() as tmp:
        for trial in range(20):
            ids=sorted(map(int,rng.choice(range(1,len(words)),12,replace=False)));ws=[words[i] for i in ids];coeff=[int(v) for v in rng.integers(-9,10,12)];coeff=[v or 1 for v in coeff];mm=[c*matrices[i] for i,c in zip(ids,coeff)];incoming=sorted(range(12),key=lambda i:(-abs(coeff[i]),i));ord=[];history=[]
            for nxt in incoming:
                old=fullW(ord,mm);scores=[]
                for pos in range(len(ord)+1):scores.append(sum(abs(v) for v in coefficients(fullW(ord[:pos]+[nxt]+ord[pos:],mm)-old,matrices)))
                where=min(range(len(scores)),key=lambda p:(scores[p],p));ord.insert(where,nxt);history.append([nxt,where,str(scores[where])])
            cp=Path(tmp)/f'{trial}.tsv';cp.touch();payload=str(len(ws))+'\n'+''.join(f'{a.masks(ws[i])[0]} {a.masks(ws[i])[1]} {coeff[i]} {i}\n' for i in incoming)
            result=subprocess.run([str(Path(__file__).with_name('v2_exact')),str(cp)],input=payload,text=True,capture_output=True,check=True)
            nativehist=[[int(t.split()[0]),int(t.split()[1]),t.split()[2]] for t in cp.read_text().splitlines()]
            if nativehist!=history or json.loads(result.stdout)!=ord:raise RuntimeError(f'Native V2 mismatch {trial}')
            # Check restart reproduces exactly without appending extra records.
            resumed=subprocess.run([str(Path(__file__).with_name('v2_exact')),str(cp)],input=payload,text=True,capture_output=True,check=True)
            if resumed.stdout!=result.stdout or len(cp.read_text().splitlines())!=12:raise RuntimeError('Resume mismatch')
            v2_checks+=1
            bo,bh=a.first_order_insertion(ws,coeff,list(map(a.masks,ws)));o=[]
            for step,nxt in enumerate(incoming):
                vals=[]
                for pos in range(len(o)+1):
                    trialorder=o[:pos]+[nxt]+o[pos:];M=sum((comm(mm[trialorder[i]],mm[trialorder[j]])/1j for i in range(len(trialorder)) for j in range(i+1,len(trialorder))),start=np.zeros_like(mm[0],dtype=complex));vals.append(sum(abs(v) for v in coefficients(M,matrices)))
                p=min(range(len(vals)),key=lambda k:(vals[k],k));o.insert(p,nxt)
                if bh[step]!=[nxt,p,str(vals[p])]:raise RuntimeError('BCH1 dense mismatch')
            if o!=bo:raise RuntimeError('BCH1 final mismatch')
            bch_checks+=1
            if trial<3:
                # Eq4's prefix indices are reversed into chronological Strang.
                oo=ord[::-1];dt=1e-4;U=np.eye(2**n,dtype=complex)
                for i in oo+oo[::-1]:U=la.expm(-1j*dt/2*mm[i])@U
                Heff=1j*la.logm(U)/dt;pred=sum(mm)+dt**2/24*fullW(ord,mm)
                err=np.linalg.norm(Heff-pred);strang_errors.append(float(err))
                if err>2e-8:raise RuntimeError('Eq4/physical Strang convention mismatch')
    report={'pauli_dense_products':product_checks,'networkx_independent_set_cases':graph_checks,'exact_dense_v2_greedy_cases':v2_checks,'exact_dense_bch1_greedy_cases':bch_checks,'checkpoint_restart_cases':v2_checks,'strang_effective_H_residuals_at_dt1e-4':strang_errors}
    Path(__file__).with_name('algorithm_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
