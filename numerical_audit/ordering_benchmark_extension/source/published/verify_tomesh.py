"""Independent explicit-set RLF and exhaustive small-clique path checks."""
import itertools,json,random
from fractions import Fraction
from pathlib import Path
import algorithms as a
import extras

def rlf_reference(adj):
    remaining=set(range(len(adj)));groups=[]
    while remaining:
        eligible=set(remaining);excluded=set();group=[]
        while eligible:
            score=lambda i:(-(adj[i]&excluded).__len__(),(adj[i]&eligible).__len__(),i) if group else (-(adj[i]&eligible).__len__(),i)
            nxt=min(eligible,key=score);group.append(nxt);eligible.remove(nxt);new=adj[nxt]&eligible;excluded|=new;eligible-=new
        groups.append(group);remaining-=set(group)
    return groups

def main():
    rng=random.Random(923712);graph_checks=0
    for n in range(3,33):
        adj=[set() for i in range(n)]
        for i in range(n):
            for j in range(i):
                if rng.random()<.4:adj[i].add(j);adj[j].add(i)
        bit=[sum(1<<j for j in s) for s in adj]
        if rlf_reference(adj)!=extras.rlf(bit):raise RuntimeError('RLF mismatch')
        graph_checks+=1
    words=[tuple((q,p) for q,p in enumerate(s) if p!='I') for s in itertools.product('IXYZ',repeat=2)];path_checks=0
    for _ in range(12):
        ws=rng.sample(words,10);coeff=[rng.choice([-5,-3,-1,1,2,4]) for w in ws];ms=list(map(a.masks,ws));adj=a.adjacency(ms);order,meta=extras.tomesh(ws,coeff,ms,adj);groups=meta['groups'];k=len(groups)
        if k>6:continue
        C=[[sum(not a.anti(ms[i],ms[j]) for i in g for j in h) for h in groups] for g in groups]
        W=[[sum(abs(coeff[i]*coeff[j]) for i in g for j in h if a.anti(ms[i],ms[j])) for h in groups] for g in groups]
        final=[]
        def expand(path,c,w):
            remain=[i for i in range(k) if i not in path]
            if not remain:final.append((path,c,w));return
            tail=path[-1]
            for nxt in sorted(remain,key=lambda j:(-C[tail][j],j))[:2]:expand(path+(nxt,),c+C[tail][nxt],w+W[tail][nxt])
        for root in range(k):expand((root,),0,0)
        win=min(final,key=lambda p:(Fraction(p[2],max(p[1],1)),-p[1],p[0]))
        if tuple(meta['clique_order'])!=win[0] or meta['connection_score']!=win[1] or int(meta['adjacent_noncommuting_pair_weight_numerator'])!=win[2]:raise RuntimeError('Beam/full small path mismatch')
        path_checks+=1
    report={'RLF_explicit_set_cases':graph_checks,'beam32_equals_exhaustive_two_successor_small_cases':path_checks,'maximum_cliques_in_exhaustive_check':6,'adaptation_not_original_implementation':True}
    Path(__file__).with_name('tomesh_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
