"""Native V2 and explicitly adapted Tomesh ordering; see PROTOCOL.md."""
import hashlib,json,math,os,subprocess,time
from fractions import Fraction
from pathlib import Path
import algorithms as a

def native_v2(design,where):
    where=Path(where);where.mkdir(parents=True,exist_ok=True);total=a.block(design['coefficients']);words=sorted(total,key=lambda w:a.dense(w,design['n_qubits']));co=[total[w] for w in words];incoming=sorted(range(len(words)),key=lambda i:(-abs(co[i]),i))
    payload=str(len(words))+'\n'+''.join(f'{a.masks(words[i])[0]} {a.masks(words[i])[1]} {co[i]} {i}\n' for i in incoming)
    ip=where/'input.txt';fingerprint={'input_sha256':hashlib.sha256(payload.encode()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).with_name('v2_exact.cpp').read_bytes()).hexdigest()}
    meta=where/'checkpoint_provenance.json'
    if meta.exists():
        if json.loads(meta.read_text())!=fingerprint:raise RuntimeError('Native resume input/source changed')
    else:
        if (where/'accepted.tsv').exists():raise RuntimeError('Unattributed checkpoint')
        ip.write_text(payload);meta.write_text(json.dumps(fingerprint,indent=2)+'\n');(where/'accepted.tsv').touch()
    start=time.monotonic()
    with (where/'progress.log').open('a') as log:
        run=subprocess.run([str(Path(__file__).with_name('v2_exact')),str(where/'accepted.tsv')],input=payload,text=True,stderr=log,capture_output=False,stdout=subprocess.PIPE,check=True)
    order=json.loads(run.stdout);history=[[int(s.split()[0]),int(s.split()[1]),s.split()[2]] for s in (where/'accepted.tsv').read_text().splitlines()]
    if sorted(order)!=list(range(len(words))) or len(history)!=len(words):raise RuntimeError('Native incomplete permutation')
    return {'schedule':[(words[i],co[i]) for i in order[::-1]],'metadata':{'source':'Tranter arXiv:1912.07555v1 Eq11 and V.3','objective':'coefficient l1 of incremental V2, not full V2 norm','arithmetic':'exact arbitrary-precision integer coefficients','incoming':'descending absolute coefficient; dense IXYZ ties','position_ties':'leftmost','eq4_order':order,'chronological_order':'reverse of Eq4 prefix order; independent matrix-log verification','history':history,'latest_native_seconds':time.monotonic()-start,**fingerprint}}

def rlf(adj):
    remaining=(1<<len(adj))-1;groups=[]
    while remaining:
        eligible=remaining;excluded=0;group=[]
        first=min(a.bits(eligible),key=lambda i:(-(adj[i]&eligible).bit_count(),i))
        while eligible:
            i=first if not group else min(a.bits(eligible),key=lambda j:(-(adj[j]&excluded).bit_count(),(adj[j]&eligible).bit_count(),j))
            group.append(i);eligible^=1<<i;new=eligible&adj[i];eligible&=~new;excluded|=new
        groups.append(group)
        for i in group:remaining^=1<<i
    return groups

def tomesh(words,coeff,ms,adj,width=32):
    import networkx as nx
    groups=rlf(adj);k=len(groups);cliques=[]
    strings=[a.dense(w,max((q for v in words for q,p in v),default=-1)+1) for w in words]
    def distance(i,j):return sum((1+(x!='I' and y!='I')) for x,y in zip(strings[i],strings[j]) if x!=y)
    for group in groups:
        if len(group)<3:cliques.append(sorted(group));continue
        G=nx.Graph();G.add_nodes_from(sorted(group));G.add_weighted_edges_from((i,j,distance(i,j)) for p,i in enumerate(sorted(group)) for j in sorted(group)[:p])
        cycle=nx.algorithms.approximation.christofides(G);cut=max(range(len(cycle)-1),key=lambda p:(distance(cycle[p],cycle[p+1]),-p));route=cycle[cut+1:-1]+cycle[:cut+1]
        if len(route)!=len(group) or set(route)!=set(group):raise RuntimeError('TSP route invalid')
        cliques.append(min(route,route[::-1]))
    connection=[[0]*k for _ in range(k)];weights=[[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i):
            count=weight=0
            for u in groups[i]:
                for v in groups[j]:
                    if anti:=bool(adj[u]&(1<<v)):weight+=abs(coeff[u]*coeff[v])
                    else:count+=1
            connection[i][j]=connection[j][i]=count;weights[i][j]=weights[j][i]=weight
    finalists=[]
    for root in range(k):
        paths=[((root,),0,0,1<<root)]
        for depth in range(1,k):
            expanded=[]
            for path,c,w,used in paths:
                tail=path[-1];choices=sorted((j for j in range(k) if not used&(1<<j)),key=lambda j:(-connection[tail][j],j))[:2]
                for j in choices:expanded.append((path+(j,),c+connection[tail][j],w+weights[tail][j],used|(1<<j)))
            paths=sorted(expanded,key=lambda p:(-p[1],p[0]))[:width]
        finalists.extend(paths)
    winner=min(finalists,key=lambda p:(Fraction(p[2],max(p[1],1)),-p[1],p[0])) if finalists else ((),0,0,0)
    order=[i for g in winner[0] for i in cliques[g]]
    return order,{'source':'Tomesh arXiv:2203.12713 section IV; dqs-term-grouping commit5797aab6cc9ee9f60bdfa0003102fc16d5a6d2aa','method':'Tomesh-inspired deterministic beam32 adaptation, not a reproduction of original implementation','groups':groups,'ordered_cliques':cliques,'clique_order':winner[0],'connection_score':winner[1],'adjacent_noncommuting_pair_weight_numerator':str(winner[2]),'beam_width_per_root':width,'changes':['Dense IXYZ input and deterministic index ties instead of randomized inputs/set iteration','Integer clique labels instead of 26-character alphabet','Beam32 per root instead of exhaustive binary path tree','Positive denominator max(connection,1) for zero connections','Retain all finalist paths rather than discard last sorted path','Christofides metric TSP (paper section IV-B) instead of public code stochastic genetic TSP']}

def build_extra_schedules(design):
    start=time.monotonic();total=a.block(design['coefficients']);words=sorted(total,key=lambda w:a.dense(w,design['n_qubits']));coeff=[total[w] for w in words];ms=list(map(a.masks,words));adj=a.adjacency(ms)
    _,cons,hmeta=a.hastings(design);out={'hastings_interleaved_consolidated':{'schedule':cons,'metadata':{**hmeta,'consolidation':'Exact duplicate sum within each commuting family; step unitary equals literal Hastings schedule'}}}
    order,meta=tomesh(words,coeff,ms,adj);out['tomesh_beam32']={'schedule':[(words[i],coeff[i]) for i in order],'metadata':meta}
    comm,reverse=a.appendix_orders(coeff,adj)
    for name,order,rule in [('tranter_commutator',comm,'Minimize commuting prefix count; then maximum coefficient magnitude'),('tranter_reverse_commutator',reverse,'Maximize commuting remaining count; then maximum coefficient magnitude')]:
        out[name]={'schedule':[(words[i],coeff[i]) for i in order],'metadata':{'source':'Tranter arXiv:1912.07555v1 AppendixA.1--2','rule':rule,'ties':'dense IXYZ index','application':'first-order common-target state-overlap benchmark'}}
    for m in out.values():m['metadata']['joint_preparation_seconds']=time.monotonic()-start
    return out
