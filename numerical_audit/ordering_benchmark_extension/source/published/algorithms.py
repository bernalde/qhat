"""Independently implemented published ordering rules; exact integer coefficients.

Tranter arXiv:1912.07555v1 IV.3--4, V.2--3, Appendix A;
Hastings arXiv:1403.1539v2 section 4. See PROTOCOL.md for adaptations.
"""
from collections import defaultdict, deque
import hashlib,json,time

def block(rows): return {tuple((int(q),p) for q,p in w):int(c) for w,c in rows}
def masks(w):
    x=z=0
    for q,p in w:
        if p in 'XY':x|=1<<q
        if p in 'YZ':z|=1<<q
    return x,z

def dense(w,n):
    a=['I']*n
    for q,p in w:a[q]=p
    return ''.join(a)
def anti(a,b):return (((a[0]&b[1])^(a[1]&b[0])).bit_count()&1)!=0
def phase(a,b):return ((a[0]&a[1]).bit_count()+(b[0]&b[1]).bit_count()-((a[0]^b[0])&(a[1]^b[1])).bit_count()+2*(a[1]&b[0]).bit_count())%4

def adjacency(ms):
    adj=[0]*len(ms)
    for i,a in enumerate(ms):
        bit=1<<i
        for j in range(i):
            if anti(a,ms[j]):adj[i]|=1<<j;adj[j]|=bit
    return adj

def bits(a):
    while a:
        b=a&-a;yield b.bit_length()-1;a^=b

def independent_groups(adj):
    # NetworkX independent_set strategy: minimum degree in the remaining induced
    # eligible graph at each inner step; explicit lowest-label tie convention.
    remain=(1<<len(adj))-1;groups=[]
    while remain:
        eligible=remain;group=[]
        while eligible:
            i=min(bits(eligible),key=lambda j:((adj[j]&eligible).bit_count(),j))
            group.append(i);eligible&=~(adj[i]|(1<<i))
        groups.append(group)
        for i in group:remain^=1<<i
    return groups

def tranter_groups(coeff,adj):
    groups=independent_groups(adj);g=[deque(sorted(s,key=lambda i:(-abs(coeff[i]),i))) for s in groups]
    deplete=[]
    while any(g):
        for s in g:
            if s:deplete.append(s.popleft())
    g=[deque(sorted(s,key=lambda i:(-abs(coeff[i]),i))) for s in groups];equal=[]
    while any(g):
        k=min((i for i,s in enumerate(g) if s),key=lambda i:(-len(g[i]),-abs(coeff[g[i][0]]),g[i][0],i));equal.append(g[k].popleft())
    return deplete,equal,groups

def appendix_orders(coeff,adj):
    n=len(coeff);remain=(1<<n)-1;prefix=[]
    while remain:
        # Minimizing commuting prefix count is maximizing anticommuting count.
        selected=((1<<n)-1)^remain
        i=min(bits(remain),key=lambda j:(-(adj[j]&selected).bit_count(),-abs(coeff[j]),j))
        prefix.append(i);remain^=1<<i
    remain=(1<<n)-1;suffix=[]
    while remain:
        # Maximizing commuting suffix count equals minimizing its anti degree.
        i=min(bits(remain),key=lambda j:((adj[j]&remain).bit_count(),-abs(coeff[j]),j))
        suffix.append(i);remain^=1<<i
    return prefix,suffix

def hastings(design):
    parents=[block(b) for b in design['blocks']['parent']];meta=design['parent_metadata'];clusters=defaultdict(list)
    for i,(b,m) in enumerate(zip(parents,meta)):
        xs={masks(w)[0] for w in b}
        if len(xs)!=1:raise ValueError('Parent does not have a single flip mask')
        x=next(iter(xs));support=tuple(sorted({int(q) for comp in m['components'] for q,_ in comp}))
        if x==0:key=(0,())
        elif x.bit_count()==2:key=(1,tuple(bits(x)))
        elif x.bit_count()==4:key=(2,tuple(bits(x)))
        else:raise ValueError('Outside real one/two-body Hastings scope')
        clusters[key].append((len(support),support,i))
    literal=[];consolidated=[];sizes=[];indices=[]
    for key,items in sorted(clusters.items()):
        part=[];pids=[]
        for _,_,i in sorted(items):
            pids.append(i);part.extend(sorted(parents[i].items()))
        # Exact internal commutation permits combining duplicate exponentials.
        ms=[masks(w) for w,_ in part]
        if any(anti(a,b) for i,a in enumerate(ms) for b in ms[:i]):raise ValueError('Noncommuting Hastings cluster')
        merged=defaultdict(int)
        for w,c in part:merged[w]+=c
        literal.extend(part);consolidated.extend((w,c) for w,c in sorted(merged.items()) if c);sizes.append(len(part));indices.append(pids)
    total=defaultdict(int)
    for w,c in literal:total[w]+=c
    if {w:c for w,c in total.items() if c}!=block(design['coefficients']):raise ValueError('Hastings changes target')
    metadata={'source':'Hastings arXiv:1403.1539v2 section 4','adaptation':'Lexicographic orbital-family completion of unspecified final parallel-depth order; exact input parents retained','cluster_parent_indices':indices,'literal_cluster_sizes':sizes,'clusters_pairwise_commuting_exact':True}
    return literal,consolidated,metadata

def first_order_insertion(words,coeff,ms):
    order=[];poly={};norm=0;history=[]
    def update(key,val):
        nonlocal norm
        old=poly.get(key,0);new=old+val;norm+=abs(new)-abs(old)
        if new:poly[key]=new
        elif old:del poly[key]
    def comm(i,j):
        if not anti(ms[i],ms[j]):return None
        # Divide all anti-Hermitian commutators by i; real integer amplitudes.
        p=phase(ms[i],ms[j]);return (ms[i][0]^ms[j][0],ms[i][1]^ms[j][1]),(2 if p==1 else -2)*coeff[i]*coeff[j]
    for nxt in sorted(range(len(words)),key=lambda i:(-abs(coeff[i]),i)):
        for old in order:
            c=comm(old,nxt)
            if c:update(*c)
        best=norm;where=len(order)
        for pos in range(len(order)-1,-1,-1):
            c=comm(order[pos],nxt)
            if c:update(c[0],-2*c[1])
            if norm<=best:best=norm;where=pos
        # Restore from the fully left trial to the selected position.
        for old in order[:where]:
            c=comm(old,nxt)
            if c:update(c[0],2*c[1])
        order.insert(where,nxt);history.append([nxt,where,str(best)])
        if norm!=best:raise RuntimeError('BCH1 insertion restore failed')
    return order,history

def build_published_schedules(design,include_first_order_insertion=True):
    started=time.monotonic();total=block(design['coefficients']);words=sorted(total,key=lambda w:dense(w,design['n_qubits']));coeff=[total[w] for w in words];ms=list(map(masks,words));adj=adjacency(ms);dep,equ,groups=tranter_groups(coeff,adj)
    def emit(order,meta):return {'schedule':[(words[i],coeff[i]) for i in order],'metadata':meta}
    common={'source':'Tranter arXiv:1912.07555v1 IV.3--4','coloring':'greedy independent_set; minimum eligible induced degree, dense IXYZ lowest-index ties','groups':groups,'application':'first-order common-target state-overlap benchmark; original study used second-order energy error'}
    out={'tranter_deplete_groups':emit(dep,{**common,'rule':'cycle groups; maximum absolute coefficient within group'}),'tranter_equalise_groups':emit(equ,{**common,'rule':'largest remaining group; then maximum coefficient magnitude and dense index'})}
    lit,cons,meta=hastings(design);out['hastings_interleaved_lex']={'schedule':lit,'metadata':meta}
    if include_first_order_insertion:
        order,hist=first_order_insertion(words,coeff,ms);out['bch1_full_l1_insertion']=emit(order,{'source':'Explicit first-order adaptation inspired by Tranter V.3; not the published V2 objective','objective':'coefficient l1 of full first-order BCH commutator sum','incoming':'descending absolute coefficient; dense IXYZ ties','position_ties':'leftmost','history':hist})
    for m in out.values():m['metadata']['joint_preparation_seconds']=time.monotonic()-started
    return out
