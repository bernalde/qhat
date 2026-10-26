"""Independent nested-commutator cubic coefficient from literal Pauli matrices."""
import math
import numpy as np
from scipy import sparse


def require(ok,message):
    if not ok:raise ValueError(message)


def literal(word,n):
    require(type(n) is int and 0<n<=8,'Literal register bound')
    require(len(dict(word))==len(word) and all(type(q) is int and 0<=q<n and p in 'XYZ' for q,p in word),'Invalid literal word')
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    out=sparse.csr_matrix([[1.+0j]]);labels=dict(word)
    for q in reversed(range(n)):out=sparse.kron(out,matrices[labels.get(q,'I')],format='csr')
    return out


class Reference:
    def __init__(self,model):
        n=model['n'];require(type(n) is int and 0<n<=8,'Independent register bound before allocation')
        require(0<len(model['blocks'])<=128 and 0<len(model['coefficients'])<=128,'Independent unit/term bound')
        require(type(model['electrons']) is int and 0<=model['electrons']<=n,'Independent occupation bound')
        den=int(model['denominator']);require(den>0,'Positive denominator')
        self.model=model;self.dimension=1<<n;self.v=(1<<model['electrons'])-1
        self.initial=np.zeros(self.dimension,complex);self.initial[self.v]=1
        self.pool={tuple((q,p) for q,p in w):literal(w,n) for w,c in model['coefficients']}
        self.blocks=[]
        for raw in model['blocks']:
            require(0<len(raw)<=128,'Independent per-unit word bound')
            G=sparse.csr_matrix((self.dimension,self.dimension),dtype=complex)
            for w,c in raw:G+=(int(c)/den)*self.pool[tuple((q,p) for q,p in w)]
            self.blocks.append(G)
        self.H=sum(self.blocks);self.lengths=[len(b) for b in model['blocks']]
        self.gaps=self.H.diagonal().real-self.H[self.v,self.v].real
        total={tuple((q,p) for q,p in w):int(c) for w,c in model['coefficients']}
        self.paulis=[(c/den)*self.pool[w] for w,c in sorted(total.items())]

    def conservation(self):
        n=self.model['n'];dimension=self.dimension
        charges=[sparse.diags([sum((v>>q)&1 for q in range(p,n,2)) for v in range(dimension)],dtype=float) for p in (0,1)]
        for G,raw in zip(self.blocks,self.model['blocks']):
            for N in charges:
                comm=G@N-N@G;require(not comm.nnz or max(abs(comm.data))<2e-12,'Unit spin-sector conservation differs')
            words=[tuple((q,p) for q,p in w) for w,c in raw]
            for at,w in enumerate(words):
                for z in words[at+1:]:require((self.pool[w]@self.pool[z]-self.pool[z]@self.pool[w]).nnz==0,'Unit flattening is not internally commuting')

    def cubic(self,order):
        require(all(type(i) is int for i in order) and sorted(order)==list(range(len(self.blocks))),'Independent integral order required')
        # Recursing through an outer A/2, inner symmetric tail B, outer A/2
        # adds -[A,[A,B]]/24 + [B,[B,A]]/12 to the local i*h^3 defect.
        B=sparse.csr_matrix(self.H.shape,dtype=complex);out=np.zeros(self.dimension,complex);v=self.initial
        for index in reversed(order):
            A=self.blocks[index]
            Av=A@v;Bv=B@v
            AAB=A@(A@Bv)-2*A@(B@Av)+B@(A@Av)
            BBA=B@(B@Av)-2*B@(A@Bv)+A@(B@Bv)
            out+=-AAB/24+BBA/12;B+=A
        require(max(abs(out.imag))<2e-12,'Real-unit cubic coefficient became complex')
        return out.real

    def score(self,order,T,rotation_cap):
        coefficient=self.cubic(order);weights=np.empty(self.dimension)
        for w,gap in enumerate(self.gaps):
            x=T*float(gap)/2;weights[w]=T*T*((math.sin(x)/x)**2 if x else 1.)
        weights[self.v]=0;JD2=float(sum(weights*coefficient**2))
        lengths=[self.lengths[i] for i in order];U=sum(lengths)
        steps=(1 if rotation_cap>=U else 0) if len(lengths)==1 else (rotation_cap-lengths[0])//(2*U-lengths[0]-lengths[-1])
        require(steps>=1,'Independent affordability failure')
        return dict(JD2=JD2,score=JD2/(steps**4),steps=steps,cubic_coefficients=coefficient.tolist())

    def first_score(self,order,T,rotation_cap,partition):
        blocks=self.blocks if partition=='conserving' else self.paulis
        require(all(type(i) is int for i in order) and sorted(order)==list(range(len(blocks))),'Independent control permutation required')
        prefix=sparse.csr_matrix(self.H.shape,dtype=complex);column=np.zeros(self.dimension,complex)
        for i in order:
            G=blocks[i];column+=(-.5j)*(prefix@(G@self.initial)-G@(prefix@self.initial));prefix+=G
        column[self.v]=0
        for w,gap in enumerate(self.gaps):
            x=T*float(gap)/2;column[w]*=T*(math.sin(x)/x if x else 1.)
        value=float(np.vdot(column,column).real);steps=rotation_cap//len(self.model['coefficients'])
        return dict(JD1=value,score=value/(steps**2),steps=steps)

    def exact_action(self,T):
        bound=float(max(np.asarray(abs(self.H).sum(axis=1)).ravel()))*T
        steps=max(1,math.ceil(bound/.25));state=self.initial.copy()
        for _ in range(steps):
            term=state.copy();result=state.copy()
            for k in range(1,25):term=(-1j*T/steps/k)*(self.H@term);result+=term
            state=result
        return state

    def replay(self,record,steps,T):
        words=[(tuple((q,p) for q,p in w),int(c)/int(self.model['denominator'])) for w,c in record['schedule']]
        offsets=np.cumsum([0]+record['unit_lengths']);a=self.initial.copy()
        def apply(unit,duration):
            nonlocal a
            for k in range(offsets[unit],offsets[unit+1]):
                word,c=words[k];theta=c*duration;a=np.cos(theta)*a-1j*np.sin(theta)*(self.pool[word]@a)
        m=len(record['unit_lengths'])
        for _ in range(steps):
            if record['formula']==1:
                for unit in range(m):apply(unit,T/steps)
            else:
                for unit in range(m-1):apply(unit,T/(2*steps))
                apply(m-1,T/steps)
                for unit in reversed(range(m-1)):apply(unit,T/(2*steps))
        return a
