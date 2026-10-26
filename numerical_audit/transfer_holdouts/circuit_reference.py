"""Independent literal transfer replay and two-pass circuit DAG audit.

Adapted from the completed compiled_accuracy/reference.py: its sparse Pauli,
Taylor reference, unfused ideal formula and independent DAG calculations are
retained. Initial-state handling is extended explicitly for external non-HF
states; producer compiler/replay helpers are not imported.
"""
from array import array
from fractions import Fraction
from functools import lru_cache
import gzip,json,math
import numpy as np
from scipy import sparse

def require(ok,message):
    if not ok:raise ValueError(message)

def literal(word,n):
    require(type(n) is int and 1<=n<=8,'Independent register bound')
    require(len({q for q,p in word})==len(word) and all(type(q) is int and 0<=q<n and p in ('X','Y','Z') for q,p in word),'Invalid literal support')
    matrices={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    labels=dict(word);out=sparse.csr_matrix([[1.+0j]])
    for q in reversed(range(n)):out=sparse.kron(out,matrices[labels.get(q,'I')],format='csr')
    return out

class Reference:
    def __init__(self,model,initial):
        require(np.finfo(np.longdouble).nmant>=63,'Independent replay requires at least63mantissa bits')
        n=model['n'];require(type(n) is int and 1<=n<=6,'Register bound before allocation')
        require(isinstance(initial,np.ndarray) and initial.shape==(1<<n,),'Initial-state shape before allocation')
        require(np.count_nonzero(initial)<=9 and np.all(np.isfinite(initial)) and abs(np.vdot(initial,initial)-1)<1e-10,'Normalized bounded-support reference required')
        require(type(model['electrons']) is int and 0<=model['electrons']<=n,'Invalid occupation')
        require(0<len(model['coefficients'])<=128 and 0<len(model['blocks'])<=128,'Term/unit storage bound')
        den=int(model['denominator']);require(den>0,'Positive denominator required')
        self.model=model;self.n=n;self.dimension=1<<n;self.v=(1<<model['electrons'])-1;self.den=den
        self.words={};self.pool={};self.H=sparse.csr_matrix((self.dimension,self.dimension),dtype=complex)
        for raw,c in model['coefficients']:
            word=tuple((q,p) for q,p in raw);require(word and word not in self.words,'Scalar-removed unique terms required')
            require(type(c) in (str,int),'Integer numerator required');coefficient=int(c)
            P=literal(word,n);self.pool[word]=P;self.words[word]=coefficient;self.H+=(coefficient/den)*P
        supplied=[]
        for block in model['blocks']:
            require(block,'Empty complete unit');supplied += [(tuple((q,p) for q,p in w),int(c)) for w,c in block]
        require(len(supplied)==len(self.words) and dict(supplied)==self.words,'Target units do not partition Pauli terms')
        for block in model['blocks']:
            operators=[self.pool[tuple((q,p) for q,p in w)] for w,c in block]
            for at,P in enumerate(operators):
                for Q in operators[at+1:]:
                    commutator=P@Q-Q@P;commutator.eliminate_zeros()
                    require(commutator.nnz==0,'Independent unit commutation required for exact fusion')
        self.initial=np.asarray(initial,dtype=np.clongdouble).copy()
        self.single={};self.cx={}
        matrices={'h':np.array([[1,1],[1,-1]],dtype=np.clongdouble)/np.sqrt(np.longdouble(2)),'s':np.diag([1,1j]),'sdg':np.diag([1,-1j])}
        for kind,M in matrices.items():
            for q in range(n):self.single[kind,q]=self.embed(M,q)
        for a in range(n):
            for b in range(n):
                if abs(a-b)!=1:continue
                columns=np.arange(self.dimension);rows=np.array([v^((1<<b) if (v>>a)&1 else 0) for v in columns])
                self.cx[a,b]=sparse.csr_matrix((np.ones(self.dimension),(rows,columns)),shape=self.H.shape,dtype=complex)

    def embed(self,matrix,q):
        out=sparse.csr_matrix([[1.+0j]])
        for bit in reversed(range(self.n)):out=sparse.kron(out,matrix if bit==q else np.eye(2),format='csr')
        return out

    @lru_cache(maxsize=128)
    def rz(self,q,numerator,denominator):
        angle=Fraction(numerator,denominator);a=np.longdouble(angle.numerator)/np.longdouble(angle.denominator)/2
        require(bool(np.isfinite(a)),'Nonfinite Rz angle')
        phase=np.clongdouble(np.cos(a)-1j*np.sin(a))
        return self.embed(np.diag([phase,phase.conjugate()]),q)

    def exact_action(self,T):
        require(math.isfinite(T) and T>0,'Positive finite time required')
        bound=float(max(np.asarray(abs(self.H).sum(axis=1)).ravel()))*T
        steps=max(1,math.ceil(bound/.25));state=self.initial.copy()
        for _ in range(steps):
            term=state.copy();out=state.copy()
            for k in range(1,25):term=(-1j*T/steps/k)*(self.H@term);out+=term
            state=out
        return state

    def ideal_formula(self,record,steps,T):
        require(type(steps) is int and 1<=steps<=128 and record['formula'] in (1,2),'Fixed-grid ideal replay required')
        schedule=[(tuple((q,p) for q,p in w),int(c)) for w,c in record['schedule']]
        require(sum(record['unit_lengths'])==len(schedule),'Unit boundaries differ')
        offsets=np.cumsum([0]+record['unit_lengths']);state=self.initial.astype(np.clongdouble)
        def unit(index,duration):
            nonlocal state
            for at in range(offsets[index],offsets[index+1]):
                word,c=schedule[at];require(c==self.words[word],'Replay coefficient differs')
                theta=np.longdouble(c)/np.longdouble(self.den)*duration;state=np.cos(theta)*state-1j*np.sin(theta)*(self.pool[word]@state)
        m=len(record['unit_lengths'])
        for _ in range(steps):
            if record['formula']==1:
                for i in range(m):unit(i,np.longdouble(T)/steps)
            else:
                for i in range(m-1):unit(i,np.longdouble(T)/(2*steps))
                unit(m-1,np.longdouble(T)/steps)
                for i in reversed(range(m-1)):unit(i,np.longdouble(T)/(2*steps))
        return state

    def audit_compiled_trace(self,path,identity,compiler_version,simulate=True):
        """Build predecessor arrays, then compute longest paths independently.

        Primitive replay uses literal sparse matrices, not compiler helpers.
        Two-million-gate bounds keep these compact integer arrays below64MiB.
        """
        previous=[-1]*self.n;left=array('i');right=array('i');twoq=bytearray();one=two=0
        if identity['preparation']=='HF_prefix':
            require(identity['state']=='HF' and abs(self.initial[self.v]-1)<1e-14 and np.count_nonzero(self.initial)==1,'HF reference differs')
            state=np.zeros(self.dimension,np.clongdouble);state[0]=1
        else:
            require(identity['preparation']=='external_state_not_compiled' and identity['state'] in ('ground','excited','coherent_real','coherent_imag'),'External preparation declaration differs')
            state=self.initial.copy()
        require(identity['compilation_status']=='complete','Censored trace cannot be simulated as complete')
        with gzip.open(path,'rt') as f:
            header=json.loads(next(f));require(header['compiler']==compiler_version and header['n']==self.n and header['identity']==identity,'Compiled target/header differs')
            for index,line in enumerate(f):
                require(index<2_000_000 and len(line)<16384,'Independent trace storage bound')
                g=json.loads(line);require(isinstance(g,list) and len(g)>=2 and type(g[1]) is int and 0<=g[1]<self.n,'Invalid primitive trace')
                kind,q=g[:2]
                if kind=='cx':
                    require(len(g)==3 and type(g[2]) is int and 0<=g[2]<self.n and abs(q-g[2])==1,'Nonlocal primitive CX')
                    t=g[2];left.append(previous[q]);right.append(previous[t]);twoq.append(1);previous[q]=previous[t]=index;two+=1
                    if simulate:state=self.cx[q,t]@state
                else:
                    require(kind in ('h','s','sdg','rz'),'Unknown primitive')
                    left.append(previous[q]);right.append(-1);twoq.append(0);previous[q]=index;one+=1
                    if kind=='rz':
                        require(len(g)==4 and isinstance(g[2],str) and isinstance(g[3],str) and int(g[3])>0,'Invalid rational angle')
                        angle=Fraction(int(g[2]),int(g[3]));require(angle and str(angle.numerator)==g[2] and str(angle.denominator)==g[3],'Noncanonical/zero angle')
                        if simulate:state=self.rz(q,angle.numerator,angle.denominator)@state
                    else:
                        require(len(g)==2,'Unexpected Clifford payload')
                        if simulate:state=self.single[kind,q]@state
        depth=array('i');cx_depth=array('i')
        for index,(a,b,cx) in enumerate(zip(left,right,twoq)):
            require(a<index and b<index,'Cyclic circuit witness')
            d=max(depth[a] if a>=0 else 0,depth[b] if b>=0 else 0)
            c=max(cx_depth[a] if a>=0 else 0,cx_depth[b] if b>=0 else 0)
            depth.append(d+1);cx_depth.append(c+cx)
        costs=dict(gates=one+two,single_qubit_gates=one,two_qubit_gates=two,
            unit_duration_ASAP_depth=max(depth,default=0),CX_dependency_depth=max(cx_depth,default=0))
        return state,costs

def metrics(exact,actual,model):
    norm=float(np.sqrt(sum(abs(actual)**2)));require(math.isfinite(norm) and norm>0,'Invalid state norm')
    a=actual/norm;z=exact/np.sqrt(sum(abs(exact)**2));inner=np.vdot(z,a);residual=a-z*inner
    loss=float(sum(abs(residual)**2));n=model['n'];hf=(1<<model['electrons'])-1
    target=[sum((hf>>q)&1 for q in range(p,n,2)) for p in (0,1)];number=joint=0.
    for v in range(len(a)):
        if v.bit_count()!=model['electrons']:number+=abs(a[v])**2
        if [sum((v>>q)&1 for q in range(p,n,2)) for p in (0,1)]!=target:joint+=abs(a[v])**2
    return dict(error=float(loss/(1+abs(inner))),infidelity=loss,norm=norm,number_leakage=float(number),joint_sector_leakage=float(joint))
