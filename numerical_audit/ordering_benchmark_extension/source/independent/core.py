"""Independent occupation-basis target and circuit calculations; no QHAT imports.

Production uses a compiled Pauli routine from QHAT. This verifier derives the
physical q0-MSB action directly, checks the affine coordinate system, then uses
in-place two-coordinate rotations. Symmetric formulas are replayed without
production's commuting-boundary merges. This is a numerical check, not a proof.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import hashlib,json,math
from pathlib import Path
import numpy as np
from numba import njit
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply


def require(condition,message):
    if not condition:raise ValueError(message)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def object_sha(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def key(word):return tuple((int(q),str(p)) for q,p in word)
def block(raw):return {key(k):int(v) for k,v in raw}
def physical_masks(word,n):
    x=z=y=0
    for q,p in word:
        require(0<=q<n and p in 'XYZ','Bad Pauli factor')
        if p in 'XY':x^=1<<(n-1-q)
        if p in 'ZY':z^=1<<(n-1-q)
        y+=p=='Y'
    return x,z,y

class OccupationAction:
    def __init__(self,design):
        self.n=int(design['n_qubits']);self.basis=np.asarray(design['affine_basis_indices'],dtype=np.int64)
        self.dim=len(self.basis);self.den=int(design['denominator']);self.total=block(design['coefficients'])
        require(self.dim>0 and self.dim&(self.dim-1)==0,'Dimension is not power of two')
        require(len(set(map(int,self.basis)))==self.dim,'Duplicate affine coordinates')
        self.translation=int(self.basis[0]);self.generators=[int(self.basis[1<<j])^self.translation for j in range(self.dim.bit_length()-1)]
        expected=np.full(self.dim,self.translation,dtype=np.int64)
        for j,g in enumerate(self.generators):expected[(np.arange(self.dim)&(1<<j))!=0]^=g
        require(np.array_equal(expected,self.basis),'Saved affine enumeration is inconsistent')
        self.positions={int(v):i for i,v in enumerate(self.basis)}
        self.parity=np.array([v.bit_count()&1 for v in range(1<<self.n)],dtype=np.int8)
        self._actions={}
    def action(self,word):
        word=key(word)
        if word not in self._actions:
            x,z,y=physical_masks(word,self.n);require(y%2==0,'Odd Y: real-only replay contract')
            require((self.translation^x) in self.positions,'Word leaves saved closure')
            flip=self.positions[self.translation^x]
            require(np.array_equal(self.basis[np.arange(self.dim)^flip],self.basis^x),'Nonuniform affine translation')
            signmask=sum(((g&z).bit_count()&1)<<j for j,g in enumerate(self.generators))
            sign0=-1 if ((self.translation&z).bit_count()+y//2)&1 else 1
            require((signmask&flip).bit_count()%2==0,'Hermitian real partner signs disagree')
            self._actions[word]=(flip,signmask,sign0)
        return self._actions[word]
    def hamiltonian(self):
        # Sum equal physical flips before sparse assembly; this avoids one
        # full set of COO coordinates per individual Pauli word on NH3.
        grouped={}
        for word,c in self.total.items():
            x,z,y=physical_masks(word,self.n);self.action(word)
            phase=(1-2*self.parity[self.basis&z])*(-1 if (y//2)&1 else 1)
            if x not in grouped:grouped[x]=np.zeros(self.dim)
            grouped[x]+=(c/self.den)*phase
        rows=[];cols=[];values=[];indices=np.arange(self.dim,dtype=np.int32)
        for x,amplitudes in grouped.items():
            keep=np.flatnonzero(amplitudes)
            if not len(keep):continue
            flip=self.positions[self.translation^x]
            rows.append((indices[keep]^flip).astype(np.int32));cols.append(indices[keep]);values.append(amplitudes[keep])
        H=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(self.dim,self.dim)).tocsr()
        delta=H-H.T
        require(not delta.nnz or np.max(np.abs(delta.data))<2e-12,'Independent target is not Hermitian')
        return H
    def reference(self,initial,T=1.):
        H=self.hamiltonian()
        if self.dim<=256:return expm(-1j*T*H.toarray())@initial,H.nnz
        return expm_multiply(-1j*T*H,initial,traceA=-1j*T*H.diagonal().sum()),H.nnz
    def prepare(self,schedule):
        actions=np.asarray([self.action(key(k)) for k,_ in schedule],dtype=np.int64)
        coefficients=np.asarray([int(v)/self.den for _,v in schedule],dtype=float)
        require(len(actions)>0,'Empty schedule')
        return coefficients,actions
    def replay(self,initial,schedule,steps,formula_order=1,T=1.):
        require(steps>0,'Cannot replay zero steps');coefficients,actions=self.prepare(schedule)
        state=np.array(initial,dtype=complex,copy=True)
        evolve_real(state,coefficients,actions,steps,formula_order,T)
        return state
    def metrics(self,exact,actual,electrons):
        a=exact/np.linalg.norm(exact);b=actual/np.linalg.norm(actual)
        overlap=np.vdot(a,b);phase=overlap/abs(overlap) if abs(overlap) else 1.
        diff=b-phase*a;error=.5*float(np.vdot(diff,diff).real)
        weights=abs(b)**2;numbers=np.array([int(v).bit_count() for v in self.basis])
        nup=np.array([sum((int(v)>>(self.n-1-q))&1 for q in range(0,self.n,2)) for v in self.basis])
        ndown=numbers-nup
        initial_nup=(int(electrons)+1)//2;initial_ndown=int(electrons)//2
        return {'one_minus_overlap':error,'infidelity':error*(2-error),
          'particle_number_leakage':float(weights[numbers!=int(electrons)].sum()),
          'spin_sector_leakage':float(weights[(nup!=initial_nup)|(ndown!=initial_ndown)].sum()),
          'norm':float(np.linalg.norm(actual)),
          'particle_number_expectation':float(np.dot(weights,numbers))}

@njit(cache=True)
def bitparity(x):
    p=0
    while x:
        p^=1;x&=x-1
    return p

@njit(cache=True)
def rotate_real(state,flip,mask,sign0,cosine,sine):
    if flip==0:
        for i in range(len(state)):
            sign=sign0*(-1 if bitparity(i&mask) else 1)
            state[i]*=cosine-1j*sine*sign
    else:
        for i in range(len(state)):
            j=i^flip
            if i<j:
                sign=sign0*(-1 if bitparity(i&mask) else 1)
                a=state[i];b=state[j]
                state[i]=cosine*a-1j*sine*sign*b
                state[j]=cosine*b-1j*sine*sign*a

@njit(cache=True)
def evolve_real(state,coefficients,actions,steps,formula_order,T):
    scale=T/steps/(2 if formula_order==2 else 1)
    cosines=np.cos(coefficients*scale);sines=np.sin(coefficients*scale)
    for _ in range(steps):
        for j in range(len(coefficients)):
            rotate_real(state,actions[j,0],actions[j,1],actions[j,2],cosines[j],sines[j])
        if formula_order==2:
            for j in range(len(coefficients)-1,-1,-1):
                rotate_real(state,actions[j,0],actions[j,1],actions[j,2],cosines[j],sines[j])


def schedule_sum(schedule):
    result={}
    for raw,v in schedule:
        word=key(raw);result[word]=result.get(word,0)+int(v)
    return {k:v for k,v in result.items() if v}


def number_commutator(raw):
    # The common factor i in [G,N] is omitted; every coefficient below is exact.
    result={}
    for word,c in block(raw).items():
        for index,(q,p) in enumerate(word):
            if p in 'XY':
                transformed=word[:index]+((q,'Y' if p=='X' else 'X'),)+word[index+1:]
                result[transformed]=result.get(transformed,0)+(c if p=='X' else -c)
    return {k:v for k,v in result.items() if v}
