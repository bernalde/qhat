"""Diagonal-reference variance for normalized complex superpositions.

No full-H diagonalization or dynamics is used in these features. Work is
bounded before allocating pair response vectors, including leaking destinations.
"""
import importlib.util
from pathlib import Path
import numpy as np

PREVIOUS=Path(__file__).resolve().parents[1]/'prospective_reference'
spec=importlib.util.spec_from_file_location('state_holdouts_determinant_engine',PREVIOUS/'engine.py')
determinant=importlib.util.module_from_spec(spec)
spec.loader.exec_module(determinant)


class Features:
    def __init__(self, blocks, diagonal, psi, T, denominator, n):
        if not 0<n<=8 or not 0<len(blocks)<=128 or not np.isfinite(T) or T<=0 or denominator<=0:
            raise ValueError('Unsupported response domain')
        dimension=1<<n
        psi=np.asarray(psi,dtype=complex)
        if psi.shape!=(dimension,) or not np.all(np.isfinite(psi)) or abs(np.vdot(psi,psi)-1)>1e-10:
            raise ValueError('Normalized full-space state required')
        for block in [diagonal]+list(blocks):
            for word,c in block.items():
                if len(set(q for q,p in word))!=len(word) or any(not 0<=q<n or p not in 'XYZ' for q,p in word):
                    raise ValueError('Invalid Pauli word')
        m=len(blocks)
        if m*(m-1)//2*dimension*16>64*1024**2:
            raise ValueError('Pair response storage exceeds 64 MiB')
        self.m=m
        self.destinations=list(range(dimension))
        self.pairs={}
        basis=np.arange(dimension,dtype=np.int64)
        descriptions=[determinant.descriptor(b) for b in blocks]
        if determinant.descriptor(diagonal)[0]!=0:
            raise ValueError('Reference must be diagonal')
        energies=np.array([determinant.value(determinant.descriptor(diagonal),int(v))/denominator
                           for v in basis])
        values=[np.array([determinant.value(d,int(v))/denominator for v in basis])
                for d in descriptions]
        for i in range(m):
            x=descriptions[i][0]
            for j in range(i+1,m):
                y=descriptions[j][0]
                destination=basis^x^y
                coefficient=-.5j*(values[j]*values[i][basis^y]
                                  -values[i]*values[j][basis^x])
                a=np.empty(dimension,complex)
                a[destination]=determinant.kernel(energies[destination]-energies,T)*coefficient*psi
                a-=psi*np.vdot(psi,a)
                ids=np.flatnonzero(a)
                if len(ids):self.pairs[i,j]=(ids,a[ids])
        self.storage_bytes=sum(ids.nbytes+a.nbytes for ids,a in self.pairs.values())

    amplitude=determinant.Features.amplitude
    crossing=determinant.Features.crossing
    score=staticmethod(determinant.Features.score)


def canonical_eigensystem(H, tolerance=1e-11):
    """Order eigenvectors deterministically within numerically degenerate spaces."""
    eigenvalues,vectors=np.linalg.eigh(H)
    start=0
    while start<len(eigenvalues):
        end=start+1
        while end<len(eigenvalues) and eigenvalues[end]-eigenvalues[start]<=tolerance:end+=1
        if end-start>1:
            space=vectors[:,start:end]
            projector=space@space.conj().T
            chosen=[]
            for column in projector.T:
                candidate=column.copy()
                # Reorthogonalize twice; tolerance here tests rank, not energy.
                for _ in range(2):
                    for old in chosen:candidate-=old*np.vdot(old,candidate)
                length=np.linalg.norm(candidate)
                if length>1e-8:chosen.append(candidate/length)
                if len(chosen)==end-start:break
            if len(chosen)!=end-start:raise ValueError('Degenerate-space rank failure')
            vectors[:,start:end]=np.column_stack(chosen)
        for at in range(start,end):
            pivot=np.argmax(abs(vectors[:,at]))
            phase=vectors[pivot,at]/abs(vectors[pivot,at])
            vectors[:,at]/=phase
        start=end
    return eigenvalues,vectors
