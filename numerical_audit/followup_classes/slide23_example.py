#!/usr/bin/env python3
"""Actual four-parent JW example; independent occupation/dense circuit checks."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-class-mpl')
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from openfermion import FermionOperator
from exact_blocks import exact_parents
from run_class_comparison import direct_fermion_matrix
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics,write,sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    pairs=[(0,1,1),(2,3,2),(0,2,3),(1,3,5)];mats=[];terms=[];expansions=[]
    for p,q,w in pairs:
        raw=FermionOperator(((p,1),(q,0)),w)+FermionOperator(((q,1),(p,0)),w)
        parents,total,den,_=exact_parents(raw,0)
        mats.append(direct_fermion_matrix(parents,den,np.arange(16),4).toarray())
        coeff={k:v/den for k,v in total.items()}
        terms.append(evolution.compile_ordered_terms(list(coeff),coeff,4,0))
        expansions.append([[k,v] for k,v in coeff.items()])
    A=mats[0]+mats[1];B=mats[2]+mats[3];H=A+B;dt=.13
    contiguous=expm(-1j*dt*mats[3])@expm(-1j*dt*mats[2])@expm(-1j*dt*mats[1])@expm(-1j*dt*mats[0])
    grouped=expm(-1j*dt*B)@expm(-1j*dt*A)
    interleaved=expm(-1j*dt*mats[3])@expm(-1j*dt*mats[1])@expm(-1j*dt*mats[2])@expm(-1j*dt*mats[0])
    E=np.eye(16,dtype=complex)
    compiled=sum(terms,[])
    pauli=np.column_stack([evolve(E[:,i],compiled,dt,1) for i in range(16)])
    checks={'within_color_commutators':[float(np.linalg.norm(mats[0]@mats[1]-mats[1]@mats[0])),float(np.linalg.norm(mats[2]@mats[3]-mats[3]@mats[2]))],
        'between_color_commutator_frobenius':float(np.linalg.norm(A@B-B@A)),
        'grouped_minus_contiguous_frobenius':float(np.linalg.norm(grouped-contiguous)),
        'pauli_minus_contiguous_frobenius':float(np.linalg.norm(pauli-contiguous)),
        'interleaved_minus_contiguous_frobenius':float(np.linalg.norm(interleaved-contiguous))}
    if max(checks['within_color_commutators']+[checks['grouped_minus_contiguous_frobenius'],checks['pauli_minus_contiguous_frobenius']])>1e-12:raise RuntimeError('Exact factorization check failed.')
    if checks['between_color_commutator_frobenius']<1 or checks['interleaved_minus_contiguous_frobenius']<1e-5:raise RuntimeError('Ordering witness degenerated.')
    # Fixed a priori one-particle input |1000>; no state selected by observed ranking.
    initial=E[:,8];nm=np.array([i.bit_count()==1 for i in range(16)]);exact=expm(-1j*H)@initial
    rows=[];states={'initial':initial,'exact':exact}
    for r in (1,10,100):
        for name,order in (('grouped_contiguous',(0,1,2,3)),('interleaved',(0,2,1,3))):
            approx=evolve(initial,[term for i in order for term in terms[i]],1.,r)
            rows.append({'ordering':name,'T':1.,'steps':r,'rotations':8*r,**physical_metrics(exact,approx,nm,nm)})
            states[f'{name}_r{r}']=approx
    write(args.output/'result.json',{'parents':pairs,'weighted_JW_expansions':expansions,'checks':checks,'runs':rows,
        'source_sha256':{str(Path(__file__)):sha(__file__),str(Path(__file__).with_name('exact_blocks.py')):sha(Path(__file__).with_name('exact_blocks.py'))},
        'state':'one particle |1000>, fixed before execution','product_convention':'chronological; later factors multiply left',
        'interpretation':'Pedagogical operator and implementation example, not molecular accuracy evidence. Grouping contiguous commuting parents changes no circuit.'})
    np.savez_compressed(args.output/'states.npz',**states)
    print(json.dumps(checks),flush=True)

if __name__=='__main__':main()
