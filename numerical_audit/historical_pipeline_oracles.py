#!/usr/bin/env python3
"""Independent sparse-matrix checks of the archived molecular study pipeline."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-provenance-mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-provenance-numba')
import argparse,hashlib,importlib.metadata,json,pathlib,sys
sys.dont_write_bytecode=True
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from openfermion import get_fermion_operator,get_sparse_operator,jordan_wigner,QubitOperator

def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=pathlib.Path);p.add_argument('output',type=pathlib.Path);a=p.parse_args();sys.path[:0]=[str(a.study/'analysis'),str(a.study)]
 import benchmark_b2_active_spaces_matrix_free as mf
 import benchmark_b2_signed_coefficient_baseline as baseline
 assert pathlib.Path(mf.__file__).resolve().is_relative_to(a.study)
 tensor=a.study/'codex_external_validation/runs/transfer-vtrgy3ds/inputs/HCN/s-1.00/6-31g/HCN_s1.00_6-31g_as-004-006.tensors.npz'
 interaction,n=baseline.load_interaction_operator(tensor);f=baseline.clean_fermion_operator(get_fermion_operator(interaction),1e-12);jw=jordan_wigner(f);jw.compress(abs_tol=1e-12);coeff={k:v for k,v in jw.terms.items() if k and abs(v)>1e-12};keys=list(coeff)
 H_ferm=get_sparse_operator(f,n_qubits=n).tocsr();cF=float(complex(f.terms.get((),0)).real);cJW=float(complex(jw.terms.get((),0)).real)
 retained=QubitOperator();retained.terms=coeff.copy();H_retained=get_sparse_operator(retained,n_qubits=n).tocsr();H_ref=H_ferm-cF*sparse.eye(2**n,format='csr');alpha=cJW-cF
 action_difference=H_ref-H_retained-alpha*sparse.eye(2**n,format='csr');matrix_res=float(sparse.linalg.norm(action_difference))
 initial=mf.build_hartree_fock_state(n,4);T=.37
 exact_retained=expm_multiply(-1j*T*H_retained,initial,traceA=-1j*T*H_retained.diagonal().sum())
 exact_full=expm_multiply(-1j*T*H_ref,initial,traceA=-1j*T*H_ref.diagonal().sum())
 sector=[]
 for spin in (False,True):
  v,indices,_,_=mf.exact_reference_state(f,n,4,T,1e-12,spin);embedded=np.zeros(2**n,complex);embedded[indices]=v
  sector.append({'spin_preserving':spin,'sector_dimension':len(v),'max_norm_drift':max(abs(np.linalg.norm(v)-1),abs(np.linalg.norm(exact_full)-1)),'sector_embedding_vs_full_space_residual':float(np.linalg.norm(embedded-exact_full)),'retained_pauli_reference_after_scalar_phase_residual':float(np.linalg.norm(np.exp(1j*alpha*T)*embedded-exact_retained)),'raw_unaligned_vector_difference':float(np.linalg.norm(embedded-exact_retained))})
 schedules=baseline.build_deterministic_orderings(f,keys,coeff,n,1e-12);checks=[]
 for name in ('fermionic_signed_coefficient_lexicographic','signed_coefficient_lexicographic','jw_magnitude_descending_lexicographic'):
  ordered=schedules[name];terms=mf.compile_ordered_terms(ordered,coeff,n,1e-12);dt=.08
  direct=initial.copy();reverse=initial.copy()
  matrices={k:get_sparse_operator(QubitOperator(k,1),n_qubits=n).tocsr() for k in ordered}
  for k in ordered:
   angle=float(complex(coeff[k]).real)*dt;direct=np.cos(angle)*direct-1j*np.sin(angle)*(matrices[k]@direct)
  for k in reversed(ordered):
   angle=float(complex(coeff[k]).real)*dt;reverse=np.cos(angle)*reverse-1j*np.sin(angle)*(matrices[k]@reverse)
  actual,count=mf.evolve_trotter_state(initial,terms,1,1,dt,2**30)
  checks.append({'ordering':name,'sparse_matrix_oracle_residual':float(np.linalg.norm(actual-direct)),'forward_vs_reverse_distance':float(np.linalg.norm(direct-reverse)),'nominal_exponentials':count})
 # General Y-phase test independently uses explicit OpenFermion matrices.
 rng=np.random.default_rng(20261009);psi=rng.normal(size=16)+1j*rng.normal(size=16);psi/=np.linalg.norm(psi)
 generic={((0,'Y'),):.31,((0,'X'),(2,'Z')):-.27,((1,'Y'),(3,'Y')):.23,((0,'Z'),(1,'X'),(3,'Y')):.19}
 ordered=list(generic);terms=mf.compile_ordered_terms(ordered,generic,4,1e-12);direct=psi.copy();T=.43;r=3
 for _ in range(r):
  for k,c in generic.items():
   P=get_sparse_operator(QubitOperator(k),n_qubits=4);direct=np.cos(c*T/r)*direct-1j*np.sin(c*T/r)*(P@direct)
 actual,_=mf.evolve_trotter_state(psi,terms,1,r,T,2**30);general_res=float(np.linalg.norm(direct-actual))
 out={'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'tensor_path':str(tensor.relative_to(a.study)),'tensor_sha256':hashlib.sha256(tensor.read_bytes()).hexdigest(),'source_files':{pathlib.Path(m.__file__).name:hashlib.sha256(pathlib.Path(m.__file__).read_bytes()).hexdigest() for m in (mf,baseline)},'environment':{name:importlib.metadata.version(name) for name in ('numpy','scipy','openfermion','numba')},'qubits':n,'retained_pauli_terms':len(coeff),'fermionic_constant':cF,'jw_identity_coefficient':cJW,'reference_minus_retained_scalar':alpha,'matrix_identity_relation_frobenius_residual':matrix_res,'sector_checks':sector,'chronological_pauli_checks':checks,'generic_complex_state_with_odd_Y_oracle_residual':general_res,'interpretation':'Independent OpenFermion sparse-matrix action agrees with historical simulator and sector reference up to the known scalar offset for this frozen HCN input. This validates these conventions for the checked case, not all archived tensors or historical output provenance.'}
 assert matrix_res<1e-9 and all(x['sector_embedding_vs_full_space_residual']<1e-10 and x['retained_pauli_reference_after_scalar_phase_residual']<1e-10 for x in sector)
 assert general_res<1e-12 and all(x['sparse_matrix_oracle_residual']<1e-12 for x in checks)
 with a.output.open('x') as dest:json.dump(out,dest,indent=2);dest.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
