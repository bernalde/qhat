#!/usr/bin/env python3
"""Frozen exploratory diagonal-reference class ordering; see PROTOCOL.md."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/ordering-response-numba')
os.environ.setdefault('MPLCONFIGDIR','/tmp/ordering-response-mpl')
import argparse,json,hashlib,shutil,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics

def require(b,msg):
 if not b:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def block(raw):return {tuple((int(q),p) for q,p in word):int(v) for word,v in raw}
def masks(word):
 x=z=y=0
 for q,p in word:
  if p in ('X','Y'):x|=1<<q
  if p in ('Z','Y'):z|=1<<q
  y+=p=='Y'
 return x,z,y
def fval(b,v):
 total=0
 for word,c in b.items():
  _,z,y=masks(word);require(y%2==0,'Non-real class')
  total += c*(-1 if ((z&v).bit_count()+y//2)%2 else 1)
 return total
def reversebits(v,n):return int(format(int(v),f'0{n}b')[::-1],2)
def pair_setup(bs,xs,hf,den):
 m=len(bs); f0=[fval(b,hf) for b in bs]
 f=[[fval(b,hf^y) for y in xs] for b in bs]
 pairs={}
 for i in range(m):
  for j in range(i+1,m):
   val=f0[j]*f[i][j]-f0[i]*f[j][i]
   if val:pairs[i,j]=(xs[i]^xs[j],val)
 diag={}
 for b,x in zip(bs,xs):
  if x==0:
   for word,value in b.items():diag[word]=diag.get(word,0)+value
 energy0=fval(diag,hf)/den
 rho={dest:float(np.sinc(((fval(diag,hf^dest)/den-energy0)/2)/np.pi)**2) for dest,_ in pairs.values()}
 return pairs,rho

def amplitudes(order,pairs):
 pos={i:p for p,i in enumerate(order)};out={}
 for (i,j),(dest,b) in pairs.items():out[dest]=out.get(dest,0)+(b if pos[i]<pos[j] else -b)
 return out

def cost(amps,rho,den):return .25*sum(rho[k]*(v/(den*den))**2 for k,v in amps.items())
def optimize(start,pairs,rho,den):
 order=list(start);a=amplitudes(order,pairs);value=cost(a,rho,den);history=[]
 for it in range(10):
  best=value;best_order=None;best_move=None
  for pos in range(len(order)):
   for direction in (-1,1):
    trial=dict(a);trial_value=value;working=list(order);at=pos
    while 0<=at+direction<len(order):
     left=min(at,at+direction);i,j=working[left:left+2];item=pairs.get((min(i,j),max(i,j)))
     if item:
      dest,b=item;b=b if i<j else -b;old=trial.get(dest,0);new=old-2*b
      trial_value += .25*rho[dest]*((new/(den*den))**2-(old/(den*den))**2);trial[dest]=new
     working[left:left+2]=[j,i];at+=direction
     if trial_value<best-max(1e-30,1e-12*best):best=trial_value;best_order=list(working);best_move=(pos,at)
  if best_order is None:break
  final=amplitudes(best_order,pairs);actual=cost(final,rho,den)
  require(abs(actual-best)<1e-10*max(value,actual,1e-20),'Insertion delta disagrees with full objective')
  require(actual<value-max(1e-30,1e-12*value),'Accepted non-improvement')
  history.append({'pass':it,'move':best_move,'old_J':value,'new_J':actual});order=best_order;a=final;value=actual
 return order,value,history

def actual_block_matrices(bs,den,basis,n):
 little=[reversebits(v,n) for v in basis];idx={v:i for i,v in enumerate(little)}
 result=[]
 for b in bs:
  flip=masks(next(iter(b)))[0];src=np.array([idx[v^flip] for v in little]);fv=np.array([fval(b,v)/den for v in little])
  result.append((src,fv))
 return result

def bch_matrix(order,bms):
 dim=len(bms[0][0]);B=np.zeros((dim,dim));cols=np.arange(dim)
 for at,i in enumerate(order):
  si,fi=bms[i]
  for j in order[at+1:]:
   sj,fj=bms[j];amp=fj*fi[sj]-fi*fj[si]
   B[si[sj],cols]+=amp
 return B

def transport(K,vals,vecs,initial):
 d=vals[:,None]-vals[None,:];integral=np.exp(.5j*d)*np.sinc(d/(2*np.pi))
 q=vecs@(((vecs.conj().T@K@vecs)*integral)@(vecs.conj().T@initial))
 projection=q-initial*np.vdot(initial,q)
 return float(np.vdot(projection,projection).real),q

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 a.output.mkdir(parents=True,exist_ok=False);began=time.monotonic()
 cases=json.loads((a.parent_run/'manifest.json').read_text())['cases']
 src=Path(__file__);protocol=src.with_name('PROTOCOL.md')
 hashes={src.name:sha(src),protocol.name:sha(protocol)}
 for p in (src,protocol):shutil.copyfile(p,a.output/p.name)
 write(a.output/'manifest.json',{'cases':cases,'source_sha256':hashes,'parent_manifest_sha256':sha(a.parent_run/'manifest.json'),'expected_rows':20,'status':'exploratory algorithm development, not held-out'})
 prepared=[]
 # Freeze every case's candidate before measuring any new circuit error.
 for case in cases:
  d=json.loads((a.parent_run/f'{case}.design.json').read_text());den=int(d['denominator']);bs=[block(b) for b in d['blocks']['class']];xs=d['class_flip_masks'];hf=(1<<d['electrons'])-1
  oldbs,oldxs=bs,xs;bs=[];xs=[]
  for oldb,oldx in zip(oldbs,oldxs):
   for newb in ([{word:value} for word,value in sorted(oldb.items())] if oldx==0 else [oldb]):bs.append(newb);xs.append(oldx)
  d['blocks']['class']=[[[word,str(value)] for word,value in sorted(b.items())] for b in bs];d['class_flip_masks']=xs
  d['orders']['class_lex']=list(range(len(bs)))
  d['unit_partition']='individual diagonal Paulis; complete nonzero flip-mask classes'
  pairs,rho=pair_setup(bs,xs,hf,den);lex=list(range(len(bs)));mag=sorted(lex,key=lambda i:(-sum(abs(v) for v in bs[i].values()),xs[i]))
  candidates=[optimize(o,pairs,rho,den) for o in (lex,mag)]
  winner=candidates[0]
  for c in candidates[1:]:
   tolerance=max(1e-30,1e-12*min(winner[1],c[1]))
   if c[1]<winner[1]-tolerance or (abs(c[1]-winner[1])<=tolerance and tuple(c[0])<tuple(winner[0])):winner=c
  order,J,history=winner
  # Exact integer covariance check; all same-destination amplitudes change by same sign.
  for spatial_signs in (1,2,3,5,10):
   sm=sum(3<<(2*p) for p in range(d['n_qubits']//2) if spatial_signs>>p&1)
   signs=[-1 if (sm&x).bit_count()%2 else 1 for x in xs]
   gb=[{w:signs[i]*v for w,v in b.items()} for i,b in enumerate(bs)]
   gp,gr=pair_setup(gb,xs,hf,den)
   require(gr==rho,'Gauge changed diagonal energies')
   require(set(gp)==set(pairs),'Gauge changed nonzero pair set')
   for ij,(dest,b) in pairs.items():require(gp[ij]==(dest,b*(-1 if (sm&dest).bit_count()%2 else 1)),'Gauge covariance failed')
   for start in (lex,mag):require(optimize(start,gp,gr,den)[:2]==optimize(start,pairs,rho,den)[:2],'Gauge changed deterministic optimized order')
  require(cost(amplitudes(list(reversed(order)),pairs),rho,den)==J,'Reversal symmetry failed')
  diag={'case':case,'parent_design_sha256':sha(a.parent_run/f'{case}.design.json'),'J_D':J,'J_D_lex':cost(amplitudes(lex,pairs),rho,den),'J_D_magnitude':cost(amplitudes(mag,pairs),rho,den),'order':order,'magnitude_order':mag,'accepted_swaps':history,'candidate_start_results':[{'order':c[0],'J':c[1],'accepted_swaps':len(c[2])} for c in candidates],'gauge_controls':5,'preparation_seconds':time.monotonic()-began}
  write(a.output/f'{case}.candidate.json',diag)
  d['orders'].update(class_refined_response=order,class_refined_magnitude=mag)
  d['planned_runs']=[{'ordering':name,'steps':r,'role':'development' if r==100 else 'convergence','T':1.} for name in ('class_refined_response','class_refined_magnitude') for r in ([100,50,200] if case=='HCN_s1.00' else [100])]
  write(a.output/f'{case}.design.json',d);prepared.append((case,d,diag,bs,pairs,rho))
 rows=[];diagnostics=[];oldrows=json.loads((a.parent_run/'results.json').read_text())
 for case,d,diag,bs,pairs,rho in prepared:
  den=int(d['denominator']);n=d['n_qubits'];total=block(d['coefficients']);keys=sorted(set(total).union(*(block(b) for b in d['blocks']['parent'])))
  unitraw=ev.compile_ordered_terms(keys,{k:1. for k in keys},n,0.);basis,reduced,_=ref.affine_space(unitraw,d['affine_basis_indices'][0]);require(list(map(int,basis))==d['affine_basis_indices'],'Changed closure');unit=dict(zip(keys,reduced))
  def compile_word(k,v):return (v/den*unit[k][0],unit[k][1],unit[k][2],unit[k][3])
  H,_=ref.assemble_hamiltonian([compile_word(k,v) for k,v in total.items()],len(basis))
  with np.load(a.parent_run/f'{case}.states.npz') as saved:initial=saved['initial'];exact=saved['exact'];dense=saved['exact_dense']
  nr=np.array([int(v).bit_count()==d['electrons'] for v in basis]);sector=set(map(int,ev.spin_sector_basis_indices(n,d['electrons'])));sr=np.array([int(v) in sector for v in basis])
  states={'initial':initial,'exact':exact,'exact_dense':dense}
  for run in d['planned_runs']:
   name,r,role=run['ordering'],run['steps'],run['role'];order=d['orders'][name];terms=[compile_word(k,v) for i in order for k,v in bs[i].items()]
   t=time.monotonic();actual=evolve(initial,terms,1.,r)
   row={'case':case,'ordering':name,'T':1.,'steps':r,'role':role,'pauli_rotations':len(terms)*r,'seconds':time.monotonic()-t,'surrogate_J':cost(amplitudes(order,pairs),rho,den),**physical_metrics(exact,actual,nr,sr)}
   for baseline in ('class_lex','owned_signed','jw_signed','jw_magnitude'):
    old=next(x for x in oldrows if x['case']==case and x['ordering']==baseline and x['steps']==r)
    row['ratio_over_'+baseline]=row['one_minus_overlap']/old['one_minus_overlap'];row[baseline+'_E']=old['one_minus_overlap']
   write(a.output/f'{case}.{name}.r{r}.{role}.json',row);rows.append(row);states[f'{name}__r{r}__{role}']=actual
   print('EVOLVE',case,name,r,row['one_minus_overlap'],row['ratio_over_jw_signed'],flush=True)
  # Diagnostic only, after all candidates have been saved, never used for selection.
  vals,vecs=eigh(H.toarray());bms=actual_block_matrices(bs,den,basis,n);index={int(v):i for i,v in enumerate(basis)};diagrow={'case':case,'orders':{}}
  for name in ('class_lex','class_refined_magnitude','class_refined_response'):
   order=d['orders'][name];B=bch_matrix(order,bms);require(np.linalg.norm(B+B.T)<1e-11,'B not skew');K=.5j*B
   expected=np.zeros(len(basis));amps=amplitudes(order,pairs)
   for mask,b in amps.items():expected[index[int(basis[0])^reversebits(mask,n)]]=b/(den*den)
   discrepancy=float(np.linalg.norm(B@initial-expected));require(discrepancy<1e-11,'Pair amplitudes differ from full matrices')
   JH,q=transport(K,vals,vecs,initial);JD=cost(amps,rho,den);J0=.25*float(np.dot(expected,expected))
   diagrow['orders'][name]={'J_H':JH,'J_D':JD,'J_initial':J0,'leading_E_r100':JH/20000,'matrix_HF_pair_discrepancy':discrepancy}
  diagnostics.append(diagrow);write(a.output/f'{case}.transport_diagnostic.json',diagrow)
  np.savez_compressed(a.output/f'{case}.states.npz',**states)
 write(a.output/'results.json',rows);write(a.output/'diagnostics.json',diagnostics)
 summary={}
 for name in ('class_refined_response','class_refined_magnitude'):
  primary=[x for x in rows if x['ordering']==name and x['steps']==100];report={}
  for baseline in ('class_lex','owned_signed','jw_signed','jw_magnitude'):
   outcomes={'win':0,'tie':0,'loss':0}
   for row in primary:outcomes['tie' if abs(row['one_minus_overlap']-row[baseline+'_E'])<=1e-12 else ('win' if row['one_minus_overlap']<row[baseline+'_E'] else 'loss')]+=1
   report[baseline]={'outcomes':outcomes,'median_ratio':float(np.median([r['ratio_over_'+baseline] for r in primary]))}
  summary[name]=report
 write(a.output/'summary.json',summary)
 require(hashes=={src.name:sha(src),protocol.name:sha(protocol)},'Source changed during run')
 write(a.output/'completion.json',{'rows':len(rows),'seconds':time.monotonic()-began,'status':'complete'})
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
