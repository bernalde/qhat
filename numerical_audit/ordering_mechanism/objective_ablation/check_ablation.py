#!/usr/bin/env python3
"""Independent exact target/schema/selection and objective checks; no QHAT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json,hashlib,math
from pathlib import Path
import numpy as np
def require(b,s):
    if not b:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def block(raw):return {tuple((int(q),p) for q,p in k):int(v) for k,v in raw}
def flip(k):return sum(1<<q for q,p in k if p in ('X','Y'))
def action(b,v):
    out=0;destination=None
    for word,c in b.items():
        w=v;real,imag=1,0
        for q,p in word:
            occupied=(v>>q)&1
            if p=='X':w^=1<<q
            elif p=='Y':
                w^=1<<q;real,imag=(imag,-real) if occupied else (-imag,real)
            elif p=='Z' and occupied:real,imag=-real,-imag
        require(imag==0,'Imaginary word action')
        require(destination is None or destination==w,'Mixed flip mask')
        destination=w;out+=c*real
    return destination,out
def objective(bs,order,hf,den,transported):
    xs=[flip(next(iter(b))) for b in bs];unique=set(xs);initial=[action(b,hf)[1] for b in bs]
    cached=[{x:action(b,hf^x)[1] for x in unique} for b in bs];amps={}
    for at,i in enumerate(order):
        for j in order[at+1:]:
            dest=xs[i]^xs[j];v=initial[j]*cached[i][xs[j]]-initial[i]*cached[j][xs[i]]
            require(dest!=0 or v==0,'Pair has nonzero diagonal HF expectation')
            amps[dest]=amps.get(dest,0)+v
    diag={k:v for b,x in zip(bs,xs) if x==0 for k,v in b.items()};energy0=action(diag,hf)[1]/den;terms=[]
    for dest,v in amps.items():
        if not v:continue
        gap=(action(diag,hf^dest)[1]/den-energy0)/2
        weight=(math.sin(gap)/gap)**2 if gap and transported else 1.
        terms.append(.25*weight*(v/(den*den))**2)
    return math.fsum(terms)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--phase2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    manifest=json.loads((a.run/'manifest.json').read_text());files=json.loads((a.run/'artifact_manifest.json').read_text())
    for rel,h in files.items():require(sha(a.run/rel)==h,'Changed artifact '+rel)
    require(sha(a.parent/'manifest.json')==manifest['parent_manifest_sha256'],'Parent manifest changed')
    require(sha(a.phase2/'manifest.json')==manifest['phase2_manifest_sha256'],'Phase2 manifest changed')
    checks=[];maxdelta=0;last_candidate=max(p.stat().st_mtime_ns for p in a.run.glob('*/*.candidate.json'));first_error=min(p.stat().st_mtime_ns for p in a.run.glob('*/*.r100.ablation.json'))
    require(last_candidate<first_error,'Candidate timestamps do not precede errors')
    for method in manifest['methods']:
        for case in manifest['cases']:
            d=json.loads((a.run/method/f'{case}.design.json').read_text());c=json.loads((a.run/method/f'{case}.candidate.json').read_text());p=json.loads((a.parent/f'{case}.design.json').read_text());fixed=json.loads((a.phase2/f'{case}.design.json').read_text())
            for key in ('coefficients','denominator','affine_basis_indices','tensor_sha256','target_sha256','initial_sha256','electrons','n_qubits'):
                require(d[key]==p[key]==fixed[key],'Unequal common target/state '+key)
            bs=[block(b) for b in d['blocks']['class']];order=d['orders'][method];total=block(d['coefficients']);seen={}
            require(set(d['orders'])=={method},'Inactive order contamination')
            require(sorted(order)==list(range(len(bs))),'Not unit permutation')
            for b in bs:
                for k,v in b.items():require(k not in seen,'Duplicate Pauli');seen[k]=v
            require(seen==total,'Wrong retained coefficients')
            if method=='class_refined_J0':require(d['blocks']['class']==fixed['blocks']['class'],'Changed refined partition')
            else:require(all(len(b)==1 for b in bs),'Pauli partition contains grouped units')
            row=json.loads((a.run/method/f'{case}.{method}.r100.ablation.json').read_text())
            require(row['T']==1 and row['steps']==100 and row['pauli_rotations']==100*len(total),'Changed time or work')
            den=int(d['denominator']);J=objective(bs,order,(1<<d['electrons'])-1,den,method=='class_pauli_JD');delta=abs(J-c['J']);maxdelta=max(maxdelta,delta)
            require(delta<1e-12*max(J,1e-20),'Independent objective mismatch')
            require(c['J']==row['surrogate_J'],'Recorded objective inconsistent')
            require(len(c['starts'])==2 and all(s['accepted_moves']<=10 for s in c['starts']),'Changed search budget')
            require(all(m['new_J']<m['old_J']-max(1e-30,1e-12*m['old_J']) for m in c['accepted_moves']),'Non-improving move')
            with np.load(a.parent/f'{case}.states.npz',allow_pickle=False) as ref,np.load(a.run/method/f'{case}.states.npz',allow_pickle=False) as saved:
                for key in ('initial','exact','exact_dense'):require(np.array_equal(ref[key],saved[key]),'Changed paired reference')
            witness_units=0
            if method=='class_pauli_JD':
                for b in bs:
                    dest,val=action(b,0)
                    if dest and val:witness_units+=1
            checks.append({'case':case,'method':method,'units':len(bs),'paulis':len(total),'independent_J':J,'J_absolute_difference':delta,
                           'exact_same_order_as_refined_JD':order==fixed['orders']['class_refined_response'] if method=='class_refined_J0' else None,
                           'pauli_units_with_vacuum_number_violation_witness':witness_units})
    result={'status':'passed','artifact_files':len(files),'checks':checks,'max_objective_absolute_difference':maxdelta,
            'all_candidates_precede_all_new_errors':True,'script_sha256':sha(Path(__file__)),'qhat_imports':False}
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
if __name__=='__main__':main()
