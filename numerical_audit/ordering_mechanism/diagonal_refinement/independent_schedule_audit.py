#!/usr/bin/env python3
"""Read-only exact schedule/number-conservation/objective audit, no QHAT imports."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json,hashlib,math
from pathlib import Path
import numpy as np

def require(b,s):
    if not b:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def parse(raw):
    out={}
    for word,value in raw:
        k=tuple((int(q),str(p)) for q,p in word)
        require(k not in out,'Repeated word in serialized block')
        out[k]=int(value)
    return out
def flip(word):return sum(1<<q for q,p in word if p in ('X','Y'))
def direct_word_action(word,v):
    dest=v;real,imag=1,0
    for q,p in word:
        occupied=(v>>q)&1
        if p=='X':dest^=1<<q
        elif p=='Y':
            dest^=1<<q
            real,imag=(imag,-real) if occupied else (-imag,real)
        elif p=='Z':
            if occupied:real,imag=-real,-imag
        else:raise RuntimeError('Invalid Pauli')
    require(imag==0,'Block contains imaginary local Pauli action')
    return dest,real
def unit_action(block,v):
    dest=None;value=0
    for word,c in block.items():
        w,s=direct_word_action(word,v)
        if dest is None:dest=w
        require(dest==w,'Unit contains different flip masks')
        value+=c*s
    return dest,value
def commute_words(a,b):
    d=dict(b)
    return sum(q in d and d[q]!=p for q,p in a)%2==0
def flatten(blocks,order):return [(k,c) for i in order for k,c in blocks[i].items()]
def direct_objective(blocks,order,hf,den):
    diag={}
    for b in blocks:
        if flip(next(iter(b)))==0:
            for k,c in b.items():diag[k]=diag.get(k,0)+c
    energy0=unit_action(diag,hf)[1]/den
    amplitudes={}
    for at,i in enumerate(order):
        desti,fi=unit_action(blocks[i],hf)
        for j in order[at+1:]:
            destj,fj=unit_action(blocks[j],hf)
            dji,fij=unit_action(blocks[i],destj)
            dij,fji=unit_action(blocks[j],desti)
            require(dji==dij,'XOR destination mismatch')
            value=fj*fij-fi*fji
            amplitudes[dji]=amplitudes.get(dji,0)+value
    require(amplitudes.get(hf,0)==0,'Nonzero diagonal commutator expectation')
    terms=[]
    for dest,c in amplitudes.items():
        if c:
            gap=(unit_action(diag,dest)[1]/den-energy0)/2
            sinc=math.sin(gap)/gap if gap else 1.
            terms.append((sinc*c/(den*den))**2/4)
    return math.fsum(terms)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    manifest=json.loads((args.run/'manifest.json').read_text())
    require(sha(args.parent/'manifest.json')==manifest['parent_manifest_sha256'],'Parent manifest mismatch')
    for rel,h in manifest['source_sha256'].items():require(sha(args.run/rel)==h,'Frozen source mismatch')
    summaries=[];objective_max=0.;totalunits=totalschedules=totalbasischecks=0
    for case in manifest['cases']:
        pd=json.loads((args.parent/f'{case}.design.json').read_text());d=json.loads((args.run/f'{case}.design.json').read_text());candidate=json.loads((args.run/f'{case}.candidate.json').read_text())
        require(sha(args.parent/f'{case}.design.json')==candidate['parent_design_sha256'],'Candidate parent-design mismatch')
        for key in ('coefficients','denominator','tensor_sha256','affine_basis_indices','initial_sha256','target_sha256','electrons','n_qubits'):
            require(d[key]==pd[key],'Changed paired target/state metadata '+key)
        with np.load(args.parent/f'{case}.states.npz',allow_pickle=False) as a,np.load(args.run/f'{case}.states.npz',allow_pickle=False) as b:
            for key in ('initial','exact','exact_dense'):require(np.array_equal(a[key],b[key]),'Changed initial/exact states')
        bs=[parse(b) for b in d['blocks']['class']];old=[parse(b) for b in pd['blocks']['class']];total=parse(d['coefficients']);den=int(d['denominator']);hf=(1<<d['electrons'])-1
        seen={}
        for b in bs:
            require(bool(b),'Empty unit')
            for k,c in b.items():require(k not in seen,'Pauli appears in two units');seen[k]=c
        require(seen==total,'Partition is not exactly the retained Hamiltonian')
        xs=[flip(next(iter(b))) for b in bs];require(xs==d['class_flip_masks'],'Saved masks mismatch')
        nonzero=[x for x in xs if x];require(len(nonzero)==len(set(nonzero)),'Repeated nonzero mask')
        for b,x in zip(bs,xs):
            if x==0:require(len(b)==1 or 'unit_partition' not in d,'Refined diagonal unit not a singleton')
            for v in range(1<<d['n_qubits']):
                dest,value=unit_action(b,v)
                require(value==0 or v.bit_count()==dest.bit_count(),'Exact particle-number violation')
                totalbasischecks+=1
            words=list(b)
            require(all(commute_words(words[i],words[j]) for i in range(len(words)) for j in range(i)),'Noncommuting internal Pauli factors')
        require(flatten(bs,d['orders']['class_lex'])==flatten(old,pd['orders']['class_lex']),'Refined fixed order changes original lex product')
        stale=[]
        for name,o in d['orders'].items():
            if name.startswith('class_') and sorted(o)!=list(range(len(bs))):stale.append(name)
        for run in d['planned_runs']:
            name=run['ordering'];o=d['orders'][name]
            require(sorted(o)==list(range(len(bs))),'Scheduled order is not a permutation')
            seq=flatten(bs,o);require(len(seq)==len(total) and dict(seq)==total,'Changed or repeated scheduled coefficient')
            row=json.loads((args.run/f"{case}.{name}.r{run['steps']}.{run['role']}.json").read_text())
            require(row['T']==run['T']==1. and row['steps']==run['steps'],'Changed time/steps')
            require(row['pauli_rotations']==len(total)*run['steps'],'Wrong equal-work claim')
            J=direct_objective(bs,o,hf,den);delta=abs(J-row['surrogate_J']);objective_max=max(objective_max,delta)
            require(delta<1e-12*max(J,1e-20),'Independent surrogate does not agree')
            totalschedules+=1
        for event in candidate['accepted_swaps']:
            require(event['new_J']<event['old_J']-max(1e-30,1e-12*event['old_J']),'Non-improving accepted move')
        require(candidate['J_D']==min(x['J'] for x in candidate['candidate_start_results']),'Candidate not minimal of starts')
        totalunits+=len(bs)
        summaries.append({'case':case,'units':len(bs),'zero_mask_units':xs.count(0),'unique_paulis':len(total),
            'all_units_conserve_N_exactly_on_full_basis':True,'all_internal_paulis_commute_exactly':True,
            'lex_product_identical_to_parent':True,'stale_inactive_orders':stale})
        print('AUDITED',case,flush=True)
    report={'status':'passed','qhat_imports':False,'scope':'Exact target, active schedule and full computational-basis block conservation; independent objective recomputation. Circuit propagation is checked separately.',
            'run_manifest_sha256':sha(args.run/'manifest.json'),'script_sha256':sha(Path(__file__)),
            'active_schedules_checked':totalschedules,'units_checked':totalunits,'exact_basis_actions_checked':totalbasischecks,
            'max_surrogate_absolute_difference':objective_max,'cases':summaries}
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='cases'}))
if __name__=='__main__':main()
