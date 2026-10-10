#!/usr/bin/env python3
"""Certify the 32 development schedules changed only inside commuting blocks.

The historical runner iterated its original dictionaries but serialized their
coefficients in sorted order. Reconstruct its original dictionaries from the
canonical components, reproduce its recorded flat hashes, and compare states.
"""
import json,sys
from pathlib import Path
from collections import defaultdict
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import *
import numpy as np
from qhat.numerical_audit.followup_classes import exact_blocks as eb
from qhat.analysis.fermionic_coefficients import _unit_transform

def anti(a,b):
 x,z,_=masks(a);xx,zz,_=masks(b);return ((x&zz)^(z&xx)).bit_count()%2

def original_blocks(d):
 parents=[]
 for meta,braw in zip(d['parent_metadata'],d['blocks']['parent']):
  comp=[tuple((int(q),int(a)) for q,a in k) for k in meta['components']];weight=int(meta['weight']);first=comp[0];partner,phase=eb.adjoint(first);components={first:weight}
  if partner!=first:components[partner]=phase*weight
  mapped=defaultdict(int)
  for component,wgt in components.items():
   for w,u in _unit_transform(component,'jw'):
    num,den=complex(u).real.as_integer_ratio();require(wgt%den==0,'Nonintegral coefficient reconstruction');mapped[w]+=wgt//den*num
  b={k:v for k,v in mapped.items() if k and v};require(b==block(braw),'Reconstructed parent differs from immutable input');parents.append(b)
 total=block(d['coefficients']);seen=set();owned=[]
 for parent in parents:
  b={k:total[k] for k in parent if k in total and k not in seen};seen.update(b)
  if b:owned.append(b)
 require([dict(b) for b in owned]==[block(b) for b in d['blocks']['owned']],'Ownership changed')
 return {'parent':parents,'owned':owned}

def main():
 oldrows=json.loads((BASE/'results.json').read_text());checks=[]
 for case in DEV_CASES:
  d,_=load_case(case);new=json.loads((HERE/'designs'/f'{case}.json').read_text());bs=original_blocks(d)
  with np.load(BASE/f'{case}.states.npz') as states:
   for name in ('parent_signed','parent_color','owned_signed','owned_color'):
    kind=name.split('_')[0];oldseq=[(w,v) for i in d['orders'][name] for w,v in bs[kind][i].items()];oldhash=object_sha(encode_schedule(oldseq));require(oldhash==d['schedule_sha256'][name],'Cannot reproduce historical actual flat hash')
    seq=decode_schedule(new['methods'][name]['schedule']);require(len(seq)==len(oldseq),'Length changed');pos=0;pairs=0
    for i in d['orders'][name]:
     b=bs[kind][i];piece=seq[pos:pos+len(b)];require(len(dict(piece))==len(piece) and dict(piece)==b,'Block contents or boundaries changed');words=list(b)
     for j,w in enumerate(words):
      for v in words[:j]:require(not anti(w,v),'Noncommuting internal rearrangement');pairs+=1
     pos+=len(b)
    require(pos==len(seq),'Block sizes incomplete')
    row=json.loads((HERE/'results'/f'{case}.{name}.b100.json').read_text());old=next(r for r in oldrows if r['case']==case and r['ordering']==name and r['steps']==row['steps']);key=f'{name}__r{old["steps"]}__{old["role"]}';before=states[key]
    with np.load(HERE/'results'/row['state_file']) as z:after=z['actual']
    delta=float(np.linalg.norm(before-after));Edelta=abs(old['one_minus_overlap']-row['one_minus_overlap']);Ndelta=abs(old['particle_number_leakage']-row['particle_number_leakage'])
    require(delta<1e-10 and Edelta<1e-12,'Historical reconstruction has material numerical discrepancy')
    checks.append({'case':case,'ordering':name,'historical_flat_sha256_reconstructed':oldhash,'new_schedule_sha256':new['methods'][name]['schedule_sha256'],'positions_reordered':sum(a!=b for a,b in zip(seq,oldseq)),'same_ordered_blocks_exact':True,'internal_pauli_pairs_checked_commuting':pairs,'steps':row['steps'],'actual_rotations':row['pauli_rotations'],'historical_state_key':key,'historical_state_npz_sha256':sha(BASE/f'{case}.states.npz'),'new_state_npz_sha256':row['state_sha256'],'state_norm_distance':delta,'one_minus_overlap_absolute_difference':Edelta,'particle_number_leakage_absolute_difference':Ndelta})
 report={'status':'pass','scope':'32 historical development flat sequences differ only by permutations inside exactly commuting parent/owned blocks; ideal product formulas are identical for every step count and time. Actual saved-state comparisons below use matched 100U budget and the historical row with the same step count.','interpretation':'The earlier 96-check compares schedules reconstructed from sorted serialized block data; it is not 96 literal historical-execution hash matches. There are 64 literal historical matches and 32 certified commuting-equivalent sequences.','checks':checks,'max_state_norm_distance':max(x['state_norm_distance'] for x in checks),'max_one_minus_overlap_absolute_difference':max(x['one_minus_overlap_absolute_difference'] for x in checks),'max_particle_number_leakage_absolute_difference':max(x['particle_number_leakage_absolute_difference'] for x in checks),'unchanged_frozen_campaign':True}
 write_atomic(HERE/'historical_internal_order_check.json',report);print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
