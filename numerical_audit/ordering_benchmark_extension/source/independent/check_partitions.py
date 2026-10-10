#!/usr/bin/env python3
"""Exact integer target, internal commutation and number-conservation checks."""
import argparse,io,json,zipfile
from pathlib import Path
from core import block,physical_masks,number_commutator,require,sha

def main():
    p=argparse.ArgumentParser();p.add_argument('--qhat-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    cases=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]+[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')];checks=[]
    for case in cases:
        if '_s' in case:db=(a.qhat_root/'numerical_audit/followup_classes/results/comparison'/f'{case}.design.json').read_bytes()
        else:
            with zipfile.ZipFile(a.qhat_root/'numerical_audit/ordering_mechanism/heldout/heldout_validation.zip') as z:db=z.read(f'results/{case}.design.json')
        d=json.loads(db);n=d['n_qubits'];total=block(d['coefficients']);groups={k:d['blocks'][k] for k in ('parent','owned','class')};refined=[]
        for b in groups['class']:
            if physical_masks(next(iter(block(b))),n)[0]==0:refined.extend([[w] for w in b])
            else:refined.append(b)
        groups['refined']=refined;counts={}
        for name,bs in groups.items():
            summed={}
            for b in bs:
                require(not number_commutator(b),'Number-nonconserving '+case+'/'+name)
                factors=[]
                for w,v in block(b).items():
                    summed[w]=summed.get(w,0)+v;x,z,y=physical_masks(w,n);require(y%2==0,'Non-real word');factors.append((x,z))
                for i,(x,z) in enumerate(factors):
                    for y,t in factors[i+1:]:require(((x&t).bit_count()+(z&y).bit_count())%2==0,'Internal noncommutation')
            require({k:v for k,v in summed.items() if v}==total,'Partition changes target')
            counts[name]={'units':len(bs),'emitted_terms':sum(len(b) for b in bs),'exact_target_sum':True,'exact_number_commutator_zero':True,'pairwise_internal_commutation':True}
        checks.append({'case':case,'checks':counts,'input_design_sha256':__import__('hashlib').sha256(db).hexdigest()});print('PARTITIONS',case,flush=True)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps({'status':'passed','source_sha256':sha(__file__),'core_sha256':sha(Path(__file__).with_name('core.py')),'cases':len(checks),'checks':checks},indent=2)+'\n')
if __name__=='__main__':main()
