#!/usr/bin/env python3
"""Exact check that complete flip classes destroy ordering information in suffix l1."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-class-mpl')
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
from exact_blocks import commutator,symplectic,suffix_objective


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise ValueError('Output exists.')
    cases=[]
    for path in sorted(a.run.glob('*.design.json')):
        d=json.loads(path.read_text());blocks=[{tuple((int(q),label) for q,label in k):int(v) for k,v in b} for b in d['blocks']['class']]
        masks=d['class_flip_masks'];comms={};seen=[set() for _ in blocks]
        for i,j in combinations(range(len(blocks)),2):
            c=commutator(blocks[i],blocks[j]);comms[i,j]=c
            if any(symplectic(k)[0] != masks[i]^masks[j] for k in c):raise RuntimeError('Wrong commutator flip mask.')
            if set(c)&seen[i] or set(c)&seen[j]:raise RuntimeError('Different partner masks produced shared Pauli words.')
            seen[i].update(c);seen[j].update(c)
        pairwise=sum(sum(abs(v) for v in c.values()) for c in comms.values())
        orders={name:suffix_objective(order,comms) for name,order in d['orders'].items() if name.startswith('class_')}
        if any(value!=pairwise for value in orders.values()):raise RuntimeError('Suffix l1 differs from order-independent pair sum.')
        cases.append({'case':d['case'],'order_independent_pairwise_l1_integer':str(pairwise),
            'suffix_l1_integers':{k:str(v) for k,v in orders.items()},'all_partner_commutator_supports_disjoint':True,
            'design_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    result={'status':'passed','cases':cases,'conclusion':'For distinct complete flip masks, coefficient-l1 suffix bound equals the unordered pairwise sum for every order. Existing test schedules verified exactly; general reason is injectivity of y -> x xor y.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS exact order-independence:',len(cases),'cases')

if __name__=='__main__':main()
