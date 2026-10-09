#!/usr/bin/env python3
"""Create clean derived views without changing frozen experiment designs."""
import argparse,hashlib,json
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def word(w):return tuple((int(q),str(p)) for q,p in w)
def block(raw):return {word(w):int(v) for w,v in raw}

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((a.run/'manifest.json').read_text());records={}
 for case in manifest['cases']:
  source=a.run/f'{case}.design.json';d=json.loads(source.read_text());names=[r['ordering'] for r in d['planned_runs']];names=list(dict.fromkeys(names));bs=[block(b) for b in d['blocks']['class']];total=block(d['coefficients']);rotations={};orders={};schedule_hashes={}
  for name in names:
   if not name.startswith('class_'):raise ValueError('This cleaner is scoped to class candidate runs.')
   order=d['orders'][name]
   if sorted(order)!=list(range(len(bs))):raise ValueError(f'Nonpermutation {case}/{name}')
   merged={};emitted=[]
   for i in order:
    for w,v in bs[i].items():
     if w in merged:raise ValueError('Duplicate Pauli')
     merged[w]=v;emitted.append([w,str(v)])
   if merged!=total:raise ValueError('Changed target')
   rotations[name]=len(emitted);orders[name]=order;schedule_hashes[name]=hashlib.sha256(json.dumps(emitted,separators=(',',':')).encode()).hexdigest()
  view={k:d[k] for k in ('case','n_qubits','electrons','tensor_sha256','denominator','coefficient_policy','coefficients','class_flip_masks','planned_runs','affine_basis_indices','dimension','target_sha256','initial_sha256','state','product_convention') if k in d}
  view.update(schema='active-class-design-v1',source_design_sha256=sha(source),source_design=str(source.name),blocks={'class':d['blocks']['class']},orders=orders,pauli_rotations_per_step=rotations,schedule_sha256=schedule_hashes,unit_partition=d.get('unit_partition','complete flip classes'))
  out=a.output/f'{case}.json';out.write_text(json.dumps(view,indent=2)+'\n');records[out.name]={'sha256':sha(out),'source_sha256':sha(source)}
 (a.output/'manifest.json').write_text(json.dumps({'scope':'Clean active views only; raw frozen designs preserved. Each active order is a permutation and emits each target Pauli exactly once. Inactive inherited graph/suffix/order/hash metadata omitted.','files':records},indent=2)+'\n')
if __name__=='__main__':main()
