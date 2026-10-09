#!/usr/bin/env python3
"""Independent saved-state, exact coefficient/tensor, and gauge invariance audit.

No QHAT imports. This rechecks recorded trajectories; it is not a second
independent propagation of all 132 trajectories. The producer separately
records direct full-space NumPy propagation for all 18 ten-qubit schedules.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
def require(x,msg):
    if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def objsha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def arrsha(a):
    a=np.asarray(a)
    return objsha({'shape':a.shape,'dtype':a.dtype.str,'raw_sha256':hashlib.sha256(a.tobytes()).hexdigest()})
def parseco(d):return {tuple((int(p),a) for p,a in k):complex(float.fromhex(r),float.fromhex(i)) for k,r,i in d['identity_free_coefficients']}
def mask(k):return sum(1<<p for p,a in k if a in ('X','Y'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    root=args.results;manifest=json.loads((root/'artifact_manifest.json').read_text())
    for rel,h in manifest.items():require(sha(root/rel)==h,'Artifact hash mismatch: '+rel)
    protocol=json.loads((root/'protocol.json').read_text());allrows=[];summaries=[]
    maxmetric=maxfixed=maxleak=0.;negative=0
    for case in protocol['cases']:
        rows=[];designs=[];archives=[]
        for g in range(6):
            gd=root/case/f'gauge{g}';d=json.loads((gd/'design.json').read_text());r=json.loads((gd/'results.json').read_text())
            require(sha(gd/'states.npz')==r['state_file_sha256'],'Changed states')
            require(sha(gd/'design.json')==r['design_sha256'],'Changed design')
            with np.load(gd/'states.npz',allow_pickle=False) as f:a={k:f[k] for k in f.files}
            for k,v in a.items():require(arrsha(v)==r['array_sha256'][k],'Changed array')
            co=parseco(d);n=d['n_qubits'];spin=d['spin_signs'];states=a['affine_basis_indices'];coord={int(s):j for j,s in enumerate(states)}
            diag=np.ones(len(states),dtype=np.int8)
            for p,s in enumerate(spin):
                if s<0:diag*=1-2*((states>>(n-1-p))&1)
            require(np.array_equal(diag,a['gauge_diagonal']),'Wrong gauge diagonal')
            with np.load(root/case/'original.tensors.npz',allow_pickle=False) as base, np.load(gd/'transformed.tensors.npz',allow_pickle=False) as got:
                sg=np.array(spin)
                require(np.array_equal(got['one_body'],base['one_body']*sg[:,None]*sg[None,:]),'Wrong 1-body transform')
                require(np.array_equal(got['two_body'],base['two_body']*sg[:,None,None,None]*sg[None,:,None,None]*sg[None,None,:,None]*sg[None,None,None,:]),'Wrong 2-body transform')
                require(np.array_equal(got['constant'],base['constant']),'Scalar transform error')
            if g==0:baseco=co;basea=a;baseclasses={}
            for k,c in co.items():
                sign=math.prod(spin[p] for p,s in k if s in ('X','Y'))
                require(c==sign*baseco[k],'Exact coefficient covariance failed')
                flip=sum(1<<(n-1-p) for p,s in k if s in ('X','Y'))
                require(all(int(s)^flip in coord for s in states),'Affine closure projects a Pauli destination')
            classes={}
            for k,c in co.items():
                classes.setdefault(mask(k),[]).append((k,Fraction(abs(c.real))))
            classes={k:sorted(v) for k,v in classes.items()}
            if g==0:baseclasses=classes
            require(classes==baseclasses,'Complete-mask coefficient magnitudes changed')
            require(np.array_equal(a['initial'],diag*basea['initial']),'Incorrect initial gauge state')
            require(np.array_equal(a['exact'],diag*basea['exact']),'Incorrect exact gauge state')
            expectedgaugeorder={m:tuple(v) for m,v in designs[0]['resorted_orders'].items()} if g else {m:tuple(v) for m,v in d['resorted_orders'].items()}
            require(tuple(d['resorted_orders']['jw_magnitude'])==expectedgaugeorder['jw_magnitude'],'Magnitude order changed')
            negative+=int(d['wrong_sign_negative_control_rejected'] is True)
            nm=np.array([int(s).bit_count()==d['n_electrons'] for s in states])
            for row in r['rows']:
                state=a[f"{row['ordering']}_{row['mode']}"];e=a['exact']/np.linalg.norm(a['exact']);v=state/np.linalg.norm(state)
                ov=np.vdot(e,v);aligned=v if abs(ov)==0 else v*ov.conjugate()/abs(ov)
                E=float(np.vdot(e-aligned,e-aligned).real/2);delta=abs(E-row['one_minus_overlap']);maxmetric=max(maxmetric,delta)
                require(delta<1e-18,'Saved overlap metric mismatch')
                leak=float(np.vdot(v[~nm],v[~nm]).real);maxleak=max(maxleak,abs(leak-row['particle_number_leakage']))
                require(abs(leak-row['particle_number_leakage'])<1e-18,'Leakage mismatch')
                if row['mode']=='fixed_control' or row['ordering']=='jw_magnitude' or g==0:
                    dist=float(np.linalg.norm(state-diag*basea[f"{row['ordering']}_resorted"]));maxfixed=max(maxfixed,dist)
                    require(dist<1e-10,'Fixed-order gauge covariance failed')
                rows.append(row);allrows.append(row)
            designs.append(d);archives.append(a)
        summary={'case':case,'n_qubits':n,'closed_dimension':len(states),'mask_classes':len(classes),
                 'all_complete_mask_coefficient_magnitudes_exactly_invariant':True,
                 'parent_order_distinct':len({tuple(d['parent_order']) for d in designs}),
                 'owner_changed_paulis':[d['owner_changed_paulis'] for d in designs],
                 'block_partition_changed':[d['block_partition_changed'] for d in designs],'methods':{}}
        for m in protocol['methods']:
            selected=[r for r in rows if r['mode']=='resorted' and r['ordering']==m];values=[r['one_minus_overlap'] for r in selected]
            summary['methods'][m]={'E_by_gauge':values,'E_min':min(values),'E_max':max(values),'E_max_over_min':max(values)/min(values),
                'distinct_word_orders':len({r['schedule_sha256'] for r in selected}),'changed_order_positions':[r['changed_order_positions'] for r in selected]}
        owned=summary['methods']['fermionic_signed']['E_by_gauge']
        for other in ('jw_magnitude','jw_signed'):
            otherE=summary['methods'][other]['E_by_gauge'];delta=[a-b for a,b in zip(owned,otherE)]
            summary['owned_vs_'+other]={'wins':sum(d<-1e-12 for d in delta),'ties':sum(abs(d)<=1e-12 for d in delta),'losses':sum(d>1e-12 for d in delta),
                'E_owned_over_other_by_gauge':[a/b for a,b in zip(owned,otherE)]}
        summaries.append(summary)
    require(len(allrows)==132,'Wrong row count')
    result={'status':'passed','rows_checked':len(allrows),'artifact_files_checked':len(manifest),'qhat_imported':False,
            'independent_propagation_performed_by_this_script':False,'exact_complete_class_magnitude_controls':24,
            'wrong_sign_negative_controls_rejected':negative,'max_metric_difference':maxmetric,'max_fixed_state_distance':maxfixed,
            'max_leakage_difference':maxleak,'cases':summaries,'producer_checks':json.loads((root/'results.json').read_text())['checks']}
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','producer_checks')}))
if __name__=='__main__':main()
