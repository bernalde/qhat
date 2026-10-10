"""Shared frozen molecular input and propagation helpers for the extension.

All inputs come from immutable previously reviewed artifacts. QHAT can be
relocated by setting QHAT_ROOT; no molecular tensors are regenerated.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import hashlib, io, json, sys, zipfile
ROOT=Path(__file__).resolve().parent
QHAT=Path(os.environ.get('QHAT_ROOT',str(ROOT.parent/'qhat'))).resolve()
sys.path.insert(0,str(QHAT.parent))
DEV_CASES=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]
PROS_CASES=[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
CASES=DEV_CASES+PROS_CASES
BUDGETS=(1,2,5,10,20,50,100)
BASE=QHAT/'numerical_audit/followup_classes/results/comparison'
MECH=QHAT/'numerical_audit/ordering_mechanism'
HELDOUT=MECH/'heldout/heldout_validation.zip'

def require(value,message):
    if not value: raise RuntimeError(message)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def object_sha(value): return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write_atomic(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def block(raw): return {tuple((int(q),str(p)) for q,p in w):int(v) for w,v in raw}
def word(raw):return tuple((int(q),str(p)) for q,p in raw)
def encode_schedule(schedule):return [[w,str(v)] for w,v in schedule]
def decode_schedule(schedule):return [(word(w),int(v)) for w,v in schedule]
def masks(word):
    x=z=y=0
    for q,p in word:
        if p in ('X','Y'):x|=1<<q
        if p in ('Y','Z'):z|=1<<q
        y+=p=='Y'
    return x,z,y

def load_case(case):
    import numpy as np
    require(case in CASES,'Unknown case '+case)
    if case in DEV_CASES:
        source=BASE/f'{case}.design.json';state=BASE/f'{case}.states.npz'
        d=json.loads(source.read_text());provenance={'design_sha256':sha(source),'states_sha256':sha(state),'design_source':str(source.relative_to(QHAT)),'state_source':str(state.relative_to(QHAT))}
        with np.load(state) as z: arrays={k:z[k].copy() for k in ('initial','exact')}
    else:
        with zipfile.ZipFile(HELDOUT) as z:
            source=f'results/{case}.design.json';state=f'results/{case}.states.npz';db=z.read(source);sb=z.read(state)
        d=json.loads(db);provenance={'design_sha256':hashlib.sha256(db).hexdigest(),'states_sha256':hashlib.sha256(sb).hexdigest(),'archive_sha256':sha(HELDOUT),'design_source':source,'state_source':state}
        with np.load(io.BytesIO(sb)) as z:arrays={k:z[k].copy() for k in ('initial','exact')}
    d['extension_input_provenance']=provenance
    return d,arrays

def normalized_case(case):
    d,states=load_case(case);bs={k:[block(b) for b in vv] for k,vv in d['blocks'].items()};total=block(d['coefficients'])
    refined=[];xs=[]
    for b,x in zip(bs['class'],d['class_flip_masks']):
        for new in ([{w:v} for w,v in sorted(b.items())] if x==0 else [b]):refined.append(new);xs.append(x)
    bs['refined']=refined
    if case in DEV_CASES:
        p=MECH/'diagonal_refinement/results'/f'{case}.design.json';candidate=json.loads(p.read_text())
        candidate_order=candidate['orders']['class_refined_response'];jw={n:[word(w) for w in vv] for n,vv in d['jw_orders'].items()}
        require([block(b) for b in candidate['blocks']['class']]==refined,'Refined archived partition differs')
    else:
        candidate_order=d['candidate_order'];jw={n:[word(w) for w,_ in d['flat_schedules'][n]] for n in ('jw_signed','jw_magnitude')}
        require(bs['candidate']==refined,'Prospective refined partition differs')
    return {'design':d,'states':states,'blocks':bs,'total':total,'denominator':int(d['denominator']),'refined_flip_masks':xs,'candidate_order':candidate_order,'jw_orders':jw}

def compile_case(design,states):
    """Compile the full recorded affine closure, retaining cancelled-parent words."""
    import numpy as np
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    from qhat.numerical_audit import reference_tools as ref
    total=block(design['coefficients']);parents=[block(b) for b in design['blocks']['parent']]
    keys=sorted(set(total).union(*parents));n=design['n_qubits'];den=int(design['denominator'])
    raw=ev.compile_ordered_terms(keys,{k:1. for k in keys},n,0.)
    basis,reduced,_=ref.affine_space(raw,design['affine_basis_indices'][0])
    require(list(map(int,basis))==design['affine_basis_indices'],'Recorded affine closure changed')
    unit=dict(zip(keys,reduced));nr=np.array([int(v).bit_count()==design['electrons'] for v in basis])
    sector=set(map(int,ev.spin_sector_basis_indices(n,design['electrons'])));sr=np.array([int(v) in sector for v in basis])
    def term(w,v):return (v/den*unit[w][0],unit[w][1],unit[w][2],unit[w][3])
    return {'basis':basis,'term':term,'number_mask':nr,'spin_mask':sr,'initial':states['initial'],'exact':states['exact']}
