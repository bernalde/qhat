#!/usr/bin/env python3
"""Frozen, bounded original-parent / owned / full flip-mask experiment.

Standalone entry point; see protocol.md. No archived experiment is overwritten.
"""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-class-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-class-numba')
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse.linalg import expm_multiply
from scipy import sparse
from openfermion import FermionOperator
import qhat
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis.fermionic_coefficients import fermion_from_interaction,jordan_wigner_molecular
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics,sha,object_sha,canonical,write
from exact_blocks import (exact_parents,decompose,add_into,number_commutator,internally_commuting,
    graph_and_color,suffix_greedy,suffix_objective,symplectic)

VERSION='exact-parent-common-target-class-comparison-v1'
ROOT=Path(baseline.__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
FIXTURES=ROOT/'numerical_audit/fixtures'
NAMES=('parent_signed','parent_color','owned_signed','owned_color','class_lex','class_color','class_suffix_l1','jw_signed','jw_magnitude')
CUTOFF=1e-12


def records(block):
    return [[key,str(value)] for key,value in sorted(block.items())]


def direct_fermion_matrix(parents,den,states,n):
    """Independent occupation action; no OpenFermion mapping or Pauli code."""
    positions={int(b):i for i,b in enumerate(states)}
    rr=[];cc=[];vv=[]
    for p in parents:
        first=p['components'][0]
        # reconstruct coefficient signs from the exact normal adjoint.
        from exact_blocks import adjoint
        partner,phase=adjoint(first)
        components=[(first,p['weight'])]
        if partner!=first: components.append((partner,phase*p['weight']))
        for key,c in components:
            for col,initial in enumerate(states):
                b=int(initial); amplitude=c/den
                for q,creation in reversed(key):
                    bit=1 << (n-1-q)
                    occupied=bool(b&bit)
                    if occupied==bool(creation): amplitude=0.;break
                    if (b>>(n-q)).bit_count()%2: amplitude=-amplitude
                    b^=bit
                if amplitude:
                    if b not in positions: raise RuntimeError('Chosen closure not invariant under a parent.')
                    rr.append(positions[b]);cc.append(col);vv.append(amplitude)
    result=sparse.coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr()
    scalar=sum(p['scalar'] for p in parents)/den
    return result-scalar*sparse.eye(len(states),format='csr')


def make_case(case,fixture,output):
    tensor=FIXTURES/case/fixture['tensor']
    if sha(tensor)!=fixture['files'][fixture['tensor']]:raise RuntimeError('Input tensor changed.')
    interaction,n=baseline.load_interaction_operator(tensor)
    raw=fermion_from_interaction(interaction)
    parents,total,den,policy=exact_parents(raw,CUTOFF)
    pb,ob,cb,masks=decompose(parents,total)
    blocks={'parent':pb,'owned':ob,'class':cb}
    for name,bs in blocks.items():
        reconstructed={}
        for b in bs:
            add_into(reconstructed,b)
            if not internally_commuting(b) or number_commutator(b):
                raise RuntimeError(f'{name} block violates exact commutation/conservation.')
        if reconstructed!=total:raise RuntimeError('Exact target reconstruction mismatch.')
    graphs={}; orders={'parent_signed':list(range(len(pb))),'owned_signed':list(range(len(ob))),
                      'class_lex':list(range(len(cb)))}
    class_comms=None
    for kind,bs in blocks.items():
        order,graph,cs=graph_and_color(bs,keep_commutators=kind=='class')
        orders[kind+'_color']=order;graphs[kind]=graph
        if kind=='class':class_comms=cs
    greedy,greedy_info=suffix_greedy(len(cb),class_comms)
    orders['class_suffix_l1']=greedy
    objectives={name:suffix_objective(order,class_comms)/(den*den)
                for name,order in orders.items() if name.startswith('class_')}
    if objectives['class_suffix_l1'] > int(greedy_info['greedy_objective_integer'])/(den*den)+1e-15:
        raise RuntimeError('Local improvement increased objective.')
    coeff={k:v/den for k,v in total.items()}
    jw_orders={'jw_signed':sorted(coeff,key=lambda k:(coeff[k],baseline.dense_pauli_string(k,n))),
               'jw_magnitude':sorted(coeff,key=lambda k:(-abs(coeff[k]),baseline.dense_pauli_string(k,n)))}
    # The parent union retains all intermediate flips, including canceled words.
    keys=sorted(set(total).union(*(set(b) for b in pb)))
    rawunit=evolution.compile_ordered_terms(keys,{k:1. for k in keys},n,0.)
    full_initial=evolution.build_hartree_fock_state(n,fixture['electrons'])
    hf=int(np.argmax(abs(full_initial)))
    states,reduced,_=ref.affine_space(rawunit,hf)
    unit=dict(zip(keys,reduced))
    def terms(block):
        return [(v/den,unit[k][1],unit[k][2],unit[k][3]) if unit[k][0]==1 else
                (v/den*unit[k][0],unit[k][1],unit[k][2],unit[k][3])
                for k,v in block.items()]
    compiled={}
    flat={}
    for name,order in orders.items():
        kind=name.split('_')[0]
        flat[name]=[(k,v) for i in order for k,v in blocks[kind][i].items()]
        compiled[name]=[term for i in order for term in terms(blocks[kind][i])]
    for name,order in jw_orders.items():
        flat[name]=[(k,total[k]) for k in order]
        compiled[name]=[terms({k:total[k]})[0] for k in order]
    H,_=ref.assemble_hamiltonian(terms(total),len(states))
    independent=direct_fermion_matrix(parents,den,states,n)
    matrix_error=float(sparse.linalg.norm(H-independent))
    if matrix_error>1e-10:raise RuntimeError('Independent occupation/JW matrices disagree.')
    initial=np.zeros(len(states),complex);initial[0]=1
    nm=np.array([int(x).bit_count()==fixture['electrons'] for x in states])
    physical=set(map(int,evolution.spin_sector_basis_indices(n,fixture['electrons'])))
    sm=np.array([int(x) in physical for x in states])
    rawjw=jordan_wigner_molecular(raw,atol=0).terms
    prior=jordan_wigner_molecular(raw,atol=1e-12).terms
    def distance(other):return math.fsum(abs(coeff.get(k,0)-other.get(k,0)) for k in set(coeff)|set(other) if k)
    # Rounding occurs only at the simulation interface; quantify schedule sum drift.
    rounding={}
    for name,schedule in flat.items():
        accum={}
        for k,v in schedule:accum.setdefault(k,[]).append(v/den)
        rounding[name]=math.fsum(abs(math.fsum(vs)-coeff.get(k,0)) for k,vs in accum.items())
    rotations={name:len(schedule) for name,schedule in flat.items()}
    matched=max(1,100*len(total)//rotations['parent_signed'])
    settings=[(100,'equal_steps')]
    if case=='HCN_s1.00':settings += [(50,'convergence'),(200,'convergence')]
    schedules=[(name,r,role) for r,role in settings for name in NAMES]
    schedules += [(name,matched,'matched_rotation_budget') for name in ('parent_signed','parent_color')]
    design={'case':case,'n_qubits':n,'electrons':fixture['electrons'],'tensor_sha256':sha(tensor),
        'version':VERSION,'denominator':str(den),'coefficient_policy':policy,
        'coefficients':records(total),'parent_metadata':[{'components':p['components'],'weight':str(p['weight']),
            'four_distinct':p['four_distinct'],'scalar':str(p['scalar'])} for p in parents],
        'blocks':{name:[records(b) for b in bs] for name,bs in blocks.items()},
        'class_flip_masks':masks,'orders':orders,'jw_orders':jw_orders,'graphs':graphs,
        'class_suffix_l1_upper_bounds':objectives,'suffix_candidate_construction':greedy_info,
        'pauli_rotations_per_step':rotations,'parent_duplication_ratio':len(flat['parent_signed'])/len(total),
        'planned_runs':[{'ordering':name,'steps':r,'role':role,'T':1.} for name,r,role in schedules],
        'matched_budget':{'consolidated_rotations':100*len(total),'parent_steps':matched,
            'parent_rotations':matched*rotations['parent_signed']},
        'affine_basis_indices':list(map(int,states)),'dimension':len(states),
        'target_sha256':object_sha(records(total)),'initial_sha256':hashlib.sha256(initial.tobytes()).hexdigest(),
        'state':'Hartree-Fock, q0 most significant bit; no sector projection',
        'comparison_to_raw_nonidentity_pauli_l1':distance(rawjw),
        'comparison_to_prior_finalcutoff_nonidentity_pauli_l1':distance(prior),
        'float_schedule_sum_target_l1':rounding,'independent_occupation_matrix_frobenius':matrix_error,
        'exact_integer_number_commutator_zero':True,'internal_pauli_commutation_exact':True,
        'schedule_sha256':{name:object_sha([[k,str(v)] for k,v in schedule]) for name,schedule in flat.items()},
        'product_convention':'chronological entries, later factors multiply on the left'}
    # Complete design, all schedules and graphs are saved before any errors are evaluated.
    write(output/f'{case}.design.json',design)
    target=output/'inputs'/case/tensor.name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(tensor,target)
    return H,initial,compiled,nm,sm,design


def summarize(rows,output):
    primary=[r for r in rows if r['role']=='equal_steps']
    bycase={}
    for row in primary:bycase.setdefault(row['case'],{})[row['ordering']]=row
    comparisons={}
    for name in NAMES:
        if name=='owned_signed':continue
        outcomes=[];ratios=[]
        for case,group in bycase.items():
            a=group[name]['one_minus_overlap'];b=group['owned_signed']['one_minus_overlap']
            outcomes.append('tie' if abs(a-b)<=1e-12 else ('win' if a<b else 'loss'))
            if b>1e-12:ratios.append(a/b)
        comparisons[name]={'candidate_vs_owned_signed':dict(Counter(outcomes)),
                           'median_error_ratio_over_owned_signed':float(np.median(ratios)) if ratios else None}
    summary={'rows':len(rows),'primary_rows':len(primary),'cases':len(bycase),
        'comparisons':comparisons,'max_completed_schedule_number_leakage':max(r['particle_number_leakage'] for r in rows if not r['ordering'].startswith('jw_')),
        'max_jw_number_leakage':max(r['particle_number_leakage'] for r in rows if r['ordering'].startswith('jw_')),
        'interpretation':'HF, fixed reconstructed canonical target, descriptive eight-case comparison. No general accuracy or graph-color advantage. Pauli rotations are not compiled gates.'}
    write(output/'summary.json',summary)
    lines=['# Conserving parent / owned / complete-class comparison','',
           'All arms share the explicitly saved canonical-parent target. Whole-parent cutoff and exact Hermitian pairing make this a separate campaign; target changes relative to earlier campaigns are recorded. These eight geometries represent four related molecular families, not eight independent populations.','',
           '| Case | Parent signed | Parent color | Owned signed | Owned color | Class lex | Class color | Class suffix l1 | JW signed | JW magnitude |',
           '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for case,group in sorted(bycase.items()):lines.append('| '+case+' | '+' | '.join(f'{group[name]["one_minus_overlap"]:.6g}' for name in NAMES)+' |')
    lines+=['','Values are stable normalized full-state E=1-|overlap| at T=1,r=100. Coloring changes the order, not the exact factorization of a contiguous commuting color. The sign-retaining suffix l1 objective is an upper-bound objective, not an accuracy certificate. Equal-step parent formulas use more rotations; saved matched_rotation_budget rows use floor(100 U/R) steps.','',
            f'Maximum completed conserving-schedule number leakage: {summary["max_completed_schedule_number_leakage"]:.6g}. All ideal integer blocks have exactly zero number commutator; that computational check is not a Lean proof.','']
    (output/'REPORT.md').write_text('\n'.join(lines))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cases',nargs='+')
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    args.output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((FIXTURES/'manifest.json').read_text())
    cases=args.cases or sorted(k for k,v in manifest.items() if v['qubits']==10)
    if any(k not in manifest for k in cases):raise ValueError('Unknown fixture.')
    sources=[Path(__file__),HERE/'exact_blocks.py',HERE/'test_exact_blocks.py',HERE/'protocol.md']
    sources += [ROOT/rel for rel in ('analysis/fermionic_coefficients.py','analysis/fermionic_ownership.py',
        'analysis/state_metrics.py','analysis/benchmark_b2_active_spaces_matrix_free.py',
        'analysis/benchmark_b2_signed_coefficient_baseline.py','analysis/benchmark_L_sweep_trotter.py',
        'analysis/benchmark_fermionic_structure_ablation.py','analysis/benchmark_b2_coloring_robustness.py',
        'numerical_audit/reference_tools.py','numerical_audit/run_frozen_audit.py')]
    hashes={str(p):sha(p) for p in sources}
    for p in sources:
        rel=p.relative_to(ROOT) if p.is_relative_to(ROOT) else Path('numerical_audit/followup_classes')/p.name
        target=args.output/'source_snapshot'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    declared={'version':VERSION,'cases':cases,'methods':NAMES,'time':1,'primary_steps':100,
        'convergence_case':'HCN_s1.00','convergence_steps':[50,200],'tie_absolute':1e-12,
        'source_sha256':hashes,'python':sys.version,'packages':{p:importlib.metadata.version(p) for p in ('numpy','scipy','openfermion','numba')},
        'git_head_before_run':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'fixture_manifest_sha256':sha(FIXTURES/'manifest.json'),'prepare_only':args.prepare_only,
        'expected_rows':11*len(cases)+(18 if 'HCN_s1.00' in cases else 0)}
    write(args.output/'manifest.json',declared)
    rows=[];start=time.monotonic()
    try:
        for case in cases:
            began=time.monotonic()
            H,initial,compiled,nm,sm,design=make_case(case,manifest[case],args.output)
            print('PREPARED',case,design['dimension'],design['pauli_rotations_per_step'],time.monotonic()-began,flush=True)
            if args.prepare_only:continue
            exact=expm_multiply(-1j*H,initial,traceA=-1j*H.diagonal().sum())
            vals,vecs=eigh(H.toarray())
            independent=vecs@(np.exp(-1j*vals)*(vecs.conj().T@initial))
            delta=float(np.linalg.norm(exact-independent))
            if delta>1e-10:raise RuntimeError('Independent exact reference disagrees.')
            saved={'initial':initial,'exact':exact,'exact_dense':independent}
            write(args.output/f'{case}.reference.json',{'sparse_dense_norm_distance':delta,'norm':float(np.linalg.norm(exact))})
            for run in design['planned_runs']:
                name,r,role=run['ordering'],run['steps'],run['role']
                began=time.monotonic();approx=evolve(initial,compiled[name],1.,r)
                row={'case':case,'ordering':name,'T':1.,'steps':r,'role':role,
                     'target_sha256':design['target_sha256'],'schedule_sha256':design['schedule_sha256'][name],
                     'pauli_rotations':r*len(compiled[name]),'seconds':time.monotonic()-began,
                     **physical_metrics(exact,approx,nm,sm)}
                saved[f'{name}__r{r}__{role}']=approx
                rows.append(row);write(args.output/f'{case}.{name}.r{r}.{role}.json',row)
                print('ROW',len(rows),case,name,r,role,row['one_minus_overlap'],flush=True)
            np.savez_compressed(args.output/f'{case}.states.npz',**saved)
        if not args.prepare_only:
            if len(rows)!=declared['expected_rows']:raise RuntimeError('Incomplete expected rows.')
            write(args.output/'results.json',rows);summarize(rows,args.output)
        for p,digest in hashes.items():
            if sha(p)!=digest:raise RuntimeError('Source changed during campaign.')
        write(args.output/'completion.json',{'status':'prepared' if args.prepare_only else 'complete','rows':len(rows),'seconds':time.monotonic()-start})
    except Exception as exc:
        write(args.output/'failure.json',{'exception':repr(exc),'rows':len(rows)});raise

if __name__=='__main__':main()
