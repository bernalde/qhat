#!/usr/bin/env python3
"""No-evolution counts on exact binary tensors, with threshold diagnostics."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/qhat-class-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/qhat-class-numba')
import argparse
from collections import Counter,defaultdict
from itertools import permutations
import json
from pathlib import Path
import resource
import shutil
import time
from openfermion import FermionOperator,hermitian_conjugated
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis.fermionic_coefficients import fermion_from_interaction
from qhat.numerical_audit.run_frozen_audit import sha,write
from exact_blocks import exact_parents,add_into,symplectic

ROOT=Path(baseline.__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
THRESHOLDS=(0.,1e-14,1e-12,1e-10)


def exact_examples():
    groups=defaultdict(list);supports=set()
    for perm in permutations(range(4)):
        a,b,c,d=perm
        m=FermionOperator(((a,1),(b,1),(c,0),(d,0)),1)
        parents,total,den,_=exact_parents(m+hermitian_conjugated(m),0)
        if len(parents)!=1 or len(total)!=8:raise RuntimeError('24-term parent check failed.')
        groups[str(parents[0]['components'])].append(perm)
        supports.add(tuple(sorted(total)))
    if sorted(map(len,groups.values()))!=[8,8,8] or len(supports)!=1:
        raise RuntimeError('24 terms are not three classes with common support.')
    physical=[]
    for weights in ((1,3,2),(0,2,2),(1,1,0),(3,5,7)):
        raw=FermionOperator()
        for perm,w in zip(((0,1,2,3),(0,2,1,3),(0,3,1,2)),weights):
            a,b,c,d=perm;m=FermionOperator(((a,1),(b,1),(c,0),(d,0)),w)
            raw+=m+hermitian_conjugated(m)
        parents,total,den,_=exact_parents(raw,0)
        R=sum(len(p['coefficients']) for p in parents);U=len(total)
        physical.append({'triad_weights':weights,'relation_a_minus_b_plus_c':weights[0]-weights[1]+weights[2],
                         'parents':len(parents),'raw_occurrences':R,'consolidated_words':U,'ratio':R/U})
        if weights[0]-weights[1]+weights[2]==0 and R!=4*U:raise RuntimeError('Physical quartet ratio differs from4.')
    return {'permutation_count':24,'hermitian_class_count':len(groups),'class_sizes':list(map(len,groups.values())),
            'common_support_count':8,'groups':dict(groups),'quartet_examples':physical,
            'arithmetic':'exact arbitrary-precision integers; no simulation or tolerance'}


def count_one(raw,threshold):
    parents,total,den,policy=exact_parents(raw,threshold)
    four=[p for p in parents if p['four_distinct']]
    rest=[p for p in parents if not p['four_distinct']]
    four_total={};rest_total={}
    for p in four:add_into(four_total,p['coefficients'])
    for p in rest:add_into(rest_total,p['coefficients'])
    if set(four_total)&set(rest_total):raise RuntimeError('Four-flip and remaining sectors overlap.')
    if any(len(p['coefficients'])!=8 for p in four):raise RuntimeError('Nonzero real quartet parent not8words.')
    R4=sum(len(p['coefficients']) for p in four);Rrest=sum(len(p['coefficients']) for p in rest)
    U4=len(four_total);Urest=len(rest_total)
    U4t=sum(abs(v/den)>threshold for v in four_total.values())
    Urestt=sum(abs(v/den)>threshold for v in rest_total.values())
    quartets=defaultdict(list)
    for p in four:quartets[tuple(sorted(q for q,_ in p['components'][0]))].append(p)
    histogram=Counter();relations=[]
    for support,ps in quartets.items():
        weights=[0,0,0];merged={}
        for p in ps:
            creators={q for q,a in p['components'][0] if a}
            if support[0] not in creators:creators=set(support)-creators
            second=next(q for q in creators if q!=support[0])
            index=support.index(second)-1
            weights[index]=p['weight']
            add_into(merged,p['coefficients'])
        relations.append((weights[0]-weights[1]+weights[2])/den)
        histogram[str((len(ps),len(merged),sum(abs(v/den)>threshold for v in merged.values())))] += 1
    ratio=lambda a,b:a/b if b else None
    return {'parent_weight_cutoff':threshold,'global_pauli_diagnostic_cutoff':threshold,
        'parent_count':len(parents),'four_distinct_parent_count':len(four),'quartet_count':len(quartets),
        'raw_occurrences':R4+Rrest,'raw_four_occurrences':R4,'raw_remaining_occurrences':Rrest,
        'exact_binary_consolidated_words':U4+Urest,'exact_binary_four_words':U4,'exact_binary_remaining_words':Urest,
        'thresholded_consolidated_words':U4t+Urestt,'thresholded_four_words':U4t,'thresholded_remaining_words':Urestt,
        'raw_over_exact_binary_total':ratio(R4+Rrest,U4+Urest),
        'raw_over_exact_binary_four':ratio(R4,U4),
        'raw_over_thresholded_total':ratio(R4+Rrest,U4t+Urestt),
        'raw_over_thresholded_four':ratio(R4,U4t),
        'quartet_grouped_over_exact_binary_total':ratio(U4+Rrest,U4+Urest),
        'globally_consolidated_duplication_ratio':1 if total else None,
        'quartet_histogram_m_exactU_thresholdU':dict(sorted(histogram.items())),
        'physical_relation_nonzero_exact_count':sum(v!=0 for v in relations),
        'physical_relation_max_absolute':max(map(abs,relations),default=0),
        'physical_relation_above_diagnostic_cutoff_count':sum(abs(v)>threshold for v in relations),
        'policy':policy}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    args.output.mkdir(parents=True,exist_ok=False)
    fixtures=json.loads((ROOT/'numerical_audit/fixtures/manifest.json').read_text())
    inputs=[]
    for case,info in sorted(fixtures.items()):
        p=ROOT/'numerical_audit/fixtures'/case/info['tensor']
        if sha(p)!=info['files'][info['tensor']]:raise RuntimeError('Changed fixture.')
        inputs.append((case,p,'frozen_fixture'))
    recovery=ROOT/'numerical_audit/h2o_recovery/data'
    hmanifest=json.loads((recovery/'manifest.json').read_text())
    for basis,info in sorted(hmanifest['bases'].items()):
        for rel,digest in sorted(info['files_sha256'].items()):
            if rel.startswith('unpruned/') and rel.endswith('.tensors.npz'):
                p=recovery/basis/rel
                if sha(p)!=digest:raise RuntimeError('Changed H2O tensor.')
                inputs.append((p.name.removesuffix('.tensors.npz'),p,'reconstructed_H2O_unpruned'))
    sources=[Path(__file__),HERE/'exact_blocks.py',HERE/'protocol.md',ROOT/'analysis/fermionic_coefficients.py']
    hashes={str(p):sha(p) for p in sources}
    write(args.output/'manifest.json',{'thresholds':THRESHOLDS,'input_count':len(inputs),
        'inputs':[{'case':c,'path':str(p),'sha256':sha(p),'provenance':kind} for c,p,kind in inputs],
        'source_sha256':hashes,'identity_convention':'exclude scalar words in both numerator and denominator',
        'warning':'Thresholded counts are diagnostics, not exact physical Coulomb certificates. Raw occurrences retain every word in each retained whole parent; tiny final coefficients may then be excluded only in the diagnostic denominator.'})
    for p in sources:shutil.copy2(p,args.output/p.name)
    write(args.output/'exact_examples.json',exact_examples())
    rows=[];start=time.monotonic()
    for case,p,kind in inputs:
        interaction,n=baseline.load_interaction_operator(p);raw=fermion_from_interaction(interaction)
        # Keep the exact source tensors, so results survive source-directory relocation.
        target=args.output/'inputs'/p.name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target)
        counts=[count_one(raw,t) for t in THRESHOLDS]
        row={'case':case,'n_qubits':n,'provenance':kind,'tensor_sha256':sha(p),'counts':counts}
        rows.append(row);write(args.output/f'{case}.json',row)
        print('COUNT',len(rows),len(inputs),case,n,counts[2]['raw_over_thresholded_four'],flush=True)
    write(args.output/'results.json',rows)
    lines=['# Molecular parent-Pauli duplication audit','',
       'Threshold-sensitive binary-input measurements are separate from the exact physical theorem. At parent cutoff 1e-12, the table shows exact-binary and thresholded (1e-12) denominators; tiny broken tensor symmetries need not cancel exactly. The ratio after global consolidation is exactly 1 for each finite retained target.','',
       '| Case | Qubits | Raw occurrences | Exact U | Threshold U | Whole raw/exact | Four raw/exact | Four raw/threshold | Quartet-grouped/exact |',
       '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    def fmt(v):return '-' if v is None else f'{v:.6g}'
    for row in rows:
        c=row['counts'][2]
        lines.append(f'| {row["case"]} | {row["n_qubits"]} | {c["raw_occurrences"]} | {c["exact_binary_consolidated_words"]} | {c["thresholded_consolidated_words"]} | '+ ' | '.join(fmt(c[k]) for k in ('raw_over_exact_binary_total','raw_over_exact_binary_four','raw_over_thresholded_four','quartet_grouped_over_exact_binary_total'))+' |')
    lines+=['','The independently enumerated 24 permutations give three classes of eight with a shared eight-word JW support. Exact physical weights (1,3,2), (0,2,2), and (1,1,0) give raw/consolidated ratio 4; unrestricted weights (3,5,7) give 3. No finite molecular series proves a universal asymptotic limit.','']
    (args.output/'REPORT.md').write_text('\n'.join(lines))
    for p,digest in hashes.items():
        if sha(p)!=digest:raise RuntimeError('Source changed while counting.')
    write(args.output/'completion.json',{'status':'complete','cases':len(rows),'threshold_rows':len(rows)*len(THRESHOLDS),'seconds':time.monotonic()-start})

if __name__=='__main__':main()
