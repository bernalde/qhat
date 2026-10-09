#!/usr/bin/env python3
"""Summarize the frozen H2O reconstruction, keeping archived data separate."""
import os
os.environ['MPLCONFIGDIR']='/tmp/qhat-h2o-mpl';os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,csv,json,math,hashlib,shutil
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
METHODS=('fermionic_signed','jw_signed','jw_magnitude')
LABELS={'fermionic_signed':'Owned, signed parent','jw_signed':'JW, signed','jw_magnitude':'JW, magnitude'}
COLORS={'fermionic_signed':'#0072B2','jw_signed':'#D55E00','jw_magnitude':'#009E73'}
def require(c,msg):
    if not c:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,o):
    with Path(p).open('w') as f:json.dump(o,f,indent=2,allow_nan=False);f.write('\n')
def outcome(x,y):return 'tie' if abs(x-y)<=1e-12 else ('win' if x<y else 'loss')
def historical(x,y):
    if min(x,y)<=1e-12:return 'below_floor'
    return 'tie' if math.isclose(x,y,rel_tol=1e-6,abs_tol=1e-14) else ('win' if x<y else 'loss')
def csvout(path,rows):
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sweep',type=Path,required=True);ap.add_argument('--low-energy',type=Path,required=True);ap.add_argument('--historical-csv',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    results=[json.loads(p.read_text()) for p in sorted(a.sweep.glob('H2O*/results.json'))];rows=[row for r in results for row in r['rows']]
    require(len(rows)==126,'Expected126 HF rows');low=json.loads((a.low_energy/'rows.json').read_text());require(len(low)==36,'Expected36 lowenergy rows')
    designs={p.parent.name:json.loads(p.read_text()) for p in a.sweep.glob('H2O*/design.json')}
    checks=[json.loads(p.read_text()) for p in a.sweep.glob('H2O*/checks.json')]
    groups={}
    for row in rows:groups.setdefault((row['case'],row['tensor_policy'],row['steps']),{})[row['ordering']]=row
    comparisons=[]
    for (case,policy,r),g in sorted(groups.items()):
        require(set(g)==set(METHODS),'Missing comparison')
        E={m:g[m]['one_minus_overlap'] for m in METHODS}
        comparisons.append({'case':case,'basis':g['fermionic_signed']['basis'],'tensor_policy':policy,'n_qubits':g['fermionic_signed']['n_qubits'],
                            'n_electrons':g['fermionic_signed']['n_electrons'],'steps':r,'errors':E,
                            'outcomes':{m:outcome(E['fermionic_signed'],E[m]) for m in METHODS[1:]},
                            'historical_convention_outcomes':{m:historical(E['fermionic_signed'],E[m]) for m in METHODS[1:]},
                            'ratios':{m:E['fermionic_signed']/E[m] if E[m]>1e-12 else None for m in METHODS[1:]}})
    policy_summary={}
    for policy in ('historical_policy','unpruned'):
        primary=[c for c in comparisons if c['tensor_policy']==policy and c['steps']==100]
        policy_summary[policy]={'cases':len(primary),'best_method_sets':dict(Counter('+'.join(sorted(m for m in METHODS if c['errors'][m]-min(c['errors'].values())<=1e-12)) for c in primary)),'new_convention':{m:dict(Counter(c['outcomes'][m] for c in primary)) for m in METHODS[1:]},
                               'historical_convention':{m:dict(Counter(c['historical_convention_outcomes'][m] for c in primary)) for m in METHODS[1:]},
                               'median_ratio':{m:float(np.median([c['ratios'][m] for c in primary if c['ratios'][m] is not None])) for m in METHODS[1:]}}
    pairdiff=[]
    for case in sorted({r['case'] for r in rows}):
        x=designs[case+'__historical_policy'];y=designs[case+'__unpruned'];gx=groups[(case,'historical_policy',100)];gy=groups[(case,'unpruned',100)]
        cx={tuple(tuple(t) for t in k):complex(float.fromhex(re),float.fromhex(im)) for k,re,im in x['identity_free_coefficients']};cy={tuple(tuple(t) for t in k):complex(float.fromhex(re),float.fromhex(im)) for k,re,im in y['identity_free_coefficients']}
        delta_l1=float(sum(abs(cx.get(k,0)-cy.get(k,0)) for k in set(cx)|set(cy)))
        pairdiff.append({'case':case,'retained_H_difference_coefficient_l1':delta_l1,'same_retained_H':x['hamiltonian_sha256']==y['hamiltonian_sha256'],
                         'same_schedule':{m:x['schedule_sha256'][m]==y['schedule_sha256'][m] for m in METHODS},
                         'max_primary_E_difference':max(abs(gx[m]['one_minus_overlap']-gy[m]['one_minus_overlap']) for m in METHODS)})
    anchors=[]
    for case,policy,r in groups:
        if r!=50:continue
        for m in METHODS:
            E={str(r):groups[(case,policy,r)][m]['one_minus_overlap'] for r in (50,100,200)}
            anchors.append({'case':case,'tensor_policy':policy,'ordering':m,'E':E,'E50_over_E100':E['50']/E['100'],'E100_over_E200':E['100']/E['200']})
    oldall=list(csv.DictReader(a.historical_csv.open()));old=[r for r in oldall if r['molecule']=='H2O'];require(len(old)==15,'Expected15 oldH2Orows')
    historical_comparison=[]
    for o in old:
        g=groups[(o['case_id'],'historical_policy',100)]
        oe={m:float(o[{'fermionic_signed':'fermionic_aware_one_minus_overlap','jw_signed':'jw_signed_one_minus_overlap','jw_magnitude':'jw_magnitude_one_minus_overlap'}[m]]) for m in METHODS}
        for m in METHODS:
            new=g[m]['one_minus_overlap'];historical_comparison.append({'case':o['case_id'],'ordering':m,'archived_E':oe[m],'reconstructed_E':new,
              'absolute_difference':abs(new-oe[m]),'relative_difference':abs(new-oe[m])/oe[m] if oe[m]>1e-12 else None,
              'interpretation':'matching configurations, missing historical tensors; no identical-input claim'})
    lg={}
    for row in low:lg.setdefault((row['case'],row['tensor_policy'],row['state']),{})[row['ordering']]=row
    lowcounts={p:{m:dict(Counter(outcome(g['fermionic_signed']['one_minus_overlap'],g[m]['one_minus_overlap']) for (c,pp,s),g in lg.items() if pp==p)) for m in METHODS[1:]} for p in ('historical_policy','unpruned')}
    lowchecks=[json.loads(p.read_text()) for p in a.low_energy.glob('H2O*/results.json')]
    lowdesigns=[json.loads(p.read_text()) for p in a.low_energy.glob('H2O*/design.json')]
    summary={'protocol':'h2o-reconstruction-sweep-v1','hf_rows':len(rows),'low_energy_rows':len(low),'primary_cases':15,'tensor_policies':2,
             'policy_summary':policy_summary,'coefficient_policy_sensitivity':pairdiff,'comparisons':comparisons,'convergence_anchors':anchors,
             'low_energy_outcomes':lowcounts,'archived_table_sha256':sha(a.historical_csv),'historical_outcomes':dict(Counter(o['outcome'] for o in old)),
             'archived_vs_reconstructed':historical_comparison,'archived_vs_reconstructed_max_absolute_difference':max(x['absolute_difference'] for x in historical_comparison),
             'archived_vs_reconstructed_max_relative_difference':max(x['relative_difference'] for x in historical_comparison if x['relative_difference'] is not None),
             'validation':{'max_reference_split_distance':max(c['reference_split_time_state_distance'] for c in checks),
                 'max_affine_full_sample_distance':max(c['sample_pauli_full_vs_affine_distance'] for c in checks),
                 'max_reference_off_spin_sector_frobenius':max(c['off_spin_sector_frobenius'] for c in checks),
                 'max_dense_pilot_distance':max(r.get('vs_independent_dense_pilot_state_distance',0) for r in rows),
                 'max_low_energy_dense_product_distance':max(r['dense_product_state_distance'] for r in low),
                 'max_low_energy_reference_dense_eigh_distance':max(v for c in lowchecks for v in c['checks'].values()),
                 'peak_RSS_GiB':max(r['peak_rss_KiB'] for r in results)/1024**2,
                 'max_norm_defect':max(max(abs(r['exact_norm']-1),abs(r['approx_norm']-1)) for r in rows+low)},
             'max_leakage_by_order':{m:{k:max(r[k] for r in rows+low if r['ordering']==m) for k in ('particle_number_leakage','spin_sector_leakage')} for m in METHODS}}
    write(a.output/'summary.json',summary);write(a.output/'rows.json',rows);write(a.output/'low_energy_rows.json',low)
    csvout(a.output/'rows.csv',rows);csvout(a.output/'low_energy_rows.csv',low);csvout(a.output/'archived_h2o_rows.csv',old);csvout(a.output/'archived_vs_reconstructed.csv',historical_comparison)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(9.0,6.0),sharex='col',layout='constrained')
    for col,basis in enumerate(('sto-6g','hgbs-5')):
        data=sorted([c for c in comparisons if c['basis']==basis and c['tensor_policy']=='historical_policy' and c['steps']==100],key=lambda c:c['n_qubits'])
        x=[c['n_qubits'] for c in data]
        for m in METHODS:
            axes[0,col].plot(x,[c['errors'][m] for c in data],marker='o',markersize=3,color=COLORS[m],label=LABELS[m])
        for m,marker in (('jw_signed','s'),('jw_magnitude','o')):
            axes[1,col].plot(x,[c['ratios'][m] if c['ratios'][m] is not None else np.nan for c in data],marker=marker,markersize=4,color=COLORS[m],label='Owned / '+LABELS[m])
        axes[0,col].axhline(1e-12,color='.5',linestyle=':',linewidth=.8);axes[0,col].set_title(basis.upper())
        axes[0,col].set_yscale('log');axes[1,col].set_yscale('log');axes[1,col].axhline(1,color='.5',linestyle=':',linewidth=.8)
        axes[1,col].set_xticks(x,[f'{c["n_qubits"]}\n({c["n_electrons"]}e)' for c in data]);axes[1,col].set_xlabel('Active qubits (active electrons)')
        for ax in axes[:,col]:ax.grid(alpha=.15)
    axes[0,0].set_ylabel(r'$E=1-|\langle\psi_{\rm exact},\psi_{\rm approx}\rangle|$');axes[1,0].set_ylabel(r'$E_{\rm owned}/E_{\rm comparator}$')
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.suptitle(r'Reconstructed H$_2$O inputs; HF state, $T=1$, $r=100$',fontsize=11)
    fig.savefig(a.output/'h2o_reconstruction.pdf');fig.savefig(a.output/'h2o_reconstruction.png',dpi=180);plt.close(fig)
    # Historical panel is intentionally an independent artifact with provenance.
    fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    for ax,basis in zip(axes,('sto-6g','hgbs-5')):
        oldbasis=sorted([o for o in old if o['basis']==basis],key=lambda o:int(o['n_qubits']))
        for m,key in (('fermionic_signed','fermionic_aware_one_minus_overlap'),('jw_magnitude','jw_magnitude_one_minus_overlap')):
            ax.plot([int(o['n_qubits']) for o in oldbasis],[float(o[key]) for o in oldbasis],marker='o',color=COLORS[m],label=LABELS[m])
        ax.set_yscale('log');ax.set_title(basis.upper());ax.set_xlabel('Active qubits');ax.grid(alpha=.15)
    axes[0].set_ylabel('Archived reported overlap error E');axes[0].legend(fontsize=8)
    fig.suptitle('Archived slide-26 source table (not new simulation)',fontsize=11)
    fig.savefig(a.output/'h2o_archived_panel.pdf');plt.close(fig)
    lines=['# H2O reconstruction sweep and small state control','','These are regenerated inputs, not a bitwise replay of unavailable historical tensors. Each within-case comparison fixes H, state, time, and Pauli exponential count. Active-space changes alter H, electron count, and frozen-core contributions.','','126 HF evolutions completed: 15 cases × two tensor policies × three methods at T=1,r=100, plus three predeclared anchors at r=50/200. A separate prescribed two-case low-energy panel adds 36 evolutions. These repeated settings are not independent molecular samples.','','## Primary outcomes','','Wins/ties/losses for the proposed owned signed schedule; absolute tie tolerance 1e-12. Ratios suppressed when the comparator error is at or below 1e-12.','','| Tensor policy | Comparator | Wins / ties / losses | Median owned/comparator ratio |','|---|---|---|---|']
    for p,s in policy_summary.items():
        for m in METHODS[1:]:
            c=s['new_convention'][m];lines.append(f'| {p} | {m} | {c.get("win",0)} / {c.get("tie",0)} / {c.get("loss",0)} | {s["median_ratio"][m]:.6g} |')
    lines+=['','Both tensor policies use corrected downstream coefficient accumulation and the same final 1e-12 Pauli cutoff. Historical tensor policy preserves the generator’s old truncation, not its downstream metric/accumulation defects. The two policies are sensitivity arms and must not be pooled as independent cases.','','## Absolute primary errors (historical tensor policy)','','| Basis | Qubits / electrons | Owned signed E | JW signed E | JW magnitude E |','|---|---|---|---|---|']
    for c in comparisons:
        if c['tensor_policy']=='historical_policy' and c['steps']==100:lines.append(f'| {c["basis"]} | {c["n_qubits"]} / {c["n_electrons"]} | {c["errors"][METHODS[0]]:.8e} | {c["errors"][METHODS[1]]:.8e} | {c["errors"][METHODS[2]]:.8e} |')
    lines+=['','## Low-energy panel','','For both eight-qubit active models, the lowest fixed-(2,2) state is a singlet and the next is an Ms=0 triplet (S²=2, negligible variance). Their equal coherent superposition has <S²>=1 and variance1. All evolutions retain the full256-dimensional space.','','| Tensor policy | Comparator | Owned wins / ties / losses across six case/state pairs |','|---|---|---|']
    for p,cs in lowcounts.items():
        for m,c in cs.items():lines.append(f'| {p} | {m} | {c.get("win",0)} / {c.get("tie",0)} / {c.get("loss",0)} |')
    lines+=['','## Checks and interpretation','',
        'The proposed schedule is uniquely best among the three in only four of fifteen HF cases. Against JW signed sorting its median error ratio is 2.36224. The selective benefits in larger HGBS-5 active spaces coexist with losses elsewhere; they do not establish a general accuracy advantage.',
        '',f'All 45 historical-policy primary errors agree with the archived table within {summary["archived_vs_reconstructed_max_absolute_difference"]:.6g} absolute. The largest relative discrepancy is {100*summary["archived_vs_reconstructed_max_relative_difference"]:.3g}% for a floor-adjacent error. This substantially supports the specific H2O table while not certifying identical missing historical tensor files.',
        '',f'Policy sensitivity: retained H and schedules are identical in 14 cases. Only the 20-qubit unpruned case adds four Pauli coefficients of about 1.15746e-12; its retained-H coefficient-l1 difference is 4.62985e-12. The maximum primary E change between policies is {max(c["max_primary_E_difference"] for c in pairdiff):.6g}; outcomes are unchanged.',
        '',f'Validation maxima: {json.dumps(summary["validation"],sort_keys=True)}.','',
        'Small HF states reproduce independently dense-verified pilot states. All36 low-energy evolutions additionally agree with independently formed dense product matrices. The full sweep exact reference uses the retained Pauli sum on the complete Pauli-invariant closure, including round-off off-sector terms. A separate checker recomputes all162 final-state metrics through the projector-residual identity and verifies input/design/state hashes and Pauli closure. This checks the large saved results but is not a second independent large-state evolution.',
        '',f'Archived tally: {json.dumps(summary["historical_outcomes"])}. Archived values are retained separately in archived_h2o_rows.csv and h2o_archived_panel.pdf. The new historical-convention classifications are {json.dumps({p:s["historical_convention"] for p,s in policy_summary.items()},sort_keys=True)}.',
        '','See summary.json for all18 anchor/method entries (two step-convergence ratios each), per-case policy equality, archived-versus-reconstructed differences, raw leakage and norms. h2o_reconstruction.pdf displays paired absolute errors and ratios; it is not a qubit-only convergence law.']
    (a.output/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'policy_summary':policy_summary,'low_energy_outcomes':lowcounts,'validation':summary['validation'],'max_leakage_by_order':summary['max_leakage_by_order']},indent=2))
if __name__=='__main__':main()
