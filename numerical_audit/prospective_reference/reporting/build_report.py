"""Render already-verified frozen evidence; never construct or select a schedule."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

ARMS=('current','pair_reference','optimized_JD','optimized_J0','pauli_JD',
      'pauli_JD_large','jw_signed','jw_magnitude','symmetric_JD')
LABELS=('Current placement','Two-level correction','Conserving J_D insertion',
        'Conserving J_0 insertion','Pauli J_D insertion','Pauli J_D larger search',
        'JW signed','JW magnitude','Symmetric J_D reuse')
CASES=('Li2_r2p95','F2_r1p55','H2O_scale1p10')


def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def csv_write(p, rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join('---' for _ in headers)+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def build(results,out):
    summary=read(results/'summary.json');decision=summary['decision']
    assert summary['status']=='passed' and summary['verification']['new_replays']==216
    rows=[read(p) for p in sorted((results/'rows').glob('*.json'))]
    assert len(rows)==216
    out.mkdir(parents=True,exist_ok=True)
    csv_write(out/'all_rows.csv',rows)
    points=[(c,s) for c in CASES for s in ('HF','near')]
    lookup={(r['case'],r['state'],r['arm'],r['budget_units']):r for r in rows}
    contexts={}
    for c in CASES:
        for p in read(results/'designs'/(c+'.json'))['points']:contexts[c,p['state']]=p['methods']
    costs=[]
    for c,s in points:
        for arm in ARMS:
            m=contexts[c,s][arm]
            costs.append(dict(case=c,state=s,arm=arm,setup_cpu_seconds=m['setup_cpu_seconds'],
                search_cpu_seconds=m['search_cpu_seconds'],
                total_cpu_seconds=m['setup_cpu_seconds']+m['search_cpu_seconds'],
                setup_wall_seconds=m['setup_seconds'],search_wall_seconds=m.get('seconds',0.),
                examined_positions=m['evaluations'],reuses_order=m.get('reuses_order','')))
    csv_write(out/'construction_costs.csv',costs)
    contrasts=[dict(budget_units=int(b),comparison=k,**v) for b,block in decision['contrasts'].items() for k,v in block.items()]
    csv_write(out/'contrasts.csv',contrasts)
    primary=[[LABELS[i]]+[f"{lookup[c,s,arm,10]['one_minus_overlap']:.5e}" for c,s in points]
             for i,arm in enumerate(ARMS)]
    primary_table=table(['Method']+[c+' '+s for c,s in points],primary)
    floor_points=', '.join(c+' '+s for c,s in points if lookup[c,s,'current',10]['one_minus_overlap']<=1e-12)
    control=decision['contrasts']['10']
    assert all(lookup[c,s,'optimized_JD',b]['one_minus_overlap']==lookup[c,s,'optimized_J0',b]['one_minus_overlap']
               for c,s in points for b in (5,10,20,50))
    contrast_table=table(['Numerator / denominator','B/U','W/T/L','Median error ratio','Excluded floors'],
        [[r['comparison'],r['budget_units'],f"{r['wins']}/{r['ties']}/{r['losses']}",
          '—' if r['median_ratio'] is None else f"{r['median_ratio']:.5g}",r['floor_exclusions']]
          for r in contrasts])
    gate_table=table(['Predeclared requirement','Observed','Pass'],[
        ['At least 4/6 wins at 10U',str(decision['primary']['wins'])+'/6',decision['conditions']['at_least_four_wins']],
        ['No strict loss exceeds 2×',f"{decision['maximum_strict_loss_ratio']:.6g}×" if decision['strict_loss_ratios'] else 'No strict losses (guard defaults to 1)',decision['conditions']['no_loss_exceeds_twofold']],
        ['Median paired total construction CPU ratio ≤2×',f"{decision['median_construction_CPU_ratio']:.6g}×",decision['conditions']['median_construction_CPU_at_most_twofold']]])
    cost_table=table(['Input','State','Current CPU s','Correction CPU s','Ratio'],[
        [r['case'],r['state'],f"{r['current_cpu']:.6g}",f"{r['corrected_cpu']:.6g}",f"{r['ratio']:.6g}"]
        for r in decision['costs']])
    absolute_max=lambda field,arms: max(r[field] for r in rows if r['arm'] in arms)
    conserving=('current','pair_reference','optimized_JD','optimized_J0','symmetric_JD')
    verification=summary['verification']
    report=f"""# Prospective ordering decision: three new molecular inputs

The fixed two-level reference correction **{'passes' if decision['correction_passes'] else 'fails'}** the predeclared development gate. At the primary rotation budget B=10U it gives {decision['primary']['wins']} wins, {decision['primary']['ties']} ties and {decision['primary']['losses']} losses against the current diagonal-placement construction. The median paired error ratio is {decision['primary']['median_ratio']:.6g}; its total construction CPU ratio is {decision['median_construction_CPU_ratio']:.6g}. These are six related determinant inputs from three molecular families, not a statistical validation of general usefulness.

{gate_table}

This experiment advances the search for an ordering by testing a fixed, inexpensive-to-describe reference correction against optimized controls on fresh physical inputs. The correction was fixed before generating these tensors; its outcome is retained without retuning or replacing it. The result does not identify an optimal ordering or isolate conservation as the cause of accuracy gains.

The correction's primary median ratio uses {decision['primary']['ratios_retained']} above-floor baselines; {floor_points} are floor-limited and supply little discrimination. Its roughly {100*(1-decision['primary']['median_ratio']):.2f}% median error reduction is small compared with the observed construction overhead. The development decision is to retain current placement and not adopt this correction on this evidence. At 10U, current placement beats conserving J_D insertion on {control['current_over_optimized_JD']['wins']} inputs with {control['current_over_optimized_JD']['ties']} ties, but loses to the larger-search Pauli and symmetric controls on those same four informative inputs. Current/control median error ratios are {control['current_over_pauli_JD_large']['median_ratio']:.2f} and {control['current_over_symmetric_JD']['median_ratio']:.2f} respectively. J_D and J_0 insertion produce the same physical errors here; this panel does not isolate a weighting advantage. These are scoped comparisons, with unequal classical work for the larger Pauli search and a different formula for the symmetric control.

## What is being ordered

The order is over complete fermionic flip-mask units after exact CAR/Jordan–Wigner consolidation, with each retained diagonal Pauli word movable separately. Complete off-diagonal units preserve particle number and spin projection. Both placement candidates fix the off-unit order by decreasing coefficient l1 norm, then choose diagonal positions with one sweep. This is an ordering of Hamiltonian evolution terms; the orbital basis and spin-orbital numbering stay fixed within every comparison. The unchanged exact one-pass diagonal step uses a newly declared coefficient-l1 initialization here, rather than the inherited optimized off-order used in the earlier 77/84 study.

The current score transports its leading commutator through diagonal D. The correction uses A=D+h(|u><v|+|v><u|), choosing the strongest retained coupling from the input determinant with fixed ties. It uses two commutator columns and a 2×2 eigensystem, without diagonalizing H for search. The reference affects scoring only. Coherent response amplitudes are added before squaring; corrected scores use binary64 complex arithmetic and do not inherit the current score's exact supplied-weight comparison guarantee.

## Panel, physical costs and controls

Neutral singlet RHF/STO-6G: Li2 at 2.95 Å, 2 active electrons/6 spin orbitals; F2 at 1.55 Å, 6/8; H2O at 1.10× the declared geometry, 4/8. Freeze lowest occupied orbitals and use the next active window. Inputs are HF and the nearest same-(Nα,Nβ) determinant (Hamming distance, then occupation-integer tie). T=1 inverse Hartree, hbar=1. Li2 and F2 also change basis/active space relative to earlier fixtures, so changes in ranking cannot be assigned to geometry alone. The full geometries, orbital metadata, target cutoff/scalar conventions and selection rules are in the frozen protocol and evidence.

U is the retained nonidentity Pauli-word count: Li2=33, F2=60, H2O=104. Budgets are 5U, 10U, 20U and 50U ideal Pauli rotations. First-order arms use floor(B/U) steps. Symmetric execution uses the largest affordable step count after exact identical-boundary fusion. It reuses the conserving J_D first-order order; it is not optimized for second order. Ideal rotations do not establish two-qubit gates, depth or hardware performance.

Conserving J_D/J_0 insertion and individual-Pauli J_D insertion use two fixed starts and the same examined-position cap E=(k+1)(m+1) as placement, counting every position, with at most ten accepted moves per start. Equal positions do not match arithmetic or CPU. The larger Pauli control has an explicitly unequal cap 20P(P−1) and the same move limit. Fixed JW signed and magnitude orders are included. All nine arms implement exactly the same retained target, and every propagation retains the entire 2^n space, including leakage.

## Absolute errors at the primary budget

E=1−|normalized overlap|, computed by a stable residual formula. Smaller is better. The table shows all nine arms and all six inputs at 10U; [all 216 rows](all_rows.csv) retain every budget, absolute error, leakage, steps, rotations, timing and state hash.

{primary_table}

![All primary absolute errors](primary_errors.png)

## Correction across budgets and construction costs

![Correction accuracy ratios and total construction cost](correction_tradeoff.png)

Wins/losses use absolute error differences beyond ±1e−12. Descriptive error ratios omit denominators ≤1e−12; the protective loss condition retains strict losses using corrected/max(current,1e−12). Both frozen implementations agree on all decision predicates and contrasts. The graph identifies excluded ratios in its legend rather than dividing by a numerical floor.

{cost_table}

CPU costs include method-specific reference/features/setup plus search; they exclude shared integral generation and physical-error validation. [All construction costs](construction_costs.csv) separate setup, search, examined positions and order reuse; [all contrasts](contrasts.csv) compare both placement candidates with every other arm at every budget. The symmetric arm reuses the J_D construction, so its cost is inherited rather than a second independent search. One timing per input under a quarter-core shared-machine quota is a development measurement, not a controlled speedup benchmark.

## Full contrast inventory

W/T/L refers to the numerator method's lower/tied/higher error. A ratio below one favors that method. Every comparator is retained, including the larger-search Pauli and symmetric controls.

{contrast_table}

## Verification and scope

All 216 circuits were independently replayed with literal NumPy Pauli rotations; symmetric circuits were replayed without boundary fusion and checked against independently counted fused budgets. Independent sparse exponential actions check exact references. Raw saved tensors are bridged through literal CAR to the retained target under the bounded declared Hermitian-adjustment/parent-cutoff allowance; orbital orthogonality, HF energy, Kronecker/bit-indexing, complete-unit conservation, partner/state choice, current exact search ledger, response scores, schedules and frozen source hashes are checked separately. The verifier imports neither QHAT nor the search/propagation engine. This complements source review; it is not formal verification of the executable.

Maximum number leakage in conserving arms: {absolute_max('number_leakage',conserving):.6g}; maximum spin-projection leakage: {absolute_max('spin_projection_leakage',conserving):.6g}. Individual-Pauli arms retain their observed leakage in the rows. Source commit before generation: `{summary['scientific_source_commit']}`. Evidence summary status: `{verification['status']}`; independent replays: {verification['new_replays']}. [Report provenance](provenance.json) binds the frozen summary, renderer and exported tables.

The original supervisor used a 4 GiB cap and 6 GiB start threshold for generation and construction. Before any physical propagation, a documented operational amendment tightened the remaining jobs to 512 MiB (448 MiB high), retained quarter-core/no-swap settings and the 5 GiB pause threshold, and lowered start/resume to 5.5 GiB. The original waiting propagation launcher was interrupted before launching a service. Scientific source, inputs, schedules and construction timings were unchanged. The retained OPERATIONAL_AMENDMENT.json and supervisor bytes disclose this deviation from the original resource protocol.

The new evidence covers small real determinant inputs at one time, three related molecular families and fixed orbitals. It does not validate correlated/coherent inputs, large spaces, a new encoding, compiled cost, second-order optimization, general timing performance or physical-error optimality. The separate JW/BK formalization is unchanged. See the paper companion for the development decision and slide-derived future work.
"""
    (out/'README.md').write_text(report)
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.size':10,'savefig.dpi':170})
    values=np.array([[lookup[c,s,a,10]['one_minus_overlap'] for c,s in points] for a in ARMS])
    fig,ax=plt.subplots(figsize=(12,6.7),layout='constrained')
    image=ax.imshow(np.maximum(values,1e-16),norm=LogNorm(vmin=1e-16,vmax=max(1e-4,float(values.max()))),cmap='viridis_r',aspect='auto')
    ax.set_xticks(range(6),[c.replace('_r','\nr').replace('_scale','\nscale')+'\n'+s for c,s in points])
    ax.set_yticks(range(9),LABELS)
    for i in range(9):
        for j in range(6):
            value=values[i,j];rgba=image.cmap(image.norm(max(value,1e-16)))
            ax.text(j,i,f'{value:.1e}',ha='center',va='center',fontsize=9,
                    color='black' if sum(rgba[:3])>1.5 else 'white')
    ax.set_title('All arms at 10U ideal rotations: absolute normalized overlap error\nColor saturates below 1e−16; annotations retain measured values')
    fig.colorbar(image,ax=ax,label='E = 1 − |normalized overlap|')
    for extension in ('png','pdf'):fig.savefig(out/('primary_errors.'+extension))
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for i,(c,s) in enumerate(points):
        budgets=(5,10,20,50);ratios=[]
        for b in budgets:
            base=lookup[c,s,'current',b]['one_minus_overlap']
            ratios.append(lookup[c,s,'pair_reference',b]['one_minus_overlap']/base if base>1e-12 else np.nan)
        excluded=sum(np.isnan(ratios))
        label=c+' '+s+(f' ({excluded}/4 floor-excluded)' if excluded else '')
        axes[0].plot(budgets,100*(1-np.array(ratios)),'o-',label=label)
    axes[0].axhline(0,color='black',ls='--');axes[0].set(xlabel='Rotation budget / U',ylabel='Relative error reduction (%)',title='Correction: modest improvements above floor')
    axes[0].set_xticks([5,10,20,50],['5','10','20','50']);axes[0].legend(fontsize=7)
    cost_ratios=[next(r['ratio'] for r in decision['costs'] if (r['case'],r['state'])==(c,s)) for c,s in points]
    axes[1].bar(range(6),cost_ratios,color='slateblue')
    axes[1].axhline(2,color='black',ls='--',label='Predeclared median limit = 2')
    axes[1].axhline(decision['median_construction_CPU_ratio'],color='darkorange',ls=':',label='Observed paired median')
    axes[1].set(yscale='log',ylabel='Corrected / current total CPU',title='Setup + search CPU: one timing per input')
    axes[1].set_xticks(range(6),[c.replace('_r','\nr').replace('_scale','\nscale')+'\n'+s for c,s in points],fontsize=8)
    axes[1].legend(fontsize=7)
    for extension in ('png','pdf'):fig.savefig(out/('correction_tradeoff.'+extension))
    plt.close(fig)
    provenance=dict(schema=1,summary_sha256=digest(results/'summary.json'),
        frozen_protocol_sha256=digest(results/'PROTOCOL_FROZEN.json'),
        schedule_freeze_sha256=digest(results/'ALL_SCHEDULES_FROZEN.json'),
        renderer_sha256=digest(Path(__file__)),scientific_source_commit=summary['scientific_source_commit'],
        artifact_sha256={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='provenance.json'})
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),gate=decision['correction_passes'],output=str(out))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.results,args.output)
