"""Build compact tables from independently verified issue #4 evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics


ARMS=('hf_current','state_current','state_JD','state_pauli_JD','hf_pauli_large',
      'hf_symmetric','state_symmetric','jw_signed','jw_magnitude')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table(path,rows):
    if not rows:raise ValueError('Empty report table')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def build(results,report):
    verified=json.loads((results/'VERIFIED.json').read_text())
    for p,h in verified['artifacts'].items():
        if sha(results/p)!=h:raise ValueError('Verified artifact changed')
    summary=json.loads((results/'SUMMARY.json').read_text())
    if summary['rows']!=432 or verified['rows']!=432 or verified['status']!='passed':raise ValueError('Incomplete report inputs')
    report.mkdir(parents=True,exist_ok=True)
    rows=[json.loads(p.read_text()) for p in sorted((results/'rows').glob('*.json'))]
    fields=('case','state','T','arm','budget_units','formula','steps','rotations','rotation_cap',
            'error','infidelity','norm','number_leakage','spin_projection_leakage','state_sha256')
    table(report/'all_rows.csv',[{k:r[k] for k in fields} for r in rows])
    table(report/'contrasts.csv',[{**r,'families':';'.join(r['families'])} for r in summary['contrasts']])
    state_rows=[];costs=[];preparation=[];campaign_costs=[]
    for p in sorted((results/'designs').glob('*.json')):
        design=json.loads(p.read_text());prep=design['state_preparation'];case=design['case']
        preparation.append(dict(case=case,cpu_seconds=prep['cpu_seconds'],wall_seconds=prep['wall_seconds']))
        for label,meta in prep['labels'].items():
            state_rows.append(dict(case=case,state=label,participation_ratio=meta['participation_ratio'],
                max_determinant_weight=meta['maximum_determinant_weight'],
                one_RDM_defect=meta['one_RDM_idempotency_defect'],
                correlated_beyond_single_Slater=meta['correlated_beyond_single_Slater'],
                eigen_residual=meta.get('eigen_residual',''),
                natural_occupations=';'.join(str(x) for x in meta['natural_occupations'])))
        for point in design['points']:
            shared_CPU=sum(c['cpu_seconds'] for c in point['setup_costs'].values())
            search_CPU=sum(point['methods'][name]['search_cost']['cpu_seconds'] for name in
                           ('state_current','state_JD','state_pauli_JD'))
            campaign_costs.append(dict(case=case,state=point['state'],T=point['T'],
                actual_feature_setup_CPU_seconds=shared_CPU,actual_search_CPU_seconds=search_CPU,
                actual_total_construction_CPU_seconds=shared_CPU+search_CPU))
            for name in ARMS:
                method=point['methods'][name];setup='pauli_features' if method['partition']=='pauli' else 'conserving_features'
                inherited=name=='state_symmetric'
                new_search=name.startswith('state_')
                costs.append(dict(case=case,state=point['state'],T=point['T'],arm=name,
                    training=method['training'],order_reuse='state_JD' if inherited else '',
                    setup_CPU_seconds=point['setup_costs'][setup]['cpu_seconds'] if new_search else 0.,
                    search_CPU_seconds=method.get('search_cost',{}).get('cpu_seconds',0.),
                    evaluations=method['evaluations'],
                    historical_HF_setup_CPU_seconds=method.get('historical_setup_cpu_seconds',''),
                    historical_HF_search_CPU_seconds=method.get('historical_search_cpu_seconds','')))
    table(report/'state_preparation.csv',preparation);table(report/'states.csv',state_rows);table(report/'construction_costs.csv',costs)
    table(report/'campaign_construction_costs.csv',campaign_costs)
    # Independent descriptive aggregation from raw physical rows.
    lookup={(r['case'],r['state'],r['T'],r['arm'],r['budget_units']):r['error'] for r in rows}
    aggregate=[]
    for budget in (10,30):
        for arm in ARMS[1:]:
            pairs=[(r['error'],lookup[r['case'],r['state'],r['T'],'hf_current',budget])
                   for r in rows if r['arm']==arm and r['budget_units']==budget]
            ratios=[a/b for a,b in pairs if b>1e-12]
            aggregate.append(dict(arm=arm,budget_units=budget,wins=sum(a-b < -1e-12 for a,b in pairs),
                ties=sum(abs(a-b)<=1e-12 for a,b in pairs),losses=sum(a-b > 1e-12 for a,b in pairs),
                median_ratio=statistics.median(ratios) if ratios else None,floor_exclusions=len(pairs)-len(ratios)))
    table(report/'aggregate.csv',aggregate)
    lines=['# Derived ordering on correlated/coherent state and time holdouts','',
        'This prospective panel compares transfer of frozen HF-trained schedules with newly state-aware diagonal-response placement and insertion. The orbital basis and retained target stay fixed within each family. All 432 circuits passed independent full-space replay. These are 24 related state/time inputs from three molecular families, not independent statistical replicates.','',
        '## State and objective scope','',
        f"The inputs are sector ground and first-excited states and (g+e)/sqrt(2), (g+i e)/sqrt(2), at T=0.37 and 1.73 inverse Hartree. One-particle density spectra classify {sum(r['correlated_beyond_single_Slater'] for r in state_rows)}/{len(state_rows)} selected states as correlated beyond a single Slater determinant under the declared 1e-8 idempotency-defect threshold. Participation ratios separately describe their expansion in the fixed orbital basis. These sector labels do not claim unrestricted global ground/excited states.",'',
        'State-aware J_D sums complex transported commutator amplitudes before the parallel-state projection and squared norm. It retains coherent interference. This is the generalized derived diagonal-response rule; earlier historical signed-parent low-energy tests do not validate it. The failed two-level reference correction is unchanged and excluded from new tuning.','',
        '## Coverage and all controls','',
        '| State | T | Informative families at 10U | Coverage |','| --- | --- | --- | --- |']
    for c in summary['coverage']:lines.append(f"| {c['state']} | {c['T']} | {', '.join(c['informative_families']) or 'none'} | {c['status']} |")
    lines+=['','Coverage requires at least two of three families with frozen-HF current-placement E>1e-12 in each state/time stratum. Failed strata remain inconclusive and are not replaced. All absolute errors and losses remain in the tables. The denominator floor is 1e-12; ratios exclude and count floor-limited baselines.','',
        '| Method / frozen-HF current | Budget | Wins / ties / losses (24 inputs) | Median error ratio | Floor exclusions |',
        '| --- | --- | --- | --- | --- |']
    for r in aggregate:
        ratio='excluded' if r['median_ratio'] is None else f"{r['median_ratio']:.6g}"
        lines.append(f"| {r['arm']} | {r['budget_units']}U | {r['wins']}/{r['ties']}/{r['losses']} | {ratio} | {r['floor_exclusions']} |")
    lines+=['','hf_current, hf_pauli_large, hf_symmetric and fixed JW controls transfer the previous HF/T=1 schedules unchanged. state_current runs generalized placement. state_JD and state_pauli_JD newly optimize the current state/time objective. Both symmetric arms reuse first-order schedules; neither is a second-order-aware optimization. State Pauli search has an explicitly larger position allowance, while the frozen HF Pauli control retains the earlier still-larger search; this is not a comparison of matched CPU budgets.','',
        '## Absolute primary errors, by family','',
        '| Family | State | T | '+' | '.join(ARMS)+' |','| --- | --- | --- | '+' | '.join(['---']*len(ARMS))+' |']
    for case,state,T in sorted({(r['case'],r['state'],r['T']) for r in rows}):
        lines.append(f'| {case} | {state} | {T} | '+' | '.join(f'{lookup[case,state,T,a,10]:.6e}' for a in ARMS)+' |')
    lines+=['','[All rows](all_rows.csv) retain both ideal-rotation budgets, formula, affordable steps, actual rotations, absolute state error, infidelity, norm and number/joint spin-resolved-sector leakage. The historical spin_projection_leakage column denotes the joint (N_alpha,N_beta) sector, not pure S_z leakage. [All contrasts](contrasts.csv) retain each family and state/time stratum. Symmetric circuits use exact central/boundary fusion and their actual affordable integer steps. Individual-Pauli controls may leak, and propagation always retains the entire Hilbert space.','',
        '## Preparation and classical work','',
        '| Family | Preparation + correlation diagnostics CPU s | Wall s |','| --- | --- | --- |']
    for r in preparation:lines.append(f"| {r['case']} | {r['cpu_seconds']:.6g} | {r['wall_seconds']:.6g} |")
    lines+=['','Preparation includes retained-H construction, sector eigensystem, explicit state vectors and density-matrix diagnostics. It is charged once per family; these exact-diagonalization states do not supply an efficient or compiled hardware preparation routine. Original molecular integral generation is reused from authenticated prior evidence.','',
        '[State diagnostics](states.csv), [stand-alone arm construction costs](construction_costs.csv) and [actual campaign construction costs](campaign_construction_costs.csv) retain the detailed costs. Conserving feature setup is shared by placement/insertion; the stand-alone table assigns the full feature cost to each method, while the campaign table counts actual setup/search exactly once. state_symmetric inherits state_JD construction, not another search. Transferred HF schedules have zero new search, with historical construction costs listed separately. These single timing observations under a shared-machine quota do not establish controlled speedups or scaling.','',
        '## Verification and limitations','',
        f"Independent literal-Pauli/BCH scoring and unfused full-space replay checked every row: maximum state-vector difference {verified['max_state_delta']:.6g}, physical metric difference {verified['max_metric_delta']:.6g}, response-score difference {verified['max_score_delta']:.6g}. Synthetic tests include complex quadrature, coherent-versus-mixture separation, determinant reduction, updates, degenerate state selection, Slater correlation, pulse fusion and resource guards.",
        '', 'The frozen protocol and all state/schedule hashes precede physical errors. The earlier pre-error correlation-diagnostic amendment left state-file bytes and all 216 schedules unchanged. The numerical source and active inputs are authenticated in the result freeze. No schedule is chosen from physical errors, and no trial is replaced after viewing errors.','',
        'This bounded panel establishes only the reported state/time and family behavior. It does not isolate conservation as an accuracy mechanism, prove J_D ranks physical errors, demonstrate matched classical efficiency, provide Lean proofs, optimize a second-order objective, or measure compiled/hardware resources. Those remain separate issues #5–#9. The completed factorial, Tranter checks and initial 8! pilot remain prior completed evidence.']
    (report/'README.md').write_text('\n'.join(lines)+'\n')
    manifest={p.name:sha(p) for p in sorted(report.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}
    (report/'MANIFEST.json').write_text(json.dumps({'files':manifest,'schedule_freeze_sha256':sha(results/'SCHEDULES_FROZEN.json'),'verification_sha256':sha(results/'VERIFIED.json'),'summary_sha256':sha(results/'SUMMARY.json')},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();build(a.results,a.report)
