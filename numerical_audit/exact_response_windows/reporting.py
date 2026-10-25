"""Deterministic reports retaining every restricted certificate and comparator."""
import argparse,csv,hashlib,json,statistics
from pathlib import Path


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table(path,rows):
    require(bool(rows),'Empty report table')
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def generate(results,report):
    require(not report.exists(),'Fresh report path required');report.mkdir(parents=True)
    verified=json.loads((results/'VERIFIED.json').read_text());require(verified['status']=='passed','Unverified study')
    protocol=json.loads((results/'PROTOCOL_FROZEN.json').read_text());certs=[];rows=[];comparisons=[]
    for path in sorted((results/'certificates').glob('*.json')):
        case=path.stem;c=json.loads(path.read_text());d=json.loads((results/'domains'/path.name).read_text())
        require(c['complete'] and c['lower_bound_integer']==c['incumbent_integer'],'Uncertified result')
        certs.append(dict(case=case,units=len(d['warm_order']),width=d['width'],offset=d['offset'],
            nonzero_internal_pairs=max(d['all_window_pair_counts']),evaluated=c['evaluated'],minimizers=c['minimizers'],
            warm_start_integer=c['warm_start_integer'],incumbent_integer=c['incumbent_integer'],lower_bound_integer=c['lower_bound_integer'],
            objective_denominator=c['objective_denominator'],absolute_gap_integer=c['absolute_gap_integer'],relative_gap=c['relative_gap'],
            initial_JD=c['initial_JD'],incumbent_JD=c['incumbent_JD'],termination=c['termination'],
            setup_CPU_seconds=d['setup_CPU_seconds'],enumeration_CPU_seconds=c['CPU_seconds'],enumeration_wall_seconds=c['wall_seconds'],
            inherited_heuristic_setup_CPU_seconds=d['prior_construction']['optimized_JD']['setup_cpu_seconds'],
            inherited_heuristic_search_CPU_seconds=d['prior_construction']['optimized_JD']['search_cpu_seconds']))
    for path in sorted((results/'rows').glob('*.json')):
        r=json.loads(path.read_text());rows.append({k:r[k] for k in ('case','record_index','arm','formula','budget_units','steps','rotations','rotation_cap','error','infidelity','norm','number_leakage','joint_sector_leakage')}|
            dict(propagation_CPU_seconds=r['propagation_cost']['CPU_seconds'],reference_CPU_seconds=r['reference_cost']['CPU_seconds']))
    require(len(certs)==len(protocol['cases']) and len(rows)==protocol['expected_circuits'],'Incomplete report grid')
    for at,a in enumerate(rows):
        for b in rows[at+1:]:
            if (a['case'],a['budget_units'])!=(b['case'],b['budget_units']):continue
            delta=a['error']-b['error'];floor=min(a['error'],b['error'])<=1e-12
            comparisons.append(dict(case=a['case'],budget_units=a['budget_units'],a_arm=a['arm'],a_formula=a['formula'],b_arm=b['arm'],b_formula=b['formula'],
                a_error=a['error'],b_error=b['error'],a_minus_b=delta,result='win' if delta < -1e-12 else 'loss' if delta>1e-12 else 'tie',
                floor_excluded_ratio=floor,a_over_b=None if floor else a['error']/b['error']))
    require(len(comparisons)==15*len(protocol['cases'])*2,'Missing all-pairs comparisons')
    table(report/'certificates.csv',certs);table(report/'all_rows.csv',rows);table(report/'comparisons.csv',comparisons)
    summary=[]
    for formula in (1,2):
        for budget in (10,30):
            matched=[]
            for case in protocol['cases']:
                a=next(r for r in rows if (r['case'],r['arm'],r['formula'],r['budget_units'])==(case,'certified_window',formula,budget))
                b=next(r for r in rows if (r['case'],r['arm'],r['formula'],r['budget_units'])==(case,'heuristic',formula,budget))
                matched.append((a,b))
            deltas=[a['error']-b['error'] for a,b in matched];ratios=[a['error']/b['error'] for a,b in matched if min(a['error'],b['error'])>1e-12]
            summary.append(dict(formula=formula,budget_units=budget,pairs=len(matched),wins=sum(d < -1e-12 for d in deltas),
                ties=sum(abs(d)<=1e-12 for d in deltas),losses=sum(d>1e-12 for d in deltas),floor_exclusions=len(matched)-len(ratios),
                median_paired_error_ratio=statistics.median(ratios) if ratios else None))
    table(report/'summary.csv',summary)
    lines=['# Exact restricted response windows','',
        f"Complete enumeration certified {len(certs)} new nine-unit windows, each over 9! = 362880 feasible orders. Outside units and their positions remain fixed. All incumbents equal the independently verified lower bounds, with zero exact-integer gap and zero tolerances.",'',
        'The objective is J_D = sum_w rho_w b_w(order)^2/(4 denominator^4), with integer commutator amplitudes and the exact rational values of supplied binary64 rho_w = sinc^2(gap_w/(2 pi denominator)) at T=1. This does not certify exact transcendental sinc weights, the full ordering domain, a molecular ground-state problem, or finite-step physical error. Complete-unit coefficients remain the authenticated rational model coefficients.','',
        'Each window is the contiguous nine-unit subset of the previously frozen HF J_D heuristic order with the most nonzero internal response pairs (earliest offset on ties). Selection uses no physical errors. Prior heuristic/current/larger-Pauli construction timings are inherited records; setup and exhaustive-enumeration costs are newly measured. No matched-classical-cost advantage is claimed.','',
        '| Case | Offset | Initial J_D | Certified J_D | Minimizers | Enumeration CPU s |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for c in certs:lines.append(f"| {c['case']} | {c['offset']} | {c['initial_JD']:.9g} | {c['incumbent_JD']:.9g} | {c['minimizers']} | {c['enumeration_CPU_seconds']:.6g} |")
    lines+=['','All six schedules per input were frozen after certification and before physical evaluation. The grid includes first-order heuristic/certified/current/larger-search individual-Pauli controls and symmetric reuse of both heuristic and certified orders, at common 10U and 30U ideal Pauli-rotation caps. Exact central and boundary fusion determines each affordable step count. Symmetric orders were not optimized for a second-order objective.','',
        '| Formula | Rotation cap/U | Certified vs heuristic W/T/L | Floor exclusions | Median paired error ratio |',
        '| --- | ---: | ---: | ---: | ---: |']
    for s in summary:lines.append(f"| {s['formula']} | {s['budget_units']} | {s['wins']}/{s['ties']}/{s['losses']} | {s['floor_exclusions']} | {s['median_paired_error_ratio']} |")
    lines+=['',
        'Physical error is the squared projective state distance 1-|<exact|actual>|, evaluated stably from the orthogonal residual in the full Hilbert space. Number and joint (N_alpha,N_beta) leakage are retained; no sector projection is applied. Independent verification uses literal Pauli matrices, unfused circuits and a scaled-Taylor reference. Every floor and pairwise win/tie/loss decision is checked independently; all losses and floor-limited rows remain in the tables. The fixed floor and absolute comparison tolerance are 1e-12. Ratios exclude pairs with either error at/below the floor. Three related molecular cases are descriptive observations, not independent statistical replications.','',
        f"Independent replay: {verified['rows']} rows, maximum state difference {verified['max_state_delta']:.9g}, maximum metric difference {verified['max_metric_delta']:.9g}. Each integer certificate was independently re-enumerated by DFS with complete Lehmer-rank coverage, and its full-order objective was checked. Independent floating checks do not constitute an interval or real-transcendental certificate.",'',
        'The earlier four eight-unit windows on CH2O_s1.00 and nh3_coeff1p00 remain completed prior restricted evidence. Their retained original result and source are included under prior/ with hashes and provenance; they were not rerun. This study extends their domain size, without converting either study into a whole-order or physical-error optimum.','']
    (report/'README.md').write_text('\n'.join(lines))
    manifest=dict(results_schedule_sha256=sha(results/'SCHEDULES_FROZEN.json'),verified_sha256=sha(results/'VERIFIED.json'),
        files={p.name:sha(p) for p in sorted(report.iterdir()) if p.is_file()})
    (report/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();generate(a.results.resolve(),a.report.resolve())
