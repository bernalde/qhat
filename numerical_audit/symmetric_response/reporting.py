"""Retain all symmetric-specific outcomes, controls and construction costs."""
import argparse,csv,hashlib,itertools,json,statistics
from pathlib import Path


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path,rows):
    require(bool(rows),'Empty report table')
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)


def generate(results,report):
    require(not report.exists(),'Fresh report required');report.mkdir(parents=True)
    read=lambda p:json.loads(p.read_text());verified=read(results/'VERIFIED.json');protocol=read(results/'PROTOCOL_FROZEN.json')
    require(verified['status']=='passed' and verified['decision_agreement'],'Unverified study')
    require(verified['schedule_freeze_sha256']==sha(results/'SCHEDULES_FROZEN.json'),'Different verified schedule freeze')
    for name,digest in verified['artifacts'].items():require(sha(results/name)==digest,'Verified artifact changed')
    physical=[];shared=[];searches=[];seeds=[];phase_costs=[]
    for path in sorted((results/'designs').glob('*.json')):
        design=read(path);case=design['case'];shared.append(dict(case=case,**design['construction']))
        executions=[]
        for tr in design['time_records']:
            T=tr['T'];executions += [('first_JD_sameT',T,10,tr['first']),('pauli_JD_large_sameT',T,10,tr['pauli'])]
            executions += [('cubic_insertion',T,cell['budget_units'],cell['cubic']) for cell in tr['cells']]
        for name,T,budget,record in executions:
            result=record['result']
            searches.append(dict(case=case,method=name,T=T,optimization_rotation_cap_units=budget,
                nominal_search_CPU_seconds=result['nominal_search_CPU_seconds'],actual_search_CPU_seconds=record['phases']['search_CPU_seconds'],
                search_CPU_overrun_seconds=max(0.,record['phases']['search_CPU_seconds']-result['nominal_search_CPU_seconds']),
                actual_local_constructor_CPU_seconds=record['actual_local_CPU_seconds'],
                attempted=result['attempted'],completed=result['completed'],evaluation_cap=result['evaluation_cap'],
                accepted_moves=sum(len(s['accepted']) for s in result['starts']),
                partial_sweeps=sum(not w['complete'] for s in result['starts'] for w in s['sweeps']),
                completed_rescores=result['accepted_rescores'],attempted_rescores=result['attempted_rescores'],
                termination=json.dumps([s['termination'] for s in result['starts']],separators=(',',':')),
                selected_score=result['record']['score'],selected_steps=result['record']['steps']))
            phase_costs.append(dict(case=case,method=name,T=T,optimization_rotation_cap_units=budget,
                actual_local_constructor_CPU_seconds=record['actual_local_CPU_seconds'],
                bookkeeping_CPU_seconds=record['bookkeeping_CPU_seconds'],cost_timing=record['cost_timing'],**record['phases']))
            for index,seed in enumerate(record['seeds']):
                seeds.append(dict(case=case,method=name,T=T,optimization_rotation_cap_units=budget,seed_index=index,
                    JD1=seed.get('JD1'),JD2=seed.get('JD2'),score=seed['score'],steps=seed['steps'],order=json.dumps(seed['order'],separators=(',',':'))))
    for path in sorted((results/'rows').glob('*.json')):
        row=read(path);flat={k:v for k,v in row.items() if not isinstance(v,dict)}
        for key,value in row.items():
            if isinstance(value,dict):flat.update({key+'_'+k:v for k,v in value.items()})
        physical.append(flat)
    physical.sort(key=lambda row:(row['case'],row['record_index']))
    require(len(physical)==protocol['expected_circuits'],'Incomplete physical grid')
    comparisons=[];classify=lambda d:'win' if d < -1e-12 else 'loss' if d>1e-12 else 'tie'
    for at,a in enumerate(physical):
        for b in physical[at+1:]:
            if (a['case'],a['T'],a['budget_units'])!=(b['case'],b['T'],b['budget_units']):continue
            difference=a['error']-b['error'];floor=min(a['error'],b['error'])<=1e-12
            comparisons.append(dict(case=a['case'],T=a['T'],budget_units=a['budget_units'],
                left_record=a['record_index'],right_record=b['record_index'],left_arm=a['arm'],right_arm=b['arm'],
                left_formula=a['formula'],right_formula=b['formula'],left_steps=a['steps'],right_steps=b['steps'],
                left_error=a['error'],right_error=b['error'],error_difference=difference,result=classify(difference),
                same_formula=a['formula']==b['formula'],floor_excluded_ratio=floor,error_ratio=None if floor else a['error']/b['error'],
                left_constructor_CPU=a['standalone_construction_CPU_charge'],right_constructor_CPU=b['standalone_construction_CPU_charge']))
    require(len(comparisons)==28*len(physical)//8,'Incomplete all-pairs comparisons')
    controls=[('current',1),('first_JD_sameT',1),('pauli_JD_large_sameT',1),('current',2),('first_JD_sameT',2),('pauli_JD_large_sameT',2),('cubic_orientation',2)]
    summary=[];cases=sorted({r['case'] for r in physical})
    for T,budget,(control,formula) in itertools.product(protocol['times'],protocol['budgets'],controls):
        pairs=[]
        for case in cases:
            a=next(r for r in physical if (r['case'],r['T'],r['budget_units'],r['arm'])==(case,T,budget,'cubic_insertion'))
            b=next(r for r in physical if (r['case'],r['T'],r['budget_units'],r['arm'],r['formula'])==(case,T,budget,control,formula))
            pairs.append((a,b))
        deltas=[a['error']-b['error'] for a,b in pairs];ratios=[a['error']/b['error'] for a,b in pairs if min(a['error'],b['error'])>1e-12]
        summary.append(dict(T=T,budget_units=budget,control=control,control_formula=formula,pairs=len(pairs),
            wins=sum(d < -1e-12 for d in deltas),ties=sum(abs(d)<=1e-12 for d in deltas),losses=sum(d>1e-12 for d in deltas),
            floor_exclusions=len(pairs)-len(ratios),median_of_paired_ratios=statistics.median(ratios) if ratios else None,
            median_of_paired_constructor_CPU_ratios=statistics.median(a['standalone_construction_CPU_charge']/b['standalone_construction_CPU_charge'] for a,b in pairs)))
    for name,rows in [('all_rows.csv',physical),('comparisons.csv',comparisons),('summary.csv',summary),
        ('shared_construction.csv',shared),('searches.csv',searches),('constructor_phases.csv',phase_costs),('seeds.csv',seeds)]:table(report/name,rows)
    lines=['# Symmetric-formula-specific response ordering','',
        f"This {'altered-target preflight' if protocol['preflight'] else 'fixed prospective panel'} retains {len(physical)} full-space circuits on {len(cases)} retained molecular targets. Times are {protocol['times']} inverse Hartree, with both10U and30U ideal Pauli-rotation caps. All eight schedules per cell were frozen before new physical errors. No input or loss is removed.",'',
        'The new objective uses the actual symmetric-step cubic defect column A3|HF>, evaluated by a degree-three Taylor recurrence. A separate literal-Pauli nested-commutator construction checks the full signed column. The diagonal-reference score is J_D2=sum(w!=HF) T^2 sinc^2(T gap_w/2) A3_w^2; optimization uses J_D2/r^4, with r the exact affordable fused symmetric step count. A positive common T^4 factor is omitted within each fixed-time cell. This is a formula-specific approximation, not an exact physical-error objective or a global optimum certificate. See DERIVATION.md and the supplied coefficients/protocol.','',
        'Three seeds are freshly optimized same-time first-order J_D, its reverse, and legacy current placement. The orientation-only control selects the lowest new score among these three; bounded insertion performs additional optimization of that new objective. Search has three equal CPU slices without unused carry, a2CPU-second nominal search allowance,1000attempted positions and at most5moves/start. All incomplete candidates/sweeps, partial improvements, rescoring and actual work remain in the frozen raw histories. Constructor setup and seed scoring are separately measured.','',
        'Strong controls include fresh same-time conserving first-order J_D (2CPU-second nominal search/1000positions/5moves per start), fresh same-time individual-Pauli J_D with a larger8CPU-second search/20m(m-1)positions/10moves per start, and separately named legacy placement whose selection reference remains T=1. Both formulas reuse all three controls. First-order orders are reused across rotation caps because their affordable first-order count is order independent. Symmetric boundary/central fusion and each actual affordable step count are retained. Formula-changing comparisons do not isolate ordering.','',
        'New-method standalone charges include upstream first-order initialization, common model/base work, current placement, cubic cache setup/cleanup, new seeds, search and complete constructor returns. These are not matched-total-CPU methods. Nominal search allowances and actual overruns are separate from complete construction charges. Shared construction and raw process resources remain distinct; standalone costs repeated across physical rows are not a sum of actual campaign work. Current placement and cache work are counted once where executed. No practical or asymptotic speedup claim follows from this small fixed-register panel.','',
        '| T | Rotation cap/U | Cubic insertion versus | Control formula | W/T/L | Floor exclusions | Median paired error ratio |',
        '| ---: | ---: | --- | ---: | ---: | ---: | ---: |']
    for row in summary:lines.append(f"| {row['T']} | {row['budget_units']} | {row['control']} | {row['control_formula']} | {row['wins']}/{row['ties']}/{row['losses']} | {row['floor_exclusions']} | {row['median_of_paired_ratios']} |")
    lines+=['',
        'All absolute projective state errors1-|overlap|,infidelities,norms,particle-number and joint(N_alpha,N_beta)leakage,actual rotations/steps and constructor costs appear in all_rows.csv. Full2^n propagation applies no symmetry-sector projection. The fixed error floor and absolute W/T/L tolerance are1e-12; ratios exclude either side at/below the floor. Paired differences/ratios are formed before descriptive aggregation over three related molecular cases. No independent statistical-replication claim is made.','',
        f"Independent replay checked every circuit, all {len(comparisons)} pairwise decisions and all floor exclusions. Maximum state difference {verified['max_state_delta']:.9g}; maximum metric difference {verified['max_metric_delta']:.9g}; full cubic-column difference {verified['max_cubic_column_delta']:.9g}. Reference propagation uses a separate scaled-Taylor approximation; these are numerical checks, not interval or kernel proofs.",'',
        'The source/protocol/input/order hashes and complete circuits/states are retained. A failed two-level correction remains rejected and was not retuned. This study supplies no conservation-mechanism identification, exact-Hamiltonian search oracle, compiled gate/depth result, hardware measurement, larger-molecule validation or whole-order optimum.','']
    (report/'README.md').write_text('\n'.join(lines))
    manifest=dict(verified_sha256=sha(results/'VERIFIED.json'),schedule_freeze_sha256=sha(results/'SCHEDULES_FROZEN.json'),
        files={p.name:sha(p) for p in sorted(report.iterdir()) if p.is_file()})
    (report/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();generate(a.results.resolve(),a.report.resolve())
