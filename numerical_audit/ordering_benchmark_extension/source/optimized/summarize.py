"""Descriptive comparison of frozen optimized and second-order controls."""
from pathlib import Path
import json,sys,statistics
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
rows=json.loads((HERE/'results.json').read_text());c.require(len(rows)==728,'Incomplete campaign');floor=1e-12
methods=list(json.loads((HERE/'designs'/f'{c.CASES[0]}.design.json').read_text())['methods']);lookup={(r['case'],r['ordering'],r['budget_units']):r for r in rows};out={'rows':len(rows),'metric':'E=1-|<exact,approx>|/(||exact|| ||approx||), evaluated stably','comparison':'candidate divided by comparator at the same requested budget100U','absolute_tie_and_ratio_denominator_floor':floor,'panels':{},'number_leakage':{},'classical_work':[],'threshold_work':{}}
for panel,cases in (('development',c.DEV_CASES),('prospective',c.PROS_CASES),('all',c.CASES)):
    report={}
    for method in methods:
        if method=='class_refined_response':continue
        counts={'win':0,'tie':0,'loss':0};ratios=[]
        for case in cases:
            candidate=lookup[case,'class_refined_response',100]['one_minus_overlap'];ref=lookup[case,method,100]['one_minus_overlap'];diff=candidate-ref
            counts['tie' if abs(diff)<=floor else 'win' if diff<0 else 'loss']+=1
            if ref>floor:ratios.append(candidate/ref)
        report[method]={'candidate_outcomes':counts,'ratios_available':len(ratios),'median_candidate_over_comparator':statistics.median(ratios) if ratios else None}
    out['panels'][panel]=report
for method in methods:out['number_leakage'][method]={'maximum_all_budgets':max(r['particle_number_leakage'] for r in rows if r['ordering']==method),'maximum_at_budget100U':max(lookup[case,method,100]['particle_number_leakage'] for case in c.CASES)}
for case in c.CASES:
    d=json.loads((HERE/'designs'/f'{case}.design.json').read_text());out['classical_work'].append({'case':case,'methods':{name:{k:r[k] for k in ('units','nonzero_pair_count','pair_setup_seconds','search_seconds','evaluations','evaluation_budget','rescored')} for name,r in d['optimizer_reports'].items()}})
    out['threshold_work'][case]={}
    for method in methods:
        candidate=[r for r in rows if r['case']==case and r['ordering']==method]
        out['threshold_work'][case][method]={str(e):min((r['pauli_rotations'] for r in candidate if r['one_minus_overlap']<=e),default=None) for e in (1e-4,1e-6,1e-8)}
out['source_sha256']={'summarize.py':c.sha(__file__),'results.json':c.sha(HERE/'results.json'),'ALL_SCHEDULES_FROZEN.json':c.sha(HERE/'ALL_SCHEDULES_FROZEN.json'),'EXECUTION_MANIFEST.json':c.sha(HERE/'EXECUTION_MANIFEST.json')}
c.write_atomic(HERE/'summary.json',out)
text=['# Optimized and Strang controls: recovered campaign','',f'Completed {len(rows)} propagations from 14 frozen molecular targets. All comparisons below use a requested budget of 100U Pauli rotations; Strang uses the largest feasible integer step count. The eight development geometries and six prospective coefficient-scaled targets are separate panels, not fourteen independent samples.','','The candidate is the conserving refined JD first-order ordering. W/T/L compares its error against each comparator with an absolute tolerance of 1e-12. Ratios whose comparator error is at or below that floor are omitted. This is a comparison within this campaign, not a conclusion against all published/fixed controls.','','| Comparator | Development W/T/L | Prospective W/T/L | All-case median candidate/comparator | Ratios |','|---|---:|---:|---:|---:|']
for method in methods:
    if method=='class_refined_response':continue
    report=out['panels']['all'][method];parts=[]
    for panel in ('development','prospective'):
        r=out['panels'][panel][method]['candidate_outcomes'];parts.append('/'.join(str(r[k]) for k in ('win','tie','loss')))
    ratio=report['median_candidate_over_comparator'];text.append('| '+method+' | '+' | '.join(parts)+' | '+('omitted' if ratio is None else f'{ratio:.6g}')+f" | {report['ratios_available']} |")
text+=['','The matched Pauli arm controls insertion-destination evaluations, not wall-clock time or setup. Complete setup/search/evaluation counts, leakage, and least tested actual rotations at errors1e-4,1e-6,1e-8 are in `summary.json`. Strang changes formula order; its comparison is not evidence for or against ordering alone. The objective, pulse-fusion, checkpoint-provenance, and literal archived-schedule checks are separate artifacts. Independent replay/objective checks and complete cross-method rankings belong to the enclosing campaign.','']
(HERE/'REPORT.md').write_text('\n'.join(text));print(json.dumps(out['panels'],indent=2),flush=True)
