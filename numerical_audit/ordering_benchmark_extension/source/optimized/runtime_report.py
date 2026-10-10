"""Describe measured classical work with explicit, heterogeneous timer scopes."""
from pathlib import Path
import json,statistics,sys,platform
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
rows=json.loads((HERE/'results.json').read_text());c.require(len(rows)==728,'Incomplete timing campaign')
designs=[json.loads((HERE/'designs'/f'{case}.design.json').read_text()) for case in c.CASES]
methods=list(designs[0]['optimizer_reports']);search={};setups={}
for method in methods:
    part=[(d['case'],d['optimizer_reports'][method]) for d in designs];times=[r['search_seconds'] for _,r in part];worst=max(part,key=lambda t:t[1]['search_seconds'])
    search[method]={'cases':len(part),'evaluations':sum(r['evaluations'] for _,r in part),'rescored_contenders':sum(r['rescored'] for _,r in part),'search_seconds_sum':sum(times),'search_seconds_median_per_case':statistics.median(times),'largest_search_case':worst[0],'largest_search_seconds':worst[1]['search_seconds']}
for name,method,users in (('refined_conserving','class_refined_response',['class_refined_response','class_refined_J0']),('individual_pauli','class_pauli_JD',['class_pauli_JD','class_pauli_JD_budget'])):
    values=[d['optimizer_reports'][method]['pair_setup_seconds'] for d in designs]
    for d in designs:c.require(all(d['optimizer_reports'][u]['pair_setup_seconds']==d['optimizer_reports'][method]['pair_setup_seconds'] for u in users),'Shared setup record differs')
    setups[name]={'shared_by':users,'actual_setup_invocations':14,'setup_seconds_sum_counted_once':sum(values),'setup_seconds_median_per_case':statistics.median(values),'setup_seconds_max':max(values)}
evolution={}
for name in designs[0]['methods']:
    rr=[r for r in rows if r['ordering']==name];tt=[r['seconds'] for r in rr]
    evolution[name]={'rows':len(rr),'actual_pauli_rotations':sum(r['pauli_rotations'] for r in rr),'recorded_row_seconds_sum':sum(tt),'recorded_row_seconds_median':statistics.median(tt),'scope':'Propagation plus pulse construction, state compression/write, and design hash; excludes metrics, state hash, row writing and setup.'}
resource=[json.loads(p.read_text()) for p in sorted((HERE/'execution_status').glob('*.json'))]
ratio=[]
for d in designs:
    a=d['optimizer_reports']['class_refined_response'];b=d['optimizer_reports']['class_pauli_JD'];ratio.append({'case':d['case'],'unrestricted_pauli_over_conserving_evaluation_ratio':b['evaluations']/a['evaluations'],'unrestricted_pauli_over_conserving_search_seconds_ratio':b['search_seconds']/a['search_seconds'],'matched_budget_evaluations_exactly_equal':d['optimizer_reports']['class_pauli_JD_budget']['evaluations']==a['evaluations']})
out={'search':search,'pair_setup_counted_once':setups,'evolution':evolution,'case_ratios':ratio,'preparation_wall_seconds_sum':sum(d['preparation_seconds'] for d in designs),'recorded_propagation_row_seconds_sum':sum(r['seconds'] for r in rows),'recorded_rotations':sum(r['pauli_rotations'] for r in rows),'largest_recorded_propagation_process_peak_rss_KiB':max(r['process_peak_rss_KiB'] for r in resource),'python':platform.python_version(),'threads_per_process':1,'per_process_address_space_cap_GiB':6,'timing_scope_sha256':c.sha(HERE/'timing_scope.json'),'source_sha256':c.sha(__file__),'caveat':'Wall times are descriptive, with concurrent jobs and variable cache/JIT state. Row intervals include I/O and cannot be compared to fixed-campaign evolve-only timers. Preparation/search totals are sums across cases, not elapsed user wait time. Shared pair setup counted once; an interrupted unfinished row is outside completed-row timers.'}
c.write_atomic(HERE/'runtime_report.json',out)
lines=['# Runtime and work of the optimized campaign','','All numbers are observations on this machine, with one thread per process and a 6 GiB address-space cap. Multiple independent jobs ran concurrently. Counts are reproducible algorithmic work measures; wall times are descriptive.','','## Search across all 14 cases','','| Arm | Examined insertion destinations | Total search seconds | Median seconds/case | Largest case seconds |','|---|---:|---:|---:|---:|']
for method,r in search.items():lines.append(f"| {method} | {r['evaluations']:,} | {r['search_seconds_sum']:.3f} | {r['search_seconds_median_per_case']:.3f} | {r['largest_search_case']}: {r['largest_search_seconds']:.3f} |")
lines+=['','The conserving JD and J0 arms share one pair cache per case; the two individual-Pauli arms share another. Counting each observed setup once:','']
for name,r in setups.items():lines.append(f"- {name}: {r['setup_seconds_sum_counted_once']:.3f}s total over 14 setups; median {r['setup_seconds_median_per_case']:.3f}s/case.")
lines+=['','Search timers include initial pair-amplitude reconstruction and between-sweep checkpoint overhead, but exclude pair setup and starting-order construction. The matched Pauli arm has exactly the conserving JD destination count on every case; it is not an equal-wall-time control.','','## Circuit measurement intervals','','| Arm | States | Pauli rotations | Sum of recorded seconds |','|---|---:|---:|---:|']
for name,r in evolution.items():lines.append(f"| {name} | {r['rows']} | {r['actual_pauli_rotations']:,} | {r['recorded_row_seconds_sum']:.3f} |")
lines+=['',f"The 728 saved rows emit {out['recorded_rotations']:,} rotations in total. Their recorded intervals sum to {out['recorded_propagation_row_seconds_sum']:.3f}s; this is not wall-clock elapsed time because jobs ran concurrently.",'','These intervals include pulse construction, state compression/output, and a design hash as well as evolution. They exclude error metrics and row writing. The fixed campaign uses evolve-only timers, so the two columns must not be treated as identical runtime measurements. See `timing_scope.json` for exact boundaries. The deliberate NH3 process split reused completed rows; work in an interrupted unfinished row is not in the totals.',f"The largest recorded propagation-process peak is {out['largest_recorded_propagation_process_peak_rss_KiB']/1024:.1f} MiB. This is a whole-process peak, possibly including earlier cases, and is not an optimizer-memory measurement.",'']
(HERE/'RUNTIME_REPORT.md').write_text('\n'.join(lines));print(json.dumps({'search':search,'setup':setups,'row_seconds_sum':out['recorded_propagation_row_seconds_sum'],'rotations':out['recorded_rotations']},indent=2))
