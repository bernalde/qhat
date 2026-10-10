#!/usr/bin/env python3
"""Describe measured timer boundaries without changing frozen scientific code."""
import datetime,json,os,platform,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import sha,write_atomic
cpuinfo=Path('/proc/cpuinfo').read_text();model=next((x.split(':',1)[1].strip() for x in cpuinfo.splitlines() if x.startswith('model name')),'unknown')
mem={k:int(v.split()[0])*1024 for k,v in (line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines()) if k in ('MemTotal','MemAvailable','SwapTotal','SwapFree')}
record={
 'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'source_sha256':sha(HERE/'run.py'),
 'interpretation':'Descriptive observed timings for this recovery campaign, not isolated microbenchmarks or proof of a speed advantage over other campaigns or implementations.',
 'row_seconds':{
  'start':'Immediately before evolve(initial, compiled_terms, T=1, steps)',
  'stop':'Immediately after evolve returns; before evaluating physical_metrics or writing any NPZ/JSON',
  'includes':['QHAT propagation routine overhead','Any phase/permutation cache setup performed inside evolve','First invocation in each worker may include Numba compilation/cache startup'],
  'excludes':['Frozen-input and design loading','Affine-basis compilation in compile_case','Converting the current exact schedule to compiled terms','Error/leakage metric computation','NPZ state compression/write','JSON result write','Schedule search and graph construction'],
  'comparison_warning':'Do not pool directly with optimized campaign row.seconds, whose timer includes state NPZ writing.'
 },
 'preparation_seconds':{
  'field':'design.extension_preparation.seconds',
  'start':'Before normalized_case input loading',
  'stop':'After all18 base,60 random and4 published schedules and target checks; before writing final verbose design JSON',
  'includes':['Input/design/state loading','Basic Pauli coloring','New order construction','Exact target-sum validation','Random permutations','Published-control schedule generation','Intermediate structural-cache JSON writes where applicable'],
  'excludes':['Writing the final verbose case design JSON','Whole-campaign source/design freeze hashing','State propagation'],
  'development_reuse':'Existing exact parent/owned/class graph color orders are read from the already-reviewed identical development designs; class commutators/suffix/BCH orders are rebuilt and matched exactly to archived sequences.',
  'positive_scale_reuse':'Parent, owned and class graph color orders and class suffix/BCH orders are reused for a coefficient-scaled case only after equality of a fingerprint containing all parent expansions and the consolidated target divided by their common integer gcd. This proves positive proportionality of the exact ordered inputs.',
  'graph_seconds_warning':'Development graph timing fields do not represent rebuilding all old graphs. A scaled cached case has no fresh graph timer. These are replay/preparation costs, not uncached algorithm complexity measurements.',
  'published_seconds':'Elapsed call to build_published_schedules; excludes adding/checking its returned schedules only insofar as the timer actually surrounds the whole loop (here validation/insertion IS included).'
 },
 'resources':{'per_worker_RLIMIT_AS_bytes':6*1024**3,'numerical_threads_per_worker':1,'fixed_workers':3,'partition':'One worker for12 non-NH3 cases; one worker for each NH3 coefficient case. The small worker has completed while the two NH3 workers continue.','other_activity':'Concurrent optimized-order search, exact native Tranter V2 preparation, published extra propagation and independent verification were running during parts of this campaign. CPU contention and shared filesystem/cache effects were not controlled.','python':sys.version,'platform':platform.platform(),'machine':platform.machine(),'processor_model':model,'logical_cpus_reported':os.cpu_count(),'memory_at_observation_bytes':mem,'memory_note':'Available memory is a contemporaneous observation, not peak resident memory. Prior live observation recorded~0.5GiB smallworker and~1.0GiB each NH3 worker; these are also observations, not measured peaks.'},
 'job_records':'job_*.json provide per-worker wall-clock start/finish times; those include imports, input loading, verification, metrics, compression, saving and checkpoint resume overhead. Their intervals overlap, so summing elapsed worker times is not total campaign wall time.',
 'pre_execution_corrections':'The exact archival owned-word alignment and reserialization after initial preparation are documented separately in baseline_alignment.json and are not retroactively added to stored preparation timers.'
}
# State the actual published timer boundary directly, without an ambiguous shorthand.
record['preparation_seconds']['published_seconds']='Elapsed loop invoking build_published_schedules and adding all4 returned schedules, including exact target-sum validation and method metadata construction; excludes writing the final design JSON.'
write_atomic(HERE/'timing_scope.json',record);print('Wrote timing_scope.json')
