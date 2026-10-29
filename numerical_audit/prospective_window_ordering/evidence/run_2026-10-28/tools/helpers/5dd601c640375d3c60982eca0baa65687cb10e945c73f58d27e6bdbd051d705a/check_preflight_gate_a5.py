"""Resource/mechanics-only decision; never use physical error values."""
from pathlib import Path
import hashlib,json,math,re,tomllib
root=Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/prospective_window_ordering_2026-10-28/a1')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
gate=read(root/'PREFLIGHT_GATE_FROZEN.a5.json');plan=read(root/'preflight-plan.a5.json');pools=read(root/'RESOURCE_POOLS.json')
case='LiH_r1p487_631g_14m4e_preflight_7over5';prefix='window-preflight-'
assert len(plan['entries'])==8 and plan['attempt_tag']=='a5'
measurements={};calls=[]
def measurement(label):
 final=read(root/'jobs'/(label+'.FINAL.json'));receipt=read(root/'logs'/(label+'.resources.json'));log=(root/'logs'/(label+'.log')).read_text()
 assert final['exit_code']==receipt['exit_code']==0 and not receipt['termination'] and not receipt['guard_error'] and not receipt['cleanup_errors']
 peaks=[]
 for sample in receipt['samples']:
  for value in re.findall(r'MemoryPeak=(\d+)',sample.get('service','')):peaks.append(int(value)/1024**2)
 factors={'K':1/1024,'M':1,'G':1024,'B':1/1024**2}
 for number,unit in re.findall(r'^Memory peak: ([0-9.]+)([KMGB])$',log,re.M):
  peaks.append((float(number)+.1)*factors[unit]) # Conservative display-rounding upper margin.
 assert peaks,'Peak memory missing'
 peak=max(peaks)
 result=dict(label=label,charged_CPU_seconds=final['charged_CPU_seconds'],wall_seconds=final['absolute_guarded_wall_seconds'],owned_cgroup_peak_MiB_upper=peak,receipt_sha256=sha(root/'logs'/(label+'.resources.json')),final_accounting_sha256=sha(root/'jobs'/(label+'.FINAL.json')))
 measurements[label]=result
 return result
for entry in plan['entries']:
 m=measurement(entry['label']);calls.append(m)
 record=read(root/entry['result']);assert record['artifact_complete'] and record['physical_error_read'] is False
 if entry['arm'].startswith('window_'):
  stats=record['constructor']['stats']
  assert stats['child_checkpoint_status']=='complete' and stats['child_checked_updates']>=0
  assert stats['window_status'] in ('complete','native_no_resolved_improvement','verified_improvement')
  if stats['child_checked_updates']:assert stats['native_verifications']>=1
  assert stats['aggregation']=='exact_real_deltas_before_transport_v2'
for kind in ['alter','plan','freeze','panel','objectives','physical','independent']:measurement(prefix+kind+'-a5')
assert all(m['owned_cgroup_peak_MiB_upper']<=400 for m in measurements.values()),'Memory-margin gate failed'
objective=tomllib.loads((root/'preflight-objectives-a5.toml').read_text());physical=tomllib.loads((root/'preflight-physical-a5/PHYSICAL.toml').read_text());independent=tomllib.loads((root/'preflight-independent-a5.toml').read_text())
assert objective['complete'] and objective['mismatched_objectives']==0 and objective['floor_censored_objectives']==0
assert physical['status']=='complete' and physical['physical_panel_complete'] and len(physical['rows'])==44
assert independent['status']=='complete' and independent['checked_logical_rows']==independent['provided_logical_rows']==44 and not independent['failures']
assert len(independent['reference_checks'])==2 and independent['distinct_successful_unfused_replays']==independent['distinct_hard_circuit_identities']
# Only input size/bit metadata and receipt costs enter these projections.
inventory=read(root/'data-a6/INPUTS_FROZEN.json');models=[read(root/'data-a6/inputs'/name/'model.json') for name in inventory['successes']]
pref=read(root/'preflight-data-a5/inputs'/case/'model.json')
def bits(model):return max([int(model['denominator']).bit_length()]+[abs(int(value)).bit_length() for _,value in model['coefficients']])
word_scale=max(1,max(len(m['coefficients']) for m in models)/len(pref['coefficients']))
unit_scale=max(1,max(len(m['blocks']) for m in models)/len(pref['blocks']))
bit_scale=max(1,max(bits(m) for m in models)/bits(pref))**2
scale=8*2*word_scale*unit_scale*bit_scale # Fourcases x2reps,2x safety, count/bit scaling.
prior={phase:[0.,0.] for phase in pools['phases']}
for p in (root/'jobs').glob('*.FINAL.json'):
 f=read(p);prior[f['phase']][0]+=f['charged_CPU_seconds'];prior[f['phase']][1]+=f['absolute_guarded_wall_seconds']
projections={}
for phase,items in [('constructors',calls),('physical',[measurements[prefix+'physical-a5']]),('independent',[measurements[prefix+'independent-a5']])]:
 cpu=sum(x['charged_CPU_seconds'] for x in items)*scale;wall=sum(x['wall_seconds'] for x in items)*scale
 remaining=[pools['phases'][phase]['CPU_seconds']-prior[phase][0],pools['phases'][phase]['wall_seconds']-prior[phase][1]]
 assert cpu<=.8*remaining[0] and wall<=.8*remaining[1],phase+' projection gate failed'
 projections[phase]=dict(projected_CPU_seconds=cpu,projected_wall_seconds=wall,remaining_CPU_wall=remaining,ceiling_fraction=.8)
v=dict(schema='prospective-window-preflight-decision-v2',decision='passed_resource_and_mechanics_gate',source_commit='722013aa6d9a0a48776e26de5f356b5a441c3f61',gate_sha256=sha(root/'PREFLIGHT_GATE_FROZEN.a5.json'),quality_values_inspected_or_used=False,windows_completed=True,objective_replay_complete=True,physical_rows_checked=44,independent_reference_times_checked=2,distinct_hard_unfused_replays=independent['distinct_successful_unfused_replays'],measurements=measurements,scaling=dict(minimum_repetition_case_factor=8,additional_safety_factor=2,word_scale=word_scale,unit_scale=unit_scale,bit_work_scale=bit_scale,total_factor=scale),projections=projections,prior_phase_actual_CPU_wall=prior,all_prior_failures_retained=True,all_resource_pools_unchanged=True,main64_calls_352_rows='ADMITTED; not yet performed',checker_sha256=sha(Path(__file__)))
with (root/'PREFLIGHT_DECISION.a5.json').open('x') as f:json.dump(v,f,indent=2);f.write('\n')
print('PREFLIGHT_GATE_PASSED','max_cgroup_MiB_upper',max(x['owned_cgroup_peak_MiB_upper'] for x in measurements.values()),'scale',scale,'projections',projections)
