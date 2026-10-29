"""SOURCE-ONLY preparation; root must guard execution after the data-a6 audit."""
import copy,datetime,hashlib,json,os,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
ROOT=R/'runs/prospective_window_ordering_2026-10-28/a1'
WT=Path('/tmp/qhat-normal-window-2026-10-28/worktree')
F='a45e836458274a4ff2f37477951c2e92eabb218d'
G='722013aa6d9a0a48776e26de5f356b5a441c3f61'
O='87da247744ed58b7a96759fd8af3a44235b6cb38'
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  while block:=stream.read(1024**2):h.update(block)
 return h.hexdigest()
read=lambda path:json.loads(path.read_text())
desc=lambda name:dict(path=name,sha256=sha(ROOT/name))
def write(name,value):
 with (ROOT/name).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
def retag(value):
 if isinstance(value,dict):return {k:retag(v) for k,v in value.items()}
 if isinstance(value,list):return [retag(v) for v in value]
 if isinstance(value,str):return 'a5' if value=='a4' else value.replace('.a4','.a5').replace('-a4','-a5')
 return value
audit=read(ROOT/'INPUT_REUSE_AUDIT.a4.json')
assert audit['status']=='passed' and audit['payload_count']==16 and audit['new_execution_source_commit']==G
assert audit['original_generation_source_commit']==O and audit['chemistry_recomputed'] is False and audit['physical_quality_inspected'] is False
assert audit['source_snapshot_Git_bindings_verified']=={'data-a2':43,'data-a6':43}
env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
git=lambda *args:subprocess.check_output(['git','-C',str(WT),*args],env=env,text=True).strip()
assert git('rev-parse','HEAD')==G and not git('status','--porcelain=v1','--untracked-files=all')
prior=read(ROOT/'CAMPAIGN_FREEZE_B.a5.json');old=read(ROOT/'data-a2/INPUT_PROTOCOL_FROZEN.json')
protocol=read(ROOT/'data-a6/INPUT_PROTOCOL_FROZEN.json')
before=read(ROOT/'data-a2/INPUTS_FROZEN.json');after=read(ROOT/'data-a6/INPUTS_FROZEN.json')
reuse=read(ROOT/'data-a6/INPUT_REUSE_PROVENANCE.json')
assert old['source_commit']==O and protocol['source_commit']==G and len(protocol['sources'])==43
assert protocol['sources'].keys()==old['sources'].keys() and protocol['configuration']==old['configuration'] and protocol['environment']==old['environment']
for name,digest in protocol['sources'].items():
 path=(WT/name).resolve();assert path.is_relative_to(WT) and sha(path)==digest
 assert sha(ROOT/'data-a6/source_snapshot'/name)==digest
assert before['successes']==after['successes'] and before['status']==after['status']=='complete' and not after['refusals']
flat={}
for case,items in before['successes'].items():
 assert set(items)=={'model.json','model.toml','tensors.npz','orbitals.npz'}
 for name,digest in items.items():
  assert sha(ROOT/'data-a2/inputs'/case/name)==sha(ROOT/'data-a6/inputs'/case/name)==digest
  flat[case+'/'+name]=digest
assert len(flat)==16 and audit['original_payload_hashes']==flat
assert audit['new_protocol_sha256']==after['protocol_sha256']==sha(ROOT/'data-a6/INPUT_PROTOCOL_FROZEN.json')
assert audit['new_inventory_sha256']==sha(ROOT/'data-a6/INPUTS_FROZEN.json') and reuse['chemistry_rerun'] is False
assert reuse['original_payload_hashes']==before['successes'] and reuse['new_execution_source_commit']==G
assert prior['original_protocol']==desc('data-a2/INPUT_PROTOCOL_FROZEN.json') and prior['original_inventory']==desc('data-a2/INPUTS_FROZEN.json')
assert sha(ROOT/'RESOURCE_POOLS.json')==prior['resource_pools_sha256']
pools=read(ROOT/'RESOURCE_POOLS.json');config=protocol['configuration']
assert pools['resource_profile']['MemoryMax_MiB']==512 and config['process_CPU_ceiling']==64 and config['process_wall_ceiling']==320
declaration=desc('HEAP_GUIDANCE_CONTINUATION_PLAN.a1.json')
assert reuse['heap_guidance_declaration_sha256']==declaration['sha256']
checker4=R/'coordination/check_preflight_gate_a4.py';checker5=R/'coordination/check_preflight_gate_a5.py'
text=checker4.read_text();sentinel='__PARENT_DATA_A6_SENTINEL__';assert sentinel not in text
expected=text.replace('data-a5',sentinel).replace('.a4','.a5').replace('-a4','-a5').replace("'a4'","'a5'").replace('"a4"','"a5"')
expected=expected.replace(sentinel,'data-a6').replace(F,G)
assert checker5.read_text()==expected and G in expected and F not in expected and str(ROOT) in expected
names=['CAMPAIGN_FREEZE_B.a6.json','PREFLIGHT_GATE_FROZEN.a5.json','PREFLIGHT_STAGE_JOB_PLAN.a5.json']
assert all(not (ROOT/name).exists() for name in names)
assert not (ROOT/'preflight-data-a5').exists() and not (ROOT/'preflight-plan.a5.json').exists()
assert not list(ROOT.glob('preflight-*a5*'))
assert not (ROOT/'ALL_ORDERS_FROZEN.preflight.a5.json').exists() and not (ROOT/'PREFLIGHT_DECISION.a5.json').exists()
assert not list((ROOT/'jobs').glob('window-preflight-*-r1-a5.PLAN.json'))
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
b=copy.deepcopy(prior);b.update(source_commit=G,input_protocol=desc('data-a6/INPUT_PROTOCOL_FROZEN.json'),
 input_inventory=desc('data-a6/INPUTS_FROZEN.json'),reuse_provenance=desc('data-a6/INPUT_REUSE_PROVENANCE.json'),
 independent_reuse_audit=desc('INPUT_REUSE_AUDIT.a4.json'),prior_freeze_sha256=sha(ROOT/'CAMPAIGN_FREEZE_B.a5.json'),
 prior_B_freeze=desc('CAMPAIGN_FREEZE_B.a5.json'),heap_guidance_declaration=declaration,frozen_UTC=now)
gate4=read(ROOT/'PREFLIGHT_GATE_FROZEN.a4.json');gate=copy.deepcopy(gate4)
failed=read(ROOT/'PREFLIGHT_DECISION.a4.json');assert failed['decision']=='failed_memory_margin' and failed['main64_352_admitted'] is False
gate.update(input_protocol_sha256=b['input_protocol']['sha256'],input_inventory_sha256=b['input_inventory']['sha256'],
 source_commit=G,attempt_tag='a5',prior_gate_sha256=sha(ROOT/'PREFLIGHT_GATE_FROZEN.a4.json'),frozen_UTC=now,
 heap_guidance_declaration=declaration,preflight_gate_checker_sha256=sha(checker5),
 source_continuation='Normal child GC guidance64M to32M only; unchanged coefficients, criteria and policies; all physical quality sealed; retain every prior failure; no pool reset')
gate['prior_failed_decisions']['PREFLIGHT_DECISION.a4.json']=sha(ROOT/'PREFLIGHT_DECISION.a4.json')
assert gate['conditions']==gate4['conditions'] and gate['preflight_max_reservations_CPU_wall']==gate4['preflight_max_reservations_CPU_wall']
stage4=read(ROOT/'PREFLIGHT_STAGE_JOB_PLAN.a4.json');stage=retag(stage4)
stage.update(scientific_source=G,frozen_UTC=now)
assert stage['root']==str(ROOT) and stage['cwd']==str(WT) and len(stage['jobs'])==5
assert [(j['phase'],j['cpu'],j['wall']) for j in stage['jobs']]==[(j['phase'],j['cpu'],j['wall']) for j in stage4['jobs']]
for name,value in zip(names,(b,gate,stage)):write(name,value)
assert sha(ROOT/'RESOURCE_POOLS.json')==prior['resource_pools_sha256']
print('HEAP32_PREFLIGHT_METADATA_FROZEN','source',G,'attempt','a5','original_payloads',16,'source_files',43,'MAIN_NOT_ADMITTED')
