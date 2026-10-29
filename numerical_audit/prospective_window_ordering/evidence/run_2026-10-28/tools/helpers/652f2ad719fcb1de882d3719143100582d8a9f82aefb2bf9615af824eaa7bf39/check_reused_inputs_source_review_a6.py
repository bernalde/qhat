"""Independent stdlib/Git authentication of chemistry reuse and source binding."""
from pathlib import Path
import hashlib,json,subprocess
root=Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/prospective_window_ordering_2026-10-28/a1');wt=Path('/tmp/qhat-normal-window-2026-10-28/worktree')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
old=root/'data-a2';new=root/'data-a6';op=read(old/'INPUT_PROTOCOL_FROZEN.json');np=read(new/'INPUT_PROTOCOL_FROZEN.json');of=read(old/'INPUTS_FROZEN.json');nf=read(new/'INPUTS_FROZEN.json');reuse=read(new/'INPUT_REUSE_PROVENANCE.json')
assert of['protocol_sha256']==sha(old/'INPUT_PROTOCOL_FROZEN.json') and nf['protocol_sha256']==sha(new/'INPUT_PROTOCOL_FROZEN.json')
assert of['status']==nf['status']=='complete' and of['expected_cases']==nf['expected_cases']==4 and not of['refusals'] and not nf['refusals']
assert of['successes']==nf['successes']==reuse['original_payload_hashes']
assert op['registry']==np['registry'] and op['configuration']==np['configuration'] and op['environment']==np['environment']
assert reuse['original_generation_source_commit']==op['source_commit'] and reuse['new_execution_source_commit']==np['source_commit']
assert reuse['original_protocol_sha256']==sha(old/'INPUT_PROTOCOL_FROZEN.json') and reuse['original_inventory_sha256']==sha(old/'INPUTS_FROZEN.json') and reuse['new_protocol_sha256']==sha(new/'INPUT_PROTOCOL_FROZEN.json') and reuse['chemistry_rerun'] is False
assert set(of['successes'])=={spec['case'] for spec in op['registry']['cases']}
payloads={}
for case,files in of['successes'].items():
 assert set(files)=={'model.json','model.toml','tensors.npz','orbitals.npz'}
 for name,expected in files.items():
  assert sha(old/'inputs'/case/name)==sha(new/'inputs'/case/name)==expected
  payloads[case+'/'+name]=expected
declaration=read(root/'SOURCE_REVIEW_CONTINUATION_PLAN.a1.json')
allowed=set(declaration['allowed_changed_scientific_paths'])
addendum=read(root/'SOURCE_REVIEW_CONTINUATION_ADDENDUM.a1.json')
allowed.update(addendum['additional_allowed_changed_scientific_paths'])
changed={name for name in set(op['sources'])|set(np['sources']) if op['sources'].get(name)!=np['sources'].get(name)}
assert changed and changed<=allowed and set(reuse['changed_sources'])==changed
for name,expected in declaration['protected_generator_conversion_registry_files'].items():
 assert op['sources'].get(name)==np['sources'].get(name)==expected
assert reuse['protected_generator_conversion_registry_files']==declaration['protected_generator_conversion_registry_files']
assert reuse['source_review_declaration_sha256']==sha(root/'SOURCE_REVIEW_CONTINUATION_PLAN.a1.json') and reuse['source_review_addendum_sha256']==sha(root/'SOURCE_REVIEW_CONTINUATION_ADDENDUM.a1.json')
assert reuse['all_chemistry_generation_conversion_registry_environment_configuration_sources_unchanged'] is True
assert reuse['heap_guidance_declaration_sha256']==sha(root/'HEAP_GUIDANCE_CONTINUATION_PLAN.a1.json')
counts={}
for directory,protocol in [(old,op),(new,np)]:
 for name,expected in protocol['sources'].items():
  assert sha(directory/'source_snapshot'/name)==expected
  blob=subprocess.check_output(['git','-C',str(wt),'show',protocol['source_commit']+':'+name])
  assert hashlib.sha256(blob).hexdigest()==expected
 counts[directory.name]=len(protocol['sources'])
out=dict(schema='independent-chemistry-reuse-audit-v1',status='passed',original_generation_source_commit=op['source_commit'],new_execution_source_commit=np['source_commit'],payload_count=len(payloads),original_payload_hashes=payloads,source_snapshot_Git_bindings_verified=counts,changed_sources=sorted(changed),original_protocol_sha256=sha(old/'INPUT_PROTOCOL_FROZEN.json'),original_inventory_sha256=sha(old/'INPUTS_FROZEN.json'),new_protocol_sha256=sha(new/'INPUT_PROTOCOL_FROZEN.json'),new_inventory_sha256=sha(new/'INPUTS_FROZEN.json'),reuse_sidecar_sha256=sha(new/'INPUT_REUSE_PROVENANCE.json'),chemistry_recomputed=False,physical_quality_inspected=False,checker_sha256=sha(Path(__file__)))
with (root/'INPUT_REUSE_AUDIT.a4.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print('REUSE_AUDIT_PASSED',len(payloads),'payloads',counts,'sourcechanges',len(changed))
