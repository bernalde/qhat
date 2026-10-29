"""Rebind byte-identical chemistry inputs after a declared runtime-only repair."""
import hashlib,json,shutil,sys
from pathlib import Path
study=Path('/tmp/qhat-normal-window-2026-10-28/worktree/numerical_audit/prospective_window_ordering');sys.path.insert(0,str(study))
import generate_inputs as g
k=g.k
root=Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/prospective_window_ordering_2026-10-28/a1')
old=root/'data-a2';new=root/'data-a6'
prior=k.read(old/'INPUT_PROTOCOL_FROZEN.json');sources=k.source_map()
declaration=k.read(root/'SOURCE_REVIEW_CONTINUATION_PLAN.a1.json')
allowed=set(declaration['allowed_changed_scientific_paths'])
addendum=k.read(root/'SOURCE_REVIEW_CONTINUATION_ADDENDUM.a1.json')
allowed.update(addendum['additional_allowed_changed_scientific_paths'])
changed={n for n in set(prior['sources'])|set(sources) if prior['sources'].get(n)!=sources.get(n)}
k.require(changed and changed<=allowed,'Only explicitly declared reviewed runtime repair paths may change')
for name,expected in declaration['protected_generator_conversion_registry_files'].items():
 k.require(prior['sources'].get(name)==sources.get(name)==expected,'Protected chemistry-generation dependency changed: '+name)
k.require(prior['configuration']==k.configuration() and prior['environment']==k.environment(),'No configuration/environment retune')
k.require(k.read(k.HERE/'REGISTRY.json')==prior['registry'],'No scientific-instance substitution')
oldfreeze=k.read(old/'INPUTS_FROZEN.json')
k.require(oldfreeze['status']=='complete' and len(oldfreeze['successes'])==4 and not oldfreeze['refusals'],'All old cases required')
k.require(oldfreeze['protocol_sha256']==k.sha(old/'INPUT_PROTOCOL_FROZEN.json'),'Prior inputfreeze/protocol linkage differs')
k.require(set(oldfreeze['successes'])=={row['case'] for row in prior['registry']['cases']},'Exact original case inventory required')
k.require(all(set(files)=={'model.json','model.toml','tensors.npz','orbitals.npz'} for files in oldfreeze['successes'].values()),'Exact original fourfile inventories required')
for case,items in oldfreeze['successes'].items():
 for name,digest in items.items():k.require(k.sha(old/'inputs'/case/name)==digest,'Prior chemistry bytes differ')
g.prepare(new)
shutil.copytree(old/'inputs',new/'inputs')
for case,items in oldfreeze['successes'].items():
 for name,digest in items.items():k.require(k.sha(new/'inputs'/case/name)==digest,'Reused chemistry bytes differ')
k.write_new(new/'INPUT_REUSE_PROVENANCE.json',dict(schema='prospective-chemistry-reuse-v1',original_generation_source_commit=prior['source_commit'],new_execution_source_commit=k.committed_sources(),original_protocol_sha256=k.sha(old/'INPUT_PROTOCOL_FROZEN.json'),original_inventory_sha256=k.sha(old/'INPUTS_FROZEN.json'),new_protocol_sha256=k.sha(new/'INPUT_PROTOCOL_FROZEN.json'),changed_sources=sorted(changed),all_chemistry_generation_conversion_registry_environment_configuration_sources_unchanged=True,chemistry_rerun=False,original_payload_hashes=oldfreeze['successes'],protected_generator_conversion_registry_files=declaration['protected_generator_conversion_registry_files'],source_review_declaration_sha256=k.sha(root/'SOURCE_REVIEW_CONTINUATION_PLAN.a1.json'),source_review_addendum_sha256=k.sha(root/'SOURCE_REVIEW_CONTINUATION_ADDENDUM.a1.json'),heap_guidance_declaration_sha256=k.sha(root/'HEAP_GUIDANCE_CONTINUATION_PLAN.a1.json'),reason='Source-review repairs and32MGCguidance after retained402.8MiBmarginfailure, before main physical observations; exact chemistry bytes and original producer retained; no SCF rerun or numerical retuning'))
k.require(g.finalize(new),'Reused input finalization must retain all original cases')
print('BYTE_IDENTICAL_CHEMISTRY_REBOUND',k.committed_sources(),k.sha(new/'INPUTS_FROZEN.json'))
