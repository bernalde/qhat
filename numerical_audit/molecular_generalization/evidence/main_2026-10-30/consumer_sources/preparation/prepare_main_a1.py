"""Freeze the accepted, committed prospective main source before new-case work."""
import hashlib,json,shutil,sys
from pathlib import Path
planpath=Path(sys.argv[1]).resolve();plan=json.loads(planpath.read_text());source=Path(plan['producer_root'])
sys.path.insert(0,str(source/'numerical_audit/molecular_generalization'))
import common as k
import runner
k.require(k.committed_sources()==plan['source_commit'],'Exact committed main source')
k.require(k.sha(Path(__file__))==plan['operator_sha256'],'Declared preparation operator')
config=k.configuration();k.require(config['protocol_state']=='frozen_main' and config['main_execution_allowed']is True,'Prospective main protocol frozen')
accepted=k.read(k.HERE/'PREFLIGHT_ACCEPTED.json')
k.require(accepted['complete']is True and accepted['numerical_relocated_replay_complete']is True and accepted['new_family_SCF_training_or_physical_errors_observed']is False,'Completed old gate precedes new case work')
manifestpath=k.HERE/accepted['evidence_manifest']['path'];k.require(k.sha(manifestpath)==accepted['evidence_manifest']['sha256'],'Accepted old manifest exact')
manifest=k.read(manifestpath)
for row in manifest['files']:
    path=k.contained(manifestpath.parent,row['path']);k.require(path.stat().st_size==row['bytes'] and k.sha(path)==row['sha256'],'Accepted evidence payload exact')
limits=Path(plan['phase_limits']);k.require(k.sha(limits)==plan['phase_limits_sha256'],'Frozen phase resource limits')
k.require(k.read(limits)==config['main_phase_limits'],'Main protocol/phase resource join')
data=Path(plan['data_root']);k.freeze(data,mode='main');runner.initialize(data,k.read(limits))
folder=data/'preparation';folder.mkdir()
for path in (planpath,Path(__file__).resolve(),limits,Path(plan['main_resource_plan'])):shutil.copyfile(path,folder/path.name)
k.write_new(folder/'PREPARED.json',dict(complete=True,mode='main',source_commit=plan['source_commit'],source_freeze_sha256=k.sha(data/'INPUT_PROTOCOL_FROZEN.json'),registry_sha256=k.sha(k.HERE/'REGISTRY.json'),phase_limits_sha256=k.sha(data/'PHASE_LIMITS.json'),old_preflight_accepted_sha256=k.sha(k.HERE/'PREFLIGHT_ACCEPTED.json'),expected_constructor_records=36,expected_logical_rows=144,SCF_training_or_physical_error_executed_in_preparation=False))
print('MAIN_SOURCE_PROTOCOL_ENVIRONMENT_CASES_AND_RESOURCES_FROZEN',flush=True)
