"""Authenticate the complete exact-circuit delivery and deterministic reports."""
import argparse,json,shutil,tempfile
from pathlib import Path
import study
from compiler import require
from reporting import report

MAX_FILES=4096
MAX_BYTES=256*1024**2

def inventory(root):
    paths=list(root.rglob('*'));require(not any(p.is_symlink() for p in paths),'Delivery symlinks refused')
    files=[p for p in paths if p.is_file()]
    require(len(files)<=MAX_FILES and sum(p.stat().st_size for p in files)<=MAX_BYTES,'Delivery storage ceiling')
    return {str(p.relative_to(root)) for p in files}

def pack(results,report_path,destination):
    require(not destination.exists(),'Fresh delivery required; preserve prior attempts')
    inventory(results);inventory(report_path);destination.mkdir(parents=True)
    shutil.copytree(results,destination/'results');shutil.copytree(report_path,destination/'report')
    manifest=dict(files={n:study.sha(destination/n) for n in sorted(inventory(destination))},sources=study.sources(),
        scope='Complete fixed-grid exact compiled circuits including censored prefixes, all physical rows/states, all target decisions and resources/cost reports')
    study.write(destination/'MANIFEST.json',manifest)

def check(delivery,replay_decisions_only=False):
    manifest=study.read(delivery/'MANIFEST.json')
    require(manifest['sources']==study.sources(),'Current delivery source versions differ')
    require(inventory(delivery)==set(manifest['files'])|{'MANIFEST.json'},'Delivery inventory differs')
    mutable={'results/VERIFIED.json','results/INDEPENDENT_ERRORS.json','results/INDEPENDENT_SELECTIONS.json'}
    for name,digest in manifest['files'].items():
        study.relative(name)
        if replay_decisions_only and name in mutable:continue
        require(study.sha(delivery/name)==digest,'Delivery artifact changed: '+name)
    output=delivery/'results';protocol,frozen=study.bind(output,independent=True);verified=study.read(output/'VERIFIED.json')
    require(verified['status']=='passed' and verified['decision_agreement'],'Unverified compiled panel')
    require(verified['verification_sources']==study.verifier_sources(),'Independent source version differs')
    require(verified['producing_source_commit']==protocol['source_commit'] and verified['producing_sources']==protocol['producer_sources'],'Producing/verification identities conflated')
    require(verified['schedule_freeze_sha256']==study.sha(output/'SCHEDULES_FROZEN.json'),'Frozen candidate identity differs')
    require(verified['candidates']==frozen['candidates']==protocol['expected_candidates'] and verified['curves']==protocol['expected_curves'],'Incomplete candidate grid')
    require(verified['feasible']==frozen['feasible'] and verified['censored']==frozen['censored'] and verified['target_decisions']==2*protocol['expected_curves'],'Incomplete feasible/censored target coverage')
    physical_names={folder+'/'+name for folder in ('rows','states') for name in inventory(output/folder)}
    physical_meta={'SELECTIONS.json','RUN_COMPLETE.json','RUN_COSTS.json'}
    require(set(verified['artifacts'])==physical_names|physical_meta,'Verified physical binding inventory differs')
    expected_results=set(frozen['artifacts'])|physical_names|physical_meta|{'SCHEDULES_FROZEN.json','INDEPENDENT_ERRORS.json','INDEPENDENT_SELECTIONS.json','VERIFIED.json'}
    require(inventory(output)==expected_results,'Unexpected or missing compiled result artifact')
    for name,digest in verified['artifacts'].items():require(study.sha(output/study.relative(name))==digest,'Verified physical artifact changed')
    require(study.read(output/'INDEPENDENT_SELECTIONS.json')==study.read(output/'SELECTIONS.json'),'Independent target decisions differ')
    independent=study.read(output/'INDEPENDENT_ERRORS.json')
    require(len(independent)==frozen['feasible'] and len({r['candidate'] for r in independent})==len(independent),'Independent physical coverage differs')
    require(verified['max_state_delta']<=2e-10 and verified['max_metric_delta']<=2e-12,'Independent numerical tolerance exceeded')
    with tempfile.TemporaryDirectory(prefix='compiled-report-check-',dir=delivery.parent) as directory:
        regenerated=Path(directory)/'report';report(output,regenerated);supplied=delivery/'report'
        require(inventory(regenerated)==inventory(supplied),'Report inventory differs')
        for name in inventory(regenerated):
            if replay_decisions_only and name in ('README.md','MANIFEST.json'):continue
            require((regenerated/name).read_bytes()==(supplied/name).read_bytes(),'Deterministic report regeneration differs: '+name)
    print('COMPILED DELIVERY PASSED',len(manifest['files']),'artifacts;',verified['candidates'],'candidates;',verified['target_decisions'],'target decisions',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('phase',choices=('pack','check'))
    for name in ('delivery','results','report'):parser.add_argument('--'+name,type=Path,required=name=='delivery')
    parser.add_argument('--replay-decisions-only',action='store_true');a=parser.parse_args()
    if a.phase=='pack':
        require(a.results is not None and a.report is not None,'Results and report required')
        pack(a.results.resolve(),a.report.resolve(),a.delivery.resolve())
    else:check(a.delivery.resolve(),a.replay_decisions_only)
