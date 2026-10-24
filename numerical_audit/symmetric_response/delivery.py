"""Bind and authenticate the complete symmetric-response delivery."""
import argparse,hashlib,json,shutil,tempfile
from pathlib import Path
from reporting import generate
from study import sources

HERE=Path(__file__).resolve().parent
def require(ok,message):
    if not ok:raise ValueError(message)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source_hashes():
    paths=[p for p in HERE.iterdir() if p.suffix in ('.py','.md','.txt')]
    paths += [HERE.parent/'prospective_reference'/n for n in ('engine.py','one_pass.py','campaign.py')]
    return {str(p.relative_to(HERE.parent)):sha(p) for p in sorted(paths)}

def pack(results,report,destination):
    require(not destination.exists(),'Fresh delivery path required');destination.mkdir(parents=True)
    shutil.copytree(results,destination/'results');shutil.copytree(report,destination/'report')
    manifest=dict(files={str(p.relative_to(destination)):sha(p) for p in sorted(destination.rglob('*')) if p.is_file()},
        sources=source_hashes(),scope='Complete prospective symmetric-specific panel, all controls and losses')
    (destination/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')

def check(delivery,replay_decisions_only=False):
    m=json.loads((delivery/'MANIFEST.json').read_text());require(m['sources']==source_hashes(),'Delivery source changed')
    actual={str(p.relative_to(delivery)) for p in delivery.rglob('*') if p.is_file() and p!=delivery/'MANIFEST.json'}
    require(actual==set(m['files']),'Delivery inventory differs')
    mutable={'results/VERIFIED.json','results/INDEPENDENT_ERRORS.json','results/INDEPENDENT_OBJECTIVES.json'}
    for path,digest in m['files'].items():
        if replay_decisions_only and path in mutable:continue
        require(sha(delivery/path)==digest,'Artifact changed: '+path)
    results=delivery/'results';protocol=json.loads((results/'PROTOCOL_FROZEN.json').read_text())
    verified=json.loads((results/'VERIFIED.json').read_text());freeze=json.loads((results/'SCHEDULES_FROZEN.json').read_text())
    require(verified['status']=='passed' and verified['decision_agreement'],'Unverified study')
    require(verified['rows']==protocol['expected_circuits'],'Incomplete verified grid')
    require(verified['sources']==freeze['sources']==protocol['sources'],'Different scientific source')
    require(verified['verification_sources']==sources(),'Independent verifier source changed')
    require(verified['schedule_freeze_sha256']==sha(results/'SCHEDULES_FROZEN.json'),'Different schedules')
    for p,h in freeze['artifacts'].items():require(sha(results/p)==h,'Frozen artifact changed')
    for p,h in verified['artifacts'].items():require(sha(results/p)==h,'Verified artifact changed')
    expected={('preflight_altered_' if protocol['preflight'] else '')+c for c in protocol['cases']}
    for folder in ('inputs','designs'):require({p.stem for p in (results/folder).glob('*.json')}==expected,'Case coverage differs')
    require(len(list((results/'rows').glob('*.json')))==len(list((results/'states').glob('*.npz')))==protocol['expected_circuits'],'Physical artifacts differ')
    require(verified['comparisons']==28*len(expected)*len(protocol['times'])*len(protocol['budgets']),'Comparison coverage differs')
    with tempfile.TemporaryDirectory() as temporary:
        regenerated=Path(temporary)/'report';generate(results,regenerated);supplied=delivery/'report'
        require({p.name for p in regenerated.iterdir()}=={p.name for p in supplied.iterdir()},'Report inventory differs')
        for p in regenerated.iterdir():
            if replay_decisions_only and p.name in ('README.md','MANIFEST.json'):continue
            require(p.read_bytes()==(supplied/p.name).read_bytes(),'Report regeneration differs: '+p.name)
    print('DELIVERY PASSED',len(m['files']),'artifacts;',verified['rows'],'circuits;',verified['comparisons'],'comparisons',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('pack','check'))
    for name in ('delivery','results','report'):p.add_argument('--'+name,type=Path,required=name=='delivery')
    p.add_argument('--replay-decisions-only',action='store_true');a=p.parse_args()
    if a.phase=='pack':pack(a.results.resolve(),a.report.resolve(),a.delivery.resolve())
    else:check(a.delivery.resolve(),a.replay_decisions_only)
