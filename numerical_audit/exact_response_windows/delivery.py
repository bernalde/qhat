"""Archive/check exact-window evidence with complete inventory and source binding."""
import argparse,hashlib,json,shutil,tempfile
from pathlib import Path
from reporting import generate

HERE=Path(__file__).resolve().parent


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def source_hashes():
    paths=[p for p in HERE.iterdir() if p.suffix in ('.py','.md','.txt')]+[HERE.parent/'prospective_reference/engine.py']
    return {str(p.relative_to(HERE.parent)):sha(p) for p in sorted(paths)}


def pack(results,report,destination,prior):
    require(not destination.exists(),'Fresh delivery path required');destination.mkdir(parents=True)
    shutil.copytree(results,destination/'results');shutil.copytree(report,destination/'report')
    require(prior.is_dir(),'Retained prior pilot required')
    shutil.copytree(prior,destination/'prior')
    manifest=dict(files={str(p.relative_to(destination)):sha(p) for p in sorted(destination.rglob('*')) if p.is_file()},
        sources=source_hashes(),policy='Complete supplied rational nine-unit certificate study; original completed eight-unit pilot retained separately')
    (destination/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')


def check(delivery,replay_decisions_only=False):
    m=json.loads((delivery/'MANIFEST.json').read_text());require(m['sources']==source_hashes(),'Active delivery source differs')
    actual={str(p.relative_to(delivery)) for p in delivery.rglob('*') if p.is_file() and p!=delivery/'MANIFEST.json'}
    require(actual==set(m['files']),'Delivery inventory differs')
    mutable={'results/VERIFIED.json','results/INDEPENDENT_ERRORS.json','results/INDEPENDENT_CERTIFICATES.json'}
    for p,h in m['files'].items():
        if replay_decisions_only and p in mutable:continue
        require(sha(delivery/p)==h,'Delivery artifact differs: '+p)
    results=delivery/'results';protocol=json.loads((results/'PROTOCOL_FROZEN.json').read_text())
    verified=json.loads((results/'VERIFIED.json').read_text())
    require(verified['status']=='passed' and verified['decision_agreement'] and verified['rows']==protocol['expected_circuits'],'Failed/incomplete verification')
    require(verified['schedule_freeze_sha256']==sha(results/'SCHEDULES_FROZEN.json'),'Verification uses different schedules')
    for p,h in verified['artifacts'].items():require(sha(results/p)==h,'Independent artifact binding differs')
    require(set(p.stem for p in (results/'inputs').glob('*.json'))==set(protocol['cases']),'Input coverage differs')
    require(len(list((results/'rows').glob('*.json')))==protocol['expected_circuits'] and len(list((results/'states').glob('*.npz')))==protocol['expected_circuits'],'Physical grid differs')
    independent=json.loads((results/'INDEPENDENT_CERTIFICATES.json').read_text())
    require({r['case'] for r in independent}==set(protocol['cases']),'Independent certificate coverage differs')
    for row in independent:
        c=json.loads((results/'certificates'/(row['case']+'.json')).read_text())
        for key in ('incumbent_integer','selected_window','minimizers','evaluated'):require(row[key]==c[key],'Certificate verifier differs')
        require(c['complete'] and c['evaluated']==362880 and c['lower_bound_integer']==c['incumbent_integer'] and c['absolute_gap_integer']=='0','Incomplete exact certificate')
    with tempfile.TemporaryDirectory() as temporary:
        regenerated=Path(temporary)/'report';generate(results,regenerated)
        supplied=delivery/'report';require({p.name for p in supplied.iterdir()}=={p.name for p in regenerated.iterdir()},'Report inventory differs')
        for p in regenerated.iterdir():
            if replay_decisions_only and p.name in ('README.md','MANIFEST.json'):continue
            require(p.read_bytes()==(supplied/p.name).read_bytes(),'Report regeneration differs: '+p.name)
    print('DELIVERY PASSED',len(m['files']),'artifacts;',verified['rows'],'physical rows;',len(independent),'complete9!certificates',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('pack','check'))
    p.add_argument('--delivery',type=Path,required=True);p.add_argument('--results',type=Path);p.add_argument('--report',type=Path);p.add_argument('--prior',type=Path)
    p.add_argument('--replay-decisions-only',action='store_true');a=p.parse_args()
    if a.phase=='pack':pack(a.results.resolve(),a.report.resolve(),a.delivery.resolve(),a.prior.resolve())
    else:check(a.delivery.resolve(),a.replay_decisions_only)
