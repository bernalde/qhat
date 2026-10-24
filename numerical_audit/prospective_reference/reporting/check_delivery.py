"""Check exported report tables and hashes against authenticated experiment rows."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)


def check(results,report,renderer):
    summary=read(results/'summary.json');provenance=read(report/'provenance.json')
    for key,name in [('summary_sha256','summary.json'),('frozen_protocol_sha256','PROTOCOL_FROZEN.json'),
                     ('schedule_freeze_sha256','ALL_SCHEDULES_FROZEN.json')]:
        require(provenance[key]==sha(results/name),'Report evidence hash: '+name)
    require(provenance['renderer_sha256']==sha(renderer),'Report renderer hash')
    for name,digest in provenance['artifact_sha256'].items():
        require(sha(report/name)==digest,'Changed report artifact: '+name)
    raw=[read(p) for p in sorted((results/'rows').glob('*.json'))]
    with (report/'all_rows.csv').open(newline='') as f:export=list(csv.DictReader(f))
    require(len(raw)==len(export)==216,'Incomplete exported row inventory')
    for original,published in zip(raw,export):
        require(set(original)==set(published),'Changed exported row fields')
        require(all(str(v)==published[k] for k,v in original.items()),'Exported row differs from raw result')
    with (report/'construction_costs.csv').open(newline='') as f:costs=list(csv.DictReader(f))
    methods={}
    for p in (results/'designs').glob('*.json'):
        design=read(p)
        for point in design['points']:
            for arm,method in point['methods'].items():methods[design['case'],point['state'],arm]=method
    require(len(costs)==len(methods)==54,'Incomplete construction inventory')
    seen=set()
    for row in costs:
        key=tuple(row[k] for k in ('case','state','arm'));require(key not in seen,'Duplicate construction row');seen.add(key)
        method=methods[key]
        for field in ('setup_cpu_seconds','search_cpu_seconds'):
            require(float(row[field])==method[field],'Construction timing mismatch')
        require(float(row['total_cpu_seconds'])==method['setup_cpu_seconds']+method['search_cpu_seconds'],
                'Total construction cost mismatch')
        require(int(row['examined_positions'])==method['evaluations'],'Examined position mismatch')
        require(row['reuses_order']==method.get('reuses_order',''),'Order reuse mismatch')
    with (report/'contrasts.csv').open(newline='') as f:export=list(csv.DictReader(f))
    expected=[dict(budget_units=int(b),comparison=k,**v) for b,block in summary['decision']['contrasts'].items() for k,v in block.items()]
    require(len(export)==len(expected)==64,'Incomplete exported contrast inventory')
    for original,published in zip(expected,export):
        require(set(original)==set(published) and all(('' if v is None else str(v))==published[k] for k,v in original.items()),
                'Exported contrast mismatch')
    manifest=results.parent/'MANIFEST.json'
    if manifest.exists():
        files=read(manifest)['file_sha256']
        actual={str(p.relative_to(results)):sha(p) for p in results.rglob('*') if p.is_file()}
        require(files==actual,'Retained result manifest mismatch')
    print(json.dumps(dict(status='passed',rows=len(raw),construction_rows=len(costs),contrasts=len(export))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True);parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--renderer',type=Path,default=Path(__file__).with_name('build_report.py'))
    args=parser.parse_args();check(args.results,args.report,args.renderer)
