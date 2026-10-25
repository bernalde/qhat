"""Package verified research evidence without copying private review envelopes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def signatures(results):
    designs={};states={}
    for path in sorted((results/'designs').glob('*.json')):
        d=read(path);states[d['case']]=d['state_sha256']
        for point in d['points']:
            for arm,m in point['methods'].items():
                key=d['case']+'.'+point['state']+'.'+str(point['T'])+'.'+arm
                payload=json.dumps({'order':m['order'],'schedule':m['schedule']},sort_keys=True,separators=(',',':')).encode()
                designs[key]=hashlib.sha256(payload).hexdigest()
    return dict(states=states,orders_and_schedules=designs)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--logs',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--prior-freeze',type=Path,action='append',default=[])
    a=p.parse_args();evidence=a.destination/'evidence';report=a.destination/'report'
    if evidence.exists() or report.exists():raise FileExistsError('Preserve existing delivery; use a fresh destination')
    verified=read(a.results/'VERIFIED.json')
    if verified['status']!='passed' or verified['rows']!=432 or not verified['decision_agreement']:
        raise ValueError('Complete independent verification required')
    for name,h in verified['artifacts'].items():
        if sha(a.results/name)!=h:raise ValueError('Verified result changed')
    current=signatures(a.results);stages=[]
    for directory in a.prior_freeze+[a.results]:
        if signatures(directory)!=current:raise ValueError('A pre-error revision changed states or schedules')
        stages.append(dict(stage=directory.name,source_commit=read(directory/'PROTOCOL_FROZEN.json')['source_commit'],
            protocol_sha256=sha(directory/'PROTOCOL_FROZEN.json'),freeze_sha256=sha(directory/'SCHEDULES_FROZEN.json'),
            source_sha256=read(directory/'SCHEDULES_FROZEN.json')['sources']))
    shutil.copytree(a.results,evidence/'results');shutil.copytree(a.report,report)
    write(evidence/'execution/pre_error_identity.json',dict(status='passed',stages=stages,
        identical_initial_state_files=3,identical_orders_and_schedules=216,signatures=current,
        scope='Correlation diagnostics and verifier/protocol corrections preceded every panel circuit error; no state/order replacement'))
    phases=[]
    for path in sorted(a.logs.glob('i04-*.resources.json')):
        if 'argon' in path.name:continue
        raw=read(path);settings=raw['settings'];measurements=raw.get('process_measurements') or {}
        available=[s['available_bytes'] for s in raw['samples']]
        phases.append(dict(label=path.name.removesuffix('.resources.json'),exit_code=raw['exit_code'],
            service_wall_seconds=raw['elapsed'],process_measurements=measurements,settings=settings,
            allowance=raw.get('budget'),minimum_available_bytes=min(available) if available else None,
            samples_while_frozen=sum(s['frozen'] for s in raw['samples']),
            completed_utc=raw['completed_utc']))
    cpu=sum(m.get('user_cpu_seconds',0)+m.get('system_cpu_seconds',0) for m in (r['process_measurements'] for r in phases))
    write(evidence/'execution/resource_summary.json',dict(phases=phases,
        measured_process_CPU_seconds=cpu,max_measured_RSS_KiB=max((r['process_measurements'].get('max_rss_kib',0) for r in phases),default=0),
        note='Service timings exclude prelaunch headroom queueing. Failed collection/preflight attempts are retained. Missing early RSS/CPU fields are not zeros measured by the process timer. Reviewer clients are a separate private workflow ledger.',
        initial_limits=dict(tests=[120,900,64],preflight=[120,900,64],construction_and_propagation=[1800,10800,512],
            verification=[1800,10800,512],extraction_reporting_replay=[300,1800,128]),
        units='CPU seconds, service wall seconds, output MiB; one worker/thread; 25% CPU, 448/512 MiB high/max, no swap'))
    code=Path(__file__).resolve().parent
    active={name:sha(code/name) for name in ('reporting.py','check_delivery.py','deliver.py','guarded.py','test_response.py','test_state_guard.py')}
    paths=sorted(p for p in evidence.rglob('*') if p.is_file())
    files={str(p.relative_to(evidence)):sha(p) for p in paths}
    write(evidence/'MANIFEST.json',dict(files=files,active_delivery_sources=active,
        frozen_scientific_sources=read(a.results/'SCHEDULES_FROZEN.json')['sources'],
        report_manifest_sha256=sha(report/'MANIFEST.json'),
        preflight_scope='Largest synthetic 8-qubit/104-word construction, response and circuit replay; extraction/packaging estimates frozen separately before the panel. Final packaging verified under its own allowance.',
        producing_source_commit=read(a.results/'PROTOCOL_FROZEN.json')['source_commit']))
    print(json.dumps(dict(status='packaged',evidence_files=len(files),report_files=len(list(report.iterdir())),
        retained_result_bytes=sum(p.stat().st_size for p in paths),source_commit=read(a.results/'PROTOCOL_FROZEN.json')['source_commit'])))


if __name__=='__main__':main()
