"""Package only authenticated scientific results, reports and execution measurements."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def package(results, report, logs, destination, labels):
    if destination.exists():
        raise ValueError('Fresh delivery directory required; preserve earlier attempts')
    verified = read(results/'VERIFIED.json')
    assert verified['status'] == 'passed' and verified['rows'] == 864 and verified['decision_agreement']
    destination.mkdir(parents=True)
    shutil.copytree(results,destination/'evidence/results')
    shutil.copytree(report,destination/'report')
    execution = destination/'evidence/execution'
    execution.mkdir()
    measurements = []
    for label in labels:
        record = read(logs/(label+'.resources.json'))
        assert label.startswith('i05-') and 'argon' not in label and 'astra' not in label
        measurements.append(dict(label=label,exit_code=record['exit_code'],settings=record['settings'],
                                 allowance=record['budget'],service_wall_seconds=record['elapsed'],
                                 process_measurements=record['process_measurements'],
                                 minimum_available_bytes=min((r['available_bytes'] for r in record['samples']),default=None),
                                 samples_while_frozen=sum(r['frozen'] for r in record['samples']),
                                 termination=record['termination']))
    (execution/'resources.json').write_text(json.dumps(dict(phases=measurements,
        scope='Explicit completed phase labels only; no reviewer-client metadata or task-private paths. Missing process fields are unmeasured, not zero.'),indent=2)+'\n')
    source = Path(__file__).resolve().parent
    active = ['reporting.py','check_delivery.py','deliver.py','guarded.py','test_bounded.py','preflight.py','requirements.txt']
    manifest = dict(files={str(p.relative_to(destination/'evidence')):sha(p) for p in sorted((destination/'evidence').rglob('*')) if p.is_file()},
                    delivery_sources={n:sha(source/n) for n in active},frozen_scientific_sources=verified['sources'],
                    report_manifest_sha256=sha(report/'MANIFEST.json'),producing_source_commit=read(results/'PROTOCOL_FROZEN.json')['source_commit'])
    (destination/'evidence/MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status='packaged',files=len(manifest['files']),bytes=sum(p.stat().st_size for p in destination.rglob('*') if p.is_file()))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--logs',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--phase-label',action='append',required=True)
    args = parser.parse_args()
    package(args.results,args.report,args.logs,args.destination,args.phase_label)
