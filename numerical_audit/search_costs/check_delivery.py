"""Authenticate the full controlled-cost bundle and independently check decisions."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import tempfile
from reporting import build


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def check(results, report, manifest, replay_decisions_only=False):
    publication = read(manifest)
    root = manifest.parent
    assert {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p != manifest} == set(publication['files'])
    for name, digest in publication['files'].items():
        assert sha(root/name) == digest, name
    assert sha(report/'MANIFEST.json') == publication['report_manifest_sha256']
    for name, digest in publication['delivery_sources'].items():
        assert sha(Path(__file__).resolve().parent/name) == digest, name
    freeze = read(results/'SCHEDULES_FROZEN.json')
    audit = Path(__file__).resolve().parent.parent
    for name, digest in freeze['sources'].items():
        assert sha(audit/name) == digest, name
    previous = audit/'prospective_reference/evidence/results/inputs'
    for case, digest in freeze['prior_models'].items():
        assert sha(previous/case/'model.json') == digest, case
    for name, digest in freeze['artifacts'].items():
        assert sha(results/name) == digest, name
    verified = read(results/'VERIFIED.json')
    assert verified['status'] == 'passed' and verified['rows'] == 864 and verified['decision_agreement']
    assert verified['sources'] == freeze['sources'] == publication['frozen_scientific_sources']
    assert verified['schedule_freeze_sha256'] == sha(results/'SCHEDULES_FROZEN.json')
    for name, digest in verified['artifacts'].items():
        assert sha(results/name) == digest, name
    rows = [read(p) for p in sorted((results/'rows').glob('*.json'))]
    key = lambda r: (r['case'],r['record_index'],r['formula'])
    physical = {key(r):r for r in rows}
    cases = {p.stem for p in (results/'designs').glob('*.json')}
    expected = set(itertools.product(cases,range(54),(1,2)))
    assert len(cases) == 8 and len(rows) == 864 and set(physical) == expected
    independent = read(results/'INDEPENDENT_ERRORS.json')
    assert independent['decision_agreement'] and not independent['ambiguities']
    alternative = {key(r):r['independent_error'] for r in independent['rows']}
    assert set(alternative) == expected
    classify = lambda d: 'win' if d < -1e-12 else 'loss' if d > 1e-12 else 'tie'
    for name, expected_count in [('paired_factorial.csv',1152),('strong_controls.csv',3456)]:
        with (report/name).open(newline='') as f:
            comparisons = list(csv.DictReader(f))
        assert len(comparisons) == expected_count
        for r in comparisons:
            case = r['case']
            if name == 'paired_factorial.csv':
                left = case,int(r['left_record']),int(r['formula'])
                right = case,int(r['right_record']),int(r['formula'])
                stored_difference = float(r['error_difference'])
            else:
                design = read(results/'designs'/(case+'.json'))
                control = [i for i,x in enumerate(design['records']) if x['arm'] == r['control'] and x['repetition'] == int(r['repetition'])]
                assert len(control) == 1
                left = case,int(r['record_index']),int(r['formula'])
                right = case,control[0],int(r['control_formula'])
                stored_difference = float(r['error_difference'])
            a,b = physical[left]['error'],physical[right]['error']
            x,y = alternative[left],alternative[right]
            assert stored_difference == a-b and r['classification'] == classify(a-b) == classify(x-y)
            assert r['denominator_above_floor'] == str(b>1e-12) == str(y>1e-12)
    if not replay_decisions_only:
        with tempfile.TemporaryDirectory(prefix='cost-report-replay-') as scratch:
            target = Path(scratch)/'report'
            build(results,target)
            assert {p.name for p in target.iterdir()} == {p.name for p in report.iterdir()}
            for path in target.iterdir():
                assert path.read_bytes() == (report/path.name).read_bytes(), path.name
    print(json.dumps(dict(status='passed',rows=864,orders=432,factorial_pairs=1152,strong_control_pairs=3456,
                          scientific_decisions_checked=True,published_manifest_checked=True,
                          report_byte_regeneration_checked=not replay_decisions_only)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--replay-decisions-only',action='store_true')
    args = parser.parse_args()
    check(args.results,args.report,args.manifest,args.replay_decisions_only)
