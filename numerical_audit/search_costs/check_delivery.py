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
    cumulative=root/'execution/cumulative.json'
    if cumulative.exists():
        ledger=read(cumulative);attempts=ledger['attempts']
        assert ledger['phase_ceilings']=={
            'tests':{'CPU_seconds':120,'wall_seconds':900},
            'preflight':{'CPU_seconds':120,'wall_seconds':900},
            'main':{'CPU_seconds':1800,'wall_seconds':10800},
            'independent':{'CPU_seconds':1800,'wall_seconds':10800},
            'reporting':{'CPU_seconds':300,'wall_seconds':1800}}
        assert {r['phase'] for r in attempts}<=set(ledger['phase_ceilings'])
        assert len({r['label'] for r in attempts})==len(attempts)
        resources=read(root/'execution/resources.json')
        assert {r['label'] for r in attempts}=={r['label'] for r in resources['phases']}
        for phase,budget in ledger['phase_ceilings'].items():
            group=[r for r in attempts if r['phase']==phase]
            assert all(r['CPU_charge_seconds']>=0 and r['controller_wall_seconds']>=0 for r in group)
            assert sum(r['CPU_charge_seconds'] for r in group)<=budget['CPU_seconds']
            assert sum(r['controller_wall_seconds'] for r in group)<=budget['wall_seconds']
        interrupted=[r for r in attempts if r['label']=='i05-delivery-check-01']
        assert len(interrupted)==1 and interrupted[0]['outcome']=='interrupted; incomplete checker'
        assert interrupted[0]['CPU_charge_seconds']>=61
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
    designs={case:read(results/'designs'/(case+'.json')) for case in cases}
    expected = set(itertools.product(cases,range(54),(1,2)))
    assert len(cases) == 8 and len(rows) == 864 and set(physical) == expected
    # Warmups are executed work, separate from the54 measured order records.
    expected_warmups=set(itertools.product(('conserving','pauli'),('JD','J0'),('cold','placement')))
    for case in cases:
        design=designs[case]
        model=read(results/'inputs'/(case+'.json'))
        warmups=design['warmups']
        assert len(warmups)==8 and {(x['partition'],x['objective'],x['initialization']) for x in warmups}==expected_warmups
        common={x['record']['common_setup_CPU_charge'] for x in warmups}
        assert len(common)==1
        for item in warmups:
            r=item['record'];m=len(model['blocks']) if item['partition']=='conserving' else len(model['coefficients'])
            assert r['CPU_budget']==.5 and (r['objective'],r['initialization'])==(item['objective'],item['initialization'])
            assert sorted(r['order'])==list(range(m))
            total=r['common_setup_CPU_charge']+r['inherited_placement_CPU_charge']+r['actual_local_CPU_seconds']
            assert abs(total-r['total_charged_CPU_seconds'])<1e-12 and r['within_CPU_budget']==(total<=.5)
            local=sum(r[k] for k in ('initializer_CPU_seconds','setup_CPU_seconds','search_CPU_seconds','accounting_overhead_CPU_seconds'))
            assert abs(local-r['actual_local_CPU_seconds'])<1e-12
            assert r['evaluation_cap']==20*m*(m-1) and 0<=r['evaluations']<=r['evaluation_cap']
            assert (r['inherited_placement_CPU_charge']==0)==(item['initialization']=='cold')
            if 'starts' in r:
                assert len(r['starts'])==2 and sum(x['evaluations'] for x in r['starts'])==r['evaluations']
                for start in r['starts']:
                    assert sum(x['evaluations'] for x in start['sweeps'])==start['evaluations'] and len(start['accepted'])<=10
                    for sweep in start['sweeps']:
                        assert 0<=sweep['evaluations']<=m*(m-1)
                        assert not sweep['complete'] or sweep['evaluations']==m*(m-1)
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
                design = designs[case]
                control = [i for i,x in enumerate(design['records']) if x['arm'] == r['control'] and x['repetition'] == int(r['repetition'])]
                assert len(control) == 1
                left = case,int(r['record_index']),int(r['formula'])
                right = case,control[0],int(r['control_formula'])
                stored_difference = float(r['error_difference'])
            a,b = physical[left]['error'],physical[right]['error']
            x,y = alternative[left],alternative[right]
            assert stored_difference == a-b and r['classification'] == classify(a-b) == classify(x-y)
            assert r['denominator_above_floor'] == str(b>1e-12) == str(y>1e-12)
    del designs
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
