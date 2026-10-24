"""Authenticate retained evidence and independently check report decisions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import tempfile
from reporting import build


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check(results,report):
    freeze=json.loads((results/'SCHEDULES_FROZEN.json').read_text())
    protocol=json.loads((results/'PROTOCOL_FROZEN.json').read_text())
    for p,h in freeze['artifacts'].items():
        if sha(results/p)!=h:raise ValueError('Frozen state/design/protocol changed')
    # Locate source/input dependencies relative to this exact checkout.
    audit=Path(__file__).resolve().parent.parent
    for p,h in freeze['sources'].items():
        if sha(audit/p)!=h:raise ValueError('Frozen scientific dependency changed')
    previous=audit/'prospective_reference/evidence/results'
    for p,h in freeze['inputs'].items():
        if sha(previous/p)!=h:raise ValueError('Prior target/HF schedule dependency changed')
    verified=json.loads((results/'VERIFIED.json').read_text())
    if verified['status']!='passed' or verified['rows']!=432:raise ValueError('Independent replay incomplete')
    if sha(results/'SCHEDULES_FROZEN.json')!=verified['schedule_freeze_sha256']:raise ValueError('Stale replay freeze')
    for p,h in verified['artifacts'].items():
        if sha(results/p)!=h:raise ValueError('Verified physical artifact changed')
    rows=[json.loads(p.read_text()) for p in sorted((results/'rows').glob('*.json'))]
    lookup={(r['case'],r['state'],r['T'],r['arm'],r['budget_units']):r['error'] for r in rows}
    if len(rows)!=432 or len(lookup)!=432:raise ValueError('Incomplete physical rows')
    summary=json.loads((results/'SUMMARY.json').read_text())
    for item in summary['coverage']:
        retained=[c for c in protocol['cases'] if lookup[c,item['state'],item['T'],'hf_current',10]>1e-12]
        if retained!=item['informative_families'] or item['sufficient']!=(len(retained)>=2):raise ValueError('Coverage classification differs')
    if len(summary['coverage'])!=8:raise ValueError('Missing state/time coverage stratum')
    for item in summary['contrasts']:
        pairs=[(lookup[c,item['state'],item['T'],item['arm'],item['budget_units']],
                lookup[c,item['state'],item['T'],'hf_current',item['budget_units']]) for c in item['families']]
        counts=(sum(a-b < -1e-12 for a,b in pairs),sum(abs(a-b)<=1e-12 for a,b in pairs),sum(a-b>1e-12 for a,b in pairs))
        if counts!=tuple(item[k] for k in ('wins','ties','losses')):raise ValueError('W/T/L classification differs')
        ratios=[a/b for a,b in pairs if b>1e-12]
        if item['floor_exclusions']!=len(pairs)-len(ratios):raise ValueError('Floor inclusion differs')
        median=statistics.median(ratios) if ratios else None
        if median!=item['median_ratio']:raise ValueError('Median ratio differs')
    with tempfile.TemporaryDirectory(prefix='state-report-replay-') as scratch:
        scratch=Path(scratch);build(results,scratch)
        if {p.name for p in scratch.iterdir()}!={p.name for p in report.iterdir()}:raise ValueError('Report file inventory differs')
        for p in scratch.iterdir():
            if p.read_bytes()!=(report/p.name).read_bytes():raise ValueError('Report regeneration differs: '+p.name)
    print('PASS: 432 rows authenticated; coverage, floors, W/T/L and report byte regeneration checked')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();check(a.results,a.report)
