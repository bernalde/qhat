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


def check(results,report,manifest=None,replay_decisions_only=False):
    freeze=json.loads((results/'SCHEDULES_FROZEN.json').read_text())
    if manifest is not None:
        delivery=json.loads(manifest.read_text());root=manifest.parent
        inventory={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=manifest}
        if inventory!=set(delivery['files']):raise ValueError('Published evidence inventory differs')
        for name,h in delivery['files'].items():
            if sha(root/name)!=h:raise ValueError('Published evidence file changed: '+name)
        if delivery['frozen_scientific_sources']!=freeze['sources']:raise ValueError('Delivery scientific-source identity differs')
        if delivery['report_manifest_sha256']!=sha(report/'MANIFEST.json'):raise ValueError('Published report manifest differs')
        for name,h in delivery['active_delivery_sources'].items():
            if sha(Path(__file__).resolve().parent/name)!=h:raise ValueError('Delivery/check source changed: '+name)
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
    independent=json.loads((results/'INDEPENDENT_ERRORS.json').read_text())
    if not independent['decisions']['agreement']:raise ValueError('Independent primary decisions are ambiguous')
    alt={(r['case'],r['state'],r['T'],r['arm'],r['budget_units']):r['error'] for r in independent['rows']}
    with (report/'candidate_contrasts.csv').open(newline='') as f:secondary=list(csv.DictReader(f))
    if len(secondary)!=64:raise ValueError('Incomplete candidate/control contrast inventory')
    classify=lambda d:'loss' if d>1e-12 else 'win' if d < -1e-12 else 'tie'
    for item in secondary:
        families=item['families'].split(';');budget=int(item['budget_units']);control=item['control']
        keys=sorted({(r['case'],r['state'],r['T']) for r in rows if r['case'] in families})
        pairs=[(lookup[c,s,T,'state_current',budget],lookup[c,s,T,control,budget]) for c,s,T in keys]
        other=[(alt[c,s,T,'state_current',budget],alt[c,s,T,control,budget]) for c,s,T in keys]
        for (a,b),(x,y) in zip(pairs,other):
            if classify(a-b)!=classify(x-y) or (b>1e-12)!=(y>1e-12):
                raise ValueError('Numerically ambiguous candidate/control decision; retain evidence without relabelling')
        counts=(sum(a-b < -1e-12 for a,b in pairs),sum(abs(a-b)<=1e-12 for a,b in pairs),sum(a-b>1e-12 for a,b in pairs))
        if counts!=tuple(int(item[k]) for k in ('wins','ties','losses')):raise ValueError('Candidate/control counts differ')
        ratios=[a/b for a,b in pairs if b>1e-12]
        if int(item['floor_exclusions'])!=len(pairs)-len(ratios):raise ValueError('Candidate/control floor inclusion differs')
        median=statistics.median(ratios) if ratios else None
        if ('' if median is None else str(median))!=item['median_ratio']:raise ValueError('Candidate/control median differs')
    with (report/'coherent_phase_pairs.csv').open(newline='') as f:phase_pairs=list(csv.DictReader(f))
    if len(phase_pairs)!=108:raise ValueError('Incomplete fixed coherent phase-pair audit')
    for item in phase_pairs:
        case,T,arm,budget=item['case'],float(item['T']),item['arm'],int(item['budget_units'])
        a=lookup[case,'coherent_real',T,arm,budget];b=lookup[case,'coherent_imag',T,arm,budget]
        x=alt[case,'coherent_real',T,arm,budget];y=alt[case,'coherent_imag',T,arm,budget]
        if classify(a-b)!=classify(x-y):raise ValueError('Numerically ambiguous coherent phase comparison')
        if float(item['error_difference'])!=a-b or item['tie_at_declared_floor']!=str(abs(a-b)<=1e-12):
            raise ValueError('Coherent phase-pair table differs')
    if not replay_decisions_only:
        with tempfile.TemporaryDirectory(prefix='state-report-replay-') as scratch:
            scratch=Path(scratch);build(results,scratch)
            if {p.name for p in scratch.iterdir()}!={p.name for p in report.iterdir()}:raise ValueError('Report file inventory differs')
            for p in scratch.iterdir():
                if p.read_bytes()!=(report/p.name).read_bytes():raise ValueError('Report regeneration differs: '+p.name)
    print(json.dumps(dict(status='passed',rows=432,scientific_decisions_checked=True,
        published_manifest_checked=manifest is not None,report_byte_regeneration_checked=not replay_decisions_only)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--manifest',type=Path)
    p.add_argument('--replay-decisions-only',action='store_true',help='Check fresh numerical replay decisions without regenerating archived diagnostic maxima')
    a=p.parse_args();check(a.results,a.report,a.manifest,a.replay_decisions_only)
