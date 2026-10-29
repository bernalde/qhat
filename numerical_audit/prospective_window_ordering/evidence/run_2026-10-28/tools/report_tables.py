"""Extract report tables from the already verified main result; no simulation."""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('root', type=Path)
p.add_argument('--report-sha256', required=True)
a = p.parse_args()
root = a.root.resolve()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
report_path = root/'INDEPENDENT_DERIVED_CHECK.main.H.a1.json'
assert sha(report_path) == a.report_sha256
d = json.loads(report_path.read_text())
assert d['status'] == 'verified_complete_derived_panel' and not d['failures']
assert d['complete_352_derived_claim_allowed'] and d['unqualified_under_ceiling_claim_allowed']
spec = json.loads((root/'INDEPENDENT_DERIVED_INDEX.main.H.a1.json').read_text())
desc = spec['analyzer_summary']
summary_path = root/desc['path']
assert sha(summary_path) == desc['sha256']
s = json.loads(summary_path.read_text())
totals = defaultdict(float)
counts = Counter()
for row in d['constructor_reconciliation']:
    totals[row['arm']] += row['final_conservative_CPU_seconds']
    counts[row['arm']] += 1
assert len(d['constructor_reconciliation']) == 64 and set(counts.values()) == {8}
barrier = json.loads((root/'MAIN_PHYSICS_BARRIER.a1.json').read_text())
assert abs(sum(totals.values())-barrier['constructor64_conservative_CPU_seconds']) <= 1e-8
arms = ['canonical','magnitude','current','first','pauli','cubic','window_sparse','window_magnitude']
costs = [dict(arm=arm, calls=counts[arm], raw_CPU_seconds=s['construction_by_arm'][arm]['CPU_seconds'],
              conservative_FINAL_CPU_seconds=totals[arm],
              absolute_wall_seconds=s['construction_by_arm'][arm]['wall_seconds']) for arm in arms]
assert len(s['comparisons']) == 12 and all(r['expected_pairs'] == 32 for r in s['comparisons'])
steps = {str(formula): sorted({r['physical']['steps'] for r in s['raw_rows'] if r['formula']==formula})
         for formula in (1,2)}
out = dict(schema='verified-main-report-tables-v1', source_report_sha256=a.report_sha256,
           source_summary=desc, extractor_sha256=sha(Path(__file__)), costs=costs,
           comparisons=s['comparisons'], physical_rows=352,
           physical_floor_rows=sum(r['physical_floor'] for r in s['raw_rows']),
           actual_validation_steps_by_formula=steps,
           sparse_over_cubic_raw_CPU_ratio=s['construction_by_arm']['window_sparse']['CPU_seconds']/s['construction_by_arm']['cubic']['CPU_seconds'],
           sparse_over_cubic_wall_ratio=s['construction_by_arm']['window_sparse']['wall_seconds']/s['construction_by_arm']['cubic']['wall_seconds'],
           scope='Tables only; paired observations are not independent chemistry samples. Post-cutoff bookkeeping is separate.')
with (root/'REPORT_NUMERIC_TABLES.a1.json').open('x') as stream:
    json.dump(out,stream,indent=2,allow_nan=False); stream.write('\n')
print('REPORT_TABLES_WRITTEN',len(costs),'arms',len(out['comparisons']),'comparisons')
