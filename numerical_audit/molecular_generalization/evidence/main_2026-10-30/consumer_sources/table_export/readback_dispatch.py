"""Bounded metadata export plus independent CSV-to-original-row readback.

Run inside the declared analysis guard. No scores, states, medians or primary
decisions are recomputed. The original incomplete verdict stays unchanged.
"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('CSV readback requires normal Python assertion checking')

plan_path = Path(sys.argv[1]).resolve()
plan = json.loads(plan_path.read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(__file__) == plan['wrapper_sha256']
data, output = Path(plan['data_root']), Path(plan['tables_output'])
consumer = Path(plan['export_consumer'])
assert sha(consumer) == plan['export_consumer_sha256']
assert sha(data/'ANALYSIS.json') == plan['analysis_sha256']
assert sha(data/'EVIDENCE_INVENTORY.json') == plan['inventory_sha256']
assert not output.exists()
subprocess.run([sys.executable, str(consumer), str(data), str(data/'ANALYSIS.json'), str(output)], check=True)

original = json.loads((data/'ANALYSIS.json').read_text())
manifest = json.loads((output/'MANIFEST.json').read_text())
assert manifest['export_delivery_complete'] is True and manifest['scientific_gate_complete'] is False
assert manifest['source_analysis_complete'] is False and manifest['source_analysis_status'] == 'incomplete'
assert json.loads((output/'PRIMARY_QUESTIONS.json').read_text()) == original['primary']
assert all(q['decision'] == 'incomplete' for q in original['primary'].values())
for name, identity in manifest['outputs'].items():
    assert (output/name).stat().st_size == identity['bytes'] and sha(output/name) == identity['sha256']


def rows(name, expected):
    with (output/name).open(newline='') as stream:
        value = list(csv.DictReader(stream))
    assert len(value) == expected == manifest['csv'][name]['rows']
    return value


def order_key(row, *, csv_row=False):
    cast = int if csv_row else lambda x: x
    return (row['case'], cast(row['repetition']), row['arm'], cast(row['cpu_budget']))


def physical_key(row, *, csv_row=False):
    return order_key(row, csv_row=csv_row)+(float(row['time']), int(row['budget_multiplier']))


costs = {order_key(row): row for row in original['constructor_costs']['constructors']}
constructor_rows = rows('constructors.csv', 36)
assert {order_key(row, csv_row=True) for row in constructor_rows} == set(costs)
for row in constructor_rows:
    key = order_key(row, csv_row=True); cost = costs[key]
    source = json.loads((data/row['constructor_path']).read_text())
    assert sha(data/row['constructor_path']) == row['constructor_sha256']
    assert float(row['selected_J']) == source['score'] and float(row['selected_cubic_response']) == source['cubic_response']
    assert int(row['selected_steps']) == source['steps'] and json.loads(row['eligible']) is source['eligible']
    for prefix, original_boundary in (('actual_own', cost['actual_own']), ('comparison_attributed', cost['comparison_attributed'])):
        for metric in ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds'):
            assert float(row[prefix+'_'+metric]) == original_boundary[metric]

physical = {physical_key(row): row for row in original['retained_physical_rows']}
physical_rows = rows('physical.csv', 144)
assert {physical_key(row, csv_row=True) for row in physical_rows} == set(physical)
for row in physical_rows:
    actual = physical[physical_key(row, csv_row=True)]
    for metric in ('projective_error', 'infidelity', 'norm_drift', 'number_leakage', 'spin_leakage'):
        assert float(row[metric]) == actual[metric]
    for metric in ('steps', 'ideal_rotations', 'rotation_budget'):
        assert int(row[metric]) == actual[metric]
    for field in ('actual_sha256', 'reference_sha256', 'circuit_sha256', 'status'):
        assert row[field] == actual[field]
    assert json.loads(row['numerical_floor']) is actual['numerical_floor']
    assert json.loads(row['scientific_gate_complete']) is False


def pair_key(row):
    return (row['case'], int(row['repetition']), row['question'], int(row['cpu_budget']),
            float(row['time']), int(row['budget_multiplier']))


pairs = {pair_key(row): row for row in original['named_pairs']}
pair_rows = rows('named_pairs.csv', 64)
assert {pair_key(row) for row in pair_rows} == set(pairs)
for row in pair_rows:
    actual = pairs[pair_key(row)]
    assert float(row['difference']) == actual['difference']
    assert json.loads(row['eligible']) is actual['eligible'] and row['outcome'] == actual['outcome']
    assert json.loads(row['scientific_gate_complete']) is False

claims = {(row['case'], row['name']): row for row in original['retained_raw_witnesses']}
raw_rows = rows('persisted_raw_claims.csv', 14)
assert {(row['case'], row['name']) for row in raw_rows} == set(claims)
for row in raw_rows:
    actual = claims[(row['case'], row['name'])]
    assert row['status'] == actual['status'] and json.loads(row['proposed_comparison']) == actual['proposed_comparison']
    if actual['status'] == 'native_refusal_censored':
        assert row['name'] == 'rep2.window.B64.raw1' and json.loads(row['native_score']) is None
        assert json.loads(row['native_comparison']) is None and row['native_refusal'] == actual['native_refusal']

cheap_rows = rows('cheap_context.csv', 192)
cheap_keys = set()
for row in cheap_rows:
    case, rep, arm, budget, duration, cap = (row['case'], int(row['repetition']), row['candidate_arm'],
        int(row['candidate_cpu_budget']), float(row['time']), int(row['budget_multiplier']))
    key = (case, rep, arm, budget, duration, cap, row['cheap_arm'])
    assert key not in cheap_keys; cheap_keys.add(key)
    candidate = physical[(case, rep, arm, budget, duration, cap)]
    baseline = physical[(case, rep, row['cheap_arm'], 0, duration, cap)]
    assert float(row['candidate_error']) == candidate['projective_error']
    assert float(row['cheap_error']) == baseline['projective_error']
    assert float(row['raw_candidate_minus_cheap']) == candidate['projective_error']-baseline['projective_error']
    assert row['primary_decision_role'] == 'none' and json.loads(row['descriptive_only']) is True
    assert json.loads(row['scientific_gate_complete']) is False

assert sha(data/'ANALYSIS.json') == plan['analysis_sha256']
assert sha(data/'EVIDENCE_INVENTORY.json') == plan['inventory_sha256']
result = dict(schema='molecular-censored-export-readback-v1', metadata_delivery_verified=True,
    scientific_gate_complete=False, source_analysis_complete=False,
    source_primary_decisions={key: value['decision'] for key, value in original['primary'].items()},
    selected_physical_rows=144, constructor_rows=36, named_pair_rows=64, persisted_raw_claim_rows=14,
    descriptive_cheap_rows=192, original_native_censored_claims=1, native_score_backfilled=False,
    analysis_sha256=plan['analysis_sha256'], inventory_sha256=plan['inventory_sha256'],
    exported_manifest_sha256=sha(output/'MANIFEST.json'), wrapper_sha256=sha(__file__),
    scope='CSV readback against original sealed JSON rows and exact decision objects; no new scientific evaluation or decision')
with Path(plan['readback_output']).open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print('CENSORED_CSV_READBACK_PASS_SCIENTIFIC_GATE_STILL_INCOMPLETE', flush=True)
