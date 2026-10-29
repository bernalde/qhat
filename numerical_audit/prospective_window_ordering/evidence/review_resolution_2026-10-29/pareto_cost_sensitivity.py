"""Append a postreview CPU-cost precision view; preserve every primary artifact."""
import argparse, hashlib, json, math
from pathlib import Path

SOURCE_SHA = 'db2e4271d53fc140168352cf70a430d4cfecfe8b939db0ed2d7738477cf387f0'
COST_EPS = 1e-12
ERROR_EPS = 1e-12
STAGE = 'postreview_numerical_cost_sensitivity_not_primary'

def need(ok, message):
    if not ok: raise ValueError(message)

def number(value):
    need(type(value) in (int, float) and math.isfinite(value) and value >= 0, 'finite nonnegative metric')
    return value

def dominates(a, b):
    return (a['cost'] <= b['cost'] + COST_EPS and
            a['projective_error'] <= b['projective_error'] + ERROR_EPS and
            (a['cost'] < b['cost'] - COST_EPS or
             a['projective_error'] < b['projective_error'] - ERROR_EPS))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    need(not a.source.is_symlink() and a.source.stat().st_size <= 32*1024**2, 'bounded regular source')
    raw = a.source.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == SOURCE_SHA, 'fixed original summary identity')
    original = json.loads(raw, parse_constant=lambda _: (_ for _ in ()).throw(ValueError('nonfinite JSON')))
    source_groups = original['pareto_memberships']
    need(len(source_groups) == 32, 'complete32group population')
    groups, changed = [], []
    for index, source in enumerate(source_groups):
        old = source['sets']['CPU_seconds']
        points = old['points']
        need(len(points) == 11 and source['expected_paths'] == 11 and old['excluded'] == [], 'all11validphysicalpaths')
        need(len({q['path'] for q in points}) == 11 and all(q['fully_validated'] is True for q in points), 'distinct validated paths')
        for point in points:
            number(point['cost']); number(point['projective_error'])
        dom = {b['path']: [q['path'] for q in points if q['path'] != b['path'] and dominates(q, b)] for b in points}
        frontier = [q['path'] for q in points if not dom[q['path']]]
        differs = frontier != old['pareto_paths']
        if differs: changed.append(index)
        groups.append(dict(index=index, case=source['case'], repetition=source['repetition'],
            time=source['time'], budget_multiplier=source['budget_multiplier'], expected_paths=source['expected_paths'],
            CPU_seconds=dict(points=points, pareto_paths=frontier, dominated_by=dom, excluded=[],
                dominance='costA<=costB+1e-12 and errorA<=errorB+1e-12; costA<costB-1e-12 or errorA<errorB-1e-12. Postreview CPU-cost sensitivity only.'),
            original_CPU_pareto_paths=old['pareto_paths'],
            wall_seconds_primary_unchanged=source['sets']['wall_seconds'], CPU_frontier_changed=differs))
    result = dict(schema='cpu-pareto-cost-sensitivity-v1', source_summary_sha256=SOURCE_SHA, stage=STAGE,
        cost_tolerance_seconds=COST_EPS, error_tolerance=ERROR_EPS, groups=groups,
        changed_CPU_frontier_indices=changed,
        scope='Postreview numerical-cost sensitivity, not a replacement or retuning of the predeclared primary analysis. '
              'All32groups/all11paths retained; CPU measurements and physical metrics unchanged. '
              'Unscored refers only to training objective; no physically validated baseline is filtered. '
              'Wall frontiers copied unchanged; no timing-repeatability guarantee. '
              'No resource-admission/ceiling tolerance is introduced.',
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output = (json.dumps(result, indent=2, allow_nan=False)+'\n').encode()
    need(len(output) <= 2*1024**2, 'bounded supplemental output')
    with a.output.open('xb') as f: f.write(output)
    print('complete32groups; changed_CPU_frontier_indices', changed)

if __name__ == '__main__': main()
