"""Independent objective, primitive/ideal action, DAG and target-choice audit.

Current source reads frozen data; archived Python is never imported or executed.
Independent mathematical paths use literal Kronecker matrices and nested
commutators, not the producing bit-mask actions or Taylor recurrence.
"""
import argparse
from fractions import Fraction
import math
from pathlib import Path
import time
import numpy as np
import common as k
import study

STATE_TOL = 2e-10
METRIC_TOL = 2e-12


def inventory(path):
    files = list(path.rglob('*'))
    k.require(not any(p.is_symlink() for p in files), 'Symlink in delivery')
    return {str(p.relative_to(path)) for p in files if p.is_file()}


def close(a, b, message, atol=1e-12):
    k.require(math.isfinite(a) and math.isfinite(b) and abs(a - b) <= atol, message)


def audit_search(search, seeds, checked, m):
    k.require(len(search['starts']) == 3 and search['evaluation_cap'] == 3000 and
              search['max_moves_per_start'] == 5 and 0 < search['nominal_search_CPU_seconds'] <= 8.,
              'Search policy changed')
    attempted = completed = accepted = 0
    for start, seed in zip(search['starts'], seeds):
        k.require(start['initial_order'] == seed['order'] and start['initial_score'] == seed['score'], 'Search seed changed')
        k.require(0 <= start['completed'] <= start['attempted'] <= 1000, 'Search allowance exceeded')
        k.require(len(start['sweeps']) <= 5 and len(start['accepted']) <= 5, 'Move cap exceeded')
        k.require(sum(s['attempted'] for s in start['sweeps']) == start['attempted'] and
                  sum(s['completed'] for s in start['sweeps']) == start['completed'], 'Sweep work differs')
        for sweep in start['sweeps']:
            k.require(0 <= sweep['completed'] <= sweep['attempted'] <= m * (m - 1), 'Invalid sweep counts')
            k.require(sweep['attempted'] - sweep['completed'] <= 1, 'Multiple partial candidates')
            k.require(not sweep['complete'] or sweep['completed'] == sweep['attempted'] == m * (m - 1),
                      'False complete-sweep label')
        previous = seed
        for move in start['accepted']:
            source, destination = move['move']
            order = list(previous['order'])
            unit = order.pop(source)
            order.insert(destination, unit)
            k.require(order == move['record']['order'], 'Accepted insertion identity differs')
            checked(move['record'])
            k.require(move['record']['score'] < previous['score'] - max(1e-30, 1e-12 * abs(previous['score'])),
                      'Nonimproving insertion accepted')
            previous = move['record']
            accepted += 1
        k.require(previous == start['record'], 'Final search start differs from accepted history')
        checked(start['record'])
        attempted += start['attempted']
        completed += start['completed']
    best = min((s['record'] for s in search['starts']), key=lambda r: (r['score'], tuple(r['order'])))
    k.require(best == search['record'] and (attempted, completed) == (search['attempted'], search['completed']),
              'Final search winner/work differs')
    k.require(accepted <= search['accepted_rescores'] <= search['attempted_rescores'] <= accepted + 3,
              'Accepted rescore accounting differs')


def audit_objectives(curve, ref):
    construction = curve['ordering']
    phases = construction['phases']
    k.require(all(math.isfinite(x) and x >= 0 for x in phases.values()), 'Invalid order phase time')
    close(sum(phases.values()) + construction['accounting_overhead_CPU_seconds'],
          construction['ordering_CPU_seconds'], 'Order timing does not close')
    close(curve['ordering_CPU_charge'], construction['ordering_CPU_seconds'], 'Standalone ordering charge differs')
    k.require(construction['nominal_ordering_CPU_seconds'] == curve['total_budget'] / 2, 'Ordering fraction changed')
    cache = {}
    count = 0
    def audit_block(block, objective):
        nonlocal count
        def checked(record):
            nonlocal count
            key = (objective, tuple(record['order']))
            if key not in cache:
                cache[key] = (ref.score(record['order'], k.TIME, 30 * len(ref.model['coefficients']))
                              if objective == 'cubic' else ref.first_score(record['order'], k.TIME,
                              30 * len(ref.model['coefficients']), objective))
            other = cache[key]
            for field in ('JD2', 'score') if objective == 'cubic' else ('JD1', 'score'):
                k.require(np.isclose(record[field], other[field], rtol=2e-9, atol=5e-18),
                          'Independent objective value differs')
            k.require(record['steps'] == other['steps'], 'Objective affordability differs')
            if objective == 'cubic':
                k.require(np.max(abs(np.array(record['cubic_coefficients']) - other['cubic_coefficients'])) < 5e-11,
                          'Signed cubic coefficient differs')
            count += 1
        for seed in block['seeds']:
            checked(seed)
        if block['search'] is not None:
            m = len(ref.paulis) if objective == 'pauli' else len(ref.blocks)
            audit_search(block['search'], block['seeds'], checked, m)
            k.require(block['best'] == block['search']['record'], 'Search result differs from chosen score')
        elif block['seeds']:
            k.require(block['best'] == min(block['seeds'], key=lambda r: (r['score'], tuple(r['order']))),
                      'Best completed seed lost')
    if construction['warm_start'] is not None:
        audit_block(construction['warm_start'], 'conserving')
    audit_block(construction['result'], 'cubic' if curve['method'].startswith('cubic') else curve['record']['partition'])
    if construction['result']['best'] is not None:
        k.require(curve['record']['order'] == construction['result']['best']['order'], 'Circuit lost optimized order')
    return count


def independent_choices(curves, candidates, errors, costs):
    choices = []
    for name, curve in sorted(curves.items()):
        panel = [candidates[name + '.r' + format(r, '03d')] for r in k.GRID]
        for target in k.TARGETS:
            good = []
            ambiguous = []
            for candidate in panel:
                if candidate['status'] != 'complete' or not candidate['eligible']:
                    continue
                label = candidate['candidate']
                error = errors[label]['error']
                if abs(error - target) <= 2e-12:
                    ambiguous.append(label)
                elif error < target:
                    cost = costs[label]
                    good.append(((cost['two_qubit_gates'], cost['unit_duration_ASAP_depth'],
                                  cost['gates'], candidate['steps']), label))
            winner = min(good)[1] if good else None
            choices.append(dict(curve=name, case=curve['case'], method=curve['method'], cohort=curve['cohort'],
                  budget=curve['total_budget'], repetition=curve['repetition'], target=target,
                  status='met' if winner else 'ambiguous' if ambiguous else 'unmet', selected_candidate=winner,
                  qualifying_candidates=len(good), ambiguous_candidates=ambiguous,
                  complete_candidates=sum(c['status'] == 'complete' for c in panel),
                  censored_candidates=sum(c['status'] != 'complete' for c in panel),
                  construction_ineligible_candidates=sum(c['status'] == 'complete' and not c['eligible'] for c in panel)))
    return choices


def verify(output, receipt_dir=None, progress_dir=None):
    receipt = output if receipt_dir is None else receipt_dir
    if receipt_dir is None:
        k.require(not (output / 'VERIFIED.json').exists(), 'Preserve original verification; use a fresh receipt directory')
    else:
        k.require(not receipt.exists(), 'Use a fresh verification receipt directory')
        receipt.mkdir(parents=True)
    if progress_dir is not None:
        k.require(not progress_dir.exists(), 'Fresh independent progress directory required')
        progress_dir.mkdir(parents=True)
    progress = []
    began, wall = time.process_time(), time.monotonic()
    protocol, frozen = study.bind(output, independent=True)
    expected_cases = k.CASES[-1:] if protocol['preflight'] else k.CASES
    repetitions = 1 if protocol['preflight'] else 3
    k.require(protocol['original_cases'] == list(expected_cases) and protocol['repetitions'] == repetitions and
              protocol['methods'] == list(k.METHODS) and protocol['matched_total_CPU_budgets'] == [2., 8.] and
              protocol['reference_arms'] == list(k.REFERENCE_ARMS), 'Study domain differs')
    version = protocol.get('controller_policy_version', 1)
    k.require(version == 2 or (version == 1 and protocol['preflight']), 'Main requires fully rotated construction dispatch')
    expected_plan = []
    upstream = output / 'upstream/evidence'
    manifest = k.read(upstream / 'MANIFEST.json')
    k.require(inventory(upstream) == set(manifest['files']) | {'MANIFEST.json'}, 'Upstream inventory differs')
    k.require(k.sha(upstream / 'MANIFEST.json') == protocol['upstream']['manifest_sha256'], 'Upstream manifest differs')
    for name, digest in manifest['files'].items():
        k.require(k.sha(upstream / k.safe(name)) == digest, 'Upstream payload differs')
    k.require(protocol['upstream']['ordering_commit'] == k.ORDERING_MERGE and
              protocol['upstream']['payload_Git_verified'], 'Actual-merge dependency was not authenticated')
    curves = {p.stem: k.read(p) for p in sorted((output / 'curves').glob('*.json'))}
    candidates = {p.stem: k.read(p) for p in sorted((output / 'candidates').glob('*.json'))}
    rows = {p.stem: k.read(p) for p in sorted((output / 'rows').glob('*.json'))}
    expected_names = set()
    max_state = max_metric = 0.
    objective_count = 0
    errors, costs = {}, {}
    for ci, original_case in enumerate(expected_cases):
        original = k.read(upstream / 'results/inputs' / (original_case + '.json'))
        if protocol['preflight']:
            original['case'] = 'preflight_altered_' + original_case
            original['coefficients'] = [[w, str(7 * int(c))] for w, c in original['coefficients']]
            original['blocks'] = [[[w, str(7 * int(c))] for w, c in block] for block in original['blocks']]
            original['denominator'] = str(5 * int(original['denominator']))
            if 'removed_scalar_exact' in original:
                original['removed_scalar_exact'] = str(Fraction(original['removed_scalar_exact']) * Fraction(7, 5))
        case = original['case']
        model_path = output / 'models' / (case + '.json')
        model = k.read(model_path)
        k.require(model == original, 'Hamiltonian changed from authenticated input')
        literal = k.reference.Reference(model)
        objective_ref = k.objective_reference.Reference(model)
        objective_ref.conservation()
        exact = literal.exact_action(k.TIME)
        producer_H = k.compiled.PrimitiveReplay(model).hamiltonian()
        diff = producer_H - literal.H
        k.require(not diff.nnz or max(abs(diff.data)) < 2e-12, 'Full Hamiltonian action differs')
        cells = [(m, b, r) for m in k.METHODS for b in k.BUDGETS for r in range(repetitions)]
        cells += [(m, None, 0) for m in k.REFERENCE_ARMS]
        expected_plan.extend(study.curve_name(case, m, b, r) for m, b, r in study.dispatch(ci, repetitions))
        expected_plan.extend(study.curve_name(case, m, None, 0) for m in k.REFERENCE_ARMS)
        if version == 2:
            for method in k.METHODS:
                warm = k.read(output / 'warmups' / (case + '.' + method + '.json'))
                compiled_warm = warm['warmup_compilation']
                trace = output / k.safe(compiled_warm['trace_file'])
                k.require(k.sha(trace) == compiled_warm['trace_sha256'] and compiled_warm['status'] == 'complete', 'Warmup trace differs')
                state, cost = literal.audit_compiled_trace(trace, compiled_warm['identity'], k.compiler.VERSION)
                ideal = literal.ideal_formula(warm['record'], 3, k.TIME)
                k.require(np.max(abs(state - ideal)) < STATE_TOL and cost == compiled_warm['compiled_resources'], 'Full compiler warmup differs')
        for method, budget, rep in cells:
            name = study.curve_name(case, method, budget, rep)
            expected_names.add(name)
            curve = curves[name]
            k.require((curve['case'], curve['method'], curve['total_budget'], curve['repetition']) ==
                      (case, method, budget, rep) and curve['model_sha256'] == k.sha(model_path), 'Curve identity differs')
            study.validate_record(model, curve['record'])
            if budget is not None:
                objective_count += audit_objectives(curve, objective_ref)
            else:
                design_path = upstream / 'results/designs' / (original_case + '.json')
                records = [r for r in k.read(design_path)['records'] if
                           (r['arm'], r['formula'], r['T'], r['budget_units']) == (method, 2, k.TIME, 30)]
                k.require(len(records) == 1, 'Reference arm missing')
                record = records[0]
                if protocol['preflight']:
                    for term in record['schedule']:
                        term[1] = str(7 * int(term[1]))
                k.require(record == curve['record'] and curve['upstream_design_sha256'] == k.sha(design_path),
                          'Frozen strong reference differs')
                close(record['standalone_construction_CPU_charge'], curve['ordering_CPU_charge'], 'Inherited order charge differs')
            for steps in k.GRID:
                label = name + '.r' + format(steps, '03d')
                candidate = candidates[label]
                k.require(candidate['curve'] == name and candidate['steps'] == steps and
                          candidate['identity']['curve_sha256'] == k.sha(output / 'curves' / (name + '.json')) and
                          candidate['identity']['model_sha256'] == k.sha(model_path) and
                          candidate['identity']['protocol_sha256'] == k.sha(output / 'PROTOCOL_FROZEN.json'), 'Candidate identity differs')
                close(candidate['ordering_CPU_charge'], curve['ordering_CPU_charge'], 'Candidate order charge differs')
                charge = candidate['ordering_CPU_charge'] + candidate['compilation_CPU_seconds']
                close(charge, candidate['total_construction_CPU_seconds'], 'Total construction sum differs')
                close(max(0., charge - budget) if budget is not None else 0., candidate['CPU_overrun_seconds'], 'CPU overrun differs')
                eligible = candidate['status'] == 'complete' and (budget is None or charge <= budget)
                k.require(candidate['eligible'] == eligible and candidate['total_CPU_budget'] == budget, 'Construction eligibility differs')
                trace = output / k.safe(candidate['trace_file'])
                k.require(k.sha(trace) == candidate['trace_sha256'], 'Primitive trace changed')
                state, cost = literal.audit_compiled_trace(trace, candidate['identity'], k.compiler.VERSION,
                                                           simulate=candidate['status'] == 'complete')
                expected_cost = candidate['compiled_resources'] if candidate['status'] == 'complete' else candidate['partial_compiled_resources']
                k.require(cost == expected_cost, 'Independent circuit DAG/resource count differs')
                if candidate['status'] != 'complete':
                    k.require(label not in rows, 'Partial circuit was physically treated as complete')
                    continue
                costs[label] = cost
                ideal = literal.ideal_formula(curve['record'], steps, k.TIME)
                delta = float(np.max(abs(state - ideal)))
                stored = np.load(output / rows[label]['state_file'], allow_pickle=False)['actual']
                delta = max(delta, float(np.max(abs(stored - state))))
                k.require(delta < STATE_TOL, 'Literal primitive/unfused/producing state disagreement')
                max_state = max(max_state, delta)
                value = k.reference.metrics(exact, state, model)
                errors[label] = value
                for field in ('error', 'infidelity', 'norm', 'number_leakage', 'joint_sector_leakage'):
                    difference = abs(value[field] - rows[label][field])
                    k.require(difference < METRIC_TOL, 'Independent physical metric differs')
                    max_metric = max(max_metric, difference)
                k.require(rows[label]['candidate_sha256'] == k.sha(output / 'candidates' / (label + '.json')) and
                          rows[label]['state_sha256'] == k.sha(output / rows[label]['state_file']) and
                          rows[label]['trace_sha256'] == candidate['trace_sha256'], 'Physical row dependency differs')
            progress.append(name)
            if progress_dir is not None:
                k.write(progress_dir / (name + '.json'), dict(status='curve_checked', curve=name,
                        curve_sha256=k.sha(output / 'curves' / (name + '.json')),
                        schedule_freeze_sha256=k.sha(output / 'SCHEDULES_FROZEN.json'),
                        current_verifier_sources=k.sources(), checked_curves=progress,
                        checked_complete_candidates=len(errors), CPU_seconds=time.process_time() - began))
            print('INDEPENDENT_CURVE_CHECKED', name, flush=True)
        print('INDEPENDENT_PRACTICAL_MODEL', case, flush=True)
    if version == 2:
        k.require(k.read(output / 'CONSTRUCTION_PLAN_FROZEN.json') == expected_plan, 'Full construction dispatch differs')
        for index, name in enumerate(expected_plan):
            k.require(curves[name]['dispatch_index'] == index, 'Curve dispatch index differs')
    k.require(set(curves) == expected_names and len(curves) == protocol['expected_curves'], 'Curve cohort incomplete')
    expected_candidates = {name + '.r' + format(r, '03d') for name in expected_names for r in k.GRID}
    k.require(set(candidates) == expected_candidates and len(candidates) == protocol['expected_candidates'], 'Grid cohort differs')
    k.require(set(rows) == set(errors) and len(rows) == frozen['complete'], 'Physical row cohort differs')
    independent = independent_choices(curves, candidates, errors, costs)
    k.require(independent == k.read(output / 'SELECTIONS.json'), 'Independent target selection differs')
    comparisons = study.compare_selections(output, independent)
    k.require(comparisons == k.read(output / 'COMPARISONS.json') and len(comparisons) == protocol['expected_comparisons'],
              'Paired comparisons differ')
    report = dict(status='passed', curves=len(curves), candidates=len(candidates), complete=len(rows),
                  censored=frozen['censored'], target_rows=len(independent), comparisons=len(comparisons),
                  objective_records=objective_count, max_state_difference=max_state, max_metric_difference=max_metric,
                  schedule_freeze_sha256=k.sha(output / 'SCHEDULES_FROZEN.json'),
                  producer_sources=protocol['sources'], current_verifier_sources=k.sources(),
                  producer_environment=protocol['environment'], current_verifier_environment=k.environment(),
                  runtime_policy='Independent replay requires identical Python release and NumPy/SciPy versions;host build/path are separately recorded',
                  decision_agreement=True, CPU_seconds=time.process_time() - began,
                  wall_seconds=time.monotonic() - wall,
                  limitation='Data/math/trace checks do not certify timing measurements or hardware execution')
    k.write(receipt / 'INDEPENDENT_ERRORS.json', errors)
    k.write(receipt / 'INDEPENDENT_SELECTIONS.json', independent)
    physical = [p for folder in ('rows', 'states') for p in (output / folder).glob('*') if p.is_file()]
    physical += [output / n for n in ('SELECTIONS.json', 'COMPARISONS.json', 'RUN_COMPLETE.json', 'RUN_COSTS.json')]
    report['physical_artifacts'] = {str(p.relative_to(output)): k.sha(p) for p in sorted(physical)}
    report['independent_errors_sha256'] = k.sha(receipt / 'INDEPENDENT_ERRORS.json')
    report['independent_selections_sha256'] = k.sha(receipt / 'INDEPENDENT_SELECTIONS.json')
    k.write(receipt / 'VERIFIED.json', report)
    print('VERIFIED_PRACTICAL', {key: report[key] for key in ('curves', 'candidates', 'target_rows', 'comparisons')}, flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt-dir', type=Path)
    parser.add_argument('--progress-dir', type=Path)
    args = parser.parse_args()
    verify(args.output.resolve(), args.receipt_dir.resolve() if args.receipt_dir else None,
           args.progress_dir.resolve() if args.progress_dir else None)
