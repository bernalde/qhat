"""Freeze all constructions, then screen accuracy on the entire retained grid."""
import argparse
import copy
import itertools
from pathlib import Path
import shutil
import time
import numpy as np
from common import (HERE, ROOT, AUDIT, CASES, TIME, GRID, TARGETS, BUDGETS, REPETITIONS,
                    METHODS, REFERENCE_ARMS, ORDERING_MERGE, AMBIGUITY, PHASE_LIMITS,
                    require, read, write, sha, safe, sources, commit_sources, environment,
                    compiled, compiler)
from construction import construct_order, construct_candidate, validate_record


def artifact_hashes(root):
    paths = sorted(root.rglob('*'))
    require(not any(p.is_symlink() for p in paths), 'Artifact symlink refused')
    return {str(p.relative_to(root)): sha(p) for p in paths if p.is_file()}


def dispatch(case_index, repetitions):
    """Rotate method order across case/repetition without inspecting outcomes."""
    cells = []
    for rep in range(repetitions):
        for budget in BUDGETS:
            offset = (case_index + rep) % len(METHODS)
            for method in METHODS[offset:] + METHODS[:offset]:
                cells.append((method, budget, rep))
    return cells


def curve_name(case, method, budget, rep):
    return case + '.' + method + ('.reference' if budget is None else '.b' + str(int(budget)) + '.rep' + str(rep))


def make_protocol(output, preflight):
    require(not output.exists(), 'Use a fresh output; preserve every failed attempt')
    source_commit = commit_sources()
    output.mkdir(parents=True)
    old = AUDIT / 'symmetric_response/evidence'
    upstream, original_sources, provenance = compiled.authenticate_upstream(old, ORDERING_MERGE, False)
    shutil.copytree(old, output / 'upstream/evidence')
    for name, data in original_sources.items():
        target = output / 'upstream/source_snapshot' / safe(name)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    write(output / 'upstream/PROVENANCE.json', provenance)
    source_map = sources()
    for name in source_map:
        target = output / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    write(output / 'source_snapshot/MANIFEST.json', dict(source_commit=source_commit, sources=source_map))
    cases = CASES[-1:] if preflight else CASES
    repetitions = 1 if preflight else REPETITIONS
    curves = len(cases) * (len(METHODS) * len(BUDGETS) * repetitions + len(REFERENCE_ARMS))
    protocol = dict(source_commit=source_commit, sources=source_map, environment=environment(),
                    preflight=preflight, original_cases=cases, T=TIME, grid=GRID, targets=TARGETS,
                    matched_total_CPU_budgets=BUDGETS, repetitions=repetitions, methods=METHODS,
                    reference_arms=REFERENCE_ARMS, upstream=provenance,
                    ambiguity=AMBIGUITY, informational_error_floor=2e-12, compiler=compiler.VERSION,
                    phase_limits=PHASE_LIMITS, expected_curves=curves,
                    expected_candidates=curves * len(GRID), expected_target_rows=curves * len(TARGETS),
                    expected_comparisons=len(cases) * len(BUDGETS) * repetitions * len(TARGETS) * 17,
                    warmup_policy='one2secondcaporderingwarmuppermethodpermodel;allworkchargedtocampaignonly',
                    reference_cost_scope='Altered-model preflight tests mechanics only' if preflight else
                    'Inherited original standalone order cost; no matched-cap claim',
                    replay_precision=dict(dtype='clongdouble', mantissa_bits=int(np.finfo(np.longdouble).nmant),
                                          per_gate_renormalization=False))
    write(output / 'PROTOCOL_FROZEN.json', protocol)
    return protocol


def prepare(output, preflight=False):
    began, wall = time.process_time(), time.monotonic()
    protocol = make_protocol(output, preflight)
    order_actual = warm_actual = compilation_actual = 0.
    for ci, original_case in enumerate(protocol['original_cases']):
        model = compiled.original_model(original_case, preflight)
        case = model['case']
        model_path = output / 'models' / (case + '.json')
        write(model_path, model)
        digest = sha(model_path)
        for method in METHODS:
            warm = construct_order(model_path, digest, method, 2.)
            warm_actual += warm['ordering_CPU_seconds']
            write(output / 'warmups' / (case + '.' + method + '.json'), warm)
        for method, budget, rep in dispatch(ci, protocol['repetitions']):
            result = construct_order(model_path, digest, method, budget)
            order_actual += result['ordering_CPU_seconds']
            label = curve_name(case, method, budget, rep)
            validate_record(model, result['record'])
            write(output / 'curves' / (label + '.json'), dict(curve=label, case=case, method=method,
                  cohort='matched', total_budget=budget, repetition=rep, record=result['record'],
                  ordering=result, ordering_CPU_charge=result['ordering_CPU_seconds'], model_sha256=digest))
        design_path = output / 'upstream/evidence/results/designs' / (original_case + '.json')
        design = read(design_path)
        for method in REFERENCE_ARMS:
            records = [r for r in design['records'] if
                       (r['arm'], r['formula'], r['T'], r['budget_units']) == (method, 2, TIME, 30)]
            require(len(records) == 1, 'Missing strong frozen reference')
            record = copy.deepcopy(records[0])
            if preflight:
                for term in record['schedule']:
                    term[1] = str(7 * int(term[1]))
            validate_record(model, record)
            label = curve_name(case, method, None, 0)
            write(output / 'curves' / (label + '.json'), dict(curve=label, case=case, method=method,
                  cohort='inherited_reference', total_budget=None, repetition=0, record=record,
                  ordering_CPU_charge=record['standalone_construction_CPU_charge'],
                  upstream_design_sha256=sha(design_path), model_sha256=digest,
                  cost_scope=protocol['reference_cost_scope']))
    orders = {str(p.relative_to(output)): sha(p) for folder in ('models', 'curves', 'warmups')
              for p in sorted((output / folder).glob('*.json'))}
    require(len(list((output / 'curves').glob('*.json'))) == protocol['expected_curves'], 'Incomplete order cohort')
    write(output / 'ORDERS_FROZEN.json', dict(artifacts=orders, policy='All orders precede all new physical errors'))
    (output / 'traces').mkdir()
    complete = censored = eligible = 0
    for path in sorted((output / 'curves').glob('*.json')):
        curve = read(path)
        model_path = output / 'models' / (curve['case'] + '.json')
        for steps in GRID:
            label = curve['curve'] + '.r' + format(steps, '03d')
            trace = output / 'traces' / (label + '.ndjson.gz')
            identity = dict(candidate=label, curve=curve['curve'], case=curve['case'],
                            method=curve['method'], steps=steps, T=TIME, formula=2,
                            model_sha256=curve['model_sha256'], curve_sha256=sha(path),
                            protocol_sha256=sha(output / 'PROTOCOL_FROZEN.json'))
            result = construct_candidate(model_path, curve['model_sha256'], curve['record'], steps,
                                         trace, identity, curve['ordering_CPU_charge'], curve['total_budget'])
            result.update(candidate=label, curve=curve['curve'], case=curve['case'], method=curve['method'],
                          steps=steps, cohort=curve['cohort'], repetition=curve['repetition'],
                          trace_file=str(trace.relative_to(output)))
            compilation_actual += result['compilation_CPU_seconds']
            complete += result['status'] == 'complete'
            censored += result['status'] != 'complete'
            eligible += result['eligible']
            write(output / 'candidates' / (label + '.json'), result)
        print('FROZE_CONSTRUCTION_CURVE', curve['curve'], flush=True)
    require(complete + censored == protocol['expected_candidates'], 'Incomplete construction panel')
    write(output / 'PREPARE_COSTS.json', dict(CPU_seconds=time.process_time() - began,
          wall_seconds=time.monotonic() - wall, actual_ordering_CPU_seconds=order_actual,
          actual_warmup_CPU_seconds=warm_actual, actual_compilation_CPU_seconds=compilation_actual,
          repeated_standalone_charges_are_campaign_work=False,
          scope='Authentication, source/input copies, all actual warmups/orders/compilation/metadata; startup/finalfreezeinprocesslogs'))
    frozen = artifact_hashes(output)
    write(output / 'SCHEDULES_FROZEN.json', dict(artifacts=frozen, complete=complete, censored=censored,
          eligible=eligible, policy='All complete/partial exact primitive circuits precede all new validation errors'))


def bind(output, independent=False):
    protocol = read(output / 'PROTOCOL_FROZEN.json')
    frozen = read(output / 'SCHEDULES_FROZEN.json')
    for name, digest in frozen['artifacts'].items():
        require(sha(output / safe(name)) == digest, 'Frozen artifact changed: ' + name)
    manifest = read(output / 'source_snapshot/MANIFEST.json')
    require(manifest['sources'] == protocol['sources'] and manifest['source_commit'] == protocol['source_commit'],
            'Producer source archive identity differs')
    require(artifact_hashes(output / 'source_snapshot') ==
            {**protocol['sources'], 'MANIFEST.json': sha(output / 'source_snapshot/MANIFEST.json')},
            'Source archive hash/inventory differs')
    if not independent:
        require(sources() == protocol['sources'], 'Current producer differs from frozen source')
    require(protocol['environment'] == environment(), 'Numerical interpreter/package identity changed')
    require(protocol['grid'] == list(GRID) and protocol['targets'] == list(TARGETS) and
            protocol['ambiguity'] == AMBIGUITY and protocol['T'] == TIME, 'Accuracy policy differs')
    return protocol, frozen


def classify(error, target):
    return 'ambiguous' if abs(error - target) <= AMBIGUITY else 'met' if error < target else 'not_met'


def select(panel, rows, target):
    complete = [c for c in panel if c['status'] == 'complete']
    eligible = [c for c in complete if c['eligible']]
    qualifying = [rows[c['candidate']] for c in eligible if classify(rows[c['candidate']]['error'], target) == 'met']
    winner = min(qualifying, key=lambda r: (r['compiled_resources']['two_qubit_gates'],
                 r['compiled_resources']['unit_duration_ASAP_depth'], r['compiled_resources']['gates'], r['steps'])) if qualifying else None
    ambiguous = [c['candidate'] for c in eligible if classify(rows[c['candidate']]['error'], target) == 'ambiguous']
    return dict(status='met' if winner else 'ambiguous' if ambiguous else 'unmet',
                selected_candidate=winner['candidate'] if winner else None,
                qualifying_candidates=len(qualifying), ambiguous_candidates=ambiguous,
                complete_candidates=len(complete), censored_candidates=len(panel) - len(complete),
                construction_ineligible_candidates=sum(not c['eligible'] for c in complete))


def choices(output, errors=None):
    rows = {p.stem: read(p) for p in sorted((output / 'rows').glob('*.json'))}
    if errors is not None:
        for label, value in errors.items():
            rows[label]['error'] = value['error']
    result = []
    for path in sorted((output / 'curves').glob('*.json')):
        curve = read(path)
        panel = [read(output / 'candidates' / (curve['curve'] + '.r' + format(r, '03d') + '.json')) for r in GRID]
        for target in TARGETS:
            result.append(dict(curve=curve['curve'], case=curve['case'], method=curve['method'],
                               cohort=curve['cohort'], budget=curve['total_budget'], repetition=curve['repetition'],
                               target=target, **select(panel, rows, target)))
    return result


def compare_selections(output, selected):
    rows = {p.stem: read(p) for p in (output / 'rows').glob('*.json')}
    comparisons = []
    matched = [s for s in selected if s['cohort'] == 'matched']
    for left in matched:
        controls = list(REFERENCE_ARMS)
        if left['method'] in ('cubic_JD', 'cubic_warm_JD'):
            controls += ['first_JD', 'pauli_JD']
        if left['method'] == 'cubic_warm_JD':
            controls += ['cubic_JD']
        for control in controls:
            inherited = control in REFERENCE_ARMS
            right = next(s for s in selected if s['case'] == left['case'] and
                         s['budget'] == (None if inherited else left['budget']) and
                         s['repetition'] == (0 if inherited else left['repetition']) and
                         s['target'] == left['target'] and s['method'] == control)
            pair = dict(case=left['case'], budget=left['budget'], repetition=left['repetition'], target=left['target'],
                        method=left['method'], control=control, control_cohort=right['cohort'],
                        method_status=left['status'], control_status=right['status'])
            if left['status'] == right['status'] == 'met':
                a, b = rows[left['selected_candidate']], rows[right['selected_candidate']]
                for field in ('two_qubit_gates', 'unit_duration_ASAP_depth'):
                    x, y = a['compiled_resources'][field], b['compiled_resources'][field]
                    pair[field] = 'win' if x < y else 'loss' if x > y else 'tie'
                pair['method_selected'] = left['selected_candidate']
                pair['control_selected'] = right['selected_candidate']
            else:
                pair['two_qubit_gates'] = pair['unit_duration_ASAP_depth'] = 'unavailable'
            comparisons.append(pair)
    return comparisons


def run(output):
    began, wall = time.process_time(), time.monotonic()
    protocol, frozen = bind(output)
    require(not (output / 'rows').exists(), 'Preserve physical-phase attempts; fresh rows required')
    (output / 'rows').mkdir()
    (output / 'states').mkdir()
    from scipy.sparse.linalg import expm_multiply
    count = 0
    for path in sorted((output / 'models').glob('*.json')):
        model = read(path)
        context = compiled.PrimitiveReplay(model)
        cpu = time.process_time()
        initial = np.zeros(context.dimension, complex)
        initial[(1 << model['electrons']) - 1] = 1
        exact = expm_multiply(-1j * TIME * context.hamiltonian(), initial)
        reference_cost = time.process_time() - cpu
        for candidate_path in sorted((output / 'candidates').glob(model['case'] + '.*.json')):
            candidate = read(candidate_path)
            if candidate['status'] != 'complete':
                continue
            cpu = time.process_time()
            state = context.replay(output / candidate['trace_file'], candidate['identity'])
            replay_cost = time.process_time() - cpu
            values = compiled.metrics(exact, state, model)
            require(abs(values['norm'] - 1) < 1e-8, 'Compiled state norm drift')
            state_path = output / 'states' / (candidate['candidate'] + '.npz')
            np.savez_compressed(state_path, actual=state)
            row = {k: candidate[k] for k in ('candidate', 'curve', 'case', 'method', 'cohort', 'repetition',
                   'steps', 'compiled_resources', 'ordering_CPU_charge', 'compilation_CPU_seconds',
                   'total_construction_CPU_seconds', 'total_CPU_budget', 'eligible', 'CPU_overrun_seconds')}
            row.update(**values, error_at_or_below_numerical_floor=values['error'] <= 2e-12, reference_CPU_seconds=reference_cost, compiled_replay_CPU_seconds=replay_cost,
                       candidate_sha256=sha(candidate_path), trace_sha256=candidate['trace_sha256'],
                       state_file=str(state_path.relative_to(output)), state_sha256=sha(state_path))
            write(output / 'rows' / (candidate['candidate'] + '.json'), row)
            count += 1
        print('SCREENED_ALL_COMPLETE', model['case'], flush=True)
    require(count == frozen['complete'], 'Incomplete physical screening cohort')
    selected = choices(output)
    require(len(selected) == protocol['expected_target_rows'], 'Incomplete target panel')
    comparisons = compare_selections(output, selected)
    write(output / 'SELECTIONS.json', selected)
    write(output / 'COMPARISONS.json', comparisons)
    write(output / 'RUN_COMPLETE.json', dict(complete_candidates=count, target_rows=len(selected),
          comparisons=len(comparisons), schedule_freeze_sha256=sha(output / 'SCHEDULES_FROZEN.json')))
    write(output / 'RUN_COSTS.json', dict(CPU_seconds=time.process_time() - began,
          wall_seconds=time.monotonic() - wall, scope='Whole offline physical screening and selection;notstandaloneconstruction'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    if args.phase == 'prepare':
        prepare(args.output.resolve(), args.preflight)
    else:
        run(args.output.resolve())
