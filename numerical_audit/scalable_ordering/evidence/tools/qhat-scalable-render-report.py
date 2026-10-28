#!/usr/bin/env python3
"""Render a completed, independently verified sparse-ordering study (stdlib only).

Example:
  python /tmp/qhat-scalable-render-report.py --study STUDY \
    --verification INDEPENDENT_DIRECTORY --destination FRESH_REPORT_DIRECTORY

This is a derived-data renderer, not another scientific verification pass. It
does not import project code, propagate states, choose cases, or alter inputs.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys


ARMS = ('canonical', 'magnitude', 'current', 'first', 'pauli', 'cubic')
CASES = ('N2_r1p10', 'N2_r2p00', 'CO_r1p13', 'CH4_r1p09')
TIMES = (0.5, 1.5)
CAPS = (10, 30)
AMBIGUITY = 1e-12
MAX_JSON_BYTES = 32 * 1024**2
RENDERER_PATH = Path(__file__).resolve()
RENDERER_BYTES = RENDERER_PATH.read_bytes()
RENDERER_IDENTITY = dict(path=str(RENDERER_PATH), bytes=len(RENDERER_BYTES),
                         sha256=hashlib.sha256(RENDERER_BYTES).hexdigest())
INTERPRETER_IDENTITY = dict(executable=str(Path(sys.executable).resolve()), version=sys.version,
    implementation=sys.implementation.name, cache_tag=sys.implementation.cache_tag, platform=sys.platform)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, label, *, nonnegative=True):
    require(type(value) in (int, float) and math.isfinite(value) and
            (not nonnegative or value >= 0), 'Invalid finite number: ' + label)
    return value


def integer(value, label, *, minimum=0):
    require(type(value) is int and value >= minimum, 'Invalid integer: ' + label)
    return value


def digest(value):
    require(isinstance(value, str) and len(value) == 64 and
            all(c in '0123456789abcdef' for c in value), 'Invalid SHA-256 digest')
    return value


def unique_object(items):
    out = {}
    for key, value in items:
        require(key not in out, 'Duplicate JSON object key: ' + key)
        out[key] = value
    return out


def reject_constant(value):
    raise ValueError('Non-finite JSON constant: ' + value)


class Reader:
    def __init__(self, study, verification):
        self.study = study.resolve(strict=True)
        self.verification = verification.resolve(strict=True)
        self.files = {}

    def read(self, relative, *, verification=False):
        root = self.verification if verification else self.study
        relative = Path(relative)
        require(not relative.is_absolute() and '..' not in relative.parts, 'Unsafe input path')
        path = (root / relative).resolve(strict=True)
        require(path.is_relative_to(root) and path.is_file(), 'Input escaped its root')
        require(path.stat().st_size <= MAX_JSON_BYTES, 'JSON input exceeds byte bound')
        with path.open('rb') as stream:
            raw = stream.read(MAX_JSON_BYTES + 1)
        require(len(raw) <= MAX_JSON_BYTES, 'JSON input grew beyond byte bound')
        key = ('verification/' if verification else 'study/') + relative.as_posix()
        self.files[key] = dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        return json.loads(raw, object_pairs_hook=unique_object, parse_constant=reject_constant)

    def sha(self, relative):
        return self.files['study/' + str(relative)]['sha256']


def close_number(left, right, label, *, absolute=2e-9, relative=0):
    finite(left, label)
    finite(right, label)
    require(abs(left - right) <= absolute + relative * max(abs(left), abs(right)), label)


def decision(left, right, error_left, error_right):
    delta = error_left - error_right
    return left if delta < -AMBIGUITY else right if delta > AMBIGUITY else 'tie'


def ratio(numerator, denominator):
    """No invented zero/zero convention; every undefined observation is retained."""
    finite(numerator, 'ratio numerator')
    finite(denominator, 'ratio denominator')
    if denominator == 0:
        return dict(value=None, reason='both_zero' if numerator == 0 else 'zero_control')
    value = numerator / denominator
    if not math.isfinite(value):
        return dict(value=None, reason='nonfinite_quotient')
    return dict(value=value, reason=None)


def aggregate_ratios(comparisons, metric, *, above_error_floor=False):
    selected = [r for r in comparisons if r['status'] == 'complete' and
                (not above_error_floor or not r['at_numerical_floor'])]
    defined = [r['ratios'][metric]['value'] for r in selected
               if r['ratios'][metric]['value'] is not None]
    undefined = Counter(r['ratios'][metric]['reason'] for r in selected
                        if r['ratios'][metric]['value'] is None)
    return dict(median=statistics.median(defined) if defined else None,
                defined_count=len(defined), undefined_counts=dict(undefined),
                eligible_complete_count=len(selected), total_count=len(comparisons))


def resource_values(value, label):
    require(isinstance(value, dict), 'Missing compiled resources: ' + label)
    for key in ('gates', 'single_qubit_gates', 'two_qubit_gates',
                'unit_duration_ASAP_depth', 'CX_dependency_depth'):
        integer(value[key], label + '/' + key)
    require(value['gates'] == value['single_qubit_gates'] + value['two_qubit_gates'],
            'Inconsistent gate counts: ' + label)
    return value


def constructor_search_work(record, arm):
    """Preserve recorded counters; an allowance is never evidence of work."""
    stats = record.get('stats')
    require(stats is None or isinstance(stats, dict), 'Invalid constructor statistics')
    fields = ('seed_attempts', 'seed_evaluations', 'incremental_candidates', 'candidate_evaluations',
              'full_move_rescores', 'full_move_rescore_attempts')
    counts = {field: None if stats is None else stats.get(field) for field in fields}
    for field, value in counts.items():
        if value is not None:
            integer(value, 'constructor stats/' + field)
    if arm in ('canonical', 'magnitude', 'current'):
        require(all(value is None or value == 0 for value in counts.values()),
                'Unsearched baseline has nonzero search counters')
        counts = dict.fromkeys(fields, 0)
    return counts


def authenticate(reader):
    """Check completion, exact inventories, hashes and cross-file result parity."""
    frozen = reader.read('ALL_ORDERS_FROZEN.json')
    prior = reader.read('PRIOR_ALL_ORDERS_FROZEN.json')
    amendment = reader.read('SOURCE_AMENDMENT.json')
    require(amendment['status'] == 'completed' and
            amendment['prior_source_commit'] == prior['source_commit'] and
            amendment['new_source_commit'] == frozen['source_commit'] and
            amendment['prior_order_freeze_sha256'] == reader.sha('PRIOR_ALL_ORDERS_FROZEN.json') and
            amendment['new_order_freeze_sha256'] == reader.sha('ALL_ORDERS_FROZEN.json'),
            'Original training/execution source lineage differs')
    require(len(prior['constructors']) == 24 and
            prior['constructors'] == frozen['constructors'] == amendment['original_constructors'] and
            amendment['same_all_24_order_hashes'] is True and
            amendment['retraining'] is False and amendment['retuning'] is False and
            amendment['ordering_errors_evaluated_in_prior_attempt'] is False and
            amendment['prior_result_rows'] == 0, 'Frozen orders changed across packaging correction')
    require(amendment['changed_scientific_sources'] == ['numerical_audit/scalable_ordering/study.py'] and
            set(prior['sources']) == set(frozen['sources']) and
            {p for p in prior['sources'] if prior['sources'][p] != frozen['sources'][p]} ==
                {'numerical_audit/scalable_ordering/study.py'} and
            all(prior[key] == frozen[key] for key in ('environment', 'models', 'arms', 'train_time',
                'times', 'caps', 'train_CPU_seconds', 'CPU_seconds', 'wall_seconds',
                'input_manifest_sha256')), 'Packaging amendment changed scientific inputs/costs')
    run = reader.read('RUN_COMPLETE.json')
    summary = reader.read('SUMMARY.json')
    pairs = reader.read('PAIRS.json')
    independent = reader.read('INDEPENDENT_VERIFIED.json', verification=True)
    require(independent['status'] == 'passed', 'Independent verification has not passed')
    require(frozen['arms'] == list(ARMS) and frozen['times'] == list(TIMES) and
            frozen['caps'] == list(CAPS) and set(frozen['models']) == set(CASES), 'Frozen panel differs')
    require(frozen['expected_circuits'] == 96 and frozen['expected_references'] == 8 and
            frozen['train_time'] == 1.0, 'Frozen counts/training time differ')
    finite(frozen['train_CPU_seconds'], 'frozen search allowance')
    require(frozen['train_CPU_seconds'] > 0, 'Search allowance must be positive')
    require(run['order_freeze_sha256'] == reader.sha('ALL_ORDERS_FROZEN.json') and
            run['source_map'] == frozen['sources'] == independent['source_map'] and
            frozen['environment'] == independent['environment'], 'Freeze/source/environment binding differs')
    require(summary['circuits'] == 96 and summary['references'] == 8 and
            summary['all_pair_count'] == 240, 'Summary completion counts differ')
    require(len(independent['candidates']) == 96 and len(independent['references']) == 8 and
            len(independent['pairs']) == 240 and len(pairs) == 240, 'Independent/producer completion counts differ')

    chemistry = reader.read('chemistry/INPUTS_FROZEN.json')
    protocol = reader.read('chemistry/INPUT_PROTOCOL_FROZEN.json')
    require(not chemistry['refusals'] and chemistry['expected_cases'] == 4 and
            chemistry['successes'] == frozen['models'], 'Incomplete or different chemistry')
    require(frozen['input_manifest_sha256'] == reader.sha('chemistry/INPUTS_FROZEN.json') and
            chemistry['protocol_sha256'] == reader.sha('chemistry/INPUT_PROTOCOL_FROZEN.json'),
            'Chemistry protocol/manifest hash differs')
    registry = protocol['registry']
    specs = {r['case']: r for r in registry['cases']}
    require(len(registry['cases']) == 4 and set(specs) == set(CASES), 'Registry cases differ')
    require(registry['basis'] == 'sto-3g' and registry['unit'] == 'Angstrom', 'Chemistry basis/unit differs')
    geometry = protocol['geometry_identity']
    require(geometry['base'] == registry['base_commit'] and
            len(geometry['cases']) == 4 and {r['case'] for r in geometry['cases']} == set(CASES) and
            all(r['matches'] == [] for r in geometry['cases']), 'Recorded baseline geometry audit did not pass')

    expected_constructors = {f'constructors/{case}.{arm}.json' for case, arm in itertools.product(CASES, ARMS)}
    require(set(frozen['constructors']) == expected_constructors, 'Exact 24-constructor inventory differs')
    models, constructors, orders, model_tables = {}, {}, [], []
    for case in CASES:
        name = f'chemistry/inputs/{case}/model.json'
        model = reader.read(name)
        require(reader.sha(name) == digest(frozen['models'][case]), 'Model hash differs: ' + case)
        require(model['case'] == case and model['n'] == 12 and model['electrons'] == 6 and
                model['family'] == specs[case]['family'] and model['scf_converged'] is True,
                'Model case/active-space/convergence differs')
        require(all(model[key] == specs[case][key] for key in
                    ('geometry', 'frozen_spatial', 'active_spatial')), 'Model geometry/window differs')
        words, units = len(model['coefficients']), len(model['blocks'])
        require(0 < words <= 4096 and 0 < units <= 1024 and
                all(len(b) > 0 for b in model['blocks']) and sum(map(len, model['blocks'])) == words,
                'Invalid model word/unit count')
        models[case] = model
        model_tables.append(dict(case=case, family=model['family'], spin_modes=model['n'],
            active_electrons=model['electrons'], Pauli_words=words, conserving_units=units,
            generation_CPU_seconds=finite(model['generation_CPU_seconds'], 'chemistry CPU'),
            generation_wall_seconds=finite(model['generation_wall_seconds'], 'chemistry wall'),
            model_sha256=frozen['models'][case]))
        for arm in ARMS:
            name = f'constructors/{case}.{arm}.json'
            record = reader.read(name)
            require(reader.sha(name) == digest(frozen['constructors'][name]), 'Constructor hash differs: ' + name)
            unit_count = words if arm == 'pauli' else units
            require(record['arm'] == arm and record['model_sha256'] == frozen['models'][case] and
                    record['formula'] == 2 and record['time'] == frozen['train_time'] and
                    record['partition'] == ('pauli' if arm == 'pauli' else 'conserving') and
                    len(record['unit_lengths']) == unit_count and
                    all(type(n) is int and n > 0 for n in record['unit_lengths']) and
                    sum(record['unit_lengths']) == words and
                    sorted(record['order']) == list(range(unit_count)), 'Constructor layout differs: ' + name)
            require(isinstance(record['status'], str) and record['status'] and
                    record['response_is_error_bound'] is False, 'Invalid constructor status/bound flag')
            integer(record['evaluations'], 'constructor evaluations')
            integer(record['accepted_moves'], 'constructor accepted moves')
            objective = reader.read(f'objectives/{case}.{arm}.json')
            require(objective['kind'] == ('first' if arm in ('first', 'pauli') else 'cubic'),
                    'Independent objective kind differs')
            finite(objective['value'], 'independent objective')
            if record['score'] is not None:
                close_number(record['score'], objective['value'], 'Frozen objective differs',
                             absolute=1e-20, relative=1e-7)
            constructors[case, arm] = record
            orders.append(dict(case=case, arm=arm, partition=record['partition'], units=unit_count,
                CPU_seconds=finite(record['CPU_seconds'], 'constructor CPU'),
                wall_seconds=finite(record['wall_seconds'], 'constructor wall'),
                nominal_CPU_seconds=record['nominal_CPU_seconds'], status=record['status'],
                evaluations=record['evaluations'], accepted_moves=record['accepted_moves'],
                search_work=constructor_search_work(record, arm), raw_stats=record.get('stats'),
                completed_training_score=record['score'], independently_audited_objective=objective['value'],
                constructor_sha256=frozen['constructors'][name]))

    expected_rows = {f'{case}.{arm}.t{ti}.c{cap}' for case, arm, ti, cap in
                     itertools.product(CASES, ARMS, range(2), CAPS)}
    require(len(run['rows']) == 96 and set(run['rows']) == expected_rows, 'Exact 96-row completion inventory differs')
    expected_refs = {f'references/{case}.t{ti}.npz' for case, ti in itertools.product(CASES, range(2))}
    require(set(run['references']) == expected_refs, 'Exact eight-reference inventory differs')
    checks = {r['candidate']: r for r in independent['candidates']}
    require(set(checks) == expected_rows and len(checks) == 96, 'Independent candidate identities differ')
    independent_refs = {(r['case'], r['T']): r for r in independent['references']}
    require(len(independent_refs) == 8 and set(independent_refs) == set(itertools.product(CASES, TIMES)),
            'Independent reference identities differ')
    for case, ti in itertools.product(CASES, range(2)):
        reference = run['references'][f'references/{case}.t{ti}.npz']
        require(reference['case'] == case and reference['T'] == TIMES[ti], 'Reference logical identity differs')
        digest(reference['sha256'])
        finite(reference['CPU_seconds'], 'reference CPU')
        finite(reference['wall_seconds'], 'reference wall')
        difference = finite(independent_refs[case, TIMES[ti]]['max_difference'], 'independent reference difference')
        require(difference < 2e-9, 'Independent reference agreement did not pass')

    rows, statuses = {}, Counter()
    for case, arm, ti, cap in itertools.product(CASES, ARMS, range(2), CAPS):
        label = f'{case}.{arm}.t{ti}.c{cap}'
        row = reader.read(f'rows/{label}.json')
        record = constructors[case, arm]
        require((row['candidate'], row['case'], row['arm'], row['T'], row['cap_multiple']) ==
                (label, case, arm, TIMES[ti], cap), 'Candidate logical identity differs')
        require(row['status'] in ('complete', 'censored') and checks[label]['status'] == row['status'],
                'Candidate/independent status differs')
        require(row['ordering_CPU_seconds'] == record['CPU_seconds'] and
                row['ordering_status'] == record['status'] and row['ordering_evaluations'] == record['evaluations'],
                'Candidate constructor cost/status differs')
        require(row['identity'] == dict(case=case, arm=arm, model_sha256=frozen['models'][case],
                constructor_sha256=frozen['constructors'][f'constructors/{case}.{arm}.json'],
                T=TIMES[ti], steps=row['steps'], cap_multiple=cap, formula=2,
                preparation='HF_prefix', compilation_status=row['status']), 'Candidate trace identity differs')
        require(row['trace_file'] == f'traces/{label}.ndjson.gz' and
                row['reference_file'] == f'references/{case}.t{ti}.npz', 'Candidate trace/reference paths differ')
        digest(row['trace_sha256'])
        finite(row['CPU_seconds'], 'compile CPU')
        finite(row['wall_seconds'], 'compile wall')
        integer(row['steps'], 'symmetric steps', minimum=1)
        integer(row['ideal_rotations'], 'ideal rotations')
        require(row['rotation_cap'] == cap * len(models[case]['coefficients']), 'Rotation cap differs')
        if row['status'] == 'complete':
            require(row['ideal_rotations'] == row['expected_ideal_rotations'] <= row['rotation_cap'],
                    'Complete circuit rotation accounting differs')
            resource_values(row['compiled_resources'], label)
            require(checks[label]['resources'] == row['compiled_resources'], 'Independent actual resources differ')
            close_number(checks[label]['error'], row['error'], 'Independent state error differs')
            for field in ('error', 'squared_infidelity', 'norm', 'number_leakage', 'joint_sector_leakage'):
                finite(row[field], label + '/' + field)
            require(abs(row['norm'] - 1) < 1e-8, 'Complete circuit state norm differs')
            if record['partition'] == 'conserving':
                require(row['joint_sector_leakage'] < 1e-16, 'Conserving circuit leaked its sector')
            require(row['state_file'] == f'states/{label}.npz', 'Candidate state path differs')
            digest(row['state_sha256'])
            finite(row['validation_CPU_seconds'], 'primitive validation CPU')
        else:
            require(row['compiled_resources'] is None, 'Censored circuit has complete resources')
            resource_values(row['partial_compiled_resources'], label + '/partial')
        rows[case, TIMES[ti], cap, arm] = row
        statuses[row['status']] += 1
    require(summary['statuses'] == dict(statuses), 'Summary status counts differ')

    expected_pairs = {(case, T, cap, left, right) for case, T, cap in
                      itertools.product(CASES, TIMES, CAPS) for left, right in itertools.combinations(ARMS, 2)}
    producer_pairs = {(r['case'], r['T'], r['cap_multiple'], r['left'], r['right']): r for r in pairs}
    auditor_pairs = {(r['case'], r['T'], r['cap_multiple'], r['left'], r['right']): r
                     for r in independent['pairs']}
    require(len(producer_pairs) == len(auditor_pairs) == 240 and
            set(producer_pairs) == set(auditor_pairs) == expected_pairs, 'Exact 240-pair inventory differs')
    cubic_decisions = {control: Counter() for control in ARMS[:-1]}
    for key in expected_pairs:
        case, T, cap, left, right = key
        a, b = rows[case, T, cap, left], rows[case, T, cap, right]
        pair = producer_pairs[key]
        complete = a['status'] == b['status'] == 'complete'
        require(pair['status'] == ('complete' if complete else 'censored'), 'Producer pair status differs')
        outcome = decision(left, right, a['error'], b['error']) if complete else 'censored'
        require(pair.get('error_decision', 'censored') == outcome and
                auditor_pairs[key]['independent_error_decision'] == outcome, 'Paired decision differs')
        if complete:
            require(pair['error_left'] == a['error'] and pair['error_right'] == b['error'] and
                    pair['error_difference'] == a['error'] - b['error'] and
                    pair['CX_left'] == a['compiled_resources']['two_qubit_gates'] and
                    pair['CX_right'] == b['compiled_resources']['two_qubit_gates'] and
                    pair['depth_left'] == a['compiled_resources']['unit_duration_ASAP_depth'] and
                    pair['depth_right'] == b['compiled_resources']['unit_duration_ASAP_depth'] and
                    pair['construction_CPU_left'] == a['ordering_CPU_seconds'] + a['CPU_seconds'] and
                    pair['construction_CPU_right'] == b['ordering_CPU_seconds'] + b['CPU_seconds'],
                    'Producer paired metrics differ')
        if right == 'cubic':
            cubic_decisions[left][outcome] += 1
    require(summary['cubic_error_comparisons'] == {k: dict(v) for k, v in cubic_decisions.items()},
            'Summary cubic decisions differ')
    return frozen, run, summary, independent, protocol, models, model_tables, constructors, orders, rows, prior, amendment


def derive(reader):
    frozen, run, summary, independent, protocol, models, model_tables, constructors, orders, rows, prior, amendment = authenticate(reader)
    settings = []
    for case, T, cap in itertools.product(CASES, TIMES, CAPS):
        panel, complete = [], []
        for arm in ARMS:
            row = rows[case, T, cap, arm]
            resources = row['compiled_resources']
            item = dict(arm=arm, status=row['status'], reason=row['reason'], steps=row['steps'],
                ideal_rotations=row['ideal_rotations'], rotation_cap=row['rotation_cap'],
                ordering_status=row['ordering_status'], ordering_CPU_seconds=row['ordering_CPU_seconds'],
                compilation_CPU_seconds=row['CPU_seconds'],
                construction_CPU_seconds=row['ordering_CPU_seconds'] + row['CPU_seconds'],
                error=row.get('error'), CX=None if resources is None else resources['two_qubit_gates'],
                ASAP_depth=None if resources is None else resources['unit_duration_ASAP_depth'],
                CX_dependency_depth=None if resources is None else resources['CX_dependency_depth'],
                partial_compiled_resources=row['partial_compiled_resources'],
                squared_infidelity=row.get('squared_infidelity'), norm=row.get('norm'),
                number_leakage=row.get('number_leakage'), joint_sector_leakage=row.get('joint_sector_leakage'),
                row_file=f"rows/{row['candidate']}.json")
            panel.append(item)
            if row['status'] == 'complete':
                complete.append(item)
        best = min((r['error'] for r in complete), default=None)
        lowest = [] if best is None else [r['arm'] for r in complete if r['error'] <= best + AMBIGUITY]
        settings.append(dict(case=case, T=T, cap_multiple=cap, rows=panel,
            lowest_error_arms_among_completed=lowest,
            censored_arms=[r['arm'] for r in panel if r['status'] != 'complete'],
            comparison_ambiguity=AMBIGUITY))

    comparisons, aggregates = [], []
    for control in ARMS[:-1]:
        group = []
        for case, T, cap in itertools.product(CASES, TIMES, CAPS):
            a, b = rows[case, T, cap, 'cubic'], rows[case, T, cap, control]
            entry = dict(case=case, T=T, cap_multiple=cap, control=control,
                status='complete' if a['status'] == b['status'] == 'complete' else 'censored',
                cubic_status=a['status'], control_status=b['status'],
                cubic_reason=a['reason'], control_reason=b['reason'])
            if entry['status'] == 'complete':
                raw_decision = decision('cubic', control, a['error'], b['error'])
                outcome = 'win' if raw_decision == 'cubic' else 'tie' if raw_decision == 'tie' else 'loss'
                entry.update(outcome=outcome, cubic_error=a['error'], control_error=b['error'],
                    at_numerical_floor=a['error'] <= AMBIGUITY or b['error'] <= AMBIGUITY,
                    cubic_CX=a['compiled_resources']['two_qubit_gates'],
                    control_CX=b['compiled_resources']['two_qubit_gates'],
                    cubic_ASAP_depth=a['compiled_resources']['unit_duration_ASAP_depth'],
                    control_ASAP_depth=b['compiled_resources']['unit_duration_ASAP_depth'],
                    cubic_construction_CPU_seconds=a['ordering_CPU_seconds'] + a['CPU_seconds'],
                    control_construction_CPU_seconds=b['ordering_CPU_seconds'] + b['CPU_seconds'])
                entry['ratios'] = dict(error=ratio(entry['cubic_error'], entry['control_error']),
                    CX=ratio(entry['cubic_CX'], entry['control_CX']),
                    ASAP_depth=ratio(entry['cubic_ASAP_depth'], entry['control_ASAP_depth']),
                    construction_CPU=ratio(entry['cubic_construction_CPU_seconds'],
                                           entry['control_construction_CPU_seconds']))
            else:
                entry.update(outcome='censored', at_numerical_floor=None, ratios=None)
            group.append(entry)
        require(len(group) == 16, 'Control comparison count differs')
        counts = Counter(r['outcome'] for r in group)
        aggregates.append(dict(control=control, comparisons=16,
            wins=counts['win'], losses=counts['loss'], ties=counts['tie'], censored=counts['censored'],
            complete_pairs_at_numerical_floor=sum(r['at_numerical_floor'] is True for r in group),
            ratios={metric: aggregate_ratios(group, metric) for metric in
                    ('error', 'CX', 'ASAP_depth', 'construction_CPU')},
            error_ratio_above_floor=aggregate_ratios(group, 'error', above_error_floor=True)))
        comparisons.extend(group)

    all_rows = list(rows.values())
    costs = dict(
        chemistry_model_generation_CPU_seconds=sum(r['generation_CPU_seconds'] for r in model_tables),
        frozen_prepare_CPU_seconds=finite(frozen['CPU_seconds'], 'frozen prepare CPU'),
        frozen_prepare_wall_seconds=finite(frozen['wall_seconds'], 'frozen prepare wall'),
        standalone_constructor_CPU_seconds=sum(r['CPU_seconds'] for r in orders),
        compilation_CPU_seconds=sum(r['CPU_seconds'] for r in all_rows),
        primitive_state_validation_CPU_seconds=sum(r.get('validation_CPU_seconds', 0) for r in all_rows),
        reference_CPU_seconds=sum(r['CPU_seconds'] for r in run['references'].values()),
        independent_reference_CPU_seconds=sum(finite(r['CPU_seconds'], 'independent reference CPU')
                                              for r in independent['references']),
        scope='Recorded successful-stage charges only; not total task CPU. Failed attempts, imports/startup, '
              'geometry scan, external review and independent replay costs beyond the recorded references '
              'require the separate task ledger. Each constructor was charged once in this total.')
    return dict(schema='qhat-scalable-ordering-derived-report-v1',
        provenance=dict(study=str(reader.study), verification=str(reader.verification),
            source_commit=frozen['source_commit'], chemistry_source_commit=protocol['source_commit'],
            original_training_source_commit=prior['source_commit'],
            original_order_freeze_sha256=reader.sha('PRIOR_ALL_ORDERS_FROZEN.json'),
            effective_order_freeze_sha256=reader.sha('ALL_ORDERS_FROZEN.json'),
            source_amendment=amendment,
            scientific_source_sha256=frozen['sources'], consumed_files=reader.files,
            renderer=dict(RENDERER_IDENTITY), interpreter=dict(INTERPRETER_IDENTITY),
            statement='Derived JSON report; no states were propagated or scientific source imported by this renderer.'),
        completion=dict(orders=24, circuits=96, references=8, pairs=240,
                        independent_status=independent['status'], circuit_statuses=summary['statuses']),
        protocol=dict(train_time=frozen['train_time'], train_rotation_cap_multiple=30,
            validation_times=frozen['times'], validation_rotation_cap_multiples=frozen['caps'],
            nominal_search_CPU_seconds=frozen['train_CPU_seconds'], maximum_candidate_evaluations=3000,
            maximum_accepted_moves_per_start=5, declared_starts=['current', 'reverse_current', 'canonical'],
            metric='1-abs(overlap)', error_comparison_absolute_ambiguity=AMBIGUITY,
            geometry_identity=protocol['geometry_identity'], basis=protocol['registry']['basis'],
            cases_are_new_geometries=True, CH4_family_previously_recorded=True,
            compiler='Exact rational H/S/Sdg/Rz/CX primitives on a nearest-neighbor line; HF prefix included.',
            resources='One worker, quarter-core CPU, one native thread, 512 MiB RAM, no swap; '
                      'start above 5.5 GiB available and pause below 5 GiB; phase ledger owns total allowances.'),
        models=model_tables, constructors=orders, setting_tables=settings,
        cubic_vs_control_aggregates=aggregates, cubic_vs_control_comparisons=comparisons,
        recorded_CPU_costs=costs,
        limits=[
            'Four 12-mode, six-active-electron cases in three families. Two N2 geometries and repeated '
            'times/budgets are dependent observations; these counts are not independent sample sizes.',
            'The geometry audit concerns the pinned recorded baseline only. CH4 was already a baseline '
            'family; new geometries do not establish family generalization.',
            'All 24 orders were frozen before validation. Exact references and full 4096-state checking '
            'are offline validation costs; only training avoids Hilbert-space-sized arrays.',
            'Equal nominal search CPU limits do not imply equal consumed CPU, useful evaluations or '
            'construction cost. CPU-limited fallback orders, missing scores and censors remain reported.',
            'The cubic response is a heuristic, not a finite-error certificate. No globally best ordering '
            'or literature novelty is claimed.',
            'CX and depth are costs of this line compiler and its declared cleanup, including HF '
            'preparation; no native hardware, synthesis-optimality, correlated-state preparation, '
            'basis-error or SCF-error claim follows.',
            'Each ratio is computed within a matched case/time/cap before taking its median. Zero '
            'control values produce explicit undefined ratios and are excluded from that median. '
            'Numerical-floor pairs remain in outcome counts and are identified separately.',
            'Independent verification is numerical, not an interval or Lean certificate. This renderer '
            'checks the saved verification result and its inventories; it is not a formal Astra review.'
        ])


def number(value):
    return '—' if value is None else f'{value:.6g}'


def candidate_work(record):
    arm = record['arm']
    if arm in ('canonical', 'magnitude', 'current'):
        return '0 (unsearched)'
    field = 'incremental_candidates' if arm == 'cubic' else 'candidate_evaluations'
    return number(record['search_work'][field]) + (' incremental' if arm == 'cubic' else ' candidates')


def escaped(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def table(headers, rows):
    lines = ['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |']
    lines += ['| ' + ' | '.join(escaped(x) for x in row) + ' |' for row in rows]
    return '\n'.join(lines)


def render(data):
    lines = ['# Sparse molecular ordering: derived completed-study report', '',
        f"Independent verification: **{data['completion']['independent_status']}**; "
        '24 frozen orders, 96 declared circuits, eight references and 240 paired comparisons. '
        f"Circuit statuses: {json.dumps(data['completion']['circuit_statuses'], sort_keys=True)}.", '',
        'This panel evaluates four new recorded-baseline geometries from three families. '
        'The CH4 family already appeared in the baseline. Results below retain each setting and all '
        'cubic/control losses, ties and censors; they do not establish a universal winner or novelty.', '',
        f"Scientific source commit: `{data['provenance']['source_commit']}`. "
        'Every consumed input has an exact byte count and SHA-256 in `report_tables.json`.', '',
        f"The orders were originally constructed and frozen at `{data['provenance']['original_training_source_commit']}`. "
        'A row-packaging correction changed only study.py before successful execution. '
        'The renderer checks the original/effective freezes and amendment: all 24 constructor hashes '
        'are identical, with no retraining or retuning and no ordering-error result from the failed attempt. '
        'That attempt and its cost remain in the evidence archive.', '',
        '## Frozen protocol and scope', '',
        'Train at time 1.0 on HF occupations with a 30U ideal-rotation cap, where U is the retained '
        'Pauli-word count. Freeze all orders before validation at times 0.5 and 1.5 with caps 10U '
        'and 30U. The searched first-response, individual-Pauli and sparse-cubic arms have the same '
        f"nominal {number(data['protocol']['nominal_search_CPU_seconds'])} CPU-second allowance, "
        'up to 3,000 candidate evaluations and five accepted moves per declared start '
        '(current placement, its reverse and canonical). Actual costs and incomplete searches appear below. '
        'Canonical and magnitude are cheap baselines. All arms use the symmetric formula.', '',
        data['protocol']['resources'], '', data['protocol']['compiler'], '',
        'Reference propagation and independent full-space checking are offline validation. '
        'The response score is not a finite-time accuracy bound.', '',
        '## Molecules', '',
        table(['Case', 'Family', 'Spin modes', 'Active electrons', 'Pauli words U', 'Conserving units', 'Chemistry CPU s'],
              [[r['case'], r['family'], r['spin_modes'], r['active_electrons'], r['Pauli_words'],
                r['conserving_units'], number(r['generation_CPU_seconds'])] for r in data['models']]), '',
        'All cases use STO-3G and the registry\'s frozen-core/active windows. The two N2 geometries '
        'and repeated time/cap settings are dependent observations.', '',
        '## Standalone order construction', '',
        'CPU includes authentication, setup, initialization, search and returned layout. '
        'A null stored score means no completed training score; the separate final-order objective '
        'audit does not retroactively make that constructor a completed search. Completed seeds, '
        'candidate work and move rescoring below come from the recorded statistics. Cubic candidate '
        'work uses incremental_candidates; first/Pauli use candidate_evaluations. A dash means the '
        'counter is unknown, including unrecorded rescore attempts; it does not mean zero. Unsearched '
        'baselines have zero search counters, while their initializer work remains in CPU charges. '
        'Configured allowances alone do not establish completed search.', '',
        table(['Case', 'Arm', 'CPU s', 'Wall s', 'Status', 'Evaluations', 'Completed seeds',
               'Candidate work', 'Rescores completed/attempted', 'Accepted moves', 'Stored score'],
              [[r['case'], r['arm'], number(r['CPU_seconds']), number(r['wall_seconds']), r['status'],
                r['evaluations'], number(r['search_work']['seed_evaluations']), candidate_work(r),
                number(r['search_work']['full_move_rescores']) + '/' +
                    number(r['search_work']['full_move_rescore_attempts']),
                r['accepted_moves'], number(r['completed_training_score'])]
               for r in data['constructors']]), '',
        '## Cubic versus controls: all 16 matched settings per control', '',
        'A win means cubic has smaller state error by more than 1e-12; losses reverse that comparison. '
        'Differences within 1e-12 are ties. Ratios are cubic/control, formed within each setting before '
        'aggregation. Values below one favor cubic for the stated metric. Undefined zero-denominator '
        'ratios are excluded from medians and counted explicitly. Error ratios near the numerical '
        'floor are descriptive only.', '',
        table(['Control', 'Wins', 'Losses', 'Ties', 'Censors', 'Median error ratio (n)',
               'Median CX ratio (n)', 'Median depth ratio (n)', 'Median construction CPU ratio (n)'],
              [[r['control'], r['wins'], r['losses'], r['ties'], r['censored']] +
               [f"{number(r['ratios'][metric]['median'])} ({r['ratios'][metric]['defined_count']})"
                for metric in ('error', 'CX', 'ASAP_depth', 'construction_CPU')]
               for r in data['cubic_vs_control_aggregates']]), '',
        table(['Control', 'Complete floor pairs', 'Error ratios undefined', 'CX ratios undefined',
               'Depth ratios undefined', 'Median error ratio above floor (n)'],
              [[r['control'], r['complete_pairs_at_numerical_floor'],
                json.dumps(r['ratios']['error']['undefined_counts'], sort_keys=True),
                json.dumps(r['ratios']['CX']['undefined_counts'], sort_keys=True),
                json.dumps(r['ratios']['ASAP_depth']['undefined_counts'], sort_keys=True),
                f"{number(r['error_ratio_above_floor']['median'])} ({r['error_ratio_above_floor']['defined_count']})"]
               for r in data['cubic_vs_control_aggregates']]), '',
        'The construction CPU ratio adds standalone ordering and that circuit\'s compilation. '
        'The same order is reused at four validation settings; the campaign totals charge it once.', '',
        '## Raw six-arm settings', '',
        'Errors are 1-|overlap|. CX and ASAP depth come from actual compiled primitive DAGs, including '
        'HF preparation. The listed lowest-error set includes every completed arm within 1e-12 of '
        'the setting\'s minimum; a censor prevents a complete six-arm ranking. Partial circuits keep '
        'their reason and partial costs in `report_tables.json`.', '']
    for setting in data['setting_tables']:
        lines += [f"### {setting['case']}: time {setting['T']}, cap {setting['cap_multiple']}U", '',
            'Lowest error among completed arms (1e-12 ambiguity): ' +
            (', '.join(setting['lowest_error_arms_among_completed']) or 'none') + '. ' +
            'Censored arms: ' + (', '.join(setting['censored_arms']) or 'none') + '.', '',
            table(['Arm', 'Status', 'Steps', 'Ideal rotations', 'Error', 'CX', 'ASAP depth',
                   'Ordering CPU s', 'Compile CPU s'],
                  [[r['arm'], r['status'], r['steps'], r['ideal_rotations'], number(r['error']),
                    number(r['CX']), number(r['ASAP_depth']), number(r['ordering_CPU_seconds']),
                    number(r['compilation_CPU_seconds'])] for r in setting['rows']]), '']
    lines += ['## Every cubic/control outcome', '',
        table(['Case', 'Time', 'Cap U', 'Control', 'Outcome', 'Error ratio', 'CX ratio', 'Depth ratio', 'Floor pair'],
              [[r['case'], r['T'], r['cap_multiple'], r['control'], r['outcome']] +
               ([number(r['ratios'][metric]['value']) for metric in ('error', 'CX', 'ASAP_depth')]
                if r['status'] == 'complete' else ['—', '—', '—']) +
               ['yes' if r['at_numerical_floor'] else 'no' if r['at_numerical_floor'] is False else '—']
               for r in data['cubic_vs_control_comparisons']]), '',
        '## Recorded CPU charges', '',
        table(['Recorded component', 'CPU s'],
              [[key, number(value)] for key, value in data['recorded_CPU_costs'].items()
               if key.endswith('_CPU_seconds')]), '', data['recorded_CPU_costs']['scope'], '',
        'Frozen preparation CPU includes its constructor work and must not be added again to the '
        'standalone-constructor sum. No whole-task speed or memory claim is derived from these totals.', '',
        '## Limits', '']
    lines += ['- ' + text for text in data['limits']]
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--verification', type=Path, required=True,
                        help='Directory containing INDEPENDENT_VERIFIED.json')
    parser.add_argument('--destination', type=Path, required=True, help='Fresh directory outside project and inputs')
    parser.add_argument('--project-root', type=Path, default=Path.cwd(),
                        help='Checkout to keep free of report output (default: current directory; '
                             'set explicitly when launched elsewhere)')
    args = parser.parse_args()
    reader = Reader(args.study, args.verification)
    target = args.destination.resolve()
    require(not target.exists(), 'Fresh report destination required; old reports are preserved')
    for root in (reader.study, reader.verification, args.project_root.resolve()):
        require(not target.is_relative_to(root) and not root.is_relative_to(target),
                'Report destination must be outside project and input trees')
    data = derive(reader)
    markdown = render(data)
    payload = json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n'
    require(RENDERER_PATH.read_bytes() == RENDERER_BYTES, 'Renderer source changed during report derivation')
    target.mkdir(parents=True, exist_ok=False)
    require(RENDERER_PATH.read_bytes() == RENDERER_BYTES, 'Renderer source changed before Markdown write')
    (target / 'REPORT.md').write_text(markdown, encoding='utf-8')
    require(RENDERER_PATH.read_bytes() == RENDERER_BYTES, 'Renderer source changed before JSON write')
    (target / 'report_tables.json').write_text(payload, encoding='utf-8')
    print(json.dumps(dict(status='rendered', destination=str(target), orders=24, circuits=96,
                          references=8, pairs=240, consumed_files=len(reader.files)), sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (KeyError, TypeError, ValueError, OSError) as error:
        sys.exit(f'Report refused: {type(error).__name__}: {error}')
