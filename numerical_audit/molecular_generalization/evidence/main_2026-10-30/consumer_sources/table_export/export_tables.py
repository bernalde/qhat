#!/usr/bin/env python3
"""Unexecuted metadata delivery for one explicitly native-censored main claim.

Usage: export_tables.py SEALED_DATA INCOMPLETE_ANALYSIS_JSON FRESH_OUTPUT
Complete CSV delivery does not complete or alter the frozen scientific gate.
Run only under root's bounded metadata allowance after every source job closes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import sys
import tomllib

PRODUCER = 'aea354e287f2277d8b11cc51695a6d49f7e84537'
FREEZE = 'f3ed635ac7831edcb3feb51d96e6116c2c95f7ad9e2b169b006ca862006776e3'
CASES = ('HF_r1p00_631g_12m8e', 'H2S_r1p40_a96_631g_12m8e')
ARMS = (('parent', 0), ('canonical', 0), ('magnitude', 0), ('full', 32),
        ('full', 64), ('window', 32), ('window', 64), ('annealing', 32), ('annealing', 64))
OPTIMIZED = ('full', 'window', 'annealing')
QUESTIONS = ('window_minus_full', 'annealing_minus_full')
COSTS = ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds')
STUDY = 'source_snapshot/numerical_audit/molecular_generalization/'
MAX_FILES, MAX_BYTES, MAX_JSON = 20000, 512*1024**2, 32*1024**2
CENSORED_CASE = CASES[0]
CENSORED_NAME = 'rep2.window.B64.raw1'
NATIVE_REFUSAL = dict(type='SearchDeadline', message='whole-call window allowance exhausted')
EXPECTED_ISSUES = [dict(phase='scientific_receipts', type='ValueError',
    message='Scientific post-teardown FINAL/PLAN/receipt chain')]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def hash_spelling(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def canonical(value):
    return json.dumps(value, ensure_ascii=True, allow_nan=False, sort_keys=True,
                      separators=(',', ':')).encode('ascii')


def unique_object(items):
    value = {}
    for key, item in items:
        require(key not in value, 'Duplicate JSON key: '+key)
        value[key] = item
    return value


def finite_float(text):
    value = float(text)
    require(math.isfinite(value), 'Nonfinite JSON number')
    return value


def load_json(raw):
    def reject(text):
        raise ValueError('Nonfinite JSON constant: '+text)
    return json.loads(raw, object_pairs_hook=unique_object,
                      parse_float=finite_float, parse_constant=reject)


def relative(name):
    require(isinstance(name, str) and 0 < len(name) <= 4096 and '\\' not in name
            and all(ord(c) >= 32 and ord(c) != 127 for c in name)
            and not name.startswith('/') and len(name.split('/')) <= 64
            and all(p not in ('', '.', '..') for p in name.split('/')), 'Canonical relative path')
    return name


def stream_file(root, name, *, maximum=MAX_BYTES, capture=False):
    """Descriptor-relative no-symlink walk; bounded regular files only."""
    parts = relative(name).split('/')
    directory = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        for component in parts[:-1]:
            child = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory)
            os.close(directory)
            directory = child
        descriptor = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory)
        with os.fdopen(descriptor, 'rb') as stream:
            before = os.fstat(stream.fileno())
            require(stat.S_ISREG(before.st_mode) and 0 <= before.st_size <= maximum, 'Bounded regular input: '+name)
            hashed, count, chunks = hashlib.sha256(), 0, []
            while True:
                chunk = stream.read(min(1024**2, maximum-count+1))
                if not chunk:
                    break
                count += len(chunk)
                require(count <= maximum, 'Growing input: '+name)
                hashed.update(chunk)
                if capture:
                    chunks.append(chunk)
            after = os.fstat(stream.fileno())
            require(count == before.st_size == after.st_size and
                    before.st_mtime_ns == after.st_mtime_ns, 'Input changed during read: '+name)
            return dict(bytes=count, sha256=hashed.hexdigest()), b''.join(chunks) if capture else None
    finally:
        os.close(directory)


class Evidence:
    def __init__(self, root):
        self.root = root
        self.inventory_identity, raw = stream_file(root, 'EVIDENCE_INVENTORY.json', maximum=MAX_JSON, capture=True)
        self.inventory = load_json(raw)
        require(self.inventory.get('schema') == 'molecular-evidence-inventory-v1'
                and self.inventory.get('source_freeze_sha256') == FREEZE, 'Sealed main inventory identity')
        self.files = self.inventory['files']
        require(isinstance(self.files, dict) and 0 < len(self.files) <= MAX_FILES, 'Inventory file bound')
        total = 0
        for name, item in self.files.items():
            relative(name)
            require(isinstance(item, dict) and set(item) == {'bytes', 'sha256'}
                    and type(item['bytes']) is int and item['bytes'] >= 0
                    and hash_spelling(item['sha256']), 'Inventory member declaration: '+name)
            total += item['bytes']
        require(total == self.inventory['uncompressed_bytes'] and total <= MAX_BYTES, 'Inventory byte bound')
        actual_names = set()
        for directory, directories, files in os.walk(root, followlinks=False):
            for name in directories+files:
                require(not (Path(directory)/name).is_symlink(), 'Symlink in sealed DATA')
            for name in files:
                path = Path(directory)/name
                require(path.is_file(), 'Nonregular sealed artifact')
                actual_names.add(path.relative_to(root).as_posix())
                require(len(actual_names) <= MAX_FILES+1, 'Observed sealed file bound')
        require(actual_names == set(self.files) | {'EVIDENCE_INVENTORY.json'}, 'Extra/missing sealed DATA artifact')
        for name, expected in self.files.items():
            actual, _ = stream_file(root, name)
            require(actual == expected, 'Sealed artifact bytes differ: '+name)

    def raw(self, name, *, maximum=MAX_JSON):
        require(name in self.files, 'Consumed artifact absent from seal: '+name)
        actual, raw = stream_file(self.root, name, maximum=maximum, capture=True)
        require(actual == self.files[name], 'Consumed artifact differs from seal: '+name)
        return raw

    def json(self, name):
        return load_json(self.raw(name))

    def toml(self, name):
        return tomllib.loads(self.raw(name).decode('utf-8'))

    def bind(self, name, expected):
        require(hash_spelling(expected) and name in self.files and self.files[name]['sha256'] == expected,
                'Actual-byte hash reference differs: '+name)


def number(value):
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0, 'Finite nonnegative recorded number')
    return value


def identity(row):
    require(isinstance(row.get('case'), str) and isinstance(row.get('arm'), str)
            and type(row.get('repetition')) is int and type(row.get('cpu_budget')) is int,
            'Typed case/repetition/arm/budget identity')
    return tuple(row[k] for k in ('case', 'repetition', 'arm', 'cpu_budget'))


def logical(row):
    require(type(row.get('time')) is float and row['time'] in (.65, 1.35)
            and type(row.get('budget_multiplier')) is int and row['budget_multiplier'] in (8, 24),
            'Typed frozen physical condition')
    return identity(row)+(row['time'], row['budget_multiplier'])


def index(rows, key, expected=None):
    require(isinstance(rows, list), 'List of records required')
    result = {}
    for row in rows:
        require(isinstance(row, dict), 'Record object required')
        name = key(row)
        require(name not in result, 'Duplicate record identity: '+str(name))
        result[name] = row
    require(expected is None or set(result) == expected, 'Complete exact record roster')
    return result


def require_same_rows(recorded, actual, key):
    left, right = index(recorded, key), index(actual, key)
    require(set(left) == set(right) and all(canonical(left[k]) == canonical(right[k]) for k in left),
            'ANALYSIS retained rows differ from authenticated report bytes')


def costs(value):
    require(isinstance(value, dict) and set(value) == set(COSTS), 'Complete cost fields')
    for item in value.values():
        number(item)
    return value


def cost_columns(**groups):
    return {group+'_'+key: value[key] for group, value in groups.items() for key in COSTS}


def order_hash(order, count, orders, *, one_based=False):
    require(isinstance(order, list) and all(type(x) is int for x in order), 'Integer original-label order')
    values = [x-1 for x in order] if one_based else list(order)
    require(sorted(values) == list(range(count)), 'Complete original-label permutation')
    raw = canonical(values)
    sha = digest(raw)
    require(sha not in orders or orders[sha] == values, 'Conflicting canonical order digest')
    orders[sha] = values
    return sha, values


def field(value, *path):
    for name in path:
        if not isinstance(value, dict) or name not in value:
            return None
        value = value[name]
    return value


def verified_comparison(value):
    require(isinstance(value, dict) and value.get('verified') is True
            and value.get('within_tolerance') is True and value.get('numerical_floor') is False,
            'Recorded noncensored proposed/native comparison must pass')


def known_native_censor(row):
    require((row['case'], row['name'], row['order_id'], row['repetition'], row['arm'], row['cpu_budget']) ==
            (CENSORED_CASE, CENSORED_NAME, 'rep2.window.B64', 2, 'window', 64)
            and row.get('selected') is False and row.get('role') == 'raw_proposed_best'
            and row.get('status') == 'native_refusal_censored', 'Only declared unpromoted native censor')
    require('native_score' not in row and 'native_comparison' not in row,
            'Original missing native result must never be backfilled')
    refusal = row.get('native_refusal')
    parsed = load_json(refusal.encode('utf-8')) if isinstance(refusal, str) else refusal
    require(parsed == NATIVE_REFUSAL, 'Exact original native deadline refusal')
    verified_comparison(row['proposed_comparison'])


def receipt_records(e, freeze):
    """Retain all six jobs, including exit2; do not relabel their source gate."""
    expected = {(phase, case) for phase in ('physical', 'independent', 'objective') for case in CASES}
    found = {}
    for name in sorted(e.files):
        if not (name.startswith('jobs/') and name.endswith('.PLAN.json')):
            continue
        plan = e.json(name)
        phase = plan.get('phase')
        if phase not in ('physical', 'independent', 'objective'):
            continue
        key = (phase, plan['identity']['case'])
        require(key in expected and key not in found, 'Exactly one original scientific phase job per case')
        label = plan['label']; final_path = 'jobs/'+label+'.FINAL.json'; receipt_path = 'logs/'+label+'.resources.json'
        final, receipt = e.json(final_path), e.json(receipt_path)
        require(name == 'jobs/'+label+'.PLAN.json' and plan['schema'] == 'molecular-job-plan-v1'
                and plan['source_freeze_sha256'] == FREEZE and plan['producer_root'] == freeze['producer_root']
                and plan['producer_data_root'] == freeze['producer_data_root'], 'Original phase PLAN identity')
        exception = key == ('objective', CENSORED_CASE)
        exit_code = 2 if exception else 0
        require(final['schema'] == 'molecular-job-final-v1' and final['phase'] == phase and final['label'] == label
                and final['plan_sha256'] == e.files[name]['sha256'] and final['receipt_path'] == receipt_path
                and final['receipt_sha256'] == e.files[receipt_path]['sha256']
                and final['guard_supervisor_exit_code'] == exit_code and receipt['exit_code'] == exit_code,
                'Only known HF objective exit2; other scientific jobs must pass')
        require(receipt['command'] == plan['command'] and receipt['budget'] == dict(cpu_seconds=plan['cpu_cap'],
                wall_seconds=plan['wall_cap'], disk_mib=plan['disk_mib'])
                and receipt.get('termination') is None and receipt.get('guard_error') is None
                and not receipt.get('cleanup_errors'), 'Closed original scientific command/cost/termination receipt')
        expected_settings = dict(CPUQuota='25%', CPUWeight='5', IOWeight='10', Nice='19',
                                 MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0')
        require(all(str(receipt['settings'][k]) == v for k, v in expected_settings.items()), 'Original scientific OS profile')
        script = freeze['producer_root']+'/numerical_audit/molecular_generalization/'+(
            'verify.jl' if phase == 'independent' else phase+'.jl')
        command = plan['command']; require(command.count(script) == 1, 'Original numerical script')
        position = command.index(script); case = key[1]; base = freeze['producer_data_root']
        arguments = [f'{base}/inputs/{case}/model.toml', f'{base}/panels/{case}.toml']
        arguments += [f'{base}/physical/{case}', f'{base}/independent/{case}'] if phase == 'independent' else [f'{base}/{phase}/{case}']
        if phase == 'objective':
            arguments.append(f'{base}/witnesses/{case}.toml')
        require(command[position+1:] == arguments and '--threads=1' in command[:position]
                and '--project='+freeze['producer_root']+'/analysis/julia_window' in command[:position]
                and command[0] == freeze['runtime']['julia']['executable'], 'Actual frozen scientific invocation')
        e.bind(STUDY+'tools/resource_guard.py', plan['guard_sha256'])
        pm = receipt['process_measurements']; raw = number(pm['user_cpu_seconds'])+number(pm['system_cpu_seconds'])
        sampled = []
        for sample in receipt.get('samples', []):
            for line in sample.get('service', '').splitlines():
                if line.startswith('CPUUsageNSec=') and line.split('=', 1)[1].isdigit():
                    sampled.append(int(line.split('=', 1)[1])/1e9)
        measured = dict(raw_CPU_seconds=raw, charged_CPU_seconds=max([raw+.011]+sampled),
                        wall_seconds=number(receipt['elapsed']))
        require(all(final[k] == v for k, v in measured.items()) and
                measured['charged_CPU_seconds'] <= number(plan['cpu_cap']) and
                measured['wall_seconds'] <= number(plan['wall_cap']), 'Original finite phase cost fields/caps')
        found[key] = dict(phase=phase, case=case, plan_path=name, final_path=final_path, receipt_path=receipt_path,
            plan_sha256=e.files[name]['sha256'], final_sha256=e.files[final_path]['sha256'],
            receipt_sha256=e.files[receipt_path]['sha256'], exit_code=exit_code,
            original_scientific_receipt_pass=not exception, known_native_censor_exit2=exception,
            source_final=final, source_receipt=receipt, **measured)
    require(set(found) == expected, 'All six original scientific jobs must be closed and retained')
    return [found[key] for key in sorted(found)]


def censored_report_scope(a, reports, case_metadata):
    require(a['issues'] == EXPECTED_ISSUES, 'Only exact source scientific_receipts issue is allowed')
    require(all(value['decision'] == 'incomplete' for value in a['primary'].values()),
            'Both original primary decisions must remain incomplete')
    raw = a['retained_raw_witnesses']
    censored = [row for row in raw if row.get('status') != 'verified']
    require(len(censored) == 1, 'Exactly one original raw native censor; no additional gaps')
    known_native_censor(censored[0])
    for case in CASES:
        is_censored = case == CENSORED_CASE
        metadata = case_metadata[case]
        require(metadata['complete'] is (not is_censored) and metadata['raw_claims_verified'] is (not is_censored)
                and metadata['issues'] == [], 'Only known case-level raw-censor incompleteness')
        physical, independent, objective = (reports[case][phase] for phase in ('physical', 'independent', 'objective'))
        require(physical['status'] == independent['status'] == 'complete' and physical['physical_panel_complete'] is True
                and physical['all_constructor_eligible'] is True and independent['all_constructor_eligible'] is True
                and physical['completed_logical_rows'] == independent['completed_logical_rows'] == 72
                and independent['reference_check_count'] == 2 and not independent.get('censored_rows')
                and not physical.get('failures') and not independent.get('failures'), 'All physical and independent rows/references must pass')
        require(physical['unique_actual_propagations'] == independent['unique_circuit_count'] ==
                independent['independent_unfused_replay_count'] == metadata['unique_circuits'], 'All distinct circuits independently replayed')
        require(objective['status'] == ('partial_or_censored' if is_censored else 'complete')
                and objective['complete'] is (not is_censored)
                and objective['all_scored_objectives_verified'] is (not is_censored)
                and objective['processed_logical_rows'] == objective['logical_rows_retained'] == 72
                and objective['witness_inventory_complete'] is True
                and objective['expected_raw_witnesses'] == objective['processed_raw_witnesses'] == len(objective['raw_witnesses'])
                and objective['mismatched_objectives'] == objective['floor_censored_objectives'] == 0
                and objective['native_refused_objectives'] == int(is_censored)
                and not objective.get('failures'), 'Only original native censor in complete objective inventory')
        require(all(check['status'] == 'replayed' for check in objective['checks']), 'Every archived objective check must replay')
        if is_censored:
            require(len(objective['checks']) == 14 and len(objective['raw_witnesses']) == 7,
                    'Exact observed HF 14-check/seven-raw-claim scope')
        for row in objective['logical_rows']:
            require(row['status'] == 'verified' and 'native_refusal' not in row, 'All 144 selected training score checks must pass')
            verified_comparison(row['proposed_comparison']); verified_comparison(row['native_comparison'])
        for row in objective['raw_witnesses']:
            verified_comparison(row['proposed_comparison'])
            if row['status'] == 'verified':
                require('native_refusal' not in row, 'Unexpected native refusal')
                verified_comparison(row['native_comparison'])
            else:
                known_native_censor(dict(case=case, **row))
        for reference in physical['reference_records']:
            require(reference['status'] == 'complete', 'Complete saved reference')
        for reference in independent['reference_checks']:
            require(reference['status'] == 'complete', 'Complete independent reference')
    return censored[0]


def export(data, analysis_path, output):
    require(not output.exists(), 'Fresh output directory required')
    require(not output.is_relative_to(data), 'Output must be outside sealed DATA')
    for ancestor in output.parents:
        require(not any((ancestor/name).exists() for name in
            ('EVIDENCE_INVENTORY.json', 'INPUT_PROTOCOL_FROZEN.json', 'EXTRACTION_VERIFIED.json')),
            'Output ancestor is a sealed/frozen evidence tree')
    e = Evidence(data)
    freeze = e.json('INPUT_PROTOCOL_FROZEN.json')
    require(e.files['INPUT_PROTOCOL_FROZEN.json']['sha256'] == FREEZE and
            freeze['source_commit'] == PRODUCER and freeze['mode'] == 'main', 'Exact producer/freeze')
    require(not output.is_relative_to(Path(freeze['producer_root']).resolve()), 'Output must be outside producer source')
    analysis_identity, raw = stream_file(analysis_path.parent, analysis_path.name, maximum=MAX_JSON, capture=True)
    analysis_sealed_members = [name for name, item in e.files.items() if item == analysis_identity]
    require(analysis_sealed_members, 'ANALYSIS must be a sealed member or its exact byte copy')
    if analysis_path.is_relative_to(data):
        e.bind(analysis_path.relative_to(data).as_posix(), analysis_identity['sha256'])
    a = load_json(raw)
    require(a.get('schema') == 'molecular-analysis-v1' and a.get('complete') is False
            and a.get('status') == 'incomplete' and a.get('issues') == EXPECTED_ISSUES and a.get('mode') == 'main'
            and a.get('cases') == list(CASES) and a.get('repetitions') == 2
            and a.get('confidence_intervals') is False, 'Exact native-censored frozen main ANALYSIS required')
    require(a['source_freeze_sha256'] == FREEZE and a['original_frozen_sources'] == freeze['sources']
            and a['current_consumer_sources'] == freeze['sources'], 'ANALYSIS source-freeze/consumer identity')
    for name, sha in freeze['sources'].items():
        e.bind('source_snapshot/'+relative(name), sha)
    require(isinstance(a['artifact_hashes'], dict) and a['artifact_hashes'], 'ANALYSIS byte graph')
    for name, sha in a['artifact_hashes'].items():
        e.bind(relative(name), sha)
    e.bind(STUDY+'PROTOCOL.json', a['protocol_sha256'])
    config = e.json(STUDY+'PROTOCOL.json')
    require(config == freeze['configuration'] and config['main_execution_allowed'] is True
            and config['protocol_state'] == 'frozen_main' and config['inventory']['total_order_records'] == 36
            and config['inventory']['logical_physical_rows'] == 144
            and config['physical']['primary']['questions'] == list(QUESTIONS), 'Frozen main configuration')
    inputs = e.json('INPUTS_FROZEN.json')
    require(inputs['status'] == 'complete' and inputs['refusals'] == {} and inputs['cases'] == list(CASES)
            and inputs['repetitions'] == 2 and inputs['source_freeze_sha256'] == FREEZE
            and inputs['protocol_sha256'] == a['protocol_sha256'] and set(inputs['successes']) == set(CASES),
            'Complete input identity')
    sizes = {}
    for case in CASES:
        for name, sha in inputs['successes'][case].items():
            e.bind('inputs/'+case+'/'+relative(name), sha)
        model = e.json('inputs/'+case+'/model.json')
        require(model['case'] == case and model['n'] == 12 and model['electrons'] == 8,
                'Frozen active space')
        sizes[case] = len(model['blocks'])
    expected = {(case, rep, arm, budget) for case in CASES for rep in (1, 2) for arm, budget in ARMS}
    grid = {key+(time, cap) for key in expected for time in (.65, 1.35) for cap in (8, 24)}
    crows = index(a['constructor_costs']['constructors'], identity, expected)
    require(a['constructor_costs']['inventory_complete'] is True and len(crows) == 36, '36 complete constructor costs')
    physical = index(a['retained_physical_rows'], logical, grid)
    independent = index(a['retained_independent_rows'], logical, grid)
    objective = index(a['retained_objective_rows'], logical, grid)
    require(len(physical) == len(independent) == len(objective) == 144, '144/144/144 row denominators')
    case_metadata = index(a['cases_metadata'], lambda row: row['case'], set(CASES))
    reports = {}
    for case in CASES:
        reports[case] = {phase: e.toml(f'{phase}/{case}/{filename}') for phase, filename in
            (('physical', 'PHYSICAL.toml'), ('independent', 'VERIFY.toml'), ('objective', 'OBJECTIVE.toml'))}
    for field_name, phase, row_name, key in (
        ('retained_physical_rows', 'physical', 'rows', logical),
        ('retained_independent_rows', 'independent', 'rows', logical),
        ('retained_objective_rows', 'objective', 'logical_rows', logical),
        ('retained_raw_witnesses', 'objective', 'raw_witnesses', lambda r: (r['case'], r['name'])),
        ('retained_reference_records', 'physical', 'reference_records', lambda r: (r['case'], r['time']))):
        actual = [dict(case=case, **row) for case in CASES for row in reports[case][phase][row_name]]
        require_same_rows(a[field_name], actual, key)
    require(len(a['retained_reference_records']) == 4, 'Four case/time reference records')
    require(set(a['primary']) == set(QUESTIONS), 'Both exact primary decision objects')
    censor = censored_report_scope(a, reports, case_metadata)
    phase_records = receipt_records(e, freeze)
    constructors, constructor_rows, orders = {}, [], {}
    for key in sorted(expected):
        case, rep, arm, budget = key
        path = f'constructors/{case}.rep{rep}.{arm}.B{budget}.json'
        require(path in a['artifact_hashes'], 'Constructor absent from ANALYSIS binding')
        record = e.json(path)
        require(identity(record) == key and record['source_freeze_sha256'] == FREEZE
                and record['protocol_sha256'] == a['protocol_sha256'] and record['score_kind'] == 'diagonal_reference_cubic_response',
                'Selected constructor identity')
        require(record['eligible'] is True and crows[key]['eligible'] is True, 'Ineligible constructor: refuse, never silently drop')
        for name in ('own_plan', 'own_receipt'):
            e.bind(record[name]['path'], record[name]['sha256'])
        if record.get('worker_output_path') is not None:
            e.bind(record['worker_output_path'], record['worker_output_sha256'])
        own = costs(crows[key]['actual_own']); attributed = costs(crows[key]['comparison_attributed'])
        prefix, continuation, total = (costs(record['costs'][name]) for name in ('prefix', 'continuation', 'total'))
        require(own == (prefix if arm == 'parent' else continuation), 'Actual own-cost join')
        require(attributed == (dict.fromkeys(COSTS, 0.) if arm == 'parent' else total), 'Attribution join')
        for score in ('score', 'cubic_response'):
            number(record[score])
        require(type(record['steps']) is int and record['steps'] > 0, 'Selected affordable steps')
        sha, order = order_hash(record['order'], sizes[case], orders)
        complete_cost = total if arm in OPTIMIZED else own
        result = record.get('result') or {}
        witness = field(result, 'raw_julia_candidate', 'witness')
        row = dict(case=case, repetition=rep, arm=arm, cpu_budget=budget,
            selected_J=record['score'], selected_cubic_response=record['cubic_response'], selected_steps=record['steps'],
            order_sha256=sha, order_zero_based=order, constructor_path=path, constructor_sha256=e.files[path]['sha256'],
            eligible=record['eligible'], status=result.get('status'), refusal=record.get('refusal'),
            case_scientific_gate_complete=case_metadata[case]['complete'], scientific_gate_complete=False,
            source_analysis_status=a['status'],
            result_refusal=result.get('refusal'), controller_fallback=record.get('controller_fallback'),
            execution_admission=record.get('execution_admission'), result_evaluations=result.get('evaluations'),
            result_accepted_moves=result.get('accepted_moves'), witness_start_scores_completed=field(witness, 'start_scores_completed'),
            witness_evaluations=field(witness, 'evaluations'), witness_proposals=field(witness, 'proposals'),
            witness_chain_proposals=field(witness, 'chain_proposals'), witness_chain_CPU_seconds=field(witness, 'chain_CPU_seconds'),
            work_stats=result.get('stats'), complete_result=result,
            **cost_columns(actual_own=own, comparison_attributed=attributed, complete_comparison=complete_cost,
                           prefix=prefix, continuation=continuation, total=total))
        constructors[key] = dict(record=record, row=row, order=order, order_sha256=sha)
        constructor_rows.append(row)
    physical_rows = []
    for key, row in sorted(physical.items()):
        ctor = constructors[key[:4]]
        sha, order = order_hash(row['order'], sizes[row['case']], orders, one_based=True)
        require(order == ctor['order'] and row['constructor_sha256'] == ctor['row']['constructor_sha256'], 'Physical/constructor join')
        require(row['eligible'] is True and row['status'] == independent[key]['status'] == 'complete'
                and objective[key]['status'] == 'verified', 'Complete physical/independent/objective eligibility')
        for name in ('projective_error', 'raw_normalized_overlap_deficit', 'infidelity', 'overlap_abs',
                     'steps', 'ideal_rotations', 'norm_drift', 'number_leakage', 'spin_leakage',
                     'reference_norm_drift', 'reference_number_leakage', 'reference_spin_leakage',
                     'actual_file', 'actual_sha256', 'reference_file', 'reference_sha256', 'circuit_sha256', 'numerical_floor'):
            require(name in row, 'Missing recorded physical field: '+name)
        number(row['projective_error'])
        e.bind(f"physical/{row['case']}/"+relative(row['actual_file']), row['actual_sha256'])
        e.bind(f"physical/{row['case']}/"+relative(row['reference_file']), row['reference_sha256'])
        require(e.files[f"physical/{row['case']}/"+row['actual_file']]['bytes'] == 16*4096
                and e.files[f"physical/{row['case']}/"+row['reference_file']]['bytes'] == 16*4096
                and independent[key]['independent_unfused_replay'] is True,
                'Full-register saved and independently replayed states')
        physical_rows.append(dict(row, order_zero_based=order, order_sha256=sha,
            selected_training_J=ctor['row']['selected_J'],
            selected_training_cubic_response=ctor['row']['selected_cubic_response'],
            selected_training_steps=ctor['row']['selected_steps'],
            independent_status=independent[key]['status'], objective_status=objective[key]['status'],
            independent_record=independent[key], objective_record=objective[key],
            constructor_status=ctor['row']['status'], constructor_refusal=ctor['row']['refusal'],
            selected_row_checks_verified=True, case_scientific_gate_complete=case_metadata[row['case']]['complete'],
            scientific_gate_complete=False, source_analysis_status=a['status'],
            **{k: v for k, v in ctor['row'].items() if any(k.startswith(prefix+'_') for prefix in
                ('actual_own', 'comparison_attributed', 'complete_comparison'))}))
    pair_expected = {(case, rep, question, budget, time, cap) for case in CASES for rep in (1, 2)
        for question in QUESTIONS for budget in (32, 64) for time in (.65, 1.35) for cap in (8, 24)}
    pairs = index(a['named_pairs'], lambda r: tuple(r[k] for k in
        ('case', 'repetition', 'question', 'cpu_budget', 'time', 'budget_multiplier')), pair_expected)
    require(len(pairs) == 64, '64 named pairs')
    pair_rows = []
    for key, row in sorted(pairs.items()):
        case, rep, question, budget, time, cap = key
        arm = question.removesuffix('_minus_full')
        left, right = physical[(case, rep, arm, budget, time, cap)], physical[(case, rep, 'full', budget, time, cap)]
        require(row['difference'] == left['projective_error']-right['projective_error']
                and row['candidate_error'] == left['projective_error'] and row['full_error'] == right['projective_error']
                and row['candidate_floor'] == left['numerical_floor'] and row['full_floor'] == right['numerical_floor']
                and row['eligible'] is case_metadata[case]['complete'], 'Copied paired raw-value and case-gate join')
        require(row['outcome'] == 'incomplete' if case == CENSORED_CASE else
                row['outcome'] in ('candidate_advantage', 'full_advantage', 'unresolved'),
                'Retain source pair outcomes; no censor-case promotion')
        primary = (budget, time, cap) == (64, 1.35, 8)
        pair_rows.append(dict(row, primary_condition=primary,
            primary_window_minus_full=primary and question == 'window_minus_full',
            primary_annealing_minus_full=primary and question == 'annealing_minus_full',
            case_scientific_gate_complete=case_metadata[case]['complete'], scientific_gate_complete=False,
            source_analysis_status=a['status']))
    raw_rows = []
    for row in a['retained_raw_witnesses']:
        key = identity(row); ctor = constructors[key]
        sha, order = order_hash(row['order'], sizes[row['case']], orders, one_based=True)
        checks = reports[row['case']]['objective']['checks']; check_index = row['check_index']
        require(type(check_index) is int and 1 <= check_index <= len(checks), 'Persisted raw objective check index')
        replay = checks[check_index-1]
        require(row['constructor_sha256'] == ctor['row']['constructor_sha256']
                and row['status'] in ('verified', 'native_refusal_censored') and replay['status'] == 'replayed'
                and replay['order'] == row['order'], 'Persisted raw claim binding')
        original_refusal = None
        if row['status'] == 'native_refusal_censored':
            known_native_censor(row)
            original = field(ctor['record'], 'result', 'stats', 'raw_julia_candidate')
            require(isinstance(original, dict) and original['order'] == order
                    and original.get('native_refusal') == NATIVE_REFUSAL and 'native_score' not in original
                    and order != ctor['order'], 'Unpromoted raw censor is the original persisted constructor refusal')
            original_refusal = original['native_refusal']
        raw_rows.append(dict(row, order_zero_based=order, order_sha256=sha,
            replayed_score=replay['replayed_score'], replayed_steps=replay['steps'],
            exact_replay_check=replay, constructor_status=ctor['row']['status'],
            original_constructor_native_refusal=original_refusal,
            case_scientific_gate_complete=case_metadata[row['case']]['complete'], scientific_gate_complete=False))
    cheap_rows = []
    for key, candidate in sorted(physical.items()):
        case, rep, arm, budget, time, cap = key
        if arm not in OPTIMIZED:
            continue
        for cheap in ('canonical', 'magnitude'):
            control = physical[(case, rep, cheap, 0, time, cap)]
            cheap_rows.append(dict(case=case, repetition=rep, candidate_arm=arm, candidate_cpu_budget=budget,
                cheap_arm=cheap, cheap_cpu_identity=0, time=time, budget_multiplier=cap,
                candidate_error=candidate['projective_error'], cheap_error=control['projective_error'],
                raw_candidate_minus_cheap=candidate['projective_error']-control['projective_error'],
                candidate_floor=candidate['numerical_floor'], cheap_floor=control['numerical_floor'],
                candidate_eligible=candidate['eligible'], cheap_eligible=control['eligible'],
                paired_eligible=candidate['eligible'] and control['eligible'],
                candidate_status=candidate['status'], cheap_status=control['status'],
                candidate_constructor_status=constructors[key[:4]]['row']['status'],
                cheap_constructor_status=constructors[(case, rep, cheap, 0)]['row']['status'],
                candidate_order_sha256=constructors[key[:4]]['order_sha256'],
                cheap_order_sha256=constructors[(case, rep, cheap, 0)]['order_sha256'],
                case_scientific_gate_complete=case_metadata[case]['complete'], scientific_gate_complete=False,
                eligibility_scope='Copied constructor/physical row eligibility; not primary scientific-gate eligibility',
                descriptive_only=True, primary_decision_role='none'))
    require(len(cheap_rows) == 192, '192 complete descriptive cheap comparisons')
    # Recheck anchor bytes before publishing derived metadata. No DATA write occurs.
    require(stream_file(data, 'EVIDENCE_INVENTORY.json', maximum=MAX_JSON)[0] == e.inventory_identity,
            'Seal changed during export')
    require(stream_file(analysis_path.parent, analysis_path.name, maximum=MAX_JSON)[0] == analysis_identity,
            'ANALYSIS changed during export')
    script_identity, script = stream_file(Path(__file__).resolve().parent, Path(__file__).name, maximum=1024**2, capture=True)
    output.mkdir(parents=True, exist_ok=False)
    emitted = {}
    def put(name, content):
        with (output/name).open('xb') as stream:
            stream.write(content)
        emitted[name] = dict(bytes=len(content), sha256=digest(content))
    def json_output(name, value):
        put(name, (json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False)+'\n').encode('ascii'))
    def csv_output(name, rows, empty_fields=()):
        import io
        columns = list(dict.fromkeys(k for row in rows for k in row)) or list(empty_fields)
        stream = io.StringIO(newline=''); writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: canonical(value).decode('ascii') if value is None or isinstance(value, (dict, list, bool))
                             else value for key in columns for value in (row.get(key),)})
        put(name, stream.getvalue().encode('utf-8'))
        return dict(rows=len(rows), fields=columns)
    counts = {
        'constructors.csv': csv_output('constructors.csv', constructor_rows),
        'physical.csv': csv_output('physical.csv', physical_rows),
        'named_pairs.csv': csv_output('named_pairs.csv', pair_rows),
        'persisted_raw_claims.csv': csv_output('persisted_raw_claims.csv', raw_rows,
            ('case', 'name', 'order_id', 'arm', 'repetition', 'cpu_budget', 'status', 'order_sha256')),
        'cheap_context.csv': csv_output('cheap_context.csv', cheap_rows)}
    json_output('PRIMARY_QUESTIONS.json', a['primary'])
    json_output('SOURCE_ANALYSIS_STATUS.json', dict(source_complete=a['complete'], source_status=a['status'],
        source_issues=a['issues'], source_cases_metadata=a['cases_metadata'],
        source_scientific_phase_costs=a.get('scientific_phase_costs'), source_primary=a['primary'],
        scientific_gate_complete=False, export_delivery_complete=True))
    json_output('CENSOR_RECORD.json', censor)
    json_output('SCIENTIFIC_PHASE_RECEIPTS.json', dict(records=phase_records,
        source_receipt_gate_issue=a['issues'], scientific_gate_complete=False,
        scope='All six original job costs/statuses retained; HF objective exit2 is not changed to scientific success'))
    put('SCIENTIFIC_SCOPE.md', (
        '# Complete table inventory; incomplete scientific gate\n\n'
        'All 36 eligible constructors and 144 physical, independent and selected-objective rows are retained. '
        'Every selected training score agrees with its archived replay/native comparison. '
        'Every persisted raw proposed score passes its archived proposed-score comparison.\n\n'
        'Exactly one unpromoted HF rep2.window.B64.raw1 proposal has the original native SearchDeadline '
        '(whole-call window allowance exhausted). Its native score remains absent and native_refusal_censored '
        'remains unchanged. No constructor was retried, raw proposal promoted or native score backfilled.\n\n'
        'Under the strict frozen rules ANALYSIS.complete is false, its status remains incomplete, '
        'and both primary decisions remain incomplete. Complete CSV inventory is delivery completeness '
        'only; no primary winner or overall scientific-goal success is claimed. Source issues and exit2 '
        'cost/receipt evidence are retained without repair.\n').encode('utf-8'))
    json_output('ORDERS.json', orders)
    put('export_consumer.py', script)
    manifest = dict(schema='molecular-main-censored-table-export-v1', export_delivery_complete=True,
        scientific_gate_complete=False, source_analysis_complete=a['complete'], source_analysis_status=a['status'],
        source_analysis_issues=a['issues'], source_cases_metadata=a['cases_metadata'],
        exact_allowed_censor=dict(case=CENSORED_CASE, name=CENSORED_NAME, original_native_refusal=NATIVE_REFUSAL),
        scientific_role='Full attempted-matrix metadata delivery with exactly one original native-censored unpromoted raw claim; scientific gate and both primaries remain incomplete; no numerical replay, repair or new decision',
        producer_commit=PRODUCER, source_freeze_sha256=FREEZE, protocol_sha256=a['protocol_sha256'],
        data_root=str(data), analysis_path=str(analysis_path), analysis_identity=analysis_identity,
        analysis_sealed_members=analysis_sealed_members,
        inventory_identity=e.inventory_identity, sealed_artifacts=e.files,
        original_producer_sources=freeze['sources'], analysis_consumer_sources=a['current_consumer_sources'],
        consumer_identity=script_identity, python=sys.version, executable=sys.executable,
        order_serialization='SHA256 of ASCII json.dumps(zero_based_original_label_integer_list, ensure_ascii=True, separators=(comma,colon), allow_nan=False); no whitespace/newline. Hamiltonian/case identity is separate.',
        primary_questions='Exact source incomplete decision objects; no recalculated medians/outcomes or equivalent winner claim',
        work_scope='Counters describe search work; only selected orders and persisted raw claims are objective-replayed, not a transient proposal trajectory',
        ratio_policy='No ratios exported. Established comparisons retain raw errors and flags; ratios are suppressed when either error <=1e-12, not produced by clipping.',
        cost_policy='actual_own and comparison_attributed copy ANALYSIS; parent aggregate attribution is zero by design, not free construction; complete_comparison uses own cost for parent/cheap and total prefix+continuation for optimized; repeated physical rows never multiply spending',
        null_policy='Literal null means field/counter not persisted; never imputed as zero',
        csv=counts, outputs=emitted)
    # Manifest excludes its own digest to avoid a self-hash cycle; write last.
    with (output/'MANIFEST.json').open('xb') as stream:
        stream.write((json.dumps(manifest, indent=2, ensure_ascii=True, allow_nan=False)+'\n').encode('ascii'))
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('analysis_json', type=Path)
    parser.add_argument('fresh_output', type=Path)
    args = parser.parse_args(argv)
    data, analysis, output = (p.absolute() for p in (args.data, args.analysis_json, args.fresh_output))
    require(data == data.resolve(strict=True) and analysis == analysis.resolve(strict=True)
            and output == output.resolve(), 'Canonical absolute arguments without symlinks')
    counts = export(data, analysis, output)
    print(json.dumps(dict(export_delivery_complete=True, scientific_gate_complete=False,
        source_analysis_status='incomplete', csv_rows={name: value['rows'] for name, value in counts.items()}),
        sort_keys=True))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError, tomllib.TOMLDecodeError) as error:
        print('REFUSED: '+str(error), file=sys.stderr)
        sys.exit(2)
