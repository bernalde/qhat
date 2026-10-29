"""Authenticate a frozen continuation run and retain every paired outcome.

This consumer uses archived source hashes, never the current checkout's freeze
validator. Numerical propagation and objective replay remain separate evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import posixpath
import re
import statistics
import tomllib

COST_FIELDS = ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds')
IDENTITY = ('name', 'order_id', 'arm', 'repetition', 'cpu_budget', 'partition',
            'formula', 'time', 'budget_multiplier', 'order', 'constructor_sha256',
            'parent_sha256', 'claimed_score', 'score_kind')
KIND = 'diagonal_reference_cubic_response'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    try:
        valid = type(value) in (int, float) and math.isfinite(value) and value >= 0
    except OverflowError:
        valid = False
    require(valid,
            'Finite nonnegative number required')
    return value


def unique(records, key):
    result = {}
    for record in records:
        identity = key(record)
        require(identity not in result, 'Duplicate identity: '+repr(identity))
        result[identity] = record
    return result


def logical(row):
    return (row['arm'], row['repetition'], row['cpu_budget'],
            row['time'], row['budget_multiplier'])


def order_key(row):
    return row['case'], row['repetition'], row['arm'], row['total_CPU_budget']


def cost_check(record, parent, config):
    """Check strict ceilings separately from floating addition reconciliation."""
    costs = record['costs']
    for boundary in ('prefix', 'continuation', 'total'):
        for field in COST_FIELDS:
            if field != 'raw_CPU_seconds' or costs[boundary][field] is not None:
                number(costs[boundary][field])
        require(costs[boundary]['raw_CPU_seconds'] is None or
                costs[boundary]['charged_CPU_seconds'] >= costs[boundary]['raw_CPU_seconds'],
                'Conservative CPU below raw CPU')
    for field in COST_FIELDS:
        require(costs['prefix'][field] == parent['costs']['total'][field],
                'Shared parent prefix differs')
        operands = [costs[part][field] for part in ('prefix', 'continuation')]
        expected = None if any(value is None for value in operands) else sum(operands)
        require((costs['total'][field] is None and expected is None) or
                (costs['total'][field] is not None and expected is not None and
                 math.isclose(costs['total'][field], expected, rel_tol=0., abs_tol=1e-9)),
                'Prefix/continuation total differs')
    cap = config['parent']['CPU_seconds'] if record['arm'] == 'parent' else record['total_CPU_budget']
    if record.get('eligible', True):
        require(costs['total']['raw_CPU_seconds'] is not None and
                costs['total']['charged_CPU_seconds'] <= cap and
                costs['total']['raw_CPU_seconds'] <= cap, 'Strict whole-construction CPU ceiling')
    if record['arm'] == 'parent':
        require(all(costs['continuation'][field] == 0 for field in COST_FIELDS), 'Parent continuation is not zero')
        if record.get('eligible', True):
            require(costs['total']['wall_seconds'] <= config['parent']['wall_seconds'], 'Parent wall ceiling')
    elif record.get('eligible', True):
        require(costs['total']['wall_seconds'] <= config['accounting']['total_wall_multiplier']*cap,
                'Strict whole-construction wall ceiling')


def json_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: '+key)
        result[key] = value
    return result


def json_float(text):
    value = float(text)
    require(math.isfinite(value), 'Nonfinite JSON number')
    return value


def parse_json(raw):
    def reject(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    # Bound parser nesting independently of the byte ceiling.
    depth, quoted, escaped = 0, False, False
    for byte in raw:
        if quoted:
            if escaped:
                escaped = False
            elif byte == 92:
                escaped = True
            elif byte == 34:
                quoted = False
        elif byte == 34:
            quoted = True
        elif byte in (91, 123):
            depth += 1
            require(depth <= 64, 'JSON nesting cap')
        elif byte in (93, 125):
            depth -= 1
    return json.loads(raw, object_pairs_hook=json_object, parse_constant=reject, parse_float=json_float)


def command_flags(command, cwd, script, value_flags, switches=()):
    """Parse one fixed Python producer grammar; never execute copied text."""
    require(isinstance(command, list) and 2 <= len(command) <= 32 and
            all(isinstance(value, str) and 0 < len(value) <= 4096 and '\x00' not in value for value in command), 'Bounded command argv')
    require(posixpath.isabs(cwd) and posixpath.normpath(cwd) == cwd and '\\' not in cwd, 'Absolute producer cwd')
    expected = posixpath.normpath(cwd+'/numerical_audit/cost_matched_continuation/'+script)
    actual = posixpath.normpath(posixpath.join(cwd, command[1]))
    require(actual == expected, 'Wrong producer script')
    result, at = {}, 2
    while at < len(command):
        flag = command[at]
        require(flag not in result and flag in set(value_flags) | set(switches), 'Duplicate/unknown producer flag')
        at += 1
        if flag in switches:
            result[flag] = True
        else:
            require(at < len(command) and not command[at].startswith('--'), 'Missing producer flag value')
            result[flag] = command[at]
            at += 1
    return result


def bind_constructor(plan, receipt, final, owner, parent, worker, config):
    require(plan['phase'] == 'constructors', 'Owning job is not a constructor')
    flags = command_flags(plan['command'], plan['cwd'], 'construct.py',
        ('--root', '--case', '--repetition', '--arm', '--budget', '--raw', '--cpu-cap', '--wall-cap', '--julia-executable'), ('--worker',))
    require(flags.get('--worker') is True, 'Constructor worker switch required')
    root = flags['--root']
    require(posixpath.isabs(root) and posixpath.normpath(root) == root and '\\' not in root, 'Absolute original data root')
    relative = owner['worker_output_path']
    require(isinstance(relative, str) and not posixpath.isabs(relative) and posixpath.normpath(relative) == relative and
            '\\' not in relative and '..' not in relative.split('/'), 'Canonical relative worker destination')
    require((flags['--case'], int(flags['--repetition']), flags['--arm'], int(flags.get('--budget', '0'))) ==
            order_key(owner), 'PLAN constructor identity differs')
    require(posixpath.normpath(flags['--raw']) == posixpath.join(root, owner['worker_output_path']),
            'PLAN worker destination differs')
    cpu, wall = number(float(flags['--cpu-cap'])), number(int(flags['--wall-cap']))
    require(cpu > 0, 'Positive constructor CPU allowance')
    require(cpu == plan['attempt_CPU_seconds'] and wall == plan['attempt_wall_seconds'], 'PLAN/CLI allowance differs')
    expected_cpu = config['parent']['CPU_seconds'] if owner['arm'] == 'parent' else (
        owner['total_CPU_budget']-parent['costs']['total']['charged_CPU_seconds'])
    expected_wall = config['parent']['wall_seconds'] if owner['arm'] == 'parent' else math.floor(
        5*owner['total_CPU_budget']-parent['costs']['total']['wall_seconds'])
    require(0 <= expected_cpu-cpu <= 1e-12 and wall == expected_wall, 'Frozen remaining allowance differs')
    require(receipt['command'] == plan['command'], 'Receipt/PLAN command differs')
    require(receipt['budget']['cpu_seconds'] == cpu and receipt['budget']['wall_seconds'] == wall,
            'Receipt guarded allowance differs')
    if worker is not None:
        require(worker['guarded_CPU_cap'] == cpu and worker['guarded_wall_cap'] == wall, 'Worker guarded allowance differs')
        require(final['charged_CPU_seconds'] >= number(worker['process_CPU_before_output']),
                'Conservative CPU below worker self clock')
    if 'cost_binding' in owner:
        require(owner['cost_binding'] == dict(original_root=root, original_worker_path=flags['--raw'],
            CPU_allowance=cpu, wall_allowance=wall, verified=True), 'Producer/consumer cost binding differs')
    return root


def derived_execution(arm, result, parent_order=None, refusal=None):
    """Independently classify raw execution metadata; do not import producers."""
    stats = result.get('stats', {}) if isinstance(result, dict) else {}
    stats = stats if isinstance(stats, dict) else {}
    ws, co, cr = [stats.get(key) for key in ('window_status', 'child_outcome', 'child_refusal')]
    def answer(yes, category, reason):
        return dict(version=1, admitted=yes, category=category, reason=reason,
            **{key: value if isinstance(value, str) else None for key, value in
               zip(('window_status', 'child_outcome', 'child_refusal'), (ws, co, cr))})
    def deny(category, reason):
        return answer(False, category, reason)
    def integer(value, low=0):
        return type(value) is int and value >= low
    if arm not in ('parent', 'full', 'window'):
        return deny('unknown_or_inconsistent', 'unknown_arm')
    if refusal is not None:
        return deny('worker_refusal', 'outer_worker_refusal')
    status = result.get('status') if isinstance(result, dict) else None
    if status == 'numerical_limit':
        return deny('numerical_refusal', 'numerical_limit')
    try:
        order, score = result['order'], result['score']
        require(isinstance(order, (list, tuple)) and order and all(type(x) is int for x in order) and
            sorted(order) == list(range(len(order))), 'Order')
        number(score['value']); number(score['cubic_response'])
        require(integer(score['steps'], 1), 'Steps')
    except (ValueError, TypeError, KeyError):
        return deny('unscored', 'missing_or_invalid_scored_incumbent')
    if arm == 'parent':
        if status in ('complete', 'cpu_limit'):
            return answer(True, 'completed' if status == 'complete' else 'budget_stop', status)
        if status == 'evaluation_limit' and integer(result.get('evaluations')) and result['evaluations'] in (2999, 3000):
            return answer(True, 'declared_parent_limit', 'declared_parent_evaluation_limit')
        return deny('unknown_or_inconsistent', 'unadmitted_parent_status_or_evaluation_count')
    if arm == 'full':
        if status not in ('complete', 'cpu_limit'):
            return deny('algorithm_control_stop' if status in ('evaluation_limit', 'support_limit') else
                'unknown_or_inconsistent', 'unadmitted_full_status')
        if not integer(stats.get('seed_evaluations')) or stats['seed_evaluations'] != 4:
            return deny('insufficient_starts', 'full_control_requires_four_scored_starts')
        return answer(True, 'completed' if status == 'complete' else 'budget_stop', status)
    if (not isinstance(parent_order, (list, tuple)) or not parent_order or len(parent_order) != len(order) or
            any(type(x) is not int for x in parent_order) or sorted(parent_order) != list(range(len(order)))):
        return deny('unknown_or_inconsistent', 'missing_or_invalid_frozen_parent_order')
    changed = list(order) != list(parent_order)
    child_stops = {'child_'+tag for tag in ('cpu_limit', 'wall_limit', 'evaluation_limit', 'action_limit', 'atom_limit')}
    if not isinstance(ws, str) or any(x is not None and not isinstance(x, str) for x in (co, cr)):
        return deny('unknown_or_inconsistent', 'missing_or_malformed_window_status')
    if status not in ('complete', 'cpu_limit', 'evaluation_limit'):
        return deny('unknown_or_inconsistent', 'unadmitted_native_window_status')
    if co is not None and co not in child_stops | {'complete'}:
        category = {'child_numerical_limit': 'numerical_refusal', 'child_admission_or_support_limit':
            'ambiguous_support_refusal', 'owned_child_cpu_unavailable': 'runtime_unavailable'}.get(co, 'backend_failure')
        return deny(category, co)
    if cr is not None and ('child_'+cr not in child_stops or co != 'child_'+cr):
        return deny({'numerical_limit': 'numerical_refusal', 'admission_or_support_limit':
            'ambiguous_support_refusal'}.get(cr, 'backend_failure'), 'refused_or_inconsistent_child:'+cr)
    failures = dict.fromkeys(('julia_unavailable', 'unsupported_toml_reader', 'unsupported_owned_cpu_accounting_platform',
        'unsupported_owned_cpu_reaper', 'owned_child_cpu_unavailable'), 'runtime_unavailable')
    failures.update(dict.fromkeys(('invalid_worker_reply', 'missing_worker_reply', 'native_julia_score_mismatch',
        'child_launch_or_transport_failure', 'child_failure', 'child_worker_failure'), 'backend_failure'))
    failures.update(dict.fromkeys(('unsupported_register_or_unit_cap', 'unsupported_word_cap', 'unsupported_coherent_reference',
        'unsupported_storage_cap', 'window_unavailable', 'unsupported_empty_or_zero_word_unit',
        'unsupported_complex_or_mixed_flip_unit'), 'domain_refusal'))
    failures.update(numerical_limit='numerical_refusal', child_numerical_limit='numerical_refusal',
        child_admission_or_support_limit='ambiguous_support_refusal')
    if ws in failures:
        return deny(failures[ws], ws)
    if co is not None and (stats.get('backend_runtime_assets_verified') is not True or
            not integer(stats.get('child_returncode'), -2**31) or co == 'complete' and stats['child_returncode'] != 0):
        return deny('unknown_or_inconsistent', 'unverified_assets_or_inconsistent_child_returncode')
    if stats.get('parent_method') != 'magnitude' or not integer(stats.get('parent_evaluations')) or stats['parent_evaluations'] != 1:
        return deny('unknown_or_inconsistent', 'supplied_parent_not_scored_once')
    capacity = ('reachable vector exceeds max_entries', 'cached reference columns exceed max_entries',
        'reachable unit column exceeds max_entries', 'combined cubic columns exceed max_entries',
        'updated cubic columns exceed max_entries', 'window rational exceeds 4096-bit admission cap', 'window request exceeds 8 MiB')
    if ws == 'support_limit':
        reason = stats.get('refusal_reason')
        if reason not in capacity:
            asset = isinstance(reason, str) and (reason.startswith(('backend ', 'malformed backend ', 'installed Julia ')) or
                reason == 'window backend requires Python 3.11+ tomllib')
            return deny('runtime_asset_failure' if asset else 'ambiguous_support_refusal',
                reason if isinstance(reason, str) else 'missing_or_malformed_support_reason')
        category = 'capacity_stop'
    elif ws in ('whole_call_limit', 'evaluation_limit', 'insufficient_evaluation_budget', 'verification_reserve_unavailable') or ws in child_stops:
        if ws in child_stops and co != ws:
            return deny('unknown_or_inconsistent', 'child_stop_outcome_differs')
        category, reason = 'budget_stop', ws
    else:
        category = None
    if category:
        return deny('unverified_checkpoint', 'stop_selected_unverified_changed_order') if changed else answer(True, category, reason)
    counts = stats.get('child_counters', {})
    completed, allocated, updates = [stats.get(key) for key in ('child_completed_evaluations', 'child_allocated_evaluations', 'child_checked_updates')]
    state = stats.get('child_checkpoint_status')
    parsed = (co is not None and stats.get('backend_runtime_assets_verified') is True and
        integer(stats.get('native_parent_revalidation_evaluations')) and stats['native_parent_revalidation_evaluations'] == 1 and
        stats.get('child_evaluation_scope') == 'checkpoint_prefix' and integer(completed, 1) and integer(allocated, 1) and
        completed <= allocated and isinstance(counts, dict) and all(integer(counts.get(key)) for key in
            ('evaluation', 'full_score', 'flat_candidate', 'exact_witness', 'action', 'atom', 'insertion_move')) and
        counts['evaluation'] == completed and counts['full_score'] >= 1 and
        sum(counts[k] for k in ('full_score', 'flat_candidate', 'exact_witness')) <= completed and integer(updates) and
        counts['full_score'] >= updates+1 and state in ('parent_checked', 'checked_update', 'complete') and
        (state != 'parent_checked' or updates == 0) and (state != 'checked_update' or updates >= 1))
    if not parsed:
        return deny('unverified_checkpoint', 'missing_or_invalid_parsed_checkpoint_metadata')
    if co == 'complete' and state != 'complete':
        return deny('unknown_or_inconsistent', 'completed_child_has_incomplete_checkpoint')
    verified = all(type(stats.get(key)) is int and stats[key] == 1 for key in ('native_verifications', 'native_verification_attempts'))
    if ws == 'complete' and co == 'complete' and not changed:
        return answer(True, 'completed', 'checked_complete_without_required_improvement')
    if ws == 'native_no_resolved_improvement' and not changed and verified:
        return answer(True, 'checked_no_improvement', ws)
    if changed and verified and updates >= 1 and ((ws == 'verified_improvement' and co == 'complete') or
            (co in child_stops and ws == 'verified_checkpoint_after_'+co)):
        return answer(True, 'verified_checkpoint', ws)
    return deny('unknown_or_inconsistent', 'unknown_or_inconsistent_window_outcome')


def terminal_accounting(receipt, final, log, reserved):
    """A polling sample alone is never a terminal CPU measurement."""
    measured, footers, partial = [], [], []
    process = receipt.get('process_measurements', {})
    try:
        process_cpu = number(number(process['user_cpu_seconds'])+number(process['system_cpu_seconds']))
        partial.append(process_cpu)
    except (ValueError, TypeError, KeyError):
        process_cpu = None
    process_valid = receipt.get('final_accounting', {}).get('process_CPU_timing_available') is True and process_cpu is not None
    if process_valid:
        measured.append(process_cpu)
    factors = dict(h=3600., min=60., s=1., ms=.001, us=1e-6, **{'µs': 1e-6, 'ns': 1e-9})
    for line in re.findall(r'^CPU time consumed: (.+)$', log, re.M):
        tokens = list(re.finditer(r'([0-9]+(?:\.[0-9]+)?)(min|ms|us|µs|ns|h|s)', line))
        residue = re.sub(r'([0-9]+(?:\.[0-9]+)?)(min|ms|us|µs|ns|h|s)', '', line)
        if tokens and not residue.strip():
            value = sum(float(match[1])*factors[match[2]] for match in tokens)
            if math.isfinite(value):
                footers.append(value)
    measured.extend(footers)
    snapshot = final.get('final_service_snapshot', {})
    service = snapshot.get('CPUUsageNSec')
    service_valid = final.get('final_service_query_exit') == 0 and isinstance(service, str) and re.fullmatch(r'[0-9]+', service)
    if service_valid:
        measured.append(int(service)/1e9)
    guard = receipt.get('final_accounting', {}).get('CPU_seconds')
    if guard is not None:
        number(guard)
        partial.append(guard)
    for sample in receipt.get('samples', []):
        partial.extend(int(value)/1e9 for value in re.findall(r'(?:^|\n)CPUUsageNSec=([0-9]+)(?:\n|$)', sample.get('service', '')))
    available = bool(measured)
    raw = max(measured+([guard] if guard is not None else [])) if available else None
    charged = max(partial+measured+[value+.001 for value in footers]+
        ([process_cpu+.02] if process_valid and not footers and not service_valid else [])+
        ([] if available else [reserved]))
    return dict(raw_terminal_CPU_seconds=raw, terminal_CPU_measurement_available=available,
                guard_reported_CPU_seconds=guard, partial_CPU_seconds=max(partial, default=None),
                minimum_charged_CPU_seconds=charged)


def cgroup_peak(receipt, final, log):
    values = []
    for text in [sample.get('service', '') for sample in receipt.get('samples', [])]:
        values.extend(float(value) for value in re.findall(r'(?:^|\n)MemoryPeak=([0-9]+)(?:\n|$)', text))
    snapshot = final.get('final_service_snapshot', {}).get('MemoryPeak')
    if final.get('final_service_query_exit') == 0 and isinstance(snapshot, str) and snapshot.isascii() and snapshot.isdigit():
        values.append(float(snapshot))
    for value, suffix in re.findall(r'^Memory peak: ([0-9]+(?:\.[0-9]+)?)([KMGT]?)(?:B)?$', log, re.M):
        unit = 1024.**('BKMGT'.index(suffix or 'B'))
        # Printed footer values are rounded; add one unit of their last digit.
        precision = len(value.split('.')[1]) if '.' in value else 0
        values.append((float(value)+10.**(-precision))*unit)
    return max(values) if values else None


def pair_grid(cases, repetitions, config, physical, constructors, verified):
    """Emit the entire grid, including unavailable or ineligible comparisons."""
    ambiguity, floor = config['physical']['ambiguity'], config['physical']['ratio_floor']
    pairs = []
    for case in cases:
        for rep in range(1, repetitions+1):
            for budget in config['total_CPU_budgets']:
                for duration in config['physical']['times']:
                    for cap in config['physical']['rotation_cap_multipliers']:
                        rows = [physical.get((case, arm, rep, budget, duration, cap)) for arm in ('full', 'window')]
                        errors = [None if row is None else number(row['projective_error']) for row in rows]
                        available = all(value is not None for value in errors)
                        delta = errors[1]-errors[0] if available else None
                        eligible = available and all(constructors.get((case, rep, arm, budget), {}).get('analysis_eligible')
                            and (case, arm, rep, budget, duration, cap) in verified for arm in ('full', 'window'))
                        pairs.append(dict(case=case, repetition=rep, total_CPU_budget=budget,
                            time=duration, rotation_cap_multiplier=cap, full_error=errors[0], window_error=errors[1],
                            difference=delta, outcome='unavailable' if delta is None else
                            'win' if delta < -ambiguity else 'loss' if delta > ambiguity else 'tie',
                            full_at_or_below_floor=errors[0] is not None and errors[0] <= floor,
                            window_at_or_below_floor=errors[1] is not None and errors[1] <= floor,
                            window_to_full_error_ratio=errors[1]/errors[0] if available and min(errors) > floor else None,
                            eligible=bool(eligible)))
    return pairs


def primary_decision(pairs, cases, repetitions, config, all_gates):
    primary = config['physical']['primary']
    selected = [row for row in pairs if (row['total_CPU_budget'], row['time'], row['rotation_cap_multiplier']) ==
                (primary['total_CPU_budget'], primary['time'], primary['rotation_cap_multiplier'])]
    medians = {case: statistics.median([row['difference'] for row in selected if row['case'] == case])
               if sum(row['case'] == case and row['difference'] is not None for row in selected) == repetitions
               else None for case in cases}
    ready = all_gates and len(selected) == len(cases)*repetitions and all(row['eligible'] for row in selected)
    ambiguity = config['physical']['ambiguity']
    decision = 'study_incomplete' if not ready else 'mixed_or_unresolved'
    if ready and all(value < -ambiguity for value in medians.values()) and not any(row['outcome'] == 'loss' for row in selected):
        decision = 'window_advantage'
    elif ready and all(value > ambiguity for value in medians.values()) and not any(row['outcome'] == 'win' for row in selected):
        decision = 'full_advantage'
    return dict(decision=decision, per_case_median_paired_difference=medians, rows=selected,
                all_primary_gates_pass=bool(ready), confidence_interval=None,
                repetition_scope=config['physical']['repetitions_scope'])


def resource_totals(constructors):
    parents = [row for row in constructors.values() if row['arm'] == 'parent']
    children = [row for row in constructors.values() if row['arm'] != 'parent']
    result = {}
    for field in COST_FIELDS:
        parent_values = [row['costs']['total'].get(field) for row in parents]
        continuation_values = [row['costs']['continuation'].get(field) for row in children]
        attributed_values = [row['costs'][boundary].get(field) for row in children for boundary in ('prefix', 'continuation')]
        summary = {}
        for name, values in dict(actual_shared_constructor_compute=parent_values+continuation_values,
                attributed_continuation_comparisons=attributed_values, actual_parent_compute=parent_values,
                actual_continuation_compute=continuation_values).items():
            known = math.fsum(number(value) for value in values if value is not None)
            unknown = sum(value is None for value in values)
            summary[name] = None if unknown else known
            summary['known_'+name] = known
            summary[name+'_unknown_components'] = unknown
        result[field] = summary
    return result


class Audit:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.inputs, self.issues = {}, []
        self.reconciled, self.original_roots = {}, set()
        self.preflight_original_root = None

    def attempt(self, scope, function, default=None):
        try:
            return function()
        except (ValueError, KeyError, TypeError, IndexError, OSError, OverflowError, RecursionError) as error:
            self.issues.append(dict(scope=scope, type=type(error).__name__, message=str(error)))
            return default

    def path(self, name):
        path = (self.root/name).resolve()
        require(path.is_relative_to(self.root), 'Path outside retained run')
        return path

    def load(self, name, digest=None, binary=False, text=False):
        path = self.path(name)
        require(path.is_file() and path.stat().st_size <= 32*1024**2, 'Missing/over-cap input: '+str(name))
        with path.open('rb') as stream:
            raw = stream.read(32*1024**2+1)
        require(len(raw) <= 32*1024**2, 'Input grew beyond byte cap')
        actual = hashlib.sha256(raw).hexdigest()
        self.inputs[str(path.relative_to(self.root))] = actual
        require(digest is None or digest == actual, 'Input digest differs: '+str(name))
        if binary:
            require(len(raw) == 16*16384, 'Full-register state byte count')
            return raw
        if text:
            return raw.decode('utf-8', errors='replace')
        return tomllib.loads(raw.decode()) if path.suffix == '.toml' else parse_json(raw)

    def source_hashes(self, sources, frozen, base=Path()):
        require(bool(sources), 'Missing producer/dependency source map')
        for name, digest in sources.items():
            canonical = str(Path('source_snapshot')/base/name)
            path = self.path(canonical)
            relative = str(path.relative_to(self.root/'source_snapshot'))
            require(frozen['sources'].get(relative) == digest, 'Producer source absent/different in source freeze: '+name)
            self.load(canonical, digest) if path.suffix in ('.json', '.toml') else self.hash_source(canonical, digest)

    def hash_source(self, name, digest):
        path = self.path(name)
        require(path.is_file() and path.stat().st_size <= 32*1024**2, 'Missing/over-cap archived source')
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        self.inputs[str(path.relative_to(self.root))] = actual
        require(actual == digest, 'Archived source digest differs: '+name)

    def worker(self, constructor):
        """Join finalized metadata to the retained worker's original witness."""
        path, digest = constructor['worker_output_path'], constructor['worker_output_sha256']
        fallback = constructor['controller_fallback']
        require(type(fallback) is bool, 'Explicit controller fallback flag')
        if digest is None:
            require(fallback and not self.path(path).exists() and constructor['result'] is None,
                    'Missing worker digest is inconsistent with retained evidence')
            return dict(verified=False, reason='missing_worker_output', path=path, worker_data=None)
        require(isinstance(digest, str), 'Worker output digest spelling')
        worker = self.load(path, digest)
        require(worker['schema'] == 'cost-matched-worker-v1', 'Worker schema')
        fields = ('case', 'repetition', 'arm', 'total_CPU_budget', 'model_json_sha256',
                  'model_toml_sha256', 'protocol_sha256', 'source_freeze_sha256', 'parent_sha256')
        require(all(worker[key] == constructor[key] for key in fields), 'Worker/canonical identity differs')
        result = worker['result']
        require(result == constructor['result'], 'Canonical native result differs from worker')
        if fallback:
            require(result is None or isinstance(result, dict) and result.get('score') is None,
                    'Controller fallback contradicts a scored worker result')
            return dict(verified=False, reason='worker_refusal' if worker.get('refusal') else
                        'missing_worker_result' if result is None else 'unscored_worker_result',
                        path=path, sha256=digest, refusal=worker.get('refusal'), worker_data=worker)
        require(isinstance(result, dict) and isinstance(result.get('score'), dict), 'Unscored worker without fallback')
        require(result['order'] == constructor['order'] and result['score'] == constructor['score'],
                'Worker/canonical selected order or score differs')
        return dict(verified=True, reason=None, path=path, sha256=digest, worker_data=worker)

    def constructor_graph(self, row, parent):
        worker = row['worker_verification'].get('worker_data')
        require(type(row['eligible']) is bool and type(row['controller_fallback']) is bool and
                type(row['guard_exit_code']) is int, 'Explicit constructor outcome flags')
        if row['controller_fallback']:
            require(row['eligible'] is False, 'Controller fallback cannot be eligible')
            if row['arm'] != 'parent':
                require(row['order'] == parent['order'] and row['score'] == parent['score'], 'Fallback differs from frozen parent')
        # Finalization adds contextual refusals after the native method returns.
        # Derive those from the actual graph, rather than trusting its saved flag.
        refusal = dict(type='missing_worker') if worker is None else worker.get('refusal')
        admission = derived_execution(row['arm'], None if worker is None else worker.get('result'),
            parent['order'] if row['arm'] == 'window' else None, refusal)
        if row['arm'] != 'parent' and parent['eligible'] is False:
            admission.update(admitted=False, category='ineligible_parent', reason='frozen_parent_not_admitted')
        if row['arm'] == 'parent' and not row['controller_fallback'] and row.get('restarts') is None:
            admission.update(admitted=False, category='missing_parent_preparation', reason='missing_measured_restart_plan')
        row['derived_execution_admission'] = admission
        if worker is not None:
            if row['arm'] == 'parent':
                require(worker.get('restarts') == row.get('restarts'), 'Parent restart witness differs')
                if not row['controller_fallback'] and row.get('restarts') is not None:
                    restarts = row['restarts']
                    expected_seed = int.from_bytes(hashlib.sha256((row['case']+':'+str(row['repetition'])+
                        ':cost-matched-v1').encode()).digest()[:8], 'big')
                    orders = restarts['orders']
                    require(restarts['seed'] == expected_seed and restarts['interior_perturbation_counts'] == [0, 1, 2, 4] and
                        len(orders) == 4 and len({tuple(value) for value in orders}) == 4 and orders[0] == row['order'] and
                        all(all(type(label) is int for label in value) and sorted(value) == list(range(len(row['order']))) and value[0] == row['order'][0] and
                            value[-1] == row['order'][-1] for value in orders), 'Frozen restart structure differs')
                require(worker.get('starts') is None, 'Native parent initialization changed')
            else:
                starts = worker.get('starts')
                if starts is None:
                    require(row['eligible'] is False and not admission['admitted'] and
                        (worker.get('result') is None or worker.get('refusal') is not None),
                        'Absent continuation starts without a failed/unstarted worker')
                else:
                    if row['arm'] == 'full':
                        require(isinstance(parent.get('restarts'), dict), 'Full continuation lacks a frozen restart graph')
                    wanted = parent['restarts']['orders'] if row['arm'] == 'full' else [parent['order']]
                    require(starts == wanted, 'Continuation starts differ from frozen parent graph')
        if 'execution_admission' in row:
            require(row['execution_admission'] == admission, 'Producer/consumer execution admission differs')
        if row['eligible']:
            require(not row['controller_fallback'] and row['guard_exit_code'] == 0 and
                row['worker_verification']['verified'] is True and admission['admitted'], 'Forged eligible execution')
        return admission

    def ledger(self, constructors, config=None, mode=None, cases=None):
        """Reconcile costs; scientific coverage additionally needs frozen context."""
        ledger = self.load('campaign/LEDGER_CUTOFF.json')
        require(ledger['schema'] == 'cost-matched-ledger-v1' and ledger['cutoff_UTC'] and ledger['scope'],
                'Explicit scientific ledger cutoff/scope')
        pools = self.load(ledger['resource_pools']['path'], ledger['resource_pools']['sha256'])
        guard = pools.get('guard', {})
        if isinstance(guard, str):
            guard = dict(path=guard)
        require(isinstance(guard, dict), 'Archived guard identity table')
        guard_path = guard.get('path', pools.get('guard_path'))
        guard_sha = guard.get('sha256', pools.get('guard_sha256'))
        require(isinstance(guard_path, str) and isinstance(guard_sha, str), 'Missing archived guard identity')
        if not guard_path.startswith('campaign/'):
            guard_path = str(Path(ledger['resource_pools']['path']).parent/guard_path)
        self.hash_source(guard_path, guard_sha)
        jobs = unique(ledger['jobs'], lambda row: row['label'])
        phases, reconciled, retained_jobs = {}, {}, []
        paths, log_paths = set(), set()
        for label, entry in jobs.items():
            retained = dict(entry, reconciled=False)
            retained_jobs.append(retained)
            def inspect(entry=entry, label=label, retained=retained):
                for kind in ('plan', 'final', 'receipt'):
                    require(entry[kind+'_path'] not in paths, 'Reused ledger artifact path')
                    paths.add(entry[kind+'_path'])
                plan = self.load(entry['plan_path'], entry['plan_sha256'])
                final = self.load(entry['final_path'], entry['final_sha256'])
                receipt = self.load(entry['receipt_path'], entry['receipt_sha256'])
                log_path = entry.get('log_path', 'campaign/logs/'+label+'.log')
                log = self.load(log_path, entry.get('log_sha256'), text=True) if self.path(log_path).exists() else ''
                if log:
                    require(log_path not in log_paths, 'Reused job log')
                    log_paths.add(log_path)
                retained.update(plan=plan, final=final, raw_final_accounting=receipt.get('final_accounting'),
                    raw_elapsed=receipt.get('elapsed'), settings=receipt.get('settings'),
                    cgroup_peak_bytes_upper=cgroup_peak(receipt, final, log),
                    process_measurements=receipt.get('process_measurements', {}))
                phase = entry['phase']
                raw = receipt.get('final_accounting', {}).get('CPU_seconds')
                charged, wall = number(final['charged_CPU_seconds']), number(final['absolute_guarded_wall_seconds'])
                cpu_cap, wall_cap = number(plan['attempt_CPU_seconds']), number(plan['attempt_wall_seconds'])
                if raw is not None:
                    number(raw)
                # Keep measured charges even when a later gate rejects the job.
                phases.setdefault(phase, []).append(dict(label=label, exit_code=final['exit_code'], raw_CPU_seconds=raw,
                    charged_CPU_seconds=charged, wall_seconds=wall, reserved_CPU_seconds=cpu_cap, reserved_wall_seconds=wall_cap))
                require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]{0,127}', label) is not None and
                    entry['plan_path'] == 'campaign/jobs/'+label+'.PLAN.json' and
                    entry['final_path'] == 'campaign/jobs/'+label+'.FINAL.json' and
                    entry['receipt_path'] == 'campaign/logs/'+label+'.resources.json' and
                    log_path == 'campaign/logs/'+label+'.log', 'Canonical job label/artifact paths differ')
                dispatcher = plan.get('dispatcher_archive_path', 'campaign/run_guarded_job.py')
                require(dispatcher == 'campaign/run_guarded_job.py' or
                    dispatcher == 'campaign/dispatcher_versions/'+plan['dispatcher_sha256']+'.py', 'Dispatcher archive location')
                self.hash_source(dispatcher, plan['dispatcher_sha256'])
                if 'dispatcher_path' in entry:
                    require(entry['dispatcher_path'] == dispatcher and entry['dispatcher_sha256'] == plan['dispatcher_sha256'],
                            'Ledger dispatcher identity differs')
                require(plan['label'] == final['label'] == label and plan['phase'] == final['phase'] == phase and
                    plan['pools_sha256'] == ledger['resource_pools']['sha256'] and
                    final['receipt_sha256'] == entry['receipt_sha256'], 'Ledger identity/pool binding')
                terminal = terminal_accounting(receipt, final, log, cpu_cap)
                retained['terminal_accounting'] = terminal
                if 'raw_terminal_CPU_seconds' in final:
                    phases[phase][-1]['raw_CPU_seconds'] = terminal['raw_terminal_CPU_seconds']
                    phases[phase][-1]['guard_reported_CPU_seconds'] = raw
                if 'raw_terminal_CPU_seconds' in final:
                    require('log_path' in entry and 'log_sha256' in entry, 'New FINAL needs authenticated terminal log')
                    require(final['raw_terminal_CPU_seconds'] == terminal['raw_terminal_CPU_seconds'] and
                        final['terminal_CPU_measurement_available'] is terminal['terminal_CPU_measurement_available'] and
                        final['guard_reported_CPU_seconds'] == terminal['guard_reported_CPU_seconds'], 'FINAL terminal CPU evidence differs')
                require(charged >= terminal['minimum_charged_CPU_seconds'], 'Unknown/raw CPU not conservatively charged')
                require(charged <= cpu_cap and wall <= wall_cap, 'Job exceeds reserved attempt ceiling')
                require(wall >= number(receipt['elapsed']), 'FINAL wall below guarded elapsed')
                expected_settings = dict(CPUQuota='25%', Nice='19', CPUWeight='5', IOWeight='10',
                    MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0')
                require(all(str(receipt['settings'][key]) == value for key, value in expected_settings.items()),
                        'Guard resource profile differs')
                require(type(final['exit_code']) is int and type(receipt['exit_code']) is int and
                        final['exit_code'] == receipt['exit_code'], 'Guard/FINAL exit code differs')
                require(receipt['command'] == plan['command'], 'Receipt/PLAN command differs')
                peak = cgroup_peak(receipt, final, log)
                retained['cgroup_peak_bytes_upper'] = peak
                retained['process_measurements'] = receipt.get('process_measurements', {})
                if mode == 'preflight':
                    require(peak is not None and 0 < peak <= 448*1024**2, 'Preflight cgroup peak missing/over448MiB')
                reconciled[label] = dict(entry=entry, plan=plan, final=final, receipt=receipt, log=log, terminal=terminal)
                retained['reconciled'] = True
            self.attempt('ledger job:'+label, inspect)
        attributed = set()
        for row in constructors.values():
            for boundary in ('prefix', 'continuation') if row['arm'] != 'parent' else ('prefix',):
                cost = row['costs'][boundary]
                def join(cost=cost, row=row, boundary=boundary):
                    label = cost['label']; job = reconciled[label]; entry = job['entry']
                    attributed.add(label)
                    require(all(cost[key] == entry[key] for key in ('final_path', 'final_sha256', 'receipt_path', 'receipt_sha256')),
                            'Constructor ledger artifact identity differs')
                    raw_cpu = job['final'].get('raw_terminal_CPU_seconds', job['receipt'].get('final_accounting', {}).get('CPU_seconds'))
                    require(cost['raw_CPU_seconds'] == raw_cpu and
                        cost['charged_CPU_seconds'] == job['final']['charged_CPU_seconds'] and
                        cost['wall_seconds'] == job['final']['absolute_guarded_wall_seconds'], 'Constructor costs differ from OS evidence')
                    if boundary == 'continuation' or row['arm'] == 'parent':
                        require(row['guard_exit_code'] == job['final']['exit_code'], 'Constructor/FINAL exit code differs')
                    if 'plan_path' in cost:
                        require(cost['plan_path'] == entry['plan_path'] and cost['plan_sha256'] == entry['plan_sha256'], 'Constructor PLAN identity differs')
                    if config is not None:
                        parent = constructors[(row['case'], row['repetition'], 'parent', 0)]
                        owner = parent if boundary == 'prefix' else row
                        require(not owner.get('eligible') or job['terminal']['terminal_CPU_measurement_available'],
                                'Eligible constructor lacks terminal CPU evidence')
                        original = bind_constructor(job['plan'], job['receipt'], job['final'], owner, parent,
                            owner.get('worker_verification', {}).get('worker_data'), config)
                        self.original_roots.add(original)
                before = len(self.issues)
                self.attempt('constructor ledger:'+str(order_key(row))+':'+boundary, join)
                if len(self.issues) != before:
                    row['analysis_eligible'] = False
        constructor_labels = {label for label, entry in jobs.items() if entry['phase'] == 'constructors'}
        self.attempt('constructor ledger inventory', lambda: require(attributed == constructor_labels,
            'Constructor ledger has missing/unattributed attempts: '+repr(sorted(attributed ^ constructor_labels))))
        observed = {str(path.relative_to(self.root)) for path in self.path('campaign/jobs').rglob('*.json')}
        observed.update(str(path.relative_to(self.root)) for path in self.path('campaign/logs').rglob('*.resources.json'))
        self.attempt('ledger artifact inventory', lambda: require(observed == paths, 'Unindexed/missing archived job JSON'))
        observed_logs = {str(path.relative_to(self.root)) for path in self.path('campaign/logs').rglob('*.log')}
        self.attempt('ledger log inventory', lambda: require(observed_logs == log_paths, 'Unindexed/missing archived job log'))
        self.reconciled = reconciled
        scientific = self.scientific_jobs(mode, cases) if config is not None else None
        totals = {}
        for phase, rows in phases.items():
            cpu = math.fsum(row['charged_CPU_seconds'] for row in rows)
            wall = math.fsum(row['wall_seconds'] for row in rows)
            self.attempt('phase ceiling:'+phase, lambda p=phase, c=cpu, w=wall: require(
                c <= pools['phases'][p]['CPU_seconds'] and w <= pools['phases'][p]['wall_seconds'], 'Strict phase pool exceeded'))
            totals[phase] = dict(jobs=rows, charged_CPU_seconds=cpu, wall_seconds=wall,
                known_raw_CPU_seconds=math.fsum(row['raw_CPU_seconds'] for row in rows if row['raw_CPU_seconds'] is not None),
                unknown_raw_CPU_jobs=sum(row['raw_CPU_seconds'] is None for row in rows))
        return dict(cutoff_UTC=ledger['cutoff_UTC'], scope=ledger['scope'], phases=totals,
            declared_jobs=len(jobs), reconciled_jobs=len(reconciled),
            jobs=retained_jobs,
            charged_CPU_seconds=math.fsum(row['charged_CPU_seconds'] for rows in phases.values() for row in rows),
            wall_seconds=math.fsum(row['wall_seconds'] for rows in phases.values() for row in rows),
            later_analysis_reporting_review_packaging_excluded=True,
            resource_pool_sha256=ledger['resource_pools']['sha256'], archived_guard_sha256=guard_sha,
            resource_pool_declaration=pools, scientific_output_guards=scientific)

    def scientific_jobs(self, mode, cases):
        """Bind selected scientific outputs to successful guarded argv records."""
        selected = {}
        def one(name, phase, matcher):
            matches = []
            for label, job in self.reconciled.items():
                try:
                    yes = job['plan']['phase'] == phase and matcher(job['plan'])
                except (ValueError, KeyError, TypeError, IndexError):
                    yes = False
                if yes and job['final']['exit_code'] == 0:
                    matches.append(label)
            require(len(matches) == 1, 'Missing/ambiguous successful scientific guard: '+name)
            selected[name] = matches[0]
        def python(plan, script, flags, switches=()):
            return command_flags(plan['command'], plan['cwd'], script, flags, switches)
        def julia(plan, script, arguments):
            argv = plan['command']
            expected = posixpath.normpath(plan['cwd']+'/numerical_audit/cost_matched_continuation/'+script)
            return (len(argv) == 4+len(arguments) and argv[1:3] == ['--startup-file=no', '--threads=1'] and
                posixpath.normpath(posixpath.join(plan['cwd'], argv[3])) == expected and argv[4:] == arguments)
        before = len(self.issues)
        self.attempt('original data root', lambda: require(len(self.original_roots) == 1, 'No unique authenticated original data root'))
        if len(self.original_roots) != 1:
            return dict(global_pass=False, verified_cases=[], selected_jobs=selected)
        root = next(iter(self.original_roots))
        if mode == 'preflight':
            self.attempt('input guard', lambda: one('inputs', 'inputs', lambda p:
                python(p, 'preflight.py', ('--output',)) == {'--output': root}))
        else:
            for action in ('--prepare', '--finalize'):
                self.attempt('input guard:'+action, lambda action=action: one('inputs'+action, 'inputs', lambda p:
                    python(p, 'generate_inputs.py', ('--output', '--case', '--preflight'), ('--prepare', '--finalize')) ==
                    ({'--output': root, action: True, '--preflight': self.preflight_original_root} if action == '--prepare' else
                     {'--output': root, action: True})))
            for case in cases:
                self.attempt('input guard:'+case, lambda case=case: one('inputs:'+case, 'inputs', lambda p:
                    python(p, 'generate_inputs.py', ('--output', '--case', '--preflight'), ('--prepare', '--finalize')) ==
                    {'--output': root, '--case': case}))
        self.attempt('panel guard', lambda: one('panels', 'packaging', lambda p:
            python(p, 'panels.py', ('--root',)) == {'--root': root}))
        global_pass = len(self.issues) == before
        verified_cases = []
        for case in cases:
            prior = len(self.issues)
            model, panel = root+'/inputs/'+case+'/model.toml', root+'/panels/'+case+'/candidates.toml'
            for phase, script, arguments in (
                ('physics', 'physical.jl', [model, panel, root+'/physical/'+case]),
                ('objective', 'objective.jl', [model, panel, root+'/objective/'+case]),
                ('independent', 'verify.jl', [model, panel, root+'/physical/'+case, root+'/independent/'+case])):
                self.attempt('scientific guard:'+case+':'+phase, lambda phase=phase, script=script, arguments=arguments:
                    one(case+':'+phase, phase, lambda p: julia(p, script, arguments)))
            if len(self.issues) == prior and global_pass:
                verified_cases.append(case)
        return dict(global_pass=global_pass, verified_cases=verified_cases, selected_jobs=selected,
                    original_data_root=root)


def preflight_proof(audit, frozen):
    proof = frozen['preflight']
    consumer = 'numerical_audit/cost_matched_continuation/analyze.py'
    digest = proof['consumer_sha256']
    require(proof['path'] == 'preflight' and
        proof['analysis_path'] == 'preflight/validation/ANALYSIS.'+digest+'.json' and
        proof['consumer_source_path'] == 'preflight/validation/analyze.'+digest+'.py', 'Preflight proof layout')
    pref = audit.load('preflight/INPUT_PROTOCOL_FROZEN.json', proof['source_freeze_sha256'])
    require(pref['sources'].keys() == frozen['sources'].keys(), 'Preflight producer source inventory differs')
    differences = {key: dict(original=pref['sources'].get(key), replayed=frozen['sources'].get(key))
        for key in set(pref['sources']) | set(frozen['sources']) if pref['sources'].get(key) != frozen['sources'].get(key)}
    require(set(differences) <= {consumer} and proof['consumer_only_source_delta'] == differences and
        frozen['sources'][consumer] == digest, 'Preflight source delta exceeds bound consumer repair')
    audit.hash_source(proof['consumer_source_path'], digest)
    require(pref.get('mode') == 'preflight' and all(pref[key] == frozen[key] for key in ('configuration', 'registry', 'environment')),
            'Main/preflight scientific identity differs')
    saved = audit.load(proof['analysis_path'], proof['analysis_sha256'])
    require(saved['study_complete'] is True and saved['mode'] == 'preflight' and
        saved['source_freeze_sha256'] == proof['source_freeze_sha256'] and
        saved['primary']['decision'] == 'preflight_only' and saved['analyzer_sha256'] == digest and
        saved['memory_gate']['passed'] is True and saved['memory_gate']['ceiling_MiB'] == 448,
        'Stored preflight admission incomplete')
    replay = analyze(audit.path('preflight'), preflight_only=True)
    for name, digest in replay['input_digests'].items():
        audit.inputs['preflight/'+name] = digest
    require(replay['study_complete'] is True and replay['primary']['decision'] == 'preflight_only' and
        replay['memory_gate']['passed'] is True and replay['memory_gate']['ceiling_MiB'] == 448,
        'Current consumer does not admit copied preflight')
    audit.preflight_original_root = replay['scientific_cutoff_ledger']['scientific_output_guards']['original_data_root']
    return dict(verified=True, analysis_sha256=proof['analysis_sha256'], source_freeze_sha256=proof['source_freeze_sha256'],
                recomputed_analyzer_sha256=replay['analyzer_sha256'])


def analyze(root, preflight_only=False):
    audit = Audit(root)
    frozen = audit.load('INPUT_PROTOCOL_FROZEN.json')
    source_sha = audit.inputs['INPUT_PROTOCOL_FROZEN.json']
    protocol_name = 'source_snapshot/numerical_audit/cost_matched_continuation/PROTOCOL.json'
    config = audit.load(protocol_name, frozen['configuration_sha256'])
    protocol_sha = audit.inputs[protocol_name]
    require(config['schema'] == 'cost-matched-continuation-v1', 'Protocol schema')
    for name, digest in frozen['sources'].items():
        audit.attempt('source:'+name, lambda n=name, d=digest: audit.hash_source('source_snapshot/'+n, d))
    audit.attempt('configuration', lambda: require(config == frozen['configuration'], 'Archived configuration differs'))
    def amendment_check():
        name = 'numerical_audit/cost_matched_continuation/'+config['preflight']['resource_amendment']
        require(name.endswith('/PREFLIGHT_RESOURCE_AMENDMENT.json'), 'Named preflight resource amendment')
        amendment = audit.load('source_snapshot/'+name, frozen['sources'][name])
        require(amendment['schema'] == 'cost-matched-preflight-resource-amendment-v1', 'Resource amendment schema')
        declaration = amendment['admission']
        require(all(declaration[key] == value for key, value in dict(ceiling_MiB=448, enforced_MemoryMax_MiB=512,
            MemoryHigh_MiB=448, MemorySwapMax=0, CPUQuota='25%', threads=1).items()) and
            'MemoryPeak' in declaration['metric'], 'Resource amendment bounds')
        return True
    amendment_verified = audit.attempt('preflight resource amendment', amendment_check, False)
    indices = {name: audit.attempt(name, lambda n=name: audit.load(n), {}) for name in
               ('PARENTS_FROZEN.json', 'CONSTRUCTORS_FROZEN.json', 'PANELS_FROZEN.json')}
    header = indices['PARENTS_FROZEN.json']
    mode = header.get('mode', 'main')
    require(not preflight_only or mode == 'preflight', 'Nested/main preflight proof refused')
    audit.attempt('freeze mode', lambda: require(mode == frozen.get('mode', 'main'), 'Freeze/index mode differs'))
    prerequisite = None if mode == 'preflight' else audit.attempt('main preflight prerequisite', lambda: preflight_proof(audit, frozen))
    repetitions = 1 if mode == 'preflight' else config['repetitions']
    cases = header.get('cases', [row['case'] for row in frozen['registry']['cases']])
    require(isinstance(cases, list) and cases and all(isinstance(case, str) for case in cases), 'Case roster')
    require(len(set(cases)) == len(cases), 'Duplicate case roster')
    expected_cases = 1 if mode == 'preflight' else config['inventory']['cases']
    audit.attempt('case roster', lambda: require(mode in ('main', 'preflight') and len(cases) == expected_cases and
        (set(cases) == {row['case'] for row in frozen['registry']['cases']} if mode == 'main' else
         cases == [config['preflight']['case']]), 'Declared case/mode grid'))
    inputs = audit.attempt('input freeze', lambda: audit.load('INPUTS_FROZEN.json'), {})
    def inputs_check():
        require(inputs['protocol_sha256'] == protocol_sha and inputs['source_freeze_sha256'] == source_sha,
                'Input freeze protocol/source binding')
        require(set(inputs['successes']) == set(cases) and inputs['status'] == 'complete' and
                not inputs.get('refusals'), 'Exact complete input generation grid')
        for case, files in inputs['successes'].items():
            require(set(files) == {'model.json', 'model.toml', 'orbitals.npz', 'tensors.npz'}, 'Complete frozen input file inventory')
            for name, digest in files.items():
                audit.hash_source('inputs/'+case+'/'+name, digest)
    audit.attempt('input freeze graph', inputs_check)
    for name, index in indices.items():
        audit.attempt(name+':header', lambda i=index: require(i['mode'] == mode and i['cases'] == cases and
            i['repetitions'] == repetitions and i['protocol_sha256'] == protocol_sha and
            i['source_freeze_sha256'] == source_sha and i['protocol_path'] == protocol_name and
            i['source_freeze_path'] == 'INPUT_PROTOCOL_FROZEN.json', 'Frozen index identities'))
    audit.attempt('panel freeze graph', lambda: require(all(indices['PANELS_FROZEN.json'][field] == audit.inputs[name]
        for field, name in (('input_freeze_sha256', 'INPUTS_FROZEN.json'), ('parents_freeze_sha256', 'PARENTS_FROZEN.json'),
                            ('constructors_freeze_sha256', 'CONSTRUCTORS_FROZEN.json'))), 'Panel freeze links differ'))
    expected = {(case, rep, arm, budget) for case in cases for rep in range(1, repetitions+1)
                for arm, budget in [('parent', 0)]+[(arm, budget) for arm in config['arms'] for budget in config['total_CPU_budgets']]}
    constructors, retained = {}, []
    for name in ('PARENTS_FROZEN.json', 'CONSTRUCTORS_FROZEN.json'):
        for entry in indices[name].get('entries', []):
            record = audit.attempt('constructor:'+str(entry.get('path')), lambda e=entry: audit.load(e['path'], e['sha256']))
            if record is None:
                continue
            record = dict(record, analysis_eligible=False, artifact_path=entry['path'], artifact_sha256=entry['sha256'])
            retained.append(record)
            record['worker_verification'] = audit.attempt('worker:'+entry['path'], lambda row=record: audit.worker(row),
                dict(verified=False, reason='invalid_worker_evidence'))
            def admit(e=entry, row=record):
                arm = 'parent' if name == 'PARENTS_FROZEN.json' else e['arm']
                budget = 0 if arm == 'parent' else e['total_CPU_budget']
                require(row['schema'] == 'cost-matched-constructor-v1' and
                    order_key(row) == (e['case'], e['repetition'], arm, budget), 'Constructor/index identity')
                key = order_key(row)
                require(key in expected and key not in constructors, 'Unknown/duplicate constructor')
                constructors[key] = row
                require(row['protocol_sha256'] == protocol_sha and row['source_freeze_sha256'] == source_sha,
                        'Constructor source/protocol binding')
                model = audit.load('inputs/'+row['case']+'/model.json', row['model_json_sha256'])
                toml = audit.load('inputs/'+row['case']+'/model.toml', row['model_toml_sha256'])
                require(model['case'] == row['case'] == toml['case'] and model['n'] == 14 and
                    model['electrons'] == 4 and toml['source_json_sha256'] == row['model_json_sha256'], 'Original model identity')
                order = row['order']
                require(isinstance(order, list) and all(type(value) is int for value in order) and
                    sorted(order) == list(range(len(model['blocks']))), 'Complete original-label order')
                parent = row if arm == 'parent' else constructors[(row['case'], row['repetition'], 'parent', 0)]
                require(arm == 'parent' or row['parent_sha256'] == parent['artifact_sha256'], 'Frozen parent hash differs')
                cost_check(row, parent, config)
                admission = audit.constructor_graph(row, parent)
                score = row['score']
                if score is None:
                    require(row['eligible'] is False, 'Unscored eligible constructor')
                    return
                number(score['value'])
                number(score['cubic_response'])
                require(type(score['steps']) is int and score['steps'] > 0, 'Objective affordable steps')
                row['analysis_eligible'] = (row['eligible'] is True and row['guard_exit_code'] == 0 and
                    row['controller_fallback'] is False and row['worker_verification']['verified'] is True and
                    admission['admitted'] and (arm == 'parent' or parent['analysis_eligible']))
            audit.attempt('constructor admission:'+entry['path'], admit)
    audit.attempt('constructor grid', lambda: require(set(constructors) == expected and len(retained) == len(expected), 'Missing/extra frozen orders'))
    physical, verified, physical_rows, check_reports = {}, set(), [], []
    panels = audit.attempt('panel roster', lambda: unique(indices['PANELS_FROZEN.json'].get('entries', []), lambda row: row['case']), {})
    audit.attempt('panel grid', lambda: require(set(panels) == set(cases), 'Missing/extra physical panels'))
    for case in cases:
        report = audit.attempt('physical:'+case, lambda: audit.load('physical/'+case+'/PHYSICAL.toml'), {})
        verify = audit.attempt('independent:'+case, lambda: audit.load('independent/'+case+'/VERIFY.toml'), {})
        objective = audit.attempt('objective:'+case, lambda: audit.load('objective/'+case+'/OBJECTIVE.toml'), {})
        check_reports.append(dict(case=case, physical=report, independent=verify, objective=objective))
        for row in report.get('rows', []):
            if not isinstance(row, dict):
                audit.issues.append(dict(scope='physical row:'+case, type='ValueError', message='Physical row must be a table'))
                continue
            physical_rows.append(dict(row, case=case))
            def retain(row=row):
                key = (case,)+logical(row)
                require(key not in physical, 'Duplicate physical row')
                number(row['projective_error'])
                physical[key] = row
            audit.attempt('physical row:'+case+':'+str(row.get('name')), retain)
        def inspect_case():
            case_verified = set()
            entry = panels[case]
            panel = audit.load(entry['path'], entry['sha256'])
            panel_sha = audit.inputs[str(audit.path(entry['path']).relative_to(audit.root))]
            require(panel['case'] == case and panel['mode'] == mode and panel['repetitions'] == repetitions and
                panel['protocol_sha256'] == protocol_sha and panel['source_freeze_sha256'] == source_sha, 'Panel header identity')
            roster = unique(panel['candidates'], logical)
            expected_rows = {(arm, rep, budget, duration, cap) for rep in range(1, repetitions+1)
                for arm, budget in [('parent', 0)]+[(arm, budget) for arm in config['arms'] for budget in config['total_CPU_budgets']]
                for duration in config['physical']['times'] for cap in config['physical']['rotation_cap_multipliers']}
            require(set(roster) == expected_rows and panel['expected_logical_rows'] == len(expected_rows), 'Panel logical grid')
            for candidate in roster.values():
                constructor = constructors[(case, candidate['repetition'], candidate['arm'], candidate['cpu_budget'])]
                parent = constructors[(case, candidate['repetition'], 'parent', 0)]
                require(candidate['constructor_sha256'] == constructor['artifact_sha256'] and
                    candidate['parent_sha256'] == parent['artifact_sha256'] and
                    candidate['order'] == [value+1 for value in constructor['order']] and candidate['formula'] == 2 and
                    candidate['partition'] == 'conserving' and candidate['score_kind'] == KIND and
                    candidate['claimed_score'] == constructor['score']['value'] and
                    candidate['constructor_score'] == candidate['claimed_score'], 'Panel constructor graph')
                require(panel['input_sha256'] == constructor['model_toml_sha256'] and
                    panel['source_json_sha256'] == constructor['model_json_sha256'], 'Panel model binding')
            for producer in (report, verify):
                require(all(producer[key] == panel[key] for key in ('case', 'mode', 'repetitions', 'expected_logical_rows',
                    'input_sha256', 'source_json_sha256', 'protocol_sha256', 'source_freeze_sha256')) and
                    producer['candidates_sha256'] == panel_sha, 'Physical/check panel binding')
            for producer in (report, verify, objective):
                require(producer['input_sha256'] == panel['input_sha256'] and producer['candidates_sha256'] == panel_sha,
                        'Objective/physical input binding')
            audit.source_hashes(report['source_hashes'], frozen)
            audit.source_hashes(verify['checker_source_hashes'], frozen)
            require(verify['primary_source_hashes'] == report['source_hashes'], 'Independent primary source binding')
            audit.source_hashes(objective['source_hashes'], frozen, Path('numerical_audit/prospective_window_ordering'))
            require(verify['physical_sha256'] == audit.inputs['physical/'+case+'/PHYSICAL.toml'], 'Independent physical file binding')
            observed = unique(report['rows'], logical)
            checked = unique(verify['rows'], logical)
            require(set(observed) == expected_rows == set(checked), 'Physical/independent logical grid')
            objective_rows = unique(objective['logical_rows'], lambda row: row['name'])
            require(set(objective_rows) == {row['name'] for row in roster.values()}, 'Objective logical grid')
            require(objective['absolute_tolerance'] == 1e-24 and objective['relative_tolerance'] == 1e-8,
                    'Fixed objective replay tolerances')
            used_checks = set()
            for key, candidate in roster.items():
                row, checked_row = observed[key], checked[key]
                require(all(row[field] == candidate[field] for field in IDENTITY), 'Physical/candidate row identity')
                require(checked_row['name'] == candidate['name'] and checked_row['order_id'] == candidate['order_id'] and
                    checked_row['status'] == 'complete' and checked_row['independent_unfused_replay'] is True and
                    checked_row['actual_sha256'] == row['actual_sha256'], 'Independent physical row not verified')
                for file_key in ('actual', 'reference'):
                    audit.load('physical/'+case+'/'+row[file_key+'_file'], row[file_key+'_sha256'], binary=True)
                info = objective_rows[candidate['name']]
                require(info['candidate_index'] == panel['candidates'].index(candidate)+1, 'Objective candidate index')
                index = info['check_index']
                require(type(index) is int and 1 <= index <= len(objective['checks']), 'Objective check index')
                used_checks.add(index)
                check = objective['checks'][index-1]
                require(check['check_index'] == index and check['verified'] is True and info['status'] == 'verified' and
                    check['order'] == candidate['order'] and check['score_kind'] == candidate['score_kind'] and
                    check['constructor_score'] == candidate['claimed_score'] and check['partition'] == 'conserving', 'Objective witness not verified')
                require(check['steps'] == constructors[(case, candidate['repetition'], candidate['arm'], candidate['cpu_budget'])]['score']['steps'],
                        'Objective affordable steps differ')
                require(abs(number(check['replayed_score'])-candidate['claimed_score']) <=
                    objective['absolute_tolerance']+objective['relative_tolerance']*max(check['replayed_score'], candidate['claimed_score']), 'Objective replay score differs')
                case_verified.add((case,)+key)
            require(used_checks == set(range(1, len(objective['checks'])+1)), 'Unlinked objective checks')
            require(report['status'] == verify['status'] == 'complete' and report['physical_panel_complete'] is True and
                not report.get('failures') and not verify.get('failures') and not verify.get('censored_rows') and
                report['completed_logical_rows'] == verify['completed_logical_rows'] == len(expected_rows), 'Physical panel/check incomplete')
            require(verify['reference_check_count'] == 2 and len(verify['reference_checks']) == 2 and
                {row['time'] for row in verify['reference_checks']} == set(config['physical']['times']) and
                all(row['status'] == 'complete' for row in verify['reference_checks']), 'Independent references incomplete')
            require(report['unique_actual_propagations'] == verify['unique_circuit_count'] ==
                verify['independent_unfused_replay_count'] and verify['state_tolerance'] == 1e-10 and
                verify['absolute_metric_tolerance'] == 1e-12, 'Independent replay count/tolerance differs')
            require(objective['complete'] is True and objective['status'] == 'complete_for_supplied_scored_orders' and
                objective['all_scored_objectives_verified'] is True and objective['mismatched_objectives'] == 0 and
                objective['processed_logical_rows'] == objective['expected_logical_rows'] == len(expected_rows), 'Objective verification incomplete')
            return case_verified
        verified.update(audit.attempt('case:'+case, inspect_case, set()))
    expected_physical = len(expected)*len(config['physical']['times'])*len(config['physical']['rotation_cap_multipliers'])
    audit.attempt('physical inventory', lambda: require(len(physical) == expected_physical, 'Missing/extra physical rows'))
    ledger = audit.attempt('scientific cutoff ledger', lambda: audit.ledger(constructors, config, mode, cases))
    guards = {} if ledger is None else ledger.get('scientific_output_guards', {})
    verified = {key for key in verified if guards.get('global_pass') and key[0] in guards.get('verified_cases', [])}
    for row in constructors.values():
        if row['arm'] != 'parent' and not constructors.get((row['case'], row['repetition'], 'parent', 0), {}).get('analysis_eligible'):
            row['analysis_eligible'] = False
    pairs = pair_grid(cases, repetitions, config, physical, constructors, verified)
    resources = audit.attempt('resource totals', lambda: resource_totals(constructors), {})
    resources['inventory_complete'] = set(constructors) == expected and len(retained) == len(expected)
    all_gates = not audit.issues and all(row['analysis_eligible'] for row in constructors.values())
    primary = primary_decision(pairs, cases, repetitions, config, all_gates)
    if mode == 'preflight':
        primary['decision'] = 'preflight_only' if all_gates else 'study_incomplete'
    secondary = []
    for budget in config['total_CPU_budgets']:
        rows = [row for row in pairs if row['total_CPU_budget'] == budget]
        secondary.append(dict(total_CPU_budget=budget, paired_rows=len(rows),
            counts={outcome: sum(row['outcome'] == outcome for row in rows) for outcome in ('win', 'tie', 'loss', 'unavailable')},
            median_paired_difference=statistics.median([row['difference'] for row in rows if row['difference'] is not None])
            if any(row['difference'] is not None for row in rows) else None))
    jobs = [] if ledger is None else ledger['jobs']
    peaks = [job.get('cgroup_peak_bytes_upper') for job in jobs]
    memory_pass = (amendment_verified and bool(jobs) and all(value is not None and 0 < value <= 448*1024**2 for value in peaks)) if mode == 'preflight' else (
        bool(jobs) and all(job.get('settings', {}).get('MemoryMax') == '512M' and
                          job.get('settings', {}).get('MemorySwapMax') == '0' for job in jobs))
    return dict(schema='cost-matched-continuation-analysis-v1', status='complete' if all_gates else 'qualified',
        study_complete=all_gates, mode=mode, cases=cases, repetitions=repetitions,
        analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), protocol_sha256=protocol_sha,
        source_freeze_sha256=source_sha, input_digests=audit.inputs, issues=audit.issues,
        expected_frozen_orders=len(expected), retained_constructor_count=len(retained), constructors=retained,
        expected_logical_physical_rows=expected_physical, physical_rows=physical_rows,
        paired_rows=pairs, primary=primary, secondary=secondary, resources=resources,
        resource_scope='Frozen constructor records: actual shared executions versus attributed comparisons; missing inventory is qualified; physical rows never multiply costs',
        verification_reports=check_reports, scientific_cutoff_ledger=ledger, input_freeze=inputs,
        preflight_prerequisite=prerequisite, memory_gate=dict(scope='preflight' if mode == 'preflight' else 'main_hard_profile',
            ceiling_MiB=448 if mode == 'preflight' else 512, passed=bool(memory_pass),
            maximum_observed_cgroup_peak_bytes_upper=max((value for value in peaks if value is not None), default=None)))


def write_new(path, result):
    path = Path(path)
    temporary = path.with_name(path.name+'.partial')
    require(not path.exists() and not temporary.exists(), 'Preserve existing output/partial')
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    import os
    os.link(temporary, path)
    temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.root)
    write_new(args.output, result)
    return 0 if result['study_complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
