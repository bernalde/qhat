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
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
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
            number(costs[boundary][field])
        require(costs[boundary]['charged_CPU_seconds'] >= costs[boundary]['raw_CPU_seconds'],
                'Conservative CPU below raw CPU')
    for field in COST_FIELDS:
        require(costs['prefix'][field] == parent['costs']['total'][field],
                'Shared parent prefix differs')
        require(math.isclose(costs['total'][field], costs['prefix'][field] +
                             costs['continuation'][field], rel_tol=0., abs_tol=1e-9),
                'Prefix/continuation total differs')
    cap = config['parent']['CPU_seconds'] if record['arm'] == 'parent' else record['total_CPU_budget']
    require(costs['total']['charged_CPU_seconds'] <= cap and
            costs['total']['raw_CPU_seconds'] <= cap, 'Strict whole-construction CPU ceiling')
    if record['arm'] == 'parent':
        require(all(costs['continuation'][field] == 0 for field in COST_FIELDS), 'Parent continuation is not zero')
        require(costs['total']['wall_seconds'] <= config['parent']['wall_seconds'], 'Parent wall ceiling')
    else:
        require(costs['total']['wall_seconds'] <= config['accounting']['total_wall_multiplier']*cap,
                'Strict whole-construction wall ceiling')


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

    def attempt(self, scope, function, default=None):
        try:
            return function()
        except (ValueError, KeyError, TypeError, IndexError, OSError, tomllib.TOMLDecodeError) as error:
            self.issues.append(dict(scope=scope, type=type(error).__name__, message=str(error)))
            return default

    def path(self, name):
        path = (self.root/name).resolve()
        require(path.is_relative_to(self.root), 'Path outside retained run')
        return path

    def load(self, name, digest=None, binary=False):
        path = self.path(name)
        require(path.is_file() and path.stat().st_size <= 32*1024**2, 'Missing/over-cap input: '+str(name))
        raw = path.read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        self.inputs[str(path.relative_to(self.root))] = actual
        require(digest is None or digest == actual, 'Input digest differs: '+str(name))
        if binary:
            require(len(raw) == 16*16384, 'Full-register state byte count')
            return raw
        return tomllib.loads(raw.decode()) if path.suffix == '.toml' else json.loads(raw)

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
            return dict(verified=False, reason='missing_worker_output', path=path)
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
                        path=path, sha256=digest, refusal=worker.get('refusal'))
        require(isinstance(result, dict) and isinstance(result.get('score'), dict), 'Unscored worker without fallback')
        require(result['order'] == constructor['order'] and result['score'] == constructor['score'],
                'Worker/canonical selected order or score differs')
        return dict(verified=True, reason=None, path=path, sha256=digest)

    def ledger(self, constructors):
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
        paths = set()
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
                retained.update(plan=plan, final=final, raw_final_accounting=receipt.get('final_accounting'),
                    raw_elapsed=receipt.get('elapsed'), settings=receipt.get('settings'))
                phase = entry['phase']
                raw = receipt.get('final_accounting', {}).get('CPU_seconds')
                charged, wall = number(final['charged_CPU_seconds']), number(final['absolute_guarded_wall_seconds'])
                cpu_cap, wall_cap = number(plan['attempt_CPU_seconds']), number(plan['attempt_wall_seconds'])
                if raw is not None:
                    number(raw)
                # Keep measured charges even when a later gate rejects the job.
                phases.setdefault(phase, []).append(dict(label=label, exit_code=final['exit_code'], raw_CPU_seconds=raw,
                    charged_CPU_seconds=charged, wall_seconds=wall, reserved_CPU_seconds=cpu_cap, reserved_wall_seconds=wall_cap))
                self.hash_source('campaign/run_guarded_job.py', plan['dispatcher_sha256'])
                require(plan['label'] == final['label'] == label and plan['phase'] == final['phase'] == phase and
                    plan['pools_sha256'] == ledger['resource_pools']['sha256'] and
                    final['receipt_sha256'] == entry['receipt_sha256'], 'Ledger identity/pool binding')
                require(charged >= (cpu_cap if raw is None else number(raw)), 'Unknown/raw CPU not conservatively charged')
                require(charged <= cpu_cap and wall <= wall_cap, 'Job exceeds reserved attempt ceiling')
                require(wall >= number(receipt['elapsed']), 'FINAL wall below guarded elapsed')
                expected_settings = dict(CPUQuota='25%', Nice='19', CPUWeight='5', IOWeight='10',
                    MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0')
                require(all(str(receipt['settings'][key]) == value for key, value in expected_settings.items()),
                        'Guard resource profile differs')
                require(final['exit_code'] == receipt['exit_code'], 'Guard/FINAL exit code differs')
                reconciled[label] = dict(entry=entry, final=final, receipt=receipt)
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
                    require(cost['raw_CPU_seconds'] == job['receipt']['final_accounting']['CPU_seconds'] and
                        cost['charged_CPU_seconds'] == job['final']['charged_CPU_seconds'] and
                        cost['wall_seconds'] == job['final']['absolute_guarded_wall_seconds'], 'Constructor costs differ from OS evidence')
                    if boundary == 'continuation' or row['arm'] == 'parent':
                        require(row['guard_exit_code'] == job['final']['exit_code'], 'Constructor/FINAL exit code differs')
                self.attempt('constructor ledger:'+str(order_key(row))+':'+boundary, join)
        constructor_labels = {label for label, entry in jobs.items() if entry['phase'] == 'constructors'}
        self.attempt('constructor ledger inventory', lambda: require(attributed == constructor_labels,
            'Constructor ledger has missing/unattributed attempts: '+repr(sorted(attributed ^ constructor_labels))))
        observed = {str(path.relative_to(self.root)) for path in self.path('campaign/jobs').rglob('*.json')}
        observed.update(str(path.relative_to(self.root)) for path in self.path('campaign/logs').rglob('*.resources.json'))
        self.attempt('ledger artifact inventory', lambda: require(observed == paths, 'Unindexed/missing archived job JSON'))
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
            resource_pool_declaration=pools)


def analyze(root):
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
    indices = {name: audit.attempt(name, lambda n=name: audit.load(n), {}) for name in
               ('PARENTS_FROZEN.json', 'CONSTRUCTORS_FROZEN.json', 'PANELS_FROZEN.json')}
    header = indices['PARENTS_FROZEN.json']
    mode = header.get('mode', 'main')
    repetitions = 1 if mode == 'preflight' else config['repetitions']
    cases = header.get('cases', [row['case'] for row in frozen['registry']['cases']])
    require(isinstance(cases, list) and cases and all(isinstance(case, str) for case in cases), 'Case roster')
    require(len(set(cases)) == len(cases), 'Duplicate case roster')
    expected_cases = 1 if mode == 'preflight' else config['inventory']['cases']
    audit.attempt('case roster', lambda: require(mode in ('main', 'preflight') and len(cases) == expected_cases and
        (mode != 'main' or set(cases) == {row['case'] for row in frozen['registry']['cases']}), 'Declared case/mode grid'))
    if audit.path('INPUTS_FROZEN.json').exists():
        def inputs_check():
            inputs = audit.load('INPUTS_FROZEN.json')
            require(inputs['protocol_sha256'] == protocol_sha and inputs['source_freeze_sha256'] == source_sha,
                    'Input freeze protocol/source binding')
            for case, files in inputs.get('successes', {}).items():
                for name, digest in files.items():
                    audit.hash_source('inputs/'+case+'/'+name, digest)
            require(inputs['status'] == 'complete' and not inputs.get('refusals'), 'Input generation incomplete/refused')
        audit.attempt('input freeze', inputs_check)
    for name, index in indices.items():
        audit.attempt(name+':header', lambda i=index: require(i['mode'] == mode and i['cases'] == cases and
            i['repetitions'] == repetitions and i['protocol_sha256'] == protocol_sha and
            i['source_freeze_sha256'] == source_sha and i['protocol_path'] == protocol_name and
            i['source_freeze_path'] == 'INPUT_PROTOCOL_FROZEN.json', 'Frozen index identities'))
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
                score = row['score']
                require(isinstance(score, dict), 'Unscored constructor')
                number(score['value'])
                number(score['cubic_response'])
                require(type(score['steps']) is int and score['steps'] > 0, 'Objective affordable steps')
                parent = row if arm == 'parent' else constructors[(row['case'], row['repetition'], 'parent', 0)]
                require(arm == 'parent' or row['parent_sha256'] == parent['artifact_sha256'], 'Frozen parent hash differs')
                cost_check(row, parent, config)
                row['analysis_eligible'] = (row['eligible'] is True and row['guard_exit_code'] == 0 and
                    row['controller_fallback'] is False and row['worker_verification']['verified'] is True)
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
            physical_rows.append(dict(row, case=case))
            def retain(row=row):
                key = (case,)+logical(row)
                require(key not in physical, 'Duplicate physical row')
                number(row['projective_error'])
                physical[key] = row
            audit.attempt('physical row:'+case+':'+str(row.get('name')), retain)
        def inspect_case():
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
                verified.add((case,)+key)
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
        audit.attempt('case:'+case, inspect_case)
    expected_physical = len(expected)*len(config['physical']['times'])*len(config['physical']['rotation_cap_multipliers'])
    audit.attempt('physical inventory', lambda: require(len(physical) == expected_physical, 'Missing/extra physical rows'))
    pairs = pair_grid(cases, repetitions, config, physical, constructors, verified)
    ledger = audit.attempt('scientific cutoff ledger', lambda: audit.ledger(constructors))
    resources = audit.attempt('resource totals', lambda: resource_totals(constructors), {})
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
    return dict(schema='cost-matched-continuation-analysis-v1', status='complete' if all_gates else 'qualified',
        study_complete=all_gates, mode=mode, cases=cases, repetitions=repetitions,
        analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), protocol_sha256=protocol_sha,
        source_freeze_sha256=source_sha, input_digests=audit.inputs, issues=audit.issues,
        expected_frozen_orders=len(expected), retained_constructor_count=len(retained), constructors=retained,
        expected_logical_physical_rows=expected_physical, physical_rows=physical_rows,
        paired_rows=pairs, primary=primary, secondary=secondary, resources=resources,
        resource_scope='Actual shared constructor executions versus attributed continuation comparisons; physical rows never multiply construction costs',
        verification_reports=check_reports, scientific_cutoff_ledger=ledger)


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
