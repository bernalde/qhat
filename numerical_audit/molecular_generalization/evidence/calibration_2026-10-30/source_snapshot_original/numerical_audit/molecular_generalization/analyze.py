"""Join every frozen row and report the two predeclared paired questions.

No numerical kernel is imported or executed. Missing, refused and floor-limited
records remain explicit. Main decisions require the complete declared evidence.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import statistics
import sys

if __package__:
    from . import common as k, panels as p
else:
    import common as k
    import panels as p

AMBIGUITY = 1e-12
require = k.require


def unique(rows, key):
    result = {}
    require(isinstance(rows, list), 'Expected row inventory')
    for row in rows:
        require(isinstance(row, dict), 'Expected row table')
        identity = key(row)
        require(identity not in result, 'Duplicate logical identity')
        result[identity] = row
    return result


def logical(row):
    return tuple(row[key] for key in ('arm', 'repetition', 'cpu_budget', 'time', 'budget_multiplier'))


def comparison(actual, claimed):
    actual, claimed = p.number(actual), p.number(claimed)
    tolerance = 1e-24+1e-8*max(actual, claimed)
    within = abs(actual-claimed) <= tolerance
    floor = max(actual, claimed) <= 1e-24
    return dict(verified=within and not floor, within_tolerance=within, numerical_floor=floor,
                tolerance=tolerance, absolute_difference=abs(actual-claimed))


def source_join(values, freeze):
    require(isinstance(values, dict) and values, 'Nonempty producer source map')
    require(all(name in freeze['sources'] and freeze['sources'][name] == p.digest(sha)
                for name, sha in values.items()), 'Original producer bytes differ from immutable freeze')


def _bound_hashes(a, info, inventory):
    names = (info['protocol_path'], info['source_freeze_path'], f"inputs/{inventory['case']}/model.json")
    paths = list(names)+[value['artifact_path'] for value in inventory['records'].values()]
    return {name: a.hashes.get(name) or k.sha(a.path(name)) for name in paths}


def _state(a, case, filename, sha, n):
    raw = a.raw(f'physical/{case}/'+filename, p.digest(sha), maximum=16*(1 << n))
    require(len(raw) == 16*(1 << n), 'Full-register complex128 state byte count')


def _objective_info(info, expected, checks, used, units):
    require(all(info.get(key) == value for key, value in expected.items()), 'Objective frozen row/witness identity')
    index = info['check_index']
    require(type(index) is int and 1 <= index <= len(checks), 'Objective numerical check index')
    check = checks[index-1]; used.add(index)
    require(check['status'] == 'replayed' and check['order'] == expected['order'] and
            check['score_kind'] == k.SCORE_KIND, 'Objective archived check identity')
    order = [v-1 for v in expected['order']]
    lengths = [len(units[i]) for i in order]; total = sum(lengths)
    steps = 1 if len(lengths) == 1 else (30*total-lengths[0])//(2*total-lengths[0]-lengths[-1])
    require(check['steps'] == steps, 'Objective original-word affordable steps')
    proposed = comparison(check['replayed_score'], expected['claimed_score'])
    require(info['proposed_comparison']['verified'] == proposed['verified'] and
            info['proposed_comparison']['within_tolerance'] == proposed['within_tolerance'] and
            info['proposed_comparison']['numerical_floor'] == proposed['numerical_floor'], 'Recomputed proposed-score comparison')
    if 'native_refusal' in expected:
        require(info['status'] == 'native_refusal_censored', 'Raw native refusal remains censored')
        return False
    native_score = expected.get('native_score', expected['claimed_score'])
    native = comparison(check['replayed_score'], native_score)
    require(info['native_score'] == native_score and info['native_comparison']['verified'] == native['verified'] and
            info['native_comparison']['within_tolerance'] == native['within_tolerance'] and
            info['native_comparison']['numerical_floor'] == native['numerical_floor'], 'Recomputed native-score comparison')
    expected_status = 'verified' if proposed['verified'] and native['verified'] else (
        'mismatch' if not proposed['within_tolerance'] or not native['within_tolerance'] else 'floor_censored')
    require(info['status'] == expected_status, 'Objective verification/floor/mismatch status')
    return expected_status == 'verified'


def inspect_case(a, freeze, info, inventory, panel, witnesses, reports):
    case = inventory['case']; physical, independent, objective = (reports[name] for name in ('physical', 'independent', 'objective'))
    panel_sha = a.hashes[f'panels/{case}.toml']; expected_bound = _bound_hashes(a, info, inventory)
    header = ('case', 'mode', 'repetitions', 'expected_logical_rows', 'n', 'electrons', 'partition',
              'input_sha256', 'source_json_sha256', 'protocol_sha256', 'source_freeze_sha256',
              'order_index_base', 'training_time', 'training_budget_multiplier', 'times', 'rotation_cap_multipliers')
    for report in reports.values():
        require(all(report[key] == panel[key] for key in header) and report['candidates_sha256'] == panel_sha,
                'Report exact candidate/model/protocol/grid headers')
        require(report['artifact_binding_status'] == 'actual_bytes_bound' and report['bound_artifact_hashes'] == expected_bound,
                'Report actual-byte artifact graph')
    source_join(physical['source_hashes'], freeze)
    source_join(independent['checker_source_hashes'], freeze)
    source_join(objective['consumer_source_hashes'], freeze)
    require(independent['primary_source_hashes'] == physical['source_hashes'] and
            independent['physical_sha256'] == a.hashes[f'physical/{case}/PHYSICAL.toml'], 'Independent exact primary bytes')
    require(physical['full_dimension'] == 1 << panel['n'] and physical['hf'] == str((1 << panel['electrons'])-1) and
            physical['reference_sector_dimension'] == math.comb(panel['n'], panel['electrons']), 'Full register and complete fixed-N reference')
    expected_eligible = all(record['eligible'] for record in inventory['records'].values())
    require(physical['all_constructor_eligible'] is expected_eligible and independent['all_constructor_eligible'] is expected_eligible,
            'Constructor eligibility summary joins')
    require(physical['ratio_floor'] == AMBIGUITY and independent['state_tolerance'] == 1e-10 and
            independent['absolute_metric_tolerance'] == 1e-12 and objective['absolute_tolerance'] == 1e-24 and
            objective['relative_tolerance'] == 1e-8, 'Unchanged physical/objective tolerance and floor policy')
    require(p.number(physical['reference_outside_coupling_max']) <= p.number(physical['reference_coupling_tolerance']) and
            p.number(physical['reference_hermiticity_residual']) <= physical['reference_coupling_tolerance'], 'Checked primary reference admission')
    roster = unique(panel['candidates'], logical)
    observed = unique(physical['rows'], logical); checked = unique(independent['rows'], logical)
    scored = unique(objective['logical_rows'], lambda row: row['name'])
    require(set(observed) == set(checked) == set(roster) and
            set(scored) == {row['name'] for row in roster.values()}, 'Every logical physical/replay/objective observation')
    reference_rows = unique(physical['reference_records'], lambda row: row['time'])
    reference_checks = unique(independent['reference_checks'], lambda row: row['time'])
    require(set(reference_rows) == set(reference_checks) == {.65, 1.35}, 'Both saved and independent reference times')
    for duration, reference in reference_rows.items():
        replay = reference_checks[duration]
        require(reference['status'] == replay['status'] == 'complete' and
                reference['sector_dimension'] == math.comb(panel['n'], panel['electrons']) and
                replay['diagnostic']['sector_dimension'] == reference['sector_dimension'] and
                all(replay[key] == reference[key] for key in ('reference_file', 'reference_sha256')),
                'Reference file, time and complete-sector replay joins')
        require(p.number(replay['state_norm_residual']) <= 1e-10 and p.number(replay['phase_aligned_residual']) <= 1e-10 and
                p.number(reference['number_leakage']) <= 1e-12 and p.number(reference['spin_leakage']) <= 1e-10 and
                p.number(reference['norm_drift']) <= 1e-9, 'Reference residual/norm/leakage checks')
        _state(a, case, reference['reference_file'], reference['reference_sha256'], panel['n'])
    circuit_artifacts = {}; file_circuits = {}; used = set(); good = set(); row_issues = []
    for identity, expected in roster.items():
        try:
            row, replay, score = observed[identity], checked[identity], scored[expected['name']]
            for consumer in (row, replay, score):
                require(all(consumer.get(key) == value for key, value in expected.items()), 'Complete logical constructor/grid/order join')
            require(row['status'] == replay['status'] == 'complete' and replay['independent_unfused_replay'] is True,
                    'Every saved physical row independently replayed')
            order = [v-1 for v in expected['order']]; lengths = [len(inventory['units'][i]) for i in order]
            total = sum(lengths); budget = expected['budget_multiplier']*total
            steps = 1 if len(lengths) == 1 else (budget-lengths[0])//(2*total-lengths[0]-lengths[-1])
            rotations = total if len(lengths) == 1 else steps*(2*total-lengths[0]-lengths[-1])+lengths[0]
            require(row['rotation_budget'] == budget and row['steps'] == steps and row['ideal_rotations'] == rotations,
                    'Unchanged original-word physical circuit accounting')
            require(all(replay[key] == row[key] for key in ('actual_sha256', 'reference_file', 'reference_sha256', 'circuit_sha256')),
                    'Independent state/circuit identity')
            reference = reference_rows[expected['time']]
            require(all(row[key] == reference[key] for key in ('reference_file', 'reference_sha256')), 'Logical row exact-reference join')
            circuit = (tuple(order), expected['time'], steps)
            artifact = (row['actual_file'], row['actual_sha256'])
            duplicate = circuit in circuit_artifacts
            require(not duplicate or circuit_artifacts[circuit] == artifact, 'Exact circuit has one saved artifact')
            require(row['actual_file'] not in file_circuits or file_circuits[row['actual_file']] == circuit,
                    'Distinct circuits cannot share a saved artifact identity')
            require(row['propagation_reused_for_identical_identity'] is duplicate, 'Exact circuit computation dedup flag')
            circuit_artifacts[circuit] = artifact; file_circuits[row['actual_file']] = circuit
            _state(a, case, row['actual_file'], row['actual_sha256'], panel['n'])
            require(all(p.number(row[key]) <= 1e-10 for key in ('number_leakage', 'spin_leakage')) and
                    p.number(row['norm_drift']) <= 1e-9 and p.number(row['projective_error']) <= 1 and
                    p.number(row['infidelity']) <= 1 and p.number(row['overlap_abs']) <= 1+1e-12 and
                    type(row['raw_normalized_overlap_deficit']) in (int, float) and math.isfinite(row['raw_normalized_overlap_deficit']) and
                    row['numerical_floor'] is (row['projective_error'] <= AMBIGUITY), 'Raw physical metrics/floors retained')
            require(p.number(replay['unfused_state_norm_residual']) <= 1e-10 and
                    p.number(replay['unfused_phase_aligned_residual']) <= 1e-10, 'Fresh unfused state residual')
            for name in ('saved_metric_discrepancies', 'independent_reference_metric_discrepancies',
                         'unfused_saved_reference_metric_discrepancies', 'unfused_independent_reference_metric_discrepancies'):
                require(set(replay[name]) == {'projective_error', 'infidelity', 'norm_drift', 'number_leakage',
                        'spin_leakage', 'reference_norm_drift'} and all(p.number(v) <= 1e-12 for v in replay[name].values()),
                        'Every metric independently joined with fixed absolute tolerance')
            objective_good = _objective_info(score, expected, objective['checks'], used, inventory['units'])
            constructor = inventory['records'][(expected['repetition'], expected['arm'], expected['cpu_budget'])]
            if objective_good and constructor['analysis_eligible']:
                good.add(identity)
        except Exception as error:
            row_issues.append(dict(name=expected['name'], type=type(error).__name__, message=str(error)))
    raw_expected = unique(witnesses['witnesses'], lambda row: row['name'])
    raw_observed = unique(objective['raw_witnesses'], lambda row: row['name'])
    require(set(raw_expected) == set(raw_observed), 'Every persisted raw proposed-best claim is represented')
    all_raw_verified = True
    for name, expected in raw_expected.items():
        all_raw_verified = _objective_info(raw_observed[name], expected, objective['checks'], used, inventory['units']) and all_raw_verified
    require(used == set(range(1, len(objective['checks'])+1)), 'No dropped or unlinked objective numerical checks')
    require(objective['witnesses_sha256'] == a.hashes[f'witnesses/{case}.toml'] and
            objective['witness_inventory_complete'] is True and objective['expected_raw_witnesses'] == len(raw_expected) and
            objective['processed_raw_witnesses'] == len(raw_expected), 'Independent raw inventory and objective witness bytes')
    require(physical['unique_actual_propagations'] == independent['unique_circuit_count'] ==
            independent['independent_unfused_replay_count'] == len(circuit_artifacts), 'Every distinct exact circuit replayed')
    expected_files = set(file_circuits) | {v['reference_file'] for v in reference_rows.values()}
    actual_files = {file.name for file in (a.root/'physical'/case).glob('*.bin')}
    require(expected_files == actual_files, 'No missing/extra saved state binaries')
    phase_complete = (not row_issues and all_raw_verified and all(report.get('status') == 'complete' for report in reports.values()) and
        physical.get('physical_panel_complete') is True and objective.get('complete') is True and
        objective.get('all_scored_objectives_verified') is True and
        all(not report.get('failures') for report in reports.values()) and not independent.get('censored_rows') and
        physical['completed_logical_rows'] == independent['completed_logical_rows'] == objective['processed_logical_rows'] == len(roster) and
        independent['reference_check_count'] == 2 and len(good) == len(roster))
    return dict(complete=phase_complete, good=good, issues=row_issues, unique_circuits=len(circuit_artifacts),
                raw_claims_verified=all_raw_verified, rows=observed)



def phase_receipts(a, freeze, info):
    """Authenticate whole-process evidence for each scientific replay phase."""
    retained = []; verified = set()
    expected = {(phase, case) for phase in ('physical', 'independent', 'objective') for case in info['cases']}
    for path in sorted((a.root/'jobs').glob('*.PLAN.json')):
        relative = str(path.relative_to(a.root)); plan = a.json(relative)
        if plan.get('phase') not in ('physical', 'independent', 'objective'):
            continue
        phase = plan['phase']; case = plan['identity']['case']; identity = (phase, case)
        require(identity in expected and identity not in verified, 'Exactly one owning scientific job per phase/case')
        require(plan['schema'] == 'molecular-job-plan-v1' and plan['source_freeze_sha256'] == info['source_freeze_sha256'] and
                plan['producer_root'] == freeze['producer_root'] and plan['producer_data_root'] == freeze['producer_data_root'],
                'Scientific owning job/source roots')
        label = plan['label']
        require(relative == 'jobs/'+label+'.PLAN.json', 'Canonical scientific PLAN label')
        receipt_name = 'logs/'+label+'.resources.json'; receipt = a.json(receipt_name)
        final = a.json('jobs/'+label+'.FINAL.json')
        require(final['schema'] == 'molecular-job-final-v1' and final['plan_sha256'] == a.hashes[relative] and
                final['receipt_path'] == receipt_name and final['receipt_sha256'] == a.hashes[receipt_name] and
                final['phase'] == phase and final['label'] == label and final['guard_supervisor_exit_code'] == 0,
                'Scientific post-teardown FINAL/PLAN/receipt chain')
        require(receipt['command'] == plan['command'] and receipt['budget'] == dict(cpu_seconds=plan['cpu_cap'],
                wall_seconds=plan['wall_cap'], disk_mib=plan['disk_mib']), 'Scientific executed command and cap joins')
        source = freeze['producer_root']; data = freeze['producer_data_root']; command = plan['command']
        script = source+'/numerical_audit/molecular_generalization/'+('verify.jl' if phase == 'independent' else phase+'.jl')
        require(script in command and command.count(script) == 1, 'Exact pinned scientific Julia script')
        position = command.index(script)
        arguments = [f'{data}/inputs/{case}/model.toml', f'{data}/panels/{case}.toml']
        if phase == 'independent':
            arguments += [f'{data}/physical/{case}', f'{data}/independent/{case}']
        else:
            arguments += [f'{data}/{phase}/{case}']
            if phase == 'objective':
                arguments += [f'{data}/witnesses/{case}.toml']
        require(command[position+1:] == arguments and '--threads=1' in command[:position] and
                '--project='+source+'/analysis/julia_window' in command[:position] and
                command[0] == freeze['runtime']['julia']['executable'], 'Actual full-register CLI and pinned Julia runtime')
        a.raw('source_snapshot/numerical_audit/molecular_generalization/tools/resource_guard.py', plan['guard_sha256'])
        require(all(str(receipt['settings'][key]) == value for key, value in dict(CPUQuota='25%', CPUWeight='5',
            IOWeight='10', Nice='19', MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0').items()), 'Scientific resource profile')
        costs = p.measured_receipt(receipt)
        require(all(final[key] == value for key, value in costs.items()) and costs['charged_CPU_seconds'] is not None and
                costs['charged_CPU_seconds'] <= p.number(plan['cpu_cap']) and costs['wall_seconds'] <= p.number(plan['wall_cap']) and
                receipt['exit_code'] == 0 and receipt.get('termination') is None and receipt.get('guard_error') is None and
                not receipt.get('cleanup_errors'), 'Scientific complete cost/termination evidence')
        retained.append(dict(phase=phase, case=case, plan=relative, receipt=receipt_name, **costs)); verified.add(identity)
    require(verified == expected, 'Missing complete scientific phase receipts')
    return retained


def paired_grid(cases, repetitions, physical, eligible):
    pairs = []
    for case in cases:
        for rep in range(1, repetitions+1):
            for budget in (32, 64):
                for duration in (.65, 1.35):
                    for cap in (8, 24):
                        for arm in ('window', 'annealing'):
                            left = (case, arm, rep, budget, duration, cap)
                            right = (case, 'full', rep, budget, duration, cap)
                            candidate, full = physical.get(left), physical.get(right)
                            raw = None
                            if candidate is not None and full is not None:
                                try:
                                    raw = p.number(candidate['projective_error'])-p.number(full['projective_error'])
                                except (KeyError, ValueError, TypeError):
                                    pass
                            admitted = raw is not None and left in eligible and right in eligible
                            pairs.append(dict(case=case, repetition=rep, question=arm+'_minus_full', cpu_budget=budget,
                                time=duration, budget_multiplier=cap, difference=raw, eligible=admitted,
                                outcome=('candidate_advantage' if raw < -AMBIGUITY else 'full_advantage' if raw > AMBIGUITY else 'unresolved')
                                    if admitted else 'incomplete',
                                candidate_error=None if candidate is None else candidate.get('projective_error'),
                                full_error=None if full is None else full.get('projective_error'),
                                candidate_floor=None if candidate is None else candidate.get('numerical_floor'),
                                full_floor=None if full is None else full.get('numerical_floor')))
    return pairs


def primary_decision(pairs, cases, *, preflight=False, complete=False):
    results = {}
    for question in ('window_minus_full', 'annealing_minus_full'):
        rows = [row for row in pairs if row['question'] == question and
                row['cpu_budget'] == 64 and row['time'] == 1.35 and row['budget_multiplier'] == 8]
        summaries = []
        for case in cases:
            within = [row for row in rows if row['case'] == case]
            values = [row['difference'] for row in within if row['difference'] is not None]
            summaries.append(dict(case=case, repetitions=within,
                median_difference=statistics.median(values) if len(values) == 2 else None))
        decision = 'incomplete'
        if preflight:
            decision = 'preflight_only' if complete else 'incomplete'
        elif complete and len(cases) == 2 and len(rows) == 4 and all(row['eligible'] for row in rows):
            medians = [summary['median_difference'] for summary in summaries]
            if all(value is not None and value < -AMBIGUITY for value in medians) and not any(row['difference'] > AMBIGUITY for row in rows):
                decision = 'candidate_advantage'
            elif all(value is not None and value > AMBIGUITY for value in medians) and not any(row['difference'] < -AMBIGUITY for row in rows):
                decision = 'full_advantage'
            else:
                decision = 'mixed_or_unresolved'
        results[question] = dict(decision=decision, case_summaries=summaries, ambiguity=AMBIGUITY,
            scope='Predeclared end-to-end comparison at B64/T1.35/8U; two molecules, paired timing/seed repetitions; no confidence interval')
    return results


def cost_totals(inventories, *, expected_count=None):
    actual = {key: 0. for key in p.COST_FIELDS}; attributed = dict(actual)
    unknown = {key: 0 for key in actual}; rows = []
    for inventory in inventories.values():
        for record in inventory['records'].values():
            arm = record['arm']; costs = record['costs']
            own = costs['prefix'] if arm == 'parent' else costs['continuation']
            # Actual parent appears once; each optimized comparison attributes it.
            charged = costs['total'] if arm in k.ARMS or arm in k.CHEAP else {key: 0. for key in actual}
            for key in actual:
                if own[key] is None:
                    actual[key] = None; unknown[key] += 1
                elif actual[key] is not None:
                    actual[key] += own[key]
                if charged[key] is None:
                    attributed[key] = None
                elif attributed[key] is not None:
                    attributed[key] += charged[key]
            rows.append(dict(case=record['case'], repetition=record['repetition'], arm=arm,
                cpu_budget=record['cpu_budget'], actual_own=own, comparison_attributed=charged, eligible=record['analysis_eligible']))
    inventory_complete = expected_count is None or len(rows) == expected_count
    known_actual, known_attributed = dict(actual), dict(attributed)
    if not inventory_complete:
        actual = {key: None for key in actual}; attributed = {key: None for key in attributed}
    return dict(actual=actual, comparison_attributed=attributed, inventory_complete=inventory_complete,
                observed_subtotal_actual=known_actual, observed_subtotal_attributed=known_attributed,
                unknown_cost_records=unknown, constructors=rows,
                scope='Constructor costs only; parent actual once per case/repetition; six attributed optimized prefixes; cheap own costs; physical grid never multiplies construction cost')


def analyze(data, *, preflight=False):
    report = dict(schema='molecular-analysis-v1', status='incomplete', complete=False, issues=[],
                  retained_physical_rows=[], retained_objective_rows=[], retained_raw_witnesses=[],
                  retained_reference_records=[], retained_independent_rows=[], cases_metadata=[], confidence_intervals=False)
    physical = {}; eligible = set(); inventories = {}; complete_cases = set()
    try:
        a, freeze, info, config, inputs = p.context(data, preflight=preflight)
        report.update(info); report['original_frozen_sources'] = freeze['sources']; report['current_consumer_sources'] = k.source_map()
        for case in info['cases']:
            reports = {}; inventory = panel = witnesses = None
            try:
                inventory, panel, witnesses = p.validate_case(a, freeze, info, config, inputs, case)
                inventories[case] = inventory
            except Exception as error:
                report['issues'].append(dict(case=case, phase='panel_consumer', type=type(error).__name__, message=str(error)))
            for phase, filename in (('physical', 'PHYSICAL.toml'), ('independent', 'VERIFY.toml'), ('objective', 'OBJECTIVE.toml')):
                try:
                    reports[phase] = a.toml(f'{phase}/{case}/{filename}')
                except Exception as error:
                    report['issues'].append(dict(case=case, phase=phase, type=type(error).__name__, message=str(error)))
            for field, phase, rows_key in (('retained_objective_rows', 'objective', 'logical_rows'),
                    ('retained_raw_witnesses', 'objective', 'raw_witnesses'),
                    ('retained_reference_records', 'physical', 'reference_records'),
                    ('retained_independent_rows', 'independent', 'rows')):
                report[field].extend(dict(case=case, **row) for row in reports.get(phase, {}).get(rows_key, []))
            for row in reports.get('physical', {}).get('rows', []):
                report['retained_physical_rows'].append(dict(case=case, **row))
                try:
                    key = (case,)+logical(row)
                    require(key not in physical, 'Duplicate retained physical row')
                    physical[key] = row
                except Exception as error:
                    report['issues'].append(dict(case=case, phase='physical_inventory', type=type(error).__name__, message=str(error)))
            if inventory is not None and len(reports) == 3:
                try:
                    joined = inspect_case(a, freeze, info, inventory, panel, witnesses, reports)
                    report['cases_metadata'].append(dict(case=case, complete=joined['complete'],
                        unique_circuits=joined['unique_circuits'], raw_claims_verified=joined['raw_claims_verified'], issues=joined['issues']))
                    if joined['complete']:
                        complete_cases.add(case)
                        eligible.update((case,)+identity for identity in joined['good'])
                except Exception as error:
                    report['issues'].append(dict(case=case, phase='artifact_joins', type=type(error).__name__, message=str(error)))
        try:
            report['scientific_phase_costs'] = phase_receipts(a, freeze, info)
        except Exception as error:
            report['issues'].append(dict(phase='scientific_receipts', type=type(error).__name__, message=str(error)))
        complete = not report['issues'] and complete_cases == set(info['cases'])
        pairs = paired_grid(info['cases'], info['repetitions'], physical, eligible)
        require(len(pairs) == (16 if preflight else 64), 'Complete declared named pair grid')
        report['named_pairs'] = pairs
        report['primary'] = primary_decision(pairs, info['cases'], preflight=preflight, complete=complete)
        report['cheap_context'] = [dict(case=case, **row) for (case, arm, *_), row in physical.items() if arm in k.CHEAP]
        report['constructor_costs'] = cost_totals(inventories, expected_count=len(info['cases'])*9*info['repetitions'])
        report['complete'] = complete; report['status'] = 'complete_preflight' if complete and preflight else 'complete' if complete else 'incomplete'
        report['artifact_hashes'] = a.hashes
    except Exception as error:
        report['issues'].append(dict(phase='frozen_context_or_inventory', type=type(error).__name__, message=str(error)))
    if 'primary' not in report:
        report['primary'] = primary_decision([], [], preflight=preflight)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args(argv)
    report = analyze(args.data, preflight=args.preflight)
    k.write_new(args.output, report)
    return 0 if report['complete'] else 2


if __name__ == '__main__':
    sys.exit(main())
