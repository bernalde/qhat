#!/usr/bin/env python3
"""Aggregate completed independent checks; never run a simulation or Lean build.

Run only after final analysis, all four campaign audits, and the final uniform
boundary sensitivity are complete. Missing/inconsistent evidence fails closed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def maximum(values):
    return max(values, default=0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root / 'analysis'
    inputs = {}

    def load(relative):
        path = root / relative
        require(path.is_file(), 'Required final evidence is missing: ' + relative)
        inputs[relative] = sha(path)
        return json.loads(path.read_text())

    def passed(record, name):
        require(record.get('status') in ('passed', 'pass', 'complete'), name + ' did not pass')

    def case_coverage(records, cases, name):
        names = [item['case'] for item in records]
        require(len(names) == len(cases) and set(names) == set(cases), name + ' needs exactly all 14 cases')

    analysis = load('analysis/summary.json')
    passed(analysis, 'Final analysis')
    require((analysis['rows'], analysis['completed_evolutions'], analysis['infeasible']) ==
            (4312, 4186, 126), 'Final analysis must contain 4312 rows, 4186 states and 126 infeasible budgets')
    cases = analysis['cases']
    require(len(cases) == 14 and len(set(cases)) == 14, 'Final analysis needs 14 distinct cases')
    raw = {}
    raw_paths = {}
    for relative, digest in analysis['input_sha256'].items():
        path = root / relative
        require(path.is_file() and sha(path) == digest, 'Analyzed raw row changed: ' + relative)
        row = json.loads(path.read_text())
        campaign = 'optimized' if relative.startswith('optimized/') else (
            'fixed' if relative.startswith('fixed/') else (
                'extra' if relative.startswith('published/extra/') else
                'v2' if relative.startswith('published/v2/') else None))
        require(campaign is not None, 'Unknown raw-row campaign: ' + relative)
        token = (campaign, row['case'], row['ordering'], int(row['budget_units']))
        require(token not in raw, 'Duplicate raw row: ' + str(token))
        require(row['case'] in cases and row['status'] in ('complete', 'infeasible'), 'Invalid raw row: ' + relative)
        raw[token] = row
        raw_paths[token] = relative
    require(len(raw) == analysis['rows'], 'Analysis raw-input manifest has incomplete coverage')
    require(sum(r['status'] == 'complete' for r in raw.values()) == 4186 and
            sum(r['status'] == 'infeasible' for r in raw.values()) == 126, 'Raw status counts disagree')

    campaigns = {}
    audited_tokens = set()
    for name in ('fixed', 'optimized', 'extra', 'v2'):
        summary = load('independent/' + name + '/summary.json')
        passed(summary, name + ' campaign summary')
        require(summary['case_count'] == 14 and len(summary['cases']) == 14 and
                set(summary['cases']) == set(cases), name + ' audit omits cases')
        check_records = []
        for case in cases:
            record = load('independent/' + name + '/' + case + '.json')
            passed(record, name + '/' + case)
            require(record['case'] == case and record['rows'] == len(record['checks']), 'Malformed case audit')
            for checked in record['checks']:
                token = (name, case, checked['ordering'], int(checked['budget_units']))
                require(token in raw and token not in audited_tokens, 'Missing/duplicate audited row: ' + str(token))
                audited_tokens.add(token)
                row = raw[token]
                require(row['extension_design_sha256'] == record['design_sha256'], 'Audit/row design digest differs')
                require(checked['row_sha256'] == analysis['input_sha256'][raw_paths[token]], 'Audit/analysis row digest differs')
                if row['status'] == 'infeasible':
                    require(checked['status'] == 'infeasible', 'Infeasible audit mismatch')
                    require(checked.get('state_distance') is None, 'Infeasible row was counted as replayed')
                else:
                    require(checked['status'] == 'passed' and checked['state_sha256'] == row['state_sha256'], 'State audit mismatch')
                    differences = checked['metric_difference']
                    require(differences and all(math.isfinite(v) and 0 <= v < 2e-12 for v in differences.values()), 'Invalid metric check')
                    distance = checked.get('state_distance')
                    require(distance is None or (math.isfinite(distance) and 0 <= distance < 2e-11), 'Invalid replay check')
                check_records.append(checked)
            require(record['replayed'] == sum(x.get('state_distance') is not None for x in record['checks']), 'Per-case replay count differs')
        count = sum(x.get('state_distance') is not None for x in check_records)
        distance = maximum(x['state_distance'] for x in check_records if x.get('state_distance') is not None)
        metric = maximum(max(x['metric_difference'].values()) for x in check_records if 'metric_difference' in x)
        require(summary['rows'] == len(check_records) and summary['replayed'] == count, 'Campaign coverage total differs')
        require(summary['max_state_distance'] == distance and summary['max_metric_difference'] == metric, 'Campaign maximum differs')
        campaigns[name] = {'cases': 14, 'audited_rows': len(check_records),
                           'saved_states_checked': sum(x['status'] == 'passed' for x in check_records),
                           'infeasible': sum(x['status'] == 'infeasible' for x in check_records),
                           'independently_replayed': count, 'maximum_state_distance': distance,
                           'maximum_metric_difference': metric, 'replay_selection': summary['replay_selection']}
    require(audited_tokens == set(raw), 'Not every analyzed row has an independent audit')

    references = load('independent/references/summary.json')
    passed(references, 'References')
    case_coverage(references['checks'], cases, 'Reference checks')
    require(references['cases'] == 14, 'Reference summary count differs')
    for item in references['checks']:
        passed(item, 'Reference ' + item['case'])
        require(0 <= item['reference_distance'] < 2e-11, 'Reference distance exceeds audit tolerance')
    reference_max = maximum(x['reference_distance'] for x in references['checks'])
    require(references['max_reference_distance'] == reference_max, 'Reference maximum differs')
    reference_limit = ('This is an aggregation of the saved independent reference checks, not a new reference solve. '
                       'Reference-checkpoint reuse verifies the checker and core source hashes but skips rereading '
                       'the original input bytes and regenerating or rehashing the reference NPZ. The prior-artifact '
                       'manifest check separately establishes original-artifact byte integrity. Campaign audit '
                       'checkpoint reuse rechecks audited raw rows and state hashes, but does not independently '
                       'revalidate a reused reference NPZ. These limitations are not hidden by the aggregate pass.')

    objectives = load('independent/objectives/summary.json')
    passed(objectives, 'Objectives')
    case_coverage(objectives['checks'], cases, 'Objective checks')
    objective_records = []
    for item in objectives['checks']:
        passed(item, 'Objective ' + item['case'])
        require(len(item['checks']) == 4 and len({x['ordering'] for x in item['checks']}) == 4, 'Need four objectives per case')
        objective_records.extend(item['checks'])
    require(objectives['cases'] == 14 and objectives['objectives'] == len(objective_records) == 56, 'Objective count differs')
    require(all(math.isfinite(x['absolute_difference']) and x['absolute_difference'] >= 0 for x in objective_records), 'Invalid objective difference')

    partitions = load('independent/partition_checks.json')
    passed(partitions, 'Partitions')
    case_coverage(partitions['checks'], cases, 'Partition checks')
    require(partitions['cases'] == 14, 'Partition count differs')
    for item in partitions['checks']:
        require(set(item['checks']) == {'parent', 'owned', 'class', 'refined'}, 'Partition family coverage differs')
        for details in item['checks'].values():
            require(all(details[k] is True for k in ('exact_target_sum', 'exact_number_commutator_zero',
                                                     'pairwise_internal_commutation')), 'A partition property failed')

    archived = load('optimized/archived_state_checks.json')
    require(archived['count'] == len(archived['checks']), 'Archived state count differs')
    require(archived['bitwise_equal_count'] == sum(x['bitwise_equal'] for x in archived['checks']), 'Bitwise count differs')
    require(archived['bitwise_equal_count'] == archived['count'] and archived['maximum_norm_difference'] == 0 and
            all(x['norm_difference'] == 0 for x in archived['checks']), 'An archived optimized state differs')
    historical = load('fixed/historical_internal_order_check.json')
    passed(historical, 'Historical internal-order checks')
    require(historical['unchanged_frozen_campaign'] is True and
            all(x['same_ordered_blocks_exact'] is True for x in historical['checks']), 'Historical equivalence failed')
    for field, entry in (('max_state_norm_distance', 'state_norm_distance'),
                         ('max_one_minus_overlap_absolute_difference', 'one_minus_overlap_absolute_difference'),
                         ('max_particle_number_leakage_absolute_difference', 'particle_number_leakage_absolute_difference')):
        require(historical[field] == maximum(x[entry] for x in historical['checks']), 'Historical maximum differs')
    historical_alignment = load('fixed/archived_schedule_check.json')
    literal_matches = sum(x['historical_flat_hash_identical'] for x in historical_alignment['checks'])
    serialized_matches = sum(x['serialized_design_schedules_identical'] for x in historical_alignment['checks'])
    require(serialized_matches == literal_matches + len(historical['checks']), 'Historical coverage differs')

    lean = load('fixed/lean_bridge_correspondence_audit.json')
    passed(lean, 'Lean manuscript correspondence')
    require(len(lean['covered_theorems']) == 19 and len(set(lean['covered_theorems'])) == 19 and not lean['issues'], 'Lean correspondence incomplete')
    prior = load('packaging/prior_artifact_check.json')
    passed(prior, 'Prior artifacts')
    require(all(prior[k] > 0 for k in ('qhat_prior_artifact_files', 'qhat_prior_embedded_files', 'paper_prior_files')), 'Prior integrity coverage absent')

    boundary = load('optimized/boundary_sensitivity.json')
    require(not boundary['pending_v2_plans'], 'Boundary sensitivity still has pending V2 plans')
    case_coverage(boundary['refined_label_correspondence'], cases, 'Boundary label correspondence')
    require(all(x['literal_flattened_refined_label_equals_complete_class_label'] is True
                for x in boundary['refined_label_correspondence']), 'Refined/complete label equivalence failed')
    symmetric = {token: row for token, row in raw.items() if row['status'] == 'complete' and row['formula_order'] == 2}
    boundary_tokens = set()
    for item in boundary['cost_rows']:
        name = 'v2' if item['campaign'] == 'published/v2' else item['campaign']
        token = (name, item['case'], item['ordering'], int(item['original_budget_units']))
        require(token in symmetric and token not in boundary_tokens, 'Boundary row absent/duplicated')
        boundary_tokens.add(token)
        row = symmetric[token]
        require(item['steps'] == row['steps'] and item['state_sha256'] == row['state_sha256'] and
                item['schedule_sha256'] == row['schedule_sha256'] and
                item['one_minus_overlap'] == row['one_minus_overlap'] and
                item['original_emitted_rotations'] == row['pauli_rotations'], 'Boundary changes the measured experiment')
        require(0 < item['uniform_endpoint_rotations'] <= item['original_emitted_rotations'], 'Invalid fused cost')
    require(boundary_tokens == set(symmetric), 'Boundary sensitivity omits saved symmetric states')
    plan_tokens = [('v2' if p['campaign'] == 'published/v2' else p['campaign'], p['case'], p['ordering'])
                   for p in boundary['plans']]
    require(len(set(plan_tokens)) == len(plan_tokens) and set(plan_tokens) == {t[:3] for t in symmetric}, 'Boundary plan coverage differs')
    dense = boundary['dense_checks']
    require(dense['dense_products'] == len(dense['records']) and dense['dense_products'] > 0 and
            dense['maximum_unitary_frobenius_difference'] == maximum(x['unitary_frobenius_difference'] for x in dense['records']) and
            dense['maximum_unitary_frobenius_difference'] < 3e-12, 'Boundary dense checks failed')

    result = {
        'status': 'passed', 'scope': 'Aggregation and cross-linking of completed evidence; no new simulation or Lean rebuild.',
        'cases': cases, 'rows': len(raw), 'saved_states_checked': sum(c['saved_states_checked'] for c in campaigns.values()),
        'infeasible': sum(c['infeasible'] for c in campaigns.values()),
        'independently_replayed': sum(c['independently_replayed'] for c in campaigns.values()),
        'maximum_replayed_state_distance': max(c['maximum_state_distance'] for c in campaigns.values()),
        'maximum_saved_metric_difference': max(c['maximum_metric_difference'] for c in campaigns.values()),
        'campaigns': campaigns,
        'references': {'cases': 14, 'maximum_reference_distance': reference_max,
                       'development_jw_magnitude_smoke_replays': sum(x.get('jw_magnitude_r100_state_distance') is not None for x in references['checks']),
                       'checkpoint_reuse_limit': reference_limit},
        'objectives': {'cases': 14, 'checked_objectives': len(objective_records),
                       'maximum_absolute_difference': maximum(x['absolute_difference'] for x in objective_records)},
        'partitions': {'cases': 14, 'families_per_case': 4,
                       'properties': ['exact target sum', 'zero exact particle-number commutator', 'pairwise internal commutation']},
        'archived_optimized_states': {'count': archived['count'], 'bitwise_equal_count': archived['bitwise_equal_count'], 'maximum_norm_difference': 0.0},
        'historical_fixed_schedules': {'literal_flat_sequence_matches': literal_matches,
                                       'commuting_equivalent_sequences': len(historical['checks']),
                                       'maximum_state_norm_distance': historical['max_state_norm_distance'],
                                       'maximum_error_difference': historical['max_one_minus_overlap_absolute_difference']},
        'lean_manuscript_correspondence': {'public_statements': len(lean['covered_theorems']), 'new_lean_rebuild': False, 'review_kind': lean['review_kind']},
        'prior_artifacts': prior,
        'boundary_sensitivity': {'plans': len(plan_tokens), 'existing_symmetric_states': len(boundary_tokens),
                                 'dense_toy_products': dense['dense_products'],
                                 'maximum_dense_unitary_difference': dense['maximum_unitary_frobenius_difference'],
                                 'new_molecular_propagations': 0,
                                 'scope': 'Uniform endpoint-fusion cost reinterpretation of unchanged saved states and step counts, not new fixed-budget error measurements.'},
        'input_sha256': inputs, 'source_sha256': sha(Path(__file__))}
    lines = [
        '# Validation coverage', '',
        f"Passed coverage gates for all {len(cases)} cases: {len(raw):,} planned rows, {result['saved_states_checked']:,} saved states and {result['infeasible']} infeasible budgets.", '',
        '| Campaign | Audited rows | Saved states checked | Independently replayed | Maximum replay distance |',
        '|---|---:|---:|---:|---:|']
    for name, item in campaigns.items():
        lines.append(f"| {name} | {item['audited_rows']:,} | {item['saved_states_checked']:,} | {item['independently_replayed']:,} | {item['maximum_state_distance']:.6g} |")
    lines += ['',
        f"All saved states received independent metric checks; {result['independently_replayed']:,} circuits were independently replayed according to the predeclared, error-independent selection. The maximum replayed state-vector distance is {result['maximum_replayed_state_distance']:.6g}; the maximum saved-metric difference is {result['maximum_saved_metric_difference']:.6g}. Replay coverage is a subset, not a claim to rerun every circuit. Counts are derived from the per-row check records.", '',
        f"Independent target reference checks cover all 14 cases (maximum distance {reference_max:.6g}). The {len(objective_records)} objective checks and four exact partition-family checks per case pass; the latter verify target sums, internal commutation and particle-number commutators. Maximum objective-score difference: {result['objectives']['maximum_absolute_difference']:.6g}.", '',
        f"All {archived['count']} checked archived optimized states agree bitwise. Historical fixed schedules comprise {literal_matches} literal flat-sequence matches and {len(historical['checks'])} exact commuting-block equivalences; the latter have maximum saved-state distance {historical['max_state_norm_distance']:.6g} and error difference {historical['max_one_minus_overlap_absolute_difference']:.6g}. Internal permutations are not described as literal sequence matches.", '',
        f"Uniform boundary sensitivity covers {len(plan_tokens)} symmetric plans and {len(boundary_tokens)} existing states, with {dense['dense_products']} dense toy-product checks (maximum difference {dense['maximum_unitary_frobenius_difference']:.6g}). It changes interpreted cost only, not a state, step count, or error, and makes no newly saturated-budget measurement claim.", '',
        f"The manuscript correspondence check covers all {len(lean['covered_theorems'])} public bridge statements. It is a source/proof-step review, with no new Lean compilation or axiom audit. Prior-artifact verification found {prior['qhat_prior_artifact_files']} QHAT files, {prior['qhat_prior_embedded_files']} embedded files and {prior['paper_prior_files']} earlier paper files unchanged; this is byte-integrity evidence, not a new scientific rerun.", '',
        '## Reference-checkpoint limitation', '', reference_limit, '',
        'This aggregator rehashes the analyzed raw JSON rows and cross-links their existing state digests and per-row audits. It does not reread all state arrays, redo matrix propagation, regenerate references, or rebuild Lean. Input evidence hashes and the exact coverage are retained in validation.json.']
    output.mkdir(parents=True, exist_ok=True)
    for name, content in (('validation.json', json.dumps(result, indent=2, allow_nan=False) + '\n'),
                          ('VALIDATION.md', '\n'.join(lines) + '\n')):
        path = output / name
        temporary = path.with_suffix(path.suffix + '.tmp')
        temporary.write_text(content)
        temporary.replace(path)
    print(json.dumps({k: result[k] for k in ('status', 'rows', 'saved_states_checked', 'infeasible', 'independently_replayed')}))


if __name__ == '__main__':
    main()
