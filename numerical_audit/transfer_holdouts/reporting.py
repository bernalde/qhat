"""Complete transfer tables, costs and losses from independently verified data."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import verify as audit

HERE = Path(__file__).resolve().parent
audit.require(Path(audit.__file__).resolve() == HERE / 'verify.py', 'Wrong report audit binding')
POOLS = {
    'chemistry.tests': (120, 900, 64), 'chemistry.preflight': (120, 900, 256),
    'chemistry.generation': (240, 1800, 256), 'reference.preflight': (120, 900, 256),
    'reference.main': (120, 900, 256), 'full.full_preflight': (600, 3600, 256),
    'full.measurement': (2400, 14400, 512), 'full.independent': (2400, 14400, 1024),
    'full.reporting': (300, 1800, 1024),
}


def resources(path):
    path = Path(path)
    audit.require(path.stat().st_size <= 32 * 1024**2, 'Resource ledger byte ceiling')
    ledger = json.loads(path.read_text())
    audit.require(set(ledger['pools']) == set(POOLS) and ledger['snapshot_scope'], 'Complete resource pool inventory required')
    used, seen = {}, set()
    for name, limits in POOLS.items():
        pool = ledger['pools'][name]
        audit.require(tuple(pool[field] for field in ('CPU_seconds', 'wall_seconds', 'disk_MiB')) == limits,
                      'Original transfer phase ceiling changed')
        cpu = wall = retained = 0.
        for attempt in pool['attempts']:
            filename = Path(attempt['receipt']).name
            audit.require(filename not in seen, 'A resource attempt is charged to more than one pool')
            seen.add(filename)
            original = ledger['receipts'][filename]
            audit.require(hashlib.sha256(original['text'].encode()).hexdigest() == original['sha256'], 'Raw resource receipt changed')
            receipt = json.loads(original['text'])
            actual_cpu, elapsed = receipt['final_accounting']['CPU_seconds'], receipt['elapsed']
            audit.require(math.isfinite(attempt['CPU_seconds']) and attempt['CPU_seconds'] + 1e-9 >= actual_cpu and
                          abs(attempt['wall_seconds'] - elapsed) <= 1e-8 and attempt['exit_code'] == receipt['exit_code'],
                          'Resource attempt undercharged or mislabeled')
            settings = receipt['settings']
            audit.require(all(settings[key] == value for key, value in
                          dict(CPUQuota='25%', Nice='19', CPUWeight='5', IOWeight='10', MemoryHigh='448M',
                               MemoryMax='512M', MemorySwapMax='0').items()), 'Numerical worker resource profile differs')
            audit.require(all(receipt['environment'][key] == '1' for key in
                          ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS')),
                          'Thread policy differs')
            cpu += attempt['CPU_seconds']
            wall += attempt['wall_seconds']
            retained += receipt['final_accounting']['output_bytes']
        audit.require(cpu <= limits[0] + 1e-8 and wall <= limits[1] + 1e-8 and
                      retained <= limits[2] * 1024**2, 'Cumulative completed phase usage exceeds its original ceiling')
        audit.require(abs(cpu - pool['consumed_CPU_seconds']) <= 1e-8 and
                      abs(wall - pool['consumed_wall_seconds']) <= 1e-8, 'Resource pool total differs from raw attempts')
        used[name] = dict(CPU_seconds=cpu, wall_seconds=wall, retained_output_bytes=retained,
                          attempts=len(pool['attempts']), failed_attempts=sum(row['exit_code'] != 0 for row in pool['attempts']))
    audit.require(seen == set(ledger['receipts']), 'Uncharged or omitted raw resource receipt')
    return ledger, used


def csv_file(path, fields, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def require_verified(output, verification):
    verified = audit.read(verification / 'VERIFIED.json')
    audit.require(verified['status'] == 'passed', 'Complete independent audit required before reporting')
    audit.require(verified['schedule_freeze_sha256'] == audit.sha(output / 'SCHEDULES_FROZEN.json') and
                  verified['independent_errors_sha256'] == audit.sha(verification / 'INDEPENDENT_ERRORS.json') and
                  verified['independent_comparisons_sha256'] == audit.sha(verification / 'INDEPENDENT_COMPARISONS.json'),
                  'Independent audit record binding differs')
    frozen = audit.read(output / 'SCHEDULES_FROZEN.json')
    for name, digest in {**frozen['artifacts'], **verified['physical_artifacts']}.items():
        audit.require(audit.sha(output / audit.safe(name)) == digest, 'Audited scientific bytes changed before reporting')
    audit.require(audit.inventory(output) == set(frozen['artifacts']) | set(verified['physical_artifacts']) | {'SCHEDULES_FROZEN.json'},
                  'Audited payload inventory differs')
    return verified


def report(output, verification, ledger_path, destination):
    output, verification, destination = Path(output), Path(verification), Path(destination)
    verified = require_verified(output, verification)
    audit.require(not destination.exists(), 'Fresh report directory required;retain prior reports')
    ledger, usage = resources(ledger_path)
    plan = audit.read(output / 'PLAN_FROZEN.json')
    protocol = audit.read(output / 'PROTOCOL_FROZEN.json')
    audit.require((verified['logical_circuits'], verified['unique_circuits'], verified['comparisons']) ==
                  (plan['expected_logical_circuits'], plan['expected_unique_circuits'], plan['expected_comparisons']),
                  'Full report cohort differs')
    destination.mkdir(parents=True)
    registry = audit.read(HERE / 'REGISTRY.json')
    families = {row['case']: row['family'] for row in registry['cases'] + [registry['chemistry_preflight']]}
    reference = {case: audit.read(output / 'reference_evidence/references' / (case + '.json')) for case in plan['cases']}
    constructor_rows, phase_rows = [], []
    def phases(name, method, row, scope):
        for phase, value in row['phases'].items():
            phase_rows.append(dict(constructor=name, method=method, scope=scope, phase=phase, CPU_seconds=value,
                                   included_in_parent=scope not in ('measured_outer', 'campaign_warmup')))
        if 'current_placement' in row:
            placed = row['current_placement']
            for phase, field in (('legacy_setup', 'setup_cpu_seconds'), ('legacy_placement', 'search_cpu_seconds')):
                phase_rows.append(dict(constructor=name, method=method, scope='included_legacy_initializer',
                                       phase=phase, CPU_seconds=placed[field], included_in_parent=True))
        if 'warm_start' in row:
            phases(name, 'first_JD', row['warm_start'], 'nested_first_order_warm')
    for specification in plan['training']:
        row = audit.read(output / 'constructors' / (specification['constructor'] + '.json'))
        search = row.get('result', {})
        placed = row.get('current_placement', row.get('warm_start', {}).get('current_placement'))
        constructor_rows.append(dict(specification, family=families[specification['case']],
            ordering_CPU_seconds=row['caller_CPU_seconds'], reference_CPU_charge=row['reference_preparation_CPU_charge'],
            standalone_reference_and_ordering_CPU_charge=row['standalone_reference_and_ordering_CPU_charge'],
            attempted=search.get('attempted'), completed=search.get('completed'),
            accepted_rescores=search.get('accepted_rescores'), nominal_search_CPU_seconds=search.get('nominal_search_CPU_seconds'),
            includes_fresh_first_order_warm=specification['method'] == 'cubic_JD',
            legacy_evaluations=placed['evaluations'], legacy_setup_CPU_seconds=placed['setup_cpu_seconds'],
            legacy_search_CPU_seconds=placed['search_cpu_seconds']))
        phases(specification['constructor'], specification['method'], row, 'measured_outer')
    for case in plan['cases']:
        for method, row in audit.read(output / 'warmups' / (case + '.json')).items():
            phases(case + '.HF.' + method, method, row, 'campaign_warmup')
    csv_file(destination / 'constructors.csv', ['constructor', 'case', 'family', 'state', 'method', 'training_time',
             'training_rotation_cap_units', 'ordering_CPU_seconds', 'reference_CPU_charge',
             'standalone_reference_and_ordering_CPU_charge', 'attempted', 'completed', 'accepted_rescores',
             'nominal_search_CPU_seconds', 'includes_fresh_first_order_warm', 'legacy_evaluations',
             'legacy_setup_CPU_seconds', 'legacy_search_CPU_seconds'], constructor_rows)
    csv_file(destination / 'constructor_phases.csv', ['constructor', 'method', 'scope', 'phase', 'CPU_seconds', 'included_in_parent'], phase_rows)
    uses = Counter(row['candidate'] for row in plan['logical'])
    unique = {}
    fields = ['candidate', 'case', 'family', 'state', 'correlation_label', 'one_RDM_idempotency_defect',
              'ground_excited_gap', 'near_degenerate', 'method', 'training_state', 'formula', 'T', 'cap_multiple',
              'steps', 'status', 'reason', 'ideal_rotations', 'expected_ideal_rotations',
              'two_qubit_gates', 'unit_duration_ASAP_depth', 'CX_dependency_depth', 'gates',
              'partial_two_qubit_gates', 'partial_unit_duration_ASAP_depth', 'partial_gates',
              'error', 'infidelity', 'norm', 'number_leakage', 'joint_sector_leakage', 'at_numerical_floor',
              'input_generation_CPU_charge', 'standalone_reference_CPU_charge', 'ordering_CPU_charge',
              'outer_compilation_CPU_seconds', 'standalone_without_input_CPU_charge', 'standalone_including_input_CPU_charge',
              'quantum_preparation_circuit_cost', 'compiled_resource_scope', 'logical_uses_of_candidate']
    for specification in plan['circuits']:
        row = audit.read(output / 'physical' / (specification['candidate'] + '.json'))
        metadata = reference[row['case']]['states'][row['state']]['metadata']
        result = dict(row, family=families[row['case']], correlation_label=metadata['correlation_label'],
                      one_RDM_idempotency_defect=metadata['one_RDM_idempotency_defect'],
                      ground_excited_gap=metadata.get('ground_excited_gap'), near_degenerate=metadata.get('near_degenerate'),
                      logical_uses_of_candidate=uses[row['candidate']])
        if row['status'] == 'complete':
            result.update(row['compiled_resources'])
            result['error'] = row['error']
        else:
            result.update({'partial_' + key: value for key, value in row['partial_compiled_resources'].items()})
        unique[row['candidate']] = result
    csv_file(destination / 'unique_circuits.csv', fields, unique.values())
    logical = [dict(unique[row['candidate']], **{key: value for key, value in row.items() if key in ('arm',)}) for row in plan['logical']]
    csv_file(destination / 'logical_observations.csv', ['arm'] + fields, logical)
    pairs = audit.read(output / 'COMPARISONS.json')
    for row in pairs:
        row['family'] = families[row['case']]
        row['correlation_label'] = reference[row['case']]['states'][row['state']]['metadata']['correlation_label']
    csv_file(destination / 'comparisons.csv', ['case', 'family', 'state', 'correlation_label', 'T', 'cap_multiple', 'a', 'b',
             'a_status', 'b_status', 'a_error', 'b_error', 'a_minus_b', 'decision', 'informative', 'both_at_floor'], pairs)
    groups = defaultdict(Counter)
    for row in pairs:
        groups[row['family'], row['state'], row['T'], row['cap_multiple'], row['a'], row['b']][row['decision']] += 1
    summaries = [dict(family=key[0], state=key[1], T=key[2], cap_multiple=key[3], a=key[4], b=key[5],
                      observations=sum(value.values()), a_wins=value['a'], b_wins=value['b'],
                      numerical_ties=value['numerical_tie'], floor_ties=value['floor_tie'], unavailable=value['unavailable'])
                 for key, value in sorted(groups.items())]
    csv_file(destination / 'grouped_pairs.csv', ['family', 'state', 'T', 'cap_multiple', 'a', 'b', 'observations',
             'a_wins', 'b_wins', 'numerical_ties', 'floor_ties', 'unavailable'], summaries)
    shutil.copyfile(ledger_path, destination / 'EXECUTION_RESOURCES_SNAPSHOT.json')
    summary = dict(status='complete_independently_audited_tables', preflight=protocol['preflight'],
        measured_outer_constructors=len(constructor_rows), unique_circuits=len(unique), logical_circuits=len(logical),
        complete=verified['complete'], censored=verified['censored'], comparisons=len(pairs),
        pair_decisions=dict(Counter(row['decision'] for row in pairs)),
        actual_preparation=audit.read(output / 'PREPARE_COSTS.json'), actual_screening=audit.read(output / 'RUN_COSTS.json'),
        independent_replay_CPU_seconds=verified['CPU_seconds'], completed_phase_usage=usage,
        resource_snapshot_scope=ledger['snapshot_scope'], producer_source_commit=protocol['source_commit'],
        schedule_freeze_sha256=verified['schedule_freeze_sha256'], verification_sha256=audit.sha(verification / 'VERIFIED.json'),
        reporter_source_sha256=audit.sha(Path(__file__)),
        limits=['Finite active-space development/transfer panel;no large-molecule or hardware scaling claim.',
                'Logical observations share molecules,states,times and constructors;counts are descriptive.',
                'Nominal search allowances differ;this is not a matched-total construction-cap experiment.',
                'Non-HF resources cover evolution only;quantum preparation circuit cost is unknown.',
                'Floating independent replay is not a rigorous interval certificate or a Lean molecular theorem.'])
    audit.write(destination / 'SUMMARY.json', summary)
    text = ('# Broader molecular, state and time transfer\n\n' +
            ('Representative LiH1.841 preflight; excluded from the main scientific panel.\n\n' if protocol['preflight'] else '') +
            f"All {len(logical)} logical observations, {len(unique)} unique circuits and {len(pairs)} arm pairs are retained. " +
            f"Complete/censored unique circuits: {verified['complete']}/{verified['censored']}.\n\n" +
            'constructors.csv retains each measured outer constructor. constructor_phases.csv also shows campaign warmups and nested first-order warm work; nested and legacy-initializer rows are included in their parent and must not be added again.\n\n' +
            'unique_circuits.csv records actual unique work. logical_observations.csv preserves every declared arm, with exact sharing counts. comparisons.csv and grouped_pairs.csv retain every loss, tie, floor and unavailable result. These dependent observations are descriptive counts.\n\n' +
            'Standalone charges include the training and evaluation reference once each when distinct, or once when shared. Model generation is separate and also shown in the from-scratch total. Actual campaign costs do not multiply shared reference/ordering charges by logical arms. Offline full-H screening and independent replay remain separate.\n\n' +
            'HF resources include its preparation prefix. Non-HF counts cover evolution only, with quantum preparation cost unknown. The search controls are nominal allowances, not matched-total construction budgets. All measurements use the frozen prototype exact-angle line policy, not hardware execution.\n\n' +
            'Resource snapshots contain every completed attempt, including failures, through the stated boundary; the current reporting worker receives its final external receipt after exit. Exact producer/verifier/source bindings remain in the complete delivery.\n')
    (destination / 'README.md').write_text(text)
    manifest = {path.name: audit.sha(path) for path in sorted(destination.iterdir()) if path.is_file()}
    audit.write(destination / 'REPORT_MANIFEST.json', manifest)
    return summary
