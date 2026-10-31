"""Private exact known-censor delivery contract; never changes scientific gates.

Only stdlib metadata joins. Actual numerical work remains the frozen Julia CLI.
The expectation file binds CLOSED producer objective artifacts, not future data.
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
HF = 'HF_r1p00_631g_12m8e'
H2S = 'H2S_r1p40_a96_631g_12m8e'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    path = Path(path); size = path.stat().st_size
    require(0 <= size <= 8*1024**2, 'Bounded exact-censor metadata')
    with path.open('rb') as stream:
        raw = stream.read(size+1)
    require(len(raw) == size, 'Censor metadata changed while read')
    return json.loads(raw)


EXPECTED = read(HERE/'EXPECTED_OUTCOMES.json')


def comparison_verified(value):
    require(value['verified'] is True and value['within_tolerance'] is True and
            value['numerical_floor'] is False, 'Existing objective comparison must remain verified above floor')


def audit_objective(case, report):
    expected = EXPECTED['reports'][case]
    require(all(report[key] == value for key,value in expected['header'].items()),
            'Exact original objective model/panel/witness/source/header bytes')
    require(report['status'] == expected['expected_status'] and report['complete'] is expected['expected_complete']
            and report['all_scored_objectives_verified'] is expected['all_scored_objectives_verified'],
            'Scientific objective completeness/censor must remain unchanged')
    require(not report['failures'] and report['mismatched_objectives'] == report['floor_censored_objectives'] == 0,
            'No additional objective mismatch/floor/failure admitted')
    require(report['witness_inventory_complete'] is True and report['processed_logical_rows'] ==
            report['logical_rows_retained'] == len(report['logical_rows']) == 72, 'All72 selected objective rows retained')
    require(report['expected_raw_witnesses'] == report['processed_raw_witnesses'] ==
            report['raw_witnesses_retained'] == len(report['raw_witnesses']) == expected['raw_count'],
            'Every exact persisted raw proposed claim retained')
    require(report['unique_objectives'] == len(report['checks']) == expected['unique_checks'] and
            all(row['status'] == 'replayed' for row in report['checks']), 'All exact objective evaluations replayed')
    require(report['native_refused_objectives'] == expected['native_refused_objectives'], 'Exact original native-refusal count')
    for row in report['logical_rows']:
        require(row['status'] == 'verified', 'All selected native and archived scores verified')
        comparison_verified(row['proposed_comparison']); comparison_verified(row['native_comparison'])
    require([row['name'] for row in report['raw_witnesses']] == expected['raw_names'], 'Exact ordered raw witness names')
    for row in report['raw_witnesses']:
        require(row['status'] == expected['raw_statuses'][row['name']], 'No added, cleared or moved raw censor')
        comparison_verified(row['proposed_comparison'])
        if row['status'] == 'native_refusal_censored':
            original = expected['known_censored_row']
            keys = set(original)-{'proposed_comparison'}
            require(all(row[key] == original[key] for key in keys), 'Exact original censored claim identity/order/refusal')
            require(row['name'] == 'rep2.window.B64.raw1' and case == HF and row['selected'] is False
                    and 'native_score' not in row and 'native_comparison' not in row,
                    'Missing original-call native score remains absent and unpromoted')
            refusal = json.loads(row['native_refusal'])
            require(refusal == dict(type='SearchDeadline', message='whole-call window allowance exhausted'),
                    'Only exact original window SearchDeadline censor admitted')
        else:
            require(row['status'] == 'verified' and 'native_score' in row, 'All other raw native claims verified')
            comparison_verified(row['native_comparison'])
    return dict(case=case, selected_verified=len(report['logical_rows']),
        raw_proposed_comparisons_verified=True, raw_claims=len(report['raw_witnesses']),
        unique_objective_evaluations=len(report['checks']),
        scientific_gate_complete=report['complete'], raw_claims_verified=case != HF,
        native_refused_raw_names=[row['name'] for row in report['raw_witnesses'] if row['status']=='native_refusal_censored'])


def audit_analysis(report):
    """Expected future source output, checked only after root actually runs it."""
    require(report['schema'] == 'molecular-analysis-v1' and report['mode'] == 'main' and
            report['source_freeze_sha256'] == EXPECTED['source_freeze_sha256'], 'Original main source analysis identity')
    require(report['complete'] is False and report['status'] == 'incomplete', 'Source analysis scientific gate MUST remain incomplete')
    require(set(report['primary']) == {'window_minus_full','annealing_minus_full'} and
            all(value['decision']=='incomplete' for value in report['primary'].values()), 'Both primary decisions MUST remain incomplete')
    require(report['issues'] == [EXPECTED['future_analysis']['expected_only_issue']],
            'Retain exact original HF scientific_receipts issue; no other issue admitted')
    require('scientific_phase_costs' not in report, 'Do not fabricate a successful original phase-cost gate')
    cases = report['cases_metadata']
    require([row['case'] for row in cases] == [HF,H2S], 'Both source case joins retained')
    for row in cases:
        require(not row['issues'] and row['complete'] is (row['case'] != HF) and
                row['raw_claims_verified'] is (row['case'] != HF), 'Only known HF raw native censor prevents completeness')
    for key in ('retained_physical_rows','retained_independent_rows','retained_objective_rows'):
        rows = report[key]
        require(len(rows)==144 and all(row['status']==('verified' if key=='retained_objective_rows' else 'complete')
                                      for row in rows), 'Full144 selected/physical/independent row grid retained')
    require(len(report['retained_reference_records'])==4 and
            all(row['status']=='complete' for row in report['retained_reference_records']), 'Both references per case retained')
    require(len(report['retained_raw_witnesses'])==14, 'Every actual raw proposed claim retained in source analysis')
    for case in (HF,H2S):
        rows=[row for row in report['retained_raw_witnesses'] if row['case']==case]
        require([row['name'] for row in rows]==EXPECTED['reports'][case]['raw_names'], 'Source raw inventory exact')
        for row in rows:
            require(row['status']==EXPECTED['reports'][case]['raw_statuses'][row['name']], 'Source native censor preserved')
            comparison_verified(row['proposed_comparison'])
    require(len(report['named_pairs'])==64 and report['constructor_costs']['inventory_complete'] is True and
            len(report['constructor_costs']['constructors'])==36, 'Full pair/cost inventories retained without changing science gate')


def audit_inspection(case, joined):
    require(not joined['issues'] and len(joined['good'])==72, 'All72 selected/state/reference row joins must pass per case')
    require(joined['complete'] is (case != HF) and joined['raw_claims_verified'] is (case != HF),
            'Unchanged inspect_case must report HF incomplete/raw-unverified and H2S complete')


def original_objective_bindings(data):
    """Authenticate actual retained original nonzero status, not a broad waiver."""
    for case, expected in EXPECTED['reports'].items():
        value=expected['producer_report']; path=data/value['path']
        require(sha(path)==value['sha256'], 'Original closed objective report bytes changed')
        audit_objective_file(case,path)
        for value in (expected['owner'][key] for key in ('plan','final','receipt')):
            require(sha(data/value['path'])==value['sha256'], 'Original objective owning receipt/PLAN/FINAL bytes changed')


def audit_objective_file(case,path):
    """Bounded post-teardown outcome check; no scientific score is recomputed."""
    import tomllib
    path=Path(path); size=path.stat().st_size
    require(0<=size<=8*1024**2,'Bounded objective report metadata')
    with path.open('rb') as stream:
        raw=stream.read(size+1)
    require(len(raw)==size,'Objective report changed during metadata read')
    return audit_objective(case,tomllib.loads(raw.decode()))


def analysis_job_binding(data, gate, freeze, analysis):
    value=gate['expected_censored_analysis_job']
    require(set(value)=={'plan','final','receipt'}, 'Root must late-bind exact actual censored analysis job')
    for item in value.values():
        require(isinstance(item['path'],str) and not Path(item['path']).is_absolute() and
                '..' not in Path(item['path']).parts and sha(data/item['path'])==item['sha256'],
                'Exact root-bound analysis job artifact')
    plan=read(data/value['plan']['path']); label=plan['label']
    require(value['plan']['path']=='jobs/'+label+'.PLAN.json' and
            value['final']['path']=='jobs/'+label+'.FINAL.json' and
            value['receipt']['path']=='logs/'+label+'.resources.json', 'Canonical exact analysis owner')
    require(plan['phase']=='analysis' and plan['source_freeze_sha256']==EXPECTED['source_freeze_sha256'] and
            plan['command']==[freeze['runtime']['executable'],freeze['producer_root']+
            '/numerical_audit/molecular_generalization/analyze.py',freeze['producer_data_root'],str(analysis)],
            'Only actual unchanged main analyzer CLI may have the expected scientific exit2')
    return value['plan']['path']
