"""Exact strictness,zero-denominator and premature-reporting checks."""
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import pytest
from reporting import winner,sufficient_epsilon_domain,authenticated_inputs


def test_ranking_keeps_shared_endpoint_unresolved():
    assert winner((F(1,7),F(2,7)),(F(3,7),F(4,7))) == 'current'
    assert winner((F(3,7),F(4,7)),(F(1,7),F(2,7))) == 'diagonal_prefix'
    assert winner((F(1,7),F(2,7)),(F(2,7),F(3,7))) == 'unresolved'


def test_sufficient_coupling_domain_keeps_strict_threshold_and_zero_case():
    result = sufficient_epsilon_domain((F(1),F(2)),(F(3),F(4)),F(3),F(5))
    assert result['epsilon_upper'] == '1/8' and result['upper_inclusive'] is False
    assert sufficient_epsilon_domain((1,2),(3,4),3,5)['epsilon_upper'] == '1/8'
    with pytest.raises(ValueError,match='Exact rational'):
        sufficient_epsilon_domain((1.,2.),(3.,4.),3.,5.)
    zero = sufficient_epsilon_domain((F(1),F(2)),(F(3),F(4)),F(0),F(0))
    assert zero['epsilon_upper'] == '1' and zero['upper_inclusive'] is True
    tied = sufficient_epsilon_domain((F(1),F(3)),(F(3),F(4)),F(0),F(0))
    assert tied['status'] == 'no_certified_slope_gap'
    with pytest.raises(ValueError,match='Nonnegative'):
        sufficient_epsilon_domain((F(1),F(2)),(F(3),F(4)),F(-1),F(0))


def test_report_cannot_authenticate_missing_independent_audit(tmp_path):
    import json
    source,delivery = tmp_path/'producer',tmp_path/'delivery'
    source.mkdir();delivery.mkdir()
    for filename in ('PRODUCER_COMPLETE.json','FROZEN.json'):
        (source/filename).write_text('{}')
    with pytest.raises(FileNotFoundError):
        authenticated_inputs(source,delivery)
    (delivery/'INDEPENDENT_VERIFIED.json').write_text(json.dumps(dict(status='running')))
    (delivery/'DELIVERY_VERIFIED.json').write_text(json.dumps(dict(status='passed')))
    with pytest.raises(ValueError,match='unfinished'):
        authenticated_inputs(source,delivery)


def test_actual_merge_is_default_and_development_requires_explicit_flag(tmp_path):
    from reporting import dependency_record
    with pytest.raises(ValueError,match='explicit flag'):
        dependency_record(None,'abc')
    assert dependency_record(None,'abc',allow_development=True)['model_sha256'] == 'abc'
    invalid = tmp_path/'dependency.json'
    invalid.write_text('{"status":"development","repository":"bernalde/qhat"}')
    with pytest.raises(ValueError,match='Actual merged'):
        dependency_record(invalid,'abc')


def test_report_comparison_cannot_hide_modified_table(tmp_path):
    import hashlib,json
    from reporting import compare_reports
    a,b = tmp_path/'a',tmp_path/'b'
    for folder in (a,b):
        folder.mkdir()
        (folder/'table.csv').write_text('error\n1/7\n')
        digest = hashlib.sha256((folder/'table.csv').read_bytes()).hexdigest()
        (folder/'REPORT_MANIFEST.json').write_text(json.dumps({'table.csv':digest}))
    (b/'table.csv').write_text('error\n0\n')
    with pytest.raises(ValueError,match='Report bytes changed'):
        compare_reports(a,b)


def resource_fixture():
    import hashlib,json
    limits = dict(tests=(120,900,64),preflight=(240,1800,256),main=(1200,7200,256),
                  independent=(1200,7200,256),reporting=(300,1800,256))
    phases = {name:dict(CPU_seconds=c,wall_seconds=w,disk_MiB=d,attempts=[],
                       consumed_CPU_seconds=0,consumed_wall_seconds=0) for name,(c,w,d) in limits.items()}
    raw = dict(exit_code=1,settings=dict(CPUQuota='25%',Nice='19',CPUWeight='5',IOWeight='10',
               MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0'),
               environment={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                             'NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS')},
               budget=dict(cpu_seconds=120,wall_seconds=900,disk_mib=64),
               final_accounting=dict(CPU_seconds=.5,elapsed_seconds=2.,output_bytes=17))
    text = json.dumps(raw,indent=2)+'\n'
    phases['tests'].update(attempts=[dict(receipt='/owned/test.resources.json',exit_code=1,
                                         CPU_seconds=.5,wall_seconds=2.)],
                           consumed_CPU_seconds=.5,consumed_wall_seconds=2.)
    return dict(phases=phases,receipts={'test.resources.json':dict(sha256=hashlib.sha256(text.encode()).hexdigest(),text=text)})


def test_resource_ledger_keeps_failed_work_and_binds_raw_costs():
    from reporting import validate_resources
    record = resource_fixture()
    validate_resources(record)
    record['phases']['tests']['attempts'][0]['CPU_seconds'] = .1
    with pytest.raises(ValueError,match='raw receipt'):
        validate_resources(record)


def test_resource_ledger_refuses_changed_ceiling_and_tampered_receipt():
    from reporting import validate_resources
    record = resource_fixture()
    record['phases']['tests']['CPU_seconds'] = 121
    with pytest.raises(ValueError,match='ceiling changed'):
        validate_resources(record)
    record = resource_fixture()
    record['receipts']['test.resources.json']['text'] += ' '
    with pytest.raises(ValueError,match='hash differs'):
        validate_resources(record)


def report_pair_fixture(tmp_path):
    import hashlib,json
    roots = tmp_path/'kept',tmp_path/'replayed'
    for folder,head in zip(roots,('a'*40,'b'*40)):
        folder.mkdir()
        for name in ('all_finite_intervals.csv','finite_comparisons.csv','response_guards.csv'):
            (folder/name).write_text('error\n1/7\n')
        (folder/'README.md').write_text('Same fixed panel and conclusions.\n')
        (folder/'EXECUTION_RESOURCES_BEFORE_REPORT.json').write_text('{}\n')
        summary = dict(circuits=132,input_dependency=dict(merge_sha='c'*40),
            source_provenance=dict(producer='d'*40,verifier=head,reporter=head,
                reporter_sha256='e'*64,verifier_sha256={
                    'numerical_audit/molecular_guard/delivery.py':'f'*64,
                    'numerical_audit/molecular_guard/independent.py':'0'*64}))
        (folder/'SUMMARY.json').write_text(json.dumps(summary))
        manifest = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir()}
        (folder/'REPORT_MANIFEST.json').write_text(json.dumps(manifest))
    return roots


def replace_summary(folder, mutate):
    import hashlib,json
    path = folder/'SUMMARY.json'; row = json.loads(path.read_text()); mutate(row)
    path.write_text(json.dumps(row))
    manifest = folder/'REPORT_MANIFEST.json'; value = json.loads(manifest.read_text())
    value['SUMMARY.json'] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.write_text(json.dumps(value))


def test_replay_accepts_new_execution_commits_with_identical_verified_sources(tmp_path):
    from reporting import compare_reports
    kept,replayed = report_pair_fixture(tmp_path)
    compare_reports(kept,replayed)


def test_replay_rejects_different_verifier_sources_despite_same_tables(tmp_path):
    from reporting import compare_reports
    kept,replayed = report_pair_fixture(tmp_path)
    replace_summary(replayed,lambda row:row['source_provenance']['verifier_sha256'].update(
        {'numerical_audit/molecular_guard/independent.py':'1'*64}))
    with pytest.raises(ValueError,match='provenance report replay differs'):
        compare_reports(kept,replayed)


def test_replay_rejects_unbound_verifier_sources(tmp_path):
    from reporting import compare_reports
    kept,replayed = report_pair_fixture(tmp_path)
    replace_summary(replayed,lambda row:row['source_provenance'].pop('verifier_sha256'))
    with pytest.raises(ValueError,match='Verifier source binding'):
        compare_reports(kept,replayed)


def test_replay_keeps_actual_input_merge_exact(tmp_path):
    from reporting import compare_reports
    kept,replayed = report_pair_fixture(tmp_path)
    replace_summary(replayed,lambda row:row['input_dependency'].update(merge_sha='9'*40))
    with pytest.raises(ValueError,match='provenance report replay differs'):
        compare_reports(kept,replayed)


def test_report_verifier_binding_rejects_forged_source_hash(tmp_path):
    import hashlib,json,subprocess
    from reporting import trusted_verifier_sources,ROOT
    # This is a source-binding fixture, not an attestation that an audit ran.
    # It follows the committed candidate when a verifier guard is repaired.
    head=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    names=('numerical_audit/molecular_guard/delivery.py',
           'numerical_audit/molecular_guard/independent.py')
    expected=dict(head=head,sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                   for name in names},scope='Test source-binding fixture only')
    (tmp_path/'TRUSTED_VERIFIER_SOURCES.json').write_text(json.dumps(expected))
    (tmp_path/'DELIVERY_VERIFIED.json').write_text(json.dumps(dict(trusted_checker_head=head)))
    assert trusted_verifier_sources(tmp_path) == expected['sha256']
    expected['sha256']['numerical_audit/molecular_guard/independent.py'] = '0'*64
    (tmp_path/'TRUSTED_VERIFIER_SOURCES.json').write_text(json.dumps(expected))
    with pytest.raises(ValueError,match='Verifier source binding'):
        trusted_verifier_sources(tmp_path)


def approved_retry_extension_fixture():
    # In-memory policy fixture, never a persisted human approval.
    import hashlib,json
    record=resource_fixture()
    approval=dict(approved=True,scope='issue18_fresh_complete_independent_replay_after_expiry',
        original_phase_CPU_seconds=1200,new_phase_CPU_seconds=2000,additional_CPU_seconds=800,
        original_phase_wall_seconds=7200,new_phase_wall_seconds=10200,additional_wall_seconds=3000,
        attempt_CPU_seconds=900,attempt_wall_seconds=4200,
        MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0',CPUQuota='25%',threads=1,disk_MiB=256)
    text=json.dumps(approval)
    record['approved_ceiling_extensions']=[dict(phase='independent',additional_CPU_seconds=800,
        additional_wall_seconds=3000,approval_record=text,
        approval_sha256=hashlib.sha256(text.encode()).hexdigest())]
    record['phases']['independent'].update(CPU_seconds=2000,wall_seconds=10200)
    return record

def test_extended_replay_cap_requires_bound_explicit_approval():
    from reporting import validate_resources
    record=resource_fixture();record['phases']['independent'].update(CPU_seconds=2000,wall_seconds=10200)
    with pytest.raises(ValueError,match='Frozen phase ceiling changed'):
        validate_resources(record)

def test_approved_retry_extension_keeps_original_failures_and_boundaries():
    from reporting import validate_resources
    record=approved_retry_extension_fixture();validate_resources(record)
    assert record['phases']['tests']['attempts'][0]['exit_code']==1
    record['phases']['independent']['wall_seconds']=10201
    with pytest.raises(ValueError,match='Frozen phase ceiling changed'):
        validate_resources(record)

def test_retry_extension_does_not_authorize_other_phase_caps():
    from reporting import validate_resources
    record=approved_retry_extension_fixture();record['phases']['tests']['CPU_seconds']=300
    with pytest.raises(ValueError,match='Frozen phase ceiling changed'):
        validate_resources(record)

def test_retry_extension_rejects_approval_record_tampering():
    from reporting import validate_resources
    record=approved_retry_extension_fixture();record['approved_ceiling_extensions'][0]['approval_record']+=' '
    with pytest.raises(ValueError,match='approval record hash differs'):
        validate_resources(record)

def reconciled_counter_fixture():
    import hashlib,json
    record=resource_fixture();key='test.resources.json';entry=record['receipts'][key]
    original=json.loads(entry['text']);raw=json.loads(entry['text'])
    raw['original_guard_receipt']=dict(entry)
    raw['post_completion_service_observation']=dict(fields=dict(
        ActiveState='failed',Result='timeout',CPUUsageNSec='600000000'))
    raw['final_accounting']['CPU_seconds']=.6
    raw['CPU_reconciliation_delta_seconds']=.6-original['final_accounting']['CPU_seconds']
    text=json.dumps(raw);record['receipts'][key]=dict(text=text,sha256=hashlib.sha256(text.encode()).hexdigest())
    pool=record['phases']['tests'];pool['attempts'][0]['CPU_seconds']=.6;pool['consumed_CPU_seconds']=.6
    return record

def test_late_final_counter_preserves_original_failed_receipt():
    import json
    from reporting import validate_resources
    record=reconciled_counter_fixture();validate_resources(record)
    raw=json.loads(record['receipts']['test.resources.json']['text'])
    original=json.loads(raw['original_guard_receipt']['text'])
    assert original['exit_code']==1 and original['final_accounting']['CPU_seconds']==.5
    assert raw['final_accounting']['CPU_seconds']==.6

def test_late_counter_cannot_rewrite_failure_as_success():
    import hashlib,json
    from reporting import validate_resources
    record=reconciled_counter_fixture();raw=json.loads(record['receipts']['test.resources.json']['text'])
    raw['exit_code']=0;text=json.dumps(raw)
    record['receipts']['test.resources.json']=dict(text=text,sha256=hashlib.sha256(text.encode()).hexdigest())
    record['phases']['tests']['attempts'][0]['exit_code']=0
    with pytest.raises(ValueError,match='Original failed resource receipt changed'):
        validate_resources(record)


def two_approval_fixture():
    # Synthetic policy input only; these fixtures never grant execution authority.
    import hashlib,json
    record=approved_retry_extension_fixture()
    approval=dict(approved=True,scope='issue18_live_replay_extension300',
        original_phase_CPU_seconds=2000,new_phase_CPU_seconds=2300,additional_CPU_seconds=300,
        original_phase_wall_seconds=10200,new_phase_wall_seconds=11100,additional_wall_seconds=900,
        original_attempt_CPU_seconds=900,new_attempt_CPU_seconds=1200,
        original_attempt_wall_seconds=4200,new_attempt_wall_seconds=5100,
        MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0',CPUQuota='25%',threads=1,disk_MiB=256)
    text=json.dumps(approval)
    record['approved_ceiling_extensions'].append(dict(phase='independent',additional_CPU_seconds=300,
        additional_wall_seconds=900,approval_record=text,
        approval_sha256=hashlib.sha256(text.encode()).hexdigest()))
    record['phases']['independent'].update(CPU_seconds=2300,wall_seconds=11100)
    return record

def test_second_approval_requires_ordered_complete_chain():
    from reporting import validate_resources
    record=two_approval_fixture();validate_resources(record)
    record['approved_ceiling_extensions'].reverse()
    with pytest.raises(ValueError,match='ordered resource approval chain'):
        validate_resources(record)
    record=two_approval_fixture();record['approved_ceiling_extensions'].pop(0)
    with pytest.raises(ValueError,match='ordered resource approval chain'):
        validate_resources(record)

def test_same_worker_extension_preserves_stricter_effective_wall_timer():
    import hashlib,json
    from reporting import validate_resources
    record=two_approval_fixture();raw=json.loads(record['receipts']['test.resources.json']['text'])
    raw['budget']=dict(cpu_seconds=1200,wall_seconds=4200,disk_mib=256)
    raw['settings']['RuntimeMaxSec']='4200'
    raw['controller_handoff']=dict(status='Approved_live_extension_applied',
        approval_sha256=record['approved_ceiling_extensions'][1]['approval_sha256'],
        original_budget=dict(cpu_seconds=900,wall_seconds=4200,disk_mib=256),
        approved_budget=dict(cpu_seconds=1200,wall_seconds=5100,disk_mib=256),
        effective_budget=dict(raw['budget']),wall_timer_extension_applied=False,
        scientific_worker_restarted=False,scientific_source_unchanged=True,
        quota_memory_swap_nice_weights_unchanged=True,whole_CPU_counter_retained=True,
        initial_service=dict(MainPID='123'),service_after_handoff=dict(MainPID='123'))
    def bind():
        text=json.dumps(raw);record['receipts']['replay.resources.json']=dict(text=text,
            sha256=hashlib.sha256(text.encode()).hexdigest())
    bind()
    record['phases']['independent'].update(attempts=[dict(receipt='/owned/replay.resources.json',
        CPU_seconds=.5,wall_seconds=2.,exit_code=1)],consumed_CPU_seconds=.5,consumed_wall_seconds=2.)
    validate_resources(record)
    raw['budget']['wall_seconds']=5100;bind()
    with pytest.raises(ValueError,match='stricter effective wall limit differs'):
        validate_resources(record)
