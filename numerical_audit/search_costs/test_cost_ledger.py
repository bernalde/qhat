"""Reject cost artifacts that remain plausible after scalar-field corruption."""
from copy import deepcopy
import pytest
from cost_ledger import audit_design_costs


def artifact():
    def optimizer(initial='cold',partition='conserving'):
        inherited=0. if initial=='cold' else .02 if partition=='conserving' else .03
        return dict(CPU_budget=.5,common_setup_CPU_charge=.01,inherited_placement_CPU_charge=inherited,
            initializer_CPU_seconds=.01,setup_CPU_seconds=.02,search_CPU_seconds=.05,
            accounting_overhead_CPU_seconds=.02,actual_local_CPU_seconds=.1,
            internal_snapshot_CPU_seconds=.09,return_time_CPU_seconds=.01,cost_timing='through run_arm return',
            total_charged_CPU_seconds=.11+inherited,CPU_overrun_seconds=0.,within_CPU_budget=True,
            initialization=initial,evaluations=0,evaluation_cap=120,order=[0,1,2],partition=partition)
    warm=dict(cpu_seconds=.03,internal_snapshot_CPU_seconds=.025,placement_CPU_seconds=.02,
        pauli_conversion_CPU_seconds=.01,cost_timing='through constructor return',charges={'conserving':.02,'pauli':.03})
    warmups=[dict(partition=p,objective=j,initialization=i,record=optimizer(i,p))
        for p in ('conserving','pauli') for j in ('JD','J0') for i in ('cold','placement')]
    records=[]
    for rep in range(3):
        records += [dict(optimizer(),arm='conserving_JD_cold',repetition=rep),
            dict(arm='current_placement',repetition=rep,CPU_budget=None,total_charged_CPU_seconds=.03,
                setup_CPU_seconds=.006,search_CPU_seconds=.008)]
    return dict(warmups=warmups,repetitions=[dict(repetition=r,common_setup_cpu_seconds=.01,placement=deepcopy(warm)) for r in range(3)],records=records)


MODEL={'blocks':[None]*3,'coefficients':[None]*3}


def test_valid_separate_baseline_and_warmup_schemas():
    result=audit_design_costs(artifact(),MODEL)
    assert result['actual_initializer_CPU_seconds']==.03 and result['placement_CPU_seconds']==.02


@pytest.mark.parametrize('field,value,reason',[
    ('return_time_CPU_seconds',-.01,'Negative return'),
    ('within_CPU_budget',False,'feasibility mislabelled')])
def test_corrupted_warmup_timing_or_cap_is_rejected(field,value,reason):
    data=artifact();data['warmups'][0]['record'][field]=value
    with pytest.raises(ValueError,match=reason):audit_design_costs(data,MODEL)


@pytest.mark.parametrize('field,value,reason',[
    ('internal_snapshot_CPU_seconds',.04,'Invalid initializer'),
    ('pauli_conversion_CPU_seconds',.015,'do not close')])
def test_corrupted_initializer_intervals_are_rejected(field,value,reason):
    data=artifact();data['repetitions'][0]['placement'][field]=value
    with pytest.raises(ValueError,match=reason):audit_design_costs(data,MODEL)


def test_uncharged_baseline_cost_is_rejected():
    data=artifact();data['records'][1]['total_charged_CPU_seconds']=.02
    with pytest.raises(ValueError,match='Current placement charge differs'):audit_design_costs(data,MODEL)


def test_inconsistent_warm_start_charges_are_rejected():
    data=artifact();data['warmups'][1]['record']['inherited_placement_CPU_charge']=.03
    with pytest.raises(ValueError,match='Warmup initializer charges differ'):audit_design_costs(data,MODEL)
