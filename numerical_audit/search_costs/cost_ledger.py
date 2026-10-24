"""Audit retained construction costs without changing frozen scientific source."""
import math


def require(ok,message):
    if not ok:raise ValueError(message)


def close(a,b,message):require(abs(a-b)<1e-12,message)


def optimizer(record,units):
    from verify_costs import audit_ledger
    audit_ledger(dict(record,arm=record.get('arm','warmup')),units)
    require(all(math.isfinite(record[k]) and record[k]>=0 for k in
        ('common_setup_CPU_charge','inherited_placement_CPU_charge','initializer_CPU_seconds',
         'setup_CPU_seconds','search_CPU_seconds','actual_local_CPU_seconds','return_time_CPU_seconds')),'Invalid optimizer CPU component')
    close(record['CPU_overrun_seconds'],max(0.,record['total_charged_CPU_seconds']-record['CPU_budget']),'Overrun amount differs')


def initializer(warm):
    require(warm['cost_timing']=='through constructor return','Initializer timing boundary differs')
    total=warm['cpu_seconds'];internal=warm['internal_snapshot_CPU_seconds'];placement=warm['placement_CPU_seconds']
    conversion=warm['pauli_conversion_CPU_seconds']
    require(all(math.isfinite(v) for v in (total,internal,placement,conversion)) and total>=internal>=placement>=0 and conversion>=0,'Invalid initializer timing interval')
    close(placement+conversion,total,'Initializer placement/conversion costs do not close')
    require(warm['charges']=={'conserving':placement,'pauli':total},'Initializer charges differ')


def warmup_shared(design):
    common={r['record']['common_setup_CPU_charge'] for r in design['warmups']}
    require(len(common)==1 and next(iter(common))>=0,'Warmup common charge differs')
    charges={}
    for partition in ('conserving','pauli'):
        values={r['record']['inherited_placement_CPU_charge'] for r in design['warmups']
            if r['partition']==partition and r['initialization']=='placement'}
        require(len(values)==1 and next(iter(values))>0,'Warmup initializer charges differ')
        charges[partition]=next(iter(values))
    require(charges['pauli']>=charges['conserving'],'Negative warmup conversion cost')
    return dict(common_setup_CPU_seconds=next(iter(common)),placement_CPU_seconds=charges['conserving'],
        pauli_conversion_CPU_seconds=charges['pauli']-charges['conserving'],actual_initializer_CPU_seconds=charges['pauli'],
        source='Retained warmup charges; shared-constructor internal snapshots were not retained')


def audit_design_costs(design,model):
    shared=warmup_shared(design)
    for item in design['warmups']:
        record=item['record'];partition=item['partition'];units=len(model['blocks']) if partition=='conserving' else len(model['coefficients'])
        require(record['CPU_budget']==.5,'Warmup allowance differs')
        optimizer(record,units)
        expected=0. if item['initialization']=='cold' else shared['placement_CPU_seconds'] if partition=='conserving' else shared['actual_initializer_CPU_seconds']
        close(record['inherited_placement_CPU_charge'],expected,'Warmup placement binding differs')
    setups={r['repetition']:r for r in design['repetitions']};require(set(setups)==set(range(3)),'Shared setup coverage differs')
    for setup in setups.values():
        initializer(setup['placement']);require(setup['common_setup_cpu_seconds']>=0,'Negative common setup CPU')
    for record in design['records']:
        setup=setups[record['repetition']];warm=setup['placement'];common=setup['common_setup_cpu_seconds']
        if record['arm']=='current_placement':
            require(record['CPU_budget'] is None,'Current placement assigned a bounded-arm cap')
            close(record['total_charged_CPU_seconds'],common+warm['placement_CPU_seconds'],'Current placement charge differs')
            require(warm['placement_CPU_seconds']+1e-12>=record['setup_CPU_seconds']+record['search_CPU_seconds'],'Current internal phases exceed caller cost')
        else:
            optimizer(record,len(record['order']))
            close(record['common_setup_CPU_charge'],common,'Measured common setup charge differs')
            expected=0. if record['initialization']=='cold' else warm['charges'][record['partition']]
            close(record['inherited_placement_CPU_charge'],expected,'Measured initializer binding differs')
    return shared
