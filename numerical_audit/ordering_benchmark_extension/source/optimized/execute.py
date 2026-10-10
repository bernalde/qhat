#!/usr/bin/env python3
"""Provenance-checked freeze, propagation, and completion entry point."""
import argparse,json,resource,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
import run

def execution_sources():return {p.name:c.sha(p) for p in (HERE/'execute.py',HERE/'checked_prepare.py',HERE/'CHECKPOINT_RECOVERY.md',HERE/'check_oracles.py',HERE/'check_archived_schedules.py')}
def qhat_sources():
    paths=('analysis/benchmark_b2_active_spaces_matrix_free.py','analysis/benchmark_b2_signed_coefficient_baseline.py','analysis/benchmark_L_sweep_trotter.py','analysis/state_metrics.py','numerical_audit/run_frozen_audit.py','numerical_audit/reference_tools.py')
    return {p:c.sha(c.QHAT/p) for p in paths}
def freeze_execution():
    import check_archived_schedules
    check_archived_schedules.main()
    run.freeze()
    bindings={}
    for case in c.CASES:
        p=HERE/'checkpoints'/f'{case}.provenance.json';q=HERE/'checkpoints'/f'{case}.completed_hashes.json'
        s=json.loads(p.read_text());completed=json.loads(q.read_text());dpath=HERE/'designs'/f'{case}.design.json';d=json.loads(dpath.read_text());raw,_=c.load_case(case)
        c.require(s['binding']['input_provenance']==raw['extension_input_provenance']==d['extension_input_provenance'],'Inputs changed before execution freeze')
        c.require(s['binding']['scientific_source_sha256']==run.source_hashes(),'Search sources changed before execution freeze')
        c.require(s['binding']['recovery_wrapper_sha256']==c.sha(HERE/'checked_prepare.py'),'Search wrapper changed')
        c.require(completed['binding_sha256']==s['binding_sha256'] and completed['design_sha256']==c.sha(dpath),'Completed preparation binding mismatch')
        for name,h in completed['checkpoints'].items():c.require(c.sha(HERE/'checkpoints'/name)==h,'Completed checkpoint changed')
        bindings[case]={'provenance_sha256':c.sha(p),'completion_sha256':c.sha(q),'input_provenance':raw['extension_input_provenance']}
    value={'execution_source_sha256':execution_sources(),'qhat_source_sha256':qhat_sources(),'schedule_freeze_sha256':c.sha(HERE/'ALL_SCHEDULES_FROZEN.json'),'case_bindings':bindings,'oracle_check_sha256':c.sha(HERE/'oracle_checks.json'),'archived_schedule_checks_sha256':c.sha(HERE/'archived_schedule_checks.json')}
    p=HERE/'EXECUTION_MANIFEST.json'
    if p.exists():c.require(json.loads(p.read_text())==value,'Execution manifest differs')
    else:c.write_atomic(p,value)

def checked_run(case):
    m=json.loads((HERE/'EXECUTION_MANIFEST.json').read_text());c.require(m['execution_source_sha256']==execution_sources(),'Execution sources changed')
    c.require(m['schedule_freeze_sha256']==c.sha(HERE/'ALL_SCHEDULES_FROZEN.json'),'Schedule freeze changed')
    raw,_=c.load_case(case);c.require(m['case_bindings'][case]['input_provenance']==raw['extension_input_provenance'],'Input design/state changed before replay')
    c.require(m['qhat_source_sha256']==qhat_sources(),'QHAT propagation source changed')
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free, state_metrics
    from qhat.numerical_audit import run_frozen_audit, reference_tools
    for module in (benchmark_b2_active_spaces_matrix_free,state_metrics,run_frozen_audit,reference_tools):
        c.require(Path(module.__file__).resolve().is_relative_to(c.QHAT),'Imported unintended QHAT checkout: '+str(module.__file__))
    began=time.monotonic();run.run_case(case)
    rows=list((HERE/'rows'/case).glob('*.json'));c.require(len(rows)==52,'Incomplete case')
    c.write_atomic(HERE/'execution_status'/f'{case}.json',{'case':case,'rows':len(rows),'seconds_this_invocation':time.monotonic()-began,'process_peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_scope':'Entire current propagation process, possibly including previous cases','address_space_limit_GiB':6})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['freeze','run','finish']);ap.add_argument('--cases',nargs='*',default=c.CASES);a=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    if a.command=='freeze':freeze_execution()
    elif a.command=='run':
        for case in a.cases:checked_run(case)
    else:
        freeze_execution();sys.argv=['run.py','finish'];run.main()
        p=HERE/'completion.json';v=json.loads(p.read_text());v['completion_process_peak_rss_KiB']=v.pop('peak_rss_KiB');v['memory_note']='This peak measures the completion process only; per-propagation-process peaks are in execution_status. Each preparation/propagation process had a 6 GiB address-space limit.';c.write_atomic(p,v)
if __name__=='__main__':main()
