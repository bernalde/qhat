"""DRAFT exact-censored MAIN package/replay dispatcher. Root must inspect before use.

No producer entry point is changed. Every real PLAN is created once immediately
before its guarded command. No retry, new pool, or writing below sealed DATA.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import censor_contract as cc

ROOT = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30')
RUN = ROOT/'runs/main-a1'
DATA = RUN/'data'
OUT = RUN/'relocated-replay-a1'
EXTRACT = OUT/'extracted'
RELOCATED = EXTRACT/'data'
CONSUMER = OUT/'consumer'
STUDY = Path('numerical_audit/molecular_generalization')
PRODUCER = DATA/'source_snapshot'
LATER = RELOCATED/'source_snapshot'
COMMIT = 'aea354e287f2277d8b11cc51695a6d49f7e84537'
FREEZE = 'f3ed635ac7831edcb3feb51d96e6116c2c95f7ad9e2b169b006ca862006776e3'
CASES = ['HF_r1p00_631g_12m8e','H2S_r1p40_a96_631g_12m8e']
FIGURE_PLAN = 'jobs/mg-main-figures-a1.PLAN.json'
SPLIT = [('pack',None,20,100),('extract',None,10,60),
         ('verify-hf',CASES[0],240,1200),('verify-h2s',CASES[1],240,1200),
         ('objective-hf',CASES[0],180,900),('objective-h2s',CASES[1],180,900),
         ('metadata',None,30,150)]


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
    require(Path(path).stat().st_size <= 8*1024**2, 'Bounded JSON metadata')
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False); stream.write('\n')


def ref(path):
    return dict(path=str(Path(path).resolve()), sha256=sha(path))


def costs(receipt):
    # Same frozen construct.cost_from_receipt policy, without importing science.
    values = [receipt.get('process_measurements',{}).get(k) for k in ('user_cpu_seconds','system_cpu_seconds')]
    raw = sum(values) if all(type(v) in (int,float) and math.isfinite(v) and v>=0 for v in values) else None
    sampled=[]
    for sample in receipt.get('samples',[]):
        fields=dict(line.split('=',1) for line in sample.get('service','').splitlines() if '=' in line)
        value=fields.get('CPUUsageNSec','')
        if value.isdigit(): sampled.append(int(value)/1e9)
    charged=None if raw is None else max([raw+.011]+sampled)
    wall=receipt.get('elapsed')
    require(type(wall) in (int,float) and math.isfinite(wall) and wall>=0,'Measured guarded wall')
    return dict(raw_CPU_seconds=raw,charged_CPU_seconds=charged,wall_seconds=wall)


def completed(plan_path, final_path, receipt_path, *, expected_exit=0):
    plan, final, receipt = read(plan_path),read(final_path),read(receipt_path)
    require(final['plan_sha256']==sha(plan_path) and final['receipt_sha256']==sha(receipt_path),'PLAN/FINAL/receipt bytes')
    require(receipt['command']==plan['command'] and receipt['budget']==dict(cpu_seconds=plan['cpu_cap'],
        wall_seconds=plan['wall_cap'],disk_mib=plan['disk_mib']),'Actual argv/caps')
    measured=costs(receipt)
    require(all(final[k]==v for k,v in measured.items()),'Actual accounting binding')
    require(expected_exit in (0,2) and final['guard_supervisor_exit_code']==expected_exit and
        receipt['exit_code']==receipt['supervisor_exit_code']==expected_exit and receipt.get('termination') is None
        and receipt.get('guard_error') is None and not receipt.get('cleanup_errors'),'Exact expected payload outcome; no guard/cap/cleanup failure')
    require(measured['charged_CPU_seconds'] is not None and measured['charged_CPU_seconds']<=plan['cpu_cap']
        and measured['wall_seconds']<=plan['wall_cap'],'Strict final CPU/wall caps')
    require(receipt['final_accounting']['CPU_accounting_available'] is True and
        receipt['final_accounting']['process_CPU_timing_available'] is True and
        receipt['final_accounting']['output_bytes']<=plan['disk_mib']*1024**2,'Terminal resource evidence')
    for key,value in dict(CPUQuota='25%',MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0',CPUWeight='5',IOWeight='10',Nice='19').items():
        require(str(receipt['settings'][key])==value,'Frozen numerical resource profile')
    return measured


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-visual-gate',type=Path,required=True)
    args=parser.parse_args()
    require(not OUT.exists(),'Fresh fixed replay identity; failed attempts cannot be replaced')
    freeze=read(DATA/'INPUT_PROTOCOL_FROZEN.json')
    require(sha(DATA/'INPUT_PROTOCOL_FROZEN.json')==FREEZE and freeze['source_commit']==COMMIT
        and freeze['mode']=='main' and freeze['configuration']['main_execution_allowed'] is True
        and freeze['configuration']['protocol_state']=='frozen_main'
        and [v['case'] for v in freeze['registry']['cases']]==CASES,'Exact frozen MAIN identity')
    python=freeze['runtime']['executable']; julia=freeze['runtime']['julia']['executable']
    require(Path(sys.executable).resolve()==Path(python).resolve(),'Pinned Python interpreter')
    require(sha(julia)==freeze['runtime']['julia']['executable_sha256'],'Pinned Julia executable')
    for name,digest in freeze['sources'].items():
        require(sha(PRODUCER/name)==digest,'Original frozen source bytes changed: '+name)
    gate=read(args.producer_visual_gate)
    require(gate['producer_pipeline_closed'] is True and gate['producer_visual_inspection_complete'] is True
        and gate['source_freeze_sha256']==FREEZE and gate['source_commit']==COMMIT,'Root producer/visual gate required')
    analysis=Path(gate['producer_analysis']['path']).resolve(strict=True)
    require(analysis.is_relative_to(DATA) and sha(analysis)==gate['producer_analysis']['sha256'],'Original analysis bytes in sealed DATA')
    cc.audit_analysis(read(analysis))
    require(gate['scientific_gate_complete'] is False,'Root gate must preserve known incomplete scientific evidence')
    cc.original_objective_bindings(DATA)
    expected_analysis=cc.analysis_job_binding(DATA,gate,freeze,analysis)
    hf_owner=cc.EXPECTED['reports'][cc.HF]['owner']['plan']['path']
    require(expected_analysis!=hf_owner,'Distinct exact expected nonzero producer owners')
    expected_nonzero={hf_owner,expected_analysis}
    require(gate['visual_artifacts'] and all(sha(v['path'])==v['sha256'] for v in gate['visual_artifacts']),'Root-inspected visual artifacts')
    ancillary=gate['post_analysis_figure_plan']
    if ancillary is not None:
        require(set(ancillary)=={'path','sha256'} and ancillary['path']==FIGURE_PLAN and
                sha(DATA/FIGURE_PLAN)==ancillary['sha256'],'Only exact root-bound main figure PLAN')

    # All producer jobs, including analysis/figures, must have closed receipts.
    prior_cpu=prior_wall=0.
    seen_nonzero=set(); original_nonzero_jobs=[]
    for path in sorted((DATA/'jobs').glob('*.PLAN.json')):
        plan=read(path); label=plan['label']
        relative=str(path.relative_to(DATA)); expected_exit=2 if relative in expected_nonzero else 0
        value=completed(path,DATA/'jobs'/(label+'.FINAL.json'),DATA/'logs'/(label+'.resources.json'),expected_exit=expected_exit)
        if expected_exit==2:
            seen_nonzero.add(relative)
            original_nonzero_jobs.append(dict(plan=ref(path),final=ref(DATA/'jobs'/(label+'.FINAL.json')),
                receipt=ref(DATA/'logs'/(label+'.resources.json')),actual_exit_code=2,**value))
        if plan['phase']=='package':
            prior_cpu+=value['charged_CPU_seconds']; prior_wall+=value['wall_seconds']
    require(seen_nonzero==expected_nonzero,'Both exact known source nonzero jobs retained; no other nonzero accepted')
    # Root inventories ANY earlier external package attempt in the same pool.
    # Failed/unknown costs remain spent; an empty list is an explicit assertion.
    earlier=gate['earlier_external_package_attempts']
    require(isinstance(earlier,list),'Cumulative earlier-package inventory')
    seen=set()
    for attempt in earlier:
        require(set(attempt)=={'plan','final','receipt'},'Earlier attempt references')
        for value in attempt.values():
            if value is not None:
                require(sha(value['path'])==value['sha256'],'Earlier attempt artifact bytes')
        p=Path(attempt['plan']['path']).resolve(); f=read(attempt['final']['path']); prior=read(p)
        require(p not in seen and not p.is_relative_to(DATA/'jobs') and prior['phase']=='package',
                'Earlier package attempt counted once')
        seen.add(p)
        require(f['plan_sha256']==sha(p) and f['receipt_sha256']==(
                None if attempt['receipt'] is None else attempt['receipt']['sha256']),
                'Earlier failed/successful cost binding')
        require(attempt['receipt'] is not None or
                f['charged_CPU_seconds'] is None and f['wall_seconds'] is None,
                'Missing earlier receipt consumes its full declared reservation')
        prior_cpu+=prior['cpu_cap'] if f['charged_CPU_seconds'] is None else f['charged_CPU_seconds']
        prior_wall+=prior['wall_cap'] if f['wall_seconds'] is None else f['wall_seconds']
    limits=read(DATA/'PHASE_LIMITS.json')['phases']['package']
    require(sum(c for _,_,c,_ in SPLIT)+prior_cpu<=limits['CPU_seconds'] and
        sum(w for _,_,_,w in SPLIT)+prior_wall<=limits['wall_seconds'],'900CPU/4510wall reservation plus all prior spending must fit original1000/5000 pool')
    require(limits==dict(CPU_seconds=1000,wall_seconds=5000,disk_MiB=1024),'Original cumulative package pool')
    guard=PRODUCER/STUDY/'tools/resource_guard.py'
    require(sha(guard)==freeze['sources'][str(STUDY/'tools/resource_guard.py')],'Frozen guard bytes')
    OUT.mkdir(); (OUT/'jobs').mkdir(); (OUT/'logs').mkdir(); CONSUMER.mkdir()
    (CONSUMER/'independent').mkdir(); (CONSUMER/'objective').mkdir()
    archive=OUT/'main-a1.data.tar.gz'
    base=[julia,'--startup-file=no','--history-file=no','--threads=1','--gcthreads=1',
          '--heap-size-hint=256M','--project='+str(LATER/'analysis/julia_window')]
    commands={
        'pack':[python,str(PRODUCER/STUDY/'package.py'),'pack',str(DATA),str(archive)],
        'extract':[python,str(PRODUCER/STUDY/'package.py'),'extract',str(archive),str(EXTRACT)],
        'metadata':[python,str(Path(__file__).with_name('metadata_consumer.py').resolve()),
            '--producer-data',str(DATA),'--data',str(RELOCATED),'--consumer-root',str(CONSUMER),
            '--numerical-root',str(CONSUMER),'--producer-analysis',str(analysis),
            '--producer-visual-gate',str(args.producer_visual_gate.resolve()),
            '--metadata-plan',str(OUT/'jobs/main-replay-metadata.PLAN.json')]}
    for stage,case,_,_ in SPLIT:
        if stage.startswith('verify-'):
            commands[stage]=base+[str(LATER/STUDY/'verify.jl'),str(RELOCATED/'inputs'/case/'model.toml'),
                str(RELOCATED/'panels'/(case+'.toml')),str(RELOCATED/'physical'/case),str(CONSUMER/'independent'/case)]
        elif stage.startswith('objective-'):
            commands[stage]=base+[str(LATER/STUDY/'objective.jl'),str(RELOCATED/'inputs'/case/'model.toml'),
                str(RELOCATED/'panels'/(case+'.toml')),str(CONSUMER/'objective'/case),str(RELOCATED/'witnesses'/(case+'.toml'))]
    write(OUT/'PACKAGE_POOL.json',dict(schema='molecular-external-package-pool-v1',CPU_seconds=1000,wall_seconds=5000,
        disk_MiB=1024,prior_CPU_seconds=prior_cpu,prior_wall_seconds=prior_wall,
        original_declaration=ref(DATA/'PHASE_LIMITS.json'),source_freeze_sha256=FREEZE,
        producer_visual_gate=ref(args.producer_visual_gate),role='later_relocated_numerical_consumer',
        operator=ref(__file__),metadata_consumer=ref(Path(__file__).with_name('metadata_consumer.py')),
        censor_contract=ref(cc.__file__),expected_outcomes=ref(Path(cc.__file__).with_name('EXPECTED_OUTCOMES.json')),
        original_expected_nonzero_jobs=original_nonzero_jobs,scientific_gate_complete=False,
        current_python=dict(version=sys.version,executable=str(Path(sys.executable).resolve()),
                            executable_sha256=sha(sys.executable))))
    actual_cpu,actual_wall=prior_cpu,prior_wall
    for stage,case,cpu,wall in SPLIT:
        require(cpu<=1000-actual_cpu and wall<=5000-actual_wall,'Remaining cumulative package allowance')
        label='main-replay-'+stage; command=commands[stage]
        plan_path=OUT/'jobs'/(label+'.PLAN.json'); receipt_path=OUT/'logs'/(label+'.resources.json')
        final_path=OUT/'jobs'/(label+'.FINAL.json')
        expected_exit=2 if stage=='objective-hf' else 0
        bound=[DATA/'INPUT_PROTOCOL_FROZEN.json',args.producer_visual_gate,analysis,guard,
               Path(cc.__file__),Path(cc.__file__).with_name('EXPECTED_OUTCOMES.json')]
        if ancillary is not None:
            bound.append(DATA/FIGURE_PLAN)
        if stage=='extract': bound += [archive,archive.with_name(archive.name+'.json')]
        if case is not None or stage=='metadata':
            for bound_case in ([case] if case else CASES):
                bound += [RELOCATED/'inputs'/bound_case/'model.toml',RELOCATED/'inputs'/bound_case/'model.json',
                    RELOCATED/'panels'/(bound_case+'.toml'),RELOCATED/'witnesses'/(bound_case+'.toml'),
                    RELOCATED/'physical'/bound_case/'PHYSICAL.toml']
                if stage=='metadata':
                    bound += [CONSUMER/'independent'/bound_case/'VERIFY.toml',
                              CONSUMER/'objective'/bound_case/'OBJECTIVE.toml']
            bound.append(RELOCATED/'EVIDENCE_INVENTORY.json')
            for name,digest in freeze['sources'].items():
                require(sha(LATER/name)==digest,'Relocated source bytes changed: '+name)
        script=Path(command[len(base)] if case is not None else command[1])
        bound.append(script)
        invocation=[python,str(guard),'--label',label,'--logs',str(OUT/'logs'),'--cpu-seconds',str(cpu),
            '--wall-seconds',str(wall),'--disk-mib','1024','--memory-profile','numerical','--output',str(RUN),'--']+command
        write(plan_path,dict(schema='molecular-external-replay-plan-v1',label=label,phase='package',subphase=stage,
            role='later_relocated_numerical_consumer',command=command,guard_command=invocation,
            cpu_cap=cpu,wall_cap=wall,disk_mib=1024,source_freeze_sha256=FREEZE,
            expected_payload_exit_code=expected_exit,scientific_gate_complete=False,
            original_producer_root=freeze['producer_root'],original_producer_data_root=freeze['producer_data_root'],
            execution_root=str(PRODUCER if stage in ('pack','extract') else LATER),monitor_root=str(RUN),
            bound_inputs=[ref(p) for p in bound],guard_sha256=sha(guard),operator=ref(__file__),
            identity=(dict(case=case,repetitions=2) if case else dict(cases=CASES,repetitions=2)),prior_package_CPU_seconds=actual_cpu,prior_package_wall_seconds=actual_wall))
        # The only scientific dispatch, executed later by root, is this guard.
        status=subprocess.run(invocation,cwd=RUN,check=False).returncode
        finalization_error=None
        receipt=None
        value=dict(raw_CPU_seconds=None,charged_CPU_seconds=None,wall_seconds=None)
        try:
            receipt=read(receipt_path) if receipt_path.is_file() else None
            if receipt is not None: value=costs(receipt)
        except Exception as error:
            finalization_error=str(error)
        write(final_path,dict(schema='molecular-external-replay-final-v1',label=label,phase='package',subphase=stage,
            role='later_relocated_numerical_consumer',plan_sha256=sha(plan_path),receipt_path=str(receipt_path),
            receipt_sha256=sha(receipt_path) if receipt_path.is_file() else None,guard_supervisor_exit_code=status,
            finalization_error=finalization_error,expected_payload_exit_code=expected_exit,
            failed_cost_policy='Unknown CPU/wall spends full attempt; preserve failed PLAN/log/receipt/FINAL',**value))
        # Only the exact HF objective payload2 is expected; all other nonzero/missing/overrun outcomes stop the fixed attempt.
        value=completed(plan_path,final_path,receipt_path,expected_exit=expected_exit)
        if stage.startswith('objective-'):
            # Stop before another dispatch if exit2 represents any different
            # refusal/mismatch/floor/shape; full joins remain in metadata guard.
            cc.audit_objective_file(case,CONSUMER/'objective'/case/'OBJECTIVE.toml')
        actual_cpu+=value['charged_CPU_seconds']; actual_wall+=value['wall_seconds']
    checks=read(CONSUMER/'REPLAY_CHECKS.json')
    require(checks['delivery_complete'] is True and checks['replay_complete'] is True
        and checks['scientific_gate_complete'] is False and checks['raw_claims_verified'] is False
        and checks['raw_proposed_comparisons_verified'] is True and checks['selected_objectives_verified']==144
        and checks['receipt_gate']=='pending_external_post_teardown_FINAL'
        and checks['source_freeze_sha256']==FREEZE and checks['logical_rows']==144
        and checks['reference_checks']==4 and [value['case'] for value in checks['cases']]==CASES
        and all(value['report_processing_complete'] is True for value in checks['cases']),
        'Exact censored numerical delivery replay; scientific gate remains incomplete')
    write(OUT/'REPLAY_FINAL.json',dict(schema='molecular-external-censored-replay-summary-v1',
        delivery_complete=True,replay_complete=True,scientific_gate_complete=False,source_analysis_complete=False,
        raw_claims_verified=False,raw_proposed_comparisons_verified=True,selected_objectives_verified=144,
        source_primary_decisions=dict(window_minus_full='incomplete',annealing_minus_full='incomplete'),
        role='later_relocated_numerical_consumer',preflight_only=False,mode='main',source_freeze_sha256=FREEZE,
        expected_outcomes=ref(Path(cc.__file__).with_name('EXPECTED_OUTCOMES.json')),
        original_expected_nonzero_jobs=original_nonzero_jobs,
        checks=ref(CONSUMER/'REPLAY_CHECKS.json'),original_inventory=ref(DATA/'EVIDENCE_INVENTORY.json'),
        archive=ref(archive),extraction=ref(EXTRACT/'EXTRACTION_VERIFIED.json'),
        package_CPU_seconds=actual_cpu,package_wall_seconds=actual_wall,
        jobs=[dict(plan=ref(OUT/'jobs'/('main-replay-'+s+'.PLAN.json')),
                   receipt=ref(OUT/'logs'/('main-replay-'+s+'.resources.json')),
                   final=ref(OUT/'jobs'/('main-replay-'+s+'.FINAL.json'))) for s,_,_,_ in SPLIT],
        qualification='Delivery/replay processing only: exact original HF unpromoted native-refused witness remains censored, source ANALYSIS and both primary decisions remain incomplete; all actual nonzero statuses/costs retained; no scientific-completeness or winner claim'))


if __name__=='__main__':
    main()
