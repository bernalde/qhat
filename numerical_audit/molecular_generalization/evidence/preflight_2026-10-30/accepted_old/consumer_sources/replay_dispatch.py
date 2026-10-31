"""DRAFT external OLD-a4 package/replay operator. Root must inspect before use.

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

ROOT = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30')
RUN = ROOT/'runs/preflight-a4'
DATA = RUN/'data'
OUT = RUN/'relocated-replay-a1'
EXTRACT = OUT/'extracted'
RELOCATED = EXTRACT/'data'
CONSUMER = OUT/'consumer'
STUDY = Path('numerical_audit/molecular_generalization')
PRODUCER = DATA/'source_snapshot'
LATER = RELOCATED/'source_snapshot'
COMMIT = '79eff7fe5d2aa7a159d718fa3bd32f77b54eda95'
FREEZE = '679df7485d8adbcce184213a21ae23199c2320d955e0ce81fa98a894fe1797fb'
CASE = 'LiH_r1p487_631g_14m4e'
SPLIT = [('pack',10,60),('extract',8,50),('relocatedverify',200,1000),
         ('relocatedobjective',120,600),('metadata',12,90)]


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


def completed(plan_path, final_path, receipt_path):
    plan, final, receipt = read(plan_path),read(final_path),read(receipt_path)
    require(final['plan_sha256']==sha(plan_path) and final['receipt_sha256']==sha(receipt_path),'PLAN/FINAL/receipt bytes')
    require(receipt['command']==plan['command'] and receipt['budget']==dict(cpu_seconds=plan['cpu_cap'],
        wall_seconds=plan['wall_cap'],disk_mib=plan['disk_mib']),'Actual argv/caps')
    measured=costs(receipt)
    require(all(final[k]==v for k,v in measured.items()),'Actual accounting binding')
    require(final['guard_supervisor_exit_code']==0 and receipt['exit_code']==0 and receipt.get('termination') is None
        and receipt.get('guard_error') is None and not receipt.get('cleanup_errors'),'Closed successful guarded command')
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
        and freeze['mode']=='preflight','Exact old-only A4 freeze')
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
    require(read(analysis)['complete'] is True and read(analysis)['status']=='complete_preflight','Complete original OLD analysis')
    require(gate['visual_artifacts'] and all(sha(v['path'])==v['sha256'] for v in gate['visual_artifacts']),'Root-inspected visual artifacts')
    # All producer jobs, including analysis/figures, must have closed receipts.
    prior_cpu=prior_wall=0.
    for path in sorted((DATA/'jobs').glob('*.PLAN.json')):
        plan=read(path); label=plan['label']
        value=completed(path,DATA/'jobs'/(label+'.FINAL.json'),DATA/'logs'/(label+'.resources.json'))
        if plan['phase']=='package':
            prior_cpu+=value['charged_CPU_seconds']; prior_wall+=value['wall_seconds']
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
    require(sum(c for _,c,_ in SPLIT)+prior_cpu<=limits['CPU_seconds'] and
        sum(w for _,_,w in SPLIT)+prior_wall<=limits['wall_seconds'],'No reset of cumulative package pool')
    require(limits['disk_MiB']==512,'Whole RUN disk declaration')
    guard=PRODUCER/STUDY/'tools/resource_guard.py'
    require(sha(guard)==freeze['sources'][str(STUDY/'tools/resource_guard.py')],'Frozen guard bytes')
    OUT.mkdir(); (OUT/'jobs').mkdir(); (OUT/'logs').mkdir(); CONSUMER.mkdir()
    (CONSUMER/'independent').mkdir(); (CONSUMER/'objective').mkdir()
    archive=OUT/'preflight-a4.data.tar.gz'
    base=[julia,'--startup-file=no','--history-file=no','--threads=1','--gcthreads=1',
          '--heap-size-hint=256M','--project='+str(LATER/'analysis/julia_window')]
    commands={
        'pack':[python,str(PRODUCER/STUDY/'package.py'),'pack',str(DATA),str(archive)],
        'extract':[python,str(PRODUCER/STUDY/'package.py'),'extract',str(archive),str(EXTRACT)],
        'relocatedverify':base+[str(LATER/STUDY/'verify.jl'),str(RELOCATED/'inputs'/CASE/'model.toml'),
            str(RELOCATED/'panels'/(CASE+'.toml')),str(RELOCATED/'physical'/CASE),str(CONSUMER/'independent'/CASE)],
        'relocatedobjective':base+[str(LATER/STUDY/'objective.jl'),str(RELOCATED/'inputs'/CASE/'model.toml'),
            str(RELOCATED/'panels'/(CASE+'.toml')),str(CONSUMER/'objective'/CASE),str(RELOCATED/'witnesses'/(CASE+'.toml'))],
        'metadata':[python,str(Path(__file__).with_name('metadata_consumer.py').resolve()),
            '--producer-data',str(DATA),'--data',str(RELOCATED),'--consumer-root',str(CONSUMER),
            '--producer-analysis',str(analysis),'--metadata-plan',str(OUT/'jobs/a4-replay-metadata.PLAN.json')]}
    write(OUT/'PACKAGE_POOL.json',dict(schema='molecular-external-package-pool-v1',CPU_seconds=350,wall_seconds=1800,
        disk_MiB=512,prior_CPU_seconds=prior_cpu,prior_wall_seconds=prior_wall,
        original_declaration=ref(DATA/'PHASE_LIMITS.json'),source_freeze_sha256=FREEZE,
        producer_visual_gate=ref(args.producer_visual_gate),role='later_relocated_numerical_consumer',
        operator=ref(__file__),metadata_consumer=ref(Path(__file__).with_name('metadata_consumer.py')),
        current_python=dict(version=sys.version,executable=str(Path(sys.executable).resolve()),
                            executable_sha256=sha(sys.executable))))
    actual_cpu,actual_wall=prior_cpu,prior_wall
    for stage,cpu,wall in SPLIT:
        require(cpu<=350-actual_cpu and wall<=1800-actual_wall,'Remaining cumulative package allowance')
        label='a4-replay-'+stage; command=commands[stage]
        plan_path=OUT/'jobs'/(label+'.PLAN.json'); receipt_path=OUT/'logs'/(label+'.resources.json')
        final_path=OUT/'jobs'/(label+'.FINAL.json')
        bound=[DATA/'INPUT_PROTOCOL_FROZEN.json',args.producer_visual_gate,analysis,guard]
        if stage=='extract': bound += [archive,archive.with_name(archive.name+'.json')]
        if stage in ('relocatedverify','relocatedobjective','metadata'):
            bound += [RELOCATED/'inputs'/CASE/'model.toml',RELOCATED/'inputs'/CASE/'model.json',
                RELOCATED/'panels'/(CASE+'.toml'),RELOCATED/'witnesses'/(CASE+'.toml'),
                RELOCATED/'physical'/CASE/'PHYSICAL.toml',RELOCATED/'EVIDENCE_INVENTORY.json']
            for name,digest in freeze['sources'].items():
                require(sha(LATER/name)==digest,'Relocated source bytes changed: '+name)
        script=Path(command[len(base)] if stage.startswith('relocated') else command[1])
        bound.append(script)
        invocation=[python,str(guard),'--label',label,'--logs',str(OUT/'logs'),'--cpu-seconds',str(cpu),
            '--wall-seconds',str(wall),'--disk-mib','512','--memory-profile','numerical','--output',str(RUN),'--']+command
        write(plan_path,dict(schema='molecular-external-replay-plan-v1',label=label,phase='package',subphase=stage,
            role='later_relocated_numerical_consumer',command=command,guard_command=invocation,
            cpu_cap=cpu,wall_cap=wall,disk_mib=512,source_freeze_sha256=FREEZE,
            original_producer_root=freeze['producer_root'],original_producer_data_root=freeze['producer_data_root'],
            execution_root=str(PRODUCER if stage in ('pack','extract') else LATER),monitor_root=str(RUN),
            bound_inputs=[ref(p) for p in bound],guard_sha256=sha(guard),operator=ref(__file__),
            identity=dict(case=CASE,repetition=1),prior_package_CPU_seconds=actual_cpu,prior_package_wall_seconds=actual_wall))
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
            finalization_error=finalization_error,
            failed_cost_policy='Unknown CPU/wall spends full attempt; preserve failed PLAN/log/receipt/FINAL',**value))
        # A nonzero/missing/overrun receipt is retained and stops this fixed attempt.
        value=completed(plan_path,final_path,receipt_path)
        actual_cpu+=value['charged_CPU_seconds']; actual_wall+=value['wall_seconds']
    checks=read(CONSUMER/'REPLAY_CHECKS.json')
    require(checks['complete'] is True and checks['receipt_gate']=='pending_external_post_teardown_FINAL','Actual numerical/report gate')
    write(OUT/'REPLAY_FINAL.json',dict(schema='molecular-external-replay-summary-v1',complete=True,
        role='later_relocated_numerical_consumer',preflight_only=True,source_freeze_sha256=FREEZE,
        checks=ref(CONSUMER/'REPLAY_CHECKS.json'),original_inventory=ref(DATA/'EVIDENCE_INVENTORY.json'),
        archive=ref(archive),extraction=ref(EXTRACT/'EXTRACTION_VERIFIED.json'),
        package_CPU_seconds=actual_cpu,package_wall_seconds=actual_wall,
        jobs=[dict(plan=ref(OUT/'jobs'/('a4-replay-'+s+'.PLAN.json')),
                   receipt=ref(OUT/'logs'/('a4-replay-'+s+'.resources.json')),
                   final=ref(OUT/'jobs'/('a4-replay-'+s+'.FINAL.json'))) for s,_,_ in SPLIT],
        qualification='Producer visual gate plus actual relocated numerical replay; root must inspect terminal receipts and reports; not main authorization'))


if __name__=='__main__':
    main()
