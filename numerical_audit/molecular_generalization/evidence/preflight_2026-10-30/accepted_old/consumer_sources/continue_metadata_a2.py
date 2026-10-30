"""One metadata-only continuation; preserve both numerical replays and all costs."""
import importlib.util,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('replay_helpers',HERE/'replay_dispatch.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
out=m.OUT; root=out/'metadata-consumer-a2'
m.require(not root.exists(),'Fresh metadata continuation identity')
prior=[]; cpu=wall=0.0
for planpath in sorted((out/'jobs').glob('*.PLAN.json')):
    plan=m.read(planpath);label=plan['label'];finalpath=out/'jobs'/(label+'.FINAL.json');receiptpath=out/'logs'/(label+'.resources.json')
    final=m.read(finalpath);receipt=m.read(receiptpath);value=m.costs(receipt)
    m.require(final['plan_sha256']==m.sha(planpath) and final['receipt_sha256']==m.sha(receiptpath),'Preserved prior identity')
    m.require(receipt['command']==plan['command'] and all(final[k]==v for k,v in value.items()),'Prior accounting exact')
    if plan['subphase']=='metadata':
        m.require(final['guard_supervisor_exit_code']==receipt['exit_code']==1 and receipt.get('termination')is None,'Retained metadata comparison failure')
    else:m.completed(planpath,finalpath,receiptpath)
    cpu+=value['charged_CPU_seconds'];wall+=value['wall_seconds'];prior.append(dict(plan=m.ref(planpath),final=m.ref(finalpath),receipt=m.ref(receiptpath),status=receipt['exit_code']))
diag=out/'metadata-diagnosis-a1';dc=m.completed(diag/'PLAN.json',diag/'FINAL.json',diag/'logs/a4-replay-metadata-diagnosis-a1.resources.json')
cpu+=dc['charged_CPU_seconds'];wall+=dc['wall_seconds']
m.require(cpu+12<=350 and wall+90<=1800,'Same cumulative350CPU1800wall pool')
root.mkdir();(root/'outputs').mkdir();label='a4-replay-metadata-a2';planpath=root/'PLAN.json';finalpath=root/'FINAL.json';receiptpath=root/'logs'/(label+'.resources.json')
script=HERE/'metadata_consumer_v2.py';analysis=m.DATA/'ANALYSIS.json';python=str(Path(sys.executable));guard=m.PRODUCER/m.STUDY/'tools/resource_guard.py'
command=[python,str(script),'--producer-data',str(m.DATA),'--data',str(m.RELOCATED),'--consumer-root',str(root/'outputs'),'--numerical-root',str(m.CONSUMER),'--producer-analysis',str(analysis),'--metadata-plan',str(planpath)]
invocation=[python,str(guard),'--label',label,'--logs',str(root/'logs'),'--cpu-seconds','12','--wall-seconds','90','--disk-mib','512','--output',str(m.RUN),'--']+command
bound=[analysis,m.DATA/'INPUT_PROTOCOL_FROZEN.json',m.RELOCATED/'EVIDENCE_INVENTORY.json',m.CONSUMER/'independent'/m.CASE/'VERIFY.toml',m.CONSUMER/'objective'/m.CASE/'OBJECTIVE.toml',diag/'DIAGNOSIS.json',script]
m.write(planpath,dict(schema='molecular-external-replay-plan-v1',label=label,phase='package',subphase='metadata-a2',role='later_relocated_metadata_consumer',command=command,guard_command=invocation,cpu_cap=12,wall_cap=90,disk_mib=512,source_freeze_sha256=m.FREEZE,bound_inputs=[m.ref(p)for p in bound],prior_package_CPU_seconds=cpu,prior_package_wall_seconds=wall,prior_jobs=prior,operator=m.ref(__file__),no_new_numerical_replay=True))
status=subprocess.run(invocation,cwd=m.RUN,check=False).returncode
receipt=m.read(receiptpath);value=m.costs(receipt)
m.write(finalpath,dict(schema='molecular-external-replay-final-v1',label=label,phase='package',subphase='metadata-a2',role='later_relocated_metadata_consumer',plan_sha256=m.sha(planpath),receipt_path=str(receiptpath),receipt_sha256=m.sha(receiptpath),guard_supervisor_exit_code=status,**value))
m.completed(planpath,finalpath,receiptpath)
checks=m.read(root/'outputs/REPLAY_CHECKS.json');m.require(checks['complete']is True and checks['logical_rows']==36 and checks['reference_checks']==2 and checks['raw_claims_verified']is True,'Complete actual numerical/report joins')
m.write(out/'REPLAY_FINAL.json',dict(schema='molecular-external-replay-summary-v2',complete=True,role='later_relocated_numerical_consumer',preflight_only=True,source_freeze_sha256=m.FREEZE,checks=m.ref(root/'outputs/REPLAY_CHECKS.json'),original_inventory=m.ref(m.DATA/'EVIDENCE_INVENTORY.json'),archive=m.ref(out/'preflight-a4.data.tar.gz'),extraction=m.ref(m.EXTRACT/'EXTRACTION_VERIFIED.json'),package_CPU_seconds=cpu+value['charged_CPU_seconds'],package_wall_seconds=wall+value['wall_seconds'],prior_jobs=prior,metadata_diagnosis=dict(plan=m.ref(diag/'PLAN.json'),receipt=m.ref(diag/'logs/a4-replay-metadata-diagnosis-a1.resources.json'),final=m.ref(diag/'FINAL.json')),metadata_continuation=dict(plan=m.ref(planpath),receipt=m.ref(receiptpath),final=m.ref(finalpath)),qualification='All prior failure/diagnosis costs retained; archived numerical replays reused without rerun; original scientific fields and digests exact, only known later plot PLAN addition authenticated separately; not main authorization'))
print('OLD_A4_ARCHIVED_NUMERICAL_REPLAY_AND_METADATA_COMPLETE',flush=True)
