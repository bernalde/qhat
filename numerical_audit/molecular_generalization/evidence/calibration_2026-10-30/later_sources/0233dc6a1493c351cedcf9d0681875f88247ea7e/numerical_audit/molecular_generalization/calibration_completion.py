"""Bounded completion of the actual closed a1 calibration stop.

Reuse original native66 and complete CO33; run unchanged calibration.jl once
for CH4 parent+32 probes using original a1 input/request bytes. No native or CO
rescoring, chemistry, physical errors or retuning. Promote/review/test/commit
this entrypoint in the numerical repo before prepare/work/finalize dispatch.
"""
from __future__ import annotations
import time
BOOT_WALL=time.monotonic()
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tomllib

ORIGINAL_COMMIT='4ba93256e11aa89038b91c907189c875771c6870'
CASES=('CO_r1p13','CH4_r1p09')
PROFILE=dict(CPUQuota='25%',CPUWeight='5',IOWeight='10',Nice='19',
             MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0')
class Refusal(ValueError):pass
def require(ok,message):
    if not ok:raise Refusal(message)
def raw(path,maximum=32*1024**2):
    path=Path(path);require(path.is_file() and path.stat().st_size<=maximum,'Bounded artifact: '+str(path))
    with path.open('rb') as stream:value=stream.read(maximum+1)
    require(len(value)<=maximum,'Growing artifact byte cap');return value
def sha(path):return hashlib.sha256(raw(path)).hexdigest()
def read(path):
    def unique(items):
        result={}
        for key,value in items:require(key not in result,'Duplicate JSON key');result[key]=value
        return result
    def nonfinite(value):raise Refusal('Nonfinite JSON: '+value)
    return json.loads(raw(path),object_pairs_hook=unique,parse_constant=nonfinite)
def ref(path):
    path=Path(path).resolve();return dict(path=str(path),sha256=sha(path))
def pin(reference):
    path=Path(reference['path']);require(sha(path)==reference['sha256'],'Pinned artifact changed: '+str(path));return path
def write_new(path,data):
    import os
    path=Path(path);temporary=path.with_name(path.name+'.partial')
    require(not path.exists() and not temporary.exists(),'Fresh artifact required')
    path.parent.mkdir(parents=True,exist_ok=True)
    encoded=(json.dumps(data,indent=2,allow_nan=False)+'\n').encode('ascii')
    require(len(encoded)<=32*1024**2,'Completion JSON byte cap')
    with temporary.open('xb') as stream:stream.write(encoded)
    os.link(temporary,path);temporary.unlink();return sha(path)
def load_c(root):
    path=Path(root)/'numerical_audit/molecular_generalization/calibration.py'
    spec=importlib.util.spec_from_file_location('_original_calibration_completion',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    require(module.ROOT==Path(root).resolve(),'Exact original scorer root');return module
def tracked_driver(tree,commit):
    tree=Path(tree).resolve();driver=Path(__file__).resolve()
    require(driver.is_relative_to(tree),'Promote/review/test the private draft into the numerical repo first')
    head=subprocess.check_output(['git','-C',str(tree),'rev-parse','HEAD'],text=True).strip()
    require(head==commit and re.fullmatch(r'[0-9a-f]{40}',commit),'Exact committed completion producer')
    name=str(driver.relative_to(tree));content=subprocess.check_output(['git','-C',str(tree),'show',commit+':'+name])
    require(hashlib.sha256(content).hexdigest()==sha(driver),'Commit completion driver before dispatch')
    return dict(commit=commit,driver=ref(driver),repository_path=name,tree=str(tree))
def duration(text):
    match=re.fullmatch(r'(?:(\d+)h\s*)?(?:(\d+)min\s*)?(\d+(?:\.\d+)?)s',text.strip())
    require(match is not None,'Unsupported terminal service CPU spelling')
    h,m,s=match.groups();digits=len(s.split('.')[1]) if '.' in s else 0
    return 3600*int(h or 0)+60*int(m or 0)+float(s),.5*10**(-digits)
def receipt_cost(path,command,*,allow_failure):
    r=read(path);require(r['command']==command and isinstance(r.get('completed_utc'),str),'Closed owning command receipt')
    require(all(r['settings'].get(k)==v for k,v in PROFILE.items()),'Same quarter-core/448-512MiB/no-swap profile')
    require(all(r['environment'].get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')),'One numerical thread')
    require(type(r['exit_code']) is int and type(r['elapsed']) in (int,float) and math.isfinite(r['elapsed']) and r['elapsed']>=0,'Terminal receipt')
    require(r.get('guard_error') is None and not r.get('cleanup_errors'),'Clean guard teardown required')
    if not allow_failure:require(r['exit_code']==0 and r.get('termination') is None,'Successful completion work')
    log=raw(r['log']).decode('utf-8');values=[];upper=[];p=r.get('process_measurements',{});a=r['final_accounting']
    if a.get('process_CPU_timing_available') is True:
        require(all(type(p.get(k)) in (int,float) and math.isfinite(p[k]) and p[k]>=0
                    for k in ('user_cpu_seconds','system_cpu_seconds')),'GNU time terminal proof')
        value=p['user_cpu_seconds']+p['system_cpu_seconds'];values.append(value);upper.append(value+.011)
    footers=re.findall(r'^CPU time consumed:\s*(.+)$',log,re.MULTILINE)
    if footers:
        value,pad=duration(footers[-1]);values.append(value);upper.append(value+pad)
    observed=a.get('CPU_seconds')
    require(observed is None or type(observed) in (int,float) and math.isfinite(observed) and observed>=0,'Finite observed CPU')
    require(values,'Unknown terminal CPU blocks admission; never charge zero')
    charged=math.ceil(max(upper+[observed or 0.])*1000)/1000
    require(type(a.get('output_bytes')) is int and a['output_bytes']>=0,'Retained disk accounting')
    return r,dict(raw_CPU_seconds=max(values),charged_CPU_seconds=charged,wall_seconds=r['elapsed'],
        output_bytes=a['output_bytes'],receipt=ref(path),log=ref(r['log']),failed_attempt=bool(r['exit_code'] or r.get('termination')))
def score_valid(score,order,model):
    require(isinstance(score,dict) and set(score)=={'value','cubic_response','steps'},'Complete native score')
    require(all(type(score[k]) in (int,float) and math.isfinite(score[k]) and score[k]>=0 for k in ('value','cubic_response')),'Finite native score')
    require(type(score['steps']) is int and score['steps']>0 and not(score['cubic_response']>0 and score['value']==0),'Native steps/underflow')
    lengths=[len(model['blocks'][i]) for i in order];U=sum(lengths)
    steps=(30*U-lengths[0])//(2*U-lengths[0]-lengths[-1])
    require(score['steps']==steps and abs(score['value']-score['cubic_response']/steps**4)<=max(1e-30,abs(score['value'])*1e-12),'Original endpoint fusion/normalization')
def validate_original(root,source_ref,original):
    source=read(pin(source_ref));require(source['source_commit']==ORIGINAL_COMMIT and source['old_cases_only'] is True
        and source['new_physical_errors_observed'] is False,'Original producer identity')
    root=Path(root).resolve()
    # Authenticate the FULL historical map against immutable original Git
    # objects, not a later whole-project tree. Unrelated prospective protocol/
    # physical-adapter amendments are allowed in the completion checkout.
    archived={}
    for name,digest in source['sources'].items():
        path=(root/name).resolve();require(path.is_relative_to(root) and '\n' not in name,'Canonical historical source path')
        content=subprocess.check_output(['git','-C',str(root),'show',ORIGINAL_COMMIT+':'+name])
        require(len(content)<=32*1024**2 and hashlib.sha256(content).hexdigest()==digest,'Original Git/source-map identity: '+name)
        if name=='numerical_audit/molecular_generalization/ENVIRONMENT.json':archived['runtime']=json.loads(content)
    name='numerical_audit/molecular_generalization/calibration.py'
    require(sha(root/name)==source['sources'][name],'Original calibration loader bytes must remain unchanged')
    c=load_c(root);numerical={name:source['sources'][name] for name in c.SOURCE_NAMES}
    for name,digest in numerical.items():require(sha(root/name)==digest,'Live numerical calibration dependency changed: '+name)
    for name in ('analysis/__init__.py','numerical_audit/__init__.py'):
        require(sha(root/name)==source['sources'][name],'Implicit original package import changed')
    original=Path(original);declaration=read(original/'DECLARATION.json');oldreport=read(original/'CALIBRATION.json')
    require(declaration['case_order']==list(CASES) and declaration['source_sha256']==numerical
        and declaration['probe_count_per_case']==32 and declaration['minimum_positive_deltas']==1
        and declaration['initial_acceptance']==.8 and declaration['final_acceptance']==.01,'Original fixed fit/probe policy')
    require(oldreport['complete'] is False and oldreport['source_sha256']==numerical and
        oldreport['declaration_sha256']==sha(original/'DECLARATION.json'),'Original incomplete report')
    jobs=oldreport['Julia_jobs'];require(len(jobs)==2 and all(job['reaped'] is True for job in jobs),'Both original children reaped')
    require(jobs[0]['case']==CASES[0] and jobs[0]['stop']=='complete' and jobs[0]['returncode']==0
        and jobs[1]['case']==CASES[1] and jobs[1]['stop']=='calibration_cpu_limit' and jobs[1]['returncode']==-9,
        'This lean contingency covers only the actual closed a1 state')
    runtime=archived['runtime']
    version=runtime['julia']['version'].removeprefix('julia version ');cases=[]
    for spec,declared in zip(c.CASES,declaration['cases']):
        case=spec['case'];model_path=original/'inputs'/(case+'.toml');model=c.read_model(model_path,spec)
        probe=c.probe_plan(case,model['parent']);require(all(declared[k]==probe[k] for k in probe),'Original deterministic probes')
        request=original/(case+'.REQUEST.toml');q=tomllib.loads(raw(request,8*1024**2).decode('ascii'))
        require(q['orders']==probe['orders'] and q['pairs']==probe['pairs'] and q['source_sha256']==numerical
            and q['model_sha256']==spec['sha256'] and q['seed']==probe['seed'],'Original request identity')
        native_refs=[];julia_refs=[]
        for index,order in enumerate(probe['orders']):
            path=original/'native'/case/('row'+str(index).zfill(3)+'.json');row=read(path)
            require(row['case']==case and type(row['index']) is int and row['index']==index and row['order']==order
                and row['model_sha256']==spec['sha256'] and row['legacy_source_model_sha256']==spec['legacy_source_model_sha256'],'Every original native identity')
            score_valid(row['score'],order,model);require(q['native_scores'][index]==row['score'],'Original request/native join')
            native_refs.append(ref(path));path=original/'julia'/case/('row'+str(index).zfill(3)+'.toml')
            if case==CASES[0]:
                checked=c._admit_child_row(path,case,index,order,row,spec['sha256'],sha(request),numerical)
                require(checked['verified'] and checked['julia_version']==version,'Original CO triple parity/runtime/floors')
                julia_refs.append(ref(path))
            else:require(not path.exists(),'Unexpected CH4 row: re-plan, do not duplicate completed work')
        cases.append(dict(spec=spec,model=ref(model_path),probe=probe,request=ref(request),native=native_refs,original_Julia=julia_refs))
    require(not list((original/'julia'/CASES[1]).glob('row*.toml')),'Exactly zero actual CH4 Julia rows')
    return c,dict(root=str(original.resolve()),scorer_root=str(root),source=source_ref,source_commit=ORIGINAL_COMMIT,
        declaration=ref(original/'DECLARATION.json'),report=ref(original/'CALIBRATION.json'),numerical_sources=numerical,runtime=runtime,cases=cases)
def prepare(args):
    producer=tracked_driver(args.completion_tree,args.completion_commit)
    c,original=validate_original(args.scorer_root,ref(args.original_source),args.original)
    runtime=original['runtime'];r,cost=receipt_cost(args.original_receipt,
        [runtime['executable'],str(Path(args.original_wrapper).resolve())],allow_failure=True)
    require(r['exit_code']==2 and r.get('termination') is None and args.a1_quiescent,'Closed a1 exit-2 and coordinator quiescence confirmation')
    original.update(wrapper=ref(args.original_wrapper),guard_receipt=ref(args.original_receipt),cost=cost,
        standalone_PLAN=None,standalone_PLAN_absence_explicit=True,prospective_task_pool_declaration=ref(args.launch_declaration),
        quiescence_basis='Closed receipt, both exact owned child wait4 records reaped, coordinator confirmation')
    for x in (args.cumulative_cpu,args.cumulative_wall,args.attempt_cpu,args.attempt_wall):require(math.isfinite(x) and x>0,'Finite amended resources')
    require(0<args.attempt_cpu<=args.cumulative_cpu-cost['charged_CPU_seconds'] and
        0<args.attempt_wall<=args.cumulative_wall-cost['wall_seconds'],'No reset: attempt fits remaining amended pool')
    disk=args.disk_mib-math.ceil(cost['output_bytes']/1024**2)
    require(type(args.disk_mib) is int and 0<disk<=1024,'Finite remaining disk allowance')
    output=Path(args.output).resolve();old=Path(original['root'])
    require(not output.exists() and not output.is_relative_to(old) and not old.is_relative_to(output),'Fresh disjoint output')
    plan=dict(schema='calibration-completion-plan-v1',producer=producer,original=original,output=str(output),
        cumulative_pool=dict(CPU_seconds=args.cumulative_cpu,wall_seconds=args.cumulative_wall,disk_MiB=args.disk_mib),
        attempt=dict(CPU_seconds=args.attempt_cpu,wall_seconds=args.attempt_wall,disk_MiB=disk),
        fixed_work='Reuse native66 and CO33; unchanged calibration.jl once for original CH4 parent+32 probes',
        work_command_prefix=[runtime['executable'],str(Path(__file__).resolve()),'work'])
    print(write_new(args.plan,plan),flush=True)
def load_plan(args):
    require(sha(args.plan)==args.plan_sha256,'Exact new prospective PLAN')
    plan=read(args.plan);require(plan['schema']=='calibration-completion-plan-v1','Completion schema');pin(plan['producer']['driver'])
    for key in ('guard_receipt','wrapper','prospective_task_pool_declaration'):pin(plan['original'][key])
    c,current=validate_original(plan['original']['scorer_root'],plan['original']['source'],plan['original']['root'])
    require(all(current[k]==plan['original'][k] for k in current),'Original cache/source/runtime changed after PLAN')
    return plan,c
def check_new_rows(plan,c):
    item=plan['original']['cases'][1];spec=item['spec'];folder=Path(plan['output'])/'julia'/CASES[1]
    rows=[];version=plan['original']['runtime']['julia']['version'].removeprefix('julia version ')
    for index,order in enumerate(item['probe']['orders']):
        parent=read(pin(item['native'][index]));path=folder/('row'+str(index).zfill(3)+'.toml')
        row=c._admit_child_row(path,CASES[1],index,order,parent,spec['sha256'],item['request']['sha256'],plan['original']['numerical_sources'])
        require(row['verified'] and row['julia_version']==version,'Every new CH4 triple parity/floor/runtime');rows.append(ref(path))
    path=folder/'SUMMARY.toml';summary=tomllib.loads(raw(path,64*1024).decode('ascii'))
    require(summary['complete'] is True and summary['processed_rows']==33 and summary['case']==CASES[1]
        and summary['request_sha256']==item['request']['sha256'] and summary['model_sha256']==spec['sha256']
        and summary['source_sha256']==plan['original']['numerical_sources'],'CH4 summary identity')
    return rows,ref(path)
def work(args):
    plan,c=load_plan(args);runtime=plan['original']['runtime'];output=Path(plan['output'])
    require(not output.exists(),'One fresh completion attempt; no automatic retry');output.mkdir(parents=True);a,native=c.adapter()
    require(sys.version==runtime['python'] and str(Path(sys.executable).resolve())==runtime['resolved_executable'],'Pinned Python')
    require(sha(runtime['julia']['executable'])==runtime['julia']['executable_sha256'],'Pinned Julia binary')
    report=dict(schema='calibration-completion-work-v1',plan_sha256=args.plan_sha256,complete_numerical=False,
        original_producer=ORIGINAL_COMMIT,completion_producer=plan['producer'],jobs=[],failures=[],
        native_scores_computed=0,CO_scores_computed=0,physical_errors_read=False)
    cpu=plan['attempt']['CPU_seconds'];wall=plan['attempt']['wall_seconds'];child_CPU=0.
    def check(live=0.):
        c.require(time.process_time()+child_CPU+live<cpu-2.,'completion_cpu_limit')
        c.require(time.monotonic()-BOOT_WALL<wall-10.,'completion_wall_limit')
    try:
        check();item=plan['original']['cases'][1];job_record=dict(case=CASES[1],raw_CPU_seconds=None,stop='dispatched');report['jobs'].append(job_record)
        job=c.run_child(runtime['julia']['executable'],pin(item['model']),pin(item['request']),output/'julia'/CASES[1],
            cpu-time.process_time()-2.,wall-(time.monotonic()-BOOT_WALL)-10.,check,a)
        job_record.update(job);child_CPU=job['raw_CPU_seconds']
        require(job['returncode']==0 and job['stop']=='complete' and job['reaped'] is True,'Completion child failed: retain outputs/costs')
        rows,summary=check_new_rows(plan,c);check();report.update(complete_numerical=True,new_CH4_rows=rows,new_CH4_summary=summary)
    except Exception as exc:report['failures'].append(dict(type=type(exc).__name__,message=str(exc)[:2048]))
    report['observed_Python_plus_Julia_CPU_before_output']=time.process_time()+child_CPU
    report['cost_scope']='External guard owns all Git verification children, Python, owned Julia startup/JIT and teardown'
    report['external_guard_receipt_required']=True
    write_new(output/'COMPLETION_WORK.json',report);return 0 if report['complete_numerical'] else 2
def finalize(args):
    plan,c=load_plan(args);output=Path(plan['output']);path=output/'COMPLETION_WORK.json';work=read(path)
    require(work['plan_sha256']==args.plan_sha256 and work['completion_producer']==plan['producer'],'Work/PLAN producer join')
    command=plan['work_command_prefix']+['--plan',str(Path(args.plan).resolve()),'--plan-sha256',args.plan_sha256]
    receipt,newcost=receipt_cost(args.work_receipt,command,allow_failure=True)
    require(receipt['budget']['cpu_seconds']==plan['attempt']['CPU_seconds'] and receipt['budget']['wall_seconds']==plan['attempt']['wall_seconds'] and
        receipt['budget']['disk_mib']==plan['attempt']['disk_MiB'],'Actual work guard matched PLAN')
    oldcost=plan['original']['cost'];costs=dict(original=oldcost,completion=newcost,
        aggregate_charged_CPU_seconds=oldcost['charged_CPU_seconds']+newcost['charged_CPU_seconds'],
        aggregate_wall_seconds=oldcost['wall_seconds']+newcost['wall_seconds'],aggregate_output_bytes=oldcost['output_bytes']+newcost['output_bytes'])
    pool=plan['cumulative_pool'];budget_ok=costs['aggregate_charged_CPU_seconds']<=pool['CPU_seconds'] and costs['aggregate_wall_seconds']<=pool['wall_seconds'] and costs['aggregate_output_bytes']<=pool['disk_MiB']*1024**2
    result=dict(schema='molecular-generalization-old-calibration-combined-v1',complete=False,
        original_report=plan['original']['report'],original_SOURCE=plan['original']['source'],new_completion_PLAN=ref(args.plan),
        new_completion_work=ref(path),costs=costs,source_sha256=plan['original']['numerical_sources'],
        source_producers=dict(original=dict(commit=ORIGINAL_COMMIT,source=plan['original']['source']),completion=plan['producer']),
        original_standalone_PLAN_absent=True,runtime=plan['original']['runtime'],cases=[],probes=[],failures=[],
        physical_errors_read=False,finalizer_metadata_guard_receipt_in_shared_study_ledger=True)
    try:
        require(work['complete_numerical'] is True and receipt['exit_code']==0 and receipt.get('termination') is None and budget_ok,'Science and cumulative resource admission')
        newrows,summary=check_new_rows(plan,c);require(newrows==work['new_CH4_rows'],'New rows unchanged')
        for item in plan['original']['cases']:
            case=item['spec']['case'];native_rows=[read(pin(value)) for value in item['native']]
            refs=item['original_Julia'] if case==CASES[0] else newrows;checked=[]
            for index,order in enumerate(item['probe']['orders']):
                row=c._admit_child_row(pin(refs[index]),case,index,order,native_rows[index],item['spec']['sha256'],item['request']['sha256'],plan['original']['numerical_sources'])
                require(row['verified'],'No floor/mismatch may be dropped');checked.append(row)
            result['cases'].append(dict(case=case,native_artifacts=item['native'],Julia_artifacts=refs,native_producer=ORIGINAL_COMMIT,
                Julia_producer=ORIGINAL_COMMIT if case==CASES[0] else plan['producer']['commit'],rows=checked,complete=True))
            parent=native_rows[0]['score']['value'];scale=max(parent,c.ENERGY_FLOOR)
            for index in range(1,33):
                trial=native_rows[index]['score']['value'];result['probes'].append(dict(case=case,index=index,
                    pair=item['probe']['pairs'][index-1],order=item['probe']['orders'][index],parent_J=parent,native_J=trial,
                    normalized_delta=trial/scale-parent/scale,positive=trial>parent,verified=True,
                    objective_at_or_below_floor=trial<=c.ATOL,energy_at_or_below_floor=trial<=c.ENERGY_FLOOR))
        result['fit']=c.fit_temperature(result['probes']);result['complete']=True
    except Exception as exc:result['failures'].append(dict(type=type(exc).__name__,message=str(exc)[:2048]))
    digest=write_new(output/'CALIBRATION.json',result)
    if result['complete']:
        a,native=c.adapter();settings=dict(a.DEFAULT_SETTINGS,calibrated=True,calibration_sha256=digest,
            theta_start=result['fit']['theta_start'],theta_end=result['fit']['theta_end']);a.validate_settings(settings)
        write_new(output/'ANNEALING_SETTINGS.json',settings)
    return 0 if result['complete'] else 2
def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);phases=parser.add_subparsers(dest='phase',required=True);p=phases.add_parser('prepare')
    for name in ('completion-tree','completion-commit','scorer-root','original','original-source','original-receipt','original-wrapper','launch-declaration','output','plan'):p.add_argument('--'+name,required=True)
    p.add_argument('--a1-quiescent',action='store_true')
    for name in ('cumulative-cpu','cumulative-wall','attempt-cpu','attempt-wall'):p.add_argument('--'+name,type=float,required=True)
    p.add_argument('--disk-mib',type=int,required=True)
    for name in ('work','finalize'):
        p=phases.add_parser(name);p.add_argument('--plan',required=True);p.add_argument('--plan-sha256',required=True)
        if name=='finalize':p.add_argument('--work-receipt',required=True)
    args=parser.parse_args(argv);return prepare(args) if args.phase=='prepare' else work(args) if args.phase=='work' else finalize(args)
if __name__=='__main__':raise SystemExit(main())
