"""Resume only missing CH4 calibration rows after two closed attempts.

prepare --prior-plan ... --prior-work ... --prior-receipt ... --resources ...
        --producer-commit ... --output FRESH_DIR --plan FRESH_PLAN.json
work --plan ... --plan-sha256 ...
finalize --plan ... --plan-sha256 ... --work-receipt ...

All native66/CO33 are inherited from 4ba; complete CH4 rows are inherited from
023. One new child scores the frozen missing indices, with unchanged kernels,
full33 request, comparisons, floors and fitting. No chemistry/physical errors.
The coordinator owns all dispatch and the prospective verification RAM profile.
"""
from __future__ import annotations
import time
BOOT_WALL=time.monotonic()
import argparse
import hashlib
import importlib.util
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tomllib
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
MEMORY=dict(MemoryHigh='640M',MemoryMax='768M',MemorySwapMax='0')
class Refusal(ValueError):pass
def require(ok,message):
    if not ok:raise Refusal(message)
def file_sha(path):
    with Path(path).open('rb') as stream:raw=stream.read(32*1024**2+1)
    require(len(raw)<=32*1024**2,'Artifact byte cap');return hashlib.sha256(raw).hexdigest()
def helper(prior_path):
    import json
    with Path(prior_path).open('rb') as stream:content=stream.read(32*1024**2+1)
    require(len(content)<=32*1024**2,'Prior PLAN byte cap');prior=json.loads(content)
    path=HERE/'calibration_completion.py'
    require(file_sha(path)==prior['producer']['driver']['sha256'],'Unchanged 023 completion helper')
    spec=importlib.util.spec_from_file_location('_resume_prior_completion',path)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
    p,c=m.load_plan(SimpleNamespace(plan=str(Path(prior_path).resolve()),plan_sha256=m.sha(prior_path)))
    return m,p,c
def producer(commit,m):
    require(re.fullmatch(r'[0-9a-f]{40}',commit),'Committed resume producer identity')
    head=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    require(head==commit,'Exact resume HEAD')
    files={}
    for path in (Path(__file__).resolve(),HERE/'calibration_resume.jl'):
        name=str(path.relative_to(ROOT));blob=subprocess.check_output(['git','-C',str(ROOT),'show',commit+':'+name])
        require(hashlib.sha256(blob).hexdigest()==m.sha(path),'Commit both resume entrypoints')
        files[name]=m.sha(path)
    return dict(commit=commit,files=files,python=m.ref(__file__),julia=m.ref(HERE/'calibration_resume.jl'))
def partition_rows(rows):
    known={}
    for row in rows:
        index=row['index'];require(type(index) is int and 0<=index<=32,'Original row index')
        require(index not in known,'Duplicate existing row is not a quality choice')
        require(row.get('verified') is True,'A mismatch/floor/refusal cannot become missing work')
        known[index]=row
    return [known[i] for i in sorted(known)],[i for i in range(33) if i not in known]
def combine_rows(reused,new):
    rows,missing=partition_rows(list(reused)+list(new))
    require(not missing and len(rows)==33,'Exact completed CH4 row roster');return rows
def resource_plan(costs,declaration):
    require(declaration.get('schema')=='calibration-resume-resources-v1' and declaration.get('declared_before_use') is True,
            'Prospective resume resource declaration')
    pool=declaration['cumulative_pool'];attempt=declaration['attempt']
    require(declaration['memory_profile']==MEMORY,'Verification-only 640/768MiB/no-swap profile')
    require(pool['CPU_seconds']==840 and pool['wall_seconds']==6000 and
            attempt['CPU_seconds']==360 and attempt['wall_seconds']==1800,'Predeclared cumulative/attempt CPU and wall')
    for values in (pool,attempt):
        require(set(values)=={'CPU_seconds','wall_seconds','disk_MiB'} and
                all(type(x) in (int,float) and math.isfinite(x) and x>0 for x in values.values()),'Finite resource pools')
    require(attempt['disk_MiB']<=pool['disk_MiB']<=1024,'Bounded disk')
    for cost in costs:
        require(all(type(cost.get(k)) in (int,float) and math.isfinite(cost[k]) and cost[k]>=0
                    for k in ('charged_CPU_seconds','wall_seconds','output_bytes')),'All prior costs known and retained')
    require(sum(c['charged_CPU_seconds'] for c in costs)+attempt['CPU_seconds']<=pool['CPU_seconds'] and
            sum(c['wall_seconds'] for c in costs)+attempt['wall_seconds']<=pool['wall_seconds'] and
            sum(c['output_bytes'] for c in costs)+attempt['disk_MiB']*1024**2<=pool['disk_MiB']*1024**2,
            'No reset: prospective attempt fits remaining amended pool')
    return pool,attempt
def read_partial(m,p,c):
    item=p['original']['cases'][1];spec=item['spec'];folder=Path(p['output'])/'julia'/'CH4_r1p09'
    version=p['original']['runtime']['julia']['version'].removeprefix('julia version ')
    rows=[]
    for path in sorted(folder.glob('row*.toml')):
        match=re.fullmatch(r'row([0-9]{3})\.toml',path.name);require(match is not None,'Unexpected prior row filename')
        index=int(match.group(1));require(0<=index<=32,'Unexpected prior row index')
        native=m.read(m.pin(item['native'][index]));order=item['probe']['orders'][index]
        checked=c._admit_child_row(path,'CH4_r1p09',index,order,native,spec['sha256'],item['request']['sha256'],p['original']['numerical_sources'])
        require(checked['verified'] and checked['julia_version']==version,'Prior row parity/floor/runtime refusal')
        rows.append(dict(index=index,verified=True,artifact=m.ref(path),producer=p['producer'],kind='completion-a1'))
    reused,missing=partition_rows(rows)
    require(reused and missing,'This resume requires an interrupted partial CH4 panel')
    return reused,missing
def prepare(args):
    m,p,c=helper(args.prior_plan);prior_work=m.read(args.prior_work)
    require(prior_work['schema']=='calibration-completion-work-v1' and prior_work['plan_sha256']==m.sha(args.prior_plan)
        and prior_work['completion_producer']==p['producer'] and prior_work['complete_numerical'] is False,'Prior PLAN/work source join')
    jobs=prior_work['jobs'];require(len(jobs)==1 and jobs[0]['case']=='CH4_r1p09' and jobs[0]['reaped'] is True
        and jobs[0]['stop']=='child_wall_limit' and jobs[0]['returncode']==-9,'Closed/reaped interrupted prior child')
    item=p['original']['cases'][1];runtime=p['original']['runtime'];root=Path(p['original']['scorer_root'])
    expected=[runtime['julia']['executable'],'--startup-file=no','--history-file=no','--threads=1','--gcthreads=1',
        '--heap-size-hint=32M','--project='+str(root/'analysis/julia_window'),
        str(root/'numerical_audit/molecular_generalization/calibration.jl'),
        str(m.pin(item['model'])),str(m.pin(item['request'])),str(Path(p['output'])/'julia'/'CH4_r1p09')]
    argv=jobs[0]['command'];require(isinstance(argv,list) and len(argv)==len(expected)+2 and argv[:-2]==expected,
        'Prior actual owned Julia command/input/output identity')
    require(all(math.isfinite(float(x)) and float(x)>0 for x in argv[-2:]),'Prior finite child caps')
    command=p['work_command_prefix']+['--plan',str(Path(args.prior_plan).resolve()),'--plan-sha256',m.sha(args.prior_plan)]
    receipt,prior_cost=m.receipt_cost(args.prior_receipt,command,allow_failure=True)
    require(receipt['exit_code']==2 and receipt.get('termination') is None,'Prior ended internally, not active')
    require(receipt['budget']['cpu_seconds']==p['attempt']['CPU_seconds'] and receipt['budget']['wall_seconds']==p['attempt']['wall_seconds'] and
        receipt['budget']['disk_mib']==p['attempt']['disk_MiB'],
            'Prior actual guard/PLAN allowance')
    reused,missing=read_partial(m,p,c);resources=m.read(args.resources)
    costs=[p['original']['cost'],prior_cost];pool,attempt=resource_plan(costs,resources)
    out=Path(args.output).resolve()
    require(not out.exists() and all(not out.is_relative_to(Path(root)) and not Path(root).is_relative_to(out)
            for root in (p['output'],p['original']['root'])),'Fresh disjoint resume output')
    ident=producer(args.producer_commit,m)
    plan=dict(schema='calibration-resume-plan-v1',producer=ident,prior_plan=m.ref(args.prior_plan),
        prior_work=m.ref(args.prior_work),prior_receipt=m.ref(args.prior_receipt),prior_costs=costs,
        resources=m.ref(args.resources),cumulative_pool=pool,attempt=attempt,
        guard_profile=dict(m.PROFILE,**MEMORY),reused_rows=reused,missing_indices=missing,output=str(out),
        original_source_commit=m.ORIGINAL_COMMIT,prior_completion_producer=p['producer'],
        work_command_prefix=[p['original']['runtime']['executable'],str(Path(__file__).resolve()),'work'],
        no_native_or_CO_rescoring=True,physical_errors_read=False)
    print(m.write_new(args.plan,plan),flush=True)
def load_plan(args):
    # Read the prospective parent reference before importing its pinned helper.
    import json
    with Path(args.plan).open('rb') as stream:content=stream.read(32*1024**2+1)
    require(len(content)<=32*1024**2,'Resume PLAN byte cap');q=json.loads(content);m,p,c=helper(q['prior_plan']['path']);q=m.read(args.plan)
    require(m.sha(args.plan)==args.plan_sha256 and q['schema']=='calibration-resume-plan-v1','Exact prospective resume PLAN')
    for key in ('prior_plan','prior_work','prior_receipt','resources'):m.pin(q[key])
    for name,digest in q['producer']['files'].items():require(m.sha(ROOT/name)==digest,'Resume source bytes changed')
    pool,attempt=resource_plan(q['prior_costs'],m.read(m.pin(q['resources'])))
    require(pool==q['cumulative_pool'] and attempt==q['attempt'] and q['guard_profile']==dict(m.PROFILE,**MEMORY),
            'Frozen PLAN/resource/profile join')
    reused,missing=read_partial(m,p,c);require(reused==q['reused_rows'] and missing==q['missing_indices'],'Prior cache changed after freeze')
    return q,m,p,c
def write_subset(path,q,m,p,a,plan_sha):
    item=p['original']['cases'][1]
    fields=dict(version=1,case='CH4_r1p09',plan_sha256=plan_sha,producer_commit=q['producer']['commit'],
        model_sha256=item['model']['sha256'],request_sha256=item['request']['sha256'],indices=q['missing_indices'],
        julia_version=p['original']['runtime']['julia']['version'].removeprefix('julia version '))
    lines=[key+'='+a._toml(value) for key,value in fields.items()]
    lines.append('[resume_source_sha256]');lines.extend(a._toml(k)+'='+a._toml(v) for k,v in q['producer']['files'].items())
    require(not path.exists(),'Fresh subset declaration')
    with path.open('xb') as stream:stream.write(('\n'.join(lines)+'\n').encode('ascii'))
    return m.ref(path)
def run_child(q,m,p,a,subset,cpu,wall,check):
    require(math.isfinite(cpu) and cpu>.2 and math.isfinite(wall) and wall>1,'Remaining bounded subset allowance')
    runtime=p['original']['runtime'];item=p['original']['cases'][1];out=Path(q['output'])/'julia'/'CH4_r1p09'
    require(not out.exists(),'Fresh resume Julia output');out.parent.mkdir(parents=True,exist_ok=True)
    argv=[runtime['julia']['executable'],'--startup-file=no','--history-file=no','--threads=1','--gcthreads=1',
        '--heap-size-hint=32M','--project='+str(ROOT/'analysis/julia_window'),str(m.pin(q['producer']['julia'])),
        str(m.pin(item['model'])),str(m.pin(item['request'])),str(m.pin(subset)),str(out),repr(cpu),repr(wall)]
    env=os.environ.copy();env.update(JULIA_NUM_THREADS='1',JULIA_NUM_GC_THREADS='1',OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',JULIA_LOAD_PATH='@:@stdlib',JULIA_PKG_PRECOMPILE_AUTO='0',JULIA_IMAGE_THREADS='1')
    log=out.parent/'CH4_r1p09.stderr.log';require(not log.exists(),'Fresh child diagnostic')
    with log.open('xb') as stream:process=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=stream,env=env,start_new_session=True)
    began=time.monotonic();reaped=False;stop=None;result=dict(command=argv,raw_CPU_seconds=None,reaped=False,stderr_path=str(log))
    try:
        while True:
            pid,status,usage=os.wait4(process.pid,os.WNOHANG)
            if pid:reaped=True;break
            live=a._proc_cpu(process.pid)
            if live is None:stop='owned_child_cpu_unavailable'
            elif live>=cpu:stop='child_cpu_limit'
            elif time.monotonic()-began>=wall:stop='child_wall_limit'
            else:
                try:check(live)
                except Refusal as exc:stop=str(exc)
            if stop:break
            time.sleep(.005)
    finally:
        if not reaped:
            try:os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            _,status,usage=os.wait4(process.pid,0)
        process.returncode=os.waitstatus_to_exitcode(status)
        result.update(reaped=True,returncode=process.returncode,raw_CPU_seconds=usage.ru_utime+usage.ru_stime,
            peak_RSS_KiB=usage.ru_maxrss,wall_seconds=time.monotonic()-began,
            stop=stop or ('complete' if process.returncode==0 else 'child_refusal'),stderr=m.ref(log))
    return result
def check_rows(q,m,p,c,subset):
    item=p['original']['cases'][1];folder=Path(q['output'])/'julia'/'CH4_r1p09';rows=[]
    for index in q['missing_indices']:
        path=folder/('row'+str(index).zfill(3)+'.toml');native=m.read(m.pin(item['native'][index]))
        checked=c._admit_child_row(path,'CH4_r1p09',index,item['probe']['orders'][index],native,item['spec']['sha256'],
            item['request']['sha256'],p['original']['numerical_sources'])
        require(checked['verified'] and checked['resume_plan_sha256']==q['_plan_sha256'] and
            checked['resume_subset_sha256']==subset['sha256'] and checked['resume_source_sha256']==q['producer']['files'] and
            checked['resume_producer_commit']==q['producer']['commit'] and
            checked['julia_version']==p['original']['runtime']['julia']['version'].removeprefix('julia version '),
            'New subset row provenance/parity/floor/runtime')
        rows.append(dict(index=index,verified=True,artifact=m.ref(path),producer=q['producer'],kind='resume'))
    actual={path.name for path in folder.glob('row*.toml')}
    require(actual=={'row'+str(i).zfill(3)+'.toml' for i in q['missing_indices']},'Only frozen missing rows were scored')
    summary_path=folder/'SUMMARY.toml';summary=tomllib.loads(m.raw(summary_path,64*1024).decode('ascii'))
    require(summary['complete'] is True and summary['indices']==q['missing_indices'] and
        summary['resume_plan_sha256']==q['_plan_sha256'] and summary['resume_subset_sha256']==subset['sha256'] and
        summary['resume_source_sha256']==q['producer']['files'] and summary['request_sha256']==item['request']['sha256'] and
        summary['model_sha256']==item['spec']['sha256'],'Complete subset summary')
    return rows,m.ref(summary_path)
def work(args):
    q,m,p,c=load_plan(args);q['_plan_sha256']=args.plan_sha256;out=Path(q['output']);require(not out.exists(),'One fresh resume identity');out.mkdir(parents=True)
    runtime=p['original']['runtime'];require(sys.version==runtime['python'] and str(Path(sys.executable).resolve())==runtime['resolved_executable'],'Pinned Python')
    require(m.sha(runtime['julia']['executable'])==runtime['julia']['executable_sha256'],'Pinned Julia')
    report=dict(schema='calibration-resume-work-v1',plan_sha256=args.plan_sha256,producer=q['producer'],complete_numerical=False,
        reused_rows=q['reused_rows'],missing_indices=q['missing_indices'],jobs=[],failures=[],native_scores_computed=0,CO_scores_computed=0,
        physical_errors_read=False)
    cpu=q['attempt']['CPU_seconds'];wall=q['attempt']['wall_seconds'];child=0.
    def check(live=0.):
        require(time.process_time()+child+live<cpu-2.,'resume_cpu_limit')
        require(time.monotonic()-BOOT_WALL<wall-10.,'resume_wall_limit')
    try:
        a,native=c.adapter();check();subset=write_subset(out/'SUBSET.toml',q,m,p,a,args.plan_sha256);report['subset']=subset
        job_record=dict(raw_CPU_seconds=None,stop='dispatched');report['jobs'].append(job_record)
        job=run_child(q,m,p,a,subset,cpu-time.process_time()-2.,wall-(time.monotonic()-BOOT_WALL)-10.,check)
        job_record.update(job);child=job['raw_CPU_seconds']
        require(job['reaped'] and job['returncode']==0 and job['stop']=='complete','Subset child refusal: preserve all outputs and costs')
        rows,summary=check_rows(q,m,p,c,subset);combine_rows(q['reused_rows'],rows);check()
        report.update(complete_numerical=True,new_rows=rows,summary=summary)
    except Exception as exc:report['failures'].append(dict(type=type(exc).__name__,message=str(exc)[:2048]))
    report['Python_plus_owned_Julia_CPU_before_output']=time.process_time()+child;report['external_OS_receipt_required']=True
    m.write_new(out/'RESUME_WORK.json',report);return 0 if report['complete_numerical'] else 2
def new_receipt(m,path,command,profile):
    r=m.read(path);require(r['command']==command and isinstance(r.get('completed_utc'),str),'Closed resume receipt command')
    require(type(r['exit_code']) is int,'Integer terminal exit code')
    require(all(r['settings'].get(k)==v for k,v in profile.items()),'Prospective verification RAM/CPU profile')
    require(all(r['environment'].get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')),'One numerical thread')
    require(r.get('guard_error') is None and not r.get('cleanup_errors'),'Clean guard teardown')
    measurements=r.get('process_measurements',{});a=r['final_accounting'];values=[];upper=[]
    if a.get('process_CPU_timing_available') is True:
        require(all(type(measurements.get(k)) in (int,float) and math.isfinite(measurements[k]) and measurements[k]>=0
            for k in ('user_cpu_seconds','system_cpu_seconds')),'Terminal GNU time')
        value=measurements['user_cpu_seconds']+measurements['system_cpu_seconds'];values.append(value);upper.append(value+.011)
    footer=re.findall(r'^CPU time consumed:\s*(.+)$',m.raw(r['log']).decode('utf-8'),re.MULTILINE)
    if footer:
        value,pad=m.duration(footer[-1]);values.append(value);upper.append(value+pad)
    observed=a.get('CPU_seconds');require(observed is None or type(observed) in (int,float) and math.isfinite(observed) and observed>=0,'Finite CPU observation')
    require(values,'Unknown terminal CPU is not zero');require(type(r['elapsed']) in (int,float) and math.isfinite(r['elapsed']) and r['elapsed']>=0,'Actual wall')
    require(type(a.get('output_bytes')) is int and a['output_bytes']>=0,'Actual disk')
    return r,dict(raw_CPU_seconds=max(values),charged_CPU_seconds=math.ceil(max(upper+[observed or 0.])*1000)/1000,
        wall_seconds=r['elapsed'],output_bytes=a['output_bytes'],receipt=m.ref(path),log=m.ref(r['log']),failed_attempt=bool(r['exit_code'] or r.get('termination')))
def finalize(args):
    q,m,p,c=load_plan(args);q['_plan_sha256']=args.plan_sha256;out=Path(q['output']);work_path=out/'RESUME_WORK.json';w=m.read(work_path)
    require(w['plan_sha256']==args.plan_sha256 and w['producer']==q['producer'],'Resume work/producer join')
    command=q['work_command_prefix']+['--plan',str(Path(args.plan).resolve()),'--plan-sha256',args.plan_sha256]
    receipt,cost=new_receipt(m,args.work_receipt,command,q['guard_profile'])
    require(receipt['budget']['cpu_seconds']==q['attempt']['CPU_seconds'] and receipt['budget']['wall_seconds']==q['attempt']['wall_seconds'] and
        receipt['budget']['disk_mib']==q['attempt']['disk_MiB'],'Resume actual guard allowance')
    costs=q['prior_costs']+[cost];pool=q['cumulative_pool']
    totals={key:sum(row[key] for row in costs) for key in ('charged_CPU_seconds','wall_seconds','output_bytes')}
    result=dict(schema='molecular-generalization-old-calibration-resumed-v1',complete=False,cases=[],probes=[],failures=[],
        original_report=p['original']['report'],prior_plan=q['prior_plan'],prior_work=q['prior_work'],resume_plan=m.ref(args.plan),
        resume_work=m.ref(work_path),costs=costs,aggregate_costs=totals,runtime=p['original']['runtime'],
        source_sha256=p['original']['numerical_sources'],source_producers=dict(original=m.ORIGINAL_COMMIT,
            completion=p['producer'],resume=q['producer']),original_standalone_PLAN_absent=True,physical_errors_read=False,
        finalizer_metadata_receipt_in_shared_study_ledger=True)
    try:
        require(w['complete_numerical'] is True and receipt['exit_code']==0 and receipt.get('termination') is None,'Complete subset work/receipt')
        require(totals['charged_CPU_seconds']<=pool['CPU_seconds'] and totals['wall_seconds']<=pool['wall_seconds'] and
            totals['output_bytes']<=pool['disk_MiB']*1024**2,'Cumulative costs include both interrupted attempts')
        rows,summary=check_rows(q,m,p,c,w['subset']);require(rows==w['new_rows'],'Subset artifacts unchanged')
        ch4=combine_rows(q['reused_rows'],rows)
        for item in p['original']['cases']:
            case=item['spec']['case'];native=[m.read(m.pin(value)) for value in item['native']]
            origins=[dict(index=i,artifact=value,producer=dict(commit=m.ORIGINAL_COMMIT),kind='original-a1')
                for i,value in enumerate(item['original_Julia'])] if case=='CO_r1p13' else ch4
            checked=[]
            for origin in origins:
                i=origin['index'];row=c._admit_child_row(m.pin(origin['artifact']),case,i,item['probe']['orders'][i],native[i],
                    item['spec']['sha256'],item['request']['sha256'],p['original']['numerical_sources'])
                require(row['verified'],'No floor/mismatch may be excluded');checked.append(dict(witness=row,origin=origin))
            result['cases'].append(dict(case=case,rows=checked,native_artifacts=item['native'],native_producer=m.ORIGINAL_COMMIT,complete=True))
            parent=native[0]['score']['value'];scale=max(parent,c.ENERGY_FLOOR)
            for i in range(1,33):
                value=native[i]['score']['value'];result['probes'].append(dict(case=case,index=i,parent_J=parent,native_J=value,
                    pair=item['probe']['pairs'][i-1],order=item['probe']['orders'][i],normalized_delta=value/scale-parent/scale,
                    positive=value>parent,verified=True,objective_at_or_below_floor=value<=c.ATOL,energy_at_or_below_floor=value<=c.ENERGY_FLOOR))
        result['fit']=c.fit_temperature(result['probes']);result['complete']=True
    except Exception as exc:result['failures'].append(dict(type=type(exc).__name__,message=str(exc)[:2048]))
    digest=m.write_new(out/'CALIBRATION.json',result)
    if result['complete']:
        a,native=c.adapter();settings=dict(a.DEFAULT_SETTINGS,calibrated=True,calibration_sha256=digest,
            theta_start=result['fit']['theta_start'],theta_end=result['fit']['theta_end']);a.validate_settings(settings)
        m.write_new(out/'ANNEALING_SETTINGS.json',settings)
    return 0 if result['complete'] else 2
def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);phases=parser.add_subparsers(dest='phase',required=True);p=phases.add_parser('prepare')
    for name in ('prior-plan','prior-work','prior-receipt','resources','producer-commit','output','plan'):p.add_argument('--'+name,required=True)
    for name in ('work','finalize'):
        p=phases.add_parser(name);p.add_argument('--plan',required=True);p.add_argument('--plan-sha256',required=True)
        if name=='finalize':p.add_argument('--work-receipt',required=True)
    args=parser.parse_args(argv);return prepare(args) if args.phase=='prepare' else work(args) if args.phase=='work' else finalize(args)
if __name__=='__main__':raise SystemExit(main())
