#!/usr/bin/env python3
"""Freeze all plans, propagate, and independently verify review Strang controls."""
import os,resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,hashlib,json,sys,time,types,shutil,gc
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'cache'/'numba'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'cache'/'mpl'))
BUDGETS=(2,5,10,20,50,100)
CASES=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]+[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
EXPECTED_MANIFEST='0f42cf3dcf15982ca1181bf0d8674d4cd7b15bc0ba7352b35ad01f7481f9bb13'
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
def objsha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def sources():return {n:sha(ROOT/n) for n in ('run.py','independent_core.py','PROTOCOL.md')}
def arch(path):
    path=path.resolve();require(sha(path/'artifact_manifest.json')==EXPECTED_MANIFEST,'Wrong original archive manifest');manifest=json.loads((path/'artifact_manifest.json').read_text())['files'];return path,manifest
def check_file(path,base,manifest):
    name=str(path.relative_to(base));require(name in manifest and sha(path)==manifest[name]['sha256'],'Archive member digest changed '+name)
def setup_archive(path):
    base,manifest=arch(path);qhat=base/'qhat';ext=base/'extension';check_file(ext/'common.py',base,manifest)
    require('qhat' not in sys.modules,'qhat imported before isolation')
    package=types.ModuleType('qhat');package.__path__=[str(qhat)];package.__package__='qhat';sys.modules['qhat']=package
    os.environ['QHAT_ROOT']=str(qhat);sys.path.insert(0,str(ext));import common
    require(Path(common.__file__).resolve()==ext/'common.py','Wrong common import')
    def imported():
        found={}
        require(list(sys.modules['qhat'].__path__)==[str(qhat)],'Namespace path merged')
        for name,mod in list(sys.modules.items()):
            if name.startswith('qhat.'):
                file=getattr(mod,'__file__',None)
                if file:
                    p=Path(file).resolve();require(p.is_relative_to(qhat),'Foreign QHAT import');check_file(p,base,manifest);found[name]={'path':str(p.relative_to(base)),'sha256':sha(p)}
                else:
                    require(all(Path(p).resolve().is_relative_to(qhat) for p in getattr(mod,'__path__',[])),'Foreign namespace path')
        return found
    return base,manifest,common,imported

def freeze(archive):
    import numpy as np
    base,manifest,c,imported=setup_archive(archive)
    require(not (ROOT/'ALL_PLANS_FROZEN.json').exists(),'Already frozen; use run/verify instead')
    designs={};provenance={}
    for case in CASES:
        dp=base/'extension/optimized/designs'/f'{case}.design.json';check_file(dp,base,manifest);d=json.loads(dp.read_text());raw,states=c.load_case(case)
        require(raw['extension_input_provenance']==d['extension_input_provenance'],'Reference/input provenance mismatch')
        (ROOT/'inputs').mkdir(exist_ok=True);sp=ROOT/'inputs'/f'{case}.npz';np.savez_compressed(sp,**states)
        methods={}
        for key in ('class_refined_response','class_pauli_JD','class_pauli_JD_budget'):
            old=d['methods'][key];units=[];offset=0
            for length in old['unit_lengths']:units.append(old['schedule'][offset:offset+length]);offset+=length
            require(offset==len(old['schedule']),'Bad unit lengths')
            for orientation in ('forward','reverse'):
                chosen=units if orientation=='forward' else list(reversed(units));schedule=[pair for unit in chosen for pair in unit];name=key+'__'+orientation
                methods[name]={'schedule':schedule,'schedule_sha256':objsha(schedule),'unit_lengths':[len(u) for u in chosen],'source_method':key,'source_schedule_sha256':old['schedule_sha256'],'orientation':orientation,'formula_order':2,'reused':key=='class_refined_response' and orientation=='forward','objective':old['objective']}
        plan={k:raw[k] for k in ('case','n_qubits','electrons','denominator','coefficients','affine_basis_indices','tensor_sha256','state','blocks')}
        plan.update(methods=methods,panel=d['panel'],T=1.,input_state_sha256=sha(sp),original_design_sha256=sha(dp),input_provenance=d['extension_input_provenance'],optimizer_reports=d['optimizer_reports'])
        out=ROOT/'designs'/f'{case}.json';write(out,plan);designs[case]=sha(out)
        for b in BUDGETS:
            oldrp=base/'extension/optimized/rows'/case/f'strang_class_refined_response.b{b}.json';check_file(oldrp,base,manifest);r=json.loads(oldrp.read_text());oldsp=base/'extension/optimized'/r['state_file'];check_file(oldsp,base,manifest)
            require(r['state_sha256']==sha(oldsp) and r['extension_design_sha256']==sha(dp),'Baseline provenance changed')
            name='class_refined_response__forward';newsp=ROOT/'states'/case/f'{name}.b{b}.npz';newsp.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(oldsp,newsp)
            row={**r,'ordering':name,'source_original_row':str(oldrp.relative_to(base)),'source_original_row_sha256':sha(oldrp),'source_original_state':str(oldsp.relative_to(base)),'source_original_state_sha256':sha(oldsp),'state_file':str(newsp.relative_to(ROOT)),'extension_design_sha256':sha(out),'reused':True,'orientation':'forward','timing_scope':'historical optimized row interval; not new propagation'}
            write(ROOT/'rows'/case/f'{name}.b{b}.json',row)
        provenance[case]={'original_design':str(dp.relative_to(base)),'sha256':sha(dp),'original_inputs':d['extension_input_provenance']}
    write(ROOT/'ALL_PLANS_FROZEN.json',{'source_sha256':sources(),'design_sha256':designs,'archive_manifest_sha256':EXPECTED_MANIFEST,'original_archive_sha256':'d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43','provenance':provenance,'expected_new_states':420,'reused_states':84,'total_rows':504,'budgets':BUDGETS,'frozen_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'imported_qhat_modules':imported()})
    print('FROZEN all 14 designs, 420 new plans, 84 reused baseline states',flush=True)

def verified_design(case):
    m=json.loads((ROOT/'ALL_PLANS_FROZEN.json').read_text());require(m['source_sha256']==sources(),'Frozen source changed');p=ROOT/'designs'/f'{case}.json';require(sha(p)==m['design_sha256'][case],'Frozen design changed');d=json.loads(p.read_text());require(sha(ROOT/'inputs'/f'{case}.npz')==d['input_state_sha256'],'Frozen inputs changed');return d

def pulses(lengths,r):
    offsets=[0]
    for length in lengths:offsets.append(offsets[-1]+length)
    units=[list(range(offsets[i],offsets[i+1])) for i in range(len(lengths))]
    yield from ((i,1) for i in units[0])
    for step in range(r):
        for unit in units[1:-1]:yield from ((i,1) for i in unit)
        yield from ((i,2) for i in units[-1])
        for unit in reversed(units[1:-1]):yield from ((i,1) for i in unit)
        yield from ((i,1 if step==r-1 else 2) for i in units[0])

def run(archive,cases):
    import numpy as np
    base,manifest,c,imported=setup_archive(archive)
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    from qhat.numerical_audit.run_frozen_audit import physical_metrics
    imports=imported();started=time.monotonic()
    for case in cases:
        d=verified_design(case)
        with np.load(ROOT/'inputs'/f'{case}.npz') as z:arrays={k:z[k].copy() for k in ('initial','exact')}
        cc=c.compile_case(d,arrays)
        for name,method in d['methods'].items():
            if method['reused']:continue
            schedule=c.decode_schedule(method['schedule']);U=len(schedule);lengths=method['unit_lengths'];terms=[cc['term'](w,v) for w,v in schedule]
            require(dict(schedule)==c.block(d['coefficients']) and len(schedule)==len(c.block(d['coefficients'])),'Wrong target')
            for b in BUDGETS:
                rp=ROOT/'rows'/case/f'{name}.b{b}.json';sp=ROOT/'states'/case/f'{name}.b{b}.npz'
                if rp.exists():
                    old=json.loads(rp.read_text());require(old['state_sha256']==sha(sp),'Existing state mismatch');continue
                B=b*U;slope=2*U-lengths[0]-lengths[-1];r=(B-lengths[0])//slope;require(r>=1,'Infeasible new plan');expected=r*slope+lengths[0];require(expected<=B and (r+1)*slope+lengths[0]>B,'Not largest feasible r')
                pp=list(pulses(lengths,r));require(len(pp)==expected,'Pulse count');durations=[0]*U
                for i,mult in pp:durations[i]+=mult
                require(durations==[2*r]*U,'Per-word time mismatch')
                before=time.monotonic();actual=cc['initial'].copy();work=np.empty_like(actual)
                for i,mult in pp:actual,work=ev.apply_one_exponential(actual,work,terms[i],mult/(2*r),2**30)
                seconds=time.monotonic()-before;sp.parent.mkdir(parents=True,exist_ok=True);tmp=sp.with_suffix('.partial.npz');np.savez_compressed(tmp,actual=actual);os.replace(tmp,sp)
                row={'case':case,'panel':d['panel'],'ordering':name,'orientation':method['orientation'],'T':1.,'formula_order':2,'budget_units':b,'rotation_budget':B,'steps':r,'status':'complete','pauli_rotations':expected,'seconds':seconds,'timing_scope':'production propagation only, after pulse assembly; includes first-call JIT if triggered','reused':False,'state_file':str(sp.relative_to(ROOT)),'state_key':'actual','state_sha256':sha(sp),'extension_design_sha256':sha(ROOT/'designs'/f'{case}.json'),'schedule_sha256':method['schedule_sha256'],'surrogate_J':method['objective'],**physical_metrics(cc['exact'],actual,cc['number_mask'],cc['spin_mask'])}
                write(rp,row);print('EVOLVE',case,name,b,'E',row['one_minus_overlap'],'seconds',seconds,flush=True)
        gc.collect()
    rows=[json.loads(p.read_text()) for p in sorted((ROOT/'rows').rglob('*.json'))];write(ROOT/'results.json',rows)
    if len(rows)==504:write(ROOT/'completion.json',{'rows':len(rows),'new_evolutions':sum(not r['reused'] for r in rows),'reused_states':sum(r['reused'] for r in rows),'source_sha256':sources(),'seconds':time.monotonic()-started,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'imported_qhat_modules':imported()})

def verify(cases):
    import numpy as np
    import independent_core as ind
    require(not any(n=='qhat' or n.startswith('qhat.') for n in sys.modules),'Independent checker imported QHAT')
    started=time.monotonic();allreports=[]
    for case in cases:
        d=verified_design(case);model=ind.OccupationAction(d)
        with np.load(ROOT/'inputs'/f'{case}.npz') as z:initial=z['initial'].copy();exact=z['exact'].copy()
        fresh,nnz=model.reference(initial);reference_distance=float(np.linalg.norm(fresh-exact));require(reference_distance<3e-11,'Independent exact reference mismatch')
        require(np.count_nonzero(initial)==1 and int(model.basis[np.argmax(abs(initial))])==(1<<int(d['electrons']))-1,'Not frozen HF')
        records=[]
        for name,method in d['methods'].items():
            require(ind.schedule_sum(method['schedule'])==model.total,'Target sum mismatch')
            schedule=method['schedule'];offset=0
            for length in method['unit_lengths']:
                unit=schedule[offset:offset+length];offset+=length
                if 'refined' in name:require(not ind.number_commutator(unit),'Conserving unit violates N')
            for b in BUDGETS:
                rp=ROOT/'rows'/case/f'{name}.b{b}.json';row=json.loads(rp.read_text());sp=ROOT/row['state_file'];require(sha(sp)==row['state_sha256'],'State hash mismatch');out=ROOT/'verification'/case/f'{name}.b{b}.json'
                if out.exists():
                    old=json.loads(out.read_text());require(old['row_sha256']==sha(rp),'Verification row changed');records.append(old);continue
                with np.load(sp) as z:actual=z['actual'].copy()
                before=time.monotonic();replay=model.replay(initial,schedule,row['steps'],2);delta=float(np.linalg.norm(replay-actual));seconds=time.monotonic()-before;require(delta<3e-11,'Unmerged independent replay mismatch')
                metrics=model.metrics(fresh,actual,d['electrons']);differences={k:abs(float(row[k])-v) for k,v in metrics.items() if k in row};require(max(differences.values(),default=0)<3e-11,'Independent metric mismatch');require(abs(metrics['norm']-1)<2e-9,'State normalization drift')
                record={'case':case,'ordering':name,'budget_units':b,'state_distance':delta,'reference_distance':reference_distance,'metric_differences':differences,'independent_metrics':metrics,'seconds':seconds,'reused':row['reused'],'row_sha256':sha(rp),'state_sha256':sha(sp),'source_sha256':sources()};write(out,record);records.append(record)
                print('VERIFY',case,name,b,delta,'seconds',seconds,flush=True)
        report={'case':case,'rows':len(records),'replays':len(records),'new_state_replays':sum(not r['reused'] for r in records),'reference_distance':reference_distance,'reference_nnz':nnz,'max_state_distance':max(r['state_distance'] for r in records),'max_metric_difference':max(max(r['metric_differences'].values(),default=0) for r in records),'records':records};write(ROOT/'verification'/f'{case}.json',report);allreports.append(report);gc.collect()
    reports=[json.loads((ROOT/'verification'/f'{case}.json').read_text()) for case in CASES if (ROOT/'verification'/f'{case}.json').exists()]
    if len(reports)==14:write(ROOT/'verification.json',{'status':'passed','cases':14,'replays':sum(r['replays'] for r in reports),'new_state_replays':sum(r['new_state_replays'] for r in reports),'max_state_distance':max(r['max_state_distance'] for r in reports),'max_metric_difference':max(r['max_metric_difference'] for r in reports),'max_reference_distance':max(r['reference_distance'] for r in reports),'seconds':time.monotonic()-started,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'All 420 new states and 84 reused baseline states replayed with independently derived occupation-basis Pauli action and unmerged Strang palindromes; all 14 reference states recomputed independently.'})

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','run','verify']);p.add_argument('--archive-root',type=Path);p.add_argument('--cases',nargs='*',default=CASES);a=p.parse_args()
    require(set(a.cases)<=set(CASES),'Unknown cases')
    if a.command in ('freeze','run'):require(a.archive_root is not None,'Explicit --archive-root required');freeze(a.archive_root) if a.command=='freeze' else run(a.archive_root,a.cases)
    else:verify(a.cases)
if __name__=='__main__':main()
