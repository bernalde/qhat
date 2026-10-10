#!/usr/bin/env python3
"""Fresh production replay of one frozen row with explicit archived import paths.

Reads frozen schedules and references but does not reuse an approximate output.
The comparison output is written to a new file chosen by the caller.
"""
import argparse,os,sys,json,time,resource
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[name]='1'
def require(ok,msg):
    if not ok:raise ValueError(msg)
def main():
    p=argparse.ArgumentParser();p.add_argument('--extension',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--qhat-root',type=Path);p.add_argument('--campaign',required=True);p.add_argument('--case',required=True);p.add_argument('--ordering',required=True);p.add_argument('--budget',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));require(not a.output.exists(),'Output already exists');ext=a.extension.resolve();qhat=(a.qhat_root or ext.parent/'qhat').resolve();require(qhat.name=='qhat','Archived package directory must be named qhat');os.environ['QHAT_ROOT']=str(qhat);os.environ.setdefault('NUMBA_CACHE_DIR',str(a.output.parent/'numba-cache'));os.environ.setdefault('MPLCONFIGDIR',str(a.output.parent/'mpl-cache'));sys.path.insert(0,str(ext));sys.path.insert(0,str(qhat.parent))
    import common as c
    require(Path(c.__file__).resolve()==ext/'common.py','Wrong extension common module')
    import numpy as np
    from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as ev
    def imported_modules():
        modules={}
        for name,mod in list(sys.modules.items()):
            if name=='qhat' or name.startswith('qhat.'):
                path=getattr(mod,'__file__',None)
                if path:
                    path=Path(path).resolve();require(path.is_relative_to(qhat),'Live-checkout import: '+str(path));modules[name]={'path':str(path.relative_to(qhat)),'sha256':c.sha(path)}
        return modules
    imported_modules()
    root=ext/a.campaign;dp=next(x for x in (root/'designs'/f'{a.case}.json',root/'designs'/f'{a.case}.design.json') if x.exists());d=json.loads(dp.read_text());method=d['methods'][a.ordering]
    rp=next(x for x in (root/'results'/f'{a.case}.{a.ordering}.b{a.budget}.json',root/'rows'/a.case/f'{a.ordering}.b{a.budget}.json') if x.exists());row=json.loads(rp.read_text());require(row['status']=='complete','Frozen row infeasible');require(row['extension_design_sha256']==c.sha(dp) and c.object_sha(method['schedule'])==row['schedule_sha256'],'Frozen row/design mismatch')
    raw,states=c.load_case(a.case);require(raw['extension_input_provenance']==d['extension_input_provenance'],'Archived input mismatch');cc=c.compile_case(raw,states);schedule=c.decode_schedule(method['schedule']);terms=[cc['term'](w,v) for w,v in schedule];r=row['steps'];started=time.monotonic()
    if row['formula_order']==1:actual=evolve(cc['initial'],terms,1.,r)
    else:
        # Unmerged palindrome: independently count equivalent emitted rotations
        # in the recorded budget; numerical replay need not merge to be valid.
        actual=cc['initial'].copy();work=np.empty_like(actual)
        for _ in range(r):
            for t in terms+list(reversed(terms)):actual,work=ev.apply_one_exponential(actual,work,t,1/(2*r),2**30)
    metrics=physical_metrics(cc['exact'],actual,cc['number_mask'],cc['spin_mask']);sp=root/row['state_file']
    if not sp.exists():sp=rp.parent/row['state_file']
    require(c.sha(sp)==row['state_sha256'],'Frozen output digest changed')
    with np.load(sp) as z:distance=float(np.linalg.norm(actual-z[row.get('state_key','actual')]))
    require(distance<2e-11,'Fresh production replay mismatch')
    result={'status':'passed','case':a.case,'ordering':a.ordering,'budget_units':a.budget,'steps':r,'state_distance':distance,'metrics':metrics,'seconds':time.monotonic()-started,'imported_qhat_modules':imported_modules(),'extension_design_sha256':c.sha(dp),'state_sha256':c.sha(sp),'scope':'fresh propagation of one frozen schedule; archive code imports explicitly checked'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('status','case','ordering','state_distance','seconds')}))
if __name__=='__main__':main()
