"""One diagnostic; unchanged scientific worker, only stderr retention changes."""
from pathlib import Path
import hashlib,json,sys
study=Path('/tmp/qhat-normal-window-2026-10-28/worktree/numerical_audit/prospective_window_ordering')
sys.path.insert(0,str(study))
import constructors as c
root=Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/prospective_window_ordering_2026-10-28/a1')
out=root/'diagnostic-window-stderr-a2';out.mkdir()
model=root/'preflight-data-a2/inputs/LiH_r1p487_631g_14m4e_preflight_7over5/model.json'
c.k.bootstrap()
from qhat.analysis import window_ordering as w
original=w.subprocess.Popen
with (out/'worker.stderr').open('x') as stream:
 def retain_stderr(argv,*args,**kwargs):
  if len(argv)>1 and Path(argv[0]).name=='julia':kwargs['stderr']=stream
  return original(argv,*args,**kwargs)
 w.subprocess.Popen=retain_stderr
 try:
  result=c.construct(model,c.k.sha(model),'window_magnitude',out/'constructor.json',repetition=1,source_freeze=root/'data-a3/INPUT_PROTOCOL_FROZEN.json',julia_executable='/home/bernalde/.julia/juliaup/julia-1.10.11+0.x64.linux.gnu/bin/julia')
 finally:w.subprocess.Popen=original
print('DIAGNOSTIC_STATUS',result['constructor']['stats'].get('window_status'))
print('UNCHANGED_WORKER_STDERR', (out/'worker.stderr').read_text()[:4096])
