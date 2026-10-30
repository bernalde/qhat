import sys,os
from pathlib import Path
sys.path.insert(0,"/tmp/qhat-molecular-generalization-2026-10-30/worktree/numerical_audit/molecular_generalization")
import common as k
head=k.committed_sources()
k.require(head=="4ba93256e11aa89038b91c907189c875771c6870","Exact declared calibration producer")
k.write_new(Path("/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30/preflight/calibration/calibration-a1.SOURCE.json"),dict(source_commit=head,sources=k.source_map(),old_cases_only=True,new_physical_errors_observed=False))
os.execv(sys.executable,[sys.executable,str(k.HERE/"calibration.py"),"--output","/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30/preflight/calibration/a1","--cpu-cap","238","--wall-cap","1190","--julia-executable","/home/bernalde/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia"])
