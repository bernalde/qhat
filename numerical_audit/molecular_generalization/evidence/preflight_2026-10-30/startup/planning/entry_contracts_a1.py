"""Affected entry/adapter checks only; unchanged numerical evidence is reused."""
import subprocess,sys
subprocess.run([sys.executable,'-m','unittest','discover','-s','numerical_audit/molecular_generalization/test','-p','test_annealing*.py','-v'],check=True)
subprocess.run(['/home/bernalde/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia','--startup-file=no','--history-file=no','--threads=1','--gcthreads=1','--heap-size-hint=256M','--project=analysis/julia_window','numerical_audit/molecular_generalization/test/test_entry_dispatch.jl'],check=True)
