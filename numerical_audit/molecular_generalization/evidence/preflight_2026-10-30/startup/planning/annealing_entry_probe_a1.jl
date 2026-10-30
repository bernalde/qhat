# Old-only entry inference experiment; numerical implementation is unchanged.
include("/tmp/qhat-molecular-generalization-2026-10-30/worktree/numerical_audit/molecular_generalization/anneal.jl")
length(ARGS)==4 || error("four literal old-request worker arguments required")
exit(Base.invokelatest(MolecularAnnealingWorker.run,ARGS[1],ARGS[2],parse(Float64,ARGS[3]),parse(Float64,ARGS[4])))
