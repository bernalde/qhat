module CostMatchedObjective

# The frozen independent objective implementation fixes training T=1 and 30U.
# The owning OS guard accounts for startup/import/JIT and all teardown costs.
using LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "numerical_audit", "prospective_window_ordering", "objective_check.jl"))

function main(args=ARGS)
    length(args) == 3 || error("usage: objective.jl INPUT_MODEL.toml CANDIDATES.toml FRESH_OUTPUT_DIRECTORY")
    input, candidates, output = args
    ispath(output) && error("fresh objective output directory required")
    BLAS.set_num_threads(1)
    Threads.nthreads() == 1 || error("one Julia numerical thread required")
    mkdir(output)
    result = ObjectiveCheck.run_check(input, candidates, joinpath(output, "OBJECTIVE.toml");
        cpu_limit=120.0, wall_limit=600.0)
    result["complete"] ? 0 : 2
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(CostMatchedObjective.main())
end
