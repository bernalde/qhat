using LinearAlgebra
include(joinpath(@__DIR__,"src","GeneralizationIndependent.jl"))
function main(args=ARGS)
    length(args)==4 || error("usage: verify.jl INPUT_MODEL.toml CANDIDATES.toml PRIMARY_DIRECTORY FRESH_OUTPUT_DIRECTORY")
    BLAS.set_num_threads(1)
    Threads.nthreads()==1 || error("one Julia numerical thread required")
    GeneralizationIndependent.run(args...)["status"]=="complete" ? 0 : 2
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(main())
end
