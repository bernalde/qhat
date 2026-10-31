using LinearAlgebra
include(joinpath(@__DIR__,"src","GeneralizationPhysical.jl"))
function main(args=ARGS)
    length(args)==3 || error("usage: physical.jl INPUT_MODEL.toml CANDIDATES.toml FRESH_OUTPUT_DIRECTORY")
    BLAS.set_num_threads(1)
    Threads.nthreads()==1 || error("one Julia numerical thread required")
    GeneralizationPhysical.run(args...)["status"]=="complete" ? 0 : 2
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(main())
end
