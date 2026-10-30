# OLD-only engineering diagnostic. Original worker/request bytes stay unchanged.
# Top-level inclusion avoids old-world global lookups and wrapper specialization.
using Profile, TOML, SHA
cpu_seconds() = ccall(:clock, Int64, ()) / 1e6
length(ARGS) == 3 || error("usage: WRAPPER ORIGINAL_ANNEAL ORIGINAL_REQUEST FRESH_OUTPUT")
source, request, output = abspath.(ARGS)
VERSION == v"1.12.6" && Threads.nthreads() == 1 || error("pinned single-thread Julia required")
!ispath(output) || error("fresh diagnostic output required")
mkdir(output)
include_began = cpu_seconds()
include(source)
include_CPU = cpu_seconds() - include_began
Profile.init(n=400_000, delay=0.05)
Profile.clear()
Base.cumulative_compile_timing(true)
compilation_began = Base.cumulative_compile_time_ns()
began = cpu_seconds()
wall_began = time_ns()
status = 2
diagnostic_error = ""
println("PROFILE ORIGINAL WORKER START"); flush(stdout)
Profile.@profile begin
    try
        global status = MolecularAnnealingWorker.run(request, joinpath(output,"CHECKPOINT.toml"), 55.0, 250.0)
    catch exception
        global diagnostic_error = first(sprint(showerror,exception),2000)
    end
end
worker_CPU = cpu_seconds() - began
worker_wall = (time_ns()-wall_began)/1e9
compiled = Base.cumulative_compile_time_ns() .- compilation_began
Base.cumulative_compile_timing(false)
println("PROFILE ORIGINAL WORKER STOP ", status); flush(stdout)
samples = Profile.fetch(limitwarn=false)
# Entire rendered reports are retained; the external guard bounds total output.
open(joinpath(output,"PROFILE.tree.txt"),"w") do io
    Profile.print(IOContext(io,:displaysize=>(80,160)),samples;format=:tree,C=true,maxdepth=40,mincount=3,recur=:flatc)
end
open(joinpath(output,"PROFILE.flat.txt"),"w") do io
    Profile.print(IOContext(io,:displaysize=>(80,160)),samples;format=:flat,C=true,sortedby=:count,mincount=3)
end
summary = Dict("schema"=>"molecular-old-annealing-profile-v2","scope"=>"Old-only engineering diagnostic, not a constructor retry or scientific budget comparison",
    "source_sha256"=>bytes2hex(sha256(read(source))),"request_sha256"=>bytes2hex(sha256(read(request))),
    "wrapper_sha256"=>bytes2hex(sha256(read(@__FILE__))),"inner_absolute_process_CPU_cap"=>55.0,"inner_wall_cap"=>250.0,
    "external_CPU_cap"=>70.0,"external_wall_cap"=>350.0,"include_CPU_seconds"=>include_CPU,
    "original_worker_CPU_seconds"=>worker_CPU,"original_worker_wall_seconds"=>worker_wall,
    "Julia_reported_compilation_seconds"=>compiled[1]/1e9,"Julia_reported_recompilation_seconds"=>compiled[2]/1e9,
    "worker_return_code"=>status,"diagnostic_error"=>diagnostic_error,"profile_IP_cap"=>400000,"profile_IPs_retained"=>length(samples),
    "profile_delay_seconds"=>0.05,"checkpoint_present"=>isfile(joinpath(output,"CHECKPOINT.toml")),
    "accounting"=>"External receipt owns startup, Profile overhead and report emission. Compilation elapsed metrics are not subtracted from CPU. No scientific budget or calibrated setting changes.")
open(joinpath(output,"PROFILE.toml"),"w") do io;TOML.print(io,summary;sorted=true);end
println("PROFILE REPORT COMPLETE");flush(stdout)
exit(status)
