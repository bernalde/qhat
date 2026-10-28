using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
include("src/ResponseCore.jl")
include("src/StructureDP.jl")
include("src/StudyIO.jl")
include("src/ArchivedWindows.jl")
using .StudyIO, .ArchivedWindows

function measured(f)
    start=cpu_seconds(); value=@timed f()
    (result=value.value,cpu_seconds=cpu_seconds()-start,wall_seconds=value.time,
     allocated_bytes=value.bytes,gc_wall_seconds=value.gctime)
end

function one_case(path;multiplier=1,repetitions=3,deadline=Inf,on_sample=(stage,value)->nothing,
                  methods=("dp","enumeration"))
    check=()->(cpu_seconds()<deadline || error("comparison CPU deadline"))
    parsed=measured(()->load_archive(path;weight_multiplier=multiplier));data=parsed.result
    function run(method)
        GC.gc()
        sample=measured(()->method=="dp" ? solve_archive_dp(data;guard=check) : enumeration(data;guard=check))
        result=sample.result
        if result.complete
            result.validation.valid || error("archive objective/witness mismatch")
            if method=="enumeration"
                result.expected_minimizers_match || error("archive minimizer count mismatch")
                result.order==data.expected_window || error("archive canonical tied optimum mismatch")
                result.evaluated==factorial(data.width) && result.coverage_complete || error("enumeration coverage mismatch")
            end
        elseif method=="enumeration"
            error("incomplete fresh enumeration")
        end
        Dict("method"=>method,"cpu_seconds"=>sample.cpu_seconds,"wall_seconds"=>sample.wall_seconds,
          "allocated_bytes"=>sample.allocated_bytes,"gc_wall_seconds"=>sample.gc_wall_seconds,
          "solver"=>result)
    end
    # First calls share this process and some JIT helpers. They are not a fair
    # independent cold-start comparison. Warm rows provide method comparisons.
    cold=Dict{String,Any}[]
    for method in methods
        sample=run(method);push!(cold,sample);on_sample("cold_"*method,sample)
    end
    samples=Dict{String,Any}[]
    for repetition in 1:repetitions
        for method in (isodd(repetition) ? ("dp","enumeration") : ("enumeration","dp"))
            sample=run(method);sample["repetition"]=repetition;push!(samples,sample)
            on_sample("warm_"*string(repetition)*"_"*method,sample)
        end
    end
    Dict("case"=>data.case,"input_sha256"=>bytes2hex(sha256(read(path))),
      "metadata"=>data.metadata,"objective_denominator"=>data.objective_denominator,
      "expected_integer_minimum"=>data.expected_minimum,
      "expected_minimizers"=>data.expected_minimizers,"original_window"=>data.window,
      "warm_integer_score"=>archived_integer_score(data,data.window),
      "width"=>data.width,"parse_cpu_seconds"=>parsed.cpu_seconds,
      "first_calls"=>cold,"warm_samples"=>samples,
      "first_call_context"=>"One shared process; later methods/cases may reuse JIT helpers. Whole-process startup is in the external receipt.",
      "scope"=>"Same archived first-order JD integer objective; original outside order fixed; no cubic or physical-error optimality claim")
end

function main(args)
    length(args)>=3 || error("usage: archive_study.jl preflight|main OUTPUT.toml INPUT.toml...")
    mode,output=args[1:2];paths=args[3:end]
    mode in ("preflight","preflight_enumeration","main") || error("unknown mode")
    preflight=mode!="main"
    preflight && length(paths)!=1 && error("one preflight case required")
    began=cpu_seconds();wall=time();deadline=began+(preflight ? 14.0 : 135.0)
    record=Dict{String,Any}("mode"=>mode,"julia_version"=>string(VERSION),
      "threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
      "source_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in
        [@__FILE__;[joinpath(@__DIR__,"src",name) for name in ("ArchivedWindows.jl","StructureDP.jl","StudyIO.jl","ResponseCore.jl")]]))
    cases=Dict{String,Any}[]
    try
        for path in paths
            callback=(stage,value)->begin
                checkpoint=output*"."*splitext(basename(path))[1]*"."*stage*".toml"
                write_result(checkpoint,value)
                println(stage," ",basename(path)," ",value["solver"].status," CPU=",value["cpu_seconds"])
                flush(stdout)
            end
            push!(cases,one_case(path;multiplier=preflight ? 7 : 1,
              repetitions=preflight ? 0 : 3,deadline=deadline,on_sample=callback,
              methods=mode=="preflight_enumeration" ? ("enumeration",) : ("dp","enumeration")))
        end
        record["status"]="complete"
    catch exception
        record["status"]="incomplete";record["reason"]=sprint(showerror,exception)
    end
    record["cases"]=cases;record["process_cpu_seconds"]=cpu_seconds()-began
    record["process_wall_seconds"]=time()-wall
    write_result(output,record)
    println(record["status"]," ",output)
    record["status"]=="complete" || exit(2)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
