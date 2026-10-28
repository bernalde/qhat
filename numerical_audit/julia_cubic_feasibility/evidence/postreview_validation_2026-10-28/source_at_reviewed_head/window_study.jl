using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
include("src/ResponseCore.jl")
include("src/StructureDP.jl")
include("src/StudyIO.jl")
include("src/WindowSearch.jl")
include("src/Diagnostics.jl")
using .ResponseCore, .StructureDP, .StudyIO, .WindowSearch, .Diagnostics

function timed(f)
    began=cpu_seconds(); value=@timed f()
    (result=value.value,cpu=cpu_seconds()-began,wall=value.time,
     allocated_bytes=value.bytes,gc_wall=value.gctime)
end

function main(args)
    length(args)==3 || error("usage: window_study.jl preflight|window|structure INPUT.toml OUTPUT.toml")
    mode,input,output=args
    mode in ("preflight","window","structure") || error("unknown mode")
    began=cpu_seconds(); wall=time(); stage="parse"
    record=Dict{String,Any}("mode"=>mode,"julia_version"=>string(VERSION),
        "threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in
          [@__FILE__; [joinpath(@__DIR__,"src",p) for p in
             ("ResponseCore.jl","StructureDP.jl","StudyIO.jl","WindowSearch.jl","Diagnostics.jl")]]))
    try
        loaded=timed(()->load_model(input;scale=mode=="preflight" ? 7//5 : 1//1))
        data=loaded.result
        merge!(record,Dict("case"=>data.case,"input_sha256"=>data.input_sha256,
          "parse_cpu_seconds"=>loaded.cpu,"words"=>data.words,"units"=>length(data.units)))
        deadline=cpu_seconds()+(mode=="structure" ? 120.0 : 32.0)
        check=()->(cpu_seconds()<deadline || error("declared method CPU allowance exhausted"))
        if mode=="structure"
            stage="structural audit"
            measured=timed(()->structural_audit(data.units,data.order,data.hf;guard=check))
            record["audit"]=measured.result;record["audit_cpu_seconds"]=measured.cpu
        else
            stage="window selection"
            selected=timed(()->select_window(data.units,data.order,data.hf,8;guard=check))
            record["selection_cpu_seconds"]=selected.cpu
            record["window_start"]=first(selected.result.positions)
            record["window_width"]=length(selected.result.positions)
            record["selection_active_pairs"]=selected.result.active_pairs
            stage="conditioned table construction"
            built=timed(()->build_window_model(data.units,data.order,data.hf,selected.result.positions))
            check();model=built.result
            flattened=timed(()->prepare_window(model));check()
            record["construction_cpu_seconds"]=built.cpu+flattened.cpu
            record["table_stats"]=model.setup_stats
            record["baseline_full_order"]=data.order
            record["baseline_window"]=model.baseline_window
            record["steps"]=affordable_steps(data.units,data.order,30data.words)
            stage="exact frozen-table enumeration"
            optimum=timed(()->enumerate_window(flattened.result;guard=check))
            stage="restricted insertion comparison"
            insertion=timed(()->insertion_window(flattened.result;guard=check))
            stage="direct full Taylor validation"
            winner=copy(data.order);winner[model.window_positions]=optimum.result.order
            direct=timed(()->score_order(data.units,winner,data.hf,30data.words))
            parent=score_order(data.units,data.order,data.hf,30data.words)
            predicted=optimum.result.score/Float64(record["steps"])^4
            abs(direct.result.value-predicted)<=max(1e-20,1e-7*abs(predicted)) || error("full Taylor differs from frozen window model")
            merge!(record,Dict("optimum"=>optimum.result,"optimum_cpu_seconds"=>optimum.cpu,
              "insertion"=>insertion.result,"insertion_cpu_seconds"=>insertion.cpu,
              "direct_validation_cpu_seconds"=>direct.cpu,"new_full_order"=>winner,
              "parent_J"=>parent.value,"new_J"=>direct.result.value,
              "frozen_table_predicted_J"=>predicted,"restricted_ratio_parent_to_optimum"=>parent.value/direct.result.value))
        end
        record["status"]="complete"
    catch exception
        record["status"]="incomplete"
        record["stage"]=stage
        record["reason"]=sprint(showerror,exception)
    end
    record["process_cpu_seconds"]=cpu_seconds()-began
    record["process_wall_seconds"]=time()-wall
    write_result(output,record)
    println(record["status"]," ",get(record,"case","unknown")," ",output)
    record["status"]=="complete" || exit(2)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
