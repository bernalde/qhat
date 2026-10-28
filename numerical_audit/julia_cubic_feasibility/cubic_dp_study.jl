using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
include("src/ResponseCore.jl")
include("src/StructureDP.jl")
include("src/StudyIO.jl")
include("src/WindowSearch.jl")
include("src/CubicDP.jl")
using .ResponseCore, .StructureDP, .StudyIO, .WindowSearch, .CubicDP

function measured(f)
    start=cpu_seconds();value=@timed f()
    (result=value.value,cpu_seconds=cpu_seconds()-start,wall_seconds=value.time,
     allocated_bytes=value.bytes,gc_wall_seconds=value.gctime)
end

function evidence_solver(result::NamedTuple)
    kept=Tuple(k for k in keys(result) if !(k in (:audit,:decomposition)) &&
        !(k==:factors && !(getproperty(result,k) isa Integer)))
    length(kept)==length(result) && return result
    NamedTuple{kept}(Tuple(getproperty(result,k) for k in kept))
end

function rational_string(value)
    pieces=split(value,"//")
    length(pieces)==2 || error("invalid archived exact rational")
    parse(BigInt,pieces[1])//parse(BigInt,pieces[2])
end

"""Canonical serialized bit-pattern identity of the supplied FlatWindow."""
function flat_digest(flat)
    io=IOBuffer()
    println(io,"flat-window-v1")
    println(io,join(flat.labels,','));println(io,join(flat.destinations,','))
    for value in flat.baseline
        println(io,bitstring(real(value)),bitstring(imag(value)))
    end
    for factor in flat.factors
        println(io,join(factor.scope,','));println(io,join(factor.destinations,','))
        println(io,size(factor.values))
        for value in factor.values
            println(io,bitstring(real(value)),bitstring(imag(value)))
        end
    end
    bytes2hex(sha256(take!(io)))
end

function source_hashes()
    paths=["cubic_dp_study.jl","src/ResponseCore.jl","src/StructureDP.jl",
        "src/StudyIO.jl","src/WindowSearch.jl","src/CubicDP.jl"]
    Dict(path=>bytes2hex(sha256(read(joinpath(@__DIR__,path)))) for path in paths)
end

function physical_binding(data,full,prior)
    full==Int.(prior["new_full_order"]) || error("prior full-order witness differs")
    path=joinpath(@__DIR__,"evidence/results/PHYSICAL_CUBIC.a1.toml")
    physical=TOML.parsefile(path)
    physical["status"]=="complete" && physical["mode"]=="cubic" || error("invalid physical archive")
    matched=only(filter(c->c["case"]==data.case,physical["cases"]))
    matched["input_sha256"]==data.input_sha256 || error("physical input differs")
    rows=filter(row->row["name"]=="window_exact",matched["rows"])
    length(rows)==4 || error("missing archived physical conditions")
    Set((row["time"],row["budget_multiplier"]) for row in rows)==Set([(0.5,10),(0.5,30),(1.5,10),(1.5,30)]) || error("physical conditions differ")
    for row in rows
        Int.(row["order"])==full && row["formula"]==2 || error("physical order/formula differs")
        row["rotation_budget"]==row["budget_multiplier"]*data.words || error("physical budget differs")
        affordable_steps(data.units,full,row["rotation_budget"])==row["steps"] || error("physical step count differs")
    end
    (source="evidence/results/PHYSICAL_CUBIC.a1.toml",sha256=bytes2hex(sha256(read(path))),
     full_order=full,rows_bound=length(rows),conditions=[(time=row["time"],
       budget_multiplier=row["budget_multiplier"],steps=row["steps"]) for row in rows],
     interpretation="Reused archived full-circuit evidence for the identical input and full order; no new physical simulation.")
end

function scheduled_runs(repetitions)
    runs=[(method="enumeration",stage="first"),(method="dp",stage="first")]
    for repetition in 1:repetitions
        for method in (isodd(repetition) ? ("dp","enumeration") : ("enumeration","dp"))
            push!(runs,(method=method,stage="warm"*string(repetition)))
        end
    end
    runs
end

"""Completion requires the whole declared schedule; a first DP refusal is explicit."""
function benchmark_status(records,width,repetitions)
    expected=scheduled_runs(repetitions)
    firstdp=filter(item->item["method"]=="dp" && item["stage"]=="first",records)
    refused=length(firstdp)==1 && !only(firstdp)["complete"] &&
        !only(firstdp)["solver"].complete &&
        string(only(firstdp)["solver"].status) in ("work_limit","memory_limit","bit_limit",
          "vector_entry_limit","factor_entry_limit","cost_bit_limit","bag_state_limit",
          "state_limit","graph_incomplete","decomposition_incomplete") &&
        !isempty(only(firstdp)["solver"].reason)
    required=refused ? filter(run->run.method!="dp" || run.stage=="first",expected) : expected
    observed=[(method=item["method"],stage=item["stage"]) for item in records]
    observed==required || return "incomplete"
    for item in records
        solver=item["solver"]
        if item["method"]=="enumeration"
            item["complete"] && string(solver.status)=="complete" &&
                solver.permutations==factorial(big(width)) || return "incomplete"
        elseif item["complete"]
            solver.complete && string(solver.status)=="complete" || return "incomplete"
        elseif !(refused && item["stage"]=="first")
            return "incomplete"
        end
    end
    refused ? "completed_with_DP_refusal" : "complete"
end

function checkpoint_sample(evaluate,method,runstage,digest,output,records)
    GC.gc();started=cpu_seconds();started_wall=time();deadline=started+60.0
    guard=()->(cpu_seconds()<deadline || error("method CPU deadline"))
    item=Dict{String,Any}("method"=>method,"stage"=>runstage,
        "complete"=>false,"flat_sha256"=>digest)
    try
        sample=measured(()->evaluate(guard))
        result=sample.result
        merge!(item,Dict("cpu_seconds"=>sample.cpu_seconds,"wall_seconds"=>sample.wall_seconds,
            "allocated_bytes"=>sample.allocated_bytes,"gc_wall_seconds"=>sample.gc_wall_seconds,
            "complete"=>(method=="dp" ? result.complete : result.status=="complete"),"solver"=>result))
    catch exception
        merge!(item,Dict("cpu_seconds"=>cpu_seconds()-started,"wall_seconds"=>time()-started_wall,
            "solver"=>(status="error",complete=false,reason=sprint(showerror,exception))))
        push!(records,item)
        write_result(output*"."*runstage*"_"*method*".toml",item)
        rethrow()
    end
    push!(records,item)
    write_result(output*"."*runstage*"_"*method*".toml",item)
    item
end

function one_case(input,priorpath,output;preflight=false,repetitions=3)
    repetitions isa Integer && repetitions>=0 || error("nonnegative warm repetition count required")
    ispath(output) && error("refuse to overwrite output $output")
    warms=preflight ? 0 : repetitions
    mode=preflight ? "preflight" : (warms==3 ? "main" : "validation")
    records=Dict{String,Any}[]
    record=Dict{String,Any}("status"=>"incomplete","mode"=>mode,
        "requested_warm_repetitions"=>warms,"scheduled_runs"=>scheduled_runs(warms),
        "records"=>records,"DP_optimum_certified"=>false,"physical_reuse"=>false,
        "julia_version"=>string(VERSION),
        "scope"=>"Exact supplied cubic response DELTA tables; fixed outside and steps. No exact-transcendental/global/physical-optimality certificate.")
    began=cpu_seconds();wall=time();stage="parse prior/input"
    try
        record["source_hashes"]=source_hashes()
        prior=TOML.parsefile(priorpath)
        prior["status"]=="complete" || error("incomplete prior window certificate")
        parsed=measured(()->load_model(input;scale=preflight ? 7//5 : 1//1));data=parsed.result
        data.input_sha256==prior["input_sha256"] || error("prior model input changed")
        data.order==Int.(prior["baseline_full_order"]) || error("prior outside/parent order changed")
        positions=prior["window_start"]:prior["window_start"]+prior["window_width"]-1
        first(positions)>1 && last(positions)<length(data.order) || error("require fixed exterior endpoints")
        stage="conditioned table construction"
        common=measured(()->begin
            model=build_window_model(data.units,data.order,data.hf,positions)
            (model=model,flat=prepare_window(model))
        end)
        model=common.result.model;flat=common.result.flat;digest=flat_digest(flat)
        steps=affordable_steps(data.units,data.order,30data.words)
        steps==prior["steps"] || error("affordable step count changed")
        expected=preflight ? nothing : rational_string(prior["optimum"]["numerator"])
        expected_order=preflight ? nothing : Int.(prior["optimum"]["order"])
        merge!(record,Dict("case"=>data.case,"input_sha256"=>data.input_sha256,
            "prior_sha256"=>bytes2hex(sha256(read(priorpath))),"flat_sha256"=>digest,
            "steps"=>steps,"width"=>length(flat.labels),"window_positions"=>collect(positions),
            "baseline_window"=>flat.labels,"common_response_setup_cpu_seconds"=>common.cpu_seconds,
            "parse_cpu_seconds"=>parsed.cpu_seconds))
        function run(method,runstage)
            stage=runstage*"_"*method
            item=checkpoint_sample(method,runstage,digest,output,records) do guard
                result=method=="dp" ? solve_cubic_dp(flat;guard=guard) : enumerate_window(flat;guard=guard)
                completed=method=="dp" ? result.complete : result.status=="complete"
                if completed
                    method=="enumeration" && result.permutations!=factorial(big(length(flat.labels))) &&
                        error("incomplete exhaustive domain")
                    value=method=="dp" ? result.numerator//result.denominator : result.numerator
                    local_order=Int[findfirst(==(label),flat.labels) for label in result.order]
                    exact_numerator(flat,local_order)==value || error("returned cost differs from direct dyadic vectors")
                    expected===nothing || value==expected || error("exact prior optimum differs")
                    expected_order===nothing || result.order==expected_order || error("prior unique optimum witness differs")
                end
                evidence_solver(result)
            end
            println(data.case," ",runstage," ",method," ",item["complete"]," CPU=",item["cpu_seconds"]);flush(stdout)
            item
        end
        # A scaled preflight needs its own exhaustive oracle.
        enumeration_first=run("enumeration","first")
        enumeration_first["complete"] || error("incomplete first enumeration")
        enum_result=enumeration_first["solver"]
        if expected===nothing
            expected=enum_result.numerator
            expected_order=enum_result.minimizers==1 ? enum_result.order : nothing
        end
        dp_first=run("dp","first")
        skipped=NamedTuple[]
        for run_spec in scheduled_runs(warms)[3:end]
            if run_spec.method=="dp" && !dp_first["complete"]
                push!(skipped,(method=run_spec.method,stage=run_spec.stage,reason="first DP call refused"))
                continue
            end
            run(run_spec.method,run_spec.stage)
        end
        record["skipped_runs"]=skipped
        status=benchmark_status(records,length(flat.labels),warms)
        status!="incomplete" || error("scheduled solver records did not all complete")
        parity=nothing;reuse=nothing
        if status=="complete"
            stage="full Taylor validation and physical binding"
            result=dp_first["solver"]
            full=copy(data.order);full[positions]=result.order
            parity=measured(()->score_order(data.units,full,data.hf,30data.words))
            predicted=Float64(expected)/Float64(steps)^4
            abs(parity.result.value-predicted)<=max(1e-20,1e-7*abs(predicted)) || error("direct Taylor parity differs")
            preflight || (reuse=physical_binding(data,full,prior))
        end
        merge!(record,Dict("expected_Q"=>expected,"expected_J"=>expected/BigInt(steps)^4,
            "expected_order"=>expected_order,"full_Taylor_validation"=>parity,"status"=>status,
            "DP_optimum_certified"=>status=="complete","physical_reuse"=>reuse!==nothing,"physical_binding"=>reuse))
    catch exception
        record["status"]="incomplete";record["stage"]=stage;record["reason"]=sprint(showerror,exception)
    end
    record["process_cpu_seconds"]=cpu_seconds()-began
    record["process_wall_seconds"]=time()-wall
    write_result(output,record)
    record
end

function main(args)
    length(args)>=4 && iszero((length(args)-1)%3) || error("usage: cubic_dp_study.jl preflight|main INPUT PRIOR OUTPUT [INPUT PRIOR OUTPUT...]")
    preflight=args[1]=="preflight";args[1] in ("preflight","main") || error("unknown mode")
    success=true
    for i in 2:3:length(args)
        result=one_case(args[i],args[i+1],args[i+2];preflight=preflight)
        success=result["status"]!="incomplete" && success
    end
    success || exit(2)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
