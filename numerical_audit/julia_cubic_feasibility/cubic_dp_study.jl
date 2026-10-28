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

function one_case(input,priorpath,output;preflight=false,repetitions=3)
    prior=TOML.parsefile(priorpath)
    prior["status"]=="complete" || error("incomplete prior window certificate")
    parsed=measured(()->load_model(input;scale=preflight ? 7//5 : 1//1));data=parsed.result
    data.input_sha256==prior["input_sha256"] || error("prior model input changed")
    data.order==Int.(prior["baseline_full_order"]) || error("prior outside/parent order changed")
    positions=prior["window_start"]:prior["window_start"]+prior["window_width"]-1
    first(positions)>1 && last(positions)<length(data.order) || error("require fixed exterior endpoints")
    common=measured(()->begin
        model=build_window_model(data.units,data.order,data.hf,positions)
        (model=model,flat=prepare_window(model))
    end)
    model=common.result.model;flat=common.result.flat;digest=flat_digest(flat)
    steps=affordable_steps(data.units,data.order,30data.words)
    steps==prior["steps"] || error("affordable step count changed")
    expected=preflight ? nothing : rational_string(prior["optimum"]["numerator"])
    expected_order=preflight ? nothing : Int.(prior["optimum"]["order"])
    records=Dict{String,Any}[]
    function run(method,stage)
        GC.gc();deadline=cpu_seconds()+60.0
        guard=()->(cpu_seconds()<deadline || error("method CPU deadline"))
        sample=measured(()->begin
            result=method=="dp" ? solve_cubic_dp(flat;guard=guard) : enumerate_window(flat;guard=guard)
            completed=method=="dp" ? result.complete : result.status=="complete"
            if completed
                value=method=="dp" ? result.numerator//result.denominator : result.numerator
                # DP returns global labels; direct vector evaluation uses1:k.
                local_order=Int[findfirst(==(label),flat.labels) for label in result.order]
                exact_numerator(flat,local_order)==value || error("returned cost differs from direct dyadic vectors")
                expected===nothing || value==expected || error("exact prior optimum differs")
                expected_order===nothing || result.order==expected_order || error("prior unique optimum witness differs")
            end
            evidence_solver(result)
        end)
        result=sample.result
        complete=method=="dp" ? result.complete : result.status=="complete"
        item=Dict("method"=>method,"stage"=>stage,"cpu_seconds"=>sample.cpu_seconds,
            "wall_seconds"=>sample.wall_seconds,"allocated_bytes"=>sample.allocated_bytes,
            "gc_wall_seconds"=>sample.gc_wall_seconds,"complete"=>complete,
            "flat_sha256"=>digest,"solver"=>result)
        push!(records,item)
        write_result(output*"."*stage*"_"*method*".toml",item)
        println(data.case," ",stage," ",method," ",complete," CPU=",sample.cpu_seconds);flush(stdout)
        item
    end
    # Complete enumeration supplies a new exact oracle for scaled preflight;
    # no homogeneity of the diagonal transport is assumed under H scaling.
    enumeration_first=run("enumeration","first")
    enum_result=enumeration_first["solver"]
    if expected===nothing
        expected=enum_result.numerator
        expected_order=enum_result.minimizers==1 ? enum_result.order : nothing
    end
    enum_result.permutations==factorial(length(flat.labels)) || error("incomplete exhaustive domain")
    dp_first=run("dp","first")
    if !preflight
        for repetition in 1:repetitions
            for method in (isodd(repetition) ? ("dp","enumeration") : ("enumeration","dp"))
                method=="dp" && !dp_first["complete"] && continue
                run(method,"warm"*string(repetition))
            end
        end
    end
    # One full Taylor check is a separate validation cost rather than a hidden
    # contribution to only one method's benchmark.
    complete_dp=[s for s in records if s["method"]=="dp" && s["complete"]]
    parity=nothing;reuse=nothing
    if !isempty(complete_dp)
        result=first(complete_dp)["solver"]
        full=copy(data.order);full[positions]=result.order
        parity=measured(()->score_order(data.units,full,data.hf,30data.words))
        predicted=Float64(expected)/Float64(steps)^4
        abs(parity.result.value-predicted)<=max(1e-20,1e-7*abs(predicted)) || error("direct Taylor parity differs")
        preflight || (reuse=physical_binding(data,full,prior))
    end
    record=Dict("case"=>data.case,"input_sha256"=>data.input_sha256,"prior_sha256"=>bytes2hex(sha256(read(priorpath))),
       "flat_sha256"=>digest,"steps"=>steps,"width"=>length(flat.labels),"window_positions"=>collect(positions),
       "baseline_window"=>flat.labels,"common_response_setup_cpu_seconds"=>common.cpu_seconds,
       "parse_cpu_seconds"=>parsed.cpu_seconds,"expected_Q"=>expected,
       "expected_J"=>expected/BigInt(steps)^4,"expected_order"=>expected_order,"records"=>records,
       "full_Taylor_validation"=>parity,
       "status"=>isempty(complete_dp) ? "completed_with_DP_refusal" : "complete",
       "source_hashes"=>source_hashes(),"julia_version"=>string(VERSION),
       "physical_reuse"=>reuse!==nothing,"physical_binding"=>reuse,
       "scope"=>"Exact supplied cubic response DELTA tables; fixed outside and steps. No exact-transcendental/global/physical-optimality certificate.")
    write_result(output,record)
    record
end

function main(args)
    length(args)>=4 && iszero((length(args)-1)%3) || error("usage: cubic_dp_study.jl preflight|main INPUT PRIOR OUTPUT [INPUT PRIOR OUTPUT...]")
    preflight=args[1]=="preflight";args[1] in ("preflight","main") || error("unknown mode")
    for i in 2:3:length(args)
        one_case(args[i],args[i+1],args[i+2];preflight=preflight)
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
