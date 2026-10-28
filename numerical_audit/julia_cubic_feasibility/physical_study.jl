using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
include("src/ResponseCore.jl")
include("src/StudyIO.jl")
include("src/PhysicalCheck.jl")
using .ResponseCore, .StudyIO, .PhysicalCheck

function archive_orders(case,data)
    orders=Dict{String,Vector{Int}}("parent"=>data.order)
    skips=String[]
    for sample in case["first_calls"]
        solver=sample["solver"]
        if solver["complete"]
            solver["validation"]["valid"] || error("uncertified archived order")
            # Archived global labels are 0-based; Unit arrays are 1-based.
            orders[sample["method"]]=Int.(solver["validation"]["full_order"]).+1
        else
            push!(skips,sample["method"]*": "*solver["reason"])
        end
    end
    orders,skips
end

function frozen_rational(value)
    value isa Integer && return BigInt(value)//big(1)
    value isa Rational && return value
    pieces=split(value,"//")
    length(pieces)==2 || error("invalid frozen window numerator")
    parse(BigInt,pieces[1])//parse(BigInt,pieces[2])
end

function cubic_orders(case,data)
    case["status"]=="complete" && case["mode"]=="window" || error("incomplete or scaled cubic window result")
    Int.(case["baseline_full_order"])==data.order || error("physical cubic parent order differs")
    at=case["window_start"];width=case["window_width"];positions=at:at+width-1
    first(positions)>1 && last(positions)<length(data.order) || error("invalid interior physical window")
    baseline=data.order[positions]
    Int.(case["baseline_window"])==baseline || error("physical window labels differ")
    optimum=case["optimum"];insertion=case["insertion"]
    optimum["status"]=="complete" && optimum["permutations"]==factorial(big(width)) &&
        optimum["minimizers"]>=1 || error("incomplete cubic optimum")
    for result in (optimum,insertion)
        order=Int.(result["order"])
        length(order)==width && Set(order)==Set(baseline) || error("physical window solver changed labels")
        score=result["score"];value=frozen_rational(result["numerator"])
        isfinite(score) && score>=0 && value>=0 && score==Float64(value) || error("physical window score mismatch")
    end
    frozen_rational(insertion["numerator"])>=frozen_rational(optimum["numerator"]) || error("insertion beat certified optimum")
    winner=copy(data.order);winner[positions]=Int.(optimum["order"])
    winner==Int.(case["new_full_order"]) || error("physical optimum/full-order mismatch")
    predicted=optimum["score"]/Float64(case["steps"])^4
    case["steps"]>0 && isfinite(case["new_J"]) && case["new_J"]>=0 &&
        abs(case["new_J"]-predicted)<=max(1e-20,1e-7*abs(predicted)) || error("physical window Taylor normalization differs")
    inserted=copy(data.order);inserted[positions]=Int.(insertion["order"])
    Dict("parent"=>copy(data.order),"window_exact"=>winner,"window_insertion"=>inserted)
end

function continuation_control!(orders,skips,resultpath,data)
    controlpath=joinpath(dirname(resultpath),data.case*".CONTINUATION.a1.toml")
    if isfile(controlpath)
        control=TOML.parsefile(controlpath)
        control["model_sha256"]==data.source["source_model_sha256"] || error("legacy control model differs")
        control["parent_sha256"]==data.source["cubic_source_sha256"] || error("legacy control parent differs")
        orders["legacy_continuation"]=Int.(control["order"])
    else
        push!(skips,"legacy_continuation: optional continuation file absent: "*controlpath)
    end
    nothing
end

function main(args)
    length(args)>=4 || error("usage: physical_study.jl preflight|archive|cubic MODEL_PATH RESULT... OUTPUT")
    mode,modelpath=args[1:2];output=args[end];results=args[3:end-1]
    began=cpu_seconds();wall=time();record=Dict{String,Any}(
      "mode"=>mode,"julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),
      "blas_threads"=>BLAS.get_num_threads(),"norm_drift_definition"=>"abs(squared norm - 1)",
      "circuit_space"=>"full Hilbert space; no sector projection",
      "reference_space"=>"fixed N only after checking discarded Hamiltonian couplings",
      "reference_check_scope"=>"Maximum-entry numerical coupling/Hermiticity checks; not operator-norm or interval certificates",
      "source_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in
       [@__FILE__;[joinpath(@__DIR__,"src",n) for n in ("PhysicalCheck.jl","ResponseCore.jl","StudyIO.jl")]]))
    cases=Dict{String,Any}[]
    try
        if mode=="preflight"
            data=load_model(modelpath);sample=TOML.parsefile(only(results))
            sample["solver"]["validation"]["valid"] || error("uncertified preflight DP order")
            orders=Dict("parent"=>data.order,"dp"=>Int.(sample["solver"]["validation"]["full_order"]).+1)
            rows=physical_compare(data.units,data.n,data.hf,orders;times=[1.0],budget_multipliers=[10],formulas=[1,2])
            push!(cases,Dict("case"=>data.case,"input_sha256"=>data.input_sha256,"rows"=>rows))
        elseif mode=="archive"
            length(results)==1 || error("one archived benchmark result required")
            benchmark=TOML.parsefile(only(results));benchmark["status"]=="complete" || error("incomplete archive benchmark")
            manifest=TOML.parsefile(joinpath(modelpath,"MANIFEST.toml"))
            for case in benchmark["cases"]
                data=load_model(joinpath(modelpath,case["case"]*".toml"))
                binding=manifest[data.case]
                data.input_sha256==binding["model_toml_sha256"] || error("physical archive model bytes changed")
                data.source["source_model_sha256"]==binding["source_model_sha256"] || error("physical archive source model differs")
                case["input_sha256"]==binding["domain_toml_sha256"] || error("physical archive domain bytes differ")
                case["metadata"]["domain_sha256"]==binding["source_domain_sha256"] || error("physical archive domain provenance differs")
                data.order==Int.(binding["baseline_order"]) || error("physical archive parent differs")
                orders,skips=archive_orders(case,data)
                rows=physical_compare(data.units,data.n,data.hf,orders)
                entry=Dict("case"=>data.case,"input_sha256"=>data.input_sha256,"rows"=>rows,"skips"=>skips)
                push!(cases,entry);write_result(output*"."*data.case*".checkpoint.toml",entry)
            end
        elseif mode=="cubic"
            for resultpath in results
                case=TOML.parsefile(resultpath);case["status"]=="complete" || error("incomplete cubic window result")
                data=load_model(joinpath(modelpath,case["case"]*".toml"))
                case["input_sha256"]==data.input_sha256 || error("physical cubic input bytes differ")
                orders=cubic_orders(case,data)
                skips=String[]
                continuation_control!(orders,skips,resultpath,data)
                rows=physical_compare(data.units,data.n,data.hf,orders;times=[0.5,1.5],formulas=[2])
                entry=Dict("case"=>data.case,"input_sha256"=>data.input_sha256,"rows"=>rows,"skips"=>skips)
                push!(cases,entry);write_result(output*"."*data.case*".checkpoint.toml",entry)
            end
        else
            error("unknown physical mode")
        end
        record["status"]="complete"
    catch exception
        record["status"]="incomplete";record["reason"]=sprint(showerror,exception)
    end
    record["cases"]=cases;record["process_cpu_seconds"]=cpu_seconds()-began
    record["process_wall_seconds"]=time()-wall
    write_result(output,record);println(record["status"]," ",output)
    record["status"]=="complete" || exit(2)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
