using LinearAlgebra, TOML, SHA, Random
BLAS.set_num_threads(1)
include("src/ResponseCore.jl")
include("src/StructureDP.jl")
include("src/StudyIO.jl")
include("src/WindowSearch.jl")
include("src/CubicDP.jl")
using .ResponseCore, .StructureDP, .StudyIO, .WindowSearch
using .CubicDP: cubic_factors

const CHECK_SEED=20261028
const CHECK_SOURCE_PATHS=("cubic_dp_study.jl","src/ResponseCore.jl","src/StructureDP.jl",
    "src/StudyIO.jl","src/WindowSearch.jl","src/CubicDP.jl")

check(condition,message)=condition || error(message)
filehash(path)=bytes2hex(sha256(read(path)))
asbig(value)=value isa Integer ? BigInt(value) : parse(BigInt,value)

function exact_rational(value)
    value isa Integer && return BigInt(value)//big(1)
    value isa Rational && return BigInt(numerator(value))//BigInt(denominator(value))
    pieces=split(value,"//")
    length(pieces)==2 || error("expected exact integer/rational serialization")
    parse(BigInt,pieces[1])//parse(BigInt,pieces[2])
end

"""Reproduce the documented bit-stream format without invoking the study driver."""
function checked_flat_digest(flat)
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

function verify_sources(main,prior)
    recorded=main["source_hashes"]
    check(Set(keys(recorded))==Set(CHECK_SOURCE_PATHS),"main source inventory differs")
    verified=Dict{String,String}()
    for path in CHECK_SOURCE_PATHS
        digest=filehash(joinpath(@__DIR__,path))
        check(digest==recorded[path],"measured source bytes changed: "*path)
        verified[path]=digest
    end
    # The old response kernel is retained separately; it need not equal the
    # current file after the documented invalid-input admission-guard changes.
    snapshot=joinpath(@__DIR__,"evidence/source_at_measurements/ResponseCore.jl")
    oldhash=filehash(snapshot)
    check(oldhash==prior["source_sha256"]["ResponseCore.jl"],"prior response snapshot binding differs")
    (current=verified,prior_response_snapshot=oldhash,checker=filehash(@__FILE__))
end

function deterministic_orders(k)
    check(k==8,"this independent decision study requires the frozen eight-label windows")
    identity=collect(1:k)
    orders=[copy(identity),reverse(identity)]
    seen=Set(Tuple(order) for order in orders)
    rng=MersenneTwister(CHECK_SEED)
    draws=0
    while length(orders)<32
        draws+=1;draws<=10000 || error("seeded order sampler failed its finite draw cap")
        candidate=shuffle(rng,identity)
        key=Tuple(candidate)
        key in seen && continue
        push!(seen,key);push!(orders,candidate)
    end
    orders,draws
end

function local_order(flat,global_order)
    proposed=Int.(global_order)
    check(length(proposed)==length(flat.labels) && Set(proposed)==Set(flat.labels),"solver changed the frozen window label set")
    mapping=Dict(label=>i for (i,label) in enumerate(flat.labels))
    result=Int[mapping[label] for label in proposed]
    check(sort(result)==collect(1:length(flat.labels)),"invalid local permutation")
    result
end

function verify_physical_reuse(data,main,full)
    check(main["physical_reuse"]===true,"main did not bind retained physical evidence")
    binding=main["physical_binding"]
    source="evidence/results/PHYSICAL_CUBIC.a1.toml"
    check(binding["source"]==source,"unexpected retained physical source")
    path=joinpath(@__DIR__,source)
    check(filehash(path)==binding["sha256"],"retained physical source hash differs")
    check(Int.(binding["full_order"])==full && binding["rows_bound"]==4,"retained physical order binding differs")
    physical=TOML.parsefile(path)
    check(physical["status"]=="complete" && physical["mode"]=="cubic","incomplete retained physical result")
    cases=filter(c->c["case"]==data.case,physical["cases"])
    check(length(cases)==1,"ambiguous retained physical case")
    case=only(cases)
    check(case["input_sha256"]==data.input_sha256,"retained physical input differs")
    rows=filter(row->row["name"]=="window_exact",case["rows"])
    expected=Set([(0.5,10),(0.5,30),(1.5,10),(1.5,30)])
    check(length(rows)==4 && Set((row["time"],row["budget_multiplier"]) for row in rows)==expected,
          "retained physical condition set differs")
    for row in rows
        check(Int.(row["order"])==full && row["formula"]==2,"retained circuit order/formula differs")
        budget=row["budget_multiplier"]*data.words
        check(row["rotation_budget"]==budget,"retained circuit budget differs")
        check(affordable_steps(data.units,full,budget)==row["steps"],"retained circuit affordable steps differ")
    end
    (rows_bound=length(rows),source=source,sha256=binding["sha256"],new_propagation_performed=false)
end

"""Derive the inventory independently; unlabelled historical main runs had three warms."""
function checked_run_inventory(main)
    labelled=haskey(main,"mode") || haskey(main,"requested_warm_repetitions")
    if labelled
        check(haskey(main,"mode") && haskey(main,"requested_warm_repetitions"),"incomplete run metadata")
        mode=main["mode"];warms=main["requested_warm_repetitions"]
        check(warms isa Integer && warms>=0,"invalid warm repetition count")
        check((mode=="main" && warms==3) || (mode=="validation" && warms!=3) ||
              (mode=="preflight" && warms==0),"mode/repetition contract differs")
    else
        mode="main";warms=3
    end
    expected=[(method="enumeration",stage="first"),(method="dp",stage="first")]
    for repetition in 1:warms
        for method in (isodd(repetition) ? ("dp","enumeration") : ("enumeration","dp"))
            push!(expected,(method=method,stage="warm"*string(repetition)))
        end
    end
    if haskey(main,"scheduled_runs")
        declared=[(method=item["method"],stage=item["stage"]) for item in main["scheduled_runs"]]
        check(declared==expected,"declared run inventory differs from mode/repetition contract")
    end
    items=main["records"]
    firstdp=filter(item->item["method"]=="dp" && item["stage"]=="first",items)
    check(length(firstdp)==1,"missing or duplicate first DP record")
    refused=!only(firstdp)["complete"]
    required=refused ? filter(run->run.method!="dp" || run.stage=="first",expected) : expected
    observed=[(method=item["method"],stage=item["stage"]) for item in items]
    check(observed==required,"missing, duplicate, reordered, or unexpected solver records")
    for item in items
        solver=item["solver"]
        if !item["complete"]
            check(refused && item["method"]=="dp" && item["stage"]=="first","incomplete scheduled solver record")
            allowed=("work_limit","memory_limit","bit_limit","vector_entry_limit","factor_entry_limit",
                "cost_bit_limit","bag_state_limit","state_limit","graph_incomplete","decomposition_incomplete")
            check(!solver["complete"] && solver["status"] in allowed && !isempty(solver["reason"]),"invalid first DP refusal")
            check(isempty(get(solver,"order",Int[])) && get(solver,"numerator","unavailable")=="unavailable" &&
                  get(solver,"score","unavailable")=="unavailable","refused DP supplied an optimum witness")
        elseif item["method"]=="dp"
            check(solver["complete"] && solver["status"]=="complete","DP completion status differs")
        else
            check(solver["status"]=="complete" && solver["permutations"]==factorial(big(main["width"])),
                  "recorded exhaustive domain differs")
        end
    end
    outcome=refused ? "completed_with_DP_refusal" : "complete"
    check(main["status"]==outcome,"producer completion flag disagrees with independent inventory")
    check(get(main,"DP_optimum_certified",!refused)==!refused,"DP certificate flag differs")
    if haskey(main,"skipped_runs")
        skipped=[(method=item["method"],stage=item["stage"]) for item in main["skipped_runs"]]
        wanted=refused ? filter(run->run.method=="dp" && run.stage!="first",expected) : NamedTuple[]
        check(skipped==wanted,"skipped run inventory differs")
    end
    (mode=mode,warms=warms,refused=refused,expected_records=length(required))
end

middle3(values)=begin
    check(length(values)==3,"warm summary requires exactly three samples")
    sort!(Float64[values...])[2]
end

function warm_summary(records,repetitions=3)
    summary=Dict{String,Any}("aggregation"=>
        "Every available field is its own sorted-middle-of-three median; ratios divide the two CPU medians, not paired ratios. No new timing is measured by this postprocessing.",
        "requested_warm_repetitions"=>repetitions)
    fields=("integerization_CPU_seconds","factor_CPU_seconds","graph_decomposition_CPU_seconds",
            "solve_CPU_seconds","witness_CPU_seconds")
    for method in ("dp","enumeration")
        warm=filter(item->item["method"]==method && startswith(item["stage"],"warm") && item["complete"],records)
        if repetitions!=3 || length(warm)!=3
            summary[method]=Dict("available"=>false,"completed_warm_samples"=>length(warm))
            continue
        end
        result=Dict{String,Any}("available"=>true,"samples"=>3,
            "CPU_median_seconds"=>middle3(item["cpu_seconds"] for item in warm),
            "wall_median_seconds"=>middle3(item["wall_seconds"] for item in warm))
        if method=="dp"
            result["phase_CPU_median_seconds"]=Dict(field=>middle3(item["solver"]["stats"][field] for item in warm) for field in fields)
        end
        summary[method]=result
    end
    if summary["dp"]["available"] && summary["enumeration"]["available"]
        summary["enumeration_over_DP_CPU_median_ratio"]=summary["enumeration"]["CPU_median_seconds"]/summary["dp"]["CPU_median_seconds"]
    end
    summary
end

function check_case(input,priorpath,resultpath,output,guard)
    began=cpu_seconds();wall=time_ns()
    record=Dict{String,Any}("status"=>"incomplete","input_sha256"=>filehash(input),
        "prior_sha256"=>filehash(priorpath),"main_result_sha256"=>filehash(resultpath),
        "julia_version"=>string(VERSION),"seed"=>CHECK_SEED,
        "factor_identities_checked"=>0,"recorded_solver_outputs_checked"=>0,
        "scope"=>"Independent nonoptimal-order checks of exact frozen-vector/scalar-factor identity and recorded witness bindings. No DP, exhaustive enumeration, or new physical propagation is performed.")
    try
        guard()
        prior=TOML.parsefile(priorpath);main=TOML.parsefile(resultpath)
        check(prior["status"]=="complete","prior window certificate incomplete")
        inventory=checked_run_inventory(main)
        check(inventory.mode!="preflight","scaled preflight is outside this unscaled independent checker")
        record["mode"]=inventory.mode;record["requested_warm_repetitions"]=inventory.warms
        record["DP_optimum_certified"]=!inventory.refused
        check(record["prior_sha256"]==main["prior_sha256"],"main bound a different prior certificate")
        record["source_checks"]=verify_sources(main,prior)
        data=load_model(input)
        check(data.case==main["case"]==prior["case"],"case identity differs")
        check(data.input_sha256==main["input_sha256"]==prior["input_sha256"],"model input hash differs")
        check(data.order==Int.(prior["baseline_full_order"]),"parent/outside order differs")
        positions=prior["window_start"]:prior["window_start"]+prior["window_width"]-1
        check(first(positions)>1 && last(positions)<length(data.order),"window changed global endpoints")
        check(collect(positions)==Int.(main["window_positions"]),"main window positions differ")
        model=build_window_model(data.units,data.order,data.hf,positions)
        flat=prepare_window(model)
        check(main["width"]==length(flat.labels),"recorded window width differs")
        check(flat.labels==Int.(main["baseline_window"])==Int.(prior["baseline_window"]),"window baseline labels differ")
        digest=checked_flat_digest(flat)
        check(digest==main["flat_sha256"],"reconstructed frozen Float64 table bits differ")
        record["case"]=data.case;record["flat_sha256"]=digest
        steps=affordable_steps(data.units,data.order,30data.words)
        check(steps==15 && steps==main["steps"]==prior["steps"],"fixed affordable step count differs")
        expected=exact_rational(prior["optimum"]["numerator"])
        check(prior["optimum"]["minimizers"]==1,"prior certificate no longer asserts a unique optimum")
        winner=Int.(prior["optimum"]["order"])
        check(expected==exact_rational(main["expected_Q"]),"main exact Q differs from prior")
        check(expected/BigInt(steps)^4==exact_rational(main["expected_J"]),"main exact J normalization differs")
        check(Int.(main["expected_order"])==winner,"main unique optimum differs from prior")
        exact_data=cubic_factors(flat;guard=guard)
        check(exact_data.complete,"independent exact factor reconstruction refused: "*exact_data.reason)
        check(exact_data.vector_denominator^2==exact_data.denominator,"D squared normalization differs")
        record["factor_stats"]=exact_data.stats
        orders,draws=deterministic_orders(length(flat.labels))
        record["sampled_local_orders"]=orders;record["random_draws"]=draws
        for order in orders
            guard()
            integer=evaluate_factors(exact_data.factors,order)
            check(denominator(integer)==1,"scalar-factor total is not integer N")
            direct=exact_numerator(flat,order)
            check(integer/exact_data.denominator==direct,"nonoptimal-order factor/vector identity differs")
            record["factor_identities_checked"]+=1
        end
        complete_dp=0;complete_enumeration=0
        for item in main["records"]
            guard();check(item["flat_sha256"]==digest,"record used different frozen table bits")
            if !item["complete"]
                # checked_run_inventory already admitted only the explicit first-DP refusal.
                check(inventory.refused && item["method"]=="dp" && item["stage"]=="first","unadmitted incomplete record")
                record["DP_refusal"]=item["solver"]
                continue
            end
            solver=item["solver"];method=item["method"]
            local_perm=local_order(flat,solver["order"])
            direct=exact_numerator(flat,local_perm)
            check(Int.(solver["order"])==winner,"record differs from prior unique order")
            if method=="dp"
                complete_dp+=1
                check(solver["complete"] && solver["status"]=="complete","DP completion status differs")
                N=asbig(solver["numerator"]);D=asbig(solver["vector_denominator"]);H=asbig(solver["denominator"])
                check(D^2==H && H==exact_data.denominator && D==exact_data.vector_denominator,"recorded integerization scale differs")
                check(Int.(solver["local_order"])==local_perm,"DP local/global mapping differs")
                check(flat.labels[local_perm]==Int.(solver["order"]),"DP global relabeling differs")
                check(exact_rational(solver["dp"]["cost"])==N,"generic DP cost differs from N")
                check(exact_rational(solver["score"])==N//H==direct==expected,"recorded DP exact objective differs")
            elseif method=="enumeration"
                complete_enumeration+=1
                check(solver["status"]=="complete" && solver["permutations"]==factorial(length(flat.labels)),"recorded exhaustive domain differs")
                check(solver["minimizers"]==1,"new exhaustive unique optimum differs")
                check(exact_rational(solver["numerator"])==direct==expected,"recorded exhaustive exact objective differs")
            else
                error("unknown recorded solver method")
            end
            record["recorded_solver_outputs_checked"]+=1
        end
        check(complete_enumeration==1+inventory.warms,"exhaustive record count differs")
        check(complete_dp==(inventory.refused ? 0 : 1+inventory.warms),"DP record count differs")
        record["complete_DP_records_checked"]=complete_dp
        record["complete_enumeration_records_checked"]=complete_enumeration
        record["warm_summary"]=warm_summary(main["records"],inventory.warms)
        if !inventory.refused
            full=copy(data.order);full[positions]=winner
            check(full==Int.(prior["new_full_order"]),"reconstructed full winner differs")
            check(first(full)==first(data.order) && last(full)==last(data.order),"winner changed fixed endpoints")
            check(full[1:first(positions)-1]==data.order[1:first(positions)-1] &&
                  full[last(positions)+1:end]==data.order[last(positions)+1:end],"winner changed outside positions")
            check(affordable_steps(data.units,full,30data.words)==steps,"winner changed affordable steps")
            record["physical_evidence_binding"]=verify_physical_reuse(data,main,full)
            guard()
            fresh=score_order(data.units,full,data.hf,30data.words)
            predicted=Float64(expected)/Float64(steps)^4
            difference=abs(fresh.value-predicted)
            check(difference<=max(1e-20,1e-7*abs(predicted)),"independent full-Taylor numerical parity differs")
            old=main["full_Taylor_validation"]["result"]["value"]
            recorded_difference=abs(old-predicted)
            check(recorded_difference<=max(1e-20,1e-7*abs(predicted)),"recorded full-Taylor numerical parity differs")
            maximum_difference=max(difference,recorded_difference)
            record["full_Taylor_numerical_parity"]=(fresh_J=fresh.value,recorded_J=old,
                predicted_J=predicted,max_absolute_difference=maximum_difference,
                max_relative_difference=iszero(predicted) ? (iszero(maximum_difference) ? 0.0 : Inf) : maximum_difference/abs(predicted),
                steps=fresh.steps,exact_transcendental_or_interval_claim=false)
            record["expected_Q"]=expected;record["expected_J"]=expected/BigInt(steps)^4
            record["status"]="passed"
        else
            check(main["physical_reuse"]===false,"refused DP claimed retained physical binding")
            check(main["full_Taylor_validation"]=="unavailable","refused DP claimed full Taylor validation")
            record["status"]="passed_with_DP_refusal"
            record["interpretation"]="Explicit DP admission refusal checked; no DP optimum or physical binding certified."
        end
    catch exception
        record["status"]="failed";record["reason"]=sprint(showerror,exception)
    end
    record["CPU_seconds"]=cpu_seconds()-began
    record["wall_seconds"]=(time_ns()-wall)/1e9
    write_result(output,record)
    println(record["status"]," ",output);flush(stdout)
    record["status"] in ("passed","passed_with_DP_refusal")
end

function main(args)
    length(args)>=4 && iszero(length(args)%4) || error("usage: independent_cubic_dp_check.jl INPUT PRIOR MAIN_RESULT OUTPUT [INPUT PRIOR MAIN_RESULT OUTPUT...]")
    deadline=cpu_seconds()+90.0
    guard=()->(cpu_seconds()<deadline || error("independent verification CPU deadline"))
    success=true
    for i in 1:4:length(args)
        success=check_case(args[i],args[i+1],args[i+2],args[i+3],guard) && success
    end
    success || exit(2)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
