module GeneralizationObjective

# Numeric replay is the archived sparse Taylor implementation, separate from
# the constructor's native SparseResponse scores. I/O admission is study-local.
using TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"src","GeneralizationPhysical.jl"))
include(joinpath(@__DIR__,"..","prospective_window_ordering","objective_check.jl"))
const IOStudy = GeneralizationPhysical
const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const MAX_WITNESSES = 256
const MAX_UNIQUE = 256
filehash(path)=bytes2hex(sha256(read(path)))
require(ok,message)=ok || throw(ArgumentError(message))
export run, admit_witnesses

function sources()
    names=("numerical_audit/molecular_generalization/objective.jl",
        "numerical_audit/molecular_generalization/src/GeneralizationPhysical.jl",
        "numerical_audit/prospective_window_ordering/objective_check.jl",
        "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
        "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl",
        "analysis/sparse_response_ordering.py","numerical_audit/scalable_ordering/training.py")
    Dict(name=>filehash(joinpath(ROOT,name)) for name in names)
end

function admit_witnesses(manifest, candidate, unit_count)
    require(manifest isa AbstractDict && length(manifest)<=32,"bounded witness manifest")
    for key in IOStudy.HEADER
        require(haskey(manifest,key) && manifest[key]==candidate[key],"witness/candidate header join: "*key)
    end
    for key in ("protocol_path","source_freeze_path","constructor_artifacts")
        require(get(manifest,key,nothing)==candidate[key],"witness artifact roster join: "*key)
    end
    witnesses=manifest["witnesses"]
    require(witnesses isa Vector && length(witnesses)<=MAX_WITNESSES &&
        IOStudy.integer(manifest["expected_witnesses"]) && length(witnesses)==manifest["expected_witnesses"] &&
        get(manifest,"witness_inventory_complete",false)===true,"complete persisted raw-witness inventory")
    originals=Dict(row["order_id"]=>row for row in candidate["candidates"]); names=Set{String}(); claims=Set{Any}()
    for row in witnesses
        require(row isa AbstractDict && length(row)<=28 && get(row,"role",nothing)=="raw_proposed_best" &&
            IOStudy.label(row["name"]) && !(row["name"] in names),"unique raw proposed-best witness")
        push!(names,row["name"])
        require(haskey(originals,row["order_id"]),"raw witness constructor order_id")
        original=originals[row["order_id"]]
        for key in ("arm","repetition","cpu_budget","constructor_sha256","partition","score_kind")
            require(row[key]==original[key],"raw witness frozen constructor join: "*key)
        end
        require(haskey(row,"parent_sha256")==haskey(original,"parent_sha256") &&
            (!haskey(row,"parent_sha256") || row["parent_sha256"]==original["parent_sha256"]),"raw witness parent graph")
        ObjectiveCheck.validate_order(row["order"],unit_count)
        claim=(row["order_id"],Tuple(row["order"]),row["claimed_score"])
        require(!(claim in claims),"duplicate exact constructor raw proposed-best claim");push!(claims,claim)
        require(IOStudy.number(row["claimed_score"]) && row["score_kind"]==IOStudy.SCORE_KIND &&
            row["selected"] isa Bool,"raw witness score and selected status")
        row["selected"] && require(row["order"]==original["order"],"selected raw witness order differs")
        hasnative=haskey(row,"native_score"); refused=haskey(row,"native_refusal")
        require(hasnative != refused,"raw witness requires exactly one native score or refusal")
        hasnative && require(IOStudy.number(row["native_score"]),"finite raw native SparseResponse score")
        refused && require((row["native_refusal"] isa String && !isempty(row["native_refusal"])) ||
            (row["native_refusal"] isa AbstractDict && !isempty(row["native_refusal"])),"explicit nonempty native refusal record")
        row["selected"] && hasnative && require(row["native_score"]==original["claimed_score"],"selected witness native score join")
    end
    witnesses
end

function objective_record(row,kind)
    record=Dict{String,Any}("role"=>kind,"status"=>"incomplete")
    for key in (IOStudy.FIELDS...,"parent_sha256","constructor_score","constructor_failure","selected","native_score","native_refusal")
        haskey(row,key) && (record[key]=row[key])
    end
    record
end

function run(input,candidates_path,output;witnesses_path=nothing,cpu_limit=180.0,wall_limit=900.0,check_budget=()->nothing)
    require(!ispath(output),"fresh objective output directory required"); mkdir(output)
    began=IOStudy.cpu_seconds(); wall=time_ns(); stage="admission"
    guard=IOStudy.phase_guard(cpu_limit,wall_limit,check_budget)
    report=Dict{String,Any}("schema"=>"molecular-generalization-objective-v1","status"=>"incomplete",
        "complete"=>false,"all_scored_objectives_verified"=>false,"checks"=>Any[],"logical_rows"=>Any[],
        "raw_witnesses"=>Any[],"failures"=>Any[],"consumer_source_hashes"=>sources(),"new_physical_propagation"=>false,
        "training_time"=>1.0,"training_budget_multiplier"=>30,"relative_tolerance"=>1e-8,"absolute_tolerance"=>1e-24,
        "json_semantic_and_raw_inventory_validation"=>"Required separate independent Python consumer gate; Julia actual-byte binding and witness assertions alone do not certify raw JSON semantics or an empty raw-claim inventory",
        "native_reference_scope"=>"Frozen native SparseResponse scores, independently replayed and semantically joined by the owning Python consumer; Julia does not call Python or select an order",
        "archived_reference_scope"=>"Unchanged ResponseCore sparse degree3 Taylor response with compensated binary64 diagonal transport; current product exact-rational diagonal-gap transport remains a separate scope",
        "cooperation_scope"=>"Callbacks before/after each unchanged cubic score; owning OS guard bounds a score in progress")
    try
        guard(); admitted=IOStudy.admit(input,candidates_path); manifest=admitted.manifest
        for key in IOStudy.HEADER;report[key]=manifest[key];end
        report["candidates_sha256"]=admitted.candidates_sha256
        report["artifact_binding_status"]=admitted.bindings.status
        report["artifact_binding_scope"]=admitted.bindings.scope
        report["bound_artifact_hashes"]=admitted.bindings.hashes
        supplied=manifest["candidates"]
        append!(report["logical_rows"],(objective_record(row,"selected") for row in supplied))
        IOStudy.write_new(joinpath(output,"OBJECTIVE.START.toml"),report)
        data=ObjectiveCheck.WindowStudyIO.load_model(input)
        require(data.input_sha256==admitted.input_sha256,"objective model changed after admission")
        raw=Any[]
        if witnesses_path===nothing
            report["witness_inventory_complete"]=false
            report["witness_inventory_reason"]="No scored-witness manifest supplied; absence of persisted raw proposals is not inferred"
        else
            parsed=IOStudy.bounded_toml(witnesses_path)
            raw=admit_witnesses(parsed.data,manifest,length(data.units))
            report["witnesses_sha256"]=parsed.sha256;report["witness_inventory_complete"]=true
        end
        append!(report["raw_witnesses"],(objective_record(row,"raw_proposed_best") for row in raw))
        budget=ObjectiveCheck.ReplayBudget(guard)
        model=ObjectiveCheck.prepare_model(data;budget=budget)
        cache=Dict{Any,Int}()
        for (kind,rows) in (("selected",supplied),("raw_proposed_best",raw))
            destination=kind=="selected" ? report["logical_rows"] : report["raw_witnesses"]
            for (index,row) in enumerate(rows)
                stage=kind*"_"*string(index);guard()
                info=destination[index]
                order=Int.(row["order"]); key=Tuple(order)
                if !haskey(cache,key)
                    require(length(cache)<MAX_UNIQUE,"objective unique-order cap")
                    record=Dict{String,Any}("order"=>order,"score_kind"=>IOStudy.SCORE_KIND,"status"=>"started")
                    push!(report["checks"],record); cache[key]=length(report["checks"])
                    try
                        replay=ObjectiveCheck.replay_objective(model,"conserving",order,IOStudy.SCORE_KIND;budget=budget)
                        guard();record["replayed_score"]=replay.value;record["steps"]=replay.steps;record["details"]=replay.details;record["status"]="replayed"
                    catch e
                        record["status"]=e isa ObjectiveCheck.ResponseCore.NumericalResolution ? "numerical_refusal" :
                            e isa ObjectiveCheck.ObjectiveLimit ? e.status : "replay_refusal"
                        record["reason"]=sprint(showerror,e)
                        # Budget callbacks remain fatal for this finite phase;
                        # local mathematical resolution refusals remain records.
                        e isa ObjectiveCheck.ResponseCore.NumericalResolution || e isa ObjectiveCheck.ObjectiveLimit || rethrow()
                    end
                end
                info["check_index"]=cache[key]; replayed=report["checks"][cache[key]]
                if replayed["status"] != "replayed"
                    info["status"]=replayed["status"];info["reason"]=get(replayed,"reason","unavailable")
                else
                    proposed=ObjectiveCheck.compare_scores(replayed["replayed_score"],row["claimed_score"])
                    info["proposed_comparison"]=Dict(string(k)=>v for (k,v) in pairs(proposed))
                    if kind=="selected" || haskey(row,"native_score")
                        native_score=kind=="selected" ? row["claimed_score"] : row["native_score"]
                        native=ObjectiveCheck.compare_scores(replayed["replayed_score"],native_score)
                        info["native_score"]=native_score
                        info["native_comparison"]=Dict(string(k)=>v for (k,v) in pairs(native))
                        info["status"]=proposed.verified && native.verified ? "verified" :
                            (!proposed.within_tolerance || !native.within_tolerance) ? "mismatch" : "floor_censored"
                    else
                        info["status"]="native_refusal_censored"
                    end
                end
                IOStudy.write_new(joinpath(output,"OBJECTIVE.PARTIAL."*kind*"."*lpad(string(index),3,'0')*".toml"),report)
            end
        end
        report["processed_logical_rows"]=count(row -> row["status"] != "incomplete",report["logical_rows"])
        report["expected_raw_witnesses"]=length(raw);report["processed_raw_witnesses"]=count(row -> row["status"] != "incomplete",report["raw_witnesses"])
        report["all_scored_objectives_verified"]=all(r->r["status"]=="verified",vcat(report["logical_rows"],report["raw_witnesses"]))
        report["complete"]=report["processed_logical_rows"]==manifest["expected_logical_rows"] &&
            report["witness_inventory_complete"] && report["all_scored_objectives_verified"]
        report["status"]=report["complete"] ? "complete" : "partial_or_censored"
        report["charged_work_compile"]=budget.work
    catch e
        push!(report["failures"],Dict("phase"=>stage,"type"=>string(typeof(e)),"message"=>sprint(showerror,e)))
    end
    report["processed_logical_rows"]=count(row -> row["status"] != "incomplete",report["logical_rows"])
    report["processed_raw_witnesses"]=count(row -> row["status"] != "incomplete",report["raw_witnesses"])
    allrows=vcat(report["logical_rows"],report["raw_witnesses"])
    report["mismatched_objectives"]=count(r->r["status"]=="mismatch",allrows)
    report["floor_censored_objectives"]=count(r->r["status"]=="floor_censored",allrows)
    report["native_refused_objectives"]=count(r->r["status"]=="native_refusal_censored",allrows)
    report["logical_rows_retained"]=length(report["logical_rows"])
    report["raw_witnesses_retained"]=length(report["raw_witnesses"])
    report["unique_objectives"]=length(report["checks"])
    report["cost_scope"]="Driver CPU only; owning OS receipt includes startup/import/JIT/teardown"
    report["CPU_seconds"]=IOStudy.cpu_seconds()-began;report["wall_seconds"]=(time_ns()-wall)/1e9
    IOStudy.write_new(joinpath(output,"OBJECTIVE.toml"),report);report
end
function main(args=ARGS)
    length(args) in (3,4) || error("usage: objective.jl INPUT_MODEL.toml CANDIDATES.toml FRESH_OUTPUT_DIRECTORY [WITNESSES.toml]")
    BLAS.set_num_threads(1);Threads.nthreads()==1 || error("one Julia numerical thread required")
    result=run(args[1],args[2],args[3];witnesses_path=length(args)==4 ? args[4] : nothing)
    result["complete"] ? 0 : 2
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(GeneralizationObjective.main())
end
