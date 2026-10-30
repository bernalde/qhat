module GeneralizationIndependent

# Only the independent raw-Pauli implementation supplies numerical operations.
using TOML, SHA, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(ROOT, "numerical_audit", "prospective_window_ordering", "independent.jl"))
using .IndependentPhysical: load_model, reference_action, affordable, unfused_state,
    state_metrics, read_state

const MAX_BYTES = 8 * 1024^2
function artifact_path(root,relative)
    require(relative isa String && 1 <= ncodeunits(relative) <= 512 && !isabspath(relative) &&
        normpath(relative) == relative && !any(==(".."),splitpath(relative)),"strict relative artifact path")
    base=realpath(root); path=realpath(joinpath(base,relative))
    require(isfile(path) && !any(==(".."),splitpath(relpath(path,base))),"artifact symlink escapes DATA")
    path
end
function data_root(input,case)
    path=realpath(input); root=realpath(dirname(dirname(dirname(path))))
    require(path == joinpath(root,"inputs",case,"model.toml") && normpath(abspath(input)) == path,
        "canonical DATA/inputs/case/model.toml required")
    root
end
function observed_read(io,observed,check=require)
    check(observed isa Integer && !(observed isa Bool) && 0 <= observed <= MAX_BYTES,"artifact byte cap")
    raw=read(io,observed+1)
    check(length(raw) == observed && eof(io),"artifact grew or shrank during read")
    raw
end
function bounded_bytes(path,check=require)
    check(isfile(path),"bounded regular artifact required")
    observed=filesize(path)
    check(0 <= observed <= MAX_BYTES,"artifact byte cap")
    open(path,"r") do io
        observed_read(io,observed,check)
    end
end
function bounded_hash(path)
    bytes2hex(sha256(bounded_bytes(path)))
end


function bounded_toml(path)
    check=IndependentPhysical.require
    raw=bounded_bytes(path,check)
    # Admit only shallow scientific tables/arrays before invoking the parser.
    quote_char=UInt8(0);escaped=false;comment=false;nesting=UInt8[]
    for (index,byte) in enumerate(raw)
        if comment
            byte in (0x0a,0x0d) && (comment=false)
        elseif quote_char!=0
            # Basic and literal strings in this schema are single-line.
            # Check this before consuming an escaped character.
            check(!(byte in (0x0a,0x0d)),"raw newline in TOML string")
            if escaped
                escaped=false
            elseif quote_char==UInt8('"') && byte==UInt8('\\')
                escaped=true
            elseif byte==quote_char
                quote_char=0
            end
        elseif byte==UInt8('#')
            # Quote and bracket bytes in a line comment are not TOML syntax.
            comment=true
        elseif byte in (UInt8('"'),UInt8('\''))
            check(!(index+2<=length(raw) && raw[index+1]==byte && raw[index+2]==byte),
                "multiline TOML strings are outside schema")
            quote_char=byte
        elseif byte==UInt8('[') || byte==UInt8('{')
            check(length(nesting)<8,"TOML nesting admission bound")
            push!(nesting,byte)
        elseif byte==UInt8(']') || byte==UInt8('}')
            opening=byte==UInt8(']') ? UInt8('[') : UInt8('{')
            check(!isempty(nesting) && last(nesting)==opening,"malformed TOML nesting")
            pop!(nesting)
        end
    end
    check(quote_char==0 && !escaped,"unterminated TOML string")
    check(isempty(nesting),"unterminated TOML nesting")
    (data=TOML.parse(String(copy(raw))),sha256=bytes2hex(sha256(raw)))
end
function bind_artifacts(input,manifest)
    root=data_root(input,manifest["case"])
    require(get(manifest,"protocol_path",nothing) == "source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json" &&
        get(manifest,"source_freeze_path",nothing) == "INPUT_PROTOCOL_FROZEN.json","explicit frozen protocol/source artifact paths")
    hashes=Dict{String,String}()
    for (pathkey,hashkey) in (("protocol_path","protocol_sha256"),("source_freeze_path","source_freeze_sha256"))
        relative=manifest[pathkey]; digest=bounded_hash(artifact_path(root,relative))
        require(digest == manifest[hashkey],"frozen protocol/source actual bytes: "*pathkey)
        hashes[relative]=digest
    end
    jsonpath="inputs/"*manifest["case"]*"/model.json"
    hashes[jsonpath]=bounded_hash(artifact_path(root,jsonpath))
    require(hashes[jsonpath] == manifest["source_json_sha256"],"actual source model JSON bytes")
    artifacts=get(manifest,"constructor_artifacts",nothing)
    require(artifacts isa Vector,"complete constructor artifact inventory required")
    representatives=Dict(row["order_id"]=>row for row in manifest["candidates"])
    require(length(artifacts) == length(representatives),"one actual artifact per constructor")
    seen=Set{String}()
    for artifact in artifacts
        require(artifact isa AbstractDict && Set(keys(artifact)) == Set(("order_id","path","sha256")),"constructor artifact fields")
        oid=artifact["order_id"]
        require(oid isa String && haskey(representatives,oid) && !(oid in seen),"unique constructor artifact order_id")
        push!(seen,oid); row=representatives[oid]
        expected="constructors/"*manifest["case"]*".rep"*string(row["repetition"])*"."*row["arm"]*".B"*string(row["cpu_budget"])*".json"
        require(artifact["path"] == expected && artifact["sha256"] == row["constructor_sha256"],"canonical constructor path/digest identity")
        digest=bounded_hash(artifact_path(root,artifact["path"]))
        require(digest == artifact["sha256"],"actual constructor bytes differ")
        hashes[artifact["path"]]=digest
    end
    require(seen == Set(keys(representatives)),"missing constructor artifact")
    (data_root=root,hashes=hashes,status="actual_bytes_bound",
        scope="Julia binds exact model/candidate TOML and actual protocol/freeze/modelJSON/constructorJSON bytes. Independent Python consumer separately checks raw JSON semantic order/score/identity joins.")
end

const TIMES = (.65, 1.35)
const CAPS = (8, 24)
const STATE_TOLERANCE = 1e-10
const METRIC_TOLERANCE = 1e-12
const FIELDS = ("name", "order_id", "arm", "repetition", "cpu_budget", "partition", "formula",
    "time", "budget_multiplier", "order", "constructor_sha256", "claimed_score", "score_kind", "eligible")
const METRICS = ("projective_error", "infidelity", "norm_drift", "number_leakage", "spin_leakage", "reference_norm_drift")
require(ok, message) = ok || throw(ArgumentError(message))
integer(value) = value isa Integer && !(value isa Bool)
hexhash(value) = value isa String && occursin(r"^[0-9a-f]{64}$", value)
label(value) = value isa String && 1 <= ncodeunits(value) <= 200 && all(c -> '!' <= c <= '~', value)
filehash(path) = bytes2hex(sha256(read(path)))
failure(error) = Dict("type" => string(typeof(error)), "message" => sprint(showerror, error))

function write_new(path, record)
    temporary = path * ".partial"
    require(!ispath(path) && !ispath(temporary), "refuse existing result/partial")
    open(temporary, "w") do io; TOML.print(io, record; sorted=true); end
    mv(temporary, path; force=false)
end

function expected_primary_sources()
    names = ("numerical_audit/molecular_generalization/physical.jl",
        "numerical_audit/molecular_generalization/src/GeneralizationPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
        "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl")
    Dict(name => filehash(joinpath(ROOT, name)) for name in names)
end

function own_sources()
    names = ("numerical_audit/molecular_generalization/verify.jl",
        "numerical_audit/molecular_generalization/src/GeneralizationIndependent.jl",
        "numerical_audit/prospective_window_ordering/independent.jl")
    Dict(name => filehash(joinpath(ROOT, name)) for name in names)
end


const CONSTRUCTORS = (("parent",0),("canonical",0),("magnitude",0),("full",32),("full",64),
    ("window",32),("window",64),("annealing",32),("annealing",64))
const SCORE_KIND = "diagonal_reference_cubic_response"
const HEADER = ("version","order_index_base","case","n","electrons","partition","input_sha256",
    "source_json_sha256","protocol_sha256","source_freeze_sha256","mode","repetitions",
    "expected_logical_rows","training_time","training_budget_multiplier","times","rotation_cap_multipliers")
number(value) = value isa Real && !(value isa Bool) && isfinite(value) && value >= 0
export run, admit_candidates, admit_manifest, expected_grid
function case_specification(mode, case)
    if mode == "preflight"
        require(case == "LiH_r1p487_631g_14m4e", "only the declared old-case preflight")
        return (14,4)
    end
    require(mode == "main" && case in ("HF_r1p00_631g_12m8e",
        "H2S_r1p40_a96_631g_12m8e","HCl_r1p30_631g_12m8e"), "unknown frozen main case")
    (12,8)
end

function expected_grid(mode)
    require(mode in ("main","preflight"), "unknown frozen mode")
    reps = mode == "main" ? 2 : 1
    Set((arm,rep,cpu,t,cap) for rep in 1:reps for (arm,cpu) in CONSTRUCTORS for t in TIMES for cap in CAPS)
end
logical_identity(row) = (row["arm"],row["repetition"],row["cpu_budget"],row["time"],row["budget_multiplier"])
constructor_identity(row) = (row["arm"],row["repetition"],row["cpu_budget"])
signature(row) = (Tuple(row[k] for k in FIELDS if !(k in ("name","time","budget_multiplier"))),
    get(row,"parent_sha256",nothing),get(row,"constructor_score",nothing),get(row,"constructor_failure",nothing))

"""Pure structural roster admission, usable on tiny fixtures without propagation.

The model identity is supplied separately by the bounded model admission. Counts
are derived from the fixed grid, never from an asserted expected-count field.
"""
function admit_manifest(m, unit_count)
    require(m isa AbstractDict && length(m) <= 32 && all(k -> haskey(m,k),HEADER), "complete bounded candidate header")
    require(m["version"] === 1 && m["order_index_base"] === 1 && m["partition"] == "conserving", "candidate boundary version/index/partition")
    require(label(m["case"]) && integer(m["n"]) && integer(m["electrons"]), "candidate model identity")
    require(1 <= m["n"] <= 14 && 0 <= m["electrons"] <= m["n"] &&
        binomial(m["n"],m["electrons"]) <= 2000, "bounded full-register/fixed-N domain")
    for key in ("input_sha256","source_json_sha256","protocol_sha256","source_freeze_sha256")
        require(hexhash(m[key]), "candidate SHA-256 spelling: "*key)
    end
    expected = expected_grid(m["mode"]); reps = m["mode"] == "main" ? 2 : 1
    require(m["repetitions"] === reps && m["expected_logical_rows"] === 36reps, "candidate roster count")
    require(m["training_time"] isa AbstractFloat && m["training_time"] == 1.0 &&
        m["training_budget_multiplier"] === 30 && m["times"] == collect(TIMES) &&
        all(v -> v isa AbstractFloat,m["times"]) && m["rotation_cap_multipliers"] == collect(CAPS) &&
        all(integer,m["rotation_cap_multipliers"]), "fixed training and physical grids")
    rows = m["candidates"]
    require(rows isa Vector && length(rows) == length(expected), "complete candidate rows")
    observed=Set{Any}(); names=Set{String}(); ids=Dict{String,Any}(); counts=Dict{String,Int}(); constructors=Dict{Any,String}()
    for row in rows
        require(row isa AbstractDict && length(row) <= 28 && all(k -> haskey(row,k),FIELDS), "candidate row fields")
        require(label(row["name"]) && label(row["order_id"]) && !(row["name"] in names), "unique bounded ASCII candidate names")
        push!(names,row["name"])
        require(all(k -> integer(row[k]),("repetition","cpu_budget","formula","budget_multiplier")), "integer row identity")
        require(row["time"] isa AbstractFloat && isfinite(row["time"]) && row["formula"] === 2 &&
            row["partition"] == "conserving" && row["eligible"] isa Bool, "physical formula/time/partition/eligibility")
        key=logical_identity(row); require(key in expected && !(key in observed), "unknown/duplicate physical identity")
        push!(observed,key)
        order=row["order"]
        require(order isa Vector && length(order) == unit_count && all(integer,order) &&
            sort(order) == collect(1:unit_count), "whole original one-based order")
        require(hexhash(row["constructor_sha256"]) && number(row["claimed_score"]) && row["score_kind"] == SCORE_KIND,
            "finite scored constructor identity")
        if haskey(row,"constructor_score")
            require(number(row["constructor_score"]) && row["constructor_score"] == row["claimed_score"], "equal optional score alias")
        end
        if row["arm"] in ("canonical","magnitude")
            require(!haskey(row,"parent_sha256"), "cheap construction has no measured parent prefix")
        else
            require(haskey(row,"parent_sha256") && hexhash(row["parent_sha256"]), "parent identity required")
        end
        oid=row["order_id"]; sig=signature(row); cid=constructor_identity(row)
        require(!haskey(ids,oid) || ids[oid] == sig, "order_id changes constructor or order")
        require(!haskey(constructors,cid) || constructors[cid] == oid, "one order_id per constructor")
        ids[oid]=sig; constructors[cid]=oid; counts[oid]=get(counts,oid,0)+1
    end
    require(observed == expected && length(ids) == 9reps && all(==(4),values(counts)), "complete constructor graph")
    for rep in 1:reps
        parent=only(filter(row -> row["arm"] == "parent" && row["repetition"] == rep &&
            row["time"] == first(TIMES) && row["budget_multiplier"] == first(CAPS),rows))
        require(parent["parent_sha256"] == parent["constructor_sha256"], "parent self-reference")
        require(all(row -> row["repetition"] != rep || row["arm"] in ("canonical","magnitude") ||
            row["parent_sha256"] == parent["constructor_sha256"],rows), "optimized common measured parent")
    end
    m
end

function admit_candidates(model,candidates_path;input)
    require(model.partition == "conserving" && binomial(model.n,model.electrons) <= 2000,"independent conserving model/fixed-N cap")
    parsed=bounded_toml(candidates_path)
    m=admit_manifest(parsed.data,length(model.units))
    require(m["input_sha256"] == model.input_sha256 && m["source_json_sha256"] == model.source_json_sha256 &&
        m["case"] == model.case && m["n"] == model.n && m["electrons"] == model.electrons,"independent candidate/model byte identity")
    require((model.n,model.electrons) == case_specification(m["mode"],model.case),"independent frozen case specification")
    bindings=bind_artifacts(input,m)
    require(realpath(candidates_path) == joinpath(bindings.data_root,"panels",m["case"]*".toml") &&
        normpath(abspath(candidates_path)) == realpath(candidates_path),"canonical DATA/panels/case.toml required")
    (data=m,sha256=parsed.sha256,bindings=bindings)
end


function numeric_metric(row, key)
    value = row[key]
    require(value isa Real && !(value isa Bool) && isfinite(value) && value >= 0, "finite nonnegative physical metric: " * key)
    Float64(value)
end

function compare_metrics(row, metrics)
    discrepancies = Dict{String,Float64}()
    for key in METRICS
        difference = abs(numeric_metric(row, key) - metrics[key])
        discrepancies[key] = difference
        require(difference <= METRIC_TOLERANCE, "physical metric mismatch: " * key)
    end
    discrepancies
end

function run(input, candidates_path, primary_dir, output; cpu_limit=600.0, wall_limit=3000.0, check_budget=()->nothing)
    require(!ispath(output), "fresh verification directory required"); mkdir(output)
    began = time_ns(); cpu_began = ccall(:clock,Int64,()) / 1e6
    require(number(cpu_limit) && cpu_limit > 0 && number(wall_limit) && wall_limit > 0,"finite positive verifier allowance")
    guard=() -> begin
        check_budget()
        require(ccall(:clock,Int64,()) / 1e6-cpu_began <= cpu_limit,"independent CPU allowance exhausted")
        require((time_ns()-began)/1e9 <= wall_limit,"independent wall allowance exhausted")
    end
    report = Dict{String,Any}("schema" => "molecular-generalization-independent-v1", "status" => "incomplete",
        "rows" => Any[], "failures" => Any[], "censored_rows" => Any[], "reference_checks" => Any[],
        "checker_source_hashes" => own_sources(), "completed_logical_rows" => 0,
        "unique_circuit_count" => 0, "independent_unfused_replay_count" => 0, "reference_check_count" => 0,
        "state_tolerance" => STATE_TOLERANCE, "absolute_metric_tolerance" => METRIC_TOLERANCE,
        "replay_scope" => "Fresh unfused full-register replay for every distinct saved circuit, and independently reconstructed references at both held-out times",
        "objective_scope" => "Candidate score metadata authenticated here; separate objective replay verifies numerical claimed scores",
        "exact_polynomial_conservation" => "Independent integer [N,H]=[Nalpha,H]=0 and original-unit conservation checks",
        "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads())
    try
        guard(); model = load_model(input); candidate = admit_candidates(model, candidates_path;input=input); candidates = candidate.data
        # Candidate identities remain visible even if the primary is incomplete.
        for row in candidates["candidates"]
            entry=Dict{String,Any}(key=>row[key] for key in FIELDS)
            for key in ("parent_sha256","constructor_score","constructor_failure")
                haskey(row,key) && (entry[key]=row[key])
            end
            entry["status"]="incomplete"; push!(report["rows"],entry)
        end
        entries_by_name=Dict(entry["name"]=>entry for entry in report["rows"])
        primary_path = joinpath(primary_dir, "PHYSICAL.toml")
        parsed = bounded_toml(primary_path); primary = parsed.data
        report["input_sha256"] = model.input_sha256; report["source_json_sha256"] = model.source_json_sha256
        report["candidates_sha256"] = candidate.sha256; report["physical_sha256"] = parsed.sha256
        report["case"] = model.case
        require(primary["artifact_binding_status"] == candidate.bindings.status && primary["bound_artifact_hashes"] == candidate.bindings.hashes,"independent actual artifact byte joins")
        report["artifact_binding_status"]=candidate.bindings.status
        report["artifact_binding_scope"]=candidate.bindings.scope
        report["bound_artifact_hashes"]=candidate.bindings.hashes
        for key in ("input_sha256", "source_json_sha256", "protocol_sha256", "source_freeze_sha256",
                "case", "mode", "repetitions", "expected_logical_rows", "n", "electrons", "partition", "order_index_base", "training_time", "training_budget_multiplier", "times", "rotation_cap_multipliers")
            require(primary[key] == candidates[key], "primary/candidate header mismatch: " * key)
            report[key] = candidates[key]
        end
        require(primary["candidates_sha256"] == candidate.sha256, "primary exact candidate bytes")
        require(primary["source_hashes"] == expected_primary_sources(), "primary direct/dependency source bytes differ")
        report["primary_source_hashes"] = primary["source_hashes"]
        require(primary["n"] == model.n && primary["electrons"] == model.electrons && primary["hf"] == string(model.hf) && primary["full_dimension"] == 1<<model.n, "primary full-register/HF identity")
        require(primary["status"] == "complete" && primary["physical_panel_complete"] === true, "primary panel incomplete/censored")
        require(primary["reference_sector_dimension"] == binomial(model.n,model.electrons),"complete fixed-N reference, never fixed-spin projection")
        coupling_tolerance=numeric_metric(primary,"reference_coupling_tolerance")
        require(numeric_metric(primary,"reference_outside_coupling_max") <= coupling_tolerance &&
            numeric_metric(primary,"reference_hermiticity_residual") <= coupling_tolerance,"saved reference sector/Hermiticity admission")
        require(primary["all_constructor_eligible"] === all(row -> row["eligible"],candidates["candidates"]),"primary constructor eligibility summary")
        rows = primary["rows"]; roster = candidates["candidates"]
        require(rows isa Vector && length(rows) == length(roster) && primary["completed_logical_rows"] == length(roster), "no missing/extra primary rows")
        byname = Dict(row["name"] => row for row in roster)
        require(all(row -> row isa AbstractDict && length(row) <= 64 && haskey(row, "name") && label(row["name"]), rows), "bounded primary rows")
        require(length(Set(row["name"] for row in rows)) == length(rows) && Set(row["name"] for row in rows) == Set(keys(byname)), "no missing/extra row names")
        for row in rows
            original = byname[row["name"]]
            require(get(row,"status",nothing) == "complete","primary logical row is incomplete")
            require(all(key -> haskey(row, key) && row[key] == original[key], FIELDS), "primary row changes candidate constructor graph")
            for key in ("parent_sha256","constructor_score","constructor_failure")
                require(haskey(row,key)==haskey(original,key) && (!haskey(row,key) || row[key]==original[key]),"optional candidate/primary join: "*key)
            end
        end
        for row in rows
            entry=entries_by_name[row["name"]]
            for key in ("rotation_budget","steps","ideal_rotations","reference_file","reference_sha256","actual_file","actual_sha256","circuit_sha256")
                haskey(row,key) && (entry[key]=row[key])
            end
        end
        savedrefs = primary["reference_records"]
        require(savedrefs isa Vector && length(savedrefs) == 2 && Set(row["time"] for row in savedrefs) == Set(TIMES), "complete unique saved references")
        refmeta = Dict(Float64(row["time"]) => row for row in savedrefs)
        references = Dict{Float64,Any}(); savedreferences = Dict{Float64,Any}()
        append!(report["reference_checks"],(Dict{String,Any}("time"=>duration,"status"=>"incomplete",
            "reference_file"=>refmeta[duration]["reference_file"],"reference_sha256"=>refmeta[duration]["reference_sha256"]) for duration in TIMES))
        for (reference_index,duration) in enumerate(TIMES)
            rec=report["reference_checks"][reference_index]
            try
                guard(); independent = reference_action(model, duration;check_budget=guard)
                require(independent.diagnostic["sector_dimension"] == binomial(model.n,model.electrons),"independent complete fixed-N basis")
                saved = refmeta[duration]
                require(get(saved,"status",nothing) == "complete" && integer(saved["sector_dimension"]) &&
                    saved["sector_dimension"] == binomial(model.n,model.electrons),"saved reference complete fixed-N dimension/status")
                reference = read_state(primary_dir, saved["reference_file"], saved["reference_sha256"], model.n)
                residual = norm(reference - independent.state)
                aligned = IndependentPhysical.phase_aligned_residual(independent.state, reference)
                merge!(rec,Dict("state_norm_residual"=>residual,"phase_aligned_residual"=>aligned,"diagnostic"=>independent.diagnostic))
                require(residual <= STATE_TOLERANCE && aligned <= STATE_TOLERANCE, "independent reference mismatch")
                metrics = state_metrics(reference, reference, model.n, model.hf)
                for key in ("norm_drift", "number_leakage", "spin_leakage")
                    require(abs(numeric_metric(saved, key) - metrics[key]) <= METRIC_TOLERANCE, "saved reference metric mismatch")
                end
                references[duration] = independent; savedreferences[duration] = reference; rec["status"] = "complete"
            catch error
                rec["failure"]=failure(error)
                push!(report["failures"], merge(Dict("phase" => "reference", "time" => duration), failure(error)))
            end
        end
        replay_cache = Dict{Any,Any}(); replay_failures = Dict{Any,Any}()
        actual_graph = Dict{Any,Any}(); file_graph = Dict{String,Any}(); reference_files = Set{String}()
        for row in rows
            guard()
            entry = entries_by_name[row["name"]]
            try
                duration = Float64(row["time"]); order = row["order"]; budget = row["budget_multiplier"] * length(model.words)
                require(integer(row["rotation_budget"]) && row["rotation_budget"] == budget, "original-word rotation budget")
                plan = affordable(model, "conserving", order, budget; formula=2)
                require(integer(row["steps"]) && row["steps"] == plan.steps && integer(row["ideal_rotations"]) &&
                    row["ideal_rotations"] == plan.ideal_rotations, "independent steps/rotation accounting")
                saved = refmeta[duration]
                require(row["reference_file"] == saved["reference_file"] && row["reference_sha256"] == saved["reference_sha256"], "row reference identity")
                push!(reference_files, row["reference_file"])
                key = (Tuple(order), 2, reinterpret(UInt64, duration), plan.steps)
                artifact = (row["actual_file"], row["actual_sha256"])
                duplicate = haskey(actual_graph, key)
                require(!duplicate || actual_graph[key] == artifact, "identical circuit maps to different saved artifact")
                require(!haskey(file_graph, row["actual_file"]) || file_graph[row["actual_file"]] == key, "saved artifact reused for different circuit")
                require(row["propagation_reused_for_identical_identity"] === duplicate, "primary circuit deduplication flag")
                actual_graph[key] = artifact; file_graph[row["actual_file"]] = key
                require(haskey(references, duration) && haskey(savedreferences, duration), "reference check incomplete")
                actual = read_state(primary_dir, row["actual_file"], row["actual_sha256"], model.n)
                reference = savedreferences[duration]; independent = references[duration].state
                metrics = state_metrics(reference, actual, model.n, model.hf)
                independent_metrics = state_metrics(independent, actual, model.n, model.hf)
                require(metrics["number_leakage"] <= 1e-10 && metrics["spin_leakage"] <= 1e-10,"conserving saved circuit leakage failure")
                rawdeficit=1-abs(dot(reference,actual))/sqrt(sum(abs2,reference)*sum(abs2,actual))
                require(row["raw_normalized_overlap_deficit"] isa Real && isfinite(row["raw_normalized_overlap_deficit"]) && abs(row["raw_normalized_overlap_deficit"]-rawdeficit)<=METRIC_TOLERANCE,"raw normalized overlap deficit mismatch")
                require(row["numerical_floor"] === (numeric_metric(row,"projective_error")<=1e-12),"numerical floor label mismatch")
                circuit_hash=bytes2hex(sha256(model.input_sha256*"|"*string(key)))
                require(row["circuit_sha256"] == circuit_hash,"exact circuit identity digest mismatch")
                entry["saved_reference_metrics"] = metrics; entry["independent_reference_metrics"] = independent_metrics
                entry["saved_metric_discrepancies"] = compare_metrics(row, metrics)
                entry["independent_reference_metric_discrepancies"] = compare_metrics(row, independent_metrics)
                overlap = min(1., abs(dot(reference, actual)) / sqrt(sum(abs2, reference) * sum(abs2, actual)))
                require(abs(numeric_metric(row, "overlap_abs") - overlap) <= METRIC_TOLERANCE, "normalized overlap mismatch")
                refmetrics = state_metrics(reference, reference, model.n, model.hf)
                for (savedkey, metrickey) in (("reference_number_leakage", "number_leakage"), ("reference_spin_leakage", "spin_leakage"))
                    require(abs(numeric_metric(row, savedkey) - refmetrics[metrickey]) <= METRIC_TOLERANCE, "reference leakage mismatch")
                end
                if !haskey(replay_cache, key) && !haskey(replay_failures, key)
                    try
                        replay_cache[key] = unfused_state(model, "conserving", order, duration, plan.steps; formula=2,check_budget=guard)
                    catch error
                        replay_failures[key] = error
                    end
                end
                haskey(replay_failures, key) && throw(replay_failures[key])
                replay = replay_cache[key]
                residual = norm(replay.state - actual)
                aligned = IndependentPhysical.phase_aligned_residual(actual, replay.state)
                require(residual <= STATE_TOLERANCE && aligned <= STATE_TOLERANCE, "fresh independent unfused circuit mismatch")
                replay_metrics = state_metrics(reference, replay.state, model.n, model.hf)
                independent_replay_metrics = state_metrics(independent, replay.state, model.n, model.hf)
                entry["unfused_saved_reference_metric_discrepancies"] = compare_metrics(row, replay_metrics)
                entry["unfused_independent_reference_metric_discrepancies"] = compare_metrics(row, independent_replay_metrics)
                entry["unfused_state_norm_residual"] = residual; entry["unfused_phase_aligned_residual"] = aligned
                entry["auxiliary_unfused_primitive_pulses"] = replay.primitive_pulses
                entry["actual_sha256"] = row["actual_sha256"]; entry["independent_unfused_replay"] = true
                entry["status"] = "complete"
            catch error
                entry["failure"] = failure(error)
                push!(report["censored_rows"], Dict("name" => row["name"], "reason" => entry["failure"]))
            end
        end
        report["unique_circuit_count"] = length(actual_graph)
        report["independent_unfused_replay_count"] = length(replay_cache)
        report["reference_check_count"] = length(savedreferences)
        require(primary["unique_actual_propagations"] == length(actual_graph), "primary unique circuit count mismatch")
        require(length(reference_files) == 2 && isempty(intersect(reference_files, Set(keys(file_graph)))), "distinct reference/actual file identities")
        expected_files = union(reference_files, Set(keys(file_graph)))
        observed_files = Set(name for name in readdir(primary_dir) if endswith(name, ".bin"))
        require(observed_files == expected_files, "unreferenced or missing saved-state file")
        report["completed_logical_rows"] = count(row -> row["status"] == "complete", report["rows"])
        if isempty(report["failures"]) && isempty(report["censored_rows"]) && length(savedreferences) == 2 &&
                length(replay_cache) == length(actual_graph) && report["completed_logical_rows"] == length(roster)
            report["status"] = "complete"
        end
    catch error
        push!(report["failures"], merge(Dict("phase" => "inventory_or_source_binding"), failure(error)))
    end
    report["completed_logical_rows"] = count(row -> row["status"] == "complete", report["rows"])
    report["wall_seconds"] = (time_ns() - began) / 1e9
    report["CPU_seconds"] = ccall(:clock,Int64,()) / 1e6 - cpu_began
    report["logical_rows_retained"] = length(report["rows"])
    report["reference_checks_retained"] = length(report["reference_checks"])
    report["all_constructor_eligible"] = !isempty(report["rows"]) && all(row -> get(row,"eligible",false),report["rows"])
    report["cost_scope"] = "Driver CPU only; owning OS receipt includes startup/import/JIT/teardown"
    report["cooperation_scope"] = "Budget callbacks inside unchanged independent reference/unfused kernels; owning OS guard is authoritative"

    write_new(joinpath(output, "VERIFY.toml"), report)
    report
end

end
