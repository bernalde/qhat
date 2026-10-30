module CostMatchedVerify

# Only the independent raw-Pauli implementation supplies numerical operations.
using TOML, SHA, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "numerical_audit", "prospective_window_ordering", "independent.jl"))
using .IndependentPhysical: load_model, reference_action, affordable, unfused_state,
    state_metrics, read_state

const TIMES = (.65, 1.35)
const CAPS = (8, 24)
const STATE_TOLERANCE = 1e-10
const METRIC_TOLERANCE = 1e-12
const FIELDS = ("name", "order_id", "arm", "repetition", "cpu_budget", "partition", "formula",
    "time", "budget_multiplier", "order", "constructor_sha256", "parent_sha256", "claimed_score", "score_kind")
const METRICS = ("projective_error", "infidelity", "norm_drift", "number_leakage", "spin_leakage", "reference_norm_drift")
require(ok, message) = ok || throw(ArgumentError(message))
integer(value) = value isa Integer && !(value isa Bool)
hexhash(value) = value isa String && occursin(r"^[0-9a-f]{64}$", value)
label(value) = value isa String && 1 <= ncodeunits(value) <= 200 && all(c -> ' ' <= c <= '~', value)
filehash(path) = bytes2hex(sha256(read(path)))
failure(error) = Dict("type" => string(typeof(error)), "message" => sprint(showerror, error))

function write_new(path, record)
    temporary = path * ".partial"
    require(!ispath(path) && !ispath(temporary), "refuse existing result/partial")
    open(temporary, "w") do io; TOML.print(io, record; sorted=true); end
    mv(temporary, path; force=false)
end

function expected_primary_sources()
    names = ("numerical_audit/cost_matched_continuation/physical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
        "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl")
    Dict(name => filehash(joinpath(ROOT, name)) for name in names)
end

function own_sources()
    names = ("numerical_audit/cost_matched_continuation/verify.jl",
        "numerical_audit/prospective_window_ordering/independent.jl")
    Dict(name => filehash(joinpath(ROOT, name)) for name in names)
end

logical_identity(row) = (row["arm"], row["repetition"], row["cpu_budget"], row["time"], row["budget_multiplier"])
signature(row) = Tuple(row[key] for key in FIELDS if !(key in ("name", "time", "budget_multiplier")))

function admit_candidates(model, candidates_path)
    parsed = IndependentPhysical.bounded_toml(candidates_path); data = parsed.data
    require(length(data) <= 24 && data["version"] === 1, "candidate header version/size")
    require(model.n == 14 && model.electrons == 4 && model.partition == "conserving", "frozen original model")
    for key in ("input_sha256", "source_json_sha256", "protocol_sha256", "source_freeze_sha256")
        require(hexhash(data[key]), "candidate SHA-256 spelling: " * key)
    end
    require(data["input_sha256"] == model.input_sha256 && data["source_json_sha256"] == model.source_json_sha256 &&
        data["case"] == model.case, "candidate input/JSON/case identity")
    mode = data["mode"]; require(mode in ("main", "preflight"), "explicit frozen mode")
    reps = mode == "main" ? 2 : 1
    require(data["repetitions"] === reps && data["expected_logical_rows"] === 20reps, "complete logical inventory count")
    expected = Set((arm, rep, cpu, t, cap) for rep in 1:reps
        for (arm, cpu) in (("parent", 0), ("full", 32), ("full", 64), ("window", 32), ("window", 64))
        for t in TIMES for cap in CAPS)
    rows = data["candidates"]
    require(rows isa Vector && length(rows) == length(expected), "complete bounded candidate rows")
    names = Set{String}(); observed = Set{Any}(); order_ids = Dict{String,Any}(); counts = Dict{String,Int}()
    for row in rows
        require(row isa AbstractDict && length(row) <= 24 && all(key -> haskey(row, key), FIELDS), "candidate fields")
        require(label(row["name"]) && label(row["order_id"]) && !(row["name"] in names), "unique bounded ASCII name")
        push!(names, row["name"])
        require(all(key -> integer(row[key]), ("repetition", "cpu_budget", "formula", "budget_multiplier")), "integer row identity")
        require(row["time"] isa AbstractFloat && isfinite(row["time"]) && row["partition"] == "conserving" && row["formula"] == 2,
            "frozen time/partition/formula")
        key = logical_identity(row); require(key in expected && !(key in observed), "unknown or duplicate logical identity")
        push!(observed, key)
        IndependentPhysical.validate_order(row["order"], length(model.units))
        require(hexhash(row["constructor_sha256"]) && hexhash(row["parent_sha256"]), "constructor and parent hash spelling")
        require(row["claimed_score"] isa AbstractFloat && isfinite(row["claimed_score"]) && row["claimed_score"] >= 0 &&
            row["score_kind"] == "diagonal_reference_cubic_response", "finite objective witness")
        require(haskey(row, "constructor_score") && row["constructor_score"] isa Real &&
            !(row["constructor_score"] isa Bool) && isfinite(row["constructor_score"]) &&
            row["constructor_score"] == row["claimed_score"], "matching constructor score aliases")
        oid = row["order_id"]; sig = signature(row)
        require(!haskey(order_ids, oid) || order_ids[oid] == sig, "order_id changes constructor/order identity")
        order_ids[oid] = sig; counts[oid] = get(counts, oid, 0) + 1
    end
    require(observed == expected && length(order_ids) == 5reps && all(==(4), values(counts)), "complete constructor graph")
    for rep in 1:reps
        parent = only(filter(row -> row["arm"] == "parent" && row["repetition"] == rep &&
            row["time"] == first(TIMES) && row["budget_multiplier"] == first(CAPS), rows))
        require(parent["parent_sha256"] == parent["constructor_sha256"], "parent self identity")
        require(all(row -> row["repetition"] != rep || row["parent_sha256"] == parent["constructor_sha256"], rows),
            "continuation frozen-parent identity differs")
    end
    (data=data, sha256=parsed.sha256)
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

function run_check(input, candidates_path, primary_dir, output)
    require(!ispath(output), "fresh verification directory required"); mkdir(output)
    began = time_ns()
    report = Dict{String,Any}("schema" => "cost-matched-continuation-independent-v1", "status" => "incomplete",
        "rows" => Any[], "failures" => Any[], "censored_rows" => Any[], "reference_checks" => Any[],
        "checker_source_hashes" => own_sources(), "completed_logical_rows" => 0,
        "unique_circuit_count" => 0, "independent_unfused_replay_count" => 0, "reference_check_count" => 0,
        "state_tolerance" => STATE_TOLERANCE, "absolute_metric_tolerance" => METRIC_TOLERANCE,
        "replay_scope" => "Fresh unfused full-register replay for every distinct saved circuit, and independently reconstructed references at both held-out times",
        "objective_scope" => "Candidate score metadata authenticated here; separate objective replay verifies numerical claimed scores",
        "exact_polynomial_conservation" => "Independent integer [N,H]=[Nalpha,H]=0 and original-unit conservation checks",
        "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads())
    try
        model = load_model(input); candidate = admit_candidates(model, candidates_path); candidates = candidate.data
        primary_path = joinpath(primary_dir, "PHYSICAL.toml")
        parsed = IndependentPhysical.bounded_toml(primary_path); primary = parsed.data
        report["input_sha256"] = model.input_sha256; report["source_json_sha256"] = model.source_json_sha256
        report["candidates_sha256"] = candidate.sha256; report["physical_sha256"] = parsed.sha256
        report["case"] = model.case
        for key in ("input_sha256", "source_json_sha256", "protocol_sha256", "source_freeze_sha256",
                "case", "mode", "repetitions", "expected_logical_rows")
            require(primary[key] == candidates[key], "primary/candidate header mismatch: " * key)
            report[key] = candidates[key]
        end
        require(primary["candidates_sha256"] == candidate.sha256, "primary exact candidate bytes")
        require(primary["source_hashes"] == expected_primary_sources(), "primary direct/dependency source bytes differ")
        report["primary_source_hashes"] = primary["source_hashes"]
        require(primary["n"] === 14 && primary["hf"] == "15", "primary full-register/HF identity")
        require(primary["status"] == "complete" && primary["physical_panel_complete"] === true, "primary panel incomplete/censored")
        rows = primary["rows"]; roster = candidates["candidates"]
        require(rows isa Vector && length(rows) == length(roster) && primary["completed_logical_rows"] == length(roster), "no missing/extra primary rows")
        byname = Dict(row["name"] => row for row in roster)
        require(all(row -> row isa AbstractDict && length(row) <= 64 && haskey(row, "name") && label(row["name"]), rows), "bounded primary rows")
        require(length(Set(row["name"] for row in rows)) == length(rows) && Set(row["name"] for row in rows) == Set(keys(byname)), "no missing/extra row names")
        for row in rows
            original = byname[row["name"]]
            require(all(key -> haskey(row, key) && row[key] == original[key], FIELDS), "primary row changes candidate constructor graph")
        end
        savedrefs = primary["reference_records"]
        require(savedrefs isa Vector && length(savedrefs) == 2 && Set(row["time"] for row in savedrefs) == Set(TIMES), "complete unique saved references")
        refmeta = Dict(Float64(row["time"]) => row for row in savedrefs)
        references = Dict{Float64,Any}(); savedreferences = Dict{Float64,Any}()
        for duration in TIMES
            try
                independent = reference_action(model, duration)
                saved = refmeta[duration]
                reference = read_state(primary_dir, saved["reference_file"], saved["reference_sha256"], model.n)
                residual = norm(reference - independent.state)
                aligned = IndependentPhysical.phase_aligned_residual(independent.state, reference)
                rec = Dict{String,Any}("time" => duration, "reference_file" => saved["reference_file"],
                    "reference_sha256" => saved["reference_sha256"], "state_norm_residual" => residual,
                    "phase_aligned_residual" => aligned, "diagnostic" => independent.diagnostic, "status" => "incomplete")
                push!(report["reference_checks"], rec)
                require(residual <= STATE_TOLERANCE && aligned <= STATE_TOLERANCE, "independent reference mismatch")
                metrics = state_metrics(reference, reference, model.n, model.hf)
                for key in ("norm_drift", "number_leakage", "spin_leakage")
                    require(abs(numeric_metric(saved, key) - metrics[key]) <= METRIC_TOLERANCE, "saved reference metric mismatch")
                end
                references[duration] = independent; savedreferences[duration] = reference; rec["status"] = "complete"
            catch error
                push!(report["failures"], merge(Dict("phase" => "reference", "time" => duration), failure(error)))
            end
        end
        replay_cache = Dict{Any,Any}(); replay_failures = Dict{Any,Any}()
        actual_graph = Dict{Any,Any}(); file_graph = Dict{String,Any}(); reference_files = Set{String}()
        for row in rows
            entry = Dict{String,Any}("name" => row["name"], "order_id" => row["order_id"], "arm" => row["arm"],
                "repetition" => row["repetition"], "cpu_budget" => row["cpu_budget"], "time" => row["time"],
                "budget_multiplier" => row["budget_multiplier"], "status" => "incomplete")
            push!(report["rows"], entry)
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
                        replay_cache[key] = unfused_state(model, "conserving", order, duration, plan.steps; formula=2)
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
    write_new(joinpath(output, "VERIFY.toml"), report)
    report
end

function main(args=ARGS)
    length(args) == 4 || error("usage: verify.jl INPUT_MODEL.toml CANDIDATES.toml PRIMARY_OUTPUT_DIRECTORY FRESH_OUTPUT_DIRECTORY")
    BLAS.set_num_threads(1); require(Threads.nthreads() == 1, "one Julia numerical thread required")
    run_check(args...)["status"] == "complete" ? 0 : 2
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(CostMatchedVerify.main())
end
