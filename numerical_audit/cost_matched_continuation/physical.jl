module CostMatchedPhysical

# New inventory driver; all numerical kernels remain in the frozen study.
using TOML, SHA, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OLD = joinpath(ROOT, "numerical_audit", "prospective_window_ordering")
include(joinpath(ROOT, "numerical_audit", "julia_cubic_feasibility", "src", "ResponseCore.jl"))
include(joinpath(OLD, "src", "WindowStudyIO.jl"))
include(joinpath(OLD, "src", "WindowPhysical.jl"))
using .WindowStudyIO, .WindowPhysical

const TIMES = (.65, 1.35)
const CAPS = (8, 24)
const MAX_BYTES = 8 * 1024^2
const ROW_FIELDS = ("name", "order_id", "arm", "repetition", "cpu_budget", "partition",
    "formula", "time", "budget_multiplier", "order", "constructor_sha256", "parent_sha256",
    "claimed_score", "score_kind")
require(ok, message) = ok || throw(ArgumentError(message))
filehash(path) = bytes2hex(sha256(read(path)))
hexhash(value) = value isa String && occursin(r"^[0-9a-f]{64}$", value)
label(value) = value isa String && 1 <= ncodeunits(value) <= 200 &&
    all(c -> ' ' <= c <= '~', value)
integer(value) = value isa Integer && !(value isa Bool)

function bounded_toml(path)
    require(isfile(path) && filesize(path) <= MAX_BYTES, "bounded TOML file required")
    raw = open(path, "r") do io; read(io, MAX_BYTES + 1); end
    require(length(raw) <= MAX_BYTES, "TOML grew beyond admission cap")
    # Bound nesting before parsing; scientific strings are single-line.
    quotechar = UInt8(0); escaped = false; comment = false; stack = UInt8[]
    for (index, byte) in enumerate(raw)
        if comment
            byte in (0x0a, 0x0d) && (comment = false)
        elseif quotechar != 0
            require(!(byte in (0x0a, 0x0d)), "multiline TOML string refused")
            if escaped; escaped = false
            elseif quotechar == UInt8('"') && byte == UInt8('\\'); escaped = true
            elseif byte == quotechar; quotechar = 0
            end
        elseif byte == UInt8('#'); comment = true
        elseif byte in (UInt8('"'), UInt8('\''))
            require(!(index + 2 <= length(raw) && raw[index + 1] == byte && raw[index + 2] == byte),
                "triple-quoted TOML strings are outside scientific schema")
            quotechar = byte
        elseif byte in (UInt8('['), UInt8('{'))
            require(length(stack) < 8, "TOML nesting cap"); push!(stack, byte)
        elseif byte in (UInt8(']'), UInt8('}'))
            opening = byte == UInt8(']') ? UInt8('[') : UInt8('{')
            require(!isempty(stack) && pop!(stack) == opening, "TOML nesting mismatch")
        end
    end
    require(quotechar == 0 && !escaped && isempty(stack), "unfinished TOML token")
    (data=TOML.parse(String(copy(raw))), sha256=bytes2hex(sha256(raw)))
end

function write_new(path, value)
    temporary = path * ".partial"
    require(!ispath(path) && !ispath(temporary), "refuse existing result/partial: " * path)
    open(temporary, "w") do io
        TOML.print(io, WindowStudyIO.clean(value); sorted=true)
    end
    mv(temporary, path; force=false)
end

function save_state(path, state)
    temporary = path * ".partial"
    require(!ispath(path) && !ispath(temporary), "refuse existing state/partial")
    open(temporary, "w") do io
        for value in state
            write(io, htol(reinterpret(UInt64, real(value))))
            write(io, htol(reinterpret(UInt64, imag(value))))
        end
    end
    mv(temporary, path; force=false)
    filehash(path)
end

function source_hashes()
    names = ("numerical_audit/cost_matched_continuation/physical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
        "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl")
    Dict(name => filehash(joinpath(ROOT, name)) for name in names)
end

logical_identity(row) = (row["arm"], row["repetition"], row["cpu_budget"], row["time"], row["budget_multiplier"])
order_identity(row) = Tuple(row[key] for key in ROW_FIELDS if !(key in ("name", "time", "budget_multiplier")))

function admit(input, candidates_path)
    parsed = bounded_toml(input); raw = parsed.data
    require(raw["n"] === 14 && raw["electrons"] === 4 && raw["partition"] == "conserving",
        "fourteen-mode/four-electron original conserving model required")
    require(label(raw["case"]) && hexhash(raw["source_json_sha256"]), "bounded model identity")
    units = raw["units"]
    require(units isa Vector && 1 <= length(units) <= 1024, "original unit cap")
    require(all(u -> u isa AbstractDict && u["words"] isa Vector && u["numerators"] isa Vector &&
        1 <= length(u["words"]) == length(u["numerators"]) <= 8192, units), "unit word arrays")
    require(sum(length(u["words"]) for u in units) <= 8192, "original word cap")
    manifest = bounded_toml(candidates_path); m = manifest.data
    require(length(m) <= 24 && m["version"] === 1, "candidate header version/size")
    for key in ("input_sha256", "source_json_sha256", "protocol_sha256", "source_freeze_sha256")
        require(hexhash(m[key]), "candidate SHA-256 spelling: " * key)
    end
    require(m["input_sha256"] == parsed.sha256 && m["source_json_sha256"] == raw["source_json_sha256"] &&
        m["case"] == raw["case"], "candidate exact input/case identity")
    mode = m["mode"]; require(mode in ("main", "preflight"), "frozen panel mode")
    reps = mode == "main" ? 2 : 1
    require(m["repetitions"] === reps && m["expected_logical_rows"] === 20reps, "candidate roster size")
    expected = Set((arm, rep, cpu, t, cap) for rep in 1:reps
        for (arm, cpu) in (("parent", 0), ("full", 32), ("full", 64), ("window", 32), ("window", 64))
        for t in TIMES for cap in CAPS)
    rows = m["candidates"]
    require(rows isa Vector && length(rows) == length(expected), "complete bounded candidate rows")
    observed = Set{Any}(); names = Set{String}(); orders = Dict{String,Any}(); counts = Dict{String,Int}()
    for row in rows
        require(row isa AbstractDict && length(row) <= 24 && all(key -> haskey(row, key), ROW_FIELDS), "candidate fields")
        require(label(row["name"]) && label(row["order_id"]) && !(row["name"] in names), "unique bounded ASCII candidate name")
        push!(names, row["name"])
        require(all(key -> integer(row[key]), ("repetition", "cpu_budget", "formula", "budget_multiplier")), "integer row metadata")
        require(row["time"] isa AbstractFloat && isfinite(row["time"]) && row["formula"] == 2 &&
            row["partition"] == "conserving", "frozen time/formula/partition")
        key = logical_identity(row); require(key in expected && !(key in observed), "unknown/duplicate logical identity")
        push!(observed, key)
        order = row["order"]
        require(order isa Vector && length(order) == length(units) && all(integer, order) &&
            sort(order) == collect(1:length(units)), "whole original-label order")
        require(hexhash(row["constructor_sha256"]) && hexhash(row["parent_sha256"]), "constructor/parent hashes")
        require(row["claimed_score"] isa AbstractFloat && isfinite(row["claimed_score"]) && row["claimed_score"] >= 0 &&
            row["score_kind"] == "diagonal_reference_cubic_response", "finite scored constructor witness")
        require(haskey(row, "constructor_score") && row["constructor_score"] isa Real &&
            !(row["constructor_score"] isa Bool) && isfinite(row["constructor_score"]) &&
            row["constructor_score"] == row["claimed_score"], "matching constructor score aliases")
        oid = row["order_id"]; signature = order_identity(row)
        require(!haskey(orders, oid) || orders[oid] == signature, "order_id changes constructor identity/order")
        orders[oid] = signature; counts[oid] = get(counts, oid, 0) + 1
    end
    require(observed == expected && length(orders) == 5reps && all(==(4), values(counts)), "complete constructor graph")
    for rep in 1:reps
        parent = only(filter(row -> row["arm"] == "parent" && row["repetition"] == rep &&
            row["time"] == first(TIMES) && row["budget_multiplier"] == first(CAPS), rows))
        require(parent["parent_sha256"] == parent["constructor_sha256"], "parent self identity")
        require(all(row -> row["repetition"] != rep || row["parent_sha256"] == parent["constructor_sha256"], rows),
            "continuations must share the exact frozen parent")
    end
    (manifest=m, candidates_sha256=manifest.sha256, input_sha256=parsed.sha256)
end

function run_case(input, candidates_path, output)
    require(!ispath(output), "fresh primary output directory required"); mkdir(output)
    began = cpu_seconds(); wall = time_ns(); stage = "input_admission"
    records = Dict{String,Any}[]
    report = Dict{String,Any}("schema" => "cost-matched-continuation-physical-v1", "status" => "incomplete",
        "physical_panel_complete" => false, "rows" => records, "source_hashes" => source_hashes(),
        "reference_records" => Any[], "failures" => Any[], "n" => 14, "hf" => "15",
        "new_physical_propagation" => true, "state_format" => "little-endian interleaved real/imag Float64;16 bytes per amplitude",
        "scope" => "Full 16384-dimensional ideal original-word circuits and checked fixed-N numerical reference; raw metrics, no interval or hardware certificate.")
    try
        admitted = admit(input, candidates_path); manifest = admitted.manifest
        for key in ("input_sha256", "source_json_sha256", "protocol_sha256", "source_freeze_sha256",
                "case", "mode", "repetitions", "expected_logical_rows")
            report[key] = manifest[key]
        end
        report["candidates_sha256"] = admitted.candidates_sha256
        data = load_model(input)
        require(data.input_sha256 == admitted.input_sha256 && data.n == 14 && data.electrons == 4,
            "model changed after admission")
        write_new(joinpath(output, "PHYSICAL.START.toml"), report)
        stage = "reference_construction"
        reference = reference_sector(data.units, data.n, data.hf; max_sector_dim=2000)
        report["reference_sector_dimension"] = length(reference.basis)
        report["reference_outside_coupling_max"] = reference.outside_coupling_max
        report["reference_coupling_tolerance"] = reference.coupling_tolerance
        report["reference_hermiticity_residual"] = reference.hermiticity_residual
        references = Dict{Float64,Any}()
        for duration in TIMES
            exact = reference_state(reference, duration)
            filename = "reference_" * bitstring(duration) * ".bin"
            digest = save_state(joinpath(output, filename), exact.state)
            leakage = state_leakage(exact.state, data.n, data.hf)
            references[duration] = (state=exact.state, file=filename, sha256=digest,
                norm_drift=exact.norm_drift, leakage=leakage)
            push!(report["reference_records"], Dict("time" => duration, "reference_file" => filename,
                "reference_sha256" => digest, "norm_drift" => exact.norm_drift,
                "number_leakage" => leakage.number_leakage, "spin_leakage" => leakage.spin_leakage))
        end
        buffers = (zeros(ComplexF64, 1 << data.n), zeros(ComplexF64, 1 << data.n))
        cache = Dict{Any,Any}()
        for (index, row) in enumerate(manifest["candidates"])
            stage = "physical_row_" * string(index)
            order = Int.(row["order"]); duration = Float64(row["time"])
            budget = row["budget_multiplier"] * data.words
            steps = WindowPhysical._step_count(data.units, order, budget, 2)
            key = (Tuple(order), 2, reinterpret(UInt64, duration), steps)
            duplicate = haskey(cache, key); exact = references[duration]
            if !duplicate
                circuit = simulate_order(data.units, data.n, data.hf, order, duration, budget; formula=2, buffers=buffers)
                metric = WindowPhysical._projective_error(exact.state, circuit.state)
                leakage = state_leakage(circuit.state, data.n, data.hf)
                filename = "actual_" * bytes2hex(sha256(string(key))) * ".bin"
                digest = save_state(joinpath(output, filename), circuit.state)
                cache[key] = (file=filename, sha256=digest, steps=circuit.steps, rotations=circuit.ideal_rotations,
                    norm_drift=circuit.norm_drift, projective_error=metric.projective_error,
                    infidelity=metric.infidelity, overlap_abs=metric.overlap_abs,
                    number_leakage=leakage.number_leakage, spin_leakage=leakage.spin_leakage)
            end
            saved = cache[key]
            record = Dict{String,Any}(key => row[key] for key in ROW_FIELDS)
            merge!(record, Dict("rotation_budget" => budget, "steps" => saved.steps, "ideal_rotations" => saved.rotations,
                "projective_error" => saved.projective_error, "infidelity" => saved.infidelity, "overlap_abs" => saved.overlap_abs,
                "norm_drift" => saved.norm_drift, "number_leakage" => saved.number_leakage, "spin_leakage" => saved.spin_leakage,
                "reference_norm_drift" => exact.norm_drift, "reference_number_leakage" => exact.leakage.number_leakage,
                "reference_spin_leakage" => exact.leakage.spin_leakage, "reference_file" => exact.file,
                "reference_sha256" => exact.sha256, "actual_file" => saved.file, "actual_sha256" => saved.sha256,
                "propagation_reused_for_identical_identity" => duplicate))
            push!(records, record)
            report["unique_actual_propagations"] = length(cache)
            write_new(joinpath(output, "PHYSICAL.PARTIAL." * lpad(string(index), 3, '0') * ".toml"), report)
            println("PHYSICAL_ROW ", data.case, " ", index, "/", length(manifest["candidates"])); flush(stdout)
        end
        report["physical_panel_complete"] = true; report["status"] = "complete"
    catch error
        push!(report["failures"], Dict("phase" => stage, "type" => string(typeof(error)), "message" => sprint(showerror, error)))
    end
    report["completed_logical_rows"] = length(records)
    report["CPU_seconds"] = cpu_seconds() - began; report["wall_seconds"] = (time_ns() - wall) / 1e9
    report["cost_scope"] = "Driver CPU only; complete startup/import/teardown cost comes from the owning OS guard receipt"
    write_new(joinpath(output, "PHYSICAL.toml"), report)
    report
end

function main(args=ARGS)
    length(args) == 3 || error("usage: physical.jl INPUT_MODEL.toml CANDIDATES.toml FRESH_OUTPUT_DIRECTORY")
    BLAS.set_num_threads(1)
    require(Threads.nthreads() == 1, "one Julia numerical thread required")
    run_case(args...)["status"] == "complete" ? 0 : 2
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(CostMatchedPhysical.main())
end
