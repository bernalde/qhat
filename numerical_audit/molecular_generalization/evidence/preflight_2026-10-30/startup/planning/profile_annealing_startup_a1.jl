# Old-only engineering diagnostic. This wrapper does not amend a constructor,
# its request, the immutable producer sources, or scientific admission.
using Profile, TOML, SHA

const INNER_CPU_CAP = 22.335336382999994
const INNER_WALL_CAP = 119.27459755403979
const PROFILE_IP_CAP = 200_000
const PROFILE_DELAY_SECONDS = 0.05
const TREE_BYTE_CAP = 96 * 1024
const FLAT_BYTE_CAP = 64 * 1024

cpu_seconds() = ccall(:clock, Int64, ()) / 1e6
require(ok, message) = ok || throw(ArgumentError(message))

struct ProfileOutputLimit <: Exception end
mutable struct LimitedProfileIO{T<:IO} <: IO
    stream::T
    remaining::Int
end
Base.isopen(io::LimitedProfileIO) = isopen(io.stream)
Base.iswritable(io::LimitedProfileIO) = true
Base.flush(io::LimitedProfileIO) = flush(io.stream)
Base.get(io::LimitedProfileIO, key::Symbol, default) = get(io.stream, key, default)
function Base.write(io::LimitedProfileIO, byte::UInt8)
    io.remaining > 0 || throw(ProfileOutputLimit())
    count = write(io.stream, byte)
    io.remaining -= count
    count
end
function Base.unsafe_write(io::LimitedProfileIO, pointer::Ptr{UInt8}, count::UInt)
    allowed = min(count, UInt(io.remaining))
    if allowed > 0
        Base.unsafe_write(io.stream, pointer, allowed)
        io.remaining -= Int(allowed)
    end
    allowed == count || throw(ProfileOutputLimit())
    Int(allowed)
end

function bounded_bytes(path, maximum)
    require(isfile(path), "regular scientific input required")
    observed = filesize(path)
    require(0 <= observed <= maximum, "scientific input byte cap")
    open(path, "r") do stream
        raw = read(stream, observed + 1)
        require(length(raw) == observed && eof(stream), "scientific input size changed")
        raw
    end
end

function write_new_toml(path, value)
    require(!ispath(path) && !ispath(path * ".partial"), "fresh diagnostic output required")
    buffer = IOBuffer()
    TOML.print(buffer, value; sorted=true)
    raw = take!(buffer)
    require(length(raw) <= 32 * 1024, "diagnostic metadata byte cap")
    open(path * ".partial", "w") do stream
        write(stream, raw)
    end
    mv(path * ".partial", path; force=false)
end

function measured_phase(operation, name, phases, output)
    println("PROFILE START ", name); flush(stdout)
    write_new_toml(joinpath(output, name * ".START.toml"),
        Dict("phase" => name, "process_CPU_seconds" => cpu_seconds()))
    began = cpu_seconds(); wall = time_ns()
    Base.cumulative_compile_timing(true)
    compilation = Base.cumulative_compile_time_ns()
    complete = false
    try
        value = operation()
        complete = true
        value
    finally
        compiled = Base.cumulative_compile_time_ns() .- compilation
        Base.cumulative_compile_timing(false)
        phases[name] = Dict("CPU_seconds" => cpu_seconds() - began,
            "wall_seconds" => (time_ns() - wall) / 1e9,
            "Julia_reported_compilation_seconds" => compiled[1] / 1e9,
            "Julia_reported_recompilation_seconds" => compiled[2] / 1e9,
            "returned_normally" => complete)
        write_new_toml(joinpath(output, name * ".END.toml"), phases[name])
        println("PROFILE END ", name); flush(stdout)
    end
end

function write_bounded_profile(path, data; format, maximum)
    truncated = false
    open(path, "w") do stream
        limited = LimitedProfileIO(stream, maximum - 64)
        context = IOContext(limited, :displaysize => (80, 150), :color => false, :limit => true)
        try
            if format == :tree
                Profile.print(context, data; format=:tree, C=true, combine=true,
                    maxdepth=40, mincount=1, noisefloor=0.0, recur=:flatc, groupby=:none)
            else
                Profile.print(context, data; format=:flat, C=true, combine=true,
                    sortedby=:count, mincount=1, groupby=:none)
            end
        catch error
            error isa ProfileOutputLimit || rethrow()
            truncated = true
            write(stream, "\n[bounded profile output truncated]\n")
        end
    end
    Dict("bytes" => filesize(path), "byte_cap" => maximum, "truncated" => truncated)
end

function main(args=ARGS)
    length(args) == 3 || error("usage: profile_annealing_startup_a1.jl IMMUTABLE_ANNEAL.jl ORIGINAL_B32_REQUEST.toml FRESH_PROFILE_DIRECTORY")
    VERSION == v"1.12.6" || error("diagnostic is pinned to Julia 1.12.6")
    Threads.nthreads() == 1 || error("one Julia numerical thread required")
    source, request = realpath(args[1]), realpath(args[2])
    require(source == normpath(abspath(args[1])) && request == normpath(abspath(args[2])), "canonical original scientific paths")
    data_root = normpath(joinpath(dirname(source), "..", "..", ".."))
    require(source == joinpath(data_root, "source_snapshot", "numerical_audit", "molecular_generalization", "anneal.jl"),
        "immutable source_snapshot worker required")
    require(request == joinpath(data_root, "workers", "mg-lih-r1p487-631g-14m4e-r1-annealing-b32.annealing.toml.request.toml"),
        "original old LiH B32 request required")
    output = normpath(abspath(args[3]))
    require(!ispath(output) && !(output == data_root || startswith(output, data_root * "/")),
        "fresh external diagnostic directory required")
    mkpath(dirname(output)); mkdir(output)
    checkpoint = joinpath(output, "DIAGNOSTIC_CHECKPOINT.toml")
    summary = Dict{String,Any}("schema" => "molecular-old-annealing-startup-profile-v1",
        "scope" => "OLD engineering diagnostic only; no constructor retry/admission or new-family error",
        "source_path" => source, "request_path" => request,
        "source_sha256" => bytes2hex(sha256(bounded_bytes(source, 1024^2))),
        "request_sha256" => bytes2hex(sha256(bounded_bytes(request, 8 * 1024^2))),
        "wrapper_sha256" => bytes2hex(sha256(bounded_bytes(@__FILE__, 1024^2))),
        "inner_absolute_process_CPU_cap" => INNER_CPU_CAP, "inner_wall_cap" => INNER_WALL_CAP,
        "profile_IP_cap" => PROFILE_IP_CAP, "profile_delay_seconds" => PROFILE_DELAY_SECONDS,
        "julia_version" => string(VERSION), "phases" => Dict{String,Any}(),
        "worker_return_code" => 2,
        "accounting_note" => "Wrapper/Profile startup consumes the unchanged absolute inner CPU cap. External 45CPU/225wall receipt is authoritative for the entire diagnostic, including report emission.",
        "phase_note" => "include/run times are observed here. Request/context times are copied only from an actual complete worker checkpoint; absence remains unavailable. Compilation fields are Julia's reported elapsed metrics, not subtracted from process CPU.")
    write_new_toml(joinpath(output, "PROFILE.START.toml"), summary)
    summary["profile_initialized"] = false
    began = cpu_seconds(); wall = time_ns()
    try
        Profile.init(n=PROFILE_IP_CAP, delay=PROFILE_DELAY_SECONDS)
        Profile.clear()
        summary["profile_initialized"] = true
        Profile.@profile begin
            measured_phase("include_original", summary["phases"], output) do
                Base.include(Main, source)
            end
            measured_phase("run_original", summary["phases"], output) do
                worker = getfield(Main, :MolecularAnnealingWorker)
                summary["worker_return_code"] = Base.invokelatest(worker.run, request, checkpoint,
                    INNER_CPU_CAP, INNER_WALL_CAP)
            end
        end
    catch error
        summary["failure_type"] = string(typeof(error))
        summary["failure_message"] = first(sprint(showerror, error), 1000)
    finally
        summary["CPU_before_profile_emission"] = cpu_seconds() - began
        summary["wall_before_profile_emission"] = (time_ns() - wall) / 1e9
        try
            require(summary["profile_initialized"], "profile sampling was not initialized")
            samples = Profile.fetch(limitwarn=false)
            summary["profile_IPs_retained"] = length(samples)
            summary["profile_buffer_may_be_full"] = length(samples) >= PROFILE_IP_CAP
            summary["tree"] = write_bounded_profile(joinpath(output, "PROFILE.tree.txt"), samples;
                format=:tree, maximum=TREE_BYTE_CAP)
            summary["flat"] = write_bounded_profile(joinpath(output, "PROFILE.flat.txt"), samples;
                format=:flat, maximum=FLAT_BYTE_CAP)
        catch error
            summary["profile_emission_failure"] = first(sprint(showerror, error), 500)
        end
        summary["checkpoint_present"] = isfile(checkpoint)
        if isfile(checkpoint)
            try
                raw = bounded_bytes(checkpoint, 64 * 1024)
                observed = TOML.parse(String(copy(raw)))
                summary["checkpoint_sha256"] = bytes2hex(sha256(raw))
                for key in ("phases", "status", "start_scores_completed", "evaluations", "proposals", "chain_proposals")
                    haskey(observed, key) && (summary["worker_" * key] = observed[key])
                end
            catch error
                summary["checkpoint_summary_failure"] = first(sprint(showerror, error), 500)
            end
        end
        summary["CPU_before_final_summary"] = cpu_seconds() - began
        summary["wall_before_final_summary"] = (time_ns() - wall) / 1e9
        write_new_toml(joinpath(output, "PROFILE.toml"), summary)
        println("PROFILE diagnostic written; worker return code ", summary["worker_return_code"])
        flush(stdout)
    end
    summary["worker_return_code"]
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
