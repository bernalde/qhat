module GeneralizationPhysical
using TOML, SHA, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OLD = joinpath(ROOT, "numerical_audit", "prospective_window_ordering")
include(joinpath(ROOT, "numerical_audit", "julia_cubic_feasibility", "src", "ResponseCore.jl"))
include(joinpath(OLD, "src", "WindowStudyIO.jl"))
include(joinpath(OLD, "src", "WindowPhysical.jl"))
using .WindowStudyIO, .WindowPhysical
export admit, admit_manifest, expected_grid, run, source_hashes, write_new, save_state
const TIMES = (.65, 1.35)
const CAPS = (8, 24)
const CONSTRUCTORS = (("parent",0),("canonical",0),("magnitude",0),("full",32),("full",64),
    ("window",32),("window",64),("annealing",32),("annealing",64))
const MAX_BYTES = 8 * 1024^2
const SCORE_KIND = "diagonal_reference_cubic_response"
const FIELDS = ("name", "order_id", "arm", "repetition", "cpu_budget", "partition", "formula",
    "time", "budget_multiplier", "order", "constructor_sha256", "claimed_score", "score_kind", "eligible")
const HEADER = ("version","order_index_base","case","n","electrons","partition","input_sha256",
    "source_json_sha256","protocol_sha256","source_freeze_sha256","mode","repetitions",
    "expected_logical_rows","training_time","training_budget_multiplier","times","rotation_cap_multipliers")
require(ok, message) = ok || throw(ArgumentError(message))
filehash(path) = bytes2hex(sha256(read(path)))
hexhash(value) = value isa String && occursin(r"^[0-9a-f]{64}$",value)
label(value) = value isa String && 1 <= ncodeunits(value) <= 200 && all(c -> '!' <= c <= '~',value)
integer(value) = value isa Integer && !(value isa Bool)
number(value) = value isa Real && !(value isa Bool) && isfinite(value) && value >= 0
failure(e) = Dict("type"=>string(typeof(e)),"message"=>sprint(showerror,e))

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
function bounded_toml(path)
    raw=bounded_bytes(path)
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
    names = ("numerical_audit/molecular_generalization/physical.jl",
        "numerical_audit/molecular_generalization/src/GeneralizationPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl",
        "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
        "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl")
    Dict(name => filehash(joinpath(ROOT,name)) for name in names)
end

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

function admit(input, candidates_path)
    parsed=bounded_toml(input); raw=parsed.data
    require(integer(raw["n"]) && integer(raw["electrons"]) && 1 <= raw["n"] <= 14 &&
        0 <= raw["electrons"] <= raw["n"] && binomial(raw["n"],raw["electrons"]) <= 2000,
        "model register/electron/fixed-N cap before allocation")
    require(raw["partition"] == "conserving" && label(raw["case"]) && hexhash(raw["source_json_sha256"]), "original conserving model identity")
    blocks=raw["units"]
    require(blocks isa Vector && 1 <= length(blocks) <= 1024 && all(u -> u isa AbstractDict &&
        u["words"] isa Vector && u["numerators"] isa Vector && 1 <= length(u["words"]) == length(u["numerators"]) <= 8192,blocks), "original unit arrays")
    require(sum(length(u["words"]) for u in blocks) <= 8192, "original word cap")
    candidate=bounded_toml(candidates_path); m=admit_manifest(candidate.data,length(blocks))
    require(m["input_sha256"] == parsed.sha256 && m["source_json_sha256"] == raw["source_json_sha256"] &&
        m["case"] == raw["case"] && m["n"] == raw["n"] && m["electrons"] == raw["electrons"], "exact candidate/model identity")
    require((raw["n"],raw["electrons"]) == case_specification(m["mode"],raw["case"]), "frozen case specification")
    bindings=bind_artifacts(input,m)
    require(realpath(candidates_path) == joinpath(bindings.data_root,"panels",m["case"]*".toml") &&
        normpath(abspath(candidates_path)) == realpath(candidates_path),"canonical DATA/panels/case.toml required")
    (manifest=m,candidates_sha256=candidate.sha256,input_sha256=parsed.sha256,bindings=bindings)
end

function phase_guard(cpu_limit,wall_limit,callback)
    require(number(cpu_limit) && cpu_limit > 0 && number(wall_limit) && wall_limit > 0, "finite positive phase allowance")
    began=cpu_seconds(); wall=time_ns()
    () -> begin
        callback()
        require(cpu_seconds()-began <= cpu_limit,"physical CPU allowance exhausted")
        require((time_ns()-wall)/1e9 <= wall_limit,"physical wall allowance exhausted")
        nothing
    end
end

function run(input,candidates_path,output;cpu_limit=600.0,wall_limit=3000.0,check_budget=()->nothing)
    require(!ispath(output),"fresh physical output directory required"); mkdir(output)
    began=cpu_seconds(); wall=time_ns(); stage="input_admission"; guard=phase_guard(cpu_limit,wall_limit,check_budget)
    rows=Dict{String,Any}[]
    report=Dict{String,Any}("schema"=>"molecular-generalization-physical-v1","status"=>"incomplete",
        "physical_panel_complete"=>false,"rows"=>rows,"reference_records"=>Any[],"failures"=>Any[],
        "source_hashes"=>source_hashes(),"consumer_source_hashes"=>source_hashes(),"new_physical_propagation"=>true,
        "ratio_floor"=>1e-12,"state_format"=>"little-endian interleaved real/imag Float64;16 bytes per amplitude",
        "reference_scope"=>"checked complete fixed-N sector; no fixed-spin projection",
        "circuit_scope"=>"full 2^n amplitudes, original words and exact unit boundary fusion; no projection")
    try
        guard(); admitted=admit(input,candidates_path); m=admitted.manifest
        for key in HEADER; report[key]=m[key]; end
        # Preserve the complete logical inventory even if a later kernel refuses.
        for row in m["candidates"]
            record=Dict{String,Any}(k=>row[k] for k in FIELDS)
            for k in ("parent_sha256","constructor_score","constructor_failure")
                haskey(row,k) && (record[k]=row[k])
            end
            record["status"]="incomplete"; push!(rows,record)
        end
        append!(report["reference_records"],(Dict{String,Any}("time"=>duration,"status"=>"incomplete") for duration in TIMES))
        report["candidates_sha256"]=admitted.candidates_sha256
        report["artifact_binding_status"]=admitted.bindings.status
        report["artifact_binding_scope"]=admitted.bindings.scope
        report["bound_artifact_hashes"]=admitted.bindings.hashes
        data=load_model(input)
        require(data.input_sha256 == admitted.input_sha256 && data.n == m["n"] && data.electrons == m["electrons"], "model changed after admission")
        report["hf"]=string(data.hf); report["full_dimension"]=1<<data.n
        report["original_words"]=data.words; report["original_units"]=length(data.units)
        report["all_constructor_eligible"]=all(row -> row["eligible"],m["candidates"])
        write_new(joinpath(output,"PHYSICAL.START.toml"),report)
        stage="fixed_N_reference_construction"; guard()
        reference=reference_sector(data.units,data.n,data.hf;max_sector_dim=2000)
        require(length(reference.basis) == binomial(data.n,data.electrons), "complete fixed-N reference dimension")
        guard()
        report["reference_sector_dimension"]=length(reference.basis)
        report["reference_outside_coupling_max"]=reference.outside_coupling_max
        report["reference_coupling_tolerance"]=reference.coupling_tolerance
        report["reference_hermiticity_residual"]=reference.hermiticity_residual
        references=Dict{Float64,Any}()
        for (reference_index,duration) in enumerate(TIMES)
            stage="reference_time_"*string(duration)
            guard(); exact=reference_state(reference,duration); guard()
            leakage=state_leakage(exact.state,data.n,data.hf)
            require(leakage.number_leakage <= 1e-12 && leakage.spin_leakage <= 1e-10,"reference number/spin leakage failure")
            filename="reference_"*bitstring(duration)*".bin"; digest=save_state(joinpath(output,filename),exact.state)
            references[duration]=(state=exact.state,file=filename,sha256=digest,norm_drift=exact.norm_drift,leakage=leakage)
            merge!(report["reference_records"][reference_index],Dict("reference_file"=>filename,"reference_sha256"=>digest,
                "norm_drift"=>exact.norm_drift,"number_leakage"=>leakage.number_leakage,"spin_leakage"=>leakage.spin_leakage,
                "sector_dimension"=>length(reference.basis),"status"=>"complete"))
        end
        buffers=(zeros(ComplexF64,1<<data.n),zeros(ComplexF64,1<<data.n)); cache=Dict{Any,Any}()
        for (index,row) in enumerate(m["candidates"])
            stage="physical_row_"*string(index); guard()
            record=rows[index]
            order=Int.(row["order"]); duration=Float64(row["time"]); budget=row["budget_multiplier"]*data.words
            steps=WindowPhysical._step_count(data.units,order,budget,2)
            key=(Tuple(order),2,reinterpret(UInt64,duration),steps); duplicate=haskey(cache,key); exact=references[duration]
            if !duplicate
                circuit=simulate_order(data.units,data.n,data.hf,order,duration,budget;formula=2,buffers=buffers); guard()
                metric=WindowPhysical._projective_error(exact.state,circuit.state)
                leakage=state_leakage(circuit.state,data.n,data.hf)
                require(leakage.number_leakage <= 1e-10 && leakage.spin_leakage <= 1e-10,"conserving circuit number/spin leakage failure")
                overlap=abs(dot(exact.state,circuit.state))/sqrt(sum(abs2,exact.state)*sum(abs2,circuit.state))
                circuit_hash=bytes2hex(sha256(data.input_sha256*"|"*string(key)))
                filename="actual_"*circuit_hash*".bin"; digest=save_state(joinpath(output,filename),circuit.state)
                cache[key]=(file=filename,sha256=digest,circuit_sha256=circuit_hash,steps=circuit.steps,rotations=circuit.ideal_rotations,
                    norm_drift=circuit.norm_drift,projective_error=metric.projective_error,infidelity=metric.infidelity,
                    overlap_abs=metric.overlap_abs,raw_normalized_overlap_deficit=1-overlap,
                    number_leakage=leakage.number_leakage,spin_leakage=leakage.spin_leakage)
            end
            saved=cache[key]
            merge!(record,Dict("rotation_budget"=>budget,"steps"=>saved.steps,"ideal_rotations"=>saved.rotations,
                "projective_error"=>saved.projective_error,"infidelity"=>saved.infidelity,"overlap_abs"=>saved.overlap_abs,
                "raw_normalized_overlap_deficit"=>saved.raw_normalized_overlap_deficit,"numerical_floor"=>saved.projective_error<=1e-12,
                "norm_drift"=>saved.norm_drift,"number_leakage"=>saved.number_leakage,"spin_leakage"=>saved.spin_leakage,
                "reference_norm_drift"=>exact.norm_drift,"reference_number_leakage"=>exact.leakage.number_leakage,
                "reference_spin_leakage"=>exact.leakage.spin_leakage,"reference_file"=>exact.file,"reference_sha256"=>exact.sha256,
                "actual_file"=>saved.file,"actual_sha256"=>saved.sha256,"circuit_sha256"=>saved.circuit_sha256,
                "propagation_reused_for_identical_identity"=>duplicate,"status"=>"complete"))
            report["unique_actual_propagations"]=length(cache)
            write_new(joinpath(output,"PHYSICAL.PARTIAL."*lpad(string(index),3,'0')*".toml"),report)
        end
        guard(); require(count(row -> row["status"] == "complete",rows) == m["expected_logical_rows"],"incomplete logical physical grid")
        report["physical_panel_complete"]=true; report["status"]="complete"
    catch e
        push!(report["failures"],merge(Dict("phase"=>stage),failure(e)))
    end
    report["completed_logical_rows"]=count(row -> row["status"] == "complete",rows)
    report["logical_rows_retained"]=length(rows)
    report["completed_references"]=count(row -> row["status"] == "complete",report["reference_records"])
    report["CPU_seconds"]=cpu_seconds()-began; report["wall_seconds"]=(time_ns()-wall)/1e9
    report["cost_scope"]="Driver CPU only; owning OS receipt includes startup/import/JIT/teardown"
    report["cooperation_scope"]="Callbacks before/after unchanged primary reference/circuit kernels; OS guard bounds a kernel in progress"
    write_new(joinpath(output,"PHYSICAL.toml"),report); report
end
end
