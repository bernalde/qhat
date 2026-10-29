# Full-space physical panel after the complete constructor/order freeze.
# Only the reference uses the checked fixed-number sector. Every distinct actual
# circuit is propagated in 2^n space and checkpointed before the next circuit.
using TOML, SHA
include(joinpath(@__DIR__,"..","julia_cubic_feasibility","src","ResponseCore.jl"))
include(joinpath(@__DIR__,"src","WindowStudyIO.jl"))
include(joinpath(@__DIR__,"src","WindowPhysical.jl"))
using .ResponseCore, .WindowStudyIO, .WindowPhysical

const ARMS=("canonical","magnitude","current","first","pauli","cubic","window_sparse","window_magnitude")
const TIMES=(0.4,1.6)
const CAPS=(8,24)
forms(arm)=arm in ("canonical","current","first") ? (1,2) : (2,)
filehash(path)=bytes2hex(sha256(read(path)))

function source_hashes()
    names=("physical.jl","src/WindowPhysical.jl","src/WindowStudyIO.jl",
           "../julia_cubic_feasibility/src/ResponseCore.jl")
    Dict(name=>filehash(joinpath(@__DIR__,name)) for name in names)
end

function save_state(path,state)
    ispath(path) && error("refuse existing state")
    temporary=path*".partial"
    ispath(temporary) && error("refuse retained partial state")
    open(temporary,"w") do io
        for value in state
            write(io,htol(reinterpret(UInt64,real(value))))
            write(io,htol(reinterpret(UInt64,imag(value))))
        end
    end
    mv(temporary,path)
    filehash(path)
end

function structured_word(word)
    isempty(word) && return ()
    Tuple((parse(Int,part[1:end-1]),last(part)) for part in split(word,','))
end

function pauli_partition(data)
    entries=[]
    for (unit,raw) in zip(data.units,data.source["units"]), (index,word) in enumerate(raw["words"])
        push!(entries,(key=structured_word(word),unit=Unit(unit.flip,
            [unit.parity[index]],[unit.signs[index]],[unit.coeff[index]])))
    end
    length(Set(entry.key for entry in entries))==data.words || error("global Pauli words must be unique")
    sort!(entries;by=entry->entry.key)
    Unit[entry.unit for entry in entries]
end

function grid(mode)
    mode in ("main","preflight") || error("unknown frozen panel mode")
    repetitions=mode=="main" ? 2 : 1
    Set((arm,rep,formula,time,cap) for arm in ARMS for rep in 1:repetitions
        for formula in forms(arm) for time in TIMES for cap in CAPS)
end

function row_identity(row)
    arm=row["arm"];rep=row["repetition"];formula=row["formula"]
    time=row["time"];cap=row["budget_multiplier"]
    arm isa String && arm in ARMS || error("unknown arm")
    all(value->typeof(value)===Int64,(rep,formula,cap)) || error("integer row identity")
    time isa AbstractFloat && isfinite(time) || error("finite physical time")
    (arm,rep,formula,Float64(time),cap)
end

function run_case(input,orders_path,output)
    ispath(output) && error("fresh physical output directory required")
    filesize(orders_path)<=8*1024^2 || error("candidate manifest byte cap")
    candidates=TOML.parsefile(orders_path)
    data=load_model(input)
    data.n==14 && data.electrons==4 && data.partition=="conserving" || error("declared fourteen-mode/four-electron conserving input")
    candidates["input_sha256"]==data.input_sha256 || error("order/TOML input identity differs")
    mode=String(candidates["mode"]);expected=grid(mode)
    rows=candidates["candidates"]
    rows isa Vector && length(rows)==length(expected) || error("complete candidate matrix required")
    keys=row_identity.(rows)
    length(Set(keys))==length(expected) && Set(keys)==expected || error("duplicate or missing candidate identity")
    names=[row["name"] for row in rows]
    all(name->name isa String && 1<=sizeof(name)<=256,names) && length(Set(names))==length(names) || error("unique bounded candidate names")
    paulies=pauli_partition(data)
    # Admit every full-order witness before allocating a sector matrix/state.
    for row in rows
        partition=row["partition"]
        partition==(row["arm"]=="pauli" ? "pauli" : "conserving") || error("arm partition differs")
        units=partition=="pauli" ? paulies : data.units
        order=row["order"]
        order isa Vector && all(value->typeof(value)===Int64,order) && sort(order)==collect(1:length(units)) || error("invalid original-label candidate")
    end
    mkdir(output)
    started=cpu_seconds();wall=time_ns();records=Dict{String,Any}[]
    header=Dict{String,Any}("status"=>"incomplete","mode"=>mode,"case"=>data.case,
        "source_hashes"=>source_hashes(),
        "input_sha256"=>data.input_sha256,"source_json_sha256"=>data.source["source_json_sha256"],
        "orders_sha256"=>filehash(orders_path),"n"=>data.n,"hf"=>string(data.hf),
        "expected_logical_rows"=>length(expected),"physical_panel_complete"=>false,
        "rows"=>records,"new_physical_propagation"=>true,
        "state_format"=>"little-endian interleaved real/imag Float64;16 bytes per amplitude",
        "scope"=>"Full-space ideal original-word circuits; checked fixed-N numerical reference only. No interval, compiled-resource or hardware claim.")
    stage="reference construction"
    try
        reference=reference_sector(data.units,data.n,data.hf;max_sector_dim=2000)
        header["reference_sector_dimension"]=length(reference.basis)
        header["reference_outside_coupling_max"]=reference.outside_coupling_max
        header["reference_coupling_tolerance"]=reference.coupling_tolerance
        header["reference_hermiticity_residual"]=reference.hermiticity_residual
        references=Dict{Float64,Any}()
        for duration in TIMES
            exact=reference_state(reference,duration)
            filename="reference_"*string(reinterpret(UInt64,duration))*".bin"
            digest=save_state(joinpath(output,filename),exact.state)
            references[duration]=(state=exact.state,file=filename,sha256=digest,norm_drift=exact.norm_drift)
        end
        buffers=(zeros(ComplexF64,1<<data.n),zeros(ComplexF64,1<<data.n))
        cache=Dict{Any,Any}()
        for (index,row) in enumerate(rows)
            stage="physical row "*string(index)
            partition=row["partition"];units=partition=="pauli" ? paulies : data.units
            order=Int.(row["order"]);duration=Float64(row["time"])
            multiplier=Int(row["budget_multiplier"]);budget=multiplier*data.words
            formula=Int(row["formula"])
            steps=WindowPhysical._step_count(units,order,budget,formula)
            key=(partition,Tuple(order),formula,reinterpret(UInt64,duration),steps)
            duplicate=haskey(cache,key)
            exact=references[duration]
            if !duplicate
                circuit=simulate_order(units,data.n,data.hf,order,duration,budget;formula=formula,buffers=buffers)
                metric=WindowPhysical._projective_error(exact.state,circuit.state)
                leakage=state_leakage(circuit.state,data.n,data.hf)
                filename="actual_"*bytes2hex(sha256(string(key)))*".bin"
                digest=save_state(joinpath(output,filename),circuit.state)
                # Keep immutable scalars/paths, never a reference to reused buffers.
                cache[key]=(file=filename,sha256=digest,steps=circuit.steps,
                    rotations=circuit.ideal_rotations,norm_drift=circuit.norm_drift,
                    projective_error=metric.projective_error,infidelity=metric.infidelity,
                    number_leakage=leakage.number_leakage,spin_leakage=leakage.spin_leakage)
            end
            saved=cache[key]
            record=Dict{String,Any}("name"=>row["name"],"arm"=>row["arm"],"repetition"=>row["repetition"],
                "partition"=>partition,"order"=>order,"formula"=>formula,"time"=>duration,
                "budget_multiplier"=>multiplier,"rotation_budget"=>budget,"steps"=>saved.steps,
                "ideal_rotations"=>saved.rotations,"projective_error"=>saved.projective_error,
                "infidelity"=>saved.infidelity,"norm_drift"=>saved.norm_drift,
                "number_leakage"=>saved.number_leakage,"spin_leakage"=>saved.spin_leakage,
                "reference_norm_drift"=>exact.norm_drift,"reference_file"=>exact.file,
                "reference_sha256"=>exact.sha256,"actual_file"=>saved.file,"actual_sha256"=>saved.sha256,
                "propagation_reused_for_identical_identity"=>duplicate)
            push!(records,record)
            write_result(joinpath(output,"ROW_"*lpad(string(index),3,'0')*".toml"),record)
            println("PHYSICAL_ROW ",data.case," ",index,"/",length(rows));flush(stdout)
        end
        length(records)==length(expected) || error("incomplete physical matrix")
        header["status"]="complete";header["physical_panel_complete"]=true
        header["unique_actual_propagations"]=length(cache)
    catch exception
        header["reason"]=sprint(showerror,exception);header["failed_stage"]=stage
    end
    header["CPU_seconds"]=cpu_seconds()-started
    header["wall_seconds"]=(time_ns()-wall)/1e9
    write_result(joinpath(output,"PHYSICAL.toml"),header)
    header
end

function main(args)
    length(args)==3 || error("usage: physical.jl INPUT_TOML FROZEN_ORDERS_TOML NEW_OUTPUT_DIRECTORY")
    result=run_case(args...)
    result["status"]=="complete" || exit(2)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
