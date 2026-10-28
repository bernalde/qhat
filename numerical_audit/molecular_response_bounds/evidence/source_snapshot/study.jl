#!/usr/bin/env julia
using TOML, SHA
include(joinpath(@__DIR__,"src","MolecularResponseBounds.jl"))
using .MolecularResponseBounds
using .MolecularResponseBounds.BoundsCommon
using .MolecularResponseBounds.ExactResponse
using .MolecularResponseBounds.StudyAPI

cpu()=ccall(:clock,Int64,())/1.0e6
function sources()
    paths=["study.jl","independent_check.jl"]
    append!(paths,["src/"*name for name in sort(readdir(joinpath(@__DIR__,"src"))) if endswith(name,".jl")])
    Dict(p=>filehash(joinpath(@__DIR__,p)) for p in paths)
end
function run_registry(registry_path,output_dir;caps_profile=nothing)
    filesize(registry_path)<=Caps().max_input_bytes || error("registry byte cap")
    registry=TOML.parsefile(registry_path)
    registry["schema"]=="molecular-bounds-registry-v1" || error("registry schema differs")
    caps,caps_declaration=select_caps(get(registry,"caps",Dict{String,Any}());cli_profile=caps_profile)
    caps_snapshot=caps_dict(caps);profile=caps_declaration["effective_profile"]
    cases=registry["cases"];1<=length(cases)<=32 || error("registry case count outside1:32")
    ids=String[c["id"] for c in cases]
    length(unique(ids))==length(ids) && all(s->occursin(r"^[A-Za-z0-9_.-]+$",s),ids) || error("invalid/duplicate case ID")
    options=get(registry,"options",Dict{String,Any}())
    allowed=Set(["iterations","max_witnesses","target_count","comparison_count","attempt_CPU_seconds","attempt_wall_seconds"])
    Set(keys(options))⊆allowed || error("unknown registry option")
    iterations=Int(get(options,"iterations",64));witnesses=Int(get(options,"max_witnesses",8))
    targetcount=Int(get(options,"target_count",8));comparisons=Int(get(options,"comparison_count",1))
    0<=iterations<=64 && 0<=witnesses<=8 && 1<=targetcount<=8 && 0<=comparisons<=3 || error("algorithm option outside frozen ceilings")
    cpu_cap=Float64(get(options,"attempt_CPU_seconds",180.0));wall_cap=Float64(get(options,"attempt_wall_seconds",900.0))
    0<cpu_cap<=240 && 0<wall_cap<=1200 || error("attempt cap outside ceiling")
    mkpath(output_dir);root=dirname(abspath(registry_path));hashes=sources()
    summary=Dict{String,Any}("schema"=>"molecular-bounds-run-v1","registry_sha256"=>filehash(registry_path),
        "status"=>"started","cases"=>Dict{String,Any}[],"source_hashes"=>hashes,"julia_version"=>string(VERSION),
        "caps"=>caps_snapshot,"caps_declaration"=>caps_declaration,"caps_profile"=>profile,
        "caps_profile_validation"=>"validated_bounded_profile","options"=>options)
    write_record(joinpath(output_dir,"RUN.START.toml"),summary;caps)
    for case in cases
        id=String(case["id"]);started=cpu();wall=time_ns()
        guard=()->begin
            cpu()-started<=cpu_cap || throw(Refusal("cpu_limit","case CPU allowance exhausted"))
            (time_ns()-wall)/1e9<=wall_cap || throw(Refusal("wall_limit","case wall allowance exhausted"))
        end
        input=abspath(joinpath(root,String(case["input"])))
        row=Dict{String,Any}("id"=>id,"complete"=>false,"status"=>"started","source_hashes"=>hashes,
            "input_file"=>String(case["input"]),
            "registry_sha256"=>summary["registry_sha256"],"caps"=>caps_snapshot,
            "caps_declaration"=>caps_declaration,"caps_profile"=>profile,
            "caps_profile_validation"=>"validated_bounded_profile","input_admission_status"=>"not_checked")
        write_record(joinpath(output_dir,id*".START.toml"),row;caps)
        try
            filesize(input)<=caps.max_input_bytes && (row["input_sha256"]=filehash(input))
            data=read_coefficients(input;caps,guard)
            row["input_file"]=String(case["input"]);row["input_sha256"]=data.input_sha256
            row["input_admission_status"]="admitted_under_declared_caps"
            row["input_n"]=data.n;row["input_units"]=length(data.units)
            row["input_original_words"]=sum(length(u.numerators) for u in data.units)
            kernelpath=if get(case,"freeze_kernel",false)
                out=joinpath(abspath(output_dir),id*".KERNEL.toml")
                freeze_kernel(input,out;time=get(case,"time",1.0),caps,guard);out
            else
                abspath(joinpath(root,String(case["kernel"])))
            end
            row["kernel_file"]=kernelpath;row["kernel_sha256"]=filehash(kernelpath)
            field=String(get(case,"order_field","cubic_order"))
            order=haskey(case,"order") ? Int.(case["order"]) : Int.(data.source[field])
            lengths=[length(u.numerators) for u in data.units]
            budget=haskey(case,"ideal_budget") ? bigparse(case["ideal_budget"],caps) : bigparse(get(case,"budget_multiplier",30),caps)*sum(lengths)
            record=molecular_record(input,kernelpath,order,budget;target_count=targetcount,
                comparison_count=comparisons,iterations,max_witnesses=witnesses,caps,guard,
                progress=(stage,value)->write_record(joinpath(output_dir,id*".PHASE."*stage*".toml"),
                    Dict("id"=>id,"phase"=>stage,"phase_complete"=>true,"full_certificate_complete"=>false,
                        "data"=>value,"source_hashes"=>hashes,"registry_sha256"=>summary["registry_sha256"],
                        "caps"=>caps_snapshot,"caps_declaration"=>caps_declaration,"caps_profile"=>profile,
                        "caps_profile_validation"=>"validated_bounded_profile");caps))
            merge!(row,record);row["input_file"]=String(case["input"]);row["kernel_file"]=kernelpath
        catch e
            row["complete"]=false;row["status"]=e isa Refusal ? e.status : "error"
            row["reason"]=sprint(showerror,e);row["positive_full_bound"]=false
            if row["input_admission_status"]=="not_checked"
                row["input_admission_status"]=e isa Refusal ? "refused_by_declared_caps" : "invalid_input"
            end
            # Errors/refusals remain distinct from completed zero certificates.
        end
        row["CPU_seconds"]=cpu()-started;row["wall_seconds"]=(time_ns()-wall)/1e9
        path=joinpath(output_dir,id*".RESULT.toml");write_record(path,row;caps)
        push!(summary["cases"],Dict("id"=>id,"complete"=>row["complete"],"status"=>row["status"],
            "result"=>basename(path),"sha256"=>filehash(path)))
    end
    summary["status"]="inventory_complete"
    summary["all_cases_completed"]=all(c->c["complete"],summary["cases"])
    write_record(joinpath(output_dir,"RUN.RESULT.toml"),summary;caps)
    summary
end
function study_arguments(args)
    profile,positional=parse_caps_cli(args)
    length(positional)==2 || error("usage: study.jl [--caps-profile default|larger14] REGISTRY.toml NEW_OUTPUT_DIRECTORY")
    (caps_profile=profile,registry=positional[1],output=positional[2])
end
function main(args=ARGS)
    options=study_arguments(args)
    result=run_registry(options.registry,options.output;caps_profile=options.caps_profile)
    result["all_cases_completed"] ? 0 : 2
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(main())
end
