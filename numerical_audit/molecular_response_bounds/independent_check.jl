#!/usr/bin/env julia
using TOML, SHA
include(joinpath(@__DIR__,"src","MolecularResponseBounds.jl"))
using .MolecularResponseBounds.BoundsCommon
using .MolecularResponseBounds.IndependentCheck

# Admission deliberately precedes every user-artifact hash, even for refusals.
# A missing artifact can support only an explicitly unbound metadata receipt.
function artifact_admission(path,caps)
    (path=="-" || !ispath(path)) && return "unavailable"
    isfile(path) || throw(ArgumentError("artifact must be a regular file"))
    filesize(path)<=caps.max_input_bytes || throw(Refusal("input_limit","input/kernel byte cap exceeded before hashing"))
    "available"
end
function artifact_binding!(out,record,key,path,admission)
    digest=get(record,key*"_sha256",nothing)
    available=digest!==nothing && digest!="unavailable"
    available && !(digest isa AbstractString && occursin(r"^[0-9a-fA-F]{64}$",digest)) &&
        throw(ArgumentError("invalid producer "*key*" digest"))
    out["producer_"*key*"_sha256"]=available ? String(digest) : "unavailable"
    actual=admission=="available" ? filehash(path) : "unavailable"
    out[key*"_sha256"]=actual
    status=if admission!="available"
        available ? "artifact_unavailable" : "producer_digest_and_artifact_unavailable"
    elseif !available
        "producer_digest_unavailable"
    elseif lowercase(String(digest))==actual
        "matched"
    else
        "mismatch"
    end
    out[key*"_binding_status"]=status
    if status=="mismatch"
        out["artifact_binding_status"]="mismatch"
        error("producer "*key*" digest differs from supplied bytes")
    end
    status=="matched"
end
function check_one(input,kernel,result,output;caps_profile=nothing)
    caps=Caps()
    began=ccall(:clock,Int64,())/1.0e6;wall=time_ns()
    guard=()->begin
        ccall(:clock,Int64,())/1.0e6-began<=100 || throw(Refusal("cpu_limit","independent replay100CPU cap"))
        (time_ns()-wall)/1e9<=500 || throw(Refusal("wall_limit","independent replay500wall cap"))
    end
    out=Dict{String,Any}("schema"=>"independent-molecular-bound-v1","complete"=>false,
        "result_sha256"=>"unavailable","input_sha256"=>"unavailable","kernel_sha256"=>"unavailable",
        "caps_validation_status"=>"not_checked","artifact_binding_status"=>"not_checked",
        "artifact_links_authenticated"=>false,"positive_certificate_verified"=>false,
        "scientific_verification_performed"=>false)
    caps_attempted=false
    try
        # Both admissions finish before either user artifact can be hashed/read.
        input_admission=artifact_admission(input,caps)
        kernel_admission=artifact_admission(kernel,caps)
        out["input_admission_status"]=input_admission;out["kernel_admission_status"]=kernel_admission
        filesize(result)<=caps.max_output_bytes || throw(Refusal("input_limit","producer result byte cap exceeded before parsing/hashing"))
        record=TOML.parsefile(result)
        out["result_sha256"]=filehash(result);out["checker_sha256"]=filehash(@__FILE__)
        # Only whitelisted bounded requests construct replay Caps. Every frozen
        # memory/work/bit/register value must equal the hard-coded profile.
        caps_attempted=true
        caps,declaration=verify_caps_declaration(record["caps_declaration"];expected_profile=caps_profile)
        snapshot=record["caps"]
        snapshot isa AbstractDict && Set(keys(snapshot))==Set(keys(caps_dict(caps))) &&
            all(v->v isa Integer && !(v isa Bool),values(snapshot)) || error("effective Caps snapshot schema differs")
        clean(snapshot)==clean(caps_dict(caps)) && record["caps_profile"]==declaration["effective_profile"] &&
            record["caps_profile_validation"]=="validated_bounded_profile" || error("recorded effective caps metadata differs")
        out["caps_declaration"]=declaration;out["caps"]=caps_dict(caps)
        out["caps_profile"]=declaration["effective_profile"]
        out["caps_validation_status"]="validated_bounded_profile"
        # Bind executable sources, not just the purported mathematical values.
        expected=["study.jl","independent_check.jl"]
        append!(expected,["src/"*name for name in sort(readdir(joinpath(@__DIR__,"src"))) if endswith(name,".jl")])
        Set(keys(record["source_hashes"]))==Set(expected) || error("source inventory differs")
        for (path,digest) in record["source_hashes"]
            isfile(joinpath(@__DIR__,path)) && filehash(joinpath(@__DIR__,path))==digest || error("source identity changed: "*path)
        end
        input_bound=artifact_binding!(out,record,"input",input,input_admission)
        kernel_bound=artifact_binding!(out,record,"kernel",kernel,kernel_admission)
        bound=input_bound && kernel_bound
        out["artifact_links_authenticated"]=bound
        out["artifact_binding_status"]=bound ? "authenticated" : "unbound"
        if record["complete"]===false
            record["positive_full_bound"]===false && record["reason"] isa AbstractString &&
                !isempty(record["reason"]) || error("invalid refusal receipt")
            out["status"]=bound ? "recorded_refusal" : "recorded_refusal_unbound"
            out["producer_refusal_status"]=record["status"];out["producer_refusal_reason"]=record["reason"]
        else
            record["complete"]===true || error("producer complete must be Boolean")
            bound || error("completed result requires available authenticated input and kernel")
            out["scientific_verification_performed"]=true
            checked=independent_full_check(input,kernel,record;caps,guard)
            merge!(out,Dict(string(k)=>v for (k,v) in pairs(checked)))
            out["complete"]=true;out["positive_certificate_verified"]=checked.L_Q>0
        end
    catch e
        caps_attempted && out["caps_validation_status"]=="not_checked" && (out["caps_validation_status"]="rejected")
        out["status"]=e isa Refusal ? e.status : "verification_failed";out["reason"]=sprint(showerror,e)
    end
    out["CPU_seconds"]=ccall(:clock,Int64,())/1.0e6-began
    out["wall_seconds"]=(time_ns()-wall)/1e9
    write_final_record(output,out;caps)
end
function main(args=ARGS)
    profile,positional=parse_caps_cli(args)
    !isempty(positional) && length(positional)%4==0 || error("usage: independent_check.jl [--caps-profile default|larger14] INPUT KERNEL RESULT OUTPUT [quadruples...]; '-' denotes an unavailable artifact for unbound refusal metadata only")
    okay=true
    for i in 1:4:length(positional)
        out=check_one(positional[i:i+3]...;caps_profile=profile)
        okay &= out["status"] in ("independently_verified","recorded_refusal")
    end
    okay ? 0 : 2
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(main())
end
