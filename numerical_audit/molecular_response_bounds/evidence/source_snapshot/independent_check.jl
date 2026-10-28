#!/usr/bin/env julia
using TOML, SHA
include(joinpath(@__DIR__,"src","MolecularResponseBounds.jl"))
using .MolecularResponseBounds.BoundsCommon
using .MolecularResponseBounds.IndependentCheck

function check_one(input,kernel,result,output;caps_profile=nothing)
    caps=Caps();filesize(result)<=caps.max_output_bytes || error("result byte cap")
    record=TOML.parsefile(result)
    began=ccall(:clock,Int64,())/1.0e6;wall=time_ns()
    guard=()->begin
        ccall(:clock,Int64,())/1.0e6-began<=100 || throw(Refusal("cpu_limit","independent replay100CPU cap"))
        (time_ns()-wall)/1e9<=500 || throw(Refusal("wall_limit","independent replay500wall cap"))
    end
    out=Dict{String,Any}("schema"=>"independent-molecular-bound-v1","complete"=>false,
        "result_sha256"=>filehash(result),"input_sha256"=>filehash(input),"kernel_sha256"=>filehash(kernel),
        "checker_sha256"=>filehash(@__FILE__),"caps_validation_status"=>"not_checked",
        "positive_certificate_verified"=>false,"scientific_verification_performed"=>false)
    try
        # Only whitelisted bounded requests construct replay Caps. Every frozen
        # memory/work/bit/register value must equal the hard-coded profile.
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
        if record["complete"]===false
            record["positive_full_bound"]===false && !isempty(record["reason"]) || error("invalid refusal receipt")
            out["status"]="recorded_refusal";out["positive_certificate_verified"]=false
            out["producer_refusal_status"]=record["status"];out["producer_refusal_reason"]=record["reason"]
        else
            out["scientific_verification_performed"]=true
            checked=independent_full_check(input,kernel,record;caps,guard)
            merge!(out,Dict(string(k)=>v for (k,v) in pairs(checked)))
            out["complete"]=true;out["positive_certificate_verified"]=checked.L_Q>0
        end
    catch e
        out["caps_validation_status"]=="not_checked" && (out["caps_validation_status"]="rejected")
        out["status"]=e isa Refusal ? e.status : "verification_failed";out["reason"]=sprint(showerror,e)
    end
    out["CPU_seconds"]=ccall(:clock,Int64,())/1.0e6-began
    out["wall_seconds"]=(time_ns()-wall)/1e9
    write_record(output,out;caps);out
end
function main(args=ARGS)
    profile,positional=parse_caps_cli(args)
    !isempty(positional) && length(positional)%4==0 || error("usage: independent_check.jl [--caps-profile default|larger14] INPUT KERNEL RESULT OUTPUT [quadruples...]")
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
