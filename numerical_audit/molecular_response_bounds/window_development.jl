#!/usr/bin/env julia
using TOML, SHA
include(joinpath(@__DIR__,"src","MolecularResponseBounds.jl"))
using .MolecularResponseBounds.BoundsCommon
using .MolecularResponseBounds.ResponseBounds
using .MolecularResponseBounds.StudyAPI
using .MolecularResponseBounds.IndependentCheck
module FrozenWindowConstruction
include(joinpath(@__DIR__,"..","julia_cubic_feasibility","src","ResponseCore.jl"))
include(joinpath(@__DIR__,"..","julia_cubic_feasibility","src","StudyIO.jl"))
include(joinpath(@__DIR__,"..","julia_cubic_feasibility","src","WindowSearch.jl"))
end
# Byte-for-byte serialization contract from the retained cubic-DP study.
# This digest binds any later comparison with its archived exhaustive oracle.
function flat_digest(flat)
    io=IOBuffer()
    println(io,"flat-window-v1")
    println(io,join(flat.labels,','));println(io,join(flat.destinations,','))
    for value in flat.baseline
        println(io,bitstring(real(value)),bitstring(imag(value)))
    end
    for factor in flat.factors
        println(io,join(factor.scope,','));println(io,join(factor.destinations,','))
        println(io,size(factor.values))
        for value in factor.values
            println(io,bitstring(real(value)),bitstring(imag(value)))
        end
    end
    bytes2hex(sha256(take!(io)))
end
function development(input,start,output;comparisons=1)
    caps=Caps();began=ccall(:clock,Int64,())/1.0e6;wall=time_ns()
    guard=()->begin
        ccall(:clock,Int64,())/1.0e6-began<=20 || throw(Refusal("cpu_limit","window development20CPU cap"))
        (time_ns()-wall)/1e9<=100 || throw(Refusal("wall_limit","window development100wall cap"))
    end
    data=FrozenWindowConstruction.StudyIO.load_model(input)
    model=FrozenWindowConstruction.ResponseCore.build_window_model(data.units,data.order,data.hf,start:start+7)
    flat=FrozenWindowConstruction.WindowSearch.prepare_window(model)
    old=joinpath(@__DIR__,"..","julia_cubic_feasibility","src","ResponseCore.jl")
    record=window_record(flat,data.order,start:start+7,length.(getfield.(data.units,:coeff)),30data.words;
        input_sha256=data.input_sha256,source_sha256=filehash(old),comparison_count=comparisons,caps,guard)
    r=encode_window(flat;binding=record["binding"],caps,guard)
    replay=independent_window_check(r,record,input;source_path=old,caps,guard)
    replay.L_Q==record["L_Q"] || error("window independent bound differs")
    record["independent_exact_arithmetic_verified"]=true
    record["flat_sha256"]=flat_digest(flat)
    record["fresh_enumeration_performed"]=false
    record["CPU_seconds"]=ccall(:clock,Int64,())/1.0e6-began
    record["wall_seconds"]=(time_ns()-wall)/1e9
    write_record(output,record;caps);record
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==3 || error("usage: window_development.jl INPUT START NEW_OUTPUT")
    development(ARGS[1],parse(Int,ARGS[2]),ARGS[3])
end
