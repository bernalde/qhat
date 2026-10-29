# Mechanics only: an analytic14-mode diagonal model, not a molecular panel.
using Test, TOML, SHA
include(joinpath(@__DIR__,"..","physical.jl"))

@testset "Primary driver admits the whole frozen matrix before physical work" begin
    mktempdir() do directory
        input=joinpath(directory,"model.toml")
        data=Dict("case"=>"mechanics_only","n"=>14,"electrons"=>4,
            "denominator"=>"1","partition"=>"conserving","source_json_sha256"=>"0"^64,
            "baseline_order"=>[1],"units"=>[Dict("words"=>["0Z"],"numerators"=>["1"])])
        open(input,"w") do io;TOML.print(io,data);end
        rows=[Dict{String,Any}("name"=>join(string.((arm,rep,formula,time,cap)),"_"),
            "arm"=>arm,"repetition"=>rep,"formula"=>formula,"time"=>time,
            "budget_multiplier"=>cap,"partition"=>arm=="pauli" ? "pauli" : "conserving",
            "order"=>[1]) for (arm,rep,formula,time,cap) in sort!(collect(grid("preflight")))]
        orders=joinpath(directory,"orders.toml")
        manifest=Dict("mode"=>"preflight","input_sha256"=>filehash(input),"candidates"=>rows)
        open(orders,"w") do io;TOML.print(io,manifest);end
        output=joinpath(directory,"physical")
        result=run_case(input,orders,output)
        @test result["status"]=="complete" && result["physical_panel_complete"]
        @test result["reference_sector_dimension"]==1001
        @test length(result["rows"])==44
        # Conserving formula 1/2 and Pauli formula 2, each at both times.
        @test result["unique_actual_propagations"]==6
        @test all(row->row["steps"]==row["ideal_rotations"]==1,result["rows"])
        @test all(row->row["projective_error"]<1e-26,result["rows"])
        @test all(row->filesize(joinpath(output,row["actual_file"]))==16*(1<<14),result["rows"])
        @test_throws ErrorException run_case(input,orders,output)
        pop!(manifest["candidates"])
        bad=joinpath(directory,"missing.toml")
        open(bad,"w") do io;TOML.print(io,manifest);end
        refused=joinpath(directory,"refused")
        @test_throws ErrorException run_case(input,bad,refused)
        @test !ispath(refused)
    end
end
