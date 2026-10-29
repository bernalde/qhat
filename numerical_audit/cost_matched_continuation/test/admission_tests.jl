using Test
include(joinpath(@__DIR__, "..", "physical.jl"))

# Admission only: no model compilation, reference, optimization or propagation.
@testset "Single-line scientific TOML boundary" begin
    mktempdir() do folder
        path = joinpath(folder, "scientific.toml")
        text = "case = \"A\"\nword = '0Z'\nnumerators = [\"1\", \"-2\"]\n"
        text *= "# Comment quotes are ignored: " * repeat("'", 3) * " " * repeat("\"", 3) * "\n"
        text *= "escaped = \"" * repeat("\\\"", 3) * "\"\n"
        write(path, text)
        parsed = CostMatchedPhysical.bounded_toml(path)
        @test parsed.data["case"] == "A"
        @test parsed.data["word"] == "0Z"
        @test parsed.data["escaped"] == repeat("\"", 3)
        for quotechar in ("\"", "'"), value in ("A", "", "\nA\n")
            write(path, "case = " * repeat(quotechar, 3) * value * repeat(quotechar, 3) * "\n")
            @test_throws ArgumentError CostMatchedPhysical.bounded_toml(path)
        end
    end
end

include(joinpath(@__DIR__, "..", "objective.jl"))
include(joinpath(@__DIR__, "..", "verify.jl"))

@testset "All-scored adapter permits objective deduplication" begin
    result = Dict{String,Any}("complete"=>true, "all_scored_objectives_verified"=>true,
        "unscored_rows"=>0, "processed_logical_rows"=>20, "scored_logical_rows"=>20,
        "expected_logical_rows"=>20, "unique_objectives"=>4)
    @test CostMatchedObjective.complete_scored_panel(result)
    for changes in (Dict("unscored_rows"=>1, "scored_logical_rows"=>19),
                    Dict("processed_logical_rows"=>19), Dict("complete"=>false),
                    Dict("all_scored_objectives_verified"=>false))
        @test !CostMatchedObjective.complete_scored_panel(merge(result, changes))
    end
end

@testset "Both candidate readers bind the legacy score alias" begin
    mktempdir() do folder
        input = joinpath(folder, "model.toml")
        candidate = joinpath(folder, "candidates.toml")
        source = repeat("a",64)
        units = [Dict("words"=>[string(q)*"Z"], "numerators"=>["1"]) for q in 0:5]
        raw = Dict("n"=>14, "electrons"=>4, "partition"=>"conserving", "case"=>"A",
                   "source_json_sha256"=>source, "units"=>units)
        open(input,"w") do io; CostMatchedPhysical.TOML.print(io,raw); end
        input_hash = CostMatchedPhysical.filehash(input)
        rows = Dict{String,Any}[]
        for (index,(arm,cpu)) in enumerate((("parent",0),("full",32),("full",64),("window",32),("window",64)))
            for t in (.65,1.35), cap in (8,24)
                oid = arm*string(cpu)
                push!(rows,Dict("name"=>oid*string(t)*string(cap), "order_id"=>oid,
                    "arm"=>arm,"repetition"=>1,"cpu_budget"=>cpu,"partition"=>"conserving",
                    "formula"=>2,"time"=>t,"budget_multiplier"=>cap,"order"=>collect(1:6),
                    "constructor_sha256"=>repeat(string(index),64),"parent_sha256"=>repeat("1",64),
                    "claimed_score"=>1.0,"constructor_score"=>1.0,"score_kind"=>"diagonal_reference_cubic_response"))
            end
        end
        panel = Dict("version"=>1,"case"=>"A","mode"=>"preflight","repetitions"=>1,
            "expected_logical_rows"=>20,"input_sha256"=>input_hash,"source_json_sha256"=>source,
            "protocol_sha256"=>source,"source_freeze_sha256"=>source,"candidates"=>rows)
        model = (n=14,electrons=4,partition="conserving",case="A",input_sha256=input_hash,
                 source_json_sha256=source,units=units)
        save(value) = open(candidate,"w") do io; CostMatchedPhysical.TOML.print(io,value); end
        save(panel)
        @test CostMatchedPhysical.admit(input,candidate).manifest["expected_logical_rows"]==20
        @test CostMatchedVerify.admit_candidates(model,candidate).data["expected_logical_rows"]==20
        for fault in (nothing,"unavailable",.5,true,Inf)
            broken=deepcopy(panel)
            fault===nothing ? delete!(broken["candidates"][1],"constructor_score") :
                (broken["candidates"][1]["constructor_score"]=fault)
            save(broken)
            @test_throws ArgumentError CostMatchedPhysical.admit(input,candidate)
            @test_throws ArgumentError CostMatchedVerify.admit_candidates(model,candidate)
        end
    end
end
