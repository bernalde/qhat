using Test
include(joinpath(@__DIR__,"..","objective.jl"))
include(joinpath(@__DIR__,"fixtures.jl"))
const Objective=GeneralizationObjective

function copy_witness_rows(manifest, indices=())
    copied=copy(manifest)
    copied["witnesses"]=copy(manifest["witnesses"])
    for index in indices
        copied["witnesses"][index]=copy(manifest["witnesses"][index])
    end
    copied
end

function check_objective_witness_admission()
    @testset "all persisted proposed-best witnesses, native refusal and score floor" begin
        candidate=fixture_manifest()
        candidate["protocol_path"]="source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json"
        candidate["source_freeze_path"]="INPUT_PROTOCOL_FROZEN.json"
        candidate["constructor_artifacts"]=Any[]
        original=first(filter(row->row["arm"]=="annealing" && row["cpu_budget"]==64,candidate["candidates"]))
        witness=Dict{String,Any}(key=>original[key] for key in ("order_id","arm","repetition","cpu_budget",
            "constructor_sha256","parent_sha256","partition","score_kind","order","claimed_score"))
        merge!(witness,Dict("name"=>"unpromoted.raw","role"=>"raw_proposed_best","selected"=>false,"native_score"=>1.1e-6))
        manifest=Dict{String,Any}(key=>candidate[key] for key in Objective.IOStudy.HEADER)
        for key in ("protocol_path","source_freeze_path","constructor_artifacts");manifest[key]=candidate[key];end
        manifest["witnesses"]=[witness];manifest["expected_witnesses"]=1;manifest["witness_inventory_complete"]=true
        @test length(Objective.admit_witnesses(manifest,candidate,2))==1
        bad=copy_witness_rows(manifest);duplicate=copy(manifest["witnesses"][1]);duplicate["name"]="same.claim.other.name"
        push!(bad["witnesses"],duplicate);bad["expected_witnesses"]=2
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad=copy_witness_rows(manifest,(1,));bad["witnesses"][1]["order"]=[0,1]
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad=copy_witness_rows(manifest,(1,));bad["witnesses"][1]["selected"]=true
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2) # Native selected score must join its frozen row.
        kept=Objective.objective_record(original,"selected")
        @test all(key -> kept[key]==original[key],Objective.IOStudy.FIELDS)
        @test kept["status"]=="incomplete"
        bad=copy_witness_rows(manifest);bad["expected_witnesses"]=0
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad=copy_witness_rows(manifest);bad["witness_inventory_complete"]=false
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad=copy_witness_rows(manifest,(1,));bad["witnesses"][1]["constructor_sha256"]=fixture_digest("different producer")
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad=copy_witness_rows(manifest,(1,));delete!(bad["witnesses"][1],"native_score")
        @test_throws ArgumentError Objective.admit_witnesses(bad,candidate,2)
        bad["witnesses"][1]["native_refusal"]="declared native support cap"
        @test length(Objective.admit_witnesses(bad,candidate,2))==1
        @test Objective.ObjectiveCheck.compare_scores(0.,0.).status=="floor_censored"
        @test !Objective.ObjectiveCheck.compare_scores(0.,0.).verified
        @test Objective.ObjectiveCheck.compare_scores(1e-6,1.000000001e-6).verified
        @test Objective.ObjectiveCheck.compare_scores(1e-6,2e-6).status=="mismatch"
    end
end

println("START objective admission: persisted witnesses, native refusal and score floor")
flush(stdout)
check_objective_witness_admission()
