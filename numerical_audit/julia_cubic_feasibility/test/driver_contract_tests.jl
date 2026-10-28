using Test, TOML
if !isdefined(@__MODULE__,:CubicCheckerEntrypointSmoke)
    include("entrypoint_tests.jl")
end

function driver_records()
    # The expected main inventory is explicit, independent of the producer helper.
    schedule=(("enumeration","first"),("dp","first"),("dp","warm1"),
        ("enumeration","warm1"),("enumeration","warm2"),("dp","warm2"),
        ("dp","warm3"),("enumeration","warm3"))
    [Dict{String,Any}("method"=>method,"stage"=>stage,"complete"=>true,
        "cpu_seconds"=>1.0,"wall_seconds"=>1.0,"flat_sha256"=>"fixture",
        "solver"=>(status="complete",complete=true,reason="",permutations=6)) for (method,stage) in schedule]
end
function checker_fixture(records;mode="main",warms=3)
    Dict{String,Any}("status"=>"complete","mode"=>mode,
        "requested_warm_repetitions"=>warms,"width"=>3,"records"=>
        CubicDPEntrypointSmoke.StudyIO.clean(records))
end

@testset "Benchmark completion is the full method/stage inventory" begin
    producer=CubicDPEntrypointSmoke;checker=CubicCheckerEntrypointSmoke
    records=driver_records()
    @test producer.benchmark_status(records,3,3)=="complete"
    @test checker.checked_run_inventory(checker_fixture(records)).expected_records==8
    for damaged in (records[1:end-1],[records;records[end]],records[[2,1,3,4,5,6,7,8]])
        @test producer.benchmark_status(damaged,3,3)=="incomplete"
        @test_throws ErrorException checker.checked_run_inventory(checker_fixture(damaged))
    end
    warmrefusal=deepcopy(records)
    warmrefusal[7]["complete"]=false
    warmrefusal[7]["solver"]=(status="bag_state_limit",complete=false,reason="cap")
    @test producer.benchmark_status(warmrefusal,3,3)=="incomplete"
    @test_throws ErrorException checker.checked_run_inventory(checker_fixture(warmrefusal))
    partialenum=deepcopy(records)
    partialenum[8]["solver"]=merge(partialenum[8]["solver"],(permutations=5,))
    @test producer.benchmark_status(partialenum,3,3)=="incomplete"
    @test_throws ErrorException checker.checked_run_inventory(checker_fixture(partialenum))
    firstonly=records[1:2]
    @test producer.benchmark_status(firstonly,3,0)=="complete"
    validation=checker_fixture(firstonly;mode="validation",warms=0)
    @test checker.checked_run_inventory(validation).warms==0
    summary=checker.warm_summary(validation["records"],0)
    @test !summary["dp"]["available"] && !summary["enumeration"]["available"]
    @test !haskey(summary,"enumeration_over_DP_CPU_median_ratio")
    legacy=checker_fixture(records);delete!(legacy,"mode");delete!(legacy,"requested_warm_repetitions")
    @test checker.checked_run_inventory(legacy).warms==3
    legacy["records"]=validation["records"]
    @test_throws ErrorException checker.checked_run_inventory(legacy)
    @test_throws ErrorException checker.checked_run_inventory(checker_fixture(firstonly;warms=0))
    missing_count=checker_fixture(records);delete!(missing_count,"requested_warm_repetitions")
    @test_throws ErrorException checker.checked_run_inventory(missing_count)
end

@testset "First DP refusal is explicit and carries no DP optimum" begin
    producer=CubicDPEntrypointSmoke;checker=CubicCheckerEntrypointSmoke
    records=filter(item->item["method"]!="dp" || item["stage"]=="first",driver_records())
    records[2]["complete"]=false
    records[2]["solver"]=(status="bag_state_limit",complete=false,reason="cap",order=Int[],numerator=nothing,score=nothing)
    @test producer.benchmark_status(records,3,3)=="completed_with_DP_refusal"
    main=checker_fixture(records);main["status"]="completed_with_DP_refusal";main["DP_optimum_certified"]=false
    @test checker.checked_run_inventory(main).refused
    claimed=deepcopy(main);claimed["status"]="complete"
    @test_throws ErrorException checker.checked_run_inventory(claimed)
    claimed=deepcopy(main);claimed["records"][2]["solver"]["order"]=[1,2,3]
    @test_throws ErrorException checker.checked_run_inventory(claimed)
    claimed=deepcopy(main);claimed["records"][2]["solver"]["status"]="error"
    @test_throws ErrorException checker.checked_run_inventory(claimed)
end

@testset "Driver errors retain checkpoints and an incomplete final record" begin
    producer=CubicDPEntrypointSmoke
    mktempdir() do directory
        output=joinpath(directory,"case.toml");records=Dict{String,Any}[]
        @test_throws ErrorException producer.checkpoint_sample("dp","warm2","fixture",output,records) do guard
            error("method CPU deadline")
        end
        checkpoint=output*".warm2_dp.toml"
        saved=read(checkpoint)
        failure=TOML.parsefile(checkpoint)
        @test !failure["complete"] && failure["solver"]["status"]=="error"
        @test occursin("method CPU deadline",failure["solver"]["reason"])
        result=producer.one_case(joinpath(directory,"absent input"),joinpath(directory,"absent prior"),output;repetitions=0)
        @test result["status"]=="incomplete"
        @test result["mode"]=="validation" && result["requested_warm_repetitions"]==0
        @test TOML.parsefile(output)["status"]=="incomplete"
        @test read(checkpoint)==saved
        final_saved=read(output)
        @test_throws ErrorException producer.one_case("absent","absent",output;repetitions=0)
        @test read(output)==final_saved
    end
end

@testset "Capped diagnostics complete execution and batch visits both cases" begin
    window=WindowEntrypointSmoke
    audit=Dict{String,Any}("status"=>"complete","scope_complete"=>true,
        "bucket_status"=>"complete","decomposition_complete"=>true,"decomposition_checked"=>true)
    @test window.diagnostic_completion(audit).scan_complete
    for (key,value) in (("status","scope_or_atom_cap"),("scope_complete",false),
            ("bucket_status","edge_work_limit"),("decomposition_complete",false))
        partial=copy(audit);partial[key]=value
        @test window.diagnostic_completion(partial).status=="diagnostic_incomplete"
    end
    visits=String[]
    run_case=(args;exit_on_failure=true)->begin
        push!(visits,args[2])
        Dict{String,Any}("execution_status"=>"complete","status"=>"diagnostic_incomplete")
    end
    rows=WindowBatchEntrypointSmoke.batch_main(["structure","first","out1","second","out2"];run_case=run_case)
    @test visits==["first","second"] && length(rows)==2
end

@testset "Window and physical consumers check the actual insertion contract" begin
    window=WindowEntrypointSmoke
    flat=window.WindowSearch.FlatWindow([2,3],UInt64[1],ComplexF64[1],window.WindowSearch.FlatFactor[],1e-14)
    optimum=(status="complete",permutations=2,minimizers=2,order=[2,3],numerator=big(1)//big(1),score=1.0)
    insertion=(order=[3,2],numerator=big(1)//big(1),score=1.0,evaluations=3,moves=0)
    @test !hasproperty(insertion,:status)
    @test isnothing(window.validate_window_solvers(flat,optimum,insertion))
    @test_throws ErrorException window.validate_window_solvers(flat,merge(optimum,(permutations=1,)),insertion)
    @test_throws ErrorException window.validate_window_solvers(flat,optimum,merge(insertion,(score=2.0,)))
    physical=PhysicalEntrypointSmoke
    case=Dict{String,Any}("status"=>"complete","mode"=>"window","baseline_full_order"=>[1,2,3,4],
        "baseline_window"=>[2,3],"window_start"=>2,"window_width"=>2,"steps"=>2,
        "new_full_order"=>[1,3,2,4],"new_J"=>0.25,
        "optimum"=>Dict("status"=>"complete","permutations"=>2,"minimizers"=>1,"order"=>[3,2],"numerator"=>"4//1","score"=>4.0),
        "insertion"=>Dict("order"=>[2,3],"numerator"=>"9//1","score"=>9.0))
    data=(order=[1,2,3,4],case="fixture",source=Dict("source_model_sha256"=>"model","cubic_source_sha256"=>"parent"))
    orders=physical.cubic_orders(case,data)
    @test orders["window_exact"]==[1,3,2,4] && orders["window_insertion"]==data.order
    @test !haskey(case["insertion"],"status")
    damaged=deepcopy(case);damaged["optimum"]["status"]="incomplete"
    @test_throws ErrorException physical.cubic_orders(damaged,data)
    damaged=deepcopy(case);damaged["new_full_order"]=[4,3,2,1]
    @test_throws ErrorException physical.cubic_orders(damaged,data)
    damaged=deepcopy(case);damaged["insertion"]["score"]=8.0
    @test_throws ErrorException physical.cubic_orders(damaged,data)
    mktempdir() do directory
        skips=String[];resultpath=joinpath(directory,"window.toml")
        physical.continuation_control!(orders,skips,resultpath,data)
        @test length(skips)==1 && occursin("optional continuation file absent",only(skips))
        @test !haskey(orders,"legacy_continuation")
        control=Dict("model_sha256"=>"model","parent_sha256"=>"parent","order"=>[1,2,3,4],"status"=>"cpu_limit")
        open(joinpath(directory,"fixture.CONTINUATION.a1.toml"),"w") do io
            TOML.print(io,control)
        end
        empty!(skips);physical.continuation_control!(orders,skips,resultpath,data)
        @test isempty(skips) && orders["legacy_continuation"]==data.order
    end
end

@testset "Archive validation already includes the declared expected minimum" begin
    archived=ArchiveEntrypointSmoke.ArchivedWindows
    raw=Dict{String,Any}("case"=>"fixture","width"=>2,"window"=>[1,2],"warm_order"=>[0,1,2,3],
        "outside_prefix"=>[0],"outside_suffix"=>[3],"objective_denominator"=>"1",
        "pairs"=>["1,2,1,1"],"integer_weights"=>["1,1"],"expected_minimum"=>"1")
    data=archived.ArchiveData(raw)
    @test archived.validate_archive_solution(data,[1,2],1).valid
    raw["expected_minimum"]="0";wrong=archived.ArchiveData(raw)
    validation=archived.validate_archive_solution(wrong,[1,2],1)
    @test validation.score_matches && !validation.expected_matches && !validation.valid
end
