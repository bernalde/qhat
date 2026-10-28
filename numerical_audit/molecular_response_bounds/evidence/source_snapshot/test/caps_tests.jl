@testset "Explicit Caps profiles have fixed non-count ceilings and strict schemas" begin
    small,default_decl=select_caps()
    @test caps_dict(small)==caps_dict(Caps())
    @test default_decl["effective_profile"]=="default"
    large,decl=select_caps(Dict("profile"=>"larger14"))
    @test large.max_units==1024 && large.max_words==8192
    for field in fieldnames(Caps)
        field in (:max_units,:max_words,:max_probes) && continue
        @test getfield(large,field)==getfield(small,field)
    end
    @test large.max_n==14 && large.max_bytes==128*1024^2
    @test small.max_probes==2_000_000 && large.max_probes==4_500_000
    # Admission arithmetic only: no4m-probe panel is performed in unit tests.
    worst_lookup_count=8*1024*(1024-1)÷2
    @test worst_lookup_count==4_190_208 && worst_lookup_count<=large.max_probes
    rebuilt,_=verify_caps_declaration(decl;expected_profile="larger14")
    @test caps_dict(rebuilt)==caps_dict(large)
    @test_throws ArgumentError verify_caps_declaration(decl;expected_profile="default")
    override,overdecl=select_caps(Dict("profile"=>"default");cli_profile="larger14")
    @test override.max_units==1024 && overdecl["selection_source"]=="cli"
    @test overdecl["registry"]["profile"]=="default" && overdecl["cli_profile"]=="larger14"
    limited,_=select_caps(Dict("profile"=>"larger14","max_units"=>700,"max_words"=>6000))
    @test limited.max_units==700 && limited.max_words==6000
    for config in (Dict("profile"=>"unknown"),Dict("profile"=>"larger14","max_units"=>1025),
            Dict("profile"=>"larger14","max_words"=>8193),Dict("max_units"=>513),
            Dict("max_units"=>0),Dict("max_words"=>true),Dict("max_units"=>512.0),
            Dict("max_units"=>"512"),Dict("max_n"=>15),Dict("max_bytes"=>1024^3),Dict("max_probes"=>4_500_001))
        @test_throws ArgumentError select_caps(config)
    end
    @test_throws ArgumentError select_caps(Dict("profile"=>"larger14","max_units"=>700);cli_profile="default")
    for field in ("max_units","max_words","max_n","max_bytes","max_work","max_bits","max_probes")
        changed=deepcopy(decl);changed["effective"][field]+=1
        @test_throws ArgumentError verify_caps_declaration(changed)
    end
    changed=deepcopy(decl);changed["hard_ceilings"]["max_bytes"]*=2
    @test_throws ArgumentError verify_caps_declaration(changed)
    changed=deepcopy(decl)
    # A Dict{String,Int} would convert the assigned Float64 back to Int and
    # accidentally test an unchanged valid declaration rather than bad input.
    changed["effective"]=Dict{String,Any}(changed["effective"])
    changed["effective"]["max_bytes"]=Float64(large.max_bytes)
    @test_throws ArgumentError verify_caps_declaration(changed)
    changed=deepcopy(decl);changed["unexpected"]="unbounded"
    @test_throws ArgumentError verify_caps_declaration(changed)
    @test parse_caps_cli(["--caps-profile","larger14","a","b"])==("larger14",["a","b"])
    @test parse_caps_cli(["a","b"])==(nothing,["a","b"])
    @test DriverFixture.study_arguments(["--caps-profile","larger14","a","b"]).caps_profile=="larger14"
    for args in (["--caps-profile"],["--caps-profile","unknown"],
            ["--caps-profile","default","--caps-profile","larger14"],["--max-bytes","999999999"])
        @test_throws ArgumentError parse_caps_cli(args)
    end
    mktempdir() do dir
        path=joinpath(dir,"caps.toml");write_record(path,decl)
        @test caps_dict(first(verify_caps_declaration(TOML.parsefile(path))))==caps_dict(large)
        # Reject an unknown resource override before any source/model dispatch.
        invalid_registry=Dict("schema"=>"molecular-bounds-registry-v1","caps"=>Dict("max_bytes"=>1024^3),
            "cases"=>[Dict("id"=>"unreached","input"=>joinpath(dir,"MUST-NOT-BE-READ.toml"))])
        invalid_path=joinpath(dir,"invalid-registry.toml");write_record(invalid_path,invalid_registry)
        @test_throws ArgumentError DriverFixture.run_registry(invalid_path,joinpath(dir,"uncreated"))
        @test !ispath(joinpath(dir,"uncreated"))
    end
end

@testset "Larger count profile admits complete source; smaller limits refuse without thinning" begin
    mktempdir() do dir
        # Admission fixture only: no1024-unit Taylor/atom numerical study is run.
        blocks=[(fill("0Z",8),ones(Int,8)) for _ in 1:1024]
        input,_=source_file(dir,"count-admission",1,0,blocks;den=1)
        before=filehash(input)
        small=first(select_caps());large=first(select_caps(Dict("profile"=>"larger14")))
        unit_error=try
            read_coefficients(input;caps=small);nothing
        catch e
            e
        end
        @test unit_error isa Refusal && unit_error.status=="unit_limit"
        admitted=read_coefficients(input;caps=large)
        @test length(admitted.units)==1024
        @test sum(length(u.numerators) for u in admitted.units)==8192
        @test all(u->length(u.numerators)==8,admitted.units)
        wordcaps=first(select_caps(Dict("profile"=>"larger14","max_words"=>4096)))
        word_error=try
            read_coefficients(input;caps=wordcaps);nothing
        catch e
            e
        end
        @test word_error isa Refusal && word_error.status=="word_limit"
        @test filehash(input)==before
        too_large,_=source_file(dir,"register-refusal",15,0,[(String["0Z"],[1])];den=1)
        @test_throws Refusal read_coefficients(too_large;caps=large)
    end
end

module CheckerFixture
include(joinpath(@__DIR__,"..","independent_check.jl"))
end
@testset "CLI replay validates bounded declarations before proof, including refusal receipts" begin
    mktempdir() do dir
        input,_=source_file(dir,"profile-driver",1,0,[(String["0X"],[1]),(String["0Z"],[1])];den=1)
        kernel=joinpath(dir,"zero.kernel.toml");freeze_kernel(input,kernel;time=0)
        registry=Dict("schema"=>"molecular-bounds-registry-v1","caps"=>Dict("profile"=>"larger14"),
            "options"=>Dict("iterations"=>0,"max_witnesses"=>0,"target_count"=>1,"comparison_count"=>0),
            "cases"=>[Dict("id"=>"zero","input"=>input,"kernel"=>kernel,"order"=>[1,2])])
        regpath=joinpath(dir,"registry.toml");write_record(regpath,registry)
        outdir=joinpath(dir,"out");summary=DriverFixture.run_registry(regpath,outdir)
        @test summary["all_cases_completed"] && summary["caps_profile"]=="larger14"
        resultpath=joinpath(outdir,"zero.RESULT.toml");record=TOML.parsefile(resultpath)
        @test record["caps"]["max_units"]==1024 && record["caps"]["max_words"]==8192
        @test record["input_admission_status"]=="admitted_under_declared_caps"
        phase=TOML.parsefile(joinpath(outdir,"zero.PHASE.model.toml"))
        @test phase["caps_declaration"]==record["caps_declaration"]
        checked=CheckerFixture.check_one(input,kernel,resultpath,joinpath(dir,"verified.toml");caps_profile="larger14")
        @test checked["status"]=="independently_verified" && checked["complete"]
        @test checked["caps_validation_status"]=="validated_bounded_profile"
        @test checked["caps"]["max_units"]==1024 && checked["scientific_verification_performed"]
        mismatch=CheckerFixture.check_one(input,kernel,resultpath,joinpath(dir,"mismatch.toml");caps_profile="default")
        @test mismatch["status"]=="verification_failed" && mismatch["caps_validation_status"]=="rejected"
        @test !mismatch["complete"] && !mismatch["scientific_verification_performed"]
        damaged=deepcopy(record);damaged["caps"]["max_bytes"]=1024^3
        damagedpath=joinpath(dir,"unbounded.toml");write_record(damagedpath,damaged)
        rejected=CheckerFixture.check_one(input,kernel,damagedpath,joinpath(dir,"rejected.toml"))
        @test rejected["status"]=="verification_failed" && !rejected["scientific_verification_performed"]
        refusal=deepcopy(record);refusal["complete"]=false;refusal["status"]="scope_limit"
        refusal["reason"]="fixture incomplete-coverage receipt";refusal["positive_full_bound"]=false
        refusalpath=joinpath(dir,"refusal.toml");write_record(refusalpath,refusal)
        receipt=CheckerFixture.check_one(input,kernel,refusalpath,joinpath(dir,"receipt.toml");caps_profile="larger14")
        @test receipt["status"]=="recorded_refusal" && !receipt["complete"]
        @test !receipt["positive_certificate_verified"] && !receipt["scientific_verification_performed"]
        @test receipt["caps_validation_status"]=="validated_bounded_profile"
    end
end
