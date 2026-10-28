module DriverFixture
include(joinpath(@__DIR__,"..","study.jl"))
end
@testset "Registry retains failed cases, phase receipts and later exact-zero outcomes" begin
    mktempdir() do dir
        input,_=source_file(dir,"driver",1,0,[(String["0X"],[1]),(String["0Z"],[1])];den=1)
        kernel=joinpath(dir,"driver.zero.kernel.toml");freeze_kernel(input,kernel;time=0.0)
        registry=Dict("schema"=>"molecular-bounds-registry-v1",
            "options"=>Dict("iterations"=>0,"max_witnesses"=>0,"target_count"=>1,"comparison_count"=>0,
                "attempt_CPU_seconds"=>30.0,"attempt_wall_seconds"=>150.0),
            "cases"=>[Dict("id"=>"bad","input"=>input,"kernel"=>kernel,"order"=>[1,1]),
                Dict("id"=>"zero","input"=>input,"kernel"=>kernel,"order"=>[1,2])])
        path=joinpath(dir,"registry.toml");write_record(path,registry)
        out=joinpath(dir,"out");summary=DriverFixture.run_registry(path,out)
        @test summary["status"]=="inventory_complete" && !summary["all_cases_completed"]
        @test length(summary["cases"])==2
        bad=TOML.parsefile(joinpath(out,"bad.RESULT.toml"))
        zero=TOML.parsefile(joinpath(out,"zero.RESULT.toml"))
        @test !bad["complete"] && bad["status"]=="error" && !bad["positive_full_bound"]
        @test isfile(joinpath(out,"bad.START.toml")) && isfile(joinpath(out,"bad.PHASE.model.toml"))
        @test zero["complete"] && zero["status"]=="zero_full_bound" && !zero["positive_full_bound"]
        @test zero["ratio"]=="unavailable"
        @test independent_full_check(input,kernel,zero).L_Q==0
        @test isfile(joinpath(out,"zero.PHASE.full_UB.toml")) && isfile(joinpath(out,"RUN.RESULT.toml"))
        @test_throws ErrorException DriverFixture.run_registry(path,out)
    end
end
