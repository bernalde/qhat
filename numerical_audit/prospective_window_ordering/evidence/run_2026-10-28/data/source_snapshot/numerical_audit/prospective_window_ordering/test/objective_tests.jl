using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"..","objective_check.jl"))
const OC=ObjectiveCheck

function oc_input(path; n=1,electrons=0,denominator="100",units=nothing)
    units=isnothing(units) ? [Dict("words"=>["0X"],"numerators"=>["100"]),
        Dict("words"=>["0Z"],"numerators"=>["100"])] : units
    data=Dict("case"=>"objective_mechanics_only","n"=>n,"electrons"=>electrons,
        "partition"=>"conserving","denominator"=>denominator,"source_json_sha256"=>"0"^64,
        "baseline_order"=>collect(1:length(units)),"units"=>units)
    open(path,"w") do io;TOML.print(io,data);end
    data
end
function oc_word_key(word)
    Tuple((parse(Int,token[1:end-1]),token[end]) for token in (isempty(word) ? String[] : split(word,',')))
end
function oc_dense_word(word,n)
    matrices=Dict('X'=>ComplexF64[0 1;1 0],'Y'=>ComplexF64[0 -im;im 0],'Z'=>ComplexF64[1 0;0 -1])
    key=Dict(oc_word_key(word));answer=ones(ComplexF64,1,1)
    for site in n-1:-1:0
        answer=kron(answer,haskey(key,site) ? matrices[key[site]] : Matrix{ComplexF64}(I,2,2))
    end
    answer
end
function oc_dense_units(source,partition)
    n=source["n"];den=parse(BigInt,source["denominator"])
    terms=[(word=word,coefficient=Float64(parse(BigInt,a)//den))
        for block in source["units"] for (word,a) in zip(block["words"],block["numerators"])]
    if partition=="pauli"
        sort!(terms;by=t->oc_word_key(t.word))
        return [t.coefficient*oc_dense_word(t.word,n) for t in terms]
    end
    [sum((Float64(parse(BigInt,a)//den)*oc_dense_word(word,n)
        for (word,a) in zip(block["words"],block["numerators"]));init=zeros(ComplexF64,1<<n,1<<n))
        for block in source["units"]]
end
function oc_dense_first(source,partition,order)
    n=source["n"];hf=(1<<source["electrons"])-1;units=oc_dense_units(source,partition)
    column=zeros(ComplexF64,1<<n);initial=zeros(ComplexF64,1<<n);initial[hf+1]=1
    for i in 1:length(order)-1,j in i+1:length(order)
        a,c=units[order[i]],units[order[j]]
        column.+=0.5.*((a*c-c*a)*initial)
    end
    column[hf+1]=0
    original=oc_dense_units(source,"conserving")
    diagonal=sum((u for u in original if isdiag(u));init=zeros(ComplexF64,1<<n,1<<n))
    energies=real.(diag(diagonal));response=copy(column)
    for state in 0:(1<<n)-1
        x=(energies[state+1]-energies[hf+1])/2
        response[state+1]*=(iszero(x) ? 1.0 : sin(x)/x)*cis(x)
    end
    (value=sum(abs2,response),column=column)
end
function oc_orders(labels)
    isempty(labels) && return [Int[]]
    [[x;tail] for x in labels for tail in oc_orders([y for y in labels if y!=x])]
end
function oc_candidate(name,kind,score;order=[1,2],partition="conserving",time=0.4,cap=8)
    row=Dict{String,Any}("name"=>name,"arm"=>"mechanics","repetition"=>1,"formula"=>2,
        "time"=>time,"budget_multiplier"=>cap,"partition"=>partition,"order"=>order,"score_kind"=>kind)
    !isnothing(score) && (row["constructor_score"]=score)
    row
end
function oc_manifest(path,input,rows;hash=nothing)
    data=Dict("mode"=>"preflight","input_sha256"=>(isnothing(hash) ? bytes2hex(sha256(read(input))) : hash),
        "candidates"=>rows)
    open(path,"w") do io;TOML.print(io,data);end
    data
end

@testset "Analytic X/Z fixes first pair sign, sinc and absence of r normalization" begin
    mktempdir() do dir
        path=joinpath(dir,"model.toml");oc_input(path)
        data=OC.WindowStudyIO.load_model(path);model=OC.prepare_model(data)
        first=OC.first_response(model,"conserving",[1,2])
        reversed=OC.first_response(model,"conserving",[2,1])
        @test first.column[UInt64(1)]//first.column_denominator==1
        @test reversed.column[UInt64(1)]//reversed.column_denominator==-1
        @test first.value≈sin(1.0)^2 rtol=2e-14
        @test reversed.value≈first.value rtol=2e-14
        @test first.steps===nothing && first.r_normalized===false
        cubic=OC.replay_objective(model,"conserving",[1,2],"diagonal_reference_cubic_response")
        @test cubic.steps==29
        @test cubic.value≈sin(1.0)^2/(9*29^4) rtol=2e-12
        @test first.value/cubic.value≈9*29^4 rtol=2e-12
        @test cubic.details["r_normalized"]===true
        @test_throws OC.ObjectiveLimit OC.first_response(model,"conserving",[1,2];budget=OC.ReplayBudget(()->nothing;max_work=1))
        @test_throws ArgumentError OC.first_response(model,"conserving",[1,1])
        @test_throws ArgumentError OC.replay_objective(model,"conserving",[1,2],"unknown")
    end
end

@testset "Independent literal dense pair expectation, conserving and structured Pauli orders" begin
    mktempdir() do dir
        path=joinpath(dir,"model.toml")
        source=oc_input(path;n=3,electrons=1,units=[
            Dict("words"=>["0X,1X","0Y,1Y"],"numerators"=>["50","50"]),
            Dict("words"=>["1X,2X","1Y,2Y"],"numerators"=>["35","35"]),
            Dict("words"=>["0Z"],"numerators"=>["20"])])
        model=OC.prepare_model(OC.WindowStudyIO.load_model(path))
        for partition in ("conserving","pauli")
            count=length(OC.partition_units(model,partition).blocks)
            for order in oc_orders(collect(1:count))
                actual=OC.first_response(model,partition,order)
                expected=oc_dense_first(source,partition,order)
                @test actual.value≈expected.value atol=2e-14 rtol=2e-12
                for state in 0:7
                    actual_amplitude=Float64(get(actual.column,UInt64(state),big(0))//actual.column_denominator)
                    @test actual_amplitude≈expected.column[state+1] atol=2e-14 rtol=2e-12
                end
            end
        end
        # This14-mode fixture tests label parsing only, never allocates a14-mode matrix.
        high=joinpath(dir,"labels.toml")
        oc_input(high;n=14,electrons=4,units=[
            Dict("words"=>["10Z"],"numerators"=>["10"]),
            Dict("words"=>["2Z"],"numerators"=>["20"]),
            Dict("words"=>["0X,13X"],"numerators"=>["30"])])
        compiled=OC.prepare_model(OC.WindowStudyIO.load_model(high));parts=OC.partition_units(compiled,"pauli")
        @test [w.key for w in compiled.sorted_words]==[((0,'X'),(13,'X')),((2,'Z'),),((10,'Z'),)]
        @test parts.floating[2].parity==UInt64[UInt64(1)<<2]
        @test parts.floating[3].parity==UInt64[UInt64(1)<<10]
    end
end

@testset "Logical physical-condition duplicates reuse only exact four-field identity" begin
    mktempdir() do dir
        input=joinpath(dir,"model.toml");oc_input(input)
        expected=sin(1.0)^2/(9*29^4)
        rows=[oc_candidate("cubic_a","diagonal_reference_cubic_response",expected),
            oc_candidate("cubic_b","diagonal_reference_cubic_response",expected;time=1.6,cap=24),
            oc_candidate("first_a","retained_first_order_response",sin(1.0)^2),
            oc_candidate("first_b","retained_first_order_response",sin(1.0)^2;time=1.6,cap=24),
            oc_candidate("unscored","unavailable",nothing),
            oc_candidate("distinct_score","diagonal_reference_cubic_response",nextfloat(expected)),
            oc_candidate("different_partition","retained_first_order_response",sin(1.0)^2;partition="pauli"),
            oc_candidate("different_order","retained_first_order_response",sin(1.0)^2;order=[2,1])]
        candidates=joinpath(dir,"candidates.toml");oc_manifest(candidates,input,rows)
        output=joinpath(dir,"check.toml");record=OC.run_check(input,candidates,output)
        @test record["complete"] && record["all_scored_objectives_verified"]
        @test record["unique_objectives"]==5 && record["scored_logical_rows"]==7 && record["unscored_rows"]==1
        @test record["logical_rows"][2]["check_index"]==record["logical_rows"][1]["check_index"]
        @test record["logical_rows"][4]["check_index"]==record["logical_rows"][3]["check_index"]
        @test record["logical_rows"][6]["check_index"]!=record["logical_rows"][1]["check_index"]
        @test record["logical_rows"][7]["check_index"]!=record["logical_rows"][3]["check_index"]
        @test record["logical_rows"][8]["check_index"]!=record["logical_rows"][3]["check_index"]
        @test record["logical_rows"][5]["status"]=="unscored_censored"
        @test isfile(output*".START.toml") && isfile(output*".CHECK.001.RESULT.toml") && isfile(output*".PARTIAL.005.toml")
        @test TOML.parsefile(output)["complete"]
        @test_throws ErrorException OC.run_check(input,candidates,output)
        wrong=deepcopy(rows[1]);wrong["constructor_score"]=2expected
        bad=joinpath(dir,"badscore.toml");oc_manifest(bad,input,[wrong])
        mismatch=OC.run_check(input,bad,joinpath(dir,"mismatch.toml"))
        @test !mismatch["complete"] && mismatch["status"]=="score_mismatch"
        @test mismatch["checks"][1]["status"]=="mismatch"
        hashed=joinpath(dir,"badhash.toml");oc_manifest(hashed,input,[rows[1]];hash="0"^64)
        @test !OC.run_check(input,hashed,joinpath(dir,"hashcheck.toml"))["complete"]
    end
end

@testset "Floors, numerical refusals and malformed candidates retain no artificial objective" begin
    mktempdir() do dir
        input=joinpath(dir,"commuting.toml")
        oc_input(input;units=[Dict("words"=>["0Z"],"numerators"=>["100"]),
            Dict("words"=>[""],"numerators"=>["70"])])
        candidate=joinpath(dir,"floor.toml");oc_manifest(candidate,input,[oc_candidate("floor","retained_first_order_response",0.0)])
        floor=OC.run_check(input,candidate,joinpath(dir,"floorcheck.toml"))
        @test !floor["complete"] && floor["floor_censored_objectives"]==1
        @test floor["checks"][1]["status"]=="floor_censored" && floor["checks"][1]["within_tolerance"]
        tiny=joinpath(dir,"tiny.toml")
        oc_input(tiny;denominator=string(big(2)^500),units=[Dict("words"=>["0X"],"numerators"=>["1"]),
            Dict("words"=>["0Z"],"numerators"=>["1"])])
        tinyrows=joinpath(dir,"tinyrows.toml");oc_manifest(tinyrows,tiny,[oc_candidate("underflow","retained_first_order_response",0.0),
            oc_candidate("cubic_underflow","diagonal_reference_cubic_response",0.0)])
        refused=OC.run_check(tiny,tinyrows,joinpath(dir,"refused.toml"))
        @test !refused["complete"] && refused["checks"][1]["status"]=="numerical_refusal"
        @test !haskey(refused["checks"][1],"replayed_score")
        @test refused["checks"][2]["status"]=="numerical_refusal" && !haskey(refused["checks"][2],"replayed_score")
        valid=joinpath(dir,"valid.toml");oc_input(valid)
        for (i,row) in enumerate([oc_candidate("missing_kind","unknown",0.7),
                oc_candidate("bad_order","retained_first_order_response",0.7;order=[1,1]),
                oc_candidate("bad_partition","retained_first_order_response",0.7;partition="unknown"),
                oc_candidate("bad_score","retained_first_order_response",NaN)])
            path=joinpath(dir,"bad_"*string(i)*".toml");oc_manifest(path,valid,[row])
            record=OC.run_check(valid,path,joinpath(dir,"badcheck_"*string(i)*".toml"))
            @test !record["complete"] && record["logical_rows"][1]["status"]=="invalid_candidate"
        end
        unscored=joinpath(dir,"unscored.toml");oc_manifest(unscored,valid,[oc_candidate("none","unavailable",nothing)])
        empty=OC.run_check(valid,unscored,joinpath(dir,"unscoredcheck.toml"))
        @test !empty["complete"] && empty["status"]=="no_scored_orders_supplied"
        callback=OC.run_check(valid,unscored,joinpath(dir,"budget.toml");
            check_budget=()->throw(OC.ObjectiveLimit("injected_refusal","fixture limit")))
        @test !callback["complete"] && callback["status"]=="injected_refusal"
        oversized=joinpath(dir,"oversized.toml")
        oc_input(oversized;units=[Dict("words"=>fill("0Z",8193),"numerators"=>fill("1",8193))])
        manifest=joinpath(dir,"oversizedrows.toml");oc_manifest(manifest,oversized,[oc_candidate("none","unavailable",nothing;order=[1])])
        @test !OC.run_check(oversized,manifest,joinpath(dir,"oversizedcheck.toml"))["complete"]
    end
    @test OC.compare_scores(1.0,1.0+5e-9).verified
    @test !OC.compare_scores(1.0,1.0+5e-7).within_tolerance
    @test OC.compare_scores(0.0,0.0).numerical_floor
    @test_throws ArgumentError OC.compare_scores(Inf,1.0)
end
