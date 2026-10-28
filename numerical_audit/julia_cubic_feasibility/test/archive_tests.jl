using Test
using TOML
if !isdefined(@__MODULE__, :StructureDP)
    include(joinpath(@__DIR__, "..", "src", "StructureDP.jl"))
end
if !isdefined(@__MODULE__, :ArchivedWindows)
    include(joinpath(@__DIR__, "..", "src", "ArchivedWindows.jl"))
end
using .StructureDP
using .ArchivedWindows

function archive_fixture()
    # Window 1,2,3 is contiguous inside 0..4. Both cross-window
    # and outside-outside amplitudes interfere with the internal response.
    return Dict{String,Any}(
        "case"=>"tiny archived first-order JD", "width"=>3,
        "window"=>[1,2,3],"warm_order"=>[0,1,2,3,4],
        "outside_prefix"=>[0],"outside_suffix"=>[4],
        "objective_denominator"=>"25",
        "pairs"=>["0,1,7,7","0,4,7,-3","1,2,7,11",
                   "2,4,7,13","1,4,7,-15","1,3,8,-5","2,3,8,2"],
        "integer_weights"=>["7,3","8,5"])
end

function archive_test_orders(labels)
    isempty(labels) && return [()]
    return [(v,rest...) for v in labels for rest in archive_test_orders([x for x in labels if x!=v])]
end

@testset "Archived JD retains outside interference and exact factorization" begin
    data=ArchiveData(archive_fixture())
    @test data.window==[1,2,3]
    @test data.fixed_outside[7]==2
    @test data.integer_weights==Dict(7=>big(3),8=>big(5))
    @test archived_integer_score(data,[1,2,3])==3*13^2+5*(-3)^2
    factors=archive_factors(data)
    @test all(length(f.scope)<=4 for f in factors)
    scores=Dict{Tuple,BigInt}()
    for order in archive_test_orders(data.window)
        expected=archived_integer_score(data,order)
        @test evaluate_factors(factors,order)==expected
        scores[order]=expected
    end
    enumerated=enumeration(data)
    @test enumerated.complete && enumerated.evaluated==6
    @test enumerated.coverage_complete && enumerated.coverage_lehmer_count==6
    @test enumerated.cost==minimum(values(scores))
    @test enumerated.minimizers==count(==(enumerated.cost),values(scores))
    @test enumerated.validation.valid
    solved=solve_archive_dp(data)
    @test solved.complete && solved.cost==enumerated.cost
    @test solved.validation.valid
    @test !solved.minimizer_count_claim
    @test archive_full_order(data,solved.order)[1]==0
    @test archive_full_order(data,solved.order)[end]==4
    @test_throws ArgumentError archive_full_order(data,[4,1,2,3,0])
    wrong=validate_archive_solution(data,solved.order,solved.cost+1)
    @test !wrong.valid && !wrong.score_matches
end

@testset "Archived first-order tables can contain genuine four-label factors" begin
    raw=Dict{String,Any}("case"=>"four labels","width"=>4,
        "window"=>[0,1,2,3],"warm_order"=>[0,1,2,3],
        "outside_prefix"=>Int[],"outside_suffix"=>Int[],
        "objective_denominator"=>"1",
        "pairs"=>["0,1,9,3","2,3,9,5"],"integer_weights"=>["9,7"])
    data=ArchiveData(raw);factors=archive_factors(data)
    @test maximum(length(f.scope) for f in factors)==4
    for order in archive_test_orders(data.window)
        @test evaluate_factors(factors,order)==archived_integer_score(data,order)
    end
    enumeration_result=enumeration(data)
    @test enumeration_result.cost==7*2^2
    @test enumeration_result.evaluated==24 && enumeration_result.coverage_complete
    dp=solve_archive_dp(data)
    @test dp.complete && dp.cost==enumeration_result.cost
end

@testset "Archived DFS folds only destinations unaffected by the window" begin
    baseline=ArchiveData(archive_fixture());old=enumeration(baseline)
    raw=archive_fixture()
    push!(raw["pairs"],"0,4,9,19")
    push!(raw["integer_weights"],"9,11")
    data=ArchiveData(raw);result=enumeration(data)
    constant=big(11)*19^2
    @test result.active_destination_count==2
    @test result.constant_destination_count==1
    @test result.original_destination_count==3
    @test result.outside_constant==constant
    @test result.cost==old.cost+constant
    @test result.minimizers==old.minimizers && result.order==old.order
    @test result.evaluated==old.evaluated && result.coverage_complete
    @test result.validation.score_matches
    for order in archive_test_orders(data.window)
        @test archived_integer_score(data,order)==archived_integer_score(baseline,order)+constant
    end
    # With no internal pair at all, every order has the same exact constant.
    inactive=archive_fixture()
    inactive["pairs"]=["0,4,9,19","0,1,10,-2"]
    inactive["integer_weights"]=["9,11","10,3"]
    all_constant=ArchiveData(inactive);fixed=enumeration(all_constant)
    @test fixed.complete && fixed.active_destination_count==0
    @test fixed.constant_destination_count==2
    @test fixed.cost==11*19^2+3*2^2
    @test fixed.outside_constant==fixed.cost
    @test fixed.minimizers==6 && fixed.evaluated==6
    @test fixed.order==sort(all_constant.window)
    @test fixed.coverage_complete && fixed.validation.score_matches
    @test fixed.validation_CPU_seconds>=0 && fixed.solve_CPU_seconds>=0
    dp=solve_archive_dp(data)
    @test dp.complete && dp.cost==result.cost
    @test dp.factor_construction_performed
    @test dp.factor_construction_CPU_seconds>=0 && dp.factor_construction_wall_seconds>=0
    @test dp.graph_decomposition_CPU_seconds>=0 && dp.graph_decomposition_wall_seconds>=0
    @test dp.validation_CPU_seconds>=0 && dp.validation_wall_seconds>=0
    prepared=archive_factors(data);cached=solve_archive_dp(data;factors=prepared)
    @test !cached.factor_construction_performed && cached.cost==dp.cost
end

@testset "TOML archive schema, weight scaling, ties, and preallocation refusal" begin
    raw=archive_fixture();plain=ArchiveData(raw);baseline=enumeration(plain)
    raw["expected_minimum"]=string(baseline.cost)
    raw["expected_minimizers"]=baseline.minimizers
    raw["expected_window"]=baseline.order
    raw["domain_sha256"]="retained-source-domain-hash"
    mktemp() do path,io
        TOML.print(io,raw);close(io)
        scaled=load_archive(path;weight_multiplier=7)
        @test scaled.expected_minimum==7*baseline.cost
        @test scaled.objective_denominator==25
        @test scaled.metadata["domain_sha256"]==raw["domain_sha256"]
        enumerated=enumeration(scaled)
        @test enumerated.cost==7*baseline.cost
        @test enumerated.validation.expected_matches
        @test enumerated.expected_minimizers_match
        @test enumerated.order==baseline.order
        dp=solve_archive_dp(scaled)
        @test dp.complete && dp.cost==7*baseline.cost
    end
    candidates=enumeration(plain;max_candidates=5)
    @test !candidates.complete && candidates.status==:candidate_limit
    memory=enumeration(plain;max_bytes=1)
    @test !memory.complete && memory.status==:memory_limit
    refused=solve_archive_dp(plain;max_bag_states=1)
    @test !refused.complete && refused.status==:bag_state_limit
    @test refused.graph_edges>=1 && refused.bag_count>=1
    @test refused.heuristic_width_upper_bound>=1
    @test_throws ArgumentError ArchiveData(raw;weight_multiplier=0)
    invalid=archive_fixture();invalid["warm_order"]=[0,2,1,3,4]
    @test_throws ArgumentError ArchiveData(invalid)
end
