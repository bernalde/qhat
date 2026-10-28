using Test
if !isdefined(@__MODULE__, :StructureDP)
    include(joinpath(@__DIR__, "..", "src", "StructureDP.jl"))
end
using .StructureDP

@testset "Integer denominator admission and rational bit limits" begin
    f=LocalFactor([1,2],Dict((1,2)=>0,(2,1)=>1))
    result=bag_order_dp([1,2],fill(f,20),[[1,2]],Tuple{Int,Int}[];max_cost_bits=16)
    @test result.complete && result.cost==0
    @test result.details.cost_bits<=16
    rational=LocalFactor([1,2],Dict((1,2)=>1//(1<<30),(2,1)=>1//(1<<29)))
    refused=bag_order_dp([1,2],[rational],[[1,2]],Tuple{Int,Int}[];max_cost_bits=16)
    @test !refused.complete && refused.status==:cost_bit_limit
end

function structure_permutations(labels)
    isempty(labels) && return [()]
    out=Tuple[]
    for label in labels
        for rest in structure_permutations([x for x in labels if x!=label])
            push!(out,(label,rest...))
        end
    end
    return out
end

function structure_indicator(scope, preferred, reward)
    return LocalFactor(scope, Dict(order => (order==preferred ? reward : 0)
                                   for order in structure_permutations(scope)))
end

@testset "Bag-order DP: true gluing, signed costs, exact tables" begin
    f1=structure_indicator([1,2,3],(1,2,3),-9)
    f2=structure_indicator([2,3,4],(2,4,3),-8)
    result=bag_order_dp(1:4,[f1,f2],[[1,2,3],[2,3,4]],[(1,2)])
    @test result.complete
    @test result.cost==-17
    @test evaluate_factors([f1,f2],result.order)==-17
    @test Tuple(x for x in result.order if x in (1,2,3))==(1,2,3)
    @test Tuple(x for x in result.order if x in (2,3,4))==(2,4,3)

    cycle=[structure_indicator([1,2],(1,2),-1),
           structure_indicator([2,3],(2,3),-1),
           structure_indicator([1,3],(3,1),-1)]
    @test_throws ArgumentError validate_decomposition(1:3,cycle,
        [[1,2],[2,3],[1,3]],[(1,2),(2,3)])
    exact=bag_order_dp(1:3,cycle,[[1,2,3]],Tuple{Int,Int}[])
    @test exact.complete && exact.cost==-2
    fixed=bag_order_dp(1:3,cycle,[[1,2,3]],Tuple{Int,Int}[];
                       first_label=3,last_label=1)
    @test fixed.order==[3,2,1]
    @test fixed.cost==evaluate_factors(cycle,[3,2,1])

    floating=LocalFactor([1,2],Dict((1,2)=>0.1,(2,1)=>0.2))
    fresult=bag_order_dp(1:2,[floating],[[1,2]],Tuple{Int,Int}[])
    @test fresult.order==[1,2]
    @test fresult.cost==Rational{BigInt}(0.1)
    @test fresult.cost!=big(1)//big(10)
    @test fresult.details.encoded_tables_exact
end

@testset "Disconnected cases and factor-scope reduction" begin
    factors=[structure_indicator([1,2],(1,2),-7),
             structure_indicator([3,4],(4,3),-3),
             LocalFactor(Int[],Dict(()=>-2))]
    decomposition=elimination_decomposition(1:5,[(1,2),(3,4)])
    @test decomposition.complete
    @test validate_decomposition(1:5,factors,decomposition.bags,decomposition.tree_edges)
    result=bag_order_dp(1:5,factors,decomposition.bags,decomposition.tree_edges)
    @test result.complete && result.cost==-12
    @test sort(result.order)==collect(1:5)
    empty=bag_order_dp(Int[],[LocalFactor(Int[],Dict(()=>3//7))],
                       [Int[]],Tuple{Int,Int}[])
    @test empty.complete && isempty(empty.order) && empty.cost==3//7
    singleton=bag_order_dp([9],LocalFactor[],[[9]],Tuple{Int,Int}[];
                           first_label=9,last_label=9)
    @test singleton.complete && singleton.order==[9]

    redundant=LocalFactor([1,2,3],Dict(order=>(findfirst(==(1),order)<findfirst(==(2),order) ? 0 : 1)
                                     for order in structure_permutations([1,2,3])))
    @test Set(redundant.scope)==Set([1,2])
    @test length(redundant.values)==2
    constant=LocalFactor([1,2],Dict((1,2)=>-4,(2,1)=>-4))
    @test isempty(constant.scope)
    @test constant.values[()]==-4
end

@testset "Structural evidence distinguishes potential from exact" begin
    f=structure_indicator([1,2,3],(1,2,3),-1)
    exact=audit_factors(1:4,[f])
    @test exact.complete && exact.graph_kind==:exact_factor_graph
    @test exact.clique_width_lower_bound==2
    @test Set(exact.clique_witness)==Set([1,2,3])
    potential=audit_destination_buckets(1:4,[0=>[[1,2],[2,4]]];max_pair_events=0)
    @test potential.complete
    @test potential.graph_kind==:support_overlap_potential_supergraph
    @test potential.clique_width_lower_bound===nothing
    @test potential.potential_pair_events==1
    @test !potential.full_pair_enumeration_admissible
    @test potential.pair_enumeration_status==:pair_event_limit
    @test potential.edges==Set([(1,2),(1,4),(2,4)])
    refused=audit_destination_buckets(1:4,[0=>[[1,2],[2,4]]];max_atom_visits=1)
    @test !refused.complete && refused.status==:atom_visit_limit
    @test !refused.pair_count_complete
    capped=audit_destination_buckets(1:4,[0=>[[1,2],[2,4]]];max_edge_visits=1)
    @test !capped.complete && capped.status==:edge_work_limit
    @test length(capped.edges)<=1
    @test_throws ArgumentError audit_destination_buckets(1:4,
        [0=>[[1,2]],0=>[[2,4]]])
    @test_throws ArgumentError LocalFactor([1,2],Dict((1,2)=>0))
end

@testset "Preallocation refusals and checked decomposition" begin
    calls=Ref(0);guard=()->(calls[]+=1)
    refused=bag_order_dp(1:100,LocalFactor[],[collect(1:100)],Tuple{Int,Int}[];
                         max_bag_states=10,guard=guard)
    @test !refused.complete && refused.status==:bag_state_limit
    @test calls[]==0 # No permutation/state generation follows refusal.
    f=structure_indicator([1,2],(1,2),-1)
    memory=bag_order_dp(1:2,[f],[[1,2]],Tuple{Int,Int}[];max_bytes=1)
    @test !memory.complete && memory.status==:memory_limit
    work=bag_order_dp(1:2,[f],[[1,2]],Tuple{Int,Int}[];max_work=0)
    @test !work.complete && work.status==:work_limit
    bitlimit=bag_order_dp(1:2,[f],[[1,2]],Tuple{Int,Int}[];max_cost_bits=1)
    @test !bitlimit.complete && bitlimit.status==:cost_bit_limit
    entries=bag_order_dp(1:2,[f],[[1,2]],Tuple{Int,Int}[];max_factor_entries=1)
    @test !entries.complete && entries.status==:factor_entry_limit
    @test_throws ArgumentError bag_order_dp(1:2,[f],[[1],[2]],[(1,2)])
    @test_throws ArgumentError bag_order_dp(1:2,[f],[[1,2]],Tuple{Int,Int}[];
                                           first_label=1,last_label=1)
    limited=elimination_decomposition(1:4,[(1,2),(2,3),(3,4)];max_work=1)
    @test !limited.complete && limited.status==:work_limit
    @test limited.width_upper_bound===nothing
end

@testset "Singleton endpoints retain normal admission and result details" begin
    constant=LocalFactor(Int[],Dict(()=>3//7))
    args=([9],[constant],[[9]],Tuple{Int,Int}[])
    result=bag_order_dp(args...;first_label=9,last_label=9)
    @test result.complete && result.order==[9] && result.cost==3//7
    @test result.details.states==1
    @test result.details.work>0 && result.details.estimated_bytes>0
    @test result.details.cost_bits>0 && !result.details.physical_response_exactness_claim
    for (name,status) in ((:max_bag_states,:bag_state_limit),(:max_states,:state_limit),
                          (:max_work,:work_limit),(:max_bytes,:memory_limit))
        caps=NamedTuple{(name,)}((0,))
        refused=bag_order_dp(args...;first_label=9,last_label=9,caps...)
        @test !refused.complete && refused.status==status
    end
    @test_throws ErrorException bag_order_dp(args...;first_label=9,last_label=9,
                                             guard=()->error("declared cancellation"))
end

@testset "Capped factor graph audit retains only complete prefix factors" begin
    pair=structure_indicator([1,2],(1,2),-1)
    triple=structure_indicator([2,3,4],(2,3,4),-1)
    partial=audit_factors(1:4,[pair,triple];max_edges=1)
    @test !partial.complete && partial.status==:edge_limit
    @test partial.factors_visited==1 && partial.entries_visited==2
    @test partial.edges==Set([(1,2)])
    @test partial.clique_witness==[1,2] && partial.clique_width_lower_bound==1
    empty=audit_factors(1:4,[triple];max_edges=2)
    @test !empty.complete && empty.status==:edge_limit
    @test empty.factors_visited==0 && empty.entries_visited==0
    @test isempty(empty.edges) && isempty(empty.clique_witness)
end
