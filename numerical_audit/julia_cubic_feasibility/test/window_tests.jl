using Test
if !isdefined(Main,:ResponseCore); include(joinpath(@__DIR__,"..","src","ResponseCore.jl")); end
if !isdefined(Main,:WindowSearch); include(joinpath(@__DIR__,"..","src","WindowSearch.jl")); end
using .ResponseCore, .WindowSearch

@testset "Exact frozen-table window certificate and enumeration" begin
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    z(a)=Unit(UInt64(0),UInt64[1],Int8[1],[a])
    units=[z(0.3),hop,z(0.5),z(0.2),hop]
    model=build_window_model(units,collect(1:5),UInt64(2),2:4;guard=false)
    flat=prepare_window(model)
    result=enumerate_window(flat)
    @test result.permutations==6
    @test sort(result.order)==[2,3,4]
    vals=[exact_numerator(flat,p) for p in ([1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1])]
    @test result.numerator==minimum(vals)
    @test result.minimizers==count(==(minimum(vals)),vals)
    insertion=insertion_window(flat)
    @test insertion.numerator>=result.numerator
    @test_throws ErrorException enumerate_window(flat;guard=()->error("deadline"))
end

@testset "Independent Heap coverage and global-label tie break" begin
    for k in 2:5
        seen=Set{Tuple{Vararg{Int}}}(); visits=Ref(0)
        WindowSearch.visit_orders(collect(1:k)) do order
            visits[]+=1;push!(seen,Tuple(order))
        end
        @test visits[]==factorial(k)
        @test length(seen)==factorial(k)
    end
    flat=WindowSearch.FlatWindow([20,10],UInt64[1],ComplexF64[1],WindowSearch.FlatFactor[],1e-14)
    result=enumerate_window(flat)
    @test result.order==[10,20]
    @test result.minimizers==2
end

@testset "Window certificates validate supplied representations and refresh bounds" begin
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    z(a)=Unit(UInt64(0),UInt64[1],Int8[1],[a])
    model=build_window_model([z(0.3),hop,z(0.5),z(0.2),hop],collect(1:5),UInt64(2),2:4;guard=false)
    flat=prepare_window(model);original=enumerate_window(flat)
    understated=WindowSearch.FlatWindow(flat.labels,flat.destinations,flat.baseline,flat.factors,0.0)
    checked=enumerate_window(understated)
    @test checked.numerator==original.numerator && checked.minimizers>=1
    @test checked.roundoff_bound>0
    changed=deepcopy(flat);changed.baseline .*= 2
    for factor in changed.factors;factor.values .*= 2;end
    rescored=enumerate_window(changed)
    @test rescored.numerator==4*original.numerator
    @test rescored.roundoff_bound>changed.error_bound
    for invalid in (-1.0,NaN,Inf)
        bad=WindowSearch.FlatWindow(flat.labels,flat.destinations,flat.baseline,flat.factors,invalid)
        @test_throws ArgumentError enumerate_window(bad)
    end
    for labels in (Int[],collect(1:10))
        bad=WindowSearch.FlatWindow(labels,UInt64[],ComplexF64[],WindowSearch.FlatFactor[],0.0)
        @test_throws ArgumentError enumerate_window(bad)
    end
    outside=WindowSearch.FlatFactor([1,2],[2],zeros(ComplexF64,1,2))
    bad=WindowSearch.FlatWindow([1,2],UInt64[0],ComplexF64[1],[outside],0.0)
    @test_throws ArgumentError enumerate_window(bad)
    duplicate=deepcopy(model);duplicate.baseline_window[2]=duplicate.baseline_window[1]
    @test_throws ArgumentError prepare_window(duplicate)
    repeated=deepcopy(model);repeated.factors[1].scope[2]=repeated.factors[1].scope[1]
    @test_throws ArgumentError prepare_window(repeated)
    incomplete=deepcopy(model);delete!(incomplete.factors[1].variants,first(keys(incomplete.factors[1].variants)))
    @test_throws ArgumentError prepare_window(incomplete)
    malformed=deepcopy(model);factor=malformed.factors[1];key=first(keys(factor.variants))
    vector=pop!(factor.variants,key);factor.variants[(first(key),first(key))]=vector
    @test_throws ArgumentError prepare_window(malformed)
end

if !isdefined(Main,:StructureDP);include(joinpath(@__DIR__,"..","src","StructureDP.jl"));end
if !isdefined(Main,:Diagnostics);include(joinpath(@__DIR__,"..","src","Diagnostics.jl"));end
@testset "Diagnostics validate parent orders and retain canonical partial traversal" begin
    z(a)=Unit(UInt64(0),UInt64[1],Int8[1],[a])
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    units=[z(0.3),hop,z(0.5),z(0.2)]
    for invalid in ([1,1,3,4],[1,2,3],[0,2,3,4],[1,2,3,5])
        @test_throws ArgumentError Diagnostics.select_window(units,invalid,UInt64(2),2)
        @test_throws ArgumentError Diagnostics.structural_audit(units,invalid,UInt64(2))
    end
    forward=Diagnostics.structural_audit(units,collect(1:4),UInt64(2);max_scopes=2)
    reverse=Diagnostics.structural_audit(units,collect(4:-1:1),UInt64(2);max_scopes=2)
    @test forward["status"]==reverse["status"]=="scope_or_atom_cap"
    @test forward["scope_count"]==reverse["scope_count"]==2
    @test forward["potential_edges_retained"]==reverse["potential_edges_retained"]
    @test forward["destinations"]==reverse["destinations"]
end
