using Test
if !isdefined(Main, :ResponseCore)
    include(joinpath(@__DIR__, "..", "src", "ResponseCore.jl"))
end
using .ResponseCore

@testset "Conditioned construction rejects oversized requests before work" begin
    z=Unit(UInt64(0),UInt64[1],Int8[1],[1.0])
    @test_throws ArgumentError build_window_model(fill(z,10),collect(1:10),UInt64(0),1:10)
    @test_throws ArgumentError build_window_model([z,z],[1,2],UInt64(0),1:2;max_factor_entries=0)
end

oneword(flip, parity, coefficient=1.0; phase=1) =
    Unit(UInt64(flip), UInt64[parity], Int8[phase], Float64[coefficient])

function close_columns(a,b; atol=2e-12)
    all(isapprox(get(a,k,0), get(b,k,0); atol=atol, rtol=2e-11) for k in union(keys(a),keys(b)))
end

function tiny_permutations(labels)
    length(labels) <= 1 && return [Int[labels...]]
    result = Vector{Int}[]
    for i in eachindex(labels)
        rest = [labels[j] for j in eachindex(labels) if i != j]
        for tail in tiny_permutations(rest)
            push!(result, [labels[i]; tail])
        end
    end
    result
end

@testset "Real equal-flip action and canonical conserving hop" begin
    hop = Unit(UInt64(3), UInt64[0,3], Int8[1,-1], Float64[0.5,0.5])
    @test unitvalue(hop,UInt64(2)) == 1.0
    @test unitvalue(hop,UInt64(1)) == 1.0
    @test unitvalue(hop,UInt64(0)) == 0.0
    @test apply_unit(hop,Dict(UInt64(2)=>1.0)) == Dict(UInt64(1)=>1.0)
    @test_throws ArgumentError Unit(UInt64(0),UInt64[0],Int8[0],Float64[1.0])
end

@testset "Independent cubic Taylor coefficients and transport" begin
    x,z = oneword(1,0),oneword(0,1)
    @test close_columns(cubic_column([x,z],[1,2],0;guard=false), Dict(UInt64(1)=>1/3,UInt64(0)=>-1/6))
    @test close_columns(cubic_column([x,z],[2,1],0;guard=false), Dict(UInt64(1)=>-1/6,UInt64(0)=>1/3))
    response = diagonal_response([x,z],cubic_column([x,z],[1,2],0;guard=false),0)
    @test !haskey(response,UInt64(0))
    @test isapprox(abs(response[UInt64(1)]),sin(1)/3;atol=2e-14)
    @test affordable_steps([x,z],[1,2],3) == 1
    @test_throws ArgumentError affordable_steps([x,z],[1,2],2)
    @test score_order([z],[1],0,1).value == 0.0
end

@testset "Pair and distinct-triple atoms reconstruct independent Taylor" begin
    units = [oneword(1,0,0.7),oneword(0,1,-0.4),oneword(2,0,0.3),oneword(0,3,0.2)]
    for order in tiny_permutations(collect(1:4))
        reconstructed = Dict{UInt64,Float64}()
        for i in 1:3, j in i+1:4
            for (k,v) in atom_column(units,(order[i],order[j]),0)
                reconstructed[k] = get(reconstructed,k,0.0)+v
            end
        end
        for i in 1:2, j in i+1:3, k in j+1:4
            for (s,v) in atom_column(units,(order[i],order[j],order[k]),0)
                reconstructed[s] = get(reconstructed,s,0.0)+v
            end
        end
        @test close_columns(reconstructed,cubic_column(units,order,0;guard=false))
    end
end

@testset "Conditioned contiguous window keeps outside interference" begin
    # The fixed prefix/suffix are noncommuting and contribute to the same
    # destinations as variable atoms. Removing the baseline would change Q.
    hop(mask,coefficient) = Unit(UInt64(mask),UInt64[0,mask],Int8[1,-1],
                                Float64[coefficient/2,coefficient/2])
    units = [oneword(0,1,0.6),hop(3,0.8),oneword(0,5,-0.3),
             hop(12,0.4),hop(6,-0.2)]
    base = [1,2,3,4,5]
    hf = UInt64(5) # Two electrons, one initially on each disjoint rail.
    model = build_window_model(units,base,hf,2:4;guard=false)
    @test model.baseline_window == [2,3,4]
    @test length(model.factors) == 4 # three conditioned pairs and one triple.
    @test window_response(model,model.baseline_window;guard=false) == model.baseline_response
    @test any(!isempty(v) for f in model.factors for v in values(f.variants))
    for local_order in tiny_permutations([2,3,4])
        whole = [1;local_order;5]
        direct = diagonal_response(units,cubic_column(units,whole,hf;guard=false),hf)
        @test close_columns(window_response(model,local_order;guard=false),direct)
        @test isapprox(window_score(model,local_order,30;guard=false).value,
                       score_order(units,whole,hf,30;guard=false).value;atol=2e-12,rtol=2e-11)
    end
    @test_throws ArgumentError window_response(model,[2,3,5];guard=false)
end

@testset "Scalar centering preserves cubic identity and ORIGINAL rotation accounting" begin
    x = oneword(1,0,0.7)
    z = Unit(UInt64(0),UInt64[1,0],Int8[1,1],Float64[0.4,3.0])
    units = [x,z,oneword(0,3,-0.2)]
    @test close_columns(cubic_column(units,[1,2,3],0;guard=false,center=false),
                        cubic_column(units,[1,2,3],0;guard=false,center=true))
    model = build_window_model(units,[1,2,3],0,1:2;guard=false)
    @test length.(getfield.(model.units,:coeff)) == [1,2,1]
    @test window_score(model,[1,2],15;guard=false).steps == affordable_steps(units,[1,2,3],15)
    @test window_score(model,[2,1],15;guard=false).steps == affordable_steps(units,[2,1,3],15)
end

@testset "Near-zero response is refused, not reported as a completed zero" begin
    # Exact selective projected cancellation: (X,Z,2X). The raw cubic column
    # has a parallel component, so the transported/projection guard catches it.
    units = [oneword(1,0,1.0),oneword(0,1,1.0),oneword(1,0,2.0)]
    @test_throws NumericalResolution score_order(units,[1,2,3],0,20;guard=true)
    @test score_order(units,[1,2,3],0,20;guard=false).cubic_response <= 1e-25
end
