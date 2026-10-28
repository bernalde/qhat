using Test
if !isdefined(@__MODULE__,:ResponseCore)
    include(joinpath(@__DIR__,"..","src","ResponseCore.jl"))
end
if !isdefined(@__MODULE__,:WindowSearch)
    include(joinpath(@__DIR__,"..","src","WindowSearch.jl"))
end
if !isdefined(@__MODULE__,:StructureDP)
    include(joinpath(@__DIR__,"..","src","StructureDP.jl"))
end
if !isdefined(@__MODULE__,:CubicDP)
    include(joinpath(@__DIR__,"..","src","CubicDP.jl"))
end
using .ResponseCore, .WindowSearch, .StructureDP, .CubicDP

function cdp_orders(labels)
    isempty(labels) && return [()]
    [(v,tail...) for v in labels for tail in cdp_orders([x for x in labels if x!=v])]
end

function cdp_factor(scope,destinations,columns)
    table=zeros(ComplexF64,length(destinations),length(scope)==2 ? 2 : 8)
    for (column,values) in columns
        table[:,column]=values
    end
    WindowSearch.FlatFactor(Int[scope...],Int[destinations...],table)
end

function cdp_check_all(flat;maximum_orders=720)
    data=cubic_factors(flat)
    @test data.complete
    orders=cdp_orders(collect(1:length(flat.labels)))
    @test length(orders)<=maximum_orders
    values=Rational{BigInt}[]
    for order in orders
        expected=exact_numerator(flat,order)
        N=evaluate_factors(data.factors,order)
        @test denominator(N)==1
        @test N/data.denominator==expected
        push!(values,expected)
    end
    solved=solve_cubic_dp(flat)
    @test solved.complete
    @test solved.score==minimum(values)
    @test solved.score==solved.numerator//solved.denominator
    @test solved.score==exact_numerator(flat,solved.local_order)
    @test solved.order==flat.labels[solved.local_order]
    @test sort(solved.order)==sort(flat.labels)
    @test solved.stats.witness_exact_score_matches
    data,solved
end

function cdp_cancellation_fixture()
    a=cdp_factor([1,2],[1,2],Dict(1=>[0im,0im],2=>[-1-im,-0.5+0.25im]))
    b=cdp_factor([2,3],[1],Dict(1=>[0im],2=>[0.25-0.5im]))
    WindowSearch.FlatWindow([91,7,20],UInt64[4,9],ComplexF64[1+im,0.5-0.25im],[a,b],0.0)
end

@testset "Exact complex cubic factors retain outside interference and cancellation" begin
    flat=cdp_cancellation_fixture()
    original_baseline=copy(flat.baseline);original_tables=[copy(f.values) for f in flat.factors]
    data,solved=cdp_check_all(flat)
    @test data.vector_denominator==4 && data.denominator==16
    @test evaluate_factors(data.factors,[1,2,3])==37
    @test solved.score==0 # Nontrivial R0 cancels against a nonzero DELTA.
    @test any(value<0 for f in data.factors for value in values(f.values))
    @test any(imag(value)!=0 for f in flat.factors for value in f.values)
    @test isequal(flat.baseline,original_baseline)
    @test all(isequal(flat.factors[i].values,original_tables[i]) for i in eachindex(flat.factors))
end

@testset "Genuine six-label factor and all six valid triple patterns" begin
    valid=(1,2,4,5,7,8)
    left=Dict{Int,Vector{ComplexF64}}();right=Dict{Int,Vector{ComplexF64}}()
    for (rank,column) in enumerate(valid)
        j=rank-1
        left[column]=ComplexF64[j/4+im*j/8]
        right[column]=ComplexF64[j/2-im*j/4]
    end
    a=cdp_factor([1,2,3],[1],left);b=cdp_factor([4,5,6],[1],right)
    # Unreachable comparison codes are not actual order variants.
    a.values[1,3]=17+9im;a.values[1,6]=-13+7im
    b.values[1,3]=11-8im;b.values[1,6]=5+3im
    flat=WindowSearch.FlatWindow([60,10,40,20,50,30],UInt64[13],ComplexF64[0.125+0.5im],[a,b],0.0)
    data,solved=cdp_check_all(flat)
    @test maximum(length(f.scope) for f in data.factors)==6
    @test data.stats.valid_triple_variants==6
    @test solved.stats.used_width<=5
    @test data.stats.gram_dot_products>0
end

@testset "Actual conditioned cubic window uses the identical frozen objective" begin
    hop(mask,c)=Unit(UInt64(mask),UInt64[0,mask],Int8[1,-1],[c/2,c/2])
    diagonal(mask,c)=Unit(UInt64(0),UInt64[mask],Int8[1],[c])
    units=[diagonal(1,0.4),hop(3,0.7),diagonal(5,-0.25),hop(12,0.3),hop(6,-0.2)]
    model=build_window_model(units,collect(1:5),UInt64(5),2:4;guard=false)
    flat=prepare_window(model)
    data,solved=cdp_check_all(flat)
    @test flat.labels==[2,3,4]
    @test solved.score==exact_numerator(flat,solved.local_order)
    @test data.stats.max_retained_scope<=3
    # Affordability belongs to the complete ORIGINAL order, not scalar tables.
    whole=copy(model.baseline_order);whole[model.window_positions]=solved.order
    r=affordable_steps(units,whole,30sum(length(u.coeff) for u in units))
    @test solved.score/r^4==exact_numerator(flat,solved.local_order)/r^4
end

@testset "All-zero response, exact pruning, and conservative admission" begin
    emptyflat=WindowSearch.FlatWindow([19,4,15],UInt64[],ComplexF64[],WindowSearch.FlatFactor[],0.0)
    data,solved=cdp_check_all(emptyflat)
    @test isempty(data.factors) && data.denominator==1
    @test solved.score==0 && solved.numerator==0
    @test data.stats.retained_factors==0
    flat=cdp_cancellation_fixture()
    @test !cubic_factors(flat;max_factor_entries=1).complete
    @test cubic_factors(flat;max_factor_entries=1).status==:factor_entry_limit
    @test cubic_factors(flat;max_exact_vector_entries=1).status==:vector_entry_limit
    @test cubic_factors(flat;max_bits=2).status==:bit_limit
    @test cubic_factors(flat;max_work=1).status==:work_limit
    @test cubic_factors(flat;max_bytes=1).status==:memory_limit
    @test_throws ErrorException cubic_factors(flat;guard=()->error("declared deadline"))
    refused=solve_cubic_dp(flat;max_bag_states=1)
    @test !refused.complete && refused.status==:bag_state_limit
    @test isempty(refused.order) && refused.score===nothing
end

@testset "Constant response factors and monotone bit-ceiling admission" begin
    factor=cdp_factor([1,2],[1],Dict(1=>[1+0im],2=>[1+0im]))
    flat=WindowSearch.FlatWindow([19,4],UInt64[0],ComplexF64[1+im],[factor],0.0)
    data,solved=cdp_check_all(flat)
    @test all(isempty(f.scope) for f in data.factors)
    @test data.stats.constant_factor_value==5 && solved.score==5
    fixture=cdp_cancellation_fixture()
    ordinary=cubic_factors(fixture;max_bits=4096,max_bytes=100_000)
    relaxed=cubic_factors(fixture;max_bits=1_000_000,max_bytes=100_000)
    @test ordinary.complete && relaxed.complete
    @test ordinary.stats.estimated_peak_bytes==relaxed.stats.estimated_peak_bytes
    malformed=WindowSearch.FlatWindow([1,2],UInt64[0,1],ComplexF64[],WindowSearch.FlatFactor[],0.0)
    @test_throws ArgumentError cubic_factors(malformed;max_exact_vector_entries=1)
    emptyfactor=WindowSearch.FlatFactor([1,2],Int[],zeros(ComplexF64,0,2))
    metadata=WindowSearch.FlatWindow([1,2],UInt64[],ComplexF64[],fill(emptyfactor,100),0.0)
    refused=cubic_factors(metadata;max_bytes=1024)
    @test !refused.complete && refused.status==:memory_limit
end
