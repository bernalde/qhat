function tiny_orders(labels)
    isempty(labels) && return [Int[]]
    [[x;tail] for x in labels for tail in tiny_orders([y for y in labels if y!=x])]
end
function bound_fixture()
    a=VectorFactor([7,91],[[7,91],[91,7]],
        [Dict{Int,BigInt}(),Dict(1=>big(-5),2=>big(2),3=>big(1))])
    orders=permutations([7,20,91])
    vs=[Dict(j=>v for (j,v) in (1=>big(i-1),2=>big(2-i),4=>big(-i)) if !iszero(v)) for i in 1:6]
    c=VectorFactor([7,20,91],orders,vs)
    IntegerResponse([91,7,20],4,big(8),Dict(1=>big(8),2=>big(-3),3=>big(11),4=>big(2)),[a,c],"fixed_window",Dict{String,Any}("fixture"=>"outside complex interference"))
end
@testset "Exact intervals, rational separators and independent branch arithmetic" begin
    r=bound_fixture();minimumQ=minimum(exact_value(r,p) for p in tiny_orders(r.labels))
    for comparisons in (Tuple{Int,Int}[],[(7,91)],[(7,91),(91,20),(20,7)])
        cert=bound_response(r;comparisons,iterations=16)
        L=verify_bound(r,cert)
        @test L==independent_bound(r,cert).L_Q
        @test 0<=L<=minimumQ
        @test all(p->L<=exact_value(r,p),tiny_orders(r.labels))
    end
    cert=bound_response(r;comparisons=[(7,91)],iterations=0)
    bad=deepcopy(cert);pop!(bad["branches"])
    @test_throws ArgumentError verify_bound(r,bad)
    @test_throws ErrorException independent_bound(r,bad)
    bad=deepcopy(cert);bad["response_sha256"]="changed"
    @test_throws ArgumentError verify_bound(r,bad)
    bad=deepcopy(cert);bad["branches"][1]["witnesses"]=[zeros(Int,4)]
    @test_throws ArgumentError verify_bound(r,bad)
    @test_throws ErrorException independent_bound(r,bad)
    cyclic=bound_response(r;comparisons=[(7,91),(91,20),(20,7)],iterations=0,max_witnesses=0)
    index=findfirst(b->b["cyclic"],cyclic["branches"])
    @test !isnothing(index)
    bad=deepcopy(cyclic);bad["branches"][index]["cyclic"]=false
    @test_throws ArgumentError verify_bound(r,bad)
    @test_throws ErrorException independent_bound(r,bad)
    @test_throws ArgumentError bound_response(r;comparisons=[(7,91),(91,7)])
    @test_throws Refusal bound_response(r;caps=Caps(max_work=1))
    @test_throws Refusal bound_response(r;caps=Caps(max_bytes=1))
end
@testset "Zero convex relaxation is retained; complete conditioning can separate" begin
    f=VectorFactor([1,2],[[1,2],[2,1]],[Dict{Int,BigInt}(),Dict(1=>big(-3))])
    r=IntegerResponse([1,2],1,big(6),Dict(1=>big(2)),[f],"fixed_window",Dict{String,Any}())
    raw=bound_response(r)
    @test verify_bound(r,raw)==0
    @test independent_bound(r,raw).L_Q==0
    split=bound_response(r;comparisons=[(1,2)])
    @test verify_bound(r,split)==big(1)//big(36)
    @test independent_bound(r,split).L_Q==big(1)//big(36)
    @test sort([rationalparse(b["L_Q"]) for b in split["branches"]])==[big(1)//36,big(1)//9]
    empty=IntegerResponse([1,2],0,big(1),Dict{Int,BigInt}(),VectorFactor[],"fixed_window",Dict{String,Any}())
    @test verify_bound(empty,bound_response(empty))==0
    constant=IntegerResponse([1,2],2,big(4),Dict(1=>big(-3),2=>big(2)),VectorFactor[],"fixed_window",Dict{String,Any}())
    @test verify_bound(constant,bound_response(constant))==big(13)//16
end

module FrozenFixtures
include(joinpath(@__DIR__,"..","..","julia_cubic_feasibility","src","ResponseCore.jl"))
include(joinpath(@__DIR__,"..","..","julia_cubic_feasibility","src","WindowSearch.jl"))
end
@testset "Frozen adapter realizes six patterns, external response and local scope" begin
    table=zeros(ComplexF64,2,8)
    valid=(1,2,4,5,7,8)
    for (i,col) in enumerate(valid);table[:,col]=[i/8-im*i/16,-i/4+im/2];end
    table[:,3].=NaN+Inf*im;table[:,6].=-Inf+NaN*im
    f=FrozenFixtures.WindowSearch.FlatFactor([1,2,3],[1,3],table)
    flat=FrozenFixtures.WindowSearch.FlatWindow([91,7,20],UInt64[2,3,8],ComplexF64[1+im,2-im,-1+0.5im],[f],0.0)
    r=encode_window(flat)
    for localorder in tiny_orders([1,2,3])
        @test exact_value(r,flat.labels[localorder])==FrozenFixtures.WindowSearch.exact_numerator(flat,localorder)
    end
    cert=bound_response(r)
    @test independent_bound(r,cert).L_Q<=minimum(exact_value(r,p) for p in tiny_orders(r.labels))
    renamed=IntegerResponse(r.labels,r.dimension,r.denominator,r.baseline,r.factors,"full_order_projection",r.binding)
    @test_throws ArgumentError verify_bound(renamed,cert)
    @test_throws ErrorException independent_bound(renamed,cert)
    @test_throws ArgumentError window_record(flat,[1,91,7,20,2],2:4,ones(Int,91),100)
    parent=[1;flat.labels;[i for i in 2:93 if !(i in flat.labels)]]
    rec=window_record(flat,parent,2:4,ones(Int,93),1000;iterations=0,max_witnesses=0)
    @test rec["domain"]=="fixed_window" && rec["positive_full_bound"]===false
    @test_throws ArgumentError window_record(flat,parent,1:3,ones(Int,93),1000)
end
