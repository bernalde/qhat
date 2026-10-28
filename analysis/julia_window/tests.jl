# Tiny independent identities and interruption checks; no molecular generation.
using Test
const EXACT_ACTION_PROBE = length(ARGS)==2 && ARGS[1]=="--exact-action-probe"
isempty(ARGS) || EXACT_ACTION_PROBE || error("usage: tests.jl [--exact-action-probe CORE_PATH]")
include(EXACT_ACTION_PROBE ? abspath(ARGS[2]) : joinpath(@__DIR__,"ResponseCore.jl"))
using .ResponseCore

function cancellation_unit(delta::Rational{BigInt})
    exact=Rational{BigInt}[big(1)//big(1)+delta,-big(1)//big(1)]
    Unit(UInt64(1),UInt64[0,0],Int8[1,1],Float64.(exact);exact=exact)
end
if EXACT_ACTION_PROBE
    # Root may point this one test set at a preserved pre-repair Core. Its
    # expected old-code failure is a NumericalResolution, never a large job.
    @testset "Exact action repair probe: supplied nonzero cancellation and exact zero" begin
        delta=big(1)//(big(1)<<80)
        @test unitvalue(cancellation_unit(delta),0)==2.0^-80
        @test unitvalue(cancellation_unit(-delta),0)==-2.0^-80
        @test unitvalue(cancellation_unit(big(0)//big(1)),0)===0.0
    end
    exit()
end
using LinearAlgebra
include(joinpath(@__DIR__,"WindowSearch.jl"))
using .WindowSearch
struct InjectedDeadline <: Exception end

function fixture()
    u(f,p,s,c)=Unit(UInt64(f),UInt64.(p),Int8.(s),Float64.(c))
    [u(0,[0],[1],[.07]),u(3,[0,3],[1,-1],[.5,.5]),u(0,[1],[1],[.3]),
     u(6,[0,6],[1,-1],[.2,.2]),u(0,[6],[1],[.13]),u(0,[1,4],[1,1],[.11,-.19])]
end
function visit_orders(f, values)
    if length(values)<=1
        f(copy(values));return
    end
    for i in eachindex(values)
        head=values[i];rest=values[setdiff(eachindex(values),[i])]
        visit_orders(rest) do tail;f([head;tail]);end
    end
end

@testset "Product tables retain full exterior interference and original costs" begin
    units=fixture();parent=[6,4,2,5,3,1];hf=UInt64(2)
    for duration in (.37,1.6), budget in (17,32,80)
        model=build_window_model(units,parent,hf,2:5;time=duration)
        @test model.setup_stats.aggregation=="exact_real_deltas_before_transport_v2"
        @test model.baseline_column==cubic_column(units,parent,hf)
        @test model.setup_stats.exact_cache_accounted_bytes<=64*1024^2
        @test model.setup_stats.cache_combined_entries<=100000
        flat=prepare_window(model)
        baseline_only=score_order(units[parent[2:5]],collect(1:4),hf,80;time=duration)
        @test !isapprox(baseline_only.cubic_response,score_order(units,parent,hf,budget;time=duration).cubic_response;rtol=1e-4)
        visit_orders(copy(parent[2:5])) do local_order
            order=copy(parent);order[2:5]=local_order
            direct=score_order(units,order,hf,budget;time=duration)
            conditioned=window_score(model,local_order,budget)
            @test conditioned.value≈direct.value rtol=1e-10 atol=1e-26
            lengths=length.(getfield.(units,:coeff));r=direct.steps
            @test r*(2sum(lengths)-lengths[first(order)]-lengths[last(order)])+lengths[first(order)]<=budget
            positions=Int[findfirst(==(label),flat.labels) for label in local_order]
            @test Float64(exact_numerator(flat,positions))≈direct.cubic_response rtol=1e-10 atol=1e-26
        end
        candidate=insertion_window(flat;max_sweeps=3)
        order=copy(parent);order[2:5]=candidate.order
        @test first(order)==first(parent) && last(order)==last(parent)
        @test score_order(units,order,hf,budget;time=duration).value<=score_order(units,parent,hf,budget;time=duration).value+1e-26
    end
end

@testset "Independent deadline callback reaches setup/action/transport/flatten/search" begin
    units=fixture();parent=[6,4,2,5,3,1];hf=UInt64(2)
    function interrupt_after(f, limit)
        calls=Ref(0)
        callback=(kind,amount)->begin
            calls[]+=1
            calls[]>=limit && throw(InjectedDeadline())
        end
        @test_throws InjectedDeadline with_budget(f,callback)
        @test calls[]==limit
    end
    interrupt_after(()->center_units(units,hf),2)
    interrupt_after(()->apply_unit(units[2],Dict(hf=>1.0)),2)
    interrupt_after(()->score_order(units,parent,hf,80;time=.37),32)
    interrupt_after(()->diagonal_response(units,Dict(UInt64(1)=>.1),hf;time=1.6),3)
    interrupt_after(()->build_window_model(units,parent,hf,2:5;time=.37),50)
    model=build_window_model(units,parent,hf,2:5;time=.37)
    interrupt_after(()->prepare_window(model),20)
    flat=prepare_window(model)
    interrupt_after(()->insertion_window(flat),25)
    interrupt_after(()->exact_numerator(flat,[1,2,3,4]),4)
    # Callback restoration, not a permanent mutation of the scientific guard.
    @test isfinite(score_order(units,parent,hf,80).value)
end

@testset "Local resolution semantics and preallocation/storage admission" begin
    @test ResponseCore._local_resolved(0.0,1.0,true,"structural zero")===nothing
    @test_throws NumericalResolution ResponseCore._local_resolved(1e-12,1.0,true,"nonzero cancellation")
    @test ResponseCore._local_resolved(1e-12,1.0,false,"guard disabled")===nothing
    units=fixture();parent=[6,4,2,5,3,1]
    @test_throws ArgumentError build_window_model(units,parent,2,1:4)
    @test_throws ArgumentError build_window_model(units,parent,2,2:5;max_factor_entries=1)
    model=build_window_model(units,parent,2,2:5)
    @test_throws ArgumentError prepare_window(model;max_flat_entries=1)
    huge=(big(1)<<16400)//big(1)
    @test_throws ArgumentError ResponseCore._check_rational_bits(huge)
end

@testset "Exact unit actions repair floating cancellation and preserve centered targets" begin
    delta=big(1)//(big(1)<<80)
    positive=cancellation_unit(delta);negative=cancellation_unit(-delta)
    # Both supplied Float64 sums vanish, while their exact source actions do not.
    @test positive.coeff==[1.0,-1.0] && negative.coeff==[1.0,-1.0]
    @test unitvalue(positive,0)==2.0^-80
    @test unitvalue(negative,0)==-2.0^-80
    @test unitvalue(positive,0;guard=false)==0.0 # Existing explicit opt-out.
    # A different rounded input has a nonzero floating residual but exact zero.
    large=big(10)^16
    zeroexact=Rational{BigInt}[large//big(1),-(large+1)//big(1),big(1)//big(1)]
    zerounit=Unit(UInt64(1),UInt64[0,0,0],Int8[1,1,1],Float64.(zeroexact);exact=zeroexact)
    @test unitvalue(zerounit,0;guard=false)==1.0
    @test unitvalue(zerounit,0)===0.0
    tinyexact=copy(zeroexact);tinyexact[2]+=delta
    tinyunit=Unit(UInt64(1),UInt64[0,0,0],Int8[1,1,1],Float64.(tinyexact);exact=tinyexact)
    @test unitvalue(tinyunit,0;guard=false)==1.0
    @test unitvalue(tinyunit,0)==2.0^-80
    original=Unit(UInt64(0),UInt64[1,2],Int8[1,1],[1.0,-1.0];exact=copy(positive.exact))
    centered=only(center_units([original],0))
    @test unitvalue(original,0)==2.0^-80
    @test unitvalue(centered,0)===0.0
    @test unitvalue(centered,3)==-2.0^-79
    @test original.exact==positive.exact && length(original.coeff)==2
    @test length(centered.coeff)==3 # Numerical scalar does not alter original counts.
    @test affordable_steps([original],[1],2)==1
    @test_throws ArgumentError affordable_steps([original],[1],1)
end

@testset "Exact midpoint cells authenticate ties and repair double-round proposals" begin
    q(x)=Rational{BigInt}(x)
    one=1.0;odd=nextfloat(one);even=nextfloat(odd)
    midpoint1=(q(one)+q(odd))/2
    midpoint2=(q(odd)+q(even))/2
    beyond256=big(1)//(big(1)<<400)
    # Expected values are adjacent binary64 literals, never Float64(rational).
    cases=[(midpoint1-beyond256,one),(midpoint1,one),(midpoint1+beyond256,odd),
           (midpoint2-beyond256,odd),(midpoint2,even),(midpoint2+beyond256,even)]
    for (exact,expected) in cases
        @test ResponseCore._rounded_exact_action(exact)===expected
        @test ResponseCore._rounded_exact_action(-exact)===-expected
        @test ResponseCore._rounding_cell_contains(exact,expected)
    end
    @test !ResponseCore._rounding_cell_contains(midpoint1,odd)
    @test !ResponseCore._rounding_cell_contains(midpoint2,odd)
    @test ResponseCore._rounded_exact_action(q(floatmin(Float64)))===floatmin(Float64)
end

@testset "Exact action rounding refuses unsupported range and bounded arithmetic" begin
    smallest=reinterpret(Float64,UInt64(1))
    second=reinterpret(Float64,UInt64(2))
    qsmall=big(1)//(big(1)<<1074)
    @test ResponseCore._rounded_exact_action(qsmall)===smallest
    @test ResponseCore._rounded_exact_action(-qsmall)===-smallest
    @test unitvalue(cancellation_unit(qsmall),0)===smallest
    @test ResponseCore._rounded_exact_action(3qsmall/2)===second # Subnormal tie-to-even.
    @test ResponseCore._rounded_exact_action(qsmall/2+big(1)//(big(1)<<1500))===smallest
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(qsmall/2)
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(qsmall/4)
    @test_throws NumericalResolution unitvalue(cancellation_unit(qsmall/4),0)
    @test_throws NumericalResolution ResponseCore._rounded_exact_action((big(1)<<1024)//big(1))
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(Rational{BigInt}(floatmax(Float64)))
    @test_throws NumericalResolution ResponseCore._rounded_exact_action(-Rational{BigInt}(floatmax(Float64)))
    overcap=big(1)//(big(1)<<16400)
    @test_throws ArgumentError ResponseCore._rounded_exact_action(overcap)
    @test_throws ArgumentError unitvalue(cancellation_unit(overcap),0)
    # Candidate/midpoint components fit individually, but an odd denominator
    # requiring oversized proof cross-products still refuses rather than grows.
    nearcap=(big(1)<<16383)+1
    @test_throws ArgumentError ResponseCore._rounding_compare(nearcap//(nearcap+2),
        ((big(1)<<52)+1)//(big(1)<<52))
    overflowing=Unit(UInt64(1),UInt64[0,0],Int8[1,1],[floatmax(Float64),floatmax(Float64)])
    @test_throws NumericalResolution unitvalue(overflowing,0)
end

@testset "Exact fallback interruption and downstream global guards remain mandatory" begin
    exact=(big(1)//big(1))+(big(1)//(big(1)<<400))
    for limit in (2,4,8,14)
        calls=Ref(0)
        callback=(kind,amount)->begin
            calls[]+=1
            calls[]>=limit && throw(InjectedDeadline())
        end
        @test_throws InjectedDeadline with_budget(()->ResponseCore._rounded_exact_action(exact),callback)
        @test calls[]==limit
    end
    @test unitvalue(cancellation_unit(big(1)//(big(1)<<80)),0)==2.0^-80 # Callback restored.
    @test_throws NumericalResolution ResponseCore._resolved(0.0,1.0,true)
    @test_throws NumericalResolution ResponseCore._resolved(1e-30,1.0,true)
    @test_throws NumericalResolution diagonal_response(Unit[],Dict(UInt64(0)=>1.0),0;
        guard=true,resolution_scale=1.0) # Genuine HF projection is still guarded.
    zpi=Unit(UInt64(0),UInt64[1],Int8[1],[Float64(pi)])
    @test_throws NumericalResolution diagonal_response([zpi],Dict(UInt64(1)=>1.0),0;
        time=1.0,guard=true,resolution_scale=1.0) # Genuine kernel suppression.
    x=Unit(UInt64(1),UInt64[0],Int8[1],[1.0])
    ztiny=Unit(UInt64(0),UInt64[1],Int8[1],[2.0^-80])
    @test_throws NumericalResolution score_order([x,ztiny],[1,2],0,64;time=1.0)
    @test_throws NumericalResolution ResponseCore._score(Dict(UInt64(1)=>ComplexF64(reinterpret(Float64,UInt64(1)))),1)
    @test_throws NumericalResolution ResponseCore._score(Dict(UInt64(1)=>ComplexF64(1e-100)),big(1)<<180)
end

const V2Q=Rational{BigInt}
function v2_u(a::V2Q=big(1)//big(1))
    exact=V2Q[a/2,a/2]
    Unit(UInt64(3),UInt64[0,3],Int8[1,-1],Float64.(exact);exact)
end
function v2_d(b::V2Q=big(1)//big(1))
    exact=V2Q[b/2,-b/2]
    Unit(UInt64(0),UInt64[2,1],Int8[1,1],Float64.(exact);exact)
end
v2_i()=Unit(UInt64(0),UInt64[0],Int8[1],[1.0])
function v2_literal_paulis()
    # Literal Pauli matrices and little-endian tensor placement; no unit action,
    # sparse atom helper or producer Taylor recurrence defines this oracle.
    one=V2Q[1 0;0 1];x=V2Q[0 1;1 0];z=V2Q[1 0;0 -1]
    y=Complex{V2Q}[0 -im;im 0]
    u=real.((kron(x,x)+kron(y,y))/2)
    d=(kron(z,one)-kron(one,z))/2
    u,d,Matrix{V2Q}(I,4,4)
end
v2_commutator(a,b)=a*b-b*a
function v2_literal_atom(matrices,order)
    a,b=matrices[order[1]],matrices[order[2]]
    if length(order)==2
        comm=v2_commutator(a,b)
        return -v2_commutator(a,comm)/24-v2_commutator(b,comm)/12
    end
    c=matrices[order[3]]
    -(v2_commutator(b,v2_commutator(a,c))+v2_commutator(c,v2_commutator(a,b)))/12
end
function v2_dense_cubic(matrices,order;physical=false)
    mats=physical ? [-im*a for a in matrices] : matrices
    T=eltype(first(mats));dim=size(first(mats),1)
    identity=Matrix{T}(I,dim,dim)
    coefficients=[identity,zeros(T,dim,dim),zeros(T,dim,dim),zeros(T,dim,dim)]
    stages=[(i,big(1)//big(2)) for i in order[1:end-1]]
    push!(stages,(last(order),big(1)//big(1)))
    append!(stages,[(i,big(1)//big(2)) for i in reverse(order[1:end-1])])
    for (label,tau) in stages
        # Independent dense convolution of literal exponential polynomials.
        exponential=[(tau^j/factorial(j))*(mats[label]^j) for j in 0:3]
        coefficients=[sum(exponential[p+1]*coefficients[j-p+1] for p in 0:j) for j in 0:3]
    end
    h=sum(mats)
    coefficients[4]-h^3/6
end
v2_column(column,dim=4)=[get(column,UInt64(j),big(0)//big(1)) for j in 0:dim-1]

@testset "V2 exact Pauli atoms: independent matrices, all orders and coincident flips" begin
    u,d,identity=v2_literal_paulis();hf=UInt64(1)
    @test u[:,2]==V2Q[0,0,1,0] && d[:,2]==V2Q[0,1,0,0]
    t=big(1)//(big(1)<<27)
    commuting=[v2_u(1+t),v2_u(1-t)]
    matrices=[(1+t)*u,(1-t)*u]
    @test v2_commutator(matrices...)==zeros(V2Q,4,4)
    for order in ((1,2),(2,1))
        @test isempty(atom_column(commuting,collect(order),hf))
        @test isempty(ResponseCore._atom_column(AmplitudeCache(commuting),order,hf,16))
        @test v2_literal_atom(matrices,order)==zeros(V2Q,4,4)
    end
    units=[v2_u(),v2_d(big(2)//big(5)),v2_d(big(3)//big(7)),v2_i()]
    mats=[u,(big(2)//big(5))*d,(big(3)//big(7))*d,identity]
    for centered in (false,true)
        cache=AmplitudeCache(centered ? center_units(units,hf) : units)
        for a in 1:3,b in a+1:4
            for order in ((a,b),(b,a))
                actual=ResponseCore._atom_column(cache,order,hf,16)
                @test v2_column(actual)==v2_literal_atom(mats,order)[:,2]
            end
        end
        visit_orders([1,2,3]) do order
            actual=ResponseCore._atom_column(cache,Tuple(order),hf,16)
            @test v2_column(actual)==v2_literal_atom(mats,order)[:,2]
        end
    end
    @test ResponseCore._exact_product3(AmplitudeCache(units),1,2,1,hf)==-big(2)//big(5)
    @test_throws ArgumentError atom_column(units,[1,2],hf;cache=AmplitudeCache(units;guard=false))
    @test_throws ArgumentError atom_column(units,[1,2],hf;guard=false,cache=AmplitudeCache(units))
    float_atom=atom_column(units,[1,2],hf;guard=false)
    @test float_atom isa Dict{UInt64,Float64}
end

@testset "V2 exact conditioned DELTA survives rounded columns and retains exterior/HF" begin
    u,d,identity=v2_literal_paulis();hf=UInt64(1);c=(big(1)<<60)//big(1)
    for b in (big(1)//big(1),big(1)//(big(1)<<30))
        units=[v2_d(c),v2_u(),v2_d(b),v2_i()]
        cache=AmplitudeCache(center_units(units,hf))
        rank=[1,2,3,4];outside=[1,4]
        base=ResponseCore._exact_pair_column(cache,(2,3),outside,rank,2,hf,16)
        changed=ResponseCore._exact_pair_column(cache,(3,2),outside,rank,2,hf,16)
        @test get(base,UInt64(2),0)==(b*b-b*c)/3
        @test get(changed,UInt64(2),0)==-(b*b+2b*c)/6
        @test Float64(base[UInt64(2)])==Float64(changed[UInt64(2)])
        delta=copy(changed);ResponseCore._exact_add!(delta,base,16;subtract=true)
        @test delta==Dict(UInt64(1)=>b/2,UInt64(2)=>-b*b/2)
        @test ResponseCore._rounded_exact_action(delta[UInt64(2)])==-Float64(b*b/2)
        # Full exact Taylor oracle in the complete invariant one-electron sector.
        mats=[c*d,u,b*d,identity]
        sector=[mat[[2,3],[2,3]] for mat in mats]
        expected=(v2_dense_cubic(sector,[1,3,2,4])-v2_dense_cubic(sector,[1,2,3,4]))[:,1]
        @test expected==V2Q[b/2,-b*b/2]
        @test V2Q[delta[UInt64(1)],delta[UInt64(2)]]==expected
        # Exterior terms individually nonzero/equal: dropping them loses the
        # conditioned columns even though they cancel in this particular delta.
        @test ResponseCore._atom_column(cache,(1,2,3),hf,16)[UInt64(2)]==-b*c/3
        @test ResponseCore._atom_column(cache,(1,3,2),hf,16)[UInt64(2)]==-b*c/3
        @test get(ResponseCore._atom_column(cache,(2,3),hf,16),UInt64(2),0)!=base[UInt64(2)]
        transported=ResponseCore._transport_exact_delta(cache,delta,hf,1.0,16)
        x=-Float64(b+c);kernel=(sin(x)/x)*cis(x) # Analytical D gap, not helper output.
        @test !haskey(transported,hf)
        @test transported[UInt64(2)]≈(-Float64(b*b/2))*kernel rtol=1e-14 atol=0
    end
    # A suffix exterior changes the DELTA itself; a dropped exterior term
    # therefore cannot pass merely because prefix contributions cancel.
    b=big(2)//big(5);c=big(3)//big(5)
    suffix_units=[v2_i(),v2_u(),v2_d(b),v2_d(c)]
    suffix_cache=AmplitudeCache(center_units(suffix_units,hf))
    suffix_base=ResponseCore._exact_pair_column(suffix_cache,(2,3),[1,4],[1,2,3,4],2,hf,16)
    suffix_changed=ResponseCore._exact_pair_column(suffix_cache,(3,2),[1,4],[1,2,3,4],2,hf,16)
    ResponseCore._exact_add!(suffix_changed,suffix_base,16;subtract=true)
    @test suffix_changed==Dict(UInt64(1)=>b/2,UInt64(2)=>-b*b/2-b*c)
    @test suffix_changed[UInt64(2)]!=-b*b/2
    suffix_mats=[identity,u,b*d,c*d]
    @test v2_column(suffix_changed)==(v2_dense_cubic(suffix_mats,[1,3,2,4])-v2_dense_cubic(suffix_mats,[1,2,3,4]))[:,2]
    # All24 full orders use independent dense polynomial algebra.
    units=[v2_d(big(3)//big(5)),v2_u(),v2_d(big(2)//big(5)),v2_i()]
    matrices=[(big(3)//big(5))*d,u,(big(2)//big(5))*d,identity]
    cache=AmplitudeCache(center_units(units,hf))
    visit_orders([1,2,3,4]) do order
        atom_sum=ResponseCore.ExactColumn()
        for ia in 1:3,ib in ia+1:4
            ResponseCore._exact_add!(atom_sum,ResponseCore._atom_column(cache,Tuple(order[[ia,ib]]),hf,16),16)
        end
        for ia in 1:2,ib in ia+1:3,ic in ib+1:4
            ResponseCore._exact_add!(atom_sum,ResponseCore._atom_column(cache,Tuple(order[[ia,ib,ic]]),hf,16),16)
        end
        @test v2_column(atom_sum)==v2_dense_cubic(matrices,order)[:,2]
    end
    @test v2_dense_cubic(matrices,[1,2,3,4];physical=true)==im*v2_dense_cubic(matrices,[1,2,3,4])
    # Nonmonotone global labels remain original labels at the window boundary.
    model=build_window_model(units,[4,3,2,1],hf,2:3;time=.37)
    @test model.baseline_window==[3,2]
    @test model.baseline_column==cubic_column(units,[4,3,2,1],hf)
    for local_order in ([3,2],[2,3])
        full=[4;local_order;1]
        @test window_score(model,local_order,40).value≈score_order(units,full,hf,40;time=.37).value rtol=1e-10 atol=1e-26
    end
end

@testset "V2 bounded exact cache, allocated capacity, transport and interruptions" begin
    unit=v2_u();hf=UInt64(1)
    cache=AmplitudeCache([unit];max_entries=1)
    ResponseCore._amplitude(cache,1,hf)
    first=cache.accounted_bytes
    ResponseCore._amplitude(cache,1,UInt64(2))
    @test length(cache.amplitudes)==1 && cache.skipped==1 && cache.accounted_bytes==first
    ResponseCore._local_kernel(cache,UInt64(2),hf,1.0)
    @test isempty(cache.kernels) && length(cache.amplitudes)+length(cache.kernels)==1
    @test cache.accounted_bytes<=cache.max_bytes==64*1024^2
    uncached=AmplitudeCache([unit];max_entries=0)
    for _ in 1:2;@test ResponseCore._amplitude(uncached,1,hf)==1;end
    @test isempty(uncached.amplitudes) && uncached.misses==2
    padded=deepcopy(big(1)//big(1))
    before=ResponseCore._cache_entry_bytes(padded)
    # Test-only capacity change preserves value; no production arithmetic uses
    # GMP mutation. Charge actual allocation, not the resulting one-bit value.
    Base.GMP.MPZ.realloc2!(numerator(padded),1_000_000)
    @test padded==1 && ResponseCore._cache_entry_bytes(padded)>before+100000
    padded_unit=Unit(UInt64(0),UInt64[0],Int8[1],[1.0];exact=V2Q[padded])
    basebytes=4096+128
    bytecache=AmplitudeCache([padded_unit];max_bytes=basebytes+2048)
    @test ResponseCore._amplitude(bytecache,1,UInt64(0))==1
    @test isempty(bytecache.amplitudes) && bytecache.skipped==1 && bytecache.accounted_bytes==basebytes
    @test_throws ArgumentError AmplitudeCache([unit];max_entries=100001)
    @test_throws ArgumentError AmplitudeCache([unit];max_bytes=64*1024^2+1)
    @test ResponseCore._qmul((big(1)<<9000)//big(3),big(3)//(big(1)<<9000))==1
    @test_throws ArgumentError ResponseCore._qmul((big(1)<<9000)//big(1),(big(1)<<9000)//big(1))
    @test_throws ArgumentError ResponseCore._qadd((big(1)<<16383)//big(1),(big(1)<<16383)//big(1))
    emptycache=AmplitudeCache(Unit[])
    tiny=big(1)//(big(1)<<1074)
    @test isempty(ResponseCore._transport_exact_delta(emptycache,Dict(hf=>tiny/4),hf,1.0,16))
    @test_throws NumericalResolution ResponseCore._transport_exact_delta(emptycache,Dict(UInt64(2)=>tiny/4),hf,1.0,16)
    @test_throws NumericalResolution ResponseCore._transport_exact_delta(emptycache,Dict(UInt64(2)=>tiny),hf,1e-300,16)
    @test_throws NumericalResolution ResponseCore._transport_exact_delta(AmplitudeCache(Unit[]),Dict(UInt64(2)=>big(1)//big(1)),hf,0.0,16)
    @test_throws NumericalResolution ResponseCore._transport_exact_delta(AmplitudeCache(Unit[]),Dict(UInt64(2)=>big(1)//big(1)),hf,Inf,16)
    for limit in (2,6,15,35)
        calls=Ref(0)
        callback=(kind,amount)->begin
            calls[]+=1;calls[]>=limit && throw(InjectedDeadline())
        end
        @test_throws InjectedDeadline with_budget(()->begin
            context=AmplitudeCache([v2_u(),v2_d()])
            column=ResponseCore._atom_column(context,(1,2),hf,16)
            ResponseCore._exact_add!(column,copy(column),16;subtract=true)
            ResponseCore._transport_exact_delta(context,column,hf,1.0,16)
        end,callback)
        @test calls[]==limit
    end
    counted=Dict(:action=>0,:atom=>0)
    with_budget((kind,amount)->haskey(counted,kind) && (counted[kind]+=amount)) do
        ResponseCore._atom_column(AmplitudeCache([v2_u(),v2_d()]),(1,2),hf,16)
    end
    @test counted[:atom]==1 && counted[:action]>0
end

using TOML, SHA
include(joinpath(@__DIR__,"worker.jl"))
@testset "Explicit 1024-label/8192-word/20-bit reader admission" begin
    # Parsing only: there is no large-domain numerical search in this test.
    mktempdir() do directory
        path=joinpath(directory,"admission.toml")
        data=Dict{String,Any}("version"=>1,"n"=>20,"reference"=>17,"time"=>"37/100","budget"=>40000,
            "order"=>collect(1023:-1:0),"word_lengths"=>fill(8,1024),"parent_score"=>1.0,
            "width"=>2,"max_windows"=>1,"max_sweeps"=>1,"max_entries"=>50000,
            "max_cached_actions"=>50000,"max_factor_entries"=>250000,"max_flat_entries"=>1000000,
            "max_actions"=>2000000,"max_atoms"=>200000,"max_evaluations"=>100,
            "units"=>[Dict("flip"=>0,"parities"=>collect(0:7),"signs"=>fill(1,8),
                           "coefficients"=>fill("1/100",8)) for _ in 1:1024])
        open(path,"w") do io;TOML.print(io,data);end
        admitted=read_request(path,()->nothing)
        @test length(admitted.units)==1024 && sum(length(u.coeff) for u in admitted.units)==8192
        @test admitted.hf==UInt64(17)
        push!(data["units"],first(data["units"]))
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
        pop!(data["units"])
        push!(data["units"][end]["parities"],8)
        push!(data["units"][end]["signs"],1)
        push!(data["units"][end]["coefficients"],"1/100")
        data["word_lengths"][end]=9
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing) # 8193 words
        data["n"]=21
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
        data["n"]=20;data["version"]=true
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws ErrorException read_request(path,()->nothing)
    end
end
