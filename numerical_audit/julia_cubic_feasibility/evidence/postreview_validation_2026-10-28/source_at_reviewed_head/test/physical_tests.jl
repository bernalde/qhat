using Test, LinearAlgebra
if !isdefined(Main,:ResponseCore)
    include(joinpath(@__DIR__,"..","src","ResponseCore.jl"))
end
if !isdefined(Main,:PhysicalCheck)
    include(joinpath(@__DIR__,"..","src","PhysicalCheck.jl"))
end
using .ResponseCore, .PhysicalCheck

@testset "Residual metric resolves near-identical states" begin
    epsilon=1e-8
    metric=PhysicalCheck._projective_error(ComplexF64[1,0],ComplexF64[1,epsilon])
    expected=setprecision(BigFloat,256) do
        e=BigFloat(epsilon)
        Float64(1-inv(sqrt(1+e*e)))
    end
    @test metric.projective_error>0
    @test isapprox(metric.projective_error,expected;rtol=1e-12)
    @test isapprox(metric.infidelity,epsilon^2;rtol=1e-12)
end

# Independent literal Kronecker matrices, little-endian site convention.
function literal_word(n,letters)
    paulis=Dict('I'=>ComplexF64[1 0;0 1],'X'=>ComplexF64[0 1;1 0],
                'Y'=>ComplexF64[0 -im;im 0],'Z'=>ComplexF64[1 0;0 -1])
    matrix=ones(ComplexF64,1,1)
    for site in n-1:-1:0
        matrix=kron(matrix,paulis[get(letters,site,'I')])
    end
    matrix
end

@testset "Full-space Pauli rotations agree with independent dense exp" begin
    state=ComplexF64[1,im,2,-1];state/=norm(state)
    output=similar(state)
    yy=literal_word(2,Dict(0=>'Y',1=>'Y'))
    apply_word_rotation!(output,state,UInt64(3),UInt64(3),-1,0.37)
    @test output ≈ exp(-im*0.37*yy)*state atol=3e-14
    @test_throws ArgumentError apply_word_rotation!(state,state,UInt64(3),UInt64(3),-1,0.37)
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    initial=zeros(ComplexF64,4);initial[2]=1
    simulated=simulate_order([hop],2,1,[1],0.6,2;formula=2)
    h=(literal_word(2,Dict(0=>'X',1=>'X'))+yy)/2
    @test simulated.state ≈ exp(-im*0.6*h)*initial atol=3e-14
    @test simulated.ideal_rotations==2 && simulated.steps==1
    @test state_leakage(simulated.state,2,1).number_leakage<1e-28
    @test state_leakage(simulated.state,2,1).spin_leakage>0.1
end

@testset "Checked sector reference matches full dense Hamiltonian" begin
    # Number-conserving hop 0<->2 includes the Jordan-Wigner Z1 string.
    hop=Unit(UInt64(5),UInt64[2,7],Int8[1,-1],[0.35,0.35])
    diagonal=Unit(UInt64(0),UInt64[1],Int8[1],[0.2])
    units=[hop,diagonal]
    reference=reference_sector(units,3,1)
    exact=reference_state(reference,0.8)
    h=0.35*(literal_word(3,Dict(0=>'X',1=>'Z',2=>'X'))+
            literal_word(3,Dict(0=>'Y',1=>'Z',2=>'Y')))+
      0.2*literal_word(3,Dict(0=>'Z'))
    initial=zeros(ComplexF64,8);initial[2]=1
    @test exact.state ≈ exp(-im*0.8*h)*initial atol=5e-14
    @test reference.outside_coupling_max==0
    @test length(reference.basis)==3
    records=physical_compare(units,3,1,Dict("forward"=>[1,2],"reverse"=>[2,1]);
                             times=[0.8],budget_multipliers=[3],formulas=[1,2])
    @test length(records)==4
    for record in records
        @test 0<=record.projective_error<=1
        @test record.ideal_rotations<=record.rotation_budget
        @test record.norm_drift<1e-12
        @test record.number_leakage<1e-25
        @test record.spin_leakage<1e-25
        @test record.reference_number_leakage==0
    end
end

@testset "Central/boundary fusion and nonconservation refusal" begin
    hop=Unit(UInt64(3),UInt64[0,3],Int8[1,-1],[0.5,0.5])
    diagonal=Unit(UInt64(0),UInt64[1],Int8[1],[0.2])
    units=[hop,diagonal]
    result=simulate_order(units,2,1,[1,2],0.5,15;formula=2)
    @test result.steps==affordable_steps(units,[1,2],15)
    @test result.ideal_rotations==result.steps*(6-2-1)+2
    # Literal dense multiplication of the unfused Strang sequence.
    xx=literal_word(2,Dict(0=>'X',1=>'X'))
    yy=literal_word(2,Dict(0=>'Y',1=>'Y'))
    zz=literal_word(2,Dict(0=>'Z'))
    a=(xx+yy)/2;b=0.2*zz;dt=0.5/result.steps
    dense=zeros(ComplexF64,4);dense[2]=1
    step=exp(-im*dt*a/2)*exp(-im*dt*b)*exp(-im*dt*a/2)
    for _ in 1:result.steps;dense=step*dense;end
    @test result.state ≈ dense atol=5e-14
    first=simulate_order(units,2,1,[1,2],0.5,15;formula=1)
    first_step=exp(-im*(0.5/first.steps)*b)*exp(-im*(0.5/first.steps)*a)
    dense_first=zeros(ComplexF64,4);dense_first[2]=1
    for _ in 1:first.steps;dense_first=first_step*dense_first;end
    @test first.state ≈ dense_first atol=5e-14
    @test first.ideal_rotations==first.steps*3
    bad=Unit(UInt64(1),UInt64[0],Int8[1],[1.0])
    @test_throws ArgumentError reference_sector([bad],2,0)
    @test_throws ArgumentError physical_compare([bad],2,0,Dict("bad"=>[1]))
    leaked=simulate_order([bad],2,0,[1],0.3,1;formula=1)
    @test isapprox(state_leakage(leaked.state,2,0).number_leakage,sin(0.3)^2;atol=2e-14)
end
