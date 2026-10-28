using Test, LinearAlgebra, TOML
const STUDY_ROOT=normpath(joinpath(@__DIR__,".."))
const PRIOR_ROOT=normpath(joinpath(STUDY_ROOT,"..","julia_cubic_feasibility"))
include(joinpath(PRIOR_ROOT,"src","ResponseCore.jl"))
include(joinpath(STUDY_ROOT,"src","WindowPhysical.jl"))
include(joinpath(STUDY_ROOT,"src","WindowStudyIO.jl"))
const PhysicalCheck=WindowPhysical
using .ResponseCore, .WindowStudyIO
# Reuse the audited low-register rotation/reference/fusion contracts against
# the new adapter, whose only numerical-domain change is the fourteen-mode cap.
include(joinpath(PRIOR_ROOT,"test","physical_tests.jl"))

@testset "Fourteen-mode high-bit action preserves the full state" begin
    n=14;bit=UInt64(1)<<13;angle=0.31
    input=zeros(ComplexF64,1<<n);input[1]=1
    output=similar(input)
    PhysicalCheck.apply_word_rotation!(output,input,bit,UInt64(0),1,angle)
    expected=zeros(ComplexF64,1<<n)
    expected[1]=cos(angle);expected[Int(bit)+1]=-im*sin(angle)
    @test norm(output-expected)<1e-14
    x=Unit(bit,UInt64[0],Int8[1],[0.7])
    circuit=PhysicalCheck.simulate_order([x],n,0,[1],0.4,1)
    leakage=PhysicalCheck.state_leakage(circuit.state,n,0)
    @test isapprox(leakage.number_leakage,sin(0.28)^2;atol=1e-14)
    @test isapprox(leakage.spin_leakage,sin(0.28)^2;atol=1e-14)
    @test circuit.ideal_rotations==1
    @test_throws ArgumentError PhysicalCheck.simulate_order([x],15,0,[1],0.4,1)
    @test_throws ArgumentError PhysicalCheck.reference_sector([x],14,0)
end

@testset "Fourteen-mode conserving reference and independently unfused circuit" begin
    n=14;high=UInt64(1)<<13;flip=high|UInt64(1)
    middle=((UInt64(1)<<13)-1) & ~UInt64(1)
    hop=Unit(flip,UInt64[middle,middle|flip],Int8[1,-1],[0.35,0.35])
    z0=Unit(UInt64(0),UInt64[1],Int8[1],[0.2])
    z13=Unit(UInt64(0),UInt64[high],Int8[1],[0.3])
    units=[hop,z0,z13];hf=UInt64(1);t=0.4
    reference=PhysicalCheck.reference_sector(units,n,hf)
    exact=PhysicalCheck.reference_state(reference,t)
    omega=sqrt(0.1^2+0.7^2)
    analytic=zeros(ComplexF64,1<<n)
    analytic[2]=cos(omega*t)-im*sin(omega*t)*0.1/omega
    analytic[Int(high)+1]=-im*sin(omega*t)*0.7/omega
    @test length(reference.basis)==14
    @test reference.outside_coupling_max==0
    @test norm(exact.state-analytic)<2e-13
    @test PhysicalCheck.state_leakage(exact.state,n,hf).number_leakage<1e-28
    # In the occupied two-state subspace, the three ORIGINAL units have these
    # literal matrices. Multiplying all unfused half/full exponentials is an
    # independent oracle for the optimized full-register fusion loop.
    matrices=[ComplexF64[0 0.7;0.7 0],ComplexF64[-0.2 0;0 0.2],
              ComplexF64[0.3 0;0 -0.3]]
    circuit=PhysicalCheck.simulate_order(units,n,hf,[1,2,3],t,120)
    vector=ComplexF64[1,0]
    for _ in 1:circuit.steps, (index,weight) in ((1,0.5),(2,0.5),(3,1.0),(2,0.5),(1,0.5))
        vector=exp(-im*(t/circuit.steps)*weight*matrices[index])*vector
    end
    expected=zeros(ComplexF64,1<<n);expected[2]=vector[1];expected[Int(high)+1]=vector[2]
    @test norm(circuit.state-expected)<2e-13
    @test circuit.steps==23
    @test circuit.ideal_rotations==117
    @test circuit.ideal_rotations<=120
    @test PhysicalCheck.state_leakage(circuit.state,n,hf).number_leakage<1e-28
    @test_throws ArgumentError PhysicalCheck.reference_sector(units,n,(UInt64(1)<<7)-1)
    @test_throws ArgumentError PhysicalCheck.reference_sector(units,n,15;max_sector_dim=1000)
    @test_throws ArgumentError PhysicalCheck.reference_sector(units,n,hf;max_sector_dim=2001)
end

@testset "Study adapter admits declared partitions and rejects excess input before work" begin
    function write_model(path;count=1,partition="conserving",denominator="1",word="13Z",numerator="1")
        data=Dict("case"=>"mechanics_only","n"=>14,"electrons"=>4,
            "denominator"=>denominator,"partition"=>partition,
            "baseline_order"=>collect(1:count),
            "units"=>[Dict("words"=>[word],"numerators"=>[numerator]) for _ in 1:count])
        open(path,"w") do io;TOML.print(io,data);end
    end
    mktempdir() do directory
        path=joinpath(directory,"model.toml")
        write_model(path)
        data=load_model(path)
        @test data.n==14 && data.electrons==4 && data.words==1
        write_model(path;count=1025)
        @test_throws ErrorException load_model(path)
        write_model(path;count=1025,partition="pauli")
        @test length(load_model(path).units)==1025
        write_model(path;count=8193,partition="pauli")
        @test_throws ErrorException load_model(path)
        write_model(path;denominator="1"^1281)
        @test_throws ErrorException load_model(path)
        write_model(path;numerator="1"^1281)
        @test_throws ErrorException load_model(path)
        write_model(path;word="0Z,"^100)
        @test_throws ErrorException load_model(path)
    end
end
