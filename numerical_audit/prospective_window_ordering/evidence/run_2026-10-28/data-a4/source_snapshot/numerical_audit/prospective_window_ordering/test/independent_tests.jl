using Test, LinearAlgebra, TOML, SHA
include(joinpath(@__DIR__,"..","independent.jl"))
using .IndependentPhysical
const IP=IndependentPhysical

function write_input(path; n=4,electrons=2,units=nothing,case="independent_mechanics")
    blocks=isnothing(units) ? [
        Dict("words"=>["0X,1Z,2X","0Y,1Z,2Y"],"numerators"=>["35","35"]),
        Dict("words"=>["0Z"],"numerators"=>["20"]),
        Dict("words"=>["2Z","1Z"],"numerators"=>["30","-11"]),
        Dict("words"=>[""],"numerators"=>["7"])] : units
    data=Dict("case"=>case,"n"=>n,"electrons"=>electrons,"partition"=>"conserving",
        "source_json_sha256"=>repeat("0",64),"denominator"=>"100",
        "baseline_order"=>collect(1:length(blocks)),"units"=>blocks)
    open(path,"w") do io;TOML.print(io,data);end
    data
end
function dense_word(key,n)
    identity=Matrix{ComplexF64}(I,2,2)
    matrices=Dict('X'=>ComplexF64[0 1;1 0], 'Y'=>ComplexF64[0 -im;im 0],
        'Z'=>ComplexF64[1 0;0 -1])
    labels=Dict(key);answer=ones(ComplexF64,1,1)
    for site in n-1:-1:0
        answer=kron(answer,haskey(labels,site) ? matrices[labels[site]] : identity)
    end
    answer
end
function dense_hamiltonian(model)
    matrix=zeros(ComplexF64,1<<model.n,1<<model.n)
    for word in model.words;matrix.+=word.coefficient*dense_word(word.key,model.n);end
    matrix
end
function dense_unfused(model,partition,order,t,steps;formula=2)
    units=IP.partition_units(model,partition)
    matrices=[sum((model.words[i].coefficient*dense_word(model.words[i].key,model.n) for i in ids);
                  init=zeros(ComplexF64,1<<model.n,1<<model.n)) for ids in units]
    state=zeros(ComplexF64,1<<model.n);state[Int(model.hf)+1]=1
    for _ in 1:steps
        for direction in (formula==1 ? (order,) : (order,reverse(order)))
            for i in direction
                state=exp(-im*(t/steps)*(formula==1 ? 1. : .5)*matrices[i])*state
            end
        end
    end
    state
end
function write_state(path,state)
    open(path,"w") do io
        for amplitude in state
            write(io,htol(reinterpret(UInt64,real(amplitude))))
            write(io,htol(reinterpret(UInt64,imag(amplitude))))
        end
    end
    bytes2hex(sha256(read(path)))
end
function oracle_metrics(reference,state,n,hf)
    nr=sum(abs2,reference);ns=sum(abs2,state)
    overlap=min(1.,abs(dot(reference,state))/sqrt(nr*ns))
    probabilities=abs2.(state)./ns
    number=0.;spin=0.
    alpha_count(b)=sum((b>>site)&1 for site in 0:2:n-1;init=0)
    for basis in 0:length(state)-1
        count_ones(UInt32(basis))==count_ones(UInt32(hf)) || (number+=probabilities[basis+1])
        (count_ones(UInt32(basis))==count_ones(UInt32(hf)) && alpha_count(basis)==alpha_count(Int(hf))) ||
            (spin+=probabilities[basis+1])
    end
    Dict("projective_error"=>1-overlap,"infidelity"=>1-overlap^2,"norm_drift"=>abs(ns-1),
        "number_leakage"=>number,"spin_leakage"=>spin)
end

@testset "Raw Pauli masks, complex Y phase and disjoint-pair rotation" begin
    for spelling in ("","0X","0Y","0Z","0Y,2Y","0X,2Y","0Y,1Z,2Y")
        word=IP.compile_word(spelling,3,big(1),big(1))
        matrix=dense_word(word.key,3)
        state=ComplexF64[complex(cos(i),sin(2i)) for i in 1:8];state./=norm(state)
        @test norm(word_action(word,state)-matrix*state)<2e-15
        for angle in (-.37,.29)
            actual=copy(state);rotate_word!(actual,word,angle)
            @test norm(actual-exp(-im*angle*matrix)*state)<3e-14
        end
    end
    word=IP.compile_word("13X",14,big(1),big(1))
    state=zeros(ComplexF64,16384);state[1]=1
    rotate_word!(state,word,.31)
    @test abs(state[1]-cos(.31))<1e-15
    @test abs(state[8193]+im*sin(.31))<1e-15
    @test count(value->!iszero(value),state)==2
    metrics=state_metrics([ComplexF64(i==1) for i in 1:16384],state,14,0)
    @test metrics["number_leakage"]≈sin(.31)^2 atol=2e-15
    @test metrics["spin_leakage"]≈sin(.31)^2 atol=2e-15
    @test_throws Refusal rotate_word!(zeros(ComplexF64,3),word,.1)
    @test_throws Refusal word_action(word,zeros(ComplexF64,8))
end

@testset "Independent sector Krylov reference versus full dense Kronecker exponential" begin
    mktempdir() do directory
        path=joinpath(directory,"model.toml");write_input(path);model=load_model(path)
        h=dense_hamiltonian(model);initial=zeros(ComplexF64,16);initial[Int(model.hf)+1]=1
        @test norm(h-adjoint(h))==0
        number=Diagonal(Float64[count_ones(UInt32(b)) for b in 0:15])
        alpha=Diagonal(Float64[count_ones(UInt32(b)&UInt32(5)) for b in 0:15])
        @test norm(number*h-h*number)<1e-14
        @test norm(alpha*h-h*alpha)<1e-14
        lone_xx=dense_word(((0,'X'),(1,'Z'),(2,'X')),4)
        @test norm(number*lone_xx-lone_xx*number)>0
        for t in (.4,1.6)
            result=reference_action(model,t)
            expected=exp(-im*t*h)*initial
            @test norm(result.state-expected)<5e-13
            @test result.diagnostic["global_phase_restored"]===true
            @test result.diagnostic["sector_dimension"]==6
            @test result.diagnostic["norm_drift"]<1e-12
        end
        @test_throws Refusal reference_action(model,1.6;max_krylov=1)
        @test_throws Refusal reference_action(model,.4;max_sector_dim=5)
        @test_throws Refusal reference_action(model,.4;max_sector_dim=2001)
        @test_throws Refusal reference_action(model,.4;check_budget=()->throw(Refusal("injected action limit")))
        # Hermitian imaginary hopping gives a separate odd-Y phase check of the
        # independent complex action, without borrowing any primary eigensystem.
        units=[Dict("words"=>["0X,1Z,2Y","0Y,1Z,2X"],"numerators"=>["30","-30"]),
               Dict("words"=>["0Z"],"numerators"=>["20"])]
        write_input(path;units=units);complex_model=load_model(path)
        expected=exp(-.4im*dense_hamiltonian(complex_model))*initial
        @test norm(reference_action(complex_model,.4).state-expected)<5e-13
    end
end

@testset "Both literal formulas, unequal lengths, original costs and unprojected circuits" begin
    mktempdir() do directory
        path=joinpath(directory,"model.toml");write_input(path);model=load_model(path)
        @test length.(model.units)==[2,1,2,1]
        for partition in ("conserving","pauli"),formula in (1,2),t in (.4,1.6)
            order=reverse(collect(1:length(IP.partition_units(model,partition))))
            result=unfused_state(model,partition,order,t,3;formula=formula)
            expected=dense_unfused(model,partition,order,t,3;formula=formula)
            @test norm(result.state-expected)<5e-13
            @test result.primitive_pulses==3length(model.words)*(formula==2 ? 2 : 1)
            plan=affordable(model,partition,order,24length(model.words);formula=formula)
            @test plan.ideal_rotations<=24length(model.words)
            @test plan.steps>=1
        end
        plan=affordable(model,"conserving",[4,2,1,3],21;formula=2)
        @test plan.steps==2 && plan.ideal_rotations==19
        @test_throws Refusal affordable(model,"conserving",[1,1,3,4],100)
        @test_throws Refusal affordable(model,"conserving",[1,2,3,4],0)
        @test_throws Refusal unfused_state(model,"conserving",[1,2,3,4],.4,3;max_primitive_pulses=1)
        # In N=4 on six modes, separate XX/YY pulses can visit N=2 even
        # though each original grouped unit and total Hamiltonian conserve N.
        write_input(path;n=6,electrons=4);leaking=load_model(path)
        xx=findfirst(w->w.key==((0,'X'),(1,'Z'),(2,'X')),leaking.words)
        yy=findfirst(w->w.key==((0,'Y'),(1,'Z'),(2,'Y')),leaking.words)
        zz=findfirst(w->w.key==((0,'Z'),),leaking.words)
        order=[xx,zz,yy,setdiff(collect(eachindex(leaking.words)),[xx,zz,yy])...]
        actual=unfused_state(leaking,"pauli",order,1.6,3;formula=1)
        reference=reference_action(leaking,1.6).state
        @test state_metrics(reference,actual.state,6,leaking.hf)["number_leakage"]>1e-10
    end
end

@testset "Exact conservation, canonical numeric labels and bounded admission" begin
    mktempdir() do directory
        path=joinpath(directory,"model.toml")
        data=write_input(path);good=load_model(path)
        @test good.input_sha256==bytes2hex(sha256(read(path)))
        @test isempty(IP.number_commutator(good.words))
        @test isempty(IP.number_commutator(good.words;alpha=true))
        # Remove one necessary YY coefficient: number conservation really fails.
        data["units"][1]=Dict("words"=>["0X,1Z,2X"],"numerators"=>["35"])
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws Refusal load_model(path)
        # Total N is conserved but interleaved alpha number is not.
        write_input(path;units=[Dict("words"=>["0X,1X","0Y,1Y"],"numerators"=>["35","35"])])
        @test_throws Refusal load_model(path)
        write_input(path;units=[Dict("words"=>["0X","0Z"],"numerators"=>["1","1"])])
        @test_throws Refusal load_model(path)
        write_input(path;units=[Dict("words"=>["0Z"],"numerators"=>["1"]),Dict("words"=>["0Z"],"numerators"=>["2"])])
        @test_throws Refusal load_model(path)
        write_input(path;n=14,electrons=4,units=[Dict("words"=>["10Z"],"numerators"=>["1"]),Dict("words"=>["2Z"],"numerators"=>["2"])])
        numeric=load_model(path)
        @test [word.key for word in numeric.words]==[((2,'Z'),),((10,'Z'),)]
        @test numeric.units==[[2],[1]]
        high_reference=reference_action(numeric,.4)
        @test high_reference.diagnostic["sector_dimension"]==1001
        @test abs(high_reference.state[16]-cis(.004))<2e-14
        @test count(value->!iszero(value),high_reference.state)==1
        @test_throws Refusal reference_action(numeric,.4;max_sector_dim=1000)
        write_input(path;n=14,electrons=7,units=[Dict("words"=>["13Z"],"numerators"=>["1"])])
        @test_throws Refusal reference_action(load_model(path),.4)
        write_input(path;n=15)
        @test_throws Refusal load_model(path)
        data=write_input(path);data["denominator"]=string(big(1)<<4096)
        open(path,"w") do io;TOML.print(io,data);end
        @test_throws Refusal load_model(path)
        write_input(path;n=14,units=[Dict("words"=>["10Z,2Z"],"numerators"=>["1"])])
        @test_throws Refusal load_model(path)
        write_input(path;units=[Dict("words"=>["0Z"],"numerators"=>["1"]) for _ in 1:1025])
        @test_throws Refusal load_model(path)
        write_input(path;units=[Dict("words"=>fill("0Z",8193),"numerators"=>fill("1",8193))])
        @test_throws Refusal load_model(path)
        write(path,"x="*repeat("[",9)*"0"*repeat("]",9))
        @test_throws Refusal IP.bounded_toml(path)
    end
end

@testset "Stable projective metric, little-endian/hash/path binding and malformed states" begin
    reference=ComplexF64[1,0,0,0];epsilon=1e-9
    state=ComplexF64[sqrt(1-epsilon^2),epsilon,0,0]
    metrics=state_metrics(reference,state,2,0)
    @test metrics["infidelity"]≈epsilon^2 rtol=1e-14
    @test metrics["projective_error"]≈epsilon^2/2 rtol=1e-14
    @test metrics["projective_error"]>0
    @test state_metrics(reference,cis(.37)*reference,2,0)["projective_error"]<1e-29
    mktempdir() do directory
        path=joinpath(directory,"state.bin");digest=write_state(path,state)
        @test read_state(directory,"state.bin",digest,2)==state
        @test_throws Refusal read_state(directory,"state.bin",repeat("0",64),2)
        @test_throws Refusal read_state(directory,"state.bin",digest,3)
        @test_throws Refusal read_state(directory,"../state.bin",digest,2)
        @test_throws Refusal read_state(directory,path,digest,2)
        @test_throws Refusal read_state(directory,"state.bin",digest,15)
        mktempdir() do outside
            foreign=joinpath(outside,"foreign.bin");foreign_hash=write_state(foreign,state)
            symlink(foreign,joinpath(directory,"escape.bin"))
            @test_throws Refusal read_state(directory,"escape.bin",foreign_hash,2)
        end
        state[1]=ComplexF64(NaN,0);digest=write_state(path,state)
        @test_throws Refusal read_state(directory,"state.bin",digest,2)
    end
end

function complete_panel(directory,model)
    h=dense_hamiltonian(model);initial=zeros(ComplexF64,1<<model.n);initial[Int(model.hf)+1]=1
    references=Dict{Float64,Any}();saved=Dict{Any,Any}();rows=Any[]
    for t in (.4,1.6)
        state=exp(-im*t*h)*initial;file="reference_"*string(t)*".bin"
        references[t]=(state=state,file=file,hash=write_state(joinpath(directory,file),state))
    end
    for identity in sort!(collect(IP.expected_inventory(1)))
        arm,repetition,formula,t,multiplier=identity
        partition=arm=="pauli" ? "pauli" : "conserving"
        order=collect(1:length(IP.partition_units(model,partition)))
        startswith(arm,"window") && reverse!(order)
        plan=affordable(model,partition,order,multiplier*length(model.words);formula=formula)
        key=(partition,Tuple(order),formula,t,plan.steps)
        if !haskey(saved,key)
            state=dense_unfused(model,partition,order,t,plan.steps;formula=formula)
            file="actual_"*string(length(saved)+1)*".bin"
            saved[key]=(state=state,file=file,hash=write_state(joinpath(directory,file),state))
        end
        actual=saved[key];reference=references[t]
        metrics=oracle_metrics(reference.state,actual.state,model.n,model.hf)
        row=Dict{String,Any}("name"=>arm,"arm"=>arm,"repetition"=>repetition,"partition"=>partition,
            "order"=>order,"formula"=>formula,"time"=>t,"budget_multiplier"=>multiplier,
            "rotation_budget"=>multiplier*length(model.words),"steps"=>plan.steps,"ideal_rotations"=>plan.ideal_rotations,
            "reference_file"=>reference.file,"reference_sha256"=>reference.hash,
            "actual_file"=>actual.file,"actual_sha256"=>actual.hash)
        for metric in ("projective_error","infidelity","norm_drift","number_leakage","spin_leakage");row[metric]=metrics[metric];end
        push!(rows,row)
    end
    Dict{String,Any}("status"=>"complete","physical_panel_complete"=>true,"mode"=>"preflight",
        "expected_logical_rows"=>44,"input_sha256"=>model.input_sha256,"case"=>model.case,
        "n"=>model.n,"hf"=>string(model.hf),"rows"=>rows)
end

@testset "Fixed matrix inventory, all saved metrics, scoped hard replay and incomplete CLI reports" begin
    mktempdir() do directory
        input=joinpath(directory,"model.toml");write_input(input);model=load_model(input)
        manifest=complete_panel(directory,model);path=joinpath(directory,"main.toml")
        open(path,"w") do io;TOML.print(io,manifest);end
        report=check_manifest(model,path)
        @test report["status"]=="complete"
        @test report["checked_logical_rows"]==44
        @test all(row->row["independent_unfused_replay"]==(row["time"]==1.6 && row["budget_multiplier"]==24),report["rows"])
        @test report["distinct_successful_unfused_replays"]==report["distinct_hard_circuit_identities"]
        output=joinpath(directory,"check.toml")
        @test IP.main([input,path,output])==0
        @test TOML.parsefile(output)["status"]=="complete"
        @test_throws Refusal IP.main([input,path,output])
        incomplete=deepcopy(manifest);pop!(incomplete["rows"]);incomplete["physical_panel_complete"]=false
        open(path,"w") do io;TOML.print(io,incomplete);end
        report=check_manifest(model,path)
        @test report["status"]=="incomplete" && report["checked_logical_rows"]==43
        output2=joinpath(directory,"incomplete.toml")
        @test IP.main([input,path,output2])==2
        duplicate=deepcopy(manifest);push!(duplicate["rows"],deepcopy(first(duplicate["rows"])))
        open(path,"w") do io;TOML.print(io,duplicate);end
        @test check_manifest(model,path)["status"]=="incomplete"
        bad=deepcopy(manifest);bad["rows"][1]["actual_sha256"]=repeat("0",64)
        open(path,"w") do io;TOML.print(io,bad);end
        report=check_manifest(model,path)
        @test report["status"]=="incomplete" && report["checked_logical_rows"]==43
        bad=deepcopy(manifest);bad["input_sha256"]=repeat("0",64)
        open(path,"w") do io;TOML.print(io,bad);end
        @test check_manifest(model,path)["status"]=="incomplete"
    end
end
