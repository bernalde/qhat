const CR=Complex{Rational{BigInt}}
function source_file(dir,name,n,electrons,blocks;den=8)
    path=joinpath(dir,name*".toml")
    data=Dict("n"=>n,"electrons"=>electrons,"denominator"=>string(den),"case"=>name,
        "cubic_order"=>collect(1:length(blocks)),
        "units"=>[Dict("words"=>words,"numerators"=>string.(nums)) for (words,nums) in blocks])
    write_record(path,data);path,data
end
function supplied_kernel(dir,name,input,n,hf)
    path=joinpath(dir,name*".kernel.toml")
    write_record(path,Dict("schema"=>"supplied-dyadic-kernel-v1","input_sha256"=>filehash(input),
        "n"=>n,"hf"=>string(hf),"time"=>1.0,"denominator"=>"16",
        "real"=>string.([7+i for i in 0:2^n-1]),"imag"=>string.([3-i for i in 0:2^n-1])))
    path
end
function literal_word(word,n)
    identity=CR[1 0;0 1];x=CR[0 1;1 0];y=CR[0 -im;im 0];z=CR[1 0;0 -1]
    letters=Dict{Int,Char}()
    for token in (isempty(word) ? String[] : split(word,','));letters[parse(Int,token[1:end-1])]=token[end];end
    out=fill(CR(1),1,1)
    for q in n-1:-1:0
        p=get(letters,q,'I');out=kron(out,p=='I' ? identity : p=='X' ? x : p=='Y' ? y : z)
    end
    out
end
function literal_units(data)
    n=data["n"];den=parse(BigInt,data["denominator"])
    [sum((parse(BigInt,num)//den)*literal_word(word,n) for (word,num) in zip(b["words"],b["numerators"])) for b in data["units"]]
end
"""Independent dense polynomial multiplication, no sparse action helper."""
function dense_defect(matrices,order;physical=false)
    d=size(first(matrices),1);identity=Matrix{CR}(I,d,d)
    generators=physical ? [-im*A for A in matrices] : matrices
    coeff=[identity,zeros(CR,d,d),zeros(CR,d,d),zeros(CR,d,d)]
    stages=[(i,big(1)//2) for i in order[1:end-1]]
    push!(stages,(last(order),big(1)//1));append!(stages,[(i,big(1)//2) for i in reverse(order[1:end-1])])
    for (label,tau) in stages
        stage=[(tau^p/factorial(p))*(generators[label]^p) for p in 0:3]
        coeff=[sum(stage[p+1]*coeff[j-p+1] for p in 0:j) for j in 0:3]
    end
    H=sum(generators);coeff[4]-(big(1)//big(6))*(H^3)
end
function atom_sum(model,order)
    b=Budget();sumcol=Dict{UInt64,BigInt}()
    for a in 1:length(order)-1,c in a+1:length(order)
        for (state,v) in atom_integer(model,order[[a,c]],b);sumcol[state]=get(sumcol,state,big(0))+v;end
    end
    for a in 1:length(order)-2,c in a+1:length(order)-1,d in c+1:length(order)
        for (state,v) in atom_integer(model,order[[a,c,d]],b);sumcol[state]=get(sumcol,state,big(0))+v;end
    end
    filter!(p->!iszero(last(p)),sumcol);sumcol
end
@testset "Exact phase/sign48 vs24 and independent dense polynomial identities" begin
    mktempdir() do dir
        fixtures=[("xz",1,0,[(String["0X"],[8]),(String["0Z"],[8])]),
            ("rational",1,0,[(String["0X"],[3]),(String["0Z"],[-2])]),
            ("conserving",3,1,[(String["0X,1X","0Y,1Y"],[2,2]),
                (String["1X,2X","1Y,2Y"],[-3,-3]),(String["0Z"],[5])])]
        for (name,n,electrons,blocks) in fixtures
            input,data=source_file(dir,name,n,electrons,blocks;den=(name=="rational" ? 5 : 8));hf=UInt64(2)^electrons-1
            kp=supplied_kernel(dir,name,input,n,hf);model=load_exact_model(input,kp);matrices=literal_units(data)
            for order in tiny_orders(collect(eachindex(model.units)))
                C=dense_defect(matrices,order);physical=dense_defect(matrices,order;physical=true)
                @test physical==im*C
                raw=atom_sum(model,order);ub=exact_full_response(model,order)
                expected=C[:,Int(hf)+1];scale=big(48)*model.coefficient_denominator^3
                @test [get(raw,UInt64(j-1),big(0))//scale for j in 1:2^n]==expected
                @test [get(ub.column,UInt64(j-1),big(0))//scale for j in 1:2^n]==expected
                @test all(iseven,values(raw)) # Natural Lie denominator is24DA³.
                Q=big(0)//big(1)
                for j in 1:2^n
                    UInt64(j-1)==hf && continue
                    re,impart=model.kernel[j];K=CR(re//model.kernel_denominator,impart//model.kernel_denominator)
                    Q+=abs2(expected[j]*K)
                end
                @test Q==ub.Q==ub.numerator//ub.denominator
            end
            if name=="xz"
                @test dense_defect(matrices,[1,2])==(big(1)//big(3))*matrices[1]-(big(1)//big(6))*matrices[2]
            end
            record=molecular_record(input,kp,collect(eachindex(model.units)),30sum(original_lengths(model));iterations=8,comparison_count=1)
            verified=independent_full_check(input,kp,record)
            @test verified.status=="independently_verified"
            @test verified.L_J<=verified.UB_J
            allJ=[exact_full_response(model,p).Q/endpoint_steps(original_lengths(model),first(p),last(p),30sum(original_lengths(model)))^4
                for p in tiny_orders(collect(eachindex(model.units)))]
            @test verified.L_J<=minimum(allJ)<=verified.UB_J
            stored=joinpath(dir,name*".record.toml");write_record(stored,record)
            parsed=TOML.parsefile(stored)
            @test independent_full_check(input,kp,parsed).L_J==verified.L_J
            for key in ("input_sha256","kernel_sha256","coverage_scope_sha256","projection_sha256")
                bad=deepcopy(record);bad[key]="tampered"
                @test_throws Exception independent_full_check(input,kp,bad)
            end
            bad=deepcopy(record);bad["UB_numerator"]+=1
            @test_throws Exception independent_full_check(input,kp,bad)
            bad=deepcopy(record);bad["r_max"]+=1
            @test_throws Exception independent_full_check(input,kp,bad)
            damaged=deepcopy(data);damaged["units"][1]["numerators"][1]=string(parse(BigInt,damaged["units"][1]["numerators"][1])+1)
            badinput=joinpath(dir,name*".badinput.toml");write_record(badinput,damaged)
            @test_throws Exception independent_full_check(badinput,kp,record)
            damagedkernel=TOML.parsefile(kp);damagedkernel["real"][1]=string(parse(BigInt,damagedkernel["real"][1])+1)
            badkernel=joinpath(dir,name*".badkernel.toml");write_record(badkernel,damagedkernel)
            @test_throws Exception independent_full_check(input,badkernel,record)
        end
    end
end
@testset "Complete projection covers both pair outputs, duplicate flips and diagonals" begin
    units=[ExactUnit(UInt64(3),UInt64[0,3],Int8[1,-1],BigInt[1,1]),
        ExactUnit(UInt64(6),UInt64[0,6],Int8[1,-1],BigInt[1,1]),
        ExactUnit(UInt64(0),UInt64[1],Int8[1],BigInt[2]),
        ExactUnit(UInt64(0),UInt64[4],Int8[1],BigInt[-1]),
        ExactUnit(UInt64(3),UInt64[0,3],Int8[1,-1],BigInt[3,3])]
    model=ExactModel(3,UInt64(2),units,big(2),big(1),[(big(1),big(0)) for _ in 1:8],"manual","manual",1.0)
    atom=atom_integer(model,[1,2],Budget())
    @test haskey(atom,UInt64(1)) && haskey(atom,UInt64(4))
    targets=UInt64[1,4,7];projection=project_full_atoms(model,targets)
    expected=Set{Tuple{Vararg{Int}}}()
    for i in 1:4,j in i+1:5
        (model.hf⊻units[i].flip in targets || model.hf⊻units[j].flip in targets) && push!(expected,(i,j))
    end
    for i in 1:3,j in i+1:4,k in j+1:5
        model.hf⊻units[i].flip⊻units[j].flip⊻units[k].flip in targets && push!(expected,(i,j,k))
    end
    @test scope_digest(expected)==projection.scope_sha256
    @test length(expected)==projection.scopes
    for order in tiny_orders(collect(1:5))
        full=exact_full_response(model,order);selected=big(0)
        for d in targets
            x,y=get(full.response,d,(big(0),big(0)));selected+=x*x+y*y
        end
        @test exact_value(projection.response,order)==selected//full.denominator
        cert=bound_response(projection.response;iterations=0,max_witnesses=0)
        @test verify_bound(projection.response,cert)<=full.Q
    end
    @test_throws Refusal project_full_atoms(model,targets;caps=Caps(max_scopes=1))
    @test_throws ArgumentError project_full_atoms(model,[model.hf])
end
@testset "Original integer coefficients and word counts survive; caps and endpoint bounds" begin
    mktempdir() do dir
        den=big(2)^60;input,_=source_file(dir,"precision",1,0,[(String["0X"],[den+1]),(String["0Z"],[den])];den)
        kp=supplied_kernel(dir,"precision",input,1,0);m=load_exact_model(input,kp)
        @test m.units[1].numerators[1]==den+1
        @test Float64((den+1)//den)==1.0
        @test original_lengths(m)==[1,1]
        @test exact_full_response(m,[1,2];caps=Caps(max_cache=1)).Q==exact_full_response(m,[1,2]).Q
        @test_throws Refusal read_coefficients(input;caps=Caps(max_input_bytes=1))
        @test_throws Refusal read_coefficients(input;caps=Caps(max_bits=4))
        @test_throws Refusal exact_full_response(m,[1,2];caps=Caps(max_work=1))
        @test_throws ErrorException exact_full_response(m,[1,2];guard=()->error("deadline"))
        @test_throws ArgumentError exact_full_response(m,[1,1])
        zero_kernel=joinpath(dir,"zero.kernel.toml");freeze_kernel(input,zero_kernel;time=0.0)
        zero_record=molecular_record(input,zero_kernel,[1,2],60;iterations=0,max_witnesses=0)
        @test zero_record["complete"] && zero_record["status"]=="zero_full_bound"
        @test zero_record["Q_inc"]==0 && zero_record["L_Q"]==0 && isnothing(zero_record["ratio"])
        @test independent_full_check(input,zero_kernel,zero_record).L_Q==0
        @test_throws ErrorException freeze_kernel(input,zero_kernel)
        broken=TOML.parsefile(kp);pop!(broken["imag"])
        partial=joinpath(dir,"partial.kernel.toml");write_record(partial,broken)
        @test_throws ArgumentError load_exact_model(input,partial)
    end
    lengths=[2,1,3,2]
    for B in (0,8,9,30,90)
        expected=maximum(endpoint_steps(lengths,first(p),last(p),B) for p in tiny_orders(collect(1:4)))
        @test maximum_steps(lengths,B)==expected
        for p in tiny_orders(collect(1:4))
            r=endpoint_steps(lengths,first(p),last(p),B)
            if r>0
                a=lengths[first(p)];c=lengths[last(p)];U=sum(lengths)
                @test 2U*r-(r-1)*a-r*c<=B
                @test 2U*(r+1)-r*a-(r+1)*c>B
            end
        end
    end
    @test maximum_steps([3],2)==0 && maximum_steps([3],3)==1
end
