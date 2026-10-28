# Driver-only interception: never read the watched sparse file, and count actual
# scientific dispatch. Unwatched paths use the original generic implementation.
@eval CheckerFixture.MolecularResponseBounds.BoundsCommon begin
    const address_hash_watch=Set{String}()
    const address_hash_calls=Ref(0)
    function filehash(path::String)
        if path in address_hash_watch
            address_hash_calls[]+=1
            error("test intercepted forbidden pre-admission hash")
        end
        invoke(filehash,Tuple{Any},path)
    end
    const address_output_watch=Set{String}()
    function write_record(path::String,value;caps=Caps())
        if path in address_output_watch && !get(value,"output_refusal_receipt",false)
            throw(Refusal("output_limit","test final-write interception"))
        end
        invoke(write_record,Tuple{Any,Any},path,value;caps)
    end
end
@eval DriverFixture.MolecularResponseBounds.BoundsCommon begin
    const address_hash_watch=Set{String}()
    const address_hash_calls=Ref(0)
    function filehash(path::String)
        if path in address_hash_watch
            address_hash_calls[]+=1
            error("test intercepted forbidden pre-admission hash")
        end
        invoke(filehash,Tuple{Any},path)
    end
    const address_output_watch=Set{String}()
    function write_record(path::String,value;caps=Caps())
        if path in address_output_watch && !get(value,"output_refusal_receipt",false)
            throw(Refusal("output_limit","test final-write interception"))
        end
        invoke(write_record,Tuple{Any,Any},path,value;caps)
    end
end
@eval CheckerFixture.MolecularResponseBounds.IndependentCheck begin
    const address_scientific_calls=Ref(0)
    function independent_full_check(input::String,kernel::String,record;caps=Caps(),guard=()->nothing)
        address_scientific_calls[]+=1
        invoke(independent_full_check,Tuple{Any,Any,Any},input,kernel,record;caps,guard)
    end
end
const ADDRESS_CHECK_COMMON=CheckerFixture.MolecularResponseBounds.BoundsCommon
const ADDRESS_DRIVER_COMMON=DriverFixture.MolecularResponseBounds.BoundsCommon
const ADDRESS_SCIENCE=CheckerFixture.MolecularResponseBounds.IndependentCheck

function address_refusal(input,kernel)
    caps,decl=select_caps()
    Dict{String,Any}("schema"=>"molecular-bound-record-v1","complete"=>false,
        "status"=>"scope_limit","reason"=>"fresh current-source metadata-only test refusal",
        "positive_full_bound"=>false,"input_sha256"=>filehash(input),"kernel_sha256"=>filehash(kernel),
        "source_hashes"=>DriverFixture.sources(),"caps"=>caps_dict(caps),"caps_declaration"=>decl,
        "caps_profile"=>decl["effective_profile"],"caps_profile_validation"=>"validated_bounded_profile")
end
function address_sparse_file(path)
    open(path,"w") do io
        seek(io,Caps().max_input_bytes);write(io,UInt8(0))
    end
    path
end
@testset "Address: checker admits both byte sizes before hashes or scientific dispatch" begin
    mktempdir() do dir
        input=joinpath(dir,"input.toml");write(input,"small input bytes\n")
        kernel=joinpath(dir,"kernel.toml");write(kernel,"small kernel bytes\n")
        result=joinpath(dir,"refusal.toml");write_record(result,address_refusal(input,kernel))
        bigfile=address_sparse_file(joinpath(dir,"sparse-over-limit.toml"))
        @test filesize(bigfile)==Caps().max_input_bytes+1
        push!(ADDRESS_CHECK_COMMON.address_hash_watch,bigfile)
        # Also watch the small input in the oversized-kernel call: admission of
        # BOTH artifacts must precede hashing either artifact.
        for (a,k,name) in ((bigfile,kernel,"input"),(input,bigfile,"kernel"))
            ADDRESS_CHECK_COMMON.address_hash_calls[]=0;ADDRESS_SCIENCE.address_scientific_calls[]=0
            a==input && push!(ADDRESS_CHECK_COMMON.address_hash_watch,input)
            checked=CheckerFixture.check_one(a,k,result,joinpath(dir,name*".CHECK.toml"))
            @test checked["status"]=="input_limit" && !checked["complete"]
            @test checked["caps_validation_status"]=="not_checked"
            @test !checked["scientific_verification_performed"] && !checked["positive_certificate_verified"]
            @test !checked["artifact_links_authenticated"]
            @test ADDRESS_CHECK_COMMON.address_hash_calls[]==0
            @test ADDRESS_SCIENCE.address_scientific_calls[]==0
            delete!(ADDRESS_CHECK_COMMON.address_hash_watch,input)
        end
        empty!(ADDRESS_CHECK_COMMON.address_hash_watch)
        # The same admission happens for complete-looking records, before any
        # unsupported claim or absent scientific data can reach numerical code.
        claimed=address_refusal(input,kernel);claimed["complete"]=true
        claimedpath=joinpath(dir,"claimed.toml");write_record(claimedpath,claimed)
        checked=CheckerFixture.check_one(bigfile,kernel,claimedpath,joinpath(dir,"claimed.CHECK.toml"))
        @test checked["status"]=="input_limit" && ADDRESS_SCIENCE.address_scientific_calls[]==0
    end
end
@testset "Address: available refusal links are enforced, unavailable links are explicit" begin
    mktempdir() do dir
        input=joinpath(dir,"input.toml");write(input,"original input\n")
        kernel=joinpath(dir,"kernel.toml");write(kernel,"original kernel\n")
        record=address_refusal(input,kernel);result=joinpath(dir,"refusal.toml");write_record(result,record)
        ADDRESS_SCIENCE.address_scientific_calls[]=0
        checked=CheckerFixture.check_one(input,kernel,result,joinpath(dir,"bound.CHECK.toml"))
        @test checked["status"]=="recorded_refusal" && checked["artifact_links_authenticated"]
        @test checked["input_binding_status"]==checked["kernel_binding_status"]=="matched"
        @test !checked["complete"] && !checked["scientific_verification_performed"] && !checked["positive_certificate_verified"]
        replacementinput=joinpath(dir,"different-input.toml");write(replacementinput,"substituted input\n")
        replacementkernel=joinpath(dir,"different-kernel.toml");write(replacementkernel,"substituted kernel\n")
        for (a,k,name,key) in ((replacementinput,kernel,"input","input"),(input,replacementkernel,"kernel","kernel"))
            bad=CheckerFixture.check_one(a,k,result,joinpath(dir,name*"-substitution.CHECK.toml"))
            @test bad["status"]=="verification_failed" && bad[key*"_binding_status"]=="mismatch"
            @test !bad["artifact_links_authenticated"] && !bad["scientific_verification_performed"] && !bad["positive_certificate_verified"]
        end
        predigest=deepcopy(record);delete!(predigest,"input_sha256");delete!(predigest,"kernel_sha256")
        predigestpath=joinpath(dir,"predigest.toml");write_record(predigestpath,predigest)
        unbound=CheckerFixture.check_one(input,kernel,predigestpath,joinpath(dir,"predigest.CHECK.toml"))
        @test unbound["status"]=="recorded_refusal_unbound" && unbound["artifact_binding_status"]=="unbound"
        @test unbound["input_binding_status"]==unbound["kernel_binding_status"]=="producer_digest_unavailable"
        @test !unbound["artifact_links_authenticated"] && !unbound["scientific_verification_performed"]
        for k in ("-",joinpath(dir,"not-created.kernel.toml"))
            name=k=="-" ? "sentinel" : "missing"
            absent=CheckerFixture.check_one(input,k,result,joinpath(dir,name*".CHECK.toml"))
            @test absent["status"]=="recorded_refusal_unbound" && absent["kernel_binding_status"]=="artifact_unavailable"
            @test !absent["artifact_links_authenticated"] && !absent["positive_certificate_verified"]
        end
        @test CheckerFixture.main([input,"-",result,joinpath(dir,"unbound-cli.CHECK.toml")])==2
        claimed=deepcopy(record);claimed["complete"]=true
        claimedpath=joinpath(dir,"claimed.toml");write_record(claimedpath,claimed)
        absent=CheckerFixture.check_one(input,"-",claimedpath,joinpath(dir,"claimed.CHECK.toml"))
        @test absent["status"]=="verification_failed" && !absent["scientific_verification_performed"]
        damaged=deepcopy(record);damaged["source_hashes"]["study.jl"]=repeat("0",64)
        damagedpath=joinpath(dir,"wrong-source.toml");write_record(damagedpath,damaged)
        rejected=CheckerFixture.check_one(input,kernel,damagedpath,joinpath(dir,"wrong-source.CHECK.toml"))
        @test rejected["status"]=="verification_failed" && !rejected["scientific_verification_performed"]
        @test ADDRESS_SCIENCE.address_scientific_calls[]==0
    end
end
@testset "Address: producer rejects oversized inputs/kernel before hashing and continues" begin
    mktempdir() do dir
        input,_=source_file(dir,"admission",1,0,[(String["0X"],[1]),(String["0Z"],[1])];den=1)
        kernel=joinpath(dir,"zero.kernel.toml");freeze_kernel(input,kernel;time=0)
        bigfile=address_sparse_file(joinpath(dir,"oversized.toml"))
        push!(ADDRESS_DRIVER_COMMON.address_hash_watch,bigfile);ADDRESS_DRIVER_COMMON.address_hash_calls[]=0
        registry=Dict("schema"=>"molecular-bounds-registry-v1",
            "options"=>Dict("iterations"=>0,"max_witnesses"=>0,"target_count"=>1,"comparison_count"=>0),
            "cases"=>[Dict("id"=>"biginput","input"=>bigfile,"kernel"=>kernel,"order"=>[1,2]),
                Dict("id"=>"bigkernel","input"=>input,"kernel"=>bigfile,"order"=>[1,2]),
                Dict("id"=>"later","input"=>input,"kernel"=>kernel,"order"=>[1,1])])
        regpath=joinpath(dir,"registry.toml");write_record(regpath,registry)
        out=joinpath(dir,"out");summary=DriverFixture.run_registry(regpath,out)
        @test length(summary["cases"])==3 && !summary["all_cases_completed"]
        for id in ("biginput","bigkernel")
            row=TOML.parsefile(joinpath(out,id*".RESULT.toml"))
            @test row["status"]=="input_limit" && !row["complete"] && !row["positive_full_bound"]
            @test !haskey(row,id=="biginput" ? "input_sha256" : "kernel_sha256")
            @test !isfile(joinpath(out,id*".PHASE.model.toml"))
        end
        @test TOML.parsefile(joinpath(out,"later.RESULT.toml"))["status"]=="error"
        @test ADDRESS_DRIVER_COMMON.address_hash_calls[]==0
        empty!(ADDRESS_DRIVER_COMMON.address_hash_watch)
    end
end
@testset "Address: output-limit receipts use the same cap and preserve I/O errors" begin
    mktempdir() do dir
        caps=Caps(max_output_bytes=4096)
        path=joinpath(dir,"limited.toml")
        payload=Dict("schema"=>"test-large","status"=>"purported_success","complete"=>true,
            "positive_full_bound"=>true,"padding"=>repeat("x",8192))
        receipt=write_final_record(path,payload;caps)
        @test receipt["status"]=="output_limit" && !receipt["complete"] && receipt["output_refusal_receipt"]
        @test !receipt["result_payload_persisted"] && !receipt["scientific_verification_performed"]
        @test !receipt["positive_full_bound"] && !receipt["positive_certificate_verified"]
        @test filesize(path)<=caps.max_output_bytes && TOML.parsefile(path)["status"]=="output_limit"
        @test_throws ErrorException write_final_record(path,payload;caps)
        tiny=joinpath(dir,"too-small.toml")
        err=try write_final_record(tiny,payload;caps=Caps(max_output_bytes=16));nothing catch e;e end
        @test err isa Refusal && err.status=="output_limit" && !ispath(tiny)
        parent=joinpath(dir,"regular-file");write(parent,"not a directory")
        err=try write_final_record(joinpath(parent,"child.toml"),Dict("small"=>true);caps);nothing catch e;e end
        @test err!==nothing && !(err isa Refusal)
    end
end
@testset "Address: both CLI batches continue after final output refusal" begin
    mktempdir() do dir
        input=joinpath(dir,"input.toml");write(input,"small input\n")
        kernel=joinpath(dir,"kernel.toml");write(kernel,"small kernel\n")
        result=joinpath(dir,"refusal.toml");write_record(result,address_refusal(input,kernel))
        firstout=joinpath(dir,"first.CHECK.toml");secondout=joinpath(dir,"second.CHECK.toml")
        push!(ADDRESS_CHECK_COMMON.address_output_watch,firstout)
        status=CheckerFixture.main([input,kernel,result,firstout,input,kernel,result,secondout])
        @test status==2 && TOML.parsefile(firstout)["status"]=="output_limit"
        @test TOML.parsefile(secondout)["status"]=="recorded_refusal"
        empty!(ADDRESS_CHECK_COMMON.address_output_watch)
        source,_=source_file(dir,"output-driver",1,0,[(String["0X"],[1]),(String["0Z"],[1])];den=1)
        kp=joinpath(dir,"zero.kernel.toml");freeze_kernel(source,kp;time=0)
        registry=Dict("schema"=>"molecular-bounds-registry-v1",
            "options"=>Dict("iterations"=>0,"max_witnesses"=>0,"target_count"=>1,"comparison_count"=>0),
            "cases"=>[Dict("id"=>"first","input"=>source,"kernel"=>kp,"order"=>[1,1]),
                Dict("id"=>"second","input"=>source,"kernel"=>kp,"order"=>[1,1])])
        regpath=joinpath(dir,"registry.toml");write_record(regpath,registry)
        out=joinpath(dir,"out");push!(ADDRESS_DRIVER_COMMON.address_output_watch,joinpath(out,"first.RESULT.toml"))
        summary=DriverFixture.run_registry(regpath,out)
        @test summary["status"]=="inventory_complete" && length(summary["cases"])==2
        @test summary["cases"][1]["status"]=="output_limit" && summary["cases"][2]["status"]=="error"
        @test TOML.parsefile(joinpath(out,"first.RESULT.toml"))["output_refusal_receipt"]
        @test isfile(joinpath(out,"second.RESULT.toml")) && isfile(joinpath(out,"RUN.RESULT.toml"))
        empty!(ADDRESS_DRIVER_COMMON.address_output_watch)
    end
end
