using Test, LinearAlgebra, TOML, SHA
include(joinpath(@__DIR__,"..","src","GeneralizationPhysical.jl"))
include(joinpath(@__DIR__,"..","src","GeneralizationIndependent.jl"))
include(joinpath(@__DIR__,"fixtures.jl"))
const Primary=GeneralizationPhysical
const Check=GeneralizationIndependent

# Keep each malformed fixture isolated without copying every nested row.
function copy_candidate_rows(manifest, indices=())
    copied=copy(manifest)
    copied["candidates"]=copy(manifest["candidates"])
    for index in indices
        copied["candidates"][index]=copy(manifest["candidates"][index])
    end
    copied
end

function physical_test_start(name)
    println("START physical admission: ",name)
    flush(stdout)
end

function check_main_roster(admit, name)
    Base.@nospecialize admit
    @testset "$name" begin
        manifest=fixture_manifest()
        @test admit(manifest,2) === manifest
        bad=copy_candidate_rows(manifest);pop!(bad["candidates"])
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest);bad["candidates"][end]=copy(bad["candidates"][1])
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(1,));bad["candidates"][1]["order"]=[0,1]
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(1,));bad["candidates"][1]["formula"]=1
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest);bad["training_budget_multiplier"]=24
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(5,));bad["candidates"][5]["parent_sha256"]=fixture_digest("invented cheap prefix")
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(13,));bad["candidates"][13]["parent_sha256"]=fixture_digest("different parent")
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(2,));bad["candidates"][2]["claimed_score"]=2e-6
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest,(1,));bad["candidates"][1]["constructor_score"]=bad["candidates"][1]["claimed_score"]
        @test_throws ArgumentError admit(bad,2) # Alias presence is constructor metadata too.
        indices=findall(row -> row["arm"]=="annealing" && row["cpu_budget"]==64,manifest["candidates"])
        fallback=copy_candidate_rows(manifest,indices)
        for index in indices
            row=fallback["candidates"][index]
            row["eligible"]=false;row["constructor_failure"]="finite CPU refusal; valid scored parent retained"
        end
        @test length(admit(fallback,2)["candidates"])==72
        # The complete fixed-N sector, not just the initial spin block, is admitted.
        @test binomial(12,8)==495
        bad=fixture_manifest(n=14,electrons=7)
        @test_throws ArgumentError admit(bad,2)
    end
end

function check_preflight_roster(admit, name)
    Base.@nospecialize admit
    @testset "$name" begin
        manifest=fixture_manifest(mode="preflight",n=14,electrons=4)
        @test admit(manifest,2) === manifest
        @test length(manifest["candidates"])==36
        bad=copy_candidate_rows(manifest);pop!(bad["candidates"])
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest);push!(bad["candidates"],manifest["candidates"][1])
        @test_throws ArgumentError admit(bad,2)
        bad=copy_candidate_rows(manifest);bad["expected_logical_rows"]=72
        @test_throws ArgumentError admit(bad,2)
    end
end

function check_artifact_bytes()
    @testset "artifact bytes and safe DATA paths" begin
        mktempdir() do root
            manifest=fixture_manifest();case=manifest["case"]
            mkpath(joinpath(root,"inputs",case));input=joinpath(root,"inputs",case,"model.toml")
            write(input,"n = 12\n");write(joinpath(dirname(input),"model.json"),"{}\n")
            manifest["input_sha256"]=Primary.filehash(input)
            manifest["source_json_sha256"]=Primary.filehash(joinpath(dirname(input),"model.json"))
            manifest["protocol_path"]="source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json"
            manifest["source_freeze_path"]="INPUT_PROTOCOL_FROZEN.json"
            for (pathkey,hashkey) in (("protocol_path","protocol_sha256"),("source_freeze_path","source_freeze_sha256"))
                path=joinpath(root,manifest[pathkey]);mkpath(dirname(path));write(path,"{}\n");manifest[hashkey]=Primary.filehash(path)
            end
            artifacts=Any[];digests=Dict{String,String}()
            for row in manifest["candidates"]
                oid=row["order_id"]
                if !haskey(digests,oid)
                    relative="constructors/$(case).rep$(row["repetition"]).$(row["arm"]).B$(row["cpu_budget"]).json"
                    path=joinpath(root,relative);mkpath(dirname(path));write(path,"{\"fixture\":\"$(oid)\"}\n")
                    digests[oid]=Primary.filehash(path)
                    push!(artifacts,Dict("order_id"=>oid,"path"=>relative,"sha256"=>digests[oid]))
                end
            end
            for row in manifest["candidates"]
                row["constructor_sha256"]=digests[row["order_id"]]
                row["arm"] in ("canonical","magnitude") || (row["parent_sha256"]=digests["rep$(row["repetition"]).parent.B0"])
            end
            manifest["constructor_artifacts"]=artifacts
            @test Primary.bind_artifacts(input,manifest).hashes==Check.bind_artifacts(input,manifest).hashes
            @test_throws ArgumentError Primary.artifact_path(root,"../outside")
            @test_throws ArgumentError Check.artifact_path(root,input)
            write(joinpath(root,first(artifacts)["path"]),"tampered\n")
            @test_throws ArgumentError Primary.bind_artifacts(input,manifest)
            @test_throws ArgumentError Check.bind_artifacts(input,manifest)
        end
    end
end

function check_full_register()
    @testset "full-register propagation, complete N reference and leakage visibility" begin
        Unit=Primary.ResponseCore.Unit
        hopping=Unit(UInt64(5),UInt64[0,5],Int8[1,-1],[.5,.5])
        diagonal=Unit(UInt64(0),UInt64[1],Int8[1],[.3])
        units=[hopping,diagonal];hf=UInt64(3)
        reference=Primary.WindowPhysical.reference_sector(units,4,hf)
        @test length(reference.basis)==6 # Fixed spin would contain four states.
        exact=Primary.WindowPhysical.reference_state(reference,.7).state
        X=ComplexF64[0 1;1 0];Y=ComplexF64[0 -im;im 0];Z=ComplexF64[1 0;0 -1];I2=Matrix{ComplexF64}(I,2,2)
        H=.5kron(I2,X,I2,X)+.5kron(I2,Y,I2,Y)+.3kron(I2,I2,I2,Z)
        initial=zeros(ComplexF64,16);initial[4]=1
        @test norm(exact-exp(-im*.7*H)*initial)<1e-12
        scratch=zeros(ComplexF64,16)
        Primary.WindowPhysical.apply_word_rotation!(scratch,initial,UInt64(1),UInt64(0),1,.3)
        leakage=Primary.WindowPhysical.state_leakage(scratch,4,hf)
        @test leakage.number_leakage≈sin(.3)^2 atol=1e-14
        @test_throws ArgumentError Primary.WindowPhysical.reference_sector([Unit(UInt64(1),UInt64[0],Int8[1],[1.])],4,hf)
        @test_throws Check.IndependentPhysical.Refusal Check.IndependentPhysical.state_metrics(initial,2initial,4,hf)
        mktempdir() do root
            path=joinpath(root,"state.bin");digest=Primary.save_state(path,scratch)
            @test Check.IndependentPhysical.read_state(root,"state.bin",digest,4)==scratch
            @test_throws ArgumentError Primary.save_state(path,scratch)
            write(path,zeros(UInt8,256))
            @test_throws Check.IndependentPhysical.Refusal Check.IndependentPhysical.read_state(root,"state.bin",digest,4)
        end
    end
end

function check_incomplete_inventory()
    @testset "old-case admission abort retains every logical row and reference" begin
        mktempdir() do root
            manifest=fixture_manifest(mode="preflight",n=14,electrons=4)
            case="LiH_r1p487_631g_14m4e";manifest["case"]=case
            input=joinpath(root,"inputs",case,"model.toml");mkpath(dirname(input))
            json=joinpath(dirname(input),"model.json");write(json,"{}\n")
            model=Dict("case"=>case,"n"=>14,"electrons"=>4,"partition"=>"conserving","denominator"=>"10",
                "source_json_sha256"=>Primary.filehash(json),"baseline_order"=>[1,2],
                "units"=>[Dict("words"=>["0Z"],"numerators"=>["3"]),Dict("words"=>["2Z"],"numerators"=>["4"])])
            open(input,"w") do io;TOML.print(io,model);end
            manifest["input_sha256"]=Primary.filehash(input);manifest["source_json_sha256"]=Primary.filehash(json)
            manifest["protocol_path"]="source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json"
            manifest["source_freeze_path"]="INPUT_PROTOCOL_FROZEN.json"
            for (pathkey,hashkey) in (("protocol_path","protocol_sha256"),("source_freeze_path","source_freeze_sha256"))
                path=joinpath(root,manifest[pathkey]);mkpath(dirname(path));write(path,"{}\n");manifest[hashkey]=Primary.filehash(path)
            end
            artifacts=Any[];digests=Dict{String,String}()
            for row in manifest["candidates"]
                oid=row["order_id"]
                if !haskey(digests,oid)
                    relative="constructors/$(case).rep$(row["repetition"]).$(row["arm"]).B$(row["cpu_budget"]).json"
                    path=joinpath(root,relative);mkpath(dirname(path));write(path,"{}\n")
                    digests[oid]=Primary.filehash(path)
                    push!(artifacts,Dict("order_id"=>oid,"path"=>relative,"sha256"=>digests[oid]))
                end
            end
            for row in manifest["candidates"]
                row["constructor_sha256"]=digests[row["order_id"]]
                row["arm"] in ("canonical","magnitude") || (row["parent_sha256"]=digests["rep1.parent.B0"])
            end
            manifest["constructor_artifacts"]=artifacts
            panel=joinpath(root,"panels",case*".toml");mkpath(dirname(panel))
            open(panel,"w") do io;TOML.print(io,manifest);end
            @test Primary.admit(input,panel).manifest["expected_logical_rows"]==36
            independent_model=Check.IndependentPhysical.load_model(input)
            @test Check.admit_candidates(independent_model,panel;input=input).data["expected_logical_rows"]==36
            moved=joinpath(root,"copied.toml");cp(panel,moved)
            @test_throws ArgumentError Primary.admit(input,moved)
            @test_throws ArgumentError Check.admit_candidates(independent_model,moved;input=input)
            calls=Ref(0)
            guard=() -> begin
                calls[]+=1
                calls[]>=2 && throw(ArgumentError("intentional pre-reference stop"))
            end
            output=joinpath(root,"physical")
            report=Primary.run(input,panel,output;check_budget=guard)
            @test report["status"]=="incomplete"
            @test report["completed_logical_rows"]==0
            @test length(report["rows"])==36
            @test all(row -> row["status"]=="incomplete",report["rows"])
            @test length(report["reference_records"])==2
            @test report["completed_references"]==0
            @test !isempty(report["failures"])
            @test isempty(filter(name -> endswith(name,".bin"),readdir(output)))
            @test_throws ArgumentError Primary.run(input,panel,output)
            verification=Check.run(input,panel,output,joinpath(root,"verify"))
            @test verification["status"]=="incomplete"
            @test verification["completed_logical_rows"]==0
            @test length(verification["rows"])==36
            @test Set(row["name"] for row in verification["rows"])==Set(row["name"] for row in manifest["candidates"])
            # JSON fixture bytes deliberately exercise binding only; the separate
            # Python semantic consumer must refuse these as scientific artifacts.
        end
    end
end

function check_bounded_reads()
    @testset "observed-size hash/TOML reads" begin
        mktempdir() do root
            path=joinpath(root,"tiny.toml");raw=Vector{UInt8}(codeunits("n = 4\n"));write(path,raw)
            digest=bytes2hex(sha256(raw))
            for adapter in (Primary,Check)
                @test adapter.bounded_hash(path)==digest
                parsed=adapter.bounded_toml(path)
                @test parsed.sha256==digest && parsed.data==Dict("n"=>4)
                @test adapter.observed_read(IOBuffer(copy(raw)),length(raw))==raw
                @test_throws ArgumentError adapter.observed_read(IOBuffer(raw[1:end-1]),length(raw))
                @test_throws ArgumentError adapter.observed_read(IOBuffer(vcat(raw,UInt8[0])),length(raw))
            end
        end
    end
end

function run_physical_admission_tests()
    physical_test_start("primary main roster")
    @test length(Primary.expected_grid("main"))==72
    check_main_roster(Primary.admit_manifest,"primary main roster and parent/cheap graph")
    physical_test_start("independent main roster")
    check_main_roster(Check.admit_manifest,"independent main roster and parent/cheap graph")
    physical_test_start("primary preflight roster/count")
    check_preflight_roster(Primary.admit_manifest,"primary preflight roster/count")
    physical_test_start("independent preflight roster/count")
    @test length(Check.expected_grid("preflight"))==36
    check_preflight_roster(Check.admit_manifest,"independent preflight roster/count")
    physical_test_start("artifact bytes and safe DATA paths")
    check_artifact_bytes()
    physical_test_start("full-register propagation, complete N reference and leakage visibility")
    check_full_register()
    physical_test_start("old-case admission abort retains every logical row and reference")
    check_incomplete_inventory()
    physical_test_start("observed-size hash/TOML reads")
    check_bounded_reads()
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_physical_admission_tests()
end
