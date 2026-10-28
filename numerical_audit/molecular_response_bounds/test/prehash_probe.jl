# Standalone red/green admission regression. The supplied directory is read only;
# it may be the immutable historical source snapshot or the active backend.
# No oversized file is read: a String-specific hash spy intercepts that path.
using Test, TOML
length(ARGS)==1 || error("usage: prehash_probe.jl BACKEND_DIRECTORY")
const PROBE_ROOT=abspath(ARGS[1])
module AdmissionProbeFixture
include(joinpath(Main.PROBE_ROOT,"independent_check.jl"))
end
@eval AdmissionProbeFixture.MolecularResponseBounds.BoundsCommon begin
    const probe_watch=Set{String}()
    const probe_calls=Ref(0)
    function filehash(path::String)
        if path in probe_watch
            probe_calls[]+=1
            error("prehash probe intercepted oversized hash; no oversized read occurred")
        end
        invoke(filehash,Tuple{Any},path)
    end
end
const PROBE_COMMON=AdmissionProbeFixture.MolecularResponseBounds.BoundsCommon
@testset "Actual checker path refuses before any input/kernel hash" begin
    mktempdir() do dir
        smallinput=joinpath(dir,"input.toml");write(smallinput,"small input\n")
        smallkernel=joinpath(dir,"kernel.toml");write(smallkernel,"small kernel\n")
        result=joinpath(dir,"refusal.toml")
        # New admission must happen before reading even this minimal record.
        open(result,"w") do io
            TOML.print(io,Dict("complete"=>false,"status"=>"scope_limit","reason"=>"probe","positive_full_bound"=>false))
        end
        large=joinpath(dir,"sparse8MiBplus1.toml")
        open(large,"w") do io
            seek(io,PROBE_COMMON.Caps().max_input_bytes);write(io,UInt8(0))
        end
        @test filesize(large)==8*1024^2+1
        push!(PROBE_COMMON.probe_watch,large)
        for (input,kernel,name) in ((large,smallkernel,"input"),(smallinput,large,"kernel"))
            PROBE_COMMON.probe_calls[]=0
            # Ensure the oversized-kernel path does not hash its small input first.
            input==smallinput && push!(PROBE_COMMON.probe_watch,smallinput)
            checked=try
                AdmissionProbeFixture.check_one(input,kernel,result,joinpath(dir,name*".CHECK.toml"))
            catch e
                e
            end
            @test checked isa AbstractDict
            if checked isa AbstractDict
                @test checked["status"]=="input_limit"
                @test !checked["complete"] && !checked["scientific_verification_performed"]
                @test !checked["positive_certificate_verified"]
            end
            @test PROBE_COMMON.probe_calls[]==0
            delete!(PROBE_COMMON.probe_watch,smallinput)
        end
    end
end
