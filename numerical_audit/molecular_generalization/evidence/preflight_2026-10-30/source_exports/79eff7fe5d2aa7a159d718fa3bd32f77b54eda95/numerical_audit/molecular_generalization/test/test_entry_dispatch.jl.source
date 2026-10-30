using Test
include(joinpath(@__DIR__,"..","anneal.jl"))
const W=MolecularAnnealingWorker

# This process owns a synthetic, more-specific worker method. It checks the
# entry boundary without loading a model, scoring an order or asserting time.
@eval W begin
    const ENTRY_CALLS=Any[]
    function run(request::String,result::String,cpu::Float64,wall::Float64)
        push!(ENTRY_CALLS,(request,result,cpu,wall))
        request=="raise" && throw(ArgumentError("injected entry refusal"))
        request=="status-two" ? 2 : 0
    end
end

@testset "Entry forwards parsed arguments, status and refusal handling" begin
    @test W.main(["request","result","20","116"])==0
    @test W.ENTRY_CALLS[end]==("request","result",20.0,116.0)
    @test W.main(["status-two","result","3.5","9"])==2
    @test W.ENTRY_CALLS[end]==("status-two","result",3.5,9.0)
    @test_throws ErrorException W.main(["incomplete"])
    mktempdir() do directory
        result=joinpath(directory,"checkpoint.toml")
        status=redirect_stderr(devnull) do
            W.main(["raise",result,"20","116"])
        end
        @test status==2
        @test !ispath(result)
        @test read(result*".refusal",String)=="admission_or_support_limit"
        write(result*".refusal","retained refusal")
        before=length(W.ENTRY_CALLS)
        status=redirect_stderr(devnull) do
            W.main(["request",result,"invalid CPU","116"])
        end
        @test status==2
        @test length(W.ENTRY_CALLS)==before
        @test read(result*".refusal",String)=="retained refusal"
    end
end
