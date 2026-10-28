include("window_study.jl")

function batch_main(args;run_case=main)
    length(args)>=3 && iseven(length(args)-1) || error("usage: window_batch.jl MODE INPUT OUTPUT [INPUT OUTPUT...]")
    records=Dict{String,Any}[]
    for i in 2:2:length(args)
        push!(records,run_case([args[1],args[i],args[i+1]];exit_on_failure=false))
    end
    records
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    records=batch_main(ARGS)
    all(record["execution_status"]=="complete" for record in records) || exit(2)
end
