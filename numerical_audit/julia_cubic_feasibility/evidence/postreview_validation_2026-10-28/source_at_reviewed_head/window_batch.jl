include("window_study.jl")
length(ARGS)>=3 && iseven(length(ARGS)-1) || error("usage: window_batch.jl MODE INPUT OUTPUT [INPUT OUTPUT...]")
for i in 2:2:length(ARGS)
    main([ARGS[1],ARGS[i],ARGS[i+1]])
end
