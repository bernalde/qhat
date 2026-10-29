# QHAT window protocol v1: one stdlib-only process, no package installation.
using TOML, SHA
if !isdefined(@__MODULE__, :ResponseCore)
    include(joinpath(@__DIR__, "ResponseCore.jl"))
end
if !isdefined(@__MODULE__, :WindowSearch)
    include(joinpath(@__DIR__, "WindowSearch.jl"))
end
using .ResponseCore, .WindowSearch

cpu_seconds() = ccall(:clock, Int64, ()) / 1.0e6
struct BudgetRefusal <: Exception
    reason::String
end
Base.showerror(io::IO, e::BudgetRefusal) = print(io, e.reason)

function bounded_integer(data, key, low, high)
    value = data[key]
    typeof(value) === Int64 && low <= value <= high || error("invalid integer $key")
    Int(value)
end
function rational_string(value; positive=false)
    value isa String && sizeof(value) <= 2500 || error("oversized rational spelling")
    occursin(r"^-?[0-9]+/[0-9]+$", value) || error("invalid rational spelling")
    a,b = split(value,'/')
    numerator, denominator = parse(BigInt,a), parse(BigInt,b)
    denominator > 0 && ndigits(abs(numerator);base=2) <= 4096 && ndigits(denominator;base=2) <= 4096 || error("rational bit cap")
    result=numerator//denominator
    positive && result<=0 && error("positive rational required")
    result
end
function schema(data, fields)
    Set(keys(data)) == Set(fields) || error("unexpected protocol fields")
end

function read_request(path, check)
    check()
    filesize(path) <= 8*1024^2 || error("request byte cap")
    raw=read(path); check()
    request_hash=bytes2hex(sha256(raw))
    data=TOML.parse(String(raw)); check()
    schema(data,["version","n","reference","time","budget","order","word_lengths","parent_score",
                 "width","max_windows","max_sweeps","max_entries","max_cached_actions",
                 "max_factor_entries","max_flat_entries","max_actions","max_atoms","max_evaluations","units"])
    typeof(data["version"]) === Int64 && data["version"] == 1 || error("unsupported protocol version")
    n=bounded_integer(data,"n",1,20)
    reference=bounded_integer(data,"reference",0,(1<<n)-1)
    duration=rational_string(data["time"];positive=true)
    time=Float64(duration)
    isfinite(time) && time>0 || error("time conversion outside resolution")
    budget=bounded_integer(data,"budget",0,typemax(Int64)÷2)
    width=bounded_integer(data,"width",2,9)
    max_windows=bounded_integer(data,"max_windows",1,32)
    max_sweeps=bounded_integer(data,"max_sweeps",1,100)
    max_entries=bounded_integer(data,"max_entries",1,50000)
    max_cached_actions=bounded_integer(data,"max_cached_actions",0,100000)
    max_factor_entries=bounded_integer(data,"max_factor_entries",1,250000)
    max_flat_entries=bounded_integer(data,"max_flat_entries",1,1000000)
    max_actions=bounded_integer(data,"max_actions",1,10000000)
    max_atoms=bounded_integer(data,"max_atoms",1,1000000)
    max_evaluations=bounded_integer(data,"max_evaluations",1,1000000)
    blocks=data["units"]
    blocks isa Vector && width+2 <= length(blocks) <= 1024 || error("unit admission cap or no interior window")
    units=Unit[]; words=0
    for block in blocks
        check()
        schema(block,["flip","parities","signs","coefficients"])
        flip=bounded_integer(block,"flip",0,(1<<n)-1)
        pars=block["parities"]; signs=block["signs"]; coefficients=block["coefficients"]
        pars isa Vector && signs isa Vector && coefficients isa Vector || error("invalid unit arrays")
        1 <= length(pars) == length(signs) == length(coefficients) <= 8192 || error("unit word count")
        words += length(pars); words <= 8192 || error("total word cap")
        parity=UInt64[]; phases=Int8[]; coeff=Float64[]; exact=Rational{BigInt}[]
        for (p,s,c) in zip(pars,signs,coefficients)
            check()
            typeof(p) === Int64 && 0<=p<(1<<n) || error("parity outside register")
            typeof(s) === Int64 && s in (-1,1) || error("invalid real phase")
            ys=count_ones(UInt64(flip)&UInt64(p))
            iseven(ys) && s==(isodd(ys÷2) ? -1 : 1) || error("phase is not a real Pauli word")
            rational=rational_string(c)
            value=Float64(rational)
            isfinite(value) && value!=0 || error("coefficient conversion outside resolution")
            # Equal-flip words commute iff x.(z_i xor z_j) is even. A common
            # parity relative to the first word avoids a quadratic word scan.
            isempty(parity) || iseven(count_ones(UInt64(flip)&(UInt64(p)⊻parity[1]))) || error("noncommuting unit")
            push!(parity,UInt64(p));push!(phases,Int8(s));push!(coeff,value);push!(exact,rational)
        end
        push!(units,Unit(UInt64(flip),parity,phases,coeff;exact=exact))
    end
    m=length(units)
    order=data["order"]
    order isa Vector && all(i->typeof(i)===Int64,order) && sort(order)==collect(0:m-1) || error("invalid global order")
    lengths=data["word_lengths"]
    lengths isa Vector && all(i->typeof(i)===Int64,lengths) &&
        lengths==[length(u.coeff) for u in units] || error("original word lengths differ")
    parent_score=data["parent_score"]
    parent_score isa AbstractFloat && isfinite(parent_score) && parent_score>=0 || error("invalid parent score")
    (data=data,units=units,order=Int.(order).+1,hf=UInt64(reference),time=time,budget=budget,width=width,
     max_windows=max_windows,max_sweeps=max_sweeps,max_entries=max_entries,max_cached_actions=max_cached_actions,
     max_factor_entries=max_factor_entries,max_flat_entries=max_flat_entries,max_actions=max_actions,
     max_atoms=max_atoms,max_evaluations=max_evaluations,request_hash=request_hash)
end

function select_window(request, order, used)
    # Declared deterministic heuristic: maximum count of nonzero pair atoms.
    # It uses supplied units/reference only, never a physical-error oracle.
    centered=center_units(request.units,request.hf)
    # Production selection counts exact preprojection pair support, including
    # HF; the one rational-only store shares its byte/count cap with kernels.
    cache=AmplitudeCache(centered;max_entries=request.max_cached_actions,guard=true)
    best=nothing; best_count=-1
    for lo in 2:length(order)-request.width
        check_budget()
        lo in used && continue
        count=0
        for a in lo:lo+request.width-2, b in a+1:lo+request.width-1
            check_budget()
            atom=ResponseCore._atom_column(cache,(order[a],order[b]),request.hf,request.max_entries)
            count += !isempty(atom)
        end
        if count>best_count
            best=lo:lo+request.width-1;best_count=count
        end
    end
    best
end

function write_checkpoint(path, request, order, score, status, counters, phases, updates, windows, history)
    check_budget()
    payload=Dict("version"=>1,"request_sha256"=>request.request_hash,"order"=>order.-1,
        "word_lengths"=>[length(u.coeff) for u in request.units],"score"=>score.value,
        "cubic_response"=>score.cubic_response,"steps"=>score.steps,"status"=>status,
        "checked_updates"=>updates,"windows"=>windows,"history"=>history,
        "counters"=>copy(counters),"phases"=>copy(phases))
    temporary=path*".tmp"
    open(temporary,"w") do io; TOML.print(io,payload;sorted=true); end
    # Only a completely checked global incumbent is ever atomically published.
    mv(temporary,path;force=true)
end

function run(request_path, result_path, cpu_cap, wall_cap)
    started=cpu_seconds(); wall_started=time_ns()
    counters=Dict("evaluation"=>0,"action"=>0,"atom"=>0,"full_score"=>0,
                  "flat_candidate"=>0,"exact_witness"=>0,"insertion_move"=>0)
    phases=Dict{String,Float64}()
    limits=Dict("evaluation"=>1000000,"action"=>10000000,"atom"=>1000000,
                "full_score"=>33,"flat_candidate"=>1000000,"exact_witness"=>1000000,"insertion_move"=>3200)
    function check(kind=:loop,amount=1)
        cpu_seconds()-started < cpu_cap || throw(BudgetRefusal("cpu_limit"))
        (time_ns()-wall_started)/1e9 < wall_cap || throw(BudgetRefusal("wall_limit"))
        key=String(kind)
        if haskey(counters,key)
            counters[key]+=amount
            counters[key]<=limits[key] || throw(BudgetRefusal(key*"_limit"))
        end
        nothing
    end
    with_budget(check) do
        request=read_request(request_path,check)
        limits["evaluation"]=request.max_evaluations
        limits["action"]=request.max_actions;limits["atom"]=request.max_atoms
        phases["request_CPU_seconds"]=cpu_seconds()-started
        order=copy(request.order); updates=0; windows=Vector{Vector{Int}}(); used=Set{Int}()
        history=Vector{Dict{String,Any}}()
        phase=cpu_seconds();check_budget(:evaluation)
        score=score_order(request.units,order,request.hf,request.budget;time=request.time,max_entries=request.max_entries)
        check_budget(:full_score)
        phases["parent_full_score_CPU_seconds"]=cpu_seconds()-phase
        # Compatibility check, never use a divergent native parent as a baseline.
        abs(score.value-request.data["parent_score"]) <= max(1e-28,abs(score.value)*1e-8) || error("parent native/Julia score mismatch")
        write_checkpoint(result_path,request,order,score,"parent_checked",counters,phases,updates,windows,history)
        for _ in 1:request.max_windows
            check_budget()
            phase=cpu_seconds(); positions=select_window(request,order,used)
            phases["selection_CPU_seconds"]=get(phases,"selection_CPU_seconds",0.0)+cpu_seconds()-phase
            isnothing(positions) && break
            push!(used,first(positions))
            phase=cpu_seconds()
            model=build_window_model(request.units,order,request.hf,positions;time=request.time,
                max_entries=request.max_entries,max_cached_actions=request.max_cached_actions,
                max_factor_entries=request.max_factor_entries)
            phases["setup_CPU_seconds"]=get(phases,"setup_CPU_seconds",0.0)+cpu_seconds()-phase
            phase=cpu_seconds();flat=prepare_window(model;max_flat_entries=request.max_flat_entries)
            phases["flatten_CPU_seconds"]=get(phases,"flatten_CPU_seconds",0.0)+cpu_seconds()-phase
            phase=cpu_seconds();proposal=insertion_window(flat;guard=check,max_sweeps=request.max_sweeps)
            phases["insertion_CPU_seconds"]=get(phases,"insertion_CPU_seconds",0.0)+cpu_seconds()-phase
            candidate=copy(order);candidate[positions]=proposal.order
            record=Dict{String,Any}("positions"=>[first(positions)-1,last(positions)-1],
                "baseline_labels"=>order[positions].-1,"proposal_labels"=>proposal.order.-1,"accepted"=>false)
            if candidate==order
                push!(history,record)
                continue
            end
            phase=cpu_seconds();check_budget(:evaluation)
            actual=score_order(request.units,candidate,request.hf,request.budget;time=request.time,max_entries=request.max_entries)
            check_budget(:full_score)
            phases["verification_CPU_seconds"]=get(phases,"verification_CPU_seconds",0.0)+cpu_seconds()-phase
            check_budget()
            actual.steps==score.steps || error("interior window changed affordable steps")
            if actual.value < score.value-max(1e-30,abs(score.value)*1e-12)
                order=candidate;score=actual;updates+=1
                record["accepted"]=true;push!(history,record)
                push!(windows,[first(positions)-1,last(positions)-1])
                write_checkpoint(result_path,request,order,score,"checked_update",counters,phases,updates,windows,history)
            else
                push!(history,record)
            end
        end
        write_checkpoint(result_path,request,order,score,"complete",counters,phases,updates,windows,history)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==4 || error("expected request, checkpoint, CPU cap, wall cap")
    try
        run(ARGS[1],ARGS[2],parse(Float64,ARGS[3]),parse(Float64,ARGS[4]))
    catch error
        # Leave the last atomic checked checkpoint in place on every refusal.
        reason=error isa BudgetRefusal ? error.reason :
            error isa NumericalResolution ? "numerical_limit" :
            error isa ArgumentError ? "admission_or_support_limit" : "worker_failure"
        # Small diagnostic only; the Python owner still authenticates/rescores
        # the checkpoint before accepting any proposal.
        open(ARGS[2]*".refusal","w") do io;write(io,reason);end
        showerror(stderr,error);println(stderr)
        exit(2)
    end
end
