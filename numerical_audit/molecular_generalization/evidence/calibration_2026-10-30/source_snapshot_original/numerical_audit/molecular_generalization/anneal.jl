module MolecularAnnealingWorker
using TOML, SHA
include(joinpath(@__DIR__,"src","FullOrderAnnealing.jl"))
using .FullOrderAnnealing

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const SOURCE_NAMES=("analysis/sparse_response_ordering.py","analysis/strang_schedule.py",
    "analysis/julia_window/ResponseCore.jl","analysis/julia_window/Project.toml",
    "numerical_audit/molecular_generalization/src/FullOrderAnnealing.jl",
    "numerical_audit/molecular_generalization/anneal.jl",
    "numerical_audit/molecular_generalization/annealing.py")
const REQUEST_FIELDS=("version","n","reference","time","budget","order","starts","word_lengths",
    "parent_score","parent_cubic_response","parent_steps","max_entries","max_cached_actions",
    "seed_material","seeds","settings_sha256","coefficients_sha256","source_sha256","settings","units")

require(ok,message)=ok || throw(ArgumentError(message))
hashspelling(x)=x isa String && occursin(r"^[0-9a-f]{64}$",x)
integer(x,name,lo,hi)=FullOrderAnnealing.integer(x,name,lo,hi)
cpu_seconds()=ccall(:clock,Int64,())/1e6
function rational_string(x)
    require(x isa String && ncodeunits(x)<=2500 && occursin(r"^-?[0-9]+/[0-9]+$",x),"bounded rational spelling")
    a,b=split(x,'/');p=parse(BigInt,a);q=parse(BigInt,b)
    require(q>0 && ndigits(abs(p);base=2)<=4096 && ndigits(q;base=2)<=4096,"rational bit cap")
    value=p//q
    require(x==string(numerator(value))*"/"*string(denominator(value)),"canonical rational spelling")
    value
end
function bounded_toml(path,maximum;max_lines=9000)
    require(isfile(path) && filesize(path)<=maximum,"request/reply byte cap")
    raw=open(path,"r") do io;read(io,maximum+1);end
    require(length(raw)<=maximum && all(byte->byte<0x80,raw),"bounded ASCII TOML")
    quotechar=UInt8(0);escaped=false;stack=UInt8[];lines=1
    for byte in raw
        if quotechar!=0
            require(!(byte in (0x0a,0x0d)),"multiline strings outside schema")
            if escaped;escaped=false
            elseif quotechar==0x22 && byte==0x5c;escaped=true
            elseif byte==quotechar;quotechar=0
            end
        elseif byte in (0x22,0x27)
            quotechar=byte
        elseif byte in (0x5b,0x7b)
            push!(stack,byte==0x5b ? 0x5d : 0x7d);require(length(stack)<=3,"TOML nesting cap")
        elseif byte in (0x5d,0x7d)
            require(!isempty(stack) && pop!(stack)==byte,"TOML bracket mismatch")
        elseif byte==0x23
            throw(ArgumentError("comments outside wire schema"))
        end
        byte==0x0a && (lines+=1);require(lines<=max_lines,"TOML line cap")
    end
    require(quotechar==0 && isempty(stack),"unterminated TOML")
    (data=TOML.parse(String(copy(raw))),sha256=bytes2hex(sha256(raw)))
end
function coefficient_hash(units)
    io=IOBuffer()
    for u in units
        print(io,u["flip"],':')
        for j in eachindex(u["parities"])
            j>1 && print(io,';')
            print(io,u["parities"][j],',',u["signs"][j],',',u["coefficients"][j])
        end
        print(io,'\n')
    end
    bytes2hex(sha256(take!(io)))
end
function read_request(path,check=FullOrderAnnealing.noop)
    check();parsed=bounded_toml(path,8*1024^2);data=parsed.data
    require(Set(keys(data))==Set(REQUEST_FIELDS),"request schema")
    require(data["version"]===1,"request version")
    n=integer(data["n"],"register",1,14)
    hf=integer(data["reference"],"HF",0,(1<<n)-1)
    require(rational_string(data["time"])==1//1,"T=1 training required")
    budget=integer(data["budget"],"training budget",0,typemax(Int)÷2)
    entries=integer(data["max_entries"],"max_entries",1,50000)
    cached=integer(data["max_cached_actions"],"max_cached_actions",0,100000)
    settings=validate_settings(data["settings"])
    require(hashspelling(data["settings_sha256"]) && hashspelling(data["coefficients_sha256"]),"request hashes")
    sources=data["source_sha256"]
    require(sources isa AbstractDict && Set(keys(sources))==Set(SOURCE_NAMES),"source inventory")
    for name in SOURCE_NAMES
        check();require(hashspelling(sources[name]),"source digest spelling")
        source=joinpath(ROOT,name)
        require(isfile(source) && filesize(source)<=1024^2,"runtime source byte cap")
        bytes=open(source,"r") do io;read(io,1024^2+1);end
        require(length(bytes)<=1024^2 && bytes2hex(sha256(bytes))==sources[name],"runtime source hash mismatch: "*name)
    end
    raw=data["units"]
    require(raw isa Vector && 1<=length(raw)<=1024,"unit admission cap")
    require(all(u->u isa AbstractDict && Set(keys(u))==Set(("flip","parities","signs","coefficients")) &&
        u["coefficients"] isa Vector && u["parities"] isa Vector && u["signs"] isa Vector,raw),"unit schema/array types")
    require(sum(length(u["coefficients"]) for u in raw)<=8192,"word admission cap")
    units=Unit[]
    for u in raw
        check();require(u isa AbstractDict && Set(keys(u))==Set(("flip","parities","signs","coefficients")),"unit schema")
        flip=integer(u["flip"],"flip",0,(1<<n)-1)
        p=u["parities"];s=u["signs"];a=u["coefficients"]
        require(p isa Vector && s isa Vector && a isa Vector && 1<=length(p)==length(s)==length(a)<=8192,"unit arrays")
        parity=UInt64[];signs=Int8[];exact=Rational{BigInt}[];coeff=Float64[]
        for (pp,ss,aa) in zip(p,s,a)
            check();push!(parity,UInt64(integer(pp,"parity",0,(1<<n)-1)))
            require(ss isa Integer && !(ss isa Bool) && ss in (-1,1),"real sign")
            push!(signs,Int8(ss));q=rational_string(aa);f=Float64(q)
            require(!iszero(q) && isfinite(f) && !iszero(f),"coefficient conversion outside resolution")
            push!(exact,q);push!(coeff,f)
        end
        push!(units,Unit(UInt64(flip),parity,signs,coeff;exact=exact))
    end
    require(coefficient_hash(raw)==data["coefficients_sha256"],"coefficient metadata hash mismatch")
    m=length(units);order=data["order"]
    valid_order(x)=x isa Vector && length(x)==m && all(i->i isa Integer && !(i isa Bool),x) && sort(x)==collect(0:m-1)
    require(valid_order(order),"complete zero-based parent order")
    starts=data["starts"]
    require(starts isa Vector && length(starts)==4 && all(valid_order,starts) && starts[1]==order &&
        length(Set(Tuple(p) for p in starts))==4,"four distinct parent-based starts")
    require(data["word_lengths"]==[length(u.coeff) for u in units] && all(x->x isa Integer && !(x isa Bool),data["word_lengths"]),"original lengths")
    require(budget==30sum(data["word_lengths"]),"original 30U budget")
    score=data["parent_score"];response=data["parent_cubic_response"]
    require(all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>=0,(score,response)),"finite parent witness")
    steps=integer(data["parent_steps"],"parent steps",1,typemax(Int)÷2)
    require(steps==FullOrderAnnealing.ResponseCore.affordable_steps(units,Int.(order).+1,budget),"parent step metadata")
    require(abs(score-response/Float64(steps)^4)<=max(1e-30,abs(score)*1e-12),"parent normalization metadata")
    material=data["seed_material"];seeds=FullOrderAnnealing._seeds(material)
    require(data["seeds"]==string.(seeds),"seed metadata mismatch")
    (data=data,units=units,n=n,hf=hf,budget=budget,order=Int.(order).+1,
        starts=[Int.(p).+1 for p in starts],settings=settings,max_entries=entries,
        max_cached_actions=cached,request_sha256=parsed.sha256)
end

function atomic_checkpoint(path,request,state,phases)
    require(state["steps"] isa Integer && !(state["steps"] isa Bool) && state["steps"]>0 &&
        state["order"] isa Vector && all(x->x isa Integer && !(x isa Bool),state["order"]) &&
        sort(state["order"])==collect(1:length(request.units)),"no partial/order-invalid checkpoint")
    require(all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>=0,
        (state["score"],state["cubic_response"])) &&
        !(state["cubic_response"]>0 && iszero(state["score"])) &&
        state["steps"]==FullOrderAnnealing.ResponseCore.affordable_steps(request.units,state["order"],request.budget) &&
        abs(state["score"]-state["cubic_response"]/Float64(state["steps"])^4)<=max(1e-30,abs(state["score"])*1e-12),
        "no incomplete/nonfinite/unnormalized score checkpoint")
    payload=deepcopy(state);payload["order"].-=1
    merge!(payload,Dict("version"=>1,"request_sha256"=>request.request_sha256,
        "settings_sha256"=>request.data["settings_sha256"],"coefficients_sha256"=>request.data["coefficients_sha256"],
        "source_sha256"=>request.data["source_sha256"],"word_lengths"=>request.data["word_lengths"],
        "time"=>1.0,"budget"=>request.budget,"julia_version"=>string(VERSION),
        "rng"=>"Random.Xoshiro","score_scope"=>"full_order_transported_projected_cubic_J_not_error_bound",
        "phases"=>copy(phases)))
    io=IOBuffer();TOML.print(io,payload;sorted=true);raw=take!(io)
    require(length(raw)<=64*1024,"checkpoint byte cap")
    temporary=path*".tmp"
    # One worker owns this path. Only the complete current witness is replaced.
    open(temporary,"w") do out;write(out,raw);end
    mv(temporary,path;force=true)
end
function run(request_path,result_path,cpu_cap,wall_cap)
    require(isfinite(cpu_cap) && cpu_cap>0 && isfinite(wall_cap) && wall_cap>0,"positive finite process caps")
    require(!ispath(result_path) && !ispath(result_path*".tmp") && !ispath(result_path*".refusal"),"fresh checkpoint required")
    Threads.nthreads()==1 || throw(ArgumentError("one Julia numerical thread required"))
    began=cpu_seconds();wall=time_ns();actions=Ref(0);action_limit=Ref(100_000_000)
    cleanup=Ref(0.05)
    check=(kind=:loop,amount=1)->begin
        cpu_seconds()<cpu_cap-cleanup[] || throw(Limit("cpu_limit","owned Julia CPU allowance"))
        (time_ns()-wall)/1e9<wall_cap || throw(Limit("wall_limit","owned Julia wall allowance"))
        if kind==:unit_action
            actions[]+=amount
            actions[]<=action_limit[] || throw(Limit("action_limit","sparse unit action ceiling"))
        end
    end
    request=read_request(request_path,check)
    action_limit[]=request.settings["max_actions"];cleanup[]=request.settings["cleanup_CPU_seconds"]
    phases=Dict("request_CPU_seconds"=>cpu_seconds()-began)
    phase=cpu_seconds()
    context=FullOrderContext(request.units,request.n,request.hf,request.budget;
        max_entries=request.max_entries,max_cached_actions=request.max_cached_actions,
        cache_bytes=request.settings["cache_bytes"],check=check)
    phases["context_CPU_seconds"]=cpu_seconds()-phase
    function checkpoint(state)
        # The first scored start is the supplied parent. Authenticate it before
        # publishing any Julia witness, without paying for a duplicate score.
        if state["start_scores_completed"]==1
            require(abs(state["score"]-request.data["parent_score"])<=max(1e-28,abs(state["score"])*1e-8) &&
                abs(state["cubic_response"]-request.data["parent_cubic_response"])<=max(1e-28,abs(state["cubic_response"])*1e-8),
                "parent native/Julia score mismatch")
        end
        atomic_checkpoint(result_path,request,state,phases)
    end
    state=anneal(context,request.starts,request.settings;seed_material=request.data["seed_material"],
        cpu_deadline=cpu_cap-cleanup[],parent_value=request.data["parent_score"],check=check,checkpoint=checkpoint)
    state["status"] in ("complete","cpu_limit","wall_limit","evaluation_limit","proposal_limit","action_limit") ? 0 : 2
end
function main(args=ARGS)
    length(args)==4 || error("usage: anneal.jl REQUEST.toml FRESH_CHECKPOINT.toml CPU_CAP WALL_CAP")
    try
        run(args[1],args[2],parse(Float64,args[3]),parse(Float64,args[4]))
    catch e
        reason=e isa Limit ? e.status : e isa FullOrderAnnealing.NumericalResolution ? "numerical_limit" :
            e isa ArgumentError ? "admission_or_support_limit" : "worker_failure"
        if !ispath(args[2]*".refusal")
            open(args[2]*".refusal","w") do io;write(io,reason);end
        end
        showerror(stderr,e);println(stderr)
        2
    end
end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(MolecularAnnealingWorker.main())
end
