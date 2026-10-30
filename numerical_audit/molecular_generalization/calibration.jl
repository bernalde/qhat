module OldAnnealingCalibration

using TOML, SHA
include(joinpath(@__DIR__,"anneal.jl"))
const W=MolecularAnnealingWorker
const A=W.FullOrderAnnealing
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"numerical_audit","prospective_window_ordering","objective_check.jl"))

const SOURCE_NAMES=("analysis/sparse_response_ordering.py","analysis/strang_schedule.py",
    "analysis/julia_window/ResponseCore.jl","analysis/julia_window/Project.toml",
    "numerical_audit/molecular_generalization/src/FullOrderAnnealing.jl",
    "numerical_audit/molecular_generalization/anneal.jl",
    "numerical_audit/molecular_generalization/annealing.py",
    "numerical_audit/molecular_generalization/calibration.py",
    "numerical_audit/molecular_generalization/calibration.jl",
    "numerical_audit/prospective_window_ordering/objective_check.jl",
    "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
    "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl")
const FIXED_MODELS=Dict(
    "CO_r1p13"=>("b78921326c04f331a0d95b8fadf1be600bb6bfd0caa65e464d24c7aa6de5fe52",
        "8a9f89166544f72f9ecae892ad664a3f110f1a2c534fceef22359aaa1fcbdbe8",195,770),
    "CH4_r1p09"=>("00a283f3736e331b2c3bfc140213e1d5d1d3a62f119b34ab1450a1b091afbc97",
        "32aeb49b2e541088efa0604e90dae046e8be58e81259e027032b65fe2d779e37",351,1722))
const FIELDS=("version","case","model_sha256","legacy_source_model_sha256","n","reference","words",
    "units","time","budget","probe_count","seed_material","seed","pairs","orders","source_sha256","native_scores")
require(ok,message)=ok || throw(ArgumentError(message))
cpu_seconds()=ccall(:clock,Int64,())/1e6

function write_new(path,record)
    require(!ispath(path) && !ispath(path*".partial"),"fresh calibration artifact")
    io=IOBuffer();TOML.print(io,ObjectiveCheck.WindowStudyIO.clean(record);sorted=true);bytes=take!(io)
    require(length(bytes)<=64*1024,"calibration artifact byte cap")
    open(path*".partial","w") do stream;write(stream,bytes);end
    mv(path*".partial",path;force=false)
end
function read_request(modelpath,requestpath,check)
    check();request=W.bounded_toml(requestpath,8*1024^2)
    raw=W.bounded_toml(modelpath,8*1024^2)
    q=request.data
    require(Set(keys(q))==Set(FIELDS) && q["version"]===1,"calibration request schema")
    case=q["case"];require(case isa String && haskey(FIXED_MODELS,case),"only fixed old calibration cases")
    expected=FIXED_MODELS[case]
    require(q["model_sha256"]==raw.sha256==expected[1] && q["legacy_source_model_sha256"]==expected[2] &&
        raw.data["case"]==case && raw.data["source_model_sha256"]==expected[2] &&
        !haskey(raw.data,"source_json_sha256"),"fixed legacy input/provenance bytes")
    require(q["n"]===12 && q["reference"]===63 && q["units"]===expected[3] && q["words"]===expected[4] &&
        q["time"]===1.0 && q["budget"]===30*expected[4] && q["probe_count"]===32,"fixed old model/training inventory")
    sources=q["source_sha256"]
    require(sources isa AbstractDict && Set(keys(sources))==Set(SOURCE_NAMES),"calibration source inventory")
    for name in SOURCE_NAMES
        check();path=joinpath(ROOT,name)
        require(isfile(path) && filesize(path)<=1024^2 && W.hashspelling(sources[name]),"source byte/digest cap")
        bytes=open(path,"r") do stream;read(stream,1024^2+1);end
        require(length(bytes)<=1024^2 && bytes2hex(sha256(bytes))==sources[name],"calibration source changed: "*name)
    end
    orders=q["orders"];pairs=q["pairs"];m=expected[3]
    validorder(order)=order isa Vector && length(order)==m && all(x->x isa Integer && !(x isa Bool),order) && sort(order)==collect(0:m-1)
    require(orders isa Vector && length(orders)==33 && all(validorder,orders),"33 complete full orders")
    require(pairs isa Vector && length(pairs)==32,"32 declared insertion pairs")
    require(raw.data["cubic_order"]==orders[1].+1,"retained parent order identity")
    for index in 1:32
        check();pair=pairs[index]
        require(pair isa Vector && length(pair)==2 && all(x->x isa Integer && !(x isa Bool) && 0<=x<m,pair) && pair[1]!=pair[2],"full insertion pair")
        require(A.insertion(orders[1],pair[1]+1,pair[2]+1)==orders[index+1],"probe changed its declared parent insertion")
    end
    require(q["seed_material"]==case*":molecular-generalization-sa-v1:calibration-v1" &&
        q["seed"] isa String && occursin(r"^[0-9]{1,20}$",q["seed"]),"predeclared seed identity")
    computed=foldl((a,b)->(a<<8)|UInt64(b),sha256(q["seed_material"])[1:8];init=UInt64(0))
    require(q["seed"]==string(computed),"seed digest mismatch")
    native=q["native_scores"]
    require(native isa Vector && length(native)==33,"33 native claims")
    for score in native
        require(score isa AbstractDict && Set(keys(score))==Set(("value","cubic_response","steps")) &&
            all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>=0,(score["value"],score["cubic_response"])) &&
            score["steps"] isa Integer && !(score["steps"] isa Bool) && score["steps"]>0,"finite complete native witness")
        require(!(score["cubic_response"]>0 && iszero(score["value"])) &&
            abs(score["value"]-score["cubic_response"]/Float64(score["steps"])^4)<=max(1e-30,abs(score["value"])*1e-12),"native normalization/underflow")
    end
    data=ObjectiveCheck.WindowStudyIO.load_model(modelpath)
    require(data.input_sha256==raw.sha256 && data.n==12 && data.electrons==6 && data.hf==63 &&
        data.words==expected[4] && length(data.units)==expected[3],"loaded legacy identity")
    denominator=parse(BigInt,data.source["denominator"])
    units=A.Unit[]
    for (unit,original) in zip(data.units,data.source["units"])
        check();exact=[parse(BigInt,a)//denominator for a in original["numerators"]]
        push!(units,A.Unit(unit.flip,unit.parity,unit.signs,unit.coeff;exact=exact))
    end
    budget=ObjectiveCheck.ReplayBudget(()->check())
    archived=ObjectiveCheck.prepare_model(data;budget=budget)
    (q=q,data=data,units=units,archived=archived,archive_budget=budget,request_sha256=request.sha256)
end

score_dict(score)=Dict{String,Any}("value"=>score.value,"cubic_response"=>score.cubic_response,"steps"=>score.steps)
function parity(row)
    checks=Dict{String,Any}();verified=true
    for (name,left,right) in (("native_full",row["native_claim"],row["full_order_score"]),
                              ("native_archived",row["native_claim"],row["archived_score"]),
                              ("full_archived",row["full_order_score"],row["archived_score"]))
        for key in ("value","cubic_response")
            result=ObjectiveCheck.compare_scores(left[key],right[key])
            checks[name*"_"*key]=result;verified&=result.verified
        end
        verified&=left["steps"]==right["steps"]
    end
    (verified=verified,comparisons=checks)
end
function run(modelpath,requestpath,output,cpu_cap,wall_cap)
    require(isfinite(cpu_cap) && cpu_cap>.1 && isfinite(wall_cap) && wall_cap>1,"finite owned calibration allowance")
    require(!ispath(output),"fresh Julia calibration output directory");mkpath(output)
    Threads.nthreads()==1 || throw(ArgumentError("one Julia thread"))
    wall=time_ns();actions=Ref(0)
    check=(kind=:loop,amount=1)->begin
        cpu_seconds()<cpu_cap-.1 || throw(A.Limit("cpu_limit","calibration Julia CPU allowance"))
        (time_ns()-wall)/1e9<wall_cap-1 || throw(A.Limit("wall_limit","calibration Julia wall allowance"))
        if kind==:unit_action
            actions[]+=amount;actions[]<=100_000_000 || throw(A.Limit("action_limit","calibration sparse action cap"))
        end
    end
    request=read_request(modelpath,requestpath,check);q=request.q
    header=Dict{String,Any}("version"=>1,"case"=>q["case"],"model_sha256"=>q["model_sha256"],
        "request_sha256"=>request.request_sha256,"source_sha256"=>q["source_sha256"],
        "julia_version"=>string(VERSION),"complete"=>false,"processed_rows"=>0,
        "physical_errors_read"=>false,"score_scope"=>"native_full_order_and_archived_Julia_cubic_parity_T1_30U")
    write_new(joinpath(output,"START.toml"),header)
    context=A.FullOrderContext(request.units,request.data.n,request.data.hf,q["budget"];check=check)
    verified=true;stopped=false
    for index in 0:32
        row=Dict{String,Any}("version"=>1,"case"=>q["case"],"index"=>index,
            "order"=>q["orders"][index+1],"model_sha256"=>q["model_sha256"],
            "legacy_source_model_sha256"=>q["legacy_source_model_sha256"],"request_sha256"=>request.request_sha256,
            "source_sha256"=>q["source_sha256"],"native_claim"=>q["native_scores"][index+1],
            "julia_version"=>string(VERSION),"status"=>"started")
        try
            check();order=Int.(row["order"]).+1
            full=A.full_score!(context,order;check=check);row["full_order_score"]=score_dict(full)
            archived=ObjectiveCheck.replay_objective(request.archived,"conserving",order,
                "diagonal_reference_cubic_response";budget=request.archive_budget)
            row["archived_score"]=Dict{String,Any}("value"=>archived.value,"cubic_response"=>archived.details["cubic_response"],"steps"=>archived.steps)
            checked=parity(row);row["parity"]=checked.comparisons;row["verified"]=checked.verified
            row["status"]="complete";verified&=checked.verified
        catch error
            verified=false;row["verified"]=false;row["status"]="refused"
            row["refusal"]=Dict("type"=>string(typeof(error)),"message"=>sprint(showerror,error))
            stopped=error isa A.Limit || error isa ObjectiveCheck.RunBudgetRefusal
        end
        write_new(joinpath(output,"row"*lpad(string(index),3,'0')*".toml"),row)
        header["processed_rows"]=index+1
        stopped && break
    end
    header["complete"]=header["processed_rows"]==33 && verified && !stopped
    header["CPU_before_summary"]=cpu_seconds();header["wall_before_summary"]=(time_ns()-wall)/1e9
    write_new(joinpath(output,"SUMMARY.toml"),header)
    header["complete"] ? 0 : 2
end
function main(args=ARGS)
    length(args)==5 || error("usage: calibration.jl MODEL.toml REQUEST.toml FRESH_OUTPUT CPU_CAP WALL_CAP")
    try
        run(args[1],args[2],args[3],parse(Float64,args[4]),parse(Float64,args[5]))
    catch error
        # Partial row files remain authoritative. No retry or input substitution.
        if isdir(args[3]) && !ispath(joinpath(args[3],"REFUSAL.toml"))
            write_new(joinpath(args[3],"REFUSAL.toml"),Dict("type"=>string(typeof(error)),"message"=>sprint(showerror,error)))
        end
        showerror(stderr,error);println(stderr);2
    end
end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(OldAnnealingCalibration.main())
end
