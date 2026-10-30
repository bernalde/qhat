module CalibrationResume
using TOML,SHA
include(joinpath(@__DIR__,"calibration.jl"))
const O=OldAnnealingCalibration
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
require(ok,message)=ok || throw(ArgumentError(message))

function run(model,request_path,subset_path,output,cpu,wall_limit)
    require(isfinite(cpu) && cpu>.1 && isfinite(wall_limit) && wall_limit>1,"finite subset allowance")
    require(!ispath(output),"fresh subset output");mkpath(output)
    Threads.nthreads()==1 || throw(ArgumentError("one Julia numerical thread"))
    wall=time_ns();actions=Ref(0)
    check=(kind=:loop,amount=1)->begin
        O.cpu_seconds()<cpu-.1 || throw(O.A.Limit("cpu_limit","subset CPU allowance"))
        (time_ns()-wall)/1e9<wall_limit-1 || throw(O.A.Limit("wall_limit","subset wall allowance"))
        if kind==:unit_action
            actions[]+=amount;actions[]<=100_000_000 || throw(O.A.Limit("action_limit","original sparse action cap"))
        end
    end
    subset=O.W.bounded_toml(subset_path,64*1024);s=subset.data
    require(Set(keys(s))==Set(("version","case","plan_sha256","producer_commit","model_sha256","request_sha256",
        "indices","julia_version","resume_source_sha256")) && s["version"]===1 && s["case"]=="CH4_r1p09","subset schema")
    indices=s["indices"]
    require(indices isa Vector && 1<=length(indices)<=33 && all(i->i isa Integer && !(i isa Bool) && 0<=i<=32,indices) &&
        indices==sort(unique(indices)),"fixed distinct sorted original row indices")
    require(O.W.hashspelling(s["plan_sha256"]) && occursin(r"^[0-9a-f]{40}$",s["producer_commit"]) &&
        s["julia_version"]==string(VERSION),"resume producer/runtime identity")
    sources=s["resume_source_sha256"]
    names=("numerical_audit/molecular_generalization/calibration_resume.py","numerical_audit/molecular_generalization/calibration_resume.jl")
    require(sources isa AbstractDict && Set(keys(sources))==Set(names),"resume source inventory")
    for name in names
        check();path=joinpath(ROOT,name);require(isfile(path) && filesize(path)<=1024^2,"resume source cap")
        bytes=open(path,"r") do stream;read(stream,1024^2+1);end
        require(length(bytes)<=1024^2 && bytes2hex(sha256(bytes))==sources[name],"resume source changed")
    end
    request=O.read_request(model,request_path,check);q=request.q
    require(q["case"]==s["case"] && q["model_sha256"]==s["model_sha256"] && request.request_sha256==s["request_sha256"],
        "unchanged original full33 model/request")
    header=Dict{String,Any}("version"=>1,"case"=>s["case"],"model_sha256"=>q["model_sha256"],
        "request_sha256"=>request.request_sha256,"source_sha256"=>q["source_sha256"],
        "resume_plan_sha256"=>s["plan_sha256"],"resume_subset_sha256"=>subset.sha256,
        "resume_producer_commit"=>s["producer_commit"],"resume_source_sha256"=>sources,
        "julia_version"=>string(VERSION),"indices"=>indices,"processed_indices"=>Int[],"complete"=>false,
        "physical_errors_read"=>false)
    O.write_new(joinpath(output,"START.toml"),header)
    context=O.A.FullOrderContext(request.units,request.data.n,request.data.hf,q["budget"];check=check)
    verified=true
    for index in indices
        row=Dict{String,Any}("version"=>1,"case"=>q["case"],"index"=>index,"order"=>q["orders"][index+1],
            "model_sha256"=>q["model_sha256"],"legacy_source_model_sha256"=>q["legacy_source_model_sha256"],
            "request_sha256"=>request.request_sha256,"source_sha256"=>q["source_sha256"],
            "native_claim"=>q["native_scores"][index+1],"julia_version"=>string(VERSION),
            "resume_plan_sha256"=>s["plan_sha256"],"resume_subset_sha256"=>subset.sha256,
            "resume_producer_commit"=>s["producer_commit"],"resume_source_sha256"=>sources,"status"=>"started")
        try
            check();order=Int.(row["order"]).+1
            row["full_order_score"]=O.score_dict(O.A.full_score!(context,order;check=check))
            archived=O.ObjectiveCheck.replay_objective(request.archived,"conserving",order,
                "diagonal_reference_cubic_response";budget=request.archive_budget)
            row["archived_score"]=Dict{String,Any}("value"=>archived.value,"cubic_response"=>archived.details["cubic_response"],"steps"=>archived.steps)
            parity=O.parity(row);row["parity"]=parity.comparisons;row["verified"]=parity.verified;row["status"]="complete"
            verified&=parity.verified
        catch error
            verified=false;row["status"]="refused";row["verified"]=false
            row["refusal"]=Dict("type"=>string(typeof(error)),"message"=>sprint(showerror,error))
        end
        O.write_new(joinpath(output,"row"*lpad(string(index),3,'0')*".toml"),row)
        push!(header["processed_indices"],index)
        verified || break # A floor/mismatch/refusal is retained, never retried.
    end
    header["complete"]=verified && header["processed_indices"]==indices
    header["CPU_before_summary"]=O.cpu_seconds();header["wall_before_summary"]=(time_ns()-wall)/1e9
    O.write_new(joinpath(output,"SUMMARY.toml"),header)
    header["complete"] ? 0 : 2
end
function main(args=ARGS)
    length(args)==6 || error("usage: calibration_resume.jl MODEL REQUEST SUBSET FRESH_OUTPUT CPU WALL")
    try
        run(args[1],args[2],args[3],args[4],parse(Float64,args[5]),parse(Float64,args[6]))
    catch error
        if isdir(args[4]) && !ispath(joinpath(args[4],"REFUSAL.toml"))
            O.write_new(joinpath(args[4],"REFUSAL.toml"),Dict("type"=>string(typeof(error)),"message"=>sprint(showerror,error)))
        end
        showerror(stderr,error);println(stderr);2
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(CalibrationResume.main())
end
