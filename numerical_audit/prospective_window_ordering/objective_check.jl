"""Numerical objective replay for every supplied scored constructor witness.

Cubic replay uses the reviewed sparse Julia Taylor action, separate from the
native Python cubic implementation. First-order replay independently forms the
exact-integer ordered pair column and applies T=1 diagonal sinc weights. Neither
comparison is an interval error bound, physical-error certificate or optimization
certificate. No primary physical driver or product WindowSearch is imported.
"""
module ObjectiveCheck
using TOML, SHA
include(joinpath(@__DIR__,"..","julia_cubic_feasibility","src","ResponseCore.jl"))
include(joinpath(@__DIR__,"src","WindowStudyIO.jl"))
using .ResponseCore: Unit
export run_check, prepare_model, partition_units, first_response, replay_objective,
       compare_scores, ObjectiveLimit

const RTOL=1e-8
const ATOL=1e-24
const MAX_WORDS=8192
const MAX_CONSERVING_UNITS=1024
const MAX_PAULI_UNITS=8192
const MAX_ENTRIES=50_000
const MAX_ROWS=1024
const MAX_UNIQUE=256
const MAX_PAIR_EVENTS=35_000_000
const MAX_INTEGER_BITS=16_384
const MAX_WORK=200_000_000
const KINDS=("diagonal_reference_cubic_response","retained_first_order_response")
cpu_seconds()=ccall(:clock,Int64,())/1.0e6
filehash(path)=bytes2hex(sha256(read(path)))
struct ObjectiveLimit <: Exception
    status::String
    reason::String
end
Base.showerror(io::IO,e::ObjectiveLimit)=print(io,e.status,": ",e.reason)
mutable struct ReplayBudget
    work::Int
    limit::Int
    check::Any
end
ReplayBudget(check=()->nothing;max_work=MAX_WORK)=begin
    max_work isa Integer && !(max_work isa Bool) && 0<max_work<=MAX_WORK || throw(ArgumentError("objective work limit outside ceiling"))
    check();ReplayBudget(0,Int(max_work),check)
end
function charge!(b::ReplayBudget,n::Int=1)
    0<=n<=b.limit-b.work || throw(ObjectiveLimit("work_refusal","objective action work limit"))
    old=b.work;b.work+=n;old÷1024!=b.work÷1024 && b.check();nothing
end
bits(x::BigInt)=iszero(x) ? 0 : ndigits(abs(x);base=2)
function checked(x::BigInt,b)
    charge!(b);bits(x)<=MAX_INTEGER_BITS || throw(ObjectiveLimit("integer_refusal","objective integer limit"));x
end
function product(a::BigInt,c::BigInt,b)
    bits(a)+bits(c)<=MAX_INTEGER_BITS || throw(ObjectiveLimit("integer_refusal","prospective objective integer product"))
    checked(a*c,b)
end
struct ObjectiveWord
    key::Tuple
    flip::UInt64
    parity::UInt64
    signed_numerator::BigInt
    original_unit::Int
    original_word::Int
end
struct ObjectiveBlock
    flip::UInt64
    words::Vector{ObjectiveWord}
end
struct ObjectiveModel
    data::Any
    denominator::BigInt
    blocks::Vector{ObjectiveBlock}
    sorted_words::Vector{ObjectiveWord}
end
function compile_word(text,n,numerator,unit_index,word_index,b)
    text isa AbstractString && ncodeunits(text)<=4n || throw(ArgumentError("bounded Pauli spelling required"))
    key=Tuple{Int,Char}[];seen=Set{Int}();x=z=UInt64(0);ys=0
    for token in (isempty(text) ? String[] : split(text,','))
        charge!(b)
        m=match(r"^([0-9]+)([XYZ])$",token);isnothing(m) && throw(ArgumentError("invalid objective Pauli token"))
        q=parse(Int,m.captures[1]);p=only(m.captures[2])
        0<=q<n && !(q in seen) || throw(ArgumentError("repeated/outside objective site"))
        push!(seen,q);push!(key,(q,p));p!='Z' && (x|=UInt64(1)<<q);p!='X' && (z|=UInt64(1)<<q);ys+=p=='Y'
    end
    issorted(key) && iseven(ys) || throw(ArgumentError("objective requires canonical real Pauli words"))
    numerator isa AbstractString && ncodeunits(numerator)<=1280 || throw(ArgumentError("bounded objective numerator spelling"))
    a=parse(BigInt,numerator);!iszero(a) && bits(a)<=4096 || throw(ArgumentError("bounded nonzero objective numerator"))
    phase=ys%4==0 ? 1 : -1
    ObjectiveWord(Tuple(key),x,z,phase==1 ? a : -a,unit_index,word_index)
end
function prepare_model(data;budget=ReplayBudget())
    1<=data.n<=14 || throw(ObjectiveLimit("register_refusal","objective pilot supports at most14 modes"))
    data.partition=="conserving" || throw(ArgumentError("objective input must retain original conserving partition"))
    1<=length(data.units)<=MAX_CONSERVING_UNITS && 0<data.words<=MAX_WORDS || throw(ObjectiveLimit("input_refusal","original objective model counts"))
    denominator=parse(BigInt,data.source["denominator"]);denominator>0 && bits(denominator)<=4096 || throw(ArgumentError("objective denominator bound"))
    blocks=ObjectiveBlock[];words=ObjectiveWord[];keys=Set{Tuple}()
    for (i,raw) in enumerate(data.source["units"])
        compiled=ObjectiveWord[]
        for (j,(text,a)) in enumerate(zip(raw["words"],raw["numerators"]))
            word=compile_word(text,data.n,a,i,j,budget)
            word.key in keys && throw(ArgumentError("duplicate global Pauli word makes singleton labels ambiguous"))
            push!(keys,word.key);push!(compiled,word);push!(words,word)
        end
        !isempty(compiled) && all(w->w.flip==compiled[1].flip,compiled) || throw(ArgumentError("objective equal-flip condition"))
        push!(blocks,ObjectiveBlock(compiled[1].flip,compiled))
    end
    length(words)==data.words || throw(ArgumentError("objective word inventory differs"))
    sort!(words;by=w->w.key) # Numeric structured labels, never raw string order.
    ObjectiveModel(data,denominator,blocks,words)
end
function partition_units(model::ObjectiveModel,partition::String)
    if partition=="conserving"
        return (blocks=model.blocks,floating=model.data.units)
    elseif partition=="pauli"
        length(model.sorted_words)<=MAX_PAULI_UNITS || throw(ObjectiveLimit("unit_refusal","Pauli objective unit count"))
        descriptors=ObjectiveBlock[];floating=Unit[]
        for w in model.sorted_words
            push!(descriptors,ObjectiveBlock(w.flip,[w]))
            original=model.data.units[w.original_unit];j=w.original_word
            push!(floating,Unit(original.flip,[original.parity[j]],[original.signs[j]],[original.coeff[j]]))
        end
        return (blocks=descriptors,floating=floating)
    end
    throw(ArgumentError("unknown objective partition"))
end
function validate_order(raw,m)
    raw isa AbstractVector && length(raw)==m && all(x->x isa Integer && !(x isa Bool),raw) || throw(ArgumentError("complete integer objective order required"))
    order=Int.(raw);sort(order)==collect(1:m) || throw(ArgumentError("objective order is not a permutation"));order
end
function value(block::ObjectiveBlock,state::UInt64,b)
    total=big(0)
    for w in block.words
        charge!(b);term=iseven(count_ones(state&w.parity)) ? w.signed_numerator : -w.signed_numerator
        total=checked(total+term,b)
    end
    total
end
function first_response(model::ObjectiveModel,partition::String,raworder;budget=ReplayBudget())
    units=partition_units(model,partition).blocks;m=length(units);order=validate_order(raworder,m)
    pair_events=big(m)*(m-1)÷2
    pair_events<=MAX_PAIR_EVENTS || throw(ObjectiveLimit("pair_refusal","ordered first-response pair limit"))
    hf=model.data.hf;at_hf=[value(u,hf,budget) for u in units]
    column=Dict{UInt64,BigInt}();nonzero_pairs=0
    for left in 1:m-1,right in left+1:m
        charge!(budget);a=order[left];c=order[right]
        destination=hf⊻units[a].flip⊻units[c].flip
        destination==hf && continue # Pperp; not a sector projection.
        iszero(at_hf[a]) && iszero(at_hf[c]) && continue
        if length(units[a].words)==length(units[c].words)==1
            wa=only(units[a].words);wc=only(units[c].words)
            iseven(count_ones(wa.flip&wc.parity)+count_ones(wc.flip&wa.parity)) && continue
        end
        # Literal Ai*Aj - Aj*Ai action, in the SUPPLIED induced order.
        ac=product(at_hf[c],value(units[a],hf⊻units[c].flip,budget),budget)
        ca=product(at_hf[a],value(units[c],hf⊻units[a].flip,budget),budget)
        delta=checked(ac-ca,budget);iszero(delta) && continue;nonzero_pairs+=1
        updated=checked(get(column,destination,big(0))+delta,budget)
        if iszero(updated)
            delete!(column,destination)
        else
            haskey(column,destination) || length(column)<MAX_ENTRIES || throw(ObjectiveLimit("support_refusal","first response support cap"))
            column[destination]=updated
        end
    end
    # column / (2DA²) is 1/2 sum_[i<j in order] [Bi,Bj]|HF>.
    scale=product(big(2),product(model.denominator,model.denominator,budget),budget)
    diagonal=[u for u in model.blocks if iszero(u.flip)]
    initial_diagonal=sum((value(u,hf,budget) for u in diagonal);init=big(0))
    norm_response=0.0
    for (destination,integer) in column
        amplitude=Float64(integer//scale)
        isfinite(amplitude) && !iszero(amplitude) || throw(ObjectiveLimit("numerical_refusal","first-column conversion is nonfinite or underflows"))
        shifted=sum((value(u,destination,budget) for u in diagonal);init=big(0))
        gap=Float64((shifted-initial_diagonal)//model.denominator)
        x=gap/2 # Fixed training T=1; no physical validation time is consumed.
        isfinite(x) || throw(ObjectiveLimit("numerical_refusal","first transport phase is nonfinite"))
        sinc_ratio=iszero(x) ? 1.0 : sin(x)/x
        transported=amplitude*sinc_ratio
        isfinite(transported) || throw(ObjectiveLimit("numerical_refusal","first transported amplitude is nonfinite"))
        norm_response=hypot(norm_response,abs(transported))
    end
    score=norm_response^2
    isfinite(score) || throw(ObjectiveLimit("numerical_refusal","first squared norm is nonfinite"))
    !isempty(column) && iszero(score) && throw(ObjectiveLimit("numerical_refusal","nonzero first response underflows after transport/square"))
    (value=score,column=column,column_denominator=scale,steps=nothing,r_normalized=false,
     nonzero_pairs=nonzero_pairs,pair_events=Int(pair_events),destinations=length(column),
     scope="Exact-integer projected ordered pair action, independently rounded T=1 diagonal sinc norm; unit-modulus phase cancels in norm; NO r normalization")
end
function replay_objective(model,partition,order,kind;budget=ReplayBudget())
    kind in KINDS || throw(ArgumentError("unsupported scored objective kind"))
    budget.check()
    if kind=="retained_first_order_response"
        result=first_response(model,partition,order;budget)
        return (value=result.value,steps="not_used",details=Dict("pair_events"=>result.pair_events,
            "nonzero_pairs"=>result.nonzero_pairs,"destinations"=>result.destinations,"r_normalized"=>false,"scope"=>result.scope))
    end
    units=partition_units(model,partition).floating;permutation=validate_order(order,length(units))
    score=ResponseCore.score_order(units,permutation,model.data.hf,30model.data.words;
        time=1.0,max_entries=MAX_ENTRIES,guard=true)
    budget.check()
    (value=score.value,steps=score.steps,details=Dict("cubic_response"=>score.cubic_response,
        "r_normalized"=>true,"scope"=>"Reviewed Julia sparse degree3 Taylor/diagonal transport, separate from native Python cubic implementation; numerical parity only"))
end
function compare_scores(actual::Real,claimed::Real)
    a=Float64(actual);c=Float64(claimed)
    all(x->isfinite(x) && x>=0,(a,c)) || throw(ArgumentError("finite nonnegative objective scores required"))
    difference=abs(a-c);tolerance=ATOL+RTOL*max(a,c)
    within=difference<=tolerance;floor=max(a,c)<=ATOL
    status=!within ? "mismatch" : floor ? "floor_censored" : "verified"
    (status=status,verified=within && !floor,within_tolerance=within,numerical_floor=floor,
     absolute_difference=difference,relative_difference=max(a,c)>0 ? difference/max(a,c) : 0.0,
     tolerance=tolerance)
end
function source_hashes()
    names=("objective_check.jl","src/WindowStudyIO.jl","../julia_cubic_feasibility/src/ResponseCore.jl",
        "../scalable_ordering/training.py","../../analysis/sparse_response_ordering.py")
    Dict(name=>filehash(joinpath(@__DIR__,name)) for name in names)
end
function write_new(path,record)
    ispath(path) && throw(ObjectiveLimit("io_refusal","refuse existing objective output"))
    partial=path*".partial";ispath(partial) && throw(ObjectiveLimit("io_refusal","refuse retained partial objective output"))
    io=IOBuffer();TOML.print(io,WindowStudyIO.clean(record);sorted=true);bytes=take!(io)
    length(bytes)<=32*1024^2 || throw(ObjectiveLimit("output_refusal","objective output byte cap"))
    mkpath(dirname(abspath(path)));open(partial,"w") do f;write(f,bytes);end
    mv(partial,path;force=false);nothing
end
function candidate_identity(row,index)
    name=get(row,"name","row_"*string(index))
    name isa String && 0<ncodeunits(name)<=256 || throw(ArgumentError("bounded candidate name required"))
    out=Dict{String,Any}("candidate_index"=>index,"name"=>name)
    for field in ("arm","repetition","formula","time","budget_multiplier")
        haskey(row,field) || continue
        value=row[field]
        valid=field=="arm" ? (value isa String && ncodeunits(value)<=64) :
            (value isa Real && !(value isa Bool) && isfinite(value) && value>=0)
        valid || throw(ArgumentError("candidate identity metadata must be bounded scalar values"))
        out[field]=value
    end
    out
end
function scored_row(row,model)
    partition=get(row,"partition",nothing);partition isa String && partition in ("conserving","pauli") || throw(ArgumentError("candidate partition required"))
    m=partition=="conserving" ? length(model.blocks) : length(model.sorted_words)
    order=validate_order(get(row,"order",nothing),m)
    raw=get(row,"constructor_score",nothing)
    if isnothing(raw) || raw=="unavailable"
        return (scored=false,partition=partition,order=order,kind="unscored",score=nothing,key="")
    end
    raw isa Real && !(raw isa Bool) && isfinite(raw) && raw>=0 || throw(ArgumentError("constructor_score must be finite nonnegative numeric or explicitly unavailable"))
    score=Float64(raw);isfinite(score) || throw(ArgumentError("constructor score Float64 range"))
    kind=get(row,"score_kind",nothing);kind isa String && kind in KINDS || throw(ArgumentError("scored row has unknown/missing score_kind"))
    # The lookup compares the whole canonical string, not only its SHA digest.
    key=partition*"|"*kind*"|"*bitstring(score)*"|"*join(order,',')
    (scored=true,partition=partition,order=order,kind=kind,score=score,key=key)
end
struct RunBudgetRefusal <: Exception
    cause::Exception
end
Base.showerror(io::IO,e::RunBudgetRefusal)=showerror(io,e.cause)
const LOCAL_CHECK_REFUSALS=("numerical_refusal","support_refusal","integer_refusal",
                          "pair_refusal","unit_refusal","register_refusal")
function run_refusal_status(e,phase)
    if e isa RunBudgetRefusal
        return e.cause isa ObjectiveLimit ? e.cause.status : "budget_refusal"
    end
    e isa ObjectiveLimit && return e.status
    e isa Union{Base.IOError,Base.SystemError} && return "io_refusal"
    phase=="replay" ? "replay_error" : phase=="checkpoint_io" ? "checkpoint_error" : "input_or_manifest_error"
end
function run_check(input,candidates_path,output;cpu_limit=180.0,wall_limit=900.0,
        check_budget=()->nothing,max_work=MAX_WORK,max_unique=MAX_UNIQUE)
    0<cpu_limit<=180 && 0<wall_limit<=900 || throw(ArgumentError("objective process allowance outside ceiling"))
    max_unique isa Integer && !(max_unique isa Bool) && 0<max_unique<=MAX_UNIQUE ||
        throw(ArgumentError("objective unique allowance outside ceiling"))
    ispath(output) && error("refuse existing objective summary")
    began=cpu_seconds();wall=time_ns()
    guard=()->begin
        try
            check_budget()
        catch e
            throw(RunBudgetRefusal(e))
        end
        cpu_seconds()-began<=cpu_limit || throw(ObjectiveLimit("cpu_refusal","whole objective replay CPU allowance"))
        (time_ns()-wall)/1e9<=wall_limit || throw(ObjectiveLimit("wall_refusal","whole objective replay wall allowance"))
    end
    header=Dict{String,Any}("schema"=>"prospective-window-objective-check-v1","status"=>"incomplete","complete"=>false,
        "all_scored_objectives_verified"=>false,"source_hashes"=>source_hashes(),
        "relative_tolerance"=>RTOL,"absolute_tolerance"=>ATOL,"training_time"=>1.0,
        "training_budget_multiplier"=>30,"numerical_comparison_only"=>true,
        "scope"=>"Every supplied scored witness, deduplicated by exact partition/order/kind/score; unscored and floor/refused checks retained. No physical/interval/optimality certificate.",
        "logical_rows"=>Dict{String,Any}[],"checks"=>Dict{String,Any}[],"new_physical_propagation"=>false,
        "input_file"=>abspath(input),"candidates_file"=>abspath(candidates_path))
    budget=ReplayBudget(()->nothing;max_work);budget.check=guard
    bykey=Dict{String,Int}();model=nothing;aborted=false;phase="checkpoint_io"
    header["unique_limit"]=Int(max_unique)
    try
        write_new(output*".START.toml",header)
        phase="input_or_manifest"
        guard()
        filesize(input)<=8*1024^2 && filesize(candidates_path)<=8*1024^2 || throw(ObjectiveLimit("input_refusal","objective input/candidate byte cap"))
        header["input_sha256"]=filehash(input);header["candidates_sha256"]=filehash(candidates_path)
        data=WindowStudyIO.load_model(input);model=prepare_model(data;budget)
        manifest=TOML.parsefile(candidates_path);rows=manifest["candidates"]
        manifest["input_sha256"]==data.input_sha256 || throw(ArgumentError("candidate/input hash mismatch"))
        haskey(manifest,"source_json_sha256") && manifest["source_json_sha256"]!=get(data.source,"source_json_sha256","") && throw(ArgumentError("candidate original JSON hash mismatch"))
        rows isa Vector && length(rows)<=MAX_ROWS || throw(ArgumentError("bounded candidate rows required"))
        header["input_sha256"]=data.input_sha256;header["candidates_sha256"]=filehash(candidates_path)
        header["case"]=data.case;header["n"]=data.n;header["hf"]=string(data.hf)
        header["original_words"]=data.words;header["training_rotation_budget"]=30data.words
        header["expected_logical_rows"]=length(rows)
        for (index,row) in enumerate(rows)
            phase="row_guard";guard()
            info=Dict{String,Any}("candidate_index"=>index,"status"=>"pending")
            push!(header["logical_rows"],info)
            parsed=nothing;phase="candidate_validation"
            try
                row isa AbstractDict || throw(ArgumentError("candidate row must be a table"))
                merge!(info,candidate_identity(row,index))
                parsed=scored_row(row,model)
                info["partition"]=parsed.partition
            catch e
                e isa Union{ArgumentError,KeyError,TypeError} || rethrow()
                info["status"]="invalid_candidate"
                info["reason"]=sprint(showerror,e)
                continue
            end
            if !parsed.scored
                info["status"]="unscored_censored";info["reason"]="No numeric constructor objective was supplied"
                continue
            end
            if haskey(bykey,parsed.key)
                info["check_index"]=bykey[parsed.key];info["deduplicated"]=true
                info["status"]=header["checks"][info["check_index"]]["status"]
                continue
            end
            phase="unique_admission"
            length(header["checks"])<max_unique || throw(ObjectiveLimit("unique_refusal","objective task count cap"))
            check_index=length(header["checks"])+1;bykey[parsed.key]=check_index
            info["check_index"]=check_index;info["deduplicated"]=false;info["status"]="started"
            rec=Dict{String,Any}("check_index"=>check_index,"status"=>"started","verified"=>false,
                "partition"=>parsed.partition,"order"=>parsed.order,"score_kind"=>parsed.kind,
                "constructor_score"=>parsed.score,"objective_sha256"=>bytes2hex(sha256(parsed.key)))
            push!(header["checks"],rec)
            checkpoint=output*".CHECK."*lpad(string(check_index),3,'0')
            phase="checkpoint_io";write_new(checkpoint*".START.toml",rec)
            phase="replay"
            try
                guard();replayed=replay_objective(model,parsed.partition,parsed.order,parsed.kind;budget)
                comparison=compare_scores(replayed.value,parsed.score)
                rec["replayed_score"]=replayed.value;rec["steps"]=replayed.steps;rec["details"]=replayed.details
                for (key,value) in pairs(comparison);rec[string(key)]=value;end
            catch e
                rec["status"]=e isa ResponseCore.NumericalResolution ? "numerical_refusal" :
                    e isa ObjectiveLimit ? e.status : run_refusal_status(e,phase)
                rec["reason"]=sprint(showerror,e);info["status"]=rec["status"]
                (e isa ResponseCore.NumericalResolution || e isa ObjectiveLimit && e.status in LOCAL_CHECK_REFUSALS) || rethrow()
            end
            info["status"]=rec["status"]
            phase="checkpoint_io"
            write_new(checkpoint*".RESULT.toml",rec)
            write_new(output*".PARTIAL."*lpad(string(check_index),3,'0')*".toml",header)
        end
    catch e
        aborted=true;header["status"]=run_refusal_status(e,phase)
        header["reason"]=sprint(showerror,e)
        header["failed_phase"]=phase
        for rec in header["checks"]
            if rec["status"]=="started";rec["status"]=header["status"];rec["reason"]=header["reason"];end
        end
        for row in header["logical_rows"]
            if row["status"] in ("pending","started");row["status"]=header["status"];row["reason"]=header["reason"];end
        end
    end
    header["run_aborted"]=aborted;header["processed_logical_rows"]=length(header["logical_rows"])
    header["unscored_rows"]=count(r->r["status"]=="unscored_censored",header["logical_rows"])
    header["scored_logical_rows"]=count(r->haskey(r,"check_index"),header["logical_rows"])
    header["unique_objectives"]=length(header["checks"])
    header["all_scored_objectives_verified"]=!aborted && !isempty(header["checks"]) && all(r->r["verified"],header["checks"])
    valid_rows=all(r->(haskey(r,"check_index") && !haskey(r,"reason")) || r["status"]=="unscored_censored",header["logical_rows"])
    header["complete"]=!aborted && valid_rows && header["all_scored_objectives_verified"]
    header["mismatched_objectives"]=count(r->r["status"]=="mismatch",header["checks"])
    header["floor_censored_objectives"]=count(r->r["status"]=="floor_censored",header["checks"])
    if !aborted
        header["status"]=header["complete"] ? "complete_for_supplied_scored_orders" : !valid_rows ? "invalid_candidates" :
            header["mismatched_objectives"]>0 ? "score_mismatch" : isempty(header["checks"]) ? "no_scored_orders_supplied" : "partial_or_censored"
    end
    header["charged_work_first_and_compile"]=budget.work
    header["work_scope"]="Own first-order/compiler counters; cubic action and order/IO parsing use bounded dimensions plus the owning process CPU/memory guard"
    header["CPU_seconds"]=cpu_seconds()-began;header["wall_seconds"]=(time_ns()-wall)/1e9
    try
        write_new(output,header)
    catch e
        header["summary_write_refusal"]=run_refusal_status(e,"checkpoint_io")
        header["summary_write_reason"]=sprint(showerror,e)
        !aborted && (header["status"]=header["summary_write_refusal"])
        header["complete"]=false;header["all_scored_objectives_verified"]=false;header["run_aborted"]=true
    end
    header
end
function main(args=ARGS)
    length(args)==3 || error("usage: objective_check.jl INPUT_TOML CANDIDATES_TOML NEW_CHECK_TOML")
    result=run_check(args[1],args[2],args[3]);result["complete"] ? 0 : 2
end
end # module
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(ObjectiveCheck.main())
end
