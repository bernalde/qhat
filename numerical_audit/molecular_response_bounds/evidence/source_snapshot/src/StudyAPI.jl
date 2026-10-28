module StudyAPI
using ..BoundsCommon, ..ResponseBounds, ..ExactResponse
export molecular_record, choose_comparisons, target_selection, window_record

function _phase_guard(parent,cpu_cap,wall_cap)
    start=ccall(:clock,Int64,())/1.0e6;wall=time_ns()
    ()->begin
        parent()
        ccall(:clock,Int64,())/1.0e6-start<=cpu_cap || throw(Refusal("cpu_limit","named phase CPU cap"))
        (time_ns()-wall)/1e9<=wall_cap || throw(Refusal("wall_limit","named phase wall cap"))
    end
end

function target_selection(ub,hf;count=8)
    1<=count<=8 || throw(ArgumentError("target count must be 1:8"))
    states=collect(keys(ub.response))
    sort!(states;by=d->(-(ub.response[d][1]^2+ub.response[d][2]^2),d))
    isempty(states) && return UInt64[hf==0 ? 1 : 0]
    sort!(states[1:min(count,length(states))])
end
"""Proposal rule only: largest summed coordinate-range mass, lexical ties."""
function choose_comparisons(r,count;caps=Caps(),guard=()->nothing)
    0<=count<=3 || throw(ArgumentError("comparison count must be 0:3"))
    b=Budget(caps;guard);scores=Dict{Tuple{Int,Int},BigInt}()
    for f in r.factors
        coords=Set(j for v in f.variants for j in keys(v));mass=big(0)
        for j in coords
            tick!(b,length(f.variants));values=[get(v,j,big(0)) for v in f.variants]
            mass=addchecked(mass,checked(maximum(values)-minimum(values),b),b)
        end
        iszero(mass) && continue
        for i in 1:length(f.scope)-1,j in i+1:length(f.scope)
            pair=(f.scope[i],f.scope[j]);scores[pair]=addchecked(get(scores,pair,big(0)),mass,b)
        end
    end
    candidates=sort!(collect(keys(scores));by=p->(-scores[p],p))
    candidates[1:min(count,length(candidates))]
end
function molecular_record(input,kernelpath,order,budget;target_count=8,comparison_count=1,
        comparisons=nothing,iterations=64,max_witnesses=8,caps=Caps(),guard=()->nothing,
        progress=(stage,value)->nothing)
    cpu()=ccall(:clock,Int64,())/1.0e6
    phase_start=cpu();m=load_exact_model(input,kernelpath;caps,guard);setup_cpu=cpu()-phase_start
    progress("model",Dict("input_sha256"=>m.input_sha256,"kernel_sha256"=>m.kernel_sha256,
        "original_lengths"=>original_lengths(m),"full_order"=>Int[order...],"ideal_budget"=>BigInt(budget)))
    permutation=Int[order...];sort(permutation)==collect(eachindex(m.units)) || throw(ArgumentError("not a full incumbent permutation"))
    endpoint_steps(original_lengths(m),first(permutation),last(permutation),budget)>0 || throw(ArgumentError("incumbent budget cannot afford one step"))
    phase_start=cpu();ub=exact_full_response(m,order;caps,guard=_phase_guard(guard,60,300));ub_cpu=cpu()-phase_start
    progress("full_UB",Dict("Q_inc"=>ub.Q,"UB_numerator"=>ub.numerator,"UB_denominator"=>ub.denominator,
        "response_denominator"=>ub.response_denominator,"work"=>ub.work,"estimated_peak_bytes"=>ub.estimated_peak_bytes))
    phase_start=cpu();targets=target_selection(ub,m.hf;count=target_count)
    projected=project_full_atoms(m,targets;caps,guard=_phase_guard(guard,60,300));r=projected.response
    projection_cpu=cpu()-phase_start
    progress("projection",Dict("targets"=>string.(targets),"projection_sha256"=>response_digest(r),
        "coverage_scope_sha256"=>projected.scope_sha256,"coverage_scopes"=>projected.scopes,
        "variants"=>projected.variants,"work"=>projected.work,"estimated_peak_bytes"=>projected.estimated_peak_bytes))
    phase_start=cpu();boundguard=_phase_guard(guard,10,50)
    pairs=isnothing(comparisons) ? choose_comparisons(r,comparison_count;caps,guard=boundguard) : comparisons
    cert=bound_response(r;comparisons=pairs,iterations,max_witnesses,caps,guard=boundguard)
    L=verify_bound(r,cert;caps,guard=boundguard);rat=full_ratio(m,order,ub.Q,L,budget;caps,guard=boundguard)
    bound_cpu=cpu()-phase_start
    progress("bound",Dict("certificate"=>cert,"L_Q"=>L,"L_J"=>rat.L_J,"UB_J"=>rat.UB_J,"ratio"=>rat.ratio))
    Dict{String,Any}("schema"=>"molecular-bound-record-v1","complete"=>true,
        "status"=>(rat.positive_full_bound ? "positive_full_bound" : "zero_full_bound"),
        "domain"=>"all complete orders; original unit partition",
        "input_sha256"=>m.input_sha256,"kernel_sha256"=>m.kernel_sha256,
        "targets"=>string.(targets),"target_selection"=>"largest exact same-model incumbent response energy; destination ties",
        "requested_target_count"=>target_count,
        "projection_sha256"=>response_digest(r),"coverage_scope_sha256"=>projected.scope_sha256,
        "coverage_scopes"=>projected.scopes,"full_order"=>Int[order...],"ideal_budget"=>BigInt(budget),
        "Q_inc"=>ub.Q,"UB_numerator"=>ub.numerator,"UB_denominator"=>ub.denominator,
        "response_denominator"=>ub.response_denominator,"certificate"=>cert,"L_Q"=>L,
        "L_J"=>rat.L_J,"UB_J"=>rat.UB_J,"ratio"=>rat.ratio,
        "r_inc"=>rat.r_inc,"r_max"=>rat.r_max,"positive_full_bound"=>rat.positive_full_bound,
        "desired_ratio_met"=>rat.desired_ratio_met,
        "UB_work"=>ub.work,"projection_work"=>projected.work,"projection_probes"=>projected.probes,
        "atom_variants"=>projected.variants,
        "setup_CPU_seconds"=>setup_cpu,"UB_CPU_seconds"=>ub_cpu,
        "projection_CPU_seconds"=>projection_cpu,"bound_CPU_seconds"=>bound_cpu,
        "estimated_peak_bytes"=>max(ub.estimated_peak_bytes,projected.estimated_peak_bytes),
        "interpretation"=>"Exact supplied rational coefficients and supplied dyadic kernel only. No exact-transcendental or physical-error guarantee.")
end

"""Fixed INTERIOR window record. Endpoint-moving windows refuse this v1 API."""
function window_record(flat,parent,positions,lengths,budget;input_sha256="",source_sha256="",
        comparison_count=1,iterations=64,max_witnesses=8,caps=Caps(),guard=()->nothing)
    full=Int[parent...];sort(full)==collect(1:length(lengths)) || throw(ArgumentError("window parent is not full"))
    pos=Int[positions...]
    !isempty(pos) && length(pos)==length(flat.labels) && pos==collect(first(pos):last(pos)) || throw(ArgumentError("window positions must be contiguous"))
    first(pos)>1 && last(pos)<length(full) || throw(ArgumentError("v1 window J record requires fixed exterior endpoints"))
    full[pos]==flat.labels || throw(ArgumentError("window/global label mapping differs"))
    binding=Dict{String,Any}("input_sha256"=>input_sha256,"response_source_sha256"=>source_sha256,
        "parent_order"=>full,"positions"=>pos,"original_lengths"=>lengths,"ideal_budget"=>BigInt(budget))
    r=encode_window(flat;binding,caps,guard)
    pairs=choose_comparisons(r,comparison_count;caps,guard)
    cert=bound_response(r;comparisons=pairs,iterations,max_witnesses,caps,guard)
    L=verify_bound(r,cert;caps,guard);Q=exact_value(r,flat.labels;caps,guard)
    steps=endpoint_steps(lengths,first(full),last(full),budget);steps>0 || throw(ArgumentError("window parent is unaffordable"))
    Dict{String,Any}("schema"=>"window-bound-record-v1","complete"=>true,"domain"=>"fixed_window",
        "response_sha256"=>response_digest(r),"certificate"=>cert,"L_Q"=>L,"Q_inc"=>Q,
        "r_inc"=>steps,"L_J"=>L/steps^4,"UB_J"=>Q/steps^4,"ratio"=>(L>0 ? Q/L : nothing),
        "positive_full_bound"=>false,"binding"=>binding)
end
end
