module BoundsCommon
using SHA, TOML
export Caps, Budget, Refusal, tick!, reserve!, checked, mulchecked, addchecked,
       bigparse, rationalparse, filehash, clean, write_record, permutations,
       VectorFactor, IntegerResponse, response_digest, validate_response,
       exact_value, factor_choice, bitlength, scope_digest, caps_dict,
       select_caps, verify_caps_declaration, parse_caps_cli

Base.@kwdef struct Caps
    max_n::Int = 14
    max_units::Int = 512
    max_words::Int = 4096
    max_entries::Int = 1_000_000
    max_states::Int = 10_000
    max_cache::Int = 100_000
    max_scopes::Int = 50_000
    max_variants::Int = 300_000
    max_probes::Int = 2_000_000
    max_work::Int = 50_000_000
    max_bits::Int = 4096
    max_bytes::Int = 128*1024^2
    max_input_bytes::Int = 8*1024^2
    max_output_bytes::Int = 32*1024^2
end
caps_dict(c::Caps)=Dict(string(f)=>getfield(c,f) for f in fieldnames(Caps))
const CAPS_PROFILES=("default","larger14")

"""Select a bounded CLI profile; only unit/word count overrides are permitted.

The larger14 profile changes512/4096 to1024/8192 and fixes4.5m probes rather
than2m. Register, work, memory, support, arithmetic, input/output and other
ceilings remain frozen. Probe count is NOT a free registry override.
An explicit CLI profile overrides the registry profile; explicit registry counts
still apply and must fit the selected profile. Both requests are retained.
"""
function select_caps(config=Dict{String,Any}();cli_profile=nothing)
    config isa AbstractDict || throw(ArgumentError("caps must be a table"))
    Set(keys(config))⊆Set(["profile","max_units","max_words"]) || throw(ArgumentError("unknown or immutable caps field"))
    raw=Dict{String,Any}(String(k)=>v for (k,v) in config)
    registry_profile=get(raw,"profile","default")
    registry_profile isa AbstractString && registry_profile in CAPS_PROFILES || throw(ArgumentError("unknown registry caps profile"))
    isnothing(cli_profile) || (cli_profile isa AbstractString && cli_profile in CAPS_PROFILES) || throw(ArgumentError("unknown CLI caps profile"))
    profile=isnothing(cli_profile) ? String(registry_profile) : String(cli_profile)
    unit_ceiling,word_ceiling=profile=="larger14" ? (1024,8192) : (512,4096)
    units=get(raw,"max_units",unit_ceiling);words=get(raw,"max_words",word_ceiling)
    units isa Integer && !(units isa Bool) && 1<=units<=unit_ceiling || throw(ArgumentError("max_units outside selected profile"))
    words isa Integer && !(words isa Bool) && 1<=words<=word_ceiling || throw(ArgumentError("max_words outside selected profile"))
    haskey(raw,"profile") && (raw["profile"]=String(registry_profile))
    haskey(raw,"max_units") && (raw["max_units"]=Int(units))
    haskey(raw,"max_words") && (raw["max_words"]=Int(words))
    caps=Caps(max_units=Int(units),max_words=Int(words),
        max_probes=(profile=="larger14" ? 4_500_000 : 2_000_000))
    declaration=Dict{String,Any}("schema"=>"bounded-caps-v1","registry"=>raw,
        "cli_profile"=>(isnothing(cli_profile) ? "none" : String(cli_profile)),
        "selection_source"=>(!isnothing(cli_profile) ? "cli" : isempty(raw) ? "default" : "registry"),
        "effective_profile"=>profile,"effective"=>caps_dict(caps),
        "hard_ceilings"=>caps_dict(Caps(max_units=1024,max_words=8192,max_probes=4_500_000)))
    caps,declaration
end

"""Rebuild from whitelisted requests and compare EVERY effective/ceiling field.
No record-supplied work/byte/bit/register value is passed into a Caps constructor.
"""
function verify_caps_declaration(declaration;expected_profile=nothing)
    declaration isa AbstractDict || throw(ArgumentError("caps declaration must be a table"))
    keys_expected=Set(["schema","registry","cli_profile","selection_source","effective_profile","effective","hard_ceilings"])
    Set(keys(declaration))==keys_expected && declaration["schema"]=="bounded-caps-v1" || throw(ArgumentError("caps declaration schema differs"))
    stored_cli=declaration["cli_profile"]
    stored_cli isa AbstractString || throw(ArgumentError("caps CLI profile type differs"))
    caps,expected=select_caps(declaration["registry"];cli_profile=(stored_cli=="none" ? nothing : stored_cli))
    for key in ("effective","hard_ceilings")
        values=declaration[key]
        values isa AbstractDict && Set(keys(values))==Set(string.(fieldnames(Caps))) || throw(ArgumentError("caps inventory differs"))
        all(v->v isa Integer && !(v isa Bool),Base.values(values)) || throw(ArgumentError("caps values must be integers"))
    end
    clean(declaration)==clean(expected) || throw(ArgumentError("declared caps exceed or differ from bounded profile"))
    isnothing(expected_profile) || (expected_profile in CAPS_PROFILES && expected_profile==expected["effective_profile"]) || throw(ArgumentError("replay caps profile does not match requested profile"))
    caps,expected
end

"""Parse only --caps-profile PROFILE plus positional arguments; no evaluation."""
function parse_caps_cli(args)
    profile=nothing;positional=String[];i=1
    while i<=length(args)
        arg=String(args[i])
        if arg=="--caps-profile"
            isnothing(profile) && i<length(args) || throw(ArgumentError("missing or repeated caps profile"))
            profile=String(args[i+1]);profile in CAPS_PROFILES || throw(ArgumentError("unknown caps profile"));i+=2
        elseif startswith(arg,"--")
            throw(ArgumentError("unknown CLI option: "*arg))
        else
            push!(positional,arg);i+=1
        end
    end
    profile,positional
end
struct Refusal <: Exception
    status::String
    message::String
end
Base.showerror(io::IO, e::Refusal) = print(io,e.status,": ",e.message)
mutable struct Budget
    caps::Caps
    guard::Any
    work::Int
    bytes::Int
    entries::Int
    peak_bytes::Int
end
function Budget(caps::Caps=Caps();guard=()->nothing)
    all(x->x>0,(getfield(caps,f) for f in fieldnames(Caps))) || throw(ArgumentError("caps must be positive"))
    guard(); Budget(caps,guard,0,0,0,0)
end
function tick!(b::Budget,n::Integer=1)
    0<=n<=b.caps.max_work-b.work || throw(Refusal("work_limit","charged work exhausted"))
    old=b.work; b.work+=Int(n)
    old÷1024!=b.work÷1024 && b.guard()
    nothing
end
function reserve!(b::Budget,entries::Integer;bytes::Integer=0)
    0<=entries<=b.caps.max_entries-b.entries || throw(Refusal("entry_limit","retained scalar entries exceed cap"))
    # Conservative BigInt/dictionary allowance including scratch/copies. Storage
    # charges are monotone within a named component, even after local deletion.
    required=BigInt(bytes)+BigInt(entries)*(256+4cld(b.caps.max_bits,8))
    0<=required<=b.caps.max_bytes-b.bytes || throw(Refusal("memory_limit","estimated allocation exceeds cap"))
    b.entries+=Int(entries); b.bytes+=Int(required); b.peak_bytes=max(b.peak_bytes,b.bytes)
    nothing
end
bitlength(x::Integer)=iszero(x) ? 0 : ndigits(abs(x);base=2)
function checked(x::BigInt,b::Budget)
    bitlength(x)<=b.caps.max_bits || throw(Refusal("bit_limit","integer exceeds cap")); x
end
function mulchecked(a::BigInt,c::BigInt,b::Budget)
    tick!(b)
    bitlength(a)+bitlength(c)<=b.caps.max_bits || throw(Refusal("bit_limit","prospective product exceeds cap"))
    checked(a*c,b)
end
addchecked(a::BigInt,c::BigInt,b::Budget)=(tick!(b);checked(a+c,b))
function bigparse(x,caps::Caps=Caps())
    s=string(x)
    length(s)<=ceil(Int,caps.max_bits*log10(2))+2 || throw(Refusal("bit_limit","integer text exceeds cap"))
    occursin(r"^-?[0-9]+$",s) || throw(ArgumentError("invalid integer text"))
    v=parse(BigInt,s); bitlength(v)<=caps.max_bits || throw(Refusal("bit_limit","integer exceeds cap"));v
end
function rationalparse(x,caps::Caps=Caps())
    x isa Rational && return BigInt(numerator(x))//BigInt(denominator(x))
    x isa Integer && return BigInt(x)//big(1)
    parts=split(String(x),"//")
    length(parts)==1 && return bigparse(parts[1],caps)//big(1)
    length(parts)==2 || throw(ArgumentError("invalid rational text"))
    den=bigparse(parts[2],caps);den>0 || throw(ArgumentError("nonpositive denominator"))
    bigparse(parts[1],caps)//den
end
filehash(path)=bytes2hex(sha256(read(path)))
function clean(x)
    x===nothing && return "unavailable"
    x isa Union{BigInt,Rational,UInt64} && return string(x)
    x isa Symbol && return String(x)
    x isa NamedTuple && return Dict(string(k)=>clean(v) for (k,v) in pairs(x))
    x isa AbstractDict && return Dict(string(k)=>clean(v) for (k,v) in x)
    x isa Union{AbstractVector,Tuple} && return [clean(v) for v in x]
    x isa AbstractFloat && !isfinite(x) && return string(x)
    x
end
function write_record(path,value;caps=Caps())
    ispath(path) && error("refuse to overwrite output")
    # Reject oversized trees before materializing integer strings / IOBuffer.
    function forecast(x)
        if x isa AbstractDict
            return sum(64+ncodeunits(string(k))+forecast(v) for (k,v) in x;init=0)
        elseif x isa NamedTuple
            return forecast(Dict(pairs(x)))
        elseif x isa Union{AbstractVector,Tuple}
            return 32+sum(16+forecast(v) for v in x;init=0)
        elseif x isa Integer
            return 8+cld(4max(1,bitlength(x)),10)
        elseif x isa Rational
            return 16+cld(4max(1,bitlength(numerator(x))+bitlength(denominator(x))),10)
        else
            return 16+ncodeunits(string(x))
        end
    end
    forecast(value)<=caps.max_output_bytes || throw(Refusal("output_limit","serialized tree forecast exceeds cap"))
    io=IOBuffer();TOML.print(io,clean(value);sorted=true); bytes=take!(io)
    length(bytes)<=caps.max_output_bytes || throw(Refusal("output_limit","serialized output exceeds cap"))
    mkpath(dirname(abspath(path)));open(path,"w") do f;write(f,bytes);end
    nothing
end
function permutations(labels)
    length(labels)<=3 || throw(ArgumentError("only pair/triple variants are generated"))
    isempty(labels) && return [Int[]]
    [[x;tail] for x in labels for tail in permutations([y for y in labels if y!=x])]
end

"""Sparse stacked real/imaginary vector variants, with EXPLICIT induced orders."""
struct VectorFactor
    scope::Vector{Int}
    orders::Vector{Vector{Int}}
    variants::Vector{Dict{Int,BigInt}}
end
struct IntegerResponse
    labels::Vector{Int}
    dimension::Int
    denominator::BigInt
    baseline::Dict{Int,BigInt}
    factors::Vector{VectorFactor}
    domain::String
    binding::Dict{String,Any}
end
function response_digest(r::IntegerResponse)
    io=IOBuffer();println(io,"integer-response-v1\n",r.domain)
    println(io,join(r.labels,','),"\n",r.dimension,"\n",r.denominator)
    for (j,v) in sort!(collect(r.baseline);by=first);println(io,j,':',v);end
    println(io,"factors:")
    for f in r.factors
        println(io,join(f.scope,','))
        for (order,vector) in zip(f.orders,f.variants)
            println(io,join(order,','))
            for (j,v) in sort!(collect(vector);by=first);println(io,j,':',v);end
            println(io,";")
        end
    end
    # All bindings are flat strings/integers/vectors; sorted TOML is canonical.
    TOML.print(io,clean(r.binding);sorted=true)
    bytes2hex(sha256(take!(io)))
end
function scope_digest(scopes)
    io=IOBuffer()
    for s in sort!(collect(scopes));println(io,join(s,','));end
    bytes2hex(sha256(take!(io)))
end
function validate_response(r::IntegerResponse,b::Budget)
    1<=length(r.labels)<=b.caps.max_units || throw(ArgumentError("invalid label count"))
    length(unique(r.labels))==length(r.labels) || throw(ArgumentError("duplicate labels"))
    0<=r.dimension<=2b.caps.max_states || throw(Refusal("dimension_limit","response dimension exceeds cap"))
    r.denominator>0 || throw(ArgumentError("nonpositive response denominator"))
    checked(r.denominator,b); mulchecked(r.denominator,r.denominator,b)
    r.domain in ("fixed_window","full_order_projection") || throw(ArgumentError("unrecognized response domain"))
    length(r.factors)<=b.caps.max_scopes || throw(Refusal("scope_limit","too many vector factors"))
    reserve!(b,0;bytes=256length(r.factors)+256sum(length(f.variants) for f in r.factors;init=0))
    seen=Set{Tuple{Vararg{Int}}}()
    for f in r.factors
        tick!(b,1+length(f.scope))
        2<=length(f.scope)<=3 && issorted(f.scope) && length(unique(f.scope))==length(f.scope) ||
            throw(ArgumentError("factor scope must be sorted distinct pair/triple"))
        all(x->x in r.labels,f.scope) || throw(ArgumentError("factor label outside domain"))
        Tuple(f.scope) in seen && throw(ArgumentError("duplicate scope"));push!(seen,Tuple(f.scope))
        length(f.orders)==length(f.variants)==factorial(length(f.scope)) || throw(ArgumentError("wrong variant count"))
        Set(Tuple(p) for p in f.orders)==Set(Tuple(p) for p in permutations(f.scope)) || throw(ArgumentError("invalid induced-order variants"))
    end
    vectors=Iterators.flatten(((r.baseline,),Iterators.flatten(f.variants for f in r.factors)))
    for vector in vectors
        reserve!(b,length(vector));tick!(b,length(vector))
        for (j,v) in vector
            1<=j<=r.dimension && !iszero(v) || throw(ArgumentError("invalid sparse coordinate"));checked(v,b)
        end
    end
    nothing
end
factor_choice(f::VectorFactor,positions)=findfirst(p->all(positions[p[j]]<positions[p[j+1]] for j in 1:length(p)-1),f.orders)
function exact_value(r::IntegerResponse,order;caps=Caps(),guard=()->nothing)
    b=Budget(caps;guard);validate_response(r,b)
    sort(collect(order))==sort(r.labels) || throw(ArgumentError("not a domain permutation"))
    pos=Dict(x=>i for (i,x) in enumerate(order));v=copy(r.baseline)
    for f in r.factors
        k=factor_choice(f,pos);isnothing(k) && error("missing order variant")
        for (j,a) in f.variants[k];v[j]=addchecked(get(v,j,big(0)),a,b);end
    end
    total=big(0)
    for a in values(v);total=addchecked(total,mulchecked(a,a,b),b);end
    total//mulchecked(r.denominator,r.denominator,b)
end
end
