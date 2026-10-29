module WindowStudyIO
using TOML, SHA
using ..ResponseCore: Unit
export load_model, write_result, cpu_seconds

cpu_seconds() = ccall(:clock, Int64, ()) / 1.0e6

function load_model(path::AbstractString; scale::Rational=1//1)
    # input_sha256 authenticates the serialized coefficient target. A scaled
    # target must therefore be written to distinct bytes before it is loaded.
    scale==1//1 || throw(ArgumentError("nonunit scale requires a distinct serialized model input"))
    filesize(path) <= 8*1024^2 || error("model input too large")
    data = TOML.parsefile(path)
    n = Int(data["n"])
    1 <= n <= 20 || error("pilot supports 1:20 qubits")
    denstr = data["denominator"]
    denstr isa AbstractString && ncodeunits(denstr) <= 1280 || error("denominator spelling bound")
    den = parse(BigInt, denstr)
    den > 0 && ndigits(den; base=2) <= 4096 || error("bounded positive denominator required")
    units = Unit[]
    partition = get(data,"partition","conserving")
    partition in ("conserving","pauli") || error("unknown study partition")
    max_units = partition=="pauli" ? 8192 : 1024
    1 <= length(data["units"]) <= max_units || error("unit count before allocation")
    sum(length(block["words"]) for block in data["units"]) <= 8192 || error("word count before allocation")
    for block in data["units"]
        length(block["words"]) == length(block["numerators"]) || error("word count mismatch")
        partition=="pauli" && length(block["words"])!=1 && error("Pauli partition requires singleton words")
        pars = UInt64[]; signs = Int8[]; coeff = Float64[]
        flip = nothing
        for (word, numerator) in zip(block["words"], block["numerators"])
            word isa AbstractString && ncodeunits(word) <= 4n || error("Pauli word spelling bound")
            numerator isa AbstractString && ncodeunits(numerator) <= 1280 || error("numerator spelling bound")
            integer = parse(BigInt, numerator)
            ndigits(abs(integer); base=2) <= 4096 || error("coefficient numerator bit bound")
            x = z = UInt64(0); ys = 0; seen = Set{Int}()
            for entry in (isempty(word) ? String[] : split(word, ','))
                q = parse(Int, entry[1:end-1]); p = entry[end]
                0 <= q < n && !(q in seen) || error("invalid Pauli site")
                push!(seen,q)
                p in ('X','Y','Z') || error("invalid Pauli letter")
                p in ('X','Y') && (x |= UInt64(1)<<q)
                p in ('Z','Y') && (z |= UInt64(1)<<q)
                ys += p == 'Y'
            end
            iseven(ys) || error("pilot requires real equal-flip units")
            flip === nothing ? (flip=x) : (flip == x || error("mixed flip unit"))
            c = Float64((integer//den)*scale)
            isfinite(c) && c != 0 || error("invalid coefficient")
            push!(pars,z); push!(signs,isodd(ys÷2) ? Int8(-1) : Int8(1)); push!(coeff,c)
        end
        isempty(coeff) && error("empty unit")
        # Same-flip Pauli words commute iff the symplectic parity vanishes.
        for i in eachindex(pars), j in 1:i-1
            iseven(count_ones(flip & (pars[i] ⊻ pars[j]))) || error("noncommuting unit")
        end
        push!(units,Unit(flip,pars,signs,coeff))
    end
    m = length(units)
    1 <= m <= max_units || error("unit count outside study cap")
    order = Int.(data[haskey(data,"baseline_order") ? "baseline_order" : "cubic_order"])
    sort(order)==collect(1:m) || error("invalid baseline permutation")
    electrons = Int(data["electrons"])
    0 <= electrons <= n || error("invalid electron count")
    hf = (UInt64(1)<<electrons)-1
    (units=units, order=order, hf=hf, n=n, electrons=electrons,
     case=String(data["case"]), words=sum(length(u.coeff) for u in units),
     partition=partition, source=data, input_sha256=bytes2hex(sha256(read(path))))
end

function clean(value)
    value === nothing && return "unavailable"
    value isa Symbol && return String(value)
    value isa UInt64 && return string(value)
    value isa BigInt && return string(value)
    value isa Rational && return string(value)
    value isa NamedTuple && return Dict(string(k)=>clean(v) for (k,v) in pairs(value))
    value isa AbstractDict && return Dict(string(k)=>clean(v) for (k,v) in value)
    value isa Tuple && return [clean(v) for v in value]
    value isa AbstractVector && return [clean(v) for v in value]
    value isa AbstractFloat && !isfinite(value) && return string(value)
    value
end

function write_result(path, value)
    ispath(path) && error("refuse to overwrite output $path")
    mkpath(dirname(path))
    open(path,"w") do io
        TOML.print(io,clean(value);sorted=true)
    end
end
end
