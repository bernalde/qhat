module StructureDP

export LocalFactor, audit_factors, audit_destination_buckets,
       elimination_decomposition, validate_decomposition, bag_order_dp,
       evaluate_factors

const ExactCost = Rational{BigInt}
const Order = Tuple{Vararg{Int}}

exact_cost(x::Integer) = BigInt(x) // BigInt(1)
exact_cost(x::Rational) = BigInt(numerator(x)) // BigInt(denominator(x))
function exact_cost(x::AbstractFloat)
    isfinite(x) || throw(ArgumentError("factor costs must be finite"))
    return Rational{BigInt}(x) # Exact dyadic value, not tolerance-based rationalize.
end
exact_cost(x::Real) = throw(ArgumentError("use integer, rational, or finite binary floating costs"))

"""A complete local relative-order table; finite floats are interpreted exactly.

Only the supplied scalar tables are optimized. Their relation to a physical
response, and any error from rounding that response, belong to the caller.
"""
struct LocalFactor
    scope::Vector{Int}
    values::Dict{Order,ExactCost}
    function LocalFactor(scope, values::AbstractDict; max_entries::Int=720)
        labels = Int[x for x in scope]
        length(unique(labels)) == length(labels) || throw(ArgumentError("duplicate factor label"))
        max_entries >= 1 || throw(ArgumentError("max_entries must be positive"))
        required = factorial(big(length(labels)))
        required <= max_entries || throw(ArgumentError("factor table exceeds entry cap"))
        length(values) == required || throw(ArgumentError("incomplete local-order table"))
        table = Dict{Order,ExactCost}()
        wanted = Set(labels)
        for (key, value) in values
            order = Tuple(Int[x for x in key])
            length(order) == length(labels) && Set(order) == wanted ||
                throw(ArgumentError("invalid local-order table key"))
            table[order] = exact_cost(value)
        end
        # Remove labels whose insertion position never changes this factor.
        # This also reduces constants to the empty scope.
        changed = true
        while changed && !isempty(labels)
            changed = false
            for label in labels
                reduced = Dict{Order,ExactCost}()
                independent = true
                for (order, cost) in table
                    key = Tuple(x for x in order if x != label)
                    if haskey(reduced, key) && reduced[key] != cost
                        independent = false
                        break
                    end
                    reduced[key] = cost
                end
                if independent
                    labels = [x for x in labels if x != label]
                    table = reduced
                    changed = true
                    break
                end
            end
        end
        new(labels, table)
    end
end

_guard(guard) = guard === nothing ? nothing : guard()
_edge(a::Int,b::Int) = a < b ? (a,b) : (b,a)
_restrict(order, scope) = Tuple(x for x in order if x in scope)
_bits(x::Integer) = ndigits(abs(x); base=2)

function _labels(labels)
    result = Int[x for x in labels]
    length(unique(result)) == length(result) || throw(ArgumentError("duplicate label"))
    return result
end

function _refusal(status::Symbol, reason::String; details=nothing)
    return (status=status, complete=false, reason=reason, order=Int[], cost=nothing,
            details=details)
end

"""Stream actual nonconstant factor scopes into their checked primal graph.

The clique bound certifies this supplied factor graph, not every possible
refactorization of the objective. No tolerance-based zero classification is used.
"""
function audit_factors(labels, factors; max_factors::Int=100_000,
                       max_edges::Int=100_000, max_entries::Int=1_000_000,
                       guard=nothing)
    vertices = _labels(labels); allowed = Set(vertices)
    edges = Set{Tuple{Int,Int}}(); clique = Int[]
    count = 0; entries = 0; status = :complete
    for factor in factors
        _guard(guard)
        count += 1
        if count > max_factors
            count -= 1; status = :factor_limit; break
        end
        all(x -> x in allowed, factor.scope) || throw(ArgumentError("factor has unknown label"))
        if entries + length(factor.values) > max_entries
            count -= 1; status = :entry_limit; break
        end
        entries += length(factor.values)
        length(factor.scope) > length(clique) && (clique = copy(factor.scope))
        for i in eachindex(factor.scope), j in 1:i-1
            e = _edge(factor.scope[i], factor.scope[j])
            if !(e in edges) && length(edges) >= max_edges
                status = :edge_limit; break
            end
            push!(edges,e)
        end
        status == :complete || break
    end
    return (status=status, complete=status==:complete, graph_kind=:exact_factor_graph,
            labels=vertices, edges=edges, factors_visited=count, entries_visited=entries,
            clique_witness=clique, clique_width_lower_bound=max(0,length(clique)-1),
            lower_bound_scope="supplied exact scalar-factor primal graph")
end

"""Stream destination => atom-scope iterators without enumerating atom pairs.

A destination's union of labels forms the same potential primal clique as all
pairwise scope unions. Sharing support is only a necessary overlap condition;
signed Gram contributions can vanish or cancel. This graph provides NO clique
lower bound for the exact coherent graph. Duplicate atom entries overestimate
pair counts and should be removed by the producer when exact counts are wanted.
"""
function audit_destination_buckets(labels, buckets;
        max_buckets::Int=100_000, max_atom_visits::Int=1_000_000,
        max_bucket_atoms::Int=100_000, max_edges::Int=100_000,
        max_pair_events::Integer=10_000_000, max_edge_visits::Int=1_000_000,
        guard=nothing)
    vertices = _labels(labels); allowed = Set(vertices)
    edges = Set{Tuple{Int,Int}}()
    nb = 0; visits = 0; max_bucket = 0; max_union = 0; edge_visits=0
    seen_destinations=Set{Any}()
    pairs = big(0); status = :complete; completed_buckets = 0
    for (destination, scopes) in buckets
        _guard(guard)
        if nb >= max_buckets
            status=:bucket_limit; break
        end
        destination in seen_destinations && throw(ArgumentError("destination bucket repeated; merge its scopes first"))
        push!(seen_destinations,destination)
        nb += 1; n=0; union_labels=Set{Int}()
        for scope in scopes
            _guard(guard)
            if visits >= max_atom_visits || n >= max_bucket_atoms
                status = visits >= max_atom_visits ? :atom_visit_limit : :bucket_atom_limit
                break
            end
            for label in scope
                label in allowed || throw(ArgumentError("bucket has unknown label"))
                push!(union_labels,label)
            end
            visits += 1; n += 1
        end
        max_bucket=max(max_bucket,n); max_union=max(max_union,length(union_labels))
        status == :complete || break
        pairs += big(n)*(n-1)÷2
        # Pair-event cap reports downstream coefficient work refusal. The
        # compressed graph remains cheap, so finish it under its own edge cap.
        scope=sort!(collect(union_labels))
        for i in eachindex(scope), j in 1:i-1
            if edge_visits>=max_edge_visits
                status=:edge_work_limit;break
            end
            edge_visits+=1
            e=_edge(scope[i],scope[j])
            if !(e in edges) && length(edges)>=max_edges
                status=:edge_limit; break
            end
            push!(edges,e)
        end
        status == :complete || break
        completed_buckets += 1
    end
    return (status=status, complete=status==:complete,
            graph_kind=:support_overlap_potential_supergraph,
            labels=vertices, edges=edges, buckets_visited=nb,
            completed_buckets=completed_buckets, atom_visits=visits,
            edge_visits=edge_visits,
            maximum_bucket_atoms=max_bucket, maximum_bucket_label_union=max_union,
            potential_pair_events=pairs, pair_count_complete=status==:complete,
            full_pair_enumeration_admissible=status==:complete && pairs<=max_pair_events,
            pair_enumeration_status=status!=:complete ? :incomplete :
                (pairs<=max_pair_events ? :within_cap : :pair_event_limit),
            clique_width_lower_bound=nothing)
end

"""Bounded deterministic minimum-degree elimination; width is an upper bound.

A large heuristic width or a cap refusal proves nothing about minimum width.
"""
function elimination_decomposition(labels, edges; max_work::Int=10_000_000,
        max_fill_edges::Int=200_000, guard=nothing)
    vertices=_labels(labels); allowed=Set(vertices)
    adj=Dict(v=>Set{Int}() for v in vertices)
    edge_count=0
    for (a,b) in edges
        a in allowed && b in allowed && a!=b || throw(ArgumentError("invalid graph edge"))
        if !(b in adj[a])
            edge_count += 1
            edge_count<=max_fill_edges || return (status=:edge_limit,complete=false,
                bags=Vector{Vector{Int}}(), tree_edges=Tuple{Int,Int}[], width_upper_bound=nothing)
            push!(adj[a],b);push!(adj[b],a)
        end
    end
    remaining=Set(vertices); bags=Vector{Vector{Int}}(); removed=Int[]; work=0
    while !isempty(remaining)
        _guard(guard)
        work += length(remaining)
        work<=max_work || return (status=:work_limit,complete=false,bags=bags,
            tree_edges=Tuple{Int,Int}[],width_upper_bound=nothing)
        v=first(remaining)
        for candidate in remaining
            (length(adj[candidate]),candidate)<(length(adj[v]),v) && (v=candidate)
        end
        neighbors=sort!(collect(adj[v])); push!(bags,sort!([v;neighbors]));push!(removed,v)
        for i in eachindex(neighbors), j in 1:i-1
            work+=1
            work<=max_work || return (status=:work_limit,complete=false,bags=bags,
                tree_edges=Tuple{Int,Int}[],width_upper_bound=nothing)
            a,b=neighbors[i],neighbors[j]
            if !(b in adj[a])
                edge_count+=1
                edge_count<=max_fill_edges || return (status=:fill_edge_limit,complete=false,
                    bags=bags,tree_edges=Tuple{Int,Int}[],width_upper_bound=nothing)
                push!(adj[a],b);push!(adj[b],a)
            end
        end
        for u in neighbors
            delete!(adj[u],v)
        end
        edge_count-=length(neighbors);empty!(adj[v]);delete!(remaining,v)
    end
    if isempty(bags)
        push!(bags,Int[])
    end
    position=Dict(v=>i for (i,v) in enumerate(removed)); tree=Tuple{Int,Int}[];roots=Int[]
    for (i,v) in enumerate(removed)
        later=[position[u] for u in bags[i] if u!=v]
        isempty(later) ? push!(roots,i) : push!(tree,(i,minimum(later)))
    end
    for i in 2:length(roots)
        push!(tree,(roots[i-1],roots[i]))
    end
    return (status=:complete,complete=true,bags=bags,tree_edges=tree,
            width_upper_bound=maximum(length.(bags))-1, work=work,
            width_is_heuristic_upper_bound=true)
end

"""Check tree shape, label coverage, running intersection, and every factor scope."""
function validate_decomposition(labels, factors, bags, tree_edges)
    vertices=_labels(labels);allowed=Set(vertices); t=length(bags)
    t>=1 || throw(ArgumentError("decomposition must contain at least one bag"))
    bagsets=Set{Int}[]
    for bag in bags
        length(unique(bag))==length(bag) || throw(ArgumentError("duplicate bag label"))
        all(x->x in allowed,bag) || throw(ArgumentError("unknown bag label"))
        push!(bagsets,Set(Int[x for x in bag]))
    end
    length(tree_edges)==t-1 || throw(ArgumentError("bag edges must form a tree"))
    tree=[Int[] for _ in 1:t];seen_edges=Set{Tuple{Int,Int}}()
    for (a,b) in tree_edges
        1<=a<=t && 1<=b<=t && a!=b || throw(ArgumentError("invalid bag edge"))
        e=_edge(a,b);e in seen_edges && throw(ArgumentError("duplicate bag edge"))
        push!(seen_edges,e);push!(tree[a],b);push!(tree[b],a)
    end
    visited=Set([1]);queue=[1]
    for node in queue, child in tree[node]
        if !(child in visited)
            push!(visited,child);push!(queue,child)
        end
    end
    length(visited)==t || throw(ArgumentError("disconnected/cyclic bag tree"))
    for label in vertices
        containing=[i for i in 1:t if label in bagsets[i]]
        isempty(containing) && throw(ArgumentError("uncovered label"))
        chosen=Set(containing);reached=Set([first(containing)]);front=[first(containing)]
        for node in front, child in tree[node]
            if child in chosen && !(child in reached)
                push!(reached,child);push!(front,child)
            end
        end
        length(reached)==length(containing) || throw(ArgumentError("running intersection fails"))
    end
    for factor in factors
        all(x->x in allowed,factor.scope) || throw(ArgumentError("factor has unknown label"))
        any(s->issubset(Set(factor.scope),s),bagsets) ||
            throw(ArgumentError("uncovered factor scope"))
    end
    return true
end

function _orders(labels::Vector{Int})
    out=Order[];prefix=Int[];remaining=sort(labels)
    function visit(rest)
        if isempty(rest)
            push!(out,Tuple(prefix));return
        end
        for i in eachindex(rest)
            push!(prefix,rest[i]);visit([rest[j] for j in eachindex(rest) if j!=i]);pop!(prefix)
        end
    end
    visit(remaining);return out
end

function _condition(factors, first_label, last_label)
    result=LocalFactor[]
    for factor in factors
        scope=[v for v in factor.scope if v!=first_label && v!=last_label]
        table=Dict{Order,ExactCost}()
        for order in _orders(scope)
            original=Int[]
            first_label in factor.scope && push!(original,first_label)
            append!(original,order)
            last_label in factor.scope && push!(original,last_label)
            table[order]=factor.values[Tuple(original)]
        end
        push!(result,LocalFactor(scope,table))
    end
    return result
end

function evaluate_factors(factors, order)
    length(unique(order))==length(order) || throw(ArgumentError("duplicate order label"))
    cost=big(0)//big(1)
    for factor in factors
        key=_restrict(order,factor.scope)
        haskey(factor.values,key) || throw(ArgumentError("order misses factor label"))
        cost+=factor.values[key]
    end
    return cost
end

function _merge_orders(left::Vector{Int}, right::Vector{Int})
    common=Set(intersect(left,right))
    _restrict(left,common)==_restrict(right,common) || error("inconsistent separator orders")
    result=Int[];i=1;j=1
    for shared in _restrict(left,common)
        while left[i]!=shared
            push!(result,left[i]);i+=1
        end
        while right[j]!=shared
            push!(result,right[j]);j+=1
        end
        push!(result,shared);i+=1;j+=1
    end
    append!(result,left[i:end]);append!(result,right[j:end]);return result
end

"""Exactly optimize supplied tables over global permutations using bag orders.

Caps are checked before any bag permutation/state allocation. Memory/work numbers
are conservative implementation estimates, not a process-RSS guarantee; the root
driver must retain its OS resource guard. Endpoints refer to original labels.
"""
function bag_order_dp(labels, factors, bags, tree_edges;
        first_label=nothing,last_label=nothing,max_states::Int=100_000,
        max_bag_states::Int=50_000,max_work::Int=10_000_000,
        max_factor_entries::Int=1_000_000,max_bytes::Int=64*1024*1024,
        max_cost_bits::Int=16_384,guard=nothing)
    vertices=_labels(labels)
    first_label===nothing || first_label in vertices || throw(ArgumentError("unknown first endpoint"))
    last_label===nothing || last_label in vertices || throw(ArgumentError("unknown last endpoint"))
    first_label!==nothing && first_label==last_label && length(vertices)!=1 &&
        throw(ArgumentError("endpoints must differ"))
    fs=LocalFactor[]
    entries=big(0);denom_bits=big(0);num_bits=1
    for factor in factors
        entries+=length(factor.values)
        entries<=max_factor_entries || return _refusal(:factor_entry_limit,"factor entry cap exceeded")
        push!(fs,factor)
        for cost in values(factor.values)
            denom_bits+=_bits(denominator(cost));num_bits=max(num_bits,_bits(numerator(cost)))
        end
    end
    cost_bits=big(num_bits)+denom_bits+_bits(max(1,length(fs)))+2
    cost_bits<=max_cost_bits || return _refusal(:cost_bit_limit,"rational arithmetic bit bound exceeded";
                                               details=(cost_bits=cost_bits,))
    validate_decomposition(vertices,fs,bags,tree_edges)
    if length(vertices)==1 && first_label==last_label && first_label!==nothing
        return (status=:complete,complete=true,reason="",order=copy(vertices),
            cost=evaluate_factors(fs,vertices),details=(states=big(1),encoded_tables_exact=true))
    end
    interior=[v for v in vertices if v!=first_label && v!=last_label]
    bs=[Int[v for v in bag if v!=first_label && v!=last_label] for bag in bags]
    counts=[factorial(big(length(bag))) for bag in bs]
    total_states=sum(counts)
    maximum(counts)<=max_bag_states || return _refusal(:bag_state_limit,"one bag exceeds factorial state cap";
        details=(states=total_states,maximum_bag_states=maximum(counts)))
    total_states<=max_states || return _refusal(:state_limit,"total factorial states exceed cap";
        details=(states=total_states,))
    t=length(bs);tree=[Int[] for _ in 1:t]
    for (a,b) in tree_edges
        push!(tree[a],b);push!(tree[b],a)
    end
    parent=zeros(Int,t);traversal=[1];parent[1]=-1
    for node in traversal,child in tree[node]
        if child!=parent[node]
            parent[child]=node;push!(traversal,child)
        end
    end
    assigned=[Int[] for _ in 1:t]
    for (fi,factor) in enumerate(fs)
        reduced_scope=[v for v in factor.scope if v!=first_label && v!=last_label]
        node=findfirst(bag->all(v->v in bag,reduced_scope),bs)
        push!(assigned[node],fi)
    end
    work=entries
    # Charge copied conditioned tables as well as state/message dictionaries.
    bytes=entries*(256+4*cld(cost_bits,8))
    for node in 1:t
        children=count(child->parent[child]==node,tree[node])
        work+=counts[node]*(1+length(assigned[node])*7+children*(length(bs[node])+1))
        bytes+=counts[node]*(512+32*length(bs[node])+4*cld(cost_bits,8))
    end
    work<=max_work || return _refusal(:work_limit,"estimated DP work exceeds cap";
        details=(states=total_states,work=work,estimated_bytes=bytes))
    bytes<=max_bytes || return _refusal(:memory_limit,"estimated state storage exceeds cap";
        details=(states=total_states,work=work,estimated_bytes=bytes))
    # Scope conditioning allocates at most the supplied factor-entry count and
    # is postponed until the state/work/storage admission checks have passed.
    cfs=_condition(fs,first_label,last_label)
    validate_decomposition(interior,cfs,bs,tree_edges)
    _guard(guard)
    orders=[_orders(bag) for bag in bs]
    costs=[Vector{ExactCost}(undef,length(orders[i])) for i in 1:t]
    messages=[Dict{Order,Int}() for _ in 1:t]
    for node in reverse(traversal)
        separator=parent[node]<1 ? Int[] : intersect(bs[node],bs[parent[node]])
        for (oi,order) in enumerate(orders[node])
            _guard(guard)
            value=big(0)//big(1)
            for fi in assigned[node]
                factor=cfs[fi];value+=factor.values[_restrict(order,factor.scope)]
            end
            for child in tree[node]
                parent[child]==node || continue
                key=_restrict(order,intersect(bs[node],bs[child]))
                value+=costs[child][messages[child][key]]
            end
            costs[node][oi]=value
            key=_restrict(order,separator)
            if !haskey(messages[node],key) || value<costs[node][messages[node][key]]
                messages[node][key]=oi
            end
        end
    end
    function recover(node,oi)
        order=collect(orders[node][oi])
        for child in tree[node]
            parent[child]==node || continue
            key=_restrict(order,intersect(bs[node],bs[child]))
            order=_merge_orders(order,recover(child,messages[child][key]))
        end
        return order
    end
    chosen=messages[1][()];order=recover(1,chosen)
    first_label===nothing || pushfirst!(order,first_label)
    last_label===nothing || push!(order,last_label)
    length(order)==length(vertices) && Set(order)==Set(vertices) || error("DP witness is not a permutation")
    actual=evaluate_factors(fs,order)
    actual==costs[1][chosen] || error("DP traceback changed factor cost")
    return (status=:complete,complete=true,reason="",order=order,cost=actual,
            details=(states=total_states,work=work,estimated_bytes=bytes,
                     cost_bits=cost_bits,encoded_tables_exact=true,
                     physical_response_exactness_claim=false))
end

end # module
