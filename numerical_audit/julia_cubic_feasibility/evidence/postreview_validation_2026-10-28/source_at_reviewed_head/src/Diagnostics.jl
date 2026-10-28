module Diagnostics
using ..ResponseCore, ..StructureDP
export select_window, structural_audit

function select_window(units,order,hf,width;guard=()->nothing)
    2<=width<=min(9,length(order)-2) || error("invalid interior window width")
    centered=center_units(units,hf); cache=AmplitudeCache(centered)
    beststart=2; bestcount=-1
    for start in 2:length(order)-width
        guard(); labels=order[start:start+width-1]; active=0
        for i in 1:width-1, j in i+1:width
            a,b=labels[i],labels[j]
            # Reference response only; physical errors never enter selection.
            nonzero=false
            for local_order in ((a,b),(b,a))
                col=atom_column(centered,local_order,hf;center=false,cache=cache)
                any(d!=hf && v!=0 for (d,v) in col) && (nonzero=true)
            end
            active+=nonzero
        end
        if active>bestcount; bestcount=active;beststart=start; end
    end
    (positions=beststart:beststart+width-1,active_pairs=bestcount,
     selection="maximum computed nonzero internal cubic pair scopes, earliest interior start on ties")
end

function structural_audit(units,order,hf;max_scopes=150_000,max_atoms=75_000,
                          max_pair_events=1_000_000,max_states=100_000,guard=()->nothing)
    m=length(units); centered=center_units(units,hf); cache=AmplitudeCache(centered)
    buckets=Dict{UInt64,Vector{NTuple{3,Int}}}()
    scopes=0; nonzero=0; stopped=false; pairscopes=0; triplescopes=0
    function record(scope,variants)
        guard(); destinations=Set{UInt64}()
        for local_order in variants
            col=atom_column(centered,local_order,hf;center=false,cache=cache)
            for (d,v) in col; d!=hf && v!=0 && push!(destinations,d); end
        end
        scopes+=1
        if !isempty(destinations)
            nonzero+=1
            for dest in destinations
                push!(get!(buckets,dest,NTuple{3,Int}[]),scope)
            end
        end
        stopped=scopes>=max_scopes || nonzero>=max_atoms
    end
    for a in 1:m-1,b in a+1:m
        record((a,b,0),((a,b),(b,a))); pairscopes+=1
        stopped && break
    end
    if !stopped
        for a in 1:m-2,b in a+1:m-1,c in b+1:m
            record((a,b,c),((a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a)))
            triplescopes+=1; stopped && break
        end
    end
    bucket_iterator=(d=>(Tuple(x for x in scope if x!=0) for scope in buckets[d]) for d in sort!(collect(keys(buckets))))
    audit=audit_destination_buckets(collect(1:m),bucket_iterator;
        max_atom_visits=2max_atoms,max_bucket_atoms=max_atoms,max_edges=100_000,
        max_pair_events=max_pair_events,max_edge_visits=2_000_000,guard=guard)
    decomposition=elimination_decomposition(collect(1:m),audit.edges;guard=guard)
    width=decomposition.width_upper_bound
    largest=decomposition.complete ? maximum(length, decomposition.bags) : 0
    states=decomposition.complete ? sum(factorial(big(length(b))) for b in decomposition.bags) : big(0)
    edge_scopes=((scope=(a,b),) for (a,b) in audit.edges)
    valid=decomposition.complete ? validate_decomposition(collect(1:m),edge_scopes,decomposition.bags,decomposition.tree_edges) : nothing
    # A small checked decomposition of an INCOMPLETE graph would not validate
    # full DP feasibility. Large widths here do not prove minimum-width bounds.
    go=!stopped && audit.complete && decomposition.complete && states<=max_states && audit.full_pair_enumeration_admissible
    Dict("status"=>stopped ? "scope_or_atom_cap" : "complete",
      "scope_count"=>scopes,"nonzero_computed_atoms"=>nonzero,"pairs_visited"=>pairscopes,
      "triples_visited"=>triplescopes,"total_pair_scopes"=>binomial(big(m),2),
      "total_triple_scopes"=>binomial(big(m),3),"destinations"=>length(buckets),
      "maximum_bucket_atoms"=>audit.maximum_bucket_atoms,
      "maximum_bucket_label_union"=>audit.maximum_bucket_label_union,
      "potential_pair_events_retained"=>audit.potential_pair_events,
      "potential_edges_retained"=>length(audit.edges),"bucket_status"=>audit.status,
      "decomposition_complete"=>decomposition.complete,"decomposition_checked"=>valid,
      "width_upper_bound_retained_potential_graph"=>width,
      "largest_bag"=>largest,"factorial_states_retained_potential_graph"=>states,
      "full_dp_go"=>go,"state_cap"=>max_states,"pair_event_cap"=>max_pair_events,
      "scope_complete"=>!stopped,
      "interpretation"=>"Floating support diagnostic; incomplete scans and potential-supergraph width do not certify the exact coherent graph or intrinsic DP impossibility.")
end
end
