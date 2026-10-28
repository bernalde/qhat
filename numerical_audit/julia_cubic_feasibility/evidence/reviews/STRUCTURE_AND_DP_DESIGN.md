# Structure audit and exact bag-order DP: implementation handoff

Owner: `/root/approximation_obstructions`; coordinator is the sole integrator and execution owner. Base supplied by coordinator: `9a19a118`. Mutable scope is confined to the new Julia study's `src/StructureDP.jl`, `test/structure_tests.jl`, and this coordination note. No Git/index/remote writes, installs, helpers, tests, or numerical commands were executed by this worker.

## Implemented API

`StructureDP.jl` uses only Julia Base. Unit labels are arbitrary distinct `Int`s.

* `LocalFactor(scope, values; max_entries=720)` takes a complete dictionary from scope permutations to signed real scalar costs. Integer/rational costs are converted to `Rational{BigInt}`; finite binary floats become their **exact supplied dyadic values**, avoiding overflow and floating accumulation ambiguities. The constructor strips a label only if every insertion position leaves the exact table unchanged, reducing constants to empty scopes. It does not infer physical zeros from a numerical tolerance.
* `audit_factors(labels, factors; max_factors, max_edges, max_entries, guard)` streams actual nonconstant scalar scopes and their clique edges. It returns exact/partial status, visited counts, and a scope clique witness with a valid treewidth lower bound **for this supplied factorization's exact primal graph**. This does not prove that every equivalent reformulation has the same width.
* `audit_destination_buckets(labels, buckets; max_buckets, max_atom_visits, max_bucket_atoms, max_edges, max_edge_visits, max_pair_events, guard)` accepts `destination => iterator_of_atom_scopes`. Each atom should occur once in each destination bucket if the pair count is to count distinct atoms. A bucket's union of labels is made a potential clique. This is equivalent to forming clique edges from all pairs of scopes in that bucket, without enumerating pairs. The routine records `sum choose(bucket_atoms,2)` using BigInt and a separate downstream-pair-enumeration admissibility flag. It can finish a compressed potential graph even when explicit Gram-pair work exceeds that cap. Atom/edge-work caps stop and return a partial graph explicitly. Repeated destination keys are rejected: splitting one destination into chunks without merging would miss cross-chunk edges.
* `elimination_decomposition(labels, edges; max_work, max_fill_edges, guard)` produces a bounded deterministic minimum-degree elimination decomposition. Disconnected components are joined by empty-separator tree edges. Its width is a **heuristic upper bound**, never a lower bound. A cap refusal or large result does not prove minimum-width infeasibility.
* `validate_decomposition(labels, factors, bags, tree_edges)` checks a connected tree, no duplicate/unknown bag labels, label coverage, running intersection, and scope coverage. It rejects cyclic fake-order decompositions before optimization.
* `bag_order_dp(labels, factors, bags, tree_edges; first_label=nothing, last_label=nothing, max_states=100000, max_bag_states=50000, max_work=10000000, max_factor_entries=1000000, max_bytes=64MiB, max_cost_bits=16384, guard=nothing)` solves the supplied exact tables with total orders on each bag, separator minima, and witness gluing. It returns a named tuple with `status`, `complete`, `order`, `cost`, and resource details. Caps refuse before scope-conditioning or bag-permutation/state allocation. Input/decomposition arrays already supplied by the caller are outside this allocation promise. Rational arithmetic bit growth and state/message storage are conservatively estimated; the coordinator's OS guard remains necessary.
* `evaluate_factors(factors, order)` independently sums the encoded scalar tables and checks the recovered DP witness's value.

The DP's endpoint keywords are **`first_label` and `last_label`**. It does not calculate rotation budgets. The root driver fixes feasible endpoint pairs/r, calls the DP on a constant denominator, and compares the resulting Q/r^4. Endpoints are removed from bag states and conditioned within each local cost. No budget feasibility is invented here.

## Correct interpretation of structural evidence

There are three distinct objects:

1. The unit noncommutation graph limits which cubic atoms may occur.
2. Destination overlap yields a **potential coherent supergraph**: an overlap permits interference but does not establish a nonzero Gram factor. A bucket clique can be large even when exact scalar cross terms cancel. Its size is not an actual-graph clique lower bound.
3. Complete, exact nonconstant scalar factor tables give an **actual supplied-factor primal graph**. A k-label scope provides a certified k-clique and width lower bound k-1 in this graph.

Only the third object's clique witness is reported as a lower bound. A potential graph decomposition can safely serve as an upper-bound decomposition for actual factors if all actual scopes are covered; a large potential width merely says that this sufficient representation did not pass the budget. Likewise the minimum-degree heuristic's width can be poor. Report `full_dp_not_admitted_under_declared_representation_and_caps`, not mathematical impossibility.

Constructing every scalar cross factor is not authorized: generic all-atom pair enumeration can be O(m^6). The root can use the compressed bucket graph and a declared state budget to stop before any such expansion. If an actual tiny-factor subproblem is available, `audit_factors` adds genuine clique evidence cheaply. The full campaign must preserve incomplete bucket counts/refusals rather than presenting a partial graph as exact.

## ResponseCore coordination

The ResponseCore worker reports a conditioned cubic `WindowModel` with full baseline response and grouped two-/three-window-label vector tables. Outside atoms remain in that baseline and therefore interfere coherently with window-dependent terms. This is the correct model for a fixed-outside cubic window; it is not the historical first-order JD window calculation.

For a supplied factor-table DP, scalar factors can be made from the baseline cross terms and pairwise vector inner products only after a bounded window has made that construction affordable. For structural scouting, destination bucket entries can instead use atom-scope unions across all variants; the resulting graph must keep its potential label. All transported/projected variants and the full coherent reference convention belong to ResponseCore, not StructureDP.

If the full graph route fails admission, the coordinator owns selection of a bounded conditioned window, its exhaustive enumeration budget, and rescoring with the full cubic response. This worker does not launch a fallback optimizer or change any benchmark protocol.

## Tests prepared, not executed

`test/structure_tests.jl` has meaningful tiny checks for:

* compatible three-label bag orders gluing to a global optimum with signed costs;
* rejection of a three-bag cyclic fake decomposition;
* a genuine cyclic pair-cost instance solved in a valid bag;
* first/last conditioning and exact dyadic float interpretation;
* disconnected components, isolated labels, empty and singleton cases;
* exact pruning of irrelevant labels/constants;
* actual versus potential graph evidence and separate pair-work refusal;
* atom/edge caps and repeated destination rejection;
* factorial, memory, work, factor-entry, and arithmetic-bit refusal before state allocation;
* a capped elimination heuristic reporting incomplete status rather than a width claim.

Coordinator should run this file only under the already established serialized quarter-core/512MiB guard, then integrate any failures. No passing result is claimed by this note.

## Scientific scope

An exact DP result is exact for the supplied rational/dyadic scalar tables. A floating ResponseCore table remains an approximation to the transcendental-kernel mathematical response; this module supplies no interval/physical-error certificate. No large-molecule feasibility follows from the tiny tests, and no conclusion about CH4/CO is available before the coordinator's guarded structural study and retained status metadata.
