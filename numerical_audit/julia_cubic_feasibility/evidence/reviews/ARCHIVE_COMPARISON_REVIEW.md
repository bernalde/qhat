# Archived nine-unit DP comparison: independent design review

Read-only scientific review, 2026-10-28. Inspected the archived protocol, all three producer certificates, domain records, `exact.py`, and the independent DFS/Lehmer verifier under `/tmp/qhat-scalable-ordering-2026-10-27/worktree/numerical_audit/exact_response_windows`. No experiment, test, source edit, or nested delegation was performed by this reviewer. This coordination note is the only write.

## Smallest sound comparison

Use the three already-frozen nine-label domains without regenerating chemistry, selecting new windows, recomputing sinc weights, or running physical propagation. First verify import/normalization on the stored warm order and stored optimum. Then compare a fresh Julia exact enumeration and the Julia bag-order DP on precisely the same integer objective, with a checked decomposition and resource refusals retained. All three cases are small enough to make independent exhaustive correctness evidence meaningful; success here establishes this restricted-domain implementation, not full-order molecular scalability or the cubic Strang result.

The archived target is the first-response objective J_D, not the cubic Strang response currently implemented in ResponseCore. Reusing the same molecule/order does not make those objectives equivalent. Do not invoke cubic columns, recompute F_D into complex Float64 amplitudes, take floating square roots of supplied weights, divide by r^4, or change the unit partition to make the adapter convenient.

## Invariants and exact score map

For each archived domain preserve byte-authenticated integer coefficient numerators, positive common coefficient denominator d, original unit labels/partition, HF occupation, T=1, full warm order, zero-based offset, selected nine-label set, outside prefix/suffix and their positions, all pair response entries, supplied dyadic weight ratios, integer weights, and objective denominator. Archive labels are zero based; a monotone +1 map to Julia IDs is safe if applied everywhere and inverted for output. A baseline-rank map is not a global-label lexicographic map.

For an order, the domain pair record [i,j,w,b] contributes +b when i precedes j and -b otherwise, using the archived pair convention. Preserve that convention consistently. Define

    b_w(order)=sum of signed archived pair numerators,
    N(order)=sum_w W_w b_w(order)^2,
    H=4*d^4*weight_scale,
    J_D(order)=N(order)/H.

Each W_w is the archived exact integer weight. All three stored weight scales are 9007199254740992=2^53. Import these numbers and H as BigInt strings, never through Float64. Stored displayed `initial_JD`/`incumbent_JD` values are summaries; the correctness oracle is the exact numerator and denominator. Do not replace stored rho values by a fresh sin/sinc call: a last-bit change defines a different optimization problem.

The frozen window is contiguous. Therefore all pairs with fewer than two window labels have fixed orientation. Let c_w contain every such pair contribution, and let s_e(order) be +/-1 for an internal pair. The exact reduced objective is

    N(order)=sum_w W_w (c_w+sum_{e internal,dest(e)=w} b_e s_e(order))^2.

The fixed outside amplitude stays INSIDE the square. Dropping it, or optimizing only internal pair norms, loses interference and changes the archived optimum.

A scalar factor expansion is exact:

    C=sum_w W_w c_w^2 + sum_e W_dest(e) b_e^2,
    pair terms 2 W_dest(e) c_dest(e) b_e s_e,
    pair-pair terms 2 W_w b_e b_f s_e s_f for e<f with shared w.

Its nonconstant factor scopes use at most four unit labels, and entries can be signed BigInts. Aggregating equal scopes and removing exactly constant tables is safe if the constants are retained. If the solver minimizes F=N-C, a returned upper or lower bound maps to N=C+F; if an additional positive integer G was divided out, use N=C+G*F. Convert to J_D only after this exact affine map, by division by H>0. A negative partial factor-sum lower bound can be valid, but calling it a lower bound for N without restoring C is incorrect. Completion should recover archived N exactly, not merely a similar displayed decimal.

## Three certificate targets

All three archives certify complete coverage of 362880=9! orders, exact incumbent=lower bound, and zero tolerances/gap. Values below are descriptive decimals only; verify the full BigInt strings in each certificate.

| Case | Frozen offset (zero based) | Warm J_D | Certified J_D | Exact minimizing orders |
| --- | ---: | ---: | ---: | ---: |
| F2_r1p55 | 30 | 0.0053979355455230065 | same | 13440 |
| H2O_scale1p10 | 35 | 0.0007001246045370596 | 0.00045499026480378383 | 44 |
| Li2_r2p95 | 12 | 0.000022640265136053056 | 0.00000037157665961366484 | 264 |

The frozen baseline window is not the certificate's selected minimizing window. For example, H2O starts from [36,37,38,39,40,41,42,43,0], while its certificate winner is [42,38,39,40,37,41,43,0,36]. Importing the winner as the warm-start order would conceal the known improvement and corrupt baseline cost accounting. F2 has unchanged optimal value but many degenerate winners; equality of objective does not imply equality of orders.

Required correctness checks:

1. Imported warm-order N and H equal the archived certificate exactly. Imported stored winner has exactly the archived incumbent N. Full order reconstructs as outside_prefix + local order + outside_suffix.
2. A new exhaustive reference covers every rank once, obtains archived optimum N and all three exact minimizer counts, and applies the same global-label lexicographic tie convention if it reports a lexicographic winner. A count of 9! callbacks alone does not independently rule out duplicate leaves.
3. DP returns a valid permutation of the same nine labels, with fixed outside order, whose independently rescored full-order N equals archived incumbent. Its restored valid lower bound equals that value on complete termination. Check the actual decomposition's coverage/running intersection and retain width, state counts and factor counts.
4. If DP does not implement global-label lexicographic selection or exact degeneracy counting, say so and compare objective plus feasible optimal witness. A different optimal order is valid in these degenerate cases. Do not equate DP table-state counts with numbers of full-order minimizers; merging independent bag suborders can allow multiple global interleavings.

The archive already has independent integer Pauli reconstruction and DFS/Lehmer coverage; reference those as retained evidence. A new Julia comparison should independently rescore its witnesses from the stored pair list, rather than merely reading back the DP scalar factor total. Full exact scalar-factor agreement at selected orders is also useful to expose sign/constant/label-transcription errors.

## Performance boundaries

Historical producer enumeration CPU times are F2 0.82174963 s, H2O 1.391734473 s, and Li2 0.829269852 s. They belong to the recorded Python implementation, host/environment, and timing boundary. They are useful historical context, not a matched timing baseline for new Julia. In particular, do not claim a language-independent DP speedup by dividing those values by new warm Julia kernel times.

Measure fresh Julia exact enumeration and Julia DP under the same environment, one-thread resource guard, arithmetic representation, domain, and completion requirements. Keep exact arithmetic on both paths. Report at least:

* input read/authentication/parse;
* objective reduction and factor-table materialization;
* graph/decomposition construction and checking;
* solver time including required witness construction;
* independent exact rescoring/checking;
* compilation/startup and whole guarded process separately.

For standalone end-to-end costs, charge each method all setup it needs. DP factor tables/decomposition cannot be charged to a shared precomputation and then omitted from a purported standalone advantage; exact enumeration does not require those objects. Materializing a nine-label factor by first evaluating all 9! orders is already exhaustive work, even if the later DP loop is quick. Prefer the at-most-four-label polynomial factor expansion above, with construction charged.

If reporting warmed solver kernels, warm both methods with matching types and declared workloads; record cold/startup totals too. Do not repeatedly run whichever method looks favorable. A paired single run supports a bounded observation; repetitions and aggregation, if chosen, must be declared and matched. Peak RSS from a shared process is a cumulative high-water mark, so label it process-wide or use separate controlled invocations before assigning memory to individual methods.

A full-bag width-eight DP has 9! states and may be slower than direct enumeration. That is an informative result. A preallocation refusal is also a result, not permission to enlarge caps silently. Graph/decomposition quality and factor count should explain the observed cost; nine-label success alone is not extrapolation evidence for 145-351 units.

No new physical-error evaluation is needed for this decision. These archived exact J_D optima were already separately propagated, and that physical study did not claim a second-order-specific objective. Keep this comparison explicitly an algorithmic exact-solver test on supplied rational first-response windows; return to cubic conditioning only as a separate subsequent study.
