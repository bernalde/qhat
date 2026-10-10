# Final published-comparator review

The recovered published-comparator extension is complete across all 14 frozen molecular cases. The four extra first-order controls produced 392 states. The V2 campaign produced 196 rows: 182 states and 14 infeasible native-Strang budget-one rows. All 28 designs, both producer freezes, all 588 collected/raw rows, and all 574 state hashes pass the final integrity check in `final_evidence_review.json`. All 11 fresh native histories reconstruct their saved prefix orders exactly. The other three V2 designs use verified positive coefficient scaling.

The shared core API also supplies the fixed campaign's Tranter deplete/equalise schedules, literal Hastings schedule, and explicit first-order BCH insertion. Its source remained identical to the fixed campaign's freeze.

## Source fidelity and mathematical checks

- The Tranter V2 control uses the coefficient-l1 norm of the incremental second-order error operator, with exact arbitrary-precision integers and leftmost prefix-position ties. Its physical schedule reverses the paper's prefix indexing. The same order is evaluated both first-order and with native Strang.
- The first-order BCH insertion control minimizes the full first-order error polynomial and is labeled as an adaptation, not as Tranter's V2 algorithm.
- Hastings's unspecified final family order is completed lexicographically. Both literal and within-family consolidated implementations are retained; dense checks verify identical one-step unitaries.
- The Tomesh control explicitly discloses the deterministic beam32 adaptation, integer clique labels, positive zero-connection denominator, retained finalists, and Christofides tours. The unlicensed third-party source inspected for comparison is not redistributed.
- Independent checks cover 256 dense Pauli products; 14 graph-grouping cases; 20 dense greedy V2 histories; 20 dense first-order BCH histories; 20 completed-checkpoint restarts; 40 dense Appendix-A heuristic cases; 30 RLF graph cases; and 12 exhaustive small-clique beam comparisons. Three matrix-log checks verify the V2/physical-Strang sign and chronological convention.
- A real 1184-term F2 search resumed from its first 592 accepted stages and reproduced every remaining exact norm and the final order. All six family/campaign uniform-scaling checks pass, including cubic V2 and quadratic Tomesh diagnostics.

## Completed timing record

The uninterrupted 5372-term NH3 native V2 search took **4952.206595 s (82.537 min)**. Including case loading and preparation before final design JSON output took 4952.863756 s. Its coefficient-1.25 companion reused the exact order: native search time is zero and cache verification/preparation took 0.597614 s before final JSON output. This reuse is not a fresh search or an algorithmic speedup. The first-order and Strang arms share the same preparation. See `timing_scope.json` for all inclusion/exclusion boundaries and `hardware_context.json` for the descriptive hardware record.

## Final uniform-compilation sensitivity

`../optimized/boundary_sensitivity.json` now covers **70 plans and 420 symmetric rows**, with no pending V2 plans. Its 249 dense exact-equivalence tests have maximum Frobenius discrepancy 5.7e-15. The analysis changes only rotation accounting at fixed measured states and step counts; it does not fill newly available budgets.

The standalone summary compares with the previously frozen conserving first-order candidate, so the campaign-wide regenerated analysis is the authoritative manuscript table. That summary records 14/14 wins over each additional first-order comparator, including the V2 order applied first-order, but only 2/14 wins against native V2 Strang. Higher-order comparators therefore prevent a general superiority claim. Under uniform boundary compilation, the conserving Strang method uses less tested work than V2 Strang in 12/14 cases at error 1e-6 and 10/14 at 1e-8; the median V2/conserving work ratios are about 1.967 and 2.002. These statements concern the finite tested grid, not a universal theorem or hardware gate-count advantage.

No new Lean verification is claimed for these numerical implementations or for the hand-derived incremental V2 update. The manuscript methods appendix was checked against the final frozen source and the primary-source locators.
