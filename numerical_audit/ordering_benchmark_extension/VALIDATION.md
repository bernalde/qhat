# Validation coverage

Passed coverage gates for all 14 cases: 4,312 planned rows, 4,186 saved states and 126 infeasible budgets.

| Campaign | Audited rows | Saved states checked | Independently replayed | Maximum replay distance |
|---|---:|---:|---:|---:|
| fixed | 2,996 | 2,884 | 2,216 | 0 |
| optimized | 728 | 728 | 588 | 1.15867e-13 |
| extra | 392 | 392 | 312 | 0 |
| v2 | 196 | 182 | 142 | 1.15824e-14 |

All saved states received independent metric checks; 3,258 circuits were independently replayed according to the predeclared, error-independent selection. The maximum replayed state-vector distance is 1.15867e-13; the maximum saved-metric difference is 3.33067e-16. Replay coverage is a subset, not a claim to rerun every circuit. Counts are derived from the per-row check records.

Independent target reference checks cover all 14 cases (maximum distance 7.88651e-14). The 56 objective checks and four exact partition-family checks per case pass; the latter verify target sums, internal commutation and particle-number commutators. Maximum objective-score difference: 0.

All 30 checked archived optimized states agree bitwise. Historical fixed schedules comprise 64 literal flat-sequence matches and 32 exact commuting-block equivalences; the latter have maximum saved-state distance 1.1999e-14 and error difference 6.05078e-18. Internal permutations are not described as literal sequence matches.

Uniform boundary sensitivity covers 70 symmetric plans and 420 existing states, with 249 dense toy-product checks (maximum difference 5.70489e-15). It changes interpreted cost only, not a state, step count, or error, and makes no newly saturated-budget measurement claim.

The manuscript correspondence check covers all 19 public bridge statements. It is a source/proof-step review, with no new Lean compilation or axiom audit. Prior-artifact verification found 766 QHAT files, 94 embedded files and 49 earlier paper files unchanged; this is byte-integrity evidence, not a new scientific rerun.

## Reference-checkpoint limitation

This is an aggregation of the saved independent reference checks, not a new reference solve. Reference-checkpoint reuse verifies the checker and core source hashes but skips rereading the original input bytes and regenerating or rehashing the reference NPZ. The prior-artifact manifest check separately establishes original-artifact byte integrity. Campaign audit checkpoint reuse rechecks audited raw rows and state hashes, but does not independently revalidate a reused reference NPZ. These limitations are not hidden by the aggregate pass.

This aggregator rehashes the analyzed raw JSON rows and cross-links their existing state digests and per-row audits. It does not reread all state arrays, redo matrix propagation, regenerate references, or rebuild Lean. Input evidence hashes and the exact coverage are retained in validation.json.
