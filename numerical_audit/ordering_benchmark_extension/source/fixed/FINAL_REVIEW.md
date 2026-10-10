# Independent fixed-campaign final review

All 2,996 planned rows are present: 2,884 completed evolutions and 112 explicitly infeasible budgets. Each of 14 targets has 206 completed and eight infeasible rows. Checked unique method/case/budget keys, valid status, time and formula order, nonnegative finite errors, cost within the common budget, and saved-state presence. State hashes and independent propagation are covered by the root audit; this script does not claim to redo them.

At budget $100U$, the fresh conserving response candidate is paired with the same target, HF initial state, evolution time and rotation budget. The comparison below uses the absolute $10^{-12}$ tie floor. These 14 targets are seven molecular families at paired geometries or coefficient scales, not 14 independent molecules.

| Fixed first-order comparator | Candidate W/T/L | Median error ratio |
|---|---:|---:|
| bch1_full_l1_insertion | 14/0/0 | 0.0048824245 |
| class_bch_l1 | 14/0/0 | 0.050658893 |
| class_color | 14/0/0 | 0.05491451 |
| class_lex | 14/0/0 | 0.050605102 |
| class_magnitude | 14/0/0 | 0.051055828 |
| class_suffix_l1 | 14/0/0 | 0.051049106 |
| hastings_interleaved_lex | 14/0/0 | 0.0045358878 |
| jw_color | 14/0/0 | 0.0062535281 |
| jw_lex | 14/0/0 | 0.0046399789 |
| jw_magnitude | 14/0/0 | 0.024459441 |
| jw_signed | 14/0/0 | 0.046485995 |
| owned_color | 14/0/0 | 0.0035538191 |
| owned_l1 | 14/0/0 | 0.0070009268 |
| owned_parent_magnitude | 14/0/0 | 0.041861213 |
| owned_reverse | 14/0/0 | 0.040248093 |
| owned_signed | 14/0/0 | 0.040243778 |
| parent_color | 14/0/0 | 6.359471e-05 |
| parent_magnitude | 14/0/0 | 0.0010344339 |
| parent_signed | 14/0/0 | 0.0014269856 |
| refined_magnitude | 14/0/0 | 0.051811955 |
| tranter_deplete_groups | 14/0/0 | 0.0029059627 |
| tranter_equalise_groups | 14/0/0 | 0.0019764642 |

The candidate wins against all 22 deterministic fixed first-order controls in all 14 targets. This is a statement about this specified panel, not a general ordering theorem. Classical preprocessing cost and second-order formulas require the separate optimized/published comparisons.

| Fixed-random family | Candidate versus sample median W/T/L | Median ratio | Strictly better than individual samples |
|---|---:|---:|---:|
| owned | 14/0/0 | 0.0034914786 | 280/280 |
| refined | 14/0/0 | 0.014395495 | 280/280 |
| pauli | 14/0/0 | 0.0027347345 | 280/280 |

Random permutations use 20 predeclared seeds per family and are reused across every time step. No measured-error selection is made. Stronger controls matter: the full-budget unconstrained Pauli response optimizer wins 10 of 14 against the conserving first-order candidate, and symmetric formulas often do better at the same rotation budget.

## Historical baseline correction

Of 96 historical schedule checks, 64 reproduce literal historical flat-sequence hashes and 32 reproduce the same ordered commuting blocks with different internal word order. Those 32 ideal products agree exactly for every time and step count. Saved historical versus fresh states differ by at most 1.1998993524231055e-14 in vector norm and 6.0507798587110945e-18 in the reported error. See HISTORICAL_BASELINE_ALIGNMENT.md and historical_internal_order_check.json. No fixed measured schedule was changed after this distinction was identified.

All frozen sources, 14 design hashes and eight QHAT dependency hashes remain unchanged at campaign completion (final_execution_check.json). Row timing measures evolution only; use timing_scope.json before comparing times across campaigns.
