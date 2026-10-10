# Uniform boundary-fusion sensitivity

This postprocessing retains every measured state, step count and error. It changes only the rotation accounting using a single compilation rule for all symmetric schedules. It does **not** create observations at the larger step counts that this compilation might permit under 100U.

For U ordered words, let a be the maximal pairwise-commuting prefix length and b the corresponding suffix length. Overlap c=max(0,a+b-U) contains words commuting with the entire list: each overlap word is in both sets, whose union covers all words. Extract those c words into one full-time rotation each. The remaining prefix and suffix are disjoint and commute internally. Fuse the suffix at each palindrome center and the prefix across successive steps. The cost is r(2U-a-b)+a, including the overlap case. If all words commute, a=b=U and the entire product costs U rotations, independent of r.

For a complete diagonal class followed by complete nonzero flip classes, refining the diagonal class into consecutive singleton words leaves the flattened list unchanged. Since those diagonal terms commute, their forward and reverse half-exponentials combine into the corresponding half-exponential of D. The complete-class label Strang formula and refined-label Strang formula therefore define exactly the same unitary at every same r. All 14 actual flattened lists agree exactly.

Dense tests include X0,Z1,Z0, where AB and BC commute but AC do not, plus an all-commuting list and 80 random lists. All integer pulse durations are checked before comparing full unitaries.

Checked 249 dense products; maximum Frobenius discrepancy 5.7e-15.

## Feasible step-count sensitivity at 100U

| Schedule | Plans | Cases permitting more steps | Largest old→possible r |
|---|---:|---:|---|
| strang_class_lex | 14 | 0 | CH2O_s1.00: 56→56 |
| strang_class_refined_response | 14 | 4 | H2S_s1.25: 50→54 |
| strang_jw_magnitude | 14 | 12 | H2S_s1.00: 50→56 |
| strang_jw_signed | 14 | 10 | H2S_s1.00: 50→56 |
| tranter_v2_l1_strang | 14 | 11 | CH2O_s1.00: 50→53 |

These possible extra steps were not simulated by this sensitivity. The following table uses only already-tested r values with recompiled costs. A ratio above 1 means the comparator needs more tested rotations than conserving-order Strang.

| Error target | Comparator | Paired cases | Conserving W/T/L | Median comparator/conserving work |
|---|---|---:|---|---:|
| 0.0001 | strang_class_lex | 14 | 5/6/3 | 1 |
| 0.0001 | strang_jw_magnitude | 14 | 3/2/9 | 0.995619 |
| 0.0001 | strang_jw_signed | 14 | 4/0/10 | 0.989389 |
| 0.0001 | tranter_v2_l1_strang | 14 | 4/0/10 | 0.993082 |
| 1e-06 | strang_class_lex | 14 | 9/0/5 | 1.43968 |
| 1e-06 | strang_jw_magnitude | 14 | 11/0/3 | 1.93026 |
| 1e-06 | strang_jw_signed | 14 | 11/0/3 | 1.99292 |
| 1e-06 | tranter_v2_l1_strang | 14 | 12/0/2 | 1.96657 |
| 1e-08 | strang_class_lex | 14 | 6/0/8 | 0.996842 |
| 1e-08 | strang_jw_magnitude | 14 | 7/0/7 | 1.4676 |
| 1e-08 | strang_jw_signed | 14 | 10/0/4 | 2.01658 |
| 1e-08 | tranter_v2_l1_strang | 14 | 10/0/4 | 2.00226 |

Pending V2 plans: none. V2 state coverage may still be incomplete until that campaign finishes; every table reports its paired count. Original emitted costs remain unchanged in the raw records. This is not a full commuting-gate optimizer or a hardware gate-count estimate.
