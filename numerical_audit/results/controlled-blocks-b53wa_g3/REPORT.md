# Controlled whole-block ordering experiment

Each Hamiltonian uses one signed-priority ownership decomposition. Its nonempty completed blocks, all retained Pauli coefficients, and each block's internal Pauli order are held fixed. Only whole-block chronology changes. The two JW comparators are external controls; they change contiguity as well as order.

The primary error is stable normalized full-state E = 1 - |<exact, approximate>|, evaluated from the phase-aligned residual. No particle/spin postselection is used. Every comparison has the same Hamiltonian, HF state, total time, and Pauli-exponential count. Twenty predeclared uniform random block permutations supply a reference distribution, not a best-of-twenty algorithm.

Completed 728 evolutions across 11 Hamiltonians and 7 molecular families. Absolute tie tolerance and ratio interpretation floor are both 1e-12. Values below this floor remain in the raw data; ratios with a denominator below the floor are omitted.

| T | Comparator | Signed wins / ties / losses |
|---|---|---|
| 0.25 | block_parent_magnitude | 5 / 2 / 4 |
| 0.25 | block_l1 | 6 / 2 / 3 |
| 0.25 | block_reverse_signed | 0 / 11 / 0 |
| 0.25 | jw_signed | 5 / 0 / 6 |
| 0.25 | jw_magnitude | 6 / 2 / 3 |
| 0.25 | random_median | 9 / 0 / 2 |
| 1.0 | block_parent_magnitude | 6 / 0 / 5 |
| 1.0 | block_l1 | 7 / 0 / 4 |
| 1.0 | block_reverse_signed | 1 / 3 / 7 |
| 1.0 | jw_signed | 4 / 0 / 7 |
| 1.0 | jw_magnitude | 8 / 0 / 3 |
| 1.0 | random_median | 9 / 0 / 2 |

Primary settings (r=100): lower E and lower signed rank are better. Random intervals are the empirical interquartile range across twenty fixed seeds. Rank is the signed method's midrank among itself and those twenty samples; ties use the declared absolute tolerance.

| Case | T | Signed E | Parent magnitude E | Block l1 E | Reverse E | Random median [Q25,Q75] | Signed rank / 21 |
|---|---|---|---|---|---|---|---|
| CH2O_s1.00 | 0.25 | 7.8285e-10 | 6.8727e-10 | 1.4748e-09 | 7.8284e-10 | 1.9547e-09 [1.2621e-09, 2.7752e-09] | 4 |
| CH2O_s1.00 | 1 | 1.8718e-07 | 1.5603e-07 | 3.642e-07 | 1.8714e-07 | 4.5825e-07 [3.1215e-07, 6.7779e-07] | 4 |
| CH2O_s1.25 | 0.25 | 5.2078e-10 | 4.4957e-10 | 1.571e-09 | 5.2078e-10 | 2.3973e-09 [1.7813e-09, 3.3389e-09] | 2 |
| CH2O_s1.25 | 1 | 1.2778e-07 | 1.0145e-07 | 3.9225e-07 | 1.278e-07 | 5.8673e-07 [4.4418e-07, 8.6587e-07] | 2 |
| H2S_s1.00 | 0.25 | 1.828e-10 | 1.8731e-10 | 1.751e-10 | 1.828e-10 | 1.4802e-10 [8.1823e-11, 1.9034e-10] | 14.5 |
| H2S_s1.00 | 1 | 4.3181e-08 | 4.4234e-08 | 4.1324e-08 | 4.318e-08 | 3.562e-08 [1.9583e-08, 4.5495e-08] | 14 |
| H2S_s1.25 | 0.25 | 5.4053e-10 | 3.6994e-10 | 2.9671e-10 | 5.4053e-10 | 8.374e-10 [4.8802e-10, 1.0054e-09] | 7 |
| H2S_s1.25 | 1 | 1.2416e-07 | 8.4355e-08 | 6.6627e-08 | 1.2415e-07 | 1.9966e-07 [1.1444e-07, 2.3542e-07] | 7 |
| HCN_s1.00 | 0.25 | 1.4661e-09 | 1.379e-09 | 9.9505e-10 | 1.4661e-09 | 2.7615e-09 [1.7773e-09, 3.782e-09] | 4 |
| HCN_s1.00 | 1 | 3.3537e-07 | 3.1799e-07 | 2.319e-07 | 3.353e-07 | 6.6228e-07 [4.2404e-07, 8.7766e-07] | 4 |
| HCN_s1.25 | 0.25 | 3.6915e-10 | 4.2131e-10 | 3.7107e-09 | 3.6915e-10 | 2.5808e-09 [1.9259e-09, 3.2386e-09] | 1 |
| HCN_s1.25 | 1 | 1.0064e-07 | 1.1223e-07 | 8.847e-07 | 1.0062e-07 | 6.4035e-07 [4.7928e-07, 7.6242e-07] | 1 |
| SiH4_s1.00 | 0.25 | 1.8739e-10 | 2.0924e-10 | 2.0981e-10 | 1.8739e-10 | 2.1098e-10 [1.4683e-10, 2.484e-10] | 10 |
| SiH4_s1.00 | 1 | 4.4069e-08 | 4.9134e-08 | 4.9312e-08 | 4.4067e-08 | 5.0346e-08 [3.5184e-08, 5.9077e-08] | 10 |
| SiH4_s1.25 | 0.25 | 2.7196e-10 | 2.7405e-10 | 2.7265e-10 | 2.7195e-10 | 2.4777e-10 [2.1227e-10, 3.2292e-10] | 13 |
| SiH4_s1.25 | 1 | 6.6682e-08 | 6.7222e-08 | 6.6872e-08 | 6.6676e-08 | 6.0999e-08 [5.2302e-08, 7.939e-08] | 13 |
| f2 | 0.25 | 4.8199e-10 | 4.8168e-10 | 4.8181e-10 | 4.8199e-10 | 2.6017e-09 [1.9581e-09, 4.0343e-09] | 2 |
| f2 | 1 | 1.0582e-07 | 1.0575e-07 | 1.0578e-07 | 1.0581e-07 | 5.7104e-07 [4.2974e-07, 8.8579e-07] | 2 |
| li2 | 0.25 | 7.2261e-12 | 8.3401e-12 | 1.304e-11 | 7.226e-12 | 3.871e-10 [2.6538e-10, 6.1271e-10] | 1 |
| li2 | 1 | 2.8553e-10 | 3.0321e-10 | 1.4966e-09 | 2.8552e-10 | 6.3564e-08 [4.4405e-08, 1.0043e-07] | 1 |
| nh3 | 0.25 | 4.057e-12 | 4.1921e-12 | 4.342e-10 | 4.057e-12 | 6.2284e-10 [4.9277e-10, 9.4185e-10] | 1 |
| nh3 | 1 | 9.6315e-10 | 9.8425e-10 | 1.0831e-07 | 9.6307e-10 | 1.5495e-07 [1.2286e-07, 2.3482e-07] | 1 |

Convergence anchors fixed before the experiment: HCN_s1.00, Li2, and F2 at T=1 and r=50,100,200. All twenty random permutations are reused at every time and step count; no schedule is selected after observing an error.

| Case | Ordering | E(r=50)/E(r=100) | E(r=100)/E(r=200) |
|---|---|---|---|
| HCN_s1.00 | block_signed | 4 | 4 |
| HCN_s1.00 | block_parent_magnitude | 4.001 | 4.001 |
| HCN_s1.00 | block_l1 | 4.004 | 4.002 |
| HCN_s1.00 | block_reverse_signed | 4 | 4 |
| HCN_s1.00 | jw_signed | 4.001 | 4 |
| HCN_s1.00 | jw_magnitude | 4.002 | 4.001 |
| li2 | block_signed | 4 | 4 |
| li2 | block_parent_magnitude | 4 | 4 |
| li2 | block_l1 | 4.001 | 4 |
| li2 | block_reverse_signed | 4 | 4 |
| li2 | jw_signed | 4.016 | 4.004 |
| li2 | jw_magnitude | 3.999 | 4 |
| f2 | block_signed | 4 | 4 |
| f2 | block_parent_magnitude | 4 | 4 |
| f2 | block_l1 | 4 | 4 |
| f2 | block_reverse_signed | 4 | 4 |
| f2 | jw_signed | 4.007 | 4.002 |
| f2 | jw_magnitude | 4 | 4 |

Maximum particle-number leakage among all whole-block treatments: 8.14238e-27; maximum spin-sector leakage: 8.14238e-27. These numerical residuals support conservation in this calculation; they are not machine-checked exact certificates.

Interpretation is conditional on the frozen signed-priority ownership partition and the Hartree–Fock initial state. This experiment isolates sorting after that partition has been chosen; it does not compare different ownership assignments, establish an advantage for all states, measure compiled circuit costs, or treat repeated geometries as independent families. Source snapshots, manifest hashes, complete block membership and schedules, exact-reference state hashes, norm diagnostics and leakage are saved beside this report.

Presentation-only follow-up: after the completed campaign passed its independent integrity check, the plotting script was made robust to a below-floor first case and zero errors, and figures were sized for manuscript readability. Explicit handles preserve legends; downward markers denote E at/below the declared display floor; ratio denominators at/below that floor receive omission markers. The original plotted-source version remains in source_snapshot, and the revised version and hashes are recorded separately in presentation_source_snapshot and presentation_provenance.json. Numerical data, thresholds, summary and experiment design are unchanged.
