# Frozen molecular ordering audit

Recomputed 264 distinct settings for 11 frozen Hamiltonians in 7 molecular families from 1 completed run(s).
The counts below describe these inputs; repeated times, step counts and orderings are not independent molecular samples.

All three methods use the same retained Hamiltonian, initial Hartree–Fock state, total time, step count and number of Pauli exponentials. Errors use normalized full states, retain leakage, and ignore a global phase. Compiled gate costs are not compared.

An error at or below 1e-12 is marked below the analysis floor. A win/tie uses relative tolerance 1e-06 and absolute tolerance 1e-14. Predictor ranking checks require an actual infidelity gap greater than 1% and both state errors above the floor.

## Anchor comparisons

At r=100; each baseline is assessed separately.

| T | Comparator | Fermionic wins | Ties | Losses | Below floor |
|---:|---|---:|---:|---:|---:|
| 0.25 | JW magnitude | 7 | 0 | 4 | 0 |
| 0.25 | JW signed | 5 | 0 | 6 | 0 |
| 1 | JW magnitude | 8 | 0 | 3 | 0 |
| 1 | JW signed | 4 | 0 | 7 | 0 |

Individual anchor errors:

| Case | T | Fermionic signed E | JW signed E | JW magnitude E | Best within tolerance |
|---|---:|---:|---:|---:|---|
| CH2O_s1.00 | 0.25 | 7.82852e-10 | 6.48613e-10 | 1.93519e-09 | JW signed |
| CH2O_s1.00 | 1 | 1.87183e-07 | 1.57621e-07 | 4.8144e-07 | JW signed |
| CH2O_s1.25 | 0.25 | 5.20776e-10 | 9.72217e-10 | 1.5149e-09 | Fermionic signed |
| CH2O_s1.25 | 1 | 1.27784e-07 | 2.41721e-07 | 3.61926e-07 | Fermionic signed |
| H2S_s1.00 | 0.25 | 1.82804e-10 | 2.19062e-11 | 1.80428e-10 | JW signed |
| H2S_s1.00 | 1 | 4.31809e-08 | 4.67899e-09 | 4.26298e-08 | JW signed |
| H2S_s1.25 | 0.25 | 5.40529e-10 | 9.09067e-11 | 3.68827e-10 | JW signed |
| H2S_s1.25 | 1 | 1.24158e-07 | 1.94141e-08 | 8.41206e-08 | JW signed |
| HCN_s1.00 | 0.25 | 1.46611e-09 | 3.68575e-10 | 1.25726e-09 | JW signed |
| HCN_s1.00 | 1 | 3.35373e-07 | 8.03305e-08 | 2.94998e-07 | JW signed |
| HCN_s1.25 | 0.25 | 3.6915e-10 | 1.86097e-09 | 4.40502e-10 | Fermionic signed |
| HCN_s1.25 | 1 | 1.00638e-07 | 4.08194e-07 | 1.10576e-07 | Fermionic signed |
| SiH4_s1.00 | 0.25 | 1.87387e-10 | 1.56496e-12 | 2.09791e-10 | JW signed |
| SiH4_s1.00 | 1 | 4.40692e-08 | 3.45201e-10 | 4.92689e-08 | JW signed |
| SiH4_s1.25 | 0.25 | 2.71956e-10 | 2.19453e-11 | 2.77535e-10 | JW signed |
| SiH4_s1.25 | 1 | 6.66818e-08 | 5.21954e-09 | 6.80857e-08 | JW signed |
| f2 | 0.25 | 4.8199e-10 | 6.05059e-09 | 5.04935e-10 | Fermionic signed |
| f2 | 1 | 1.05815e-07 | 4.34574e-08 | 1.05928e-07 | JW signed |
| li2 | 0.25 | 7.2261e-12 | 1.32235e-09 | 6.61988e-12 | JW magnitude |
| li2 | 1 | 2.8553e-10 | 1.14261e-08 | 3.72798e-10 | Fermionic signed |
| nh3 | 0.25 | 4.05703e-12 | 5.83223e-12 | 4.12859e-12 | Fermionic signed |
| nh3 | 1 | 9.63148e-10 | 1.03958e-09 | 9.94819e-10 | Fermionic signed |

## Convergence and error prediction

The last two above-floor points in each case/time/ordering series give log–log slopes from -2.00154 to -1.99927 (median -2.00001); first-order state infidelity is expected to scale as r⁻² in the asymptotic regime. Exact series and successive step-error ratios are in summary.json.

Dynamic predictions use the time-propagated projected leading defect, computed before the corresponding finite-step states. The static diagnostic uses the initial-HF BCH norm. Relative prediction errors below refer to infidelity I, not E.

| r | Dynamic prediction rows | Max relative error | Static proxy max relative error | Dynamic correct rank pairs | Static correct on same pairs |
|---:|---:|---:|---:|---:|---:|
| 25 | 60 | 2.072% | 7855% | 59/59 | 56/59 |
| 50 | 60 | 0.5481% | 7979% | 59/59 | 56/59 |
| 100 | 60 | 0.1508% | 8011% | 59/59 | 56/59 |
| 200 | 59 | 0.05057% | 8020% | 57/57 | 54/57 |

Dynamic prediction is absent for 24 rows (the runner has a matrix-size limit). Omitted rows are not counted as successful predictions. The table compares both predictors on the same available rows and rank pairs; additional static-only results remain in summary.json. Pairwise ranking comparisons share Hamiltonians and are not independent statistical trials.

Anchor rankings change across sampled times for: f2 versus JW signed; li2 versus JW magnitude. The exact ratios and times are recorded in summary.json.

F₂ reverses its JW-signed/fermionic ranking between the sampled total times while step convergence remains regular. The initial-HF ratio is time independent and therefore cannot describe that reversal; the transported prediction can be checked against it in the figure. F₂ is a previously identified diagnostic example, not a newly prospective discovery.

## Numerical accuracy and old metric effects

Maximum state-norm errors: exact 5.66e-14, approximate 1.76e-12. Maximum full-space/closed-space slice residual: 0. Exact-reference phase-aligned discrepancy: 7.82e-12; the phase-sensitive discrepancy reaches 1.99 because the historical scalar-energy convention differs.

1 rows fall at or below the analysis floor. Replacing the stable metric by the old clipped unnormalized-overlap subtraction changes 0 floor labels and 1 pair classifications (including floor labels). The maximum above-floor relative metric change is 1.708. These counts use the stated tie threshold; they do not estimate floating-point uncertainty independently.

216 rows have saved prior frozen-input results. The largest above-floor relative difference is 1.793e-10. Agreement with these September frozen campaigns does not repair missing input provenance in the older August aggregate.

Observed classification changes:

- nh3, T=0.25, r=100 versus JW magnitude: stable E values 4.05703e-12 / 4.12859e-12, old values 4.85789e-12 / 4.82836e-12 (fermionic / baseline)

## Symmetry and fairness

| Ordering | Maximum particle-number leakage | Maximum spin-sector leakage |
|---|---:|---:|
| Fermionic signed | 1.6186e-34 | 1.6186e-34 |
| JW signed | 1.4008e-06 | 1.4008e-06 |
| JW magnitude | 5.7696e-09 | 5.7696e-09 |

Maximum owned-block number-commutator coefficient-l1 residual: 1.06e-16. Symmetry conservation and smaller total error are distinct properties: a leaking schedule can still have less total infidelity. Sector-conditional errors in the raw data are diagnostic and are not a free postselection benefit.

The BCH pair denominator is common to permutations of each fixed Pauli list. Its relative initial-state ratio is meaningful as a cancellation diagnostic; it is not the fixed-time fidelity coefficient. The present predictions use that time-dependent coefficient explicitly.

## Repeated-run and artifact record

Repeated-run rows: 0; disagreements beyond rtol=1e-8, atol=1e-20: 0. Repeats are not included in unique-case performance counts.

summary.json contains the thresholds, per-row diagnostics, pair comparisons, input hashes, run manifests and generator hash. normalized_results.csv is a flat export of the accepted unique rows. Source run files were read only.

Figures (convergence panels use separate vertical ranges; a horizontal dashed line marks the analysis floor when visible):

- [convergence_T0p25.pdf](convergence_T0p25.pdf)
- [convergence_T1.pdf](convergence_T1.pdf)
- [f2_time_reversal.pdf](f2_time_reversal.pdf)
