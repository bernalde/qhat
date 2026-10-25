# Symmetric-formula-specific response ordering

This fixed prospective panel retains 96 full-space circuits on 3 retained molecular targets. Times are [0.83, 2.11] inverse Hartree, with both10U and30U ideal Pauli-rotation caps. All eight schedules per cell were frozen before new physical errors. No input or loss is removed.

The new objective uses the actual symmetric-step cubic defect column A3|HF>, evaluated by a degree-three Taylor recurrence. A separate literal-Pauli nested-commutator construction checks the full signed column. The diagonal-reference score is J_D2=sum(w!=HF) T^2 sinc^2(T gap_w/2) A3_w^2; optimization uses J_D2/r^4, with r the exact affordable fused symmetric step count. A positive common T^4 factor is omitted within each fixed-time cell. This is a formula-specific approximation, not an exact physical-error objective or a global optimum certificate. See DERIVATION.md and the supplied coefficients/protocol.

Three seeds are freshly optimized same-time first-order J_D, its reverse, and legacy current placement. The orientation-only control selects the lowest new score among these three; bounded insertion performs additional optimization of that new objective. Search has three equal CPU slices without unused carry, a2CPU-second nominal search allowance,1000attempted positions and at most5moves/start. All incomplete candidates/sweeps, partial improvements, rescoring and actual work remain in the frozen raw histories. Constructor setup and seed scoring are separately measured.

Strong controls include fresh same-time conserving first-order J_D (2CPU-second nominal search/1000positions/5moves per start), fresh same-time individual-Pauli J_D with a larger8CPU-second search/20m(m-1)positions/10moves per start, and separately named legacy placement whose selection reference remains T=1. Both formulas reuse all three controls. First-order orders are reused across rotation caps because their affordable first-order count is order independent. Symmetric boundary/central fusion and each actual affordable step count are retained. Formula-changing comparisons do not isolate ordering.

New-method standalone charges include upstream first-order initialization, common model/base work, current placement, cubic cache setup/cleanup, new seeds, search and complete constructor returns. These are not matched-total-CPU methods. Nominal search allowances and actual overruns are separate from complete construction charges. Shared construction and raw process resources remain distinct; standalone costs repeated across physical rows are not a sum of actual campaign work. Current placement and cache work are counted once where executed. No practical or asymptotic speedup claim follows from this small fixed-register panel.

| T | Rotation cap/U | Cubic insertion versus | Control formula | W/T/L | Floor exclusions | Median paired error ratio |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 0.83 | 10 | current | 1 | 3/0/0 | 1 | 0.000432696007407994 |
| 0.83 | 10 | first_JD_sameT | 1 | 3/0/0 | 1 | 0.000432696007407994 |
| 0.83 | 10 | pauli_JD_large_sameT | 1 | 3/0/0 | 1 | 0.0011439227796448208 |
| 0.83 | 10 | current | 2 | 2/1/0 | 1 | 0.404959554188998 |
| 0.83 | 10 | first_JD_sameT | 2 | 2/1/0 | 1 | 0.404959554188998 |
| 0.83 | 10 | pauli_JD_large_sameT | 2 | 3/0/0 | 1 | 0.19719244332073665 |
| 0.83 | 10 | cubic_orientation | 2 | 2/1/0 | 1 | 0.6564223850646479 |
| 0.83 | 30 | current | 1 | 3/0/0 | 1 | 4.623210867375916e-05 |
| 0.83 | 30 | first_JD_sameT | 1 | 3/0/0 | 1 | 4.623210867375916e-05 |
| 0.83 | 30 | pauli_JD_large_sameT | 1 | 3/0/0 | 1 | 0.00018598006492819501 |
| 0.83 | 30 | current | 2 | 1/2/0 | 1 | 0.3542330941084619 |
| 0.83 | 30 | first_JD_sameT | 2 | 1/2/0 | 1 | 0.3542330941084619 |
| 0.83 | 30 | pauli_JD_large_sameT | 2 | 2/1/0 | 1 | 0.16941465868210195 |
| 0.83 | 30 | cubic_orientation | 2 | 1/2/0 | 1 | 0.6057477422746327 |
| 2.11 | 10 | current | 1 | 3/0/0 | 0 | 0.0018131989609446067 |
| 2.11 | 10 | first_JD_sameT | 1 | 3/0/0 | 0 | 0.0018131989609446067 |
| 2.11 | 10 | pauli_JD_large_sameT | 1 | 3/0/0 | 0 | 0.0041215648843178875 |
| 2.11 | 10 | current | 2 | 2/1/0 | 0 | 0.7124183865311027 |
| 2.11 | 10 | first_JD_sameT | 2 | 2/1/0 | 0 | 0.7124183865311027 |
| 2.11 | 10 | pauli_JD_large_sameT | 2 | 3/0/0 | 0 | 0.03392473388386707 |
| 2.11 | 10 | cubic_orientation | 2 | 2/1/0 | 0 | 0.9090942872439721 |
| 2.11 | 30 | current | 1 | 3/0/0 | 1 | 0.0007678244408131816 |
| 2.11 | 30 | first_JD_sameT | 1 | 3/0/0 | 1 | 0.0007678244408131816 |
| 2.11 | 30 | pauli_JD_large_sameT | 1 | 3/0/0 | 1 | 0.0016826037301057501 |
| 2.11 | 30 | current | 2 | 2/1/0 | 1 | 0.3480905218133608 |
| 2.11 | 30 | first_JD_sameT | 2 | 2/1/0 | 1 | 0.3480905218133608 |
| 2.11 | 30 | pauli_JD_large_sameT | 2 | 3/0/0 | 1 | 0.27577857058389305 |
| 2.11 | 30 | cubic_orientation | 2 | 2/1/0 | 1 | 0.7253356477593289 |

All absolute projective state errors1-|overlap|,infidelities,norms,particle-number and joint(N_alpha,N_beta)leakage,actual rotations/steps and constructor costs appear in all_rows.csv. Full2^n propagation applies no symmetry-sector projection. The fixed error floor and absolute W/T/L tolerance are1e-12; ratios exclude either side at/below the floor. Paired differences/ratios are formed before descriptive aggregation over three related molecular cases. No independent statistical-replication claim is made.

Independent replay checked every circuit, all 336 pairwise decisions and all floor exclusions. Maximum state difference 7.6915002e-15; maximum metric difference 7.66053887e-15; full cubic-column difference 2.64762171e-15. Reference propagation uses a separate scaled-Taylor approximation; these are numerical checks, not interval or kernel proofs.

The source/protocol/input/order hashes and complete circuits/states are retained. A failed two-level correction remains rejected and was not retuned. This study supplies no conservation-mechanism identification, exact-Hamiltonian search oracle, compiled gate/depth result, hardware measurement, larger-molecule validation or whole-order optimum.
