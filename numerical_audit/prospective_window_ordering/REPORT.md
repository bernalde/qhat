# Prospective molecular window-ordering study

Sparse-parent window insertion improved physical error over full-order cubic training in 14 of 32 paired observations, tied in 18 and lost in none, at substantially greater construction cost. The magnitude-parent variant was weaker than cubic training (4 wins, 2 ties, 26 losses), although it improved its own magnitude baseline in 12 observations (18 ties, 2 losses). These are descriptive results on four fixed molecular instances, not a universal ordering advantage.

The [independent derived audit](evidence/run_2026-10-28/INDEPENDENT_DERIVED_CHECK.main.H.a1.json) verified all 352 physical rows, 384 paired observations, 192 paired-repetition medians, 12 comparison summaries and 32 Pareto condition groups. All 64 constructor charges were reconciled. The retained result has no failures, censored cases, unverified claimed objectives or accounting veto. Artifact review, publication and CI are separate gates from this numerical audit.

## Methods and scope

The [fixed registry](REGISTRY.json) contains four 14-mode, four-active-electron instances in the previously studied LiH and BeH₂ families, using 6-31G and seven spatial active orbitals. These are new geometry/basis/active-space instances, without a new-family, more-atoms or universally-more-terms claim.

| Case; distance in Å | Original Pauli words U | Conserving units |
| --- | ---: | ---: |
| LiH, 1.487 | 1,289 | 286 |
| LiH, 2.113 | 1,289 | 286 |
| BeH₂, 1.369 | 665 | 198 |
| BeH₂, 1.907 | 757 | 218 |

The [protocol](PROTOCOL.md) compares canonical, magnitude, current placement, conserving first-response insertion, individual-Pauli insertion, full-order sparse cubic training, and window insertion with fresh sparse or magnitude parents. Two fresh deterministic process repetitions, with reversed arm launch order, produce 64 whole construction calls. They are not independent chemistry samples and support no confidence interval.

Training uses HF, T=1 and 30U ideal rotations. Cubic training minimizes transported cubic-response norm squared divided by the fourth power of affordable steps. Window insertion uses one strictly interior width-eight window, at most ten sweeps, exact-real local differences before transport, full cubic Taylor checks and authoritative native final acceptance. Fresh parent construction is charged. First/Pauli retain their first-response objectives; the native 3,000-evaluation setting does not equalize arithmetic across arms, and their seed/incremental counters exclude separately recorded full move rescoring.

All 64 orders/checkpoints were frozen before the four PANEL authentication jobs and physical evaluation. Canonical, current and conserving-first use formulas 1 and 2; other arms use formula 2. Eleven paths, four cases, two repetitions, T∈{0.4,1.6} and budgets {8U,24U} give 352 rows. Each order recomputes affordable steps and central/boundary fusion from original words. Observed validation-step sets were {8,24} for formula 1 and {4,12} for formula 2.

Circuits propagate the full 16,384-dimensional register. The numerical exact-H reference uses the 1,001-dimensional N=4 sector, without constraining circuit leakage. Independent checks reconstruct references at both times and recompute projective error, infidelity, norm drift and number/spin leakage for every row. Projective error is 1−|normalized overlap|. Fresh unfused propagation additionally checks distinct hard-setting circuits at T=1.6,24U against both saved and independently reconstructed references. These are floating numerical checks, not interval-certified physical-error bounds.

## Paired physical outcomes

Each row below covers 32 formula-two observations: 16 case/time/budget settings × two repetitions. For Δ=window error−baseline error, wins require Δ<−10⁻¹², losses Δ>10⁻¹²; the remainder are ties. The last column is the median of the 32 already-paired differences, not a difference of method medians. The separate 192 within-setting two-repetition medians remain in the [complete analysis](evidence/run_2026-10-28/ANALYSIS_SUMMARY.main.H.a1.json).

| Window parent | Baseline | Wins | Ties | Losses | Median paired Δ |
| --- | --- | ---: | ---: | ---: | ---: |
| Sparse | Cubic | 14 | 18 | 0 | −4.90123e−13 |
| Sparse | Current | 22 | 8 | 2 | −1.04542e−10 |
| Sparse | Magnitude | 26 | 2 | 4 | −1.49091e−9 |
| Sparse | Canonical | 26 | 2 | 4 | −1.47439e−9 |
| Sparse | First | 22 | 8 | 2 | −1.04542e−10 |
| Sparse | Pauli | 22 | 8 | 2 | −9.33686e−11 |
| Magnitude | Cubic | 4 | 2 | 26 | +1.48768e−9 |
| Magnitude | Current | 6 | 2 | 24 | +1.45358e−9 |
| Magnitude | Magnitude | 12 | 18 | 2 | −1.75223e−14 |
| Magnitude | Canonical | 4 | 10 | 18 | +1.50986e−11 |
| Magnitude | First | 6 | 2 | 24 | +1.45565e−9 |
| Magnitude | Pauli | 6 | 2 | 24 | +1.45393e−9 |

There are no censored pairs. Sparse-versus-cubic’s median lies inside the ambiguity band despite its 14 resolved wins. Ten physical rows are at/below the analysis floor; floor operands suppress error ratios. There are no objective-floor refusals. No global optimum, general approximation ratio or universal winner follows.

## Whole-call construction costs

Each arm total covers eight calls, not repeated charges across physical rows. Raw guard CPU and absolute wall include imports, parent search, Julia startup/JIT, validation and cleanup. Each call had 64 CPU/320 wall seconds under serialized quarter-core, single-thread execution and a 512 MiB hard memory limit. The child's 32M heap hint is separate from its 64 MiB exact-local cache allowance.

| Arm | Raw CPU, s | Absolute wall, s | Conservative FINAL CPU, s |
| --- | ---: | ---: | ---: |
| Canonical | 0.920 | 24.454 | 1.011 |
| Magnitude | 0.870 | 24.468 | 0.969 |
| Current | 2.110 | 24.422 | 2.216 |
| First | 2.950 | 24.455 | 3.044 |
| Pauli | 13.370 | 72.803 | 13.464 |
| Cubic | 28.467 | 127.120 | 28.561 |
| Window, sparse parent | 177.660 | 722.041 | 177.767 |
| Window, magnitude parent | 151.330 | 613.384 | 151.439 |

The 64-call totals are 377.677049 raw CPU seconds, 378.471 conservative CPU seconds and 1,633.146128 wall seconds. Cost-table displays round to three decimals; full precision and table provenance remain in the [numeric tables](evidence/run_2026-10-28/REPORT_NUMERIC_TABLES.a1.json) and their [extractor](evidence/run_2026-10-28/tools/report_tables.py). Sparse-parent window construction costs approximately 6.24× cubic CPU and 5.68× cubic wall: ratios of eight-call sums, not median speed ratios. Conservative charges remain distinct from raw measurements. These timings are neither compiled-gate counts nor hardware-performance guarantees.

## Fixed-condition figures and full data

The predeclared hard-condition figures show all four cases, eleven paths and both repetitions: 88 expected points at T=1.6,24U. [CPU PNG](evidence/run_2026-10-28/plots/main-hard.a2/cpu_vs_error.png) / [SVG](evidence/run_2026-10-28/plots/main-hard.a2/cpu_vs_error.svg) and [wall PNG](evidence/run_2026-10-28/plots/main-hard.a2/wall_vs_error.png) / [SVG](evidence/run_2026-10-28/plots/main-hard.a2/wall_vs_error.svg) use a logarithmic error axis with an explicitly marked 10⁻¹² display floor. Actual values remain in the [CSV](evidence/run_2026-10-28/plots/main-hard.a2/plotted_data.csv) and [JSON](evidence/run_2026-10-28/plots/main-hard.a2/plotted_data.json). Revision a2 moves annotations only; its 88-row CSV is byte-identical to a1. Repeated points may overlap; no jitter or confidence intervals are added. The tables cover both times and budgets. Pareto memberships are descriptive within matching case/repetition/condition and cost axis.

The [Pareto interpretation notes](PARETO_INTERPRETATION.md) distinguish training-objective status from physical validation, explain the absolute error band and short-run timing noise, and link a separately labelled CPU-cost sensitivity view. That view removes one-ULP cost tie-breaking using a fixed 10⁻¹²-second CPU-cost tolerance across all 32 groups, with independent checking; the original primary frontiers, measured costs, errors and resource ceilings remain unchanged. These axes measure ordering construction, excluding later propagation and quantum execution.

## Provenance and retained accounting

Validation cohorts retain their original revision and scope. The [F installed-wheel API run](evidence/run_2026-10-28/validation/integration_2026-10-28/PRODUCT_API_CHECK.a5.json) passed 163 tests. The [G five-test runtime check](evidence/run_2026-10-28/validation/integration_2026-10-28/PRODUCT_API_CHECK.a6.json) completed its functional tests but [exceeded its 35-CPU-second guard](evidence/run_2026-10-28/validation/integration_2026-10-28/HEAP32_RUNTIME_CHECK_QUALIFICATION.a1.json); the subsequent [G installed CLI check](evidence/run_2026-10-28/validation/integration_2026-10-28/PRODUCT_API_CHECK.a7.json) passed one test. These are narrower checks, despite inherited broader historical `scope` strings. The [final consumer metadata cohort](evidence/run_2026-10-28/validation/integration_2026-10-28/SOURCE_REVIEW_METADATA_CHECK.a4.green.json) passed 120 test methods plus 170 subtests on the exact later H consumer bytes; earlier F binding records remain F records. Frozen `DRAFT/UNRUN` headers describe the scripts before their retained executions.

The [43-file scientific snapshot](evidence/run_2026-10-28/data-a6/INPUT_PROTOCOL_FROZEN.json) is bound to G `722013aa6d9a0a48776e26de5f356b5a441c3f61`. H `c27fcd7337a965b879e41dc094ba8fb0f8b5bd8e` preserves the scientific and 19-file product identities and supplies the separately bound analysis/checker. [Input reuse](evidence/run_2026-10-28/data-a6/INPUT_REUSE_PROVENANCE.json) retains original chemistry producer `87da247744ed58b7a96759fd8af3a44235b6cb38`; it does not describe regenerated SCF data. The [analysis index](evidence/run_2026-10-28/ANALYSIS_INDEX.main.H.a1.json) and [check index](evidence/run_2026-10-28/INDEPENDENT_DERIVED_INDEX.main.H.a1.json) bind the actual sources and outputs.

The [validation cohort notes](VALIDATION_COHORTS.md) distinguish call reports, test methods and subtests, retain the expected red-run failures, and link the exact metadata drivers and their unchanged test-selection logic.

Historical validation/build wrappers have mixed provenance: 22 helper versions have retained hash evidence; 16 have [retrospective exact-file readbacks and filesystem continuity observations](evidence/run_2026-10-28/HISTORICAL_HELPER_READBACK.a1.json), rather than original runtime-byte attestations. This limitation concerns those historical wrappers; the scientific/product snapshots and complete main numerical audit remain separately hash-bound.

Three source versions from the earliest integration freezes could not be recovered: the initial CI workflow and the initial panel driver and its test. Those early attempts therefore cannot be fully reproduced from their exact declared source. Their declarations and the failed closure audit remain retained, and historical source closure remains explicitly incomplete. The [archive scope](evidence/run_2026-10-28/ARCHIVE_SCOPE.a1.md) lists the exact missing hashes and distinguishes this limitation from the complete final scientific snapshots and main numerical checks.

The demonstrated portability covers the analyzer and independent derived checker. Historical dispatchers retain machine-specific guard paths, and the fourth integration resource extension lacks an explicit declaration timestamp. The [future campaign-template requirements](REPRODUCIBILITY_FOLLOWUPS.md) record the requested improvements while preserving the frozen historical records.

The complete six-phase scientific/output ledger through **2026-10-29 03:15:13.811928 UTC** contains 184 jobs, including earlier failed attempts. These six phases had complete raw accounting and stayed within their cumulative ceilings; no main-call conservative-cap veto was triggered. The separate [integration/build/test ledger](evidence/run_2026-10-28/validation/integration_2026-10-28/) and its [versioned allowances](evidence/run_2026-10-28/validation/integration_2026-10-28/RESOURCE_EXTENSION.a4.json) are excluded from these counts and sums; the main independent audit does not claim to reconcile that separate ledger:

| Phase | Jobs | Raw CPU, s | Conservative CPU, s | Wall, s |
| --- | ---: | ---: | ---: | ---: |
| Identity/input | 16 | 58.178 | 58.378 | 262.377 |
| Preflight | 80 | 513.936 | 515.062 | 2,210.886 |
| Constructors, including setup | 69 | 384.677 | 385.532 | 1,672.600 |
| Physical | 5 | 131.090 | 131.159 | 535.333 |
| Independent | 8 | 144.850 | 144.970 | 592.962 |
| Packaging | 6 | 1.800 | 1.881 | 18.208 |

The first LiH physical launch failed on a missing parent output directory before producing outputs: 7.650 raw CPU seconds, **7.656 conservative CPU seconds**, and 33.288 wall seconds. Its [a1 receipt](evidence/run_2026-10-28/logs/window-main-physical-LiH-r1p487-631g-14m4e-a1.resources.json) and [bounded a2 retry](evidence/run_2026-10-28/logs/window-main-physical-LiH-r1p487-631g-14m4e-a2.resources.json) remain retained; no algorithm changed. Earlier preflight refusals/gate failures remain alongside [G/a5 admission](evidence/run_2026-10-28/PREFLIGHT_DECISION.a5.json).

Postcutoff indexing, checks, plotting, packaging, guard proofs and review costs require separate later accounting. The table is a cutoff ledger, not all future campaign work. Targeted [actual-data guard proofs](evidence/run_2026-10-28/main-negative-guards.a1/OUTCOME.json) retained a passing mirror baseline, precise conservative-charge and paired-delta rejections, and unchanged originals. No training source, input, configuration or budget was retuned after inspecting completed main outcomes. Publication and live CI/merge status are tracked separately from these frozen numerical results.

The [artifact closeout](evidence/closeout_2026-10-29/CLOSEOUT_MANIFEST.a1.json) records 15 later jobs, bringing the six-phase scientific/output ledger to 199 closed jobs, excluding the separate integration/build/test ledger. Analyzer and independent-checker replays from the copied archive exactly reproduce their original JSON outputs; staged Git bytes and report links were checked separately. One archive-audit launch failed before worker execution because its sandbox could not access the systemd bus: raw CPU was unavailable, so its full 20-CPU-second reservation remains conservatively charged. The bounded continuation used the established resource-guard launch context. These closeout checks preserve the historical provenance limitations above.
