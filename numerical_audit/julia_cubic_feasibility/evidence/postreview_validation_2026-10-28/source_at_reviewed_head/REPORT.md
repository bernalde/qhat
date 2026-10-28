# Julia ordering DP and cubic-window decision study

The Julia dynamic program reproduced all three archived nine-unit optima exactly and reduced warm computation time relative to a fresh run of the unchanged Python enumerator. The gain varied substantially: about 38.5 times for F2, 1.15 times for H2O, and 3.83 times for Li2. These are restricted first-order response certificates, with every outside position fixed. They establish a useful exact solver for these domains, not a full molecular ordering algorithm.

The separate cubic Strang pilot gave a more cautious practical result. Cheap insertion within an eight-unit window found the same frozen-table optimum as complete enumeration on both development cases. CH4 improved its surrogate only slightly. CO improved its surrogate by about 24%, but its physical error became worse at the longer validation time. The existing whole-order optimizer, continued for an additional 32 CPU seconds, produced better physical errors than either parent or window order on every tested setting. Full cubic molecular DP was not established. The immediate candidate for further practical study is inexpensive conditioned-window insertion, with physical validation and its costs retained; these results do not support a universal replacement of the existing optimizer.

## Fixed questions, inputs, and arithmetic

The [protocol](PROTOCOL.md) separates two objectives. The archived comparison uses F2_r1p55, H2O_scale1p10, and Li2_r2p95 from the completed `exact_response_windows` study. Only their frozen contiguous nine-unit window varies: 9!=362,880 permutations per case. Units, integer pair amplitudes, supplied binary64 sinc-squared weights, outside order, and normalization remain unchanged. The exact objective is first-order `J_D=N/H`, where N is an integer weighted coherent square and H is the retained positive denominator. No cubic coefficient or Strang step-count factor is inserted.

The adapter expands this square into signed factors on at most four labels, retaining the fixed outside amplitude inside the square and preserving the objective constant. Bag-order DP uses exact integer/rational arithmetic, checks its decomposition, and reconstructs a feasible whole-order witness. Direct pair-amplitude rescoring checks that witness independently of the factor expansion. Fresh Julia DFS enumeration verifies complete Lehmer-rank coverage, optimum, degeneracy, and the original global-label tie rule. DP is not claimed to count all minimizers or return the same tied witness.

The cubic pilot instead uses the retained twelve-mode HF inputs CH4_r1p09 (351 units, 1,722 words) and CO_r1p13 (195 units, 770 words), with frozen cubic parent orders. Training time is T=1 and the ideal rotation cap is 30U, with U the original word count. Its objective is the squared transported, projected cubic response divided by the affordable fused Strang step count to the fourth power. Both selected windows leave global endpoints fixed; r=15 throughout each restricted search.

All new numerical kernels are Julia with standard-library dependencies. The existing Python optimizer and archived Python enumerator remain unchanged controls. Inputs are direct TOML serializations of retained data, with source/input bindings recorded under [evidence](evidence/). These are known development cases, not new holdouts.

## Archived DP: exact agreement and measured cost

All three DP minima equal the archived integer optimum exactly. Fresh exhaustive checks recovered 13,440 minimizing orders for F2, 44 for H2O, and 264 for Li2. Their denominator mappings, fixed outside positions, and full-order witness scores also match. Degeneracy explains why a valid DP witness can differ from the archive's canonical winner; the changed witnesses received physical checks rather than inheriting another order's errors.

The table reports median process CPU seconds over three warm repetitions. Julia DP charges include factor construction, graph/decomposition, solving, and witness validation. Parsing and process startup are separate. Python timings are fresh executions of the unchanged retained solver, not the historical archive times. The [comparison table](evidence/ARCHIVE_COMPARISON_TABLE.json) and raw samples retain the exact values.

| Frozen case | Julia DP CPU s | Fresh Python enumeration CPU s | Fresh Julia DFS CPU s | DP bag states | Used width |
| --- | ---: | ---: | ---: | ---: | ---: |
| F2_r1p55 | 0.019660 | 0.756453 | 1.954923 | 1,737 | 5 |
| H2O_scale1p10 | 1.164137 | 1.337233 | 8.072798 | 46,353 | 7 |
| Li2_r2p95 | 0.199409 | 0.764457 | 1.926618 | 6,753 | 6 |

Widths describe the checked decompositions actually used; they are not certified minimum treewidths. State counts are bag permutations, not numbers of feasible global orders. The Julia DFS was slower than the Python control here, so language choice alone does not explain the DP result. The comparison supports these implementations and these three domains, without a scaling claim.

First invocations share one Julia process and compiled helpers. They are retained as first-call observations, not independent cold-start comparisons. Whole guarded process receipts include startup/JIT costs; warm repetitions alternate method order. Setup and verification remain visible instead of being omitted as shared preprocessing. The weight-times-seven H2O preflight preserved the minimizers and checked the full adapter before main measurements.

## Cubic structure and restricted optima

The full-order structural scans each stopped at 150,000 examined pair/triple scopes. They retained 26,192 computed nonzero atoms for CH4 and 21,207 for CO. Destination-overlap diagnostics produced heuristic width bounds of 346 and 192 for the retained potential graphs. These are partial floating support scans and potential interaction graphs, not the exact coherent factor graphs. The large widths are not lower bounds on the best possible decomposition, and do not prove intrinsic molecular hardness. They supplied no affordable route to a full bag-order DP within this pilot's admission limits; none was attempted.

The bounded alternative selected an interior eight-unit window using computed nonzero internal cubic pair scopes, with earliest-start tie breaking and no physical-error input. Starts are Julia positions 157 for CH4 and 168 for CO. The conditioned representation retains the entire baseline response and adds pair/triple changes; outside interference is preserved. Complete enumeration covered all 40,320 window permutations. An explicit roundoff envelope screened candidates, followed by exact dyadic evaluation of every overlapping contender.

| Case | Parent J | Restricted optimum J | Same-window insertion |
| --- | ---: | ---: | --- |
| CH4 | 2.666715076e-9 | 2.665840143e-9 | Same order and exact table optimum |
| CO | 1.913095251e-8 | 1.454631048e-8 | Same order and exact table optimum |

Each frozen-table optimum was unique, with one exact contender after screening. CH4's parent/optimum ratio was 1.000328; CO's was 1.315176. These certificates concern supplied frozen binary64 response DELTA tables. They do not certify exact transcendental kernel weights, the full permutation domain, or physical error. Independent literal nested-commutator checks support agreement with the intended cubic objective at the returned orders.

Construction costs matter when the improvement is this small:

| Case | Selection CPU s | Response-table setup CPU s | Insertion CPU s | Additional exact enumeration CPU s |
| --- | ---: | ---: | ---: | ---: |
| CH4 | 0.234623 | 0.245528 | 0.006655 | 0.356271 |
| CO | 0.012066 | 0.048967 | 0.002292 | 0.280780 |

Parsing and direct Taylor validation are additional; shared JIT/startup cannot be assigned independently to each method. The complete window batch used about 8.2 service CPU seconds including JIT. These are additional-ordering costs from an existing parent. Historical parent construction cost was about 29.985 CPU seconds for CH4 and 6.405 for CO; a standalone end-to-end proposal cannot treat those parents as free.

## Physical results and continuation control

Orders were frozen before physical evaluation. All cubic circuits use the symmetric formula at T=0.5 and 1.5 with 10U and 30U caps. Original Pauli-word rotations are propagated in the full 4,096-dimensional space with central/boundary fusion; circuits are never projected into a number sector. Only the exact-H numerical reference uses a checked fixed-number sector. Recorded coupling/Hermiticity checks are numerical diagnostics, not interval or operator-norm certificates. Error is stable projective distance `1-|overlap|`; norm drift and number/spin leakage are retained.

Window insertion and exact enumeration return identical orders, so their physical results coincide. Approximate window/parent error ratios expose the limitation of surrogate optimization:

| Case | T | Ratio at 10U | Ratio at 30U |
| --- | ---: | ---: | ---: |
| CH4 | 0.5 | 0.99972 | 0.99972 |
| CH4 | 1.5 | 0.999748 | 0.999748 |
| CO | 0.5 | 0.783 | 0.784 |
| CO | 1.5 | 1.0209 | 1.0244 |

Using the fixed absolute 1e-12 ambiguity, both CH4 T=0.5 comparisons are ties despite their ratios below one. CH4 improves at T=1.5. CO improves at the shorter time and loses at both longer-time budgets. For example, at T=1.5 and 10U, CO error rises from 3.401064e-6 to 3.472126e-6. Its 24% lower training surrogate therefore does not justify an unconditional physical upgrade.

The prospectively recorded [control amendment](evidence/CONTROL_AMENDMENT.json) admitted an unchanged legacy whole-order continuation from the same parent, with 32 additional CPU seconds. It returned J=2.359617128e-9 for CH4 and 1.010389778e-8 for CO, with constructor charges 30.712 and 31.729 seconds. Both terminated with `cpu_limit`; neither is a local/global optimality claim. Their physical errors improved over the parent and window orders in every tested case/time/cap setting. This broader search used much more additional construction time and a different domain than the subsecond window search, so the comparison is a cost/quality observation rather than an equal-work win for either method.

## Verification, provenance, and resource limits

The complete 263-test suite passed, followed by 66 affected checks and four public-package smoke assertions. Tests include independent dense exponentials, Taylor/atom identities, conserving-window outside interference, exact archived factorization/coverage, and admission refusals. [Independent result checking](evidence/results/INDEPENDENT_CHECK.a1.json) verified 24 archived method outputs and compared 24 archived physical rows against retained independently validated results.

Fresh independent work additionally checked four literal nested-commutator cubic objectives and four unfused full-space physical replays: window and continuation arms for both twelve-mode cases at T=1.5, 10U. The largest replay error difference was below 1.4e-18. This is not fresh independent replay of all 68 physical rows; other new settings use the validated Julia kernel. Full and partial outputs, completed-method checkpoints, failures, and receipts remain retained.

The measured source is bound to the baseline scientific tree `9a19a118247aeca74f968084c74533c1fd529d85`, previously recorded as matching merge `aca1cc33e8e2bf8c61b99775f342800e7791b988`. Measured `ResponseCore.jl` hash `24736a2eb1bbdd9e45bf13c3c2c4e747e9a9672a4239a65f48c2c6df80260d3c` is preserved under `evidence/source_at_measurements/`. The current core adds window-width and cumulative-entry admission guards after timing; these affect invalid/over-cap inputs and do not change the measured eight-unit numerical path. The snapshot makes that distinction reviewable rather than silently attributing old timings to new bytes.

Julia 1.10.11 ran with one Julia/BLAS/GC thread, quarter-core quota, 512 MiB memory, no swap, and memory-headroom monitoring. [Resource accounting](evidence/RESOURCE_SUMMARY.json) includes failed attempts: total 314.11 CPU seconds and 1,291.0 guarded wall seconds, within every cumulative phase ceiling. Preflight used 55.68/60 CPU seconds, so failures were consequential and retained. Language-model usage and unguarded packaging are separate. No GitHub publication, existing optimizer replacement, source CI change, or new general approximation-ratio claim follows from this local study.

The packaged portable legacy verifier was also replayed successfully using repository-relative paths and bundled original model JSON. Original recording adapters remain in `evidence/checks/` as historical source; `tools/verify_legacy_oracles.py` is the portable entry point.
