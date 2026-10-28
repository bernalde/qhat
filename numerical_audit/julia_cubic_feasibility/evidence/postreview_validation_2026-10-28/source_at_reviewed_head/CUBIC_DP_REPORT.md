# Cubic Strang DP integration: correctness and cost

2026-10-28. Follow-up to the first-order DP and separate cubic-window pilot at local commit `5b9e57f5873ec6b6948b70671a5b00d0a6cffb34`. This report addresses the previously unvalidated cubic integration. The original `REPORT.md`, input copies, certificates, and physical results remain intact.

## Outcome

The Julia DP now optimizes the coherent **cubic Strang response objective** on the two existing fixed eight-unit windows. It exactly reproduces both previously exhaustive optima and their unique orders. Fresh complete enumeration of all 40,320 orders per window agrees, including every first and warm repetition.

The performance result is negative for these windows: the cubic DP is about **31–34 times slower** than the existing certified enumerator. Both exact retained factor graphs are complete eight-vertex graphs, so this decomposition requires an eight-label root bag. Its 40,320 root states already cover every window order; smaller bags bring the total to 46,233 states. Factor construction and exact signed-cost evaluation add work without reducing that domain.

Keep insertion search for routine proposals and the certified enumerator for these small-window certificates. Retain cubic DP as an optional exact backend for cases with a demonstrated smaller decomposition. These results do not support making it the default insertion method or claiming a cubic speedup from the old first-order timings.

## Exact target and implementation

For the conditioned window, the supplied frozen complex response is

\[
R(\pi)=b+\sum_S v_S(\pi|_S),\qquad Q(\pi)=\|R(\pi)\|^2,
\qquad J(\pi)=Q(\pi)/r^4.
\]

The baseline b retains the full contribution from the fixed exterior. Each supplied real/imaginary binary64 component is interpreted as its exact dyadic value. A common denominator D gives integer vectors B and V. `CubicDP.jl` constructs the complete signed factorization

\[
N=\|B\|^2+\sum_S\{2B\cdot V_S+\|V_S\|^2\}
+2\sum_{S<T}V_S\cdot V_T,
\qquad Q=N/D^2.
\]

The dot product includes real and imaginary coordinates, equivalently the real part of the Hermitian inner product. The variant-dependent self norm is retained. Pair/triple response factors produce scalar scopes of at most six labels. Only the six realizable triple comparison patterns enter the objective. Equal scopes are aggregated, and exact zeros, constants, and irrelevant labels are reduced without a floating tolerance.

The adapter converts local indices back through the original, generally unsorted global window labels exactly once. It verifies the returned integer DP score against direct exact frozen-vector evaluation. The caller restores the fixed r^4 denominator. A conservative admission-accounting fix in `StructureDP` makes integer denominator one contribute zero denominator-growth bits; regression tests retain refusal for genuinely large rational denominators. This changes admission accounting, not the encoded objective.

The implementation uses Julia 1.10.11 and standard libraries. Factor conversion/construction, graph audit/decomposition, DP state/work/storage, and arithmetic bits all have finite limits. Refusals remain explicit and do not become optimum certificates. The cubic adapter accepts windows of width 1–9; the default bag/state caps may refuse a dense nine-label domain. It is not a full-molecule ordering solver.

## Frozen comparison

`CUBIC_DP_PROTOCOL.md` was written before the new errors or timings. The original CH4_r1p09 and CO_r1p13 input bytes, parent orders, eight-unit windows, HF state, T=1, and budget 30U were retained. Both windows are strictly interior; global first/last units stay fixed and r=15. No window selection, retraining, or case substitution was performed.

| Case | Total units | Window positions | Returned global window order | Exact J as decimal |
|---|---:|---|---|---:|
| CH4_r1p09 | 351 | 157:164 | 152,154,293,331,174,178,221,333 | 2.6658401425367236e-9 |
| CO_r1p13 | 195 | 168:175 | 99,118,90,102,84,131,108,61 | 1.4546310479658738e-8 |

The authoritative values are exact rational strings in `evidence/cubic_dp_validation/results/*MAIN.a1.toml`, not the rounded decimals above. Each result binds its input/prior/source hashes, full order and prepared response-table bit patterns. `provenance/MAIN_FREEZE.sha256` authenticates the pre-main source/input inventory.

A representative preflight first scaled CO coefficients by 7/5 and used fresh complete enumeration as the exact oracle. It matched the DP exactly. Scaling changes diagonal transport, so no simple score-homogeneity assumption was made. The unscaled main cases then reproduced the archived exact rational values and unique witnesses.

## Measured cost

CPU seconds below are medians of three warm calls. Each method starts from the **same prepared FlatWindow**. DP timing includes integerization, scalar-factor construction, graph/decomposition, solving, and witness validation. Enumeration includes its floating screening and exact contender checks. Both paths certify the same exact supplied-table objective using different arithmetic strategies.

| Case | Cubic DP | Complete enumeration | DP / enumeration |
|---|---:|---:|---:|
| CH4_r1p09 | 8.174066 | 0.240049 | 34.1 |
| CO_r1p13 | 5.585027 | 0.182267 | 30.6 |

First calls were 12.849393/0.516070 seconds for CH4 DP/enumeration and 5.946597/0.190036 for CO. They share one process and JIT helpers, so they are not independent cold-process comparisons. Warm calls alternate method order. Explicit pre-call `GC.gc()` is outside each method timer; collection during a method is included. Whole-process resource receipts also charge startup, forced collections, serialization and validation.

Common response preparation is separate: 0.431 seconds for CH4 and 0.029 seconds for CO in this process (see the exact retained fields). Parsing and process/JIT costs are also separate. No already-built-factor-only result is presented as the complete DP cost. Historical same-window insertion also found these unique optima, but its old timings are not a fresh controlled benchmark in this study.

The warm solve dominates the DP cost. Separate phase medians (which need not sum to the median of the total) are:

| Case | Integerization | Scalar factors | Graph/decomposition | DP solve | Witness |
|---|---:|---:|---:|---:|---:|
| CH4 | 0.012428 | 0.023132 | 0.000192 | 8.122419 | 0.005184 |
| CO | 0.004613 | 0.030117 | 0.000120 | 5.544073 | 0.001830 |

Thus avoiding repeated setup would not remove the observed slowdown. Evaluating the exact signed scalar costs on a root bag that already enumerates every order is the principal cost in this implementation.

Both cases retain six-label interactions and 28 exact graph edges. CH4 has 159 retained scalar factors/3,271 entries; unscaled CO has 103/2,759. Their DP estimates admit 46,233 bag states. These describe this particular exact representation and decomposition; they are not lower bounds on every alternative representation or an intrinsic hardness proof. The timing comparison supports an implementation choice on these development cases, not a universal complexity or Julia-versus-Python claim.

## Validation and physical interpretation

Integration tests check complex cancellation and outside interference, negative cross costs, variant-dependent norms, nonmonotone labels, impossible triple columns, zero responses, admission refusals, and a real conditioned Pauli-unit window. A genuine six-label interaction is checked against direct exact-vector evaluation for every one of its 720 permutations.

The affected integration suite passed 1,591 assertions. The final full Julia regression suite passed **1,815 assertions across 27 test sets**, including the existing response, archived first-order, physical-check and entry-point tests. All five guarded attempts exited successfully.

The separate Julia checker reconstructs the molecular tables and factors, compares their exact scalar sums with direct exact-vector scores at 32 deterministic permutations per case, and verifies all eight recorded solver outputs per case. It checks D² normalization, N, archived Q and J, label mapping, fixed exterior/endpoints, actual affordable r, source/input/table hashes, and full Taylor numerical parity. Both cases passed; the recorded/fresh Taylor values and predicted Float64 J had zero difference at the checked precision. It does not call either optimization method. Test/check outcomes and original receipts are retained under `evidence/cubic_dp_validation/`.

The returned complete molecular orders are identical to those in the previous physical study. Both the driver and independent checker bind the original physical archive by SHA-256, case/input, complete order, Strang formula, all four declared time/budget pairs, and affordable step counts. The existing mixed physical outcomes therefore apply unchanged: the CO surrogate gain helps at T=0.5 but loses by about 2% at T=1.5; methane's gains are tiny. No physical propagation was rerun and no new physical advantage is claimed.

The exact certificate concerns the supplied dyadic response tables over these fixed windows. Full Taylor parity is a separate numerical check. Neither gives a real-transcendental, global molecular, physical-error optimum, or approximation-ratio theorem. Larger-window usefulness remains unestablished; a further DP study should first demonstrate a sufficiently small exact factor decomposition.

## Resources and review

All numerical jobs were serialized under the frozen quarter-core, 512 MiB/no-swap, one-thread policy, with finite phase ceilings of 120/100/240/100 CPU seconds for implementation tests/preflight/main/independent checks and corresponding 600/500/1,200/500 wall seconds. Output ceiling: 128 MiB. Failed and successful attempts are retained; none was silently replaced by a fresh allocation.

| Phase | Actual CPU seconds | Actual guarded wall seconds | CPU / wall ceiling |
|---|---:|---:|---|
| Integration tests | 14.88 | 60.17 | 120 / 600 |
| Representative preflight | 20.92 | 84.23 | 100 / 500 |
| Main comparison | 75.50 | 303.73 | 240 / 1,200 |
| Independent checks and full suite | 42.71 | 174.47 | 100 / 500 |
| Total | **154.01** | **622.60** | 560 / 2,800 |

The main job's maximum process RSS was 449,376 KiB, below the 512 MiB limit. `provenance/ACCOUNTING.json` aggregates the unmodified per-attempt receipts. A separate read-only agent reviewed the algebra, mapping, admission guards, tests and comparison boundaries and found no blocking defect. Its original note and dispositions are retained; its final pass caught a stale scaled-preflight factor count in the draft report, which was corrected to the unscaled CO counts. This is local numerical validation and review, not a Lean proof or a new external Claude review.
