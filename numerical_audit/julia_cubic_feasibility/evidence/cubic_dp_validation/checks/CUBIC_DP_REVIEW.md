# Independent review: integrating exact cubic response tables with bag-order DP

Read-only scientific audit, 2026-10-28. Coordinator supplied base `5b9e57` and owns all execution/integration. This review reads the existing Julia `WindowSearch.FlatWindow`, `ResponseCore`, and `StructureDP` interfaces. No Julia/numerical execution, source edit, Git operation, external write, or helper was used. This note is the only write by this reviewer.

## Exact mapping that the integration must preserve

For one supplied conditioned contiguous window, the frozen response is

\[
R(\pi)=b+\sum_S v_S(\pi|_S),\qquad |S|\in\{2,3\}.
\]

The complex vector b is the **full baseline response**, including all fixed-outside interference. Each delta is the exact dyadic value of a supplied ComplexF64 table component; this is not exact transcendental transport. Both real and imaginary coordinates must be retained, including tiny nonzero components whose complex magnitude is dominated by the other coordinate.

Every finite binary floating component is rational p/2^e. Select the common denominator D as the largest component denominator (or their exact LCM); the maximum suffices only because all denominators are powers of two. Integerize real and imaginary coordinates separately. For the resulting integer vectors B and V_S,

\[
N(\pi)=\left\|B+\sum_S V_S(\pi|_S)\right\|_2^2,
\quad Q(\pi)=N(\pi)/D^2.
\]

Expand N into signed local scalar tables:

\[
c=\|B\|^2,\quad
h_S(\sigma)=2B\cdot V_S(\sigma)+\|V_S(\sigma)\|^2,
\quad
h_{ST}(\tau)=2V_S(\tau|_S)\cdot V_T(\tau|_T).
\]

The dot product is on integer real/imaginary coordinates, equivalently the real part of the Hermitian complex inner product. A bilinear complex product without conjugation would be wrong. The self norm in h_S is essential: unlike the archived first-order +/- pair-vector family, cubic variant norms need not be constant. Cross scopes have size at most six. Tables with the same scope may be summed exactly; exact constants and any zero/irrelevant-label reductions must preserve c.

If a positive common factor G is divided out to reduce integer sizes, the returned objective must restore it. For example a global gcd of **all response coordinates** gives Q=G^2 N_reduced/D^2, not G N_reduced/D^2. A gcd of the already expanded scalar costs instead scales N linearly. A constant offset C removed before optimization must also be restored. The simplest initial implementation should retain unscaled N and D explicitly.

For the existing interior eight-label domains, the original global first/last labels are fixed, hence the affordable Strang count r is constant under the declared 30U cap. Therefore

\[
J(\pi)=N(\pi)/(D^2r^4).
\]

The physical/reference time T is already inside b and V; do not insert another T power. The coefficient denominator of the molecular Hamiltonian has already been used when forming the frozen vectors and must not be divided out a second time. In particular this is not the archived first-order JD denominator mapping.

## Label and table conventions

`FlatWindow.labels` stores original global unit labels in baseline-window order. `FlatFactor.scope` and the order accepted by `WindowSearch.exact_numerator` use local labels 1:k. These need an explicit bijection; the baseline window is generally not sorted. A witness cannot be mapped by global numeric rank unless every factor and output applies that same map. Preserve global labels for reconstruction and old-winner comparison.

Three-label factors use eight pair-comparison bit slots, of which only six are realizable total orders. Only realizable columns participate in scalar tables. Zero-filled impossible slots do not represent extra allowed states. A complete local order table must evaluate each union permutation consistently on both subscopes.

Destination coordinates are sparse, but global baseline support can be much larger than each atom support. Dropping coordinates untouched by all deltas is safe only after retaining their exact contribution sum(|B_d|^2) as a constant. A coordinate touched by even one delta keeps its full baseline coordinate inside the coherent square. Overlap buckets can avoid zero dot products, but potential support overlap is not an exact nonzero-factor certificate.

## Primary validation domain and acceptance criteria

The two existing CH4_r1p09 and CO_r1p13 eight-label windows are sufficient to validate this integration. They already have unique exact frozen-table optima, full 8! coverage, fixed outside order, T=1 and r=15. A larger nine-label cubic study is not necessary merely to make the work appear broader.

For each case:

1. Authenticate the existing model, baseline parent, window positions, time, original rotation budget/count, and old certificate. Reuse the same FlatWindow instance for DP and fresh enumeration. Retain a canonical digest of its baseline and each realizable delta component.
2. Check exact equality between the old rational optimum and fresh `exact_numerator` at the old witness. Match the old unique global witness and exact optimum again with a fresh bounded complete enumerator.
3. At a declared small set of deterministic orders (baseline, old winner, reverse window, and seeded examples), compare the expanded scalar factor sum against D^2 times `exact_numerator`, using exact integers/rationals with zero tolerance.
4. For a completed DP, verify a genuine global permutation of precisely the same local labels, fixed outside prefix/suffix, exact restored N, exact N/D^2 equality to direct frozen-vector rescoring, and exact J mapping. Because the previous optimum is unique, the returned global witness should equal it; no new DP minimizer-count claim is needed.
5. Keep the independent full Taylor/literal cubic check as a secondary numerical parity check with its prior tolerance. Such a check does not upgrade the dyadic frozen-table certificate to a real-transcendental or physical-error certificate.

Tiny tests should include complex interference (a nonzero imaginary component), varying self norms, signed cross costs, disjoint supports, a genuine six-label union, nonmonotone global labels, all-zero/constant cases, and explicit preallocation refusal. These expose errors the two molecular winners alone may miss.

## Resource and performance interpretation

Integerization is linear in the stored response component count, but scalar expansion can be quadratic in the factor count and require up to 6!=720 local orders per union. Width eight gives 84 vector factors, hence 3,486 unordered factor pairs before sparsity. This is bounded but not automatically inexpensive. Charge construction and storage before appending large scalar tables, and check numerator/denominator bit limits before state allocation. Common denominator exponent can reach 1,074 for subnormal binary64 components; squaring doubles that scale, and conversion through Float64 would lose exactness.

A correct generic DP bit bound should charge denominator one zero bits of denominator growth. For arbitrary unrelated rational denominators, a bound based on the sum of ceil(log2 denominator) is safe for the common denominator product; integer costs should not acquire one fake denominator bit per table entry. Exact integer scalar factors allow a particularly simple bound: largest cost numerator bits plus log2(number of selected factors) and sign margin. Admission estimates remain conservative forecasts, not process-RSS guarantees.

Compare complete elapsed/CPU work from the same prepared FlatWindow: integerization, scalar factor construction, graph/decomposition, DP solve and direct validation versus enumerator preparation, all-candidate floating screening, exact contender checks and validation. Either charge shared conditioning to both complete methods or present it as shared prerequisite plus each incremental charge. Warm both paths symmetrically and retain first-call/JIT costs separately. A refusal must retain structural/resource diagnostics and cannot be reported as zero score or completed optimization.

The existing Float-screened enumerator and integer DP can certify the same exact frozen optimum, but their arithmetic workloads differ. A timing comparison is still meaningful as an end-to-end implementation decision. It cannot isolate algorithmic gains under identical arithmetic without an all-integer enumeration control. Such a control is optional unless that stronger attribution is claimed. A slower or refused DP on a dense eight-label coherent graph is a useful negative feasibility result; it does not invalidate the exact mapping or imply that every factorization/decomposition is inefficient.

## Scope of any successful outcome

Success establishes that supplied cubic pair/triple vector tables can be integrated with the exact scalar-table bag DP on these fixed windows. It does not establish an unrestricted molecular approximation ratio, a global full-order optimum, usefulness on larger molecules, exact transcendental kernel values, or improved finite-time physical error. The existing mixed physical results remain applicable to the identical returned orders; there is no need to repeat physical propagation merely because a different exact optimizer reproduces an authenticated same-order witness.

## New-module audit status

Read `src/CubicDP.jl` after it arrived. **No blocking algebraic, normalization, or label-mapping defect found.** This is a static audit, not a passing test result.

The implementation's valid triple columns `(1,2,4,5,7,8)` agree with the pair-comparison encoding; the other two slots are not used. It includes baseline coordinates and only realizable factor variants in the common denominator. `divrem` checks exact divisibility during integerization. The real/imaginary sparse dot product is correct, the variant-dependent self norms are included, and cross terms are doubled exactly. Scope grouping and later constant/irrelevant-label reductions preserve the full signed factor sum. Direct integer DP cost is required to be integral and nonnegative, then N/D^2 is checked exactly against `WindowSearch.exact_numerator` before local indices are mapped once through `flat.labels`. Division by fixed r^4 remains explicitly outside the adapter.

The constructor has factor-entry, exact-vector entry, integer-bit, estimated-memory, and work budgets; the generic bag DP retains its own state/work/memory admission. It returns explicit refusals rather than a zero objective or a certificate on incomplete work. For the study's 84-factor/eight-label inputs, scalar table scopes are at most six and each new local table has at most 720 entries. Cross construction caches at most 36 scalar dots for one factor pair at a time rather than a full dense Gram tensor.

One minor public-input resource point was sent to the author: the first version of `_validate_flat` charges work after scanning factor headers and their destination uniqueness checks, while an arbitrary list of zero-destination factors does not contribute to its scalar-entry count. A per-header guard/work charge or explicit response-factor count bound closes this gap. It does not alter the authenticated 84-factor scientific domain. The author/coordinator owns any patch and tests.

The coordinator clarified that its representative preflight scales physical Hamiltonian coefficients by 7/5 and uses a fresh exact enumeration on the resulting same FlatWindow as oracle. It assumes **no** score homogeneity. That is sound: fixed-T diagonal transport changes when H is scaled. Main unscaled cases still require exact old rational optimum and unique old witness reproduction.

## Follow-up review after successful scaled-CO preflight

Re-read the current adapter and `cubic_dp_study.jl` after the coordinator reported the scaled-CO preflight: complete 40,320-order enumeration, exact DP agreement, width seven, 46,233 bag states, roughly 10.44 DP CPU seconds versus 0.43 enumeration CPU seconds. These numbers are coordinator-reported preflight observations, not computations run by this reviewer.

The adapter's prior header-loop admission issue is now closed: factor headers, destination uniqueness inputs, and finite-component scans incur work/memory admission before conversion. The graph/decomposition refusal reasons are also explicit. The mathematical mapping remains sound. Tiny tests now cover imaginary cancellation, a true six-label cross scope over all 720 orders, nonmonotone global labels, unreachable triple slots, an actual conditioned Pauli window, zero response, and resource refusals. I inspected those tests but did not execute them.

The driver uses one prepared FlatWindow for both methods, preserves a canonical raw-bit component digest, obtains a fresh exact preflight oracle without a scaling assumption, checks direct exact dyadic scores, and requires exact old rational optimum/global witness equality on main inputs. Both methods receive post-construction witness rescoring inside their measured call. Shared conditioning and a full Taylor parity check are reported separately. Per-method explicit `GC.gc()` is outside the measured call; whole-process resource receipts must retain that cleanup charge, and timing descriptions should identify this boundary. The current preflight already supports an honest negative speed result rather than a claim that exact DP is generally faster.

Before main, the coordinator plans to strengthen full-order and physical-reuse bindings. A local-window witness match alone is insufficient to label old physical results reused: require the authenticated same model bytes/source, entire reconstructed order, original word-count cost/r, and the same physical times/caps. Preflight physical H scaling must never inherit unscaled physical errors.

### Lightweight independent main-result check

A separate standard-library TOML/Fraction/SHA pass is enough to check the new integration's immutable evidence without another propagation study:

* Verify input/prior hashes and one common FlatWindow digest across all completed records of a case.
* Verify each local witness is a permutation of 1:8, map through the recorded baseline-window labels exactly once, and reconstruct the whole order with unchanged outside positions. Match the old unique complete-order witness exactly.
* Parse DP integer N and vector denominator D directly from strings; assert recorded denominator equals D^2, DP table cost is integral and equals N, N>=0, N/D^2 equals the old exact rational Q, and recorded J equals Q/r^4.
* Recompute r from the authenticated original unit word counts and unchanged global endpoints under 30U; check the main cases' r=15. Recheck fresh enumerator rational objective and 40,320 coverage against the same old certificate. Do not infer a new DP minimizer count.
* Bind any physical-reuse statement to the old physical file hashes and matching complete orders/parameters. Existing authenticated literal-cubic and full-space physical checks at the identical witnesses can be reused; no new quantum-state propagation is required by the optimizer change.

This verifies evidence transformations independently of the new scalar factor implementation. Exhaustive tiny-table identities and the derivation above supply the complementary implementation check. Optional component snapshots would permit a separate direct integer-vector scorer at selected orders, but they are not necessary to justify the two already-certified primary domains or to avoid a larger nine-unit study.

## Final main/check/report review and dispositions

Inspected the final `CUBIC_DP_REPORT.md`, `cubic_dp_study.jl`, `independent_cubic_dp_check.jl`, both `*MAIN.a1.toml` results and both `*CHECK.a1.toml` results. **No blocking scientific or implementation finding remains.** This reviewer performed passive source/evidence reads and file hashing, but no numerical execution.

Resolved findings:

* **Header admission guard:** closed. The current `_validate_flat` charges header/destination/finite-component work and admits relevant storage before integer conversion; the formerly uncharged empty-factor header scan is no longer present.
* **Physical reuse binding:** closed. Both study driver and checker require the authenticated retained physical file, exact input/case, complete reconstructed winner, formula two, all four declared time/budget pairs, original rotation budgets and affordable step counts. The main results record `physical_reuse=true`, the full order, four bound rows, and the physical file hash. Preflight scaling cannot claim this reuse.
* **Source binding:** closed. Current file hashes match the independent-check inventories: adapter `d1ad1e9f...`, generic DP `1ececadd...`, study driver `5b0d603c...`, checker `b914779e...`. The old response snapshot hash `24736a2e...` is checked against each prior certificate, distinct from current guarded ResponseCore `636ce3a0...`. The generic denominator-one accounting fix retains the previous conservative charge for other denominators and does not alter encoded costs.
* **Report scalar-count mismatch:** one low-severity stale preflight count was found and reported. Unscaled CO main/check results have **103 retained scalar factors and 2,759 entries**, not 105/2,763. The coordinator reports that the report has been corrected. CH4's 159/3,271 counts were already consistent.

Evidence and interpretation checked:

* Both main cases have four completed DP and four completed enumeration records, the same per-case FlatWindow digest, the exact old rational Q/J and unique global witness, width seven, 28 actual factor-graph edges, and 46,233 bag states. D^2 restoration and direct dyadic witness checks are explicit.
* The separate checker performs 32 deterministic exact scalar/vector identity checks per case, checks all eight recorded solver outputs, reconstructs the full fixed-outside winner, and verifies physical reuse. It calls neither optimization method. It shares the factor constructor under test, but compares its sums to the independent direct frozen-vector evaluator at nonoptimized orders; this limitation is correctly described rather than presented as a separately reimplemented optimizer.
* The report's warm medians and individual phase medians match the checker summaries: CH4 DP/enumeration 8.174066/0.240049 CPU seconds, CO 5.585027/0.182267. The stated roughly 34.1x/30.6x slowdown is valid for incremental method calls from the same prepared FlatWindow. Shared preparation, forced GC outside each timer, shared JIT, per-phase medians, and different arithmetic strategies are disclosed. The result does not imply a language or universal algorithmic speed claim.
* Complete eight-vertex factor graphs explain the eight-label root state space for this representation. The report does not promote this into an intrinsic lower bound for every alternative representation or full-molecule hardness.
* Prior mixed physical findings are preserved: identical CO window orders help at the shorter time and worsen error by about 2% at the longer time. No new physical propagation or advantage is asserted.

The coordinator reports the final full suite passed **1,815 assertions in 27 testsets**, all five jobs exited zero, and total retained study usage was **154.01 CPU / 622.5989 guarded wall seconds**. The independent check's 64 exact sample identities passed and its reported Taylor differences were zero at Float64 precision. These final execution outcomes were supplied by the coordinator; this reviewer did not rerun them. Resource/test-count entries in the report are being finalized from those receipts.

Remaining limitations are material but explicitly scoped: exactness is for the supplied frozen dyadic tables and fixed eight-unit domains; real-transcendental transport, global molecular optimum, approximation ratio, larger-window usefulness and new physical accuracy are not established. The current cubic DP is substantially slower than the certified enumerator on both primary cases, which is a valid negative efficiency result.
