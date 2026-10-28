# Independent implementation and protocol review — round 1

Reviewed `bernalde/qhat` at `e5f6b1fe68c419cf4095a489beed2b5ea9b50c35`, against merge-base `d8ff197fde4c140eae07dafa500f92b8f5a56d7a`, in `/tmp/qhat-scalable-ordering-2026-10-27/worktree`. Local HEAD, merge-base and the clean worktree agree with the supplied inventory. No PR exists; this is a local, no-post review under the gh-review-pr skill.

**Disposition: changes required.** The inventory contains **1 Blocking, 2 Nonblocking and 1 Question** finding. The blocking finding concerns numerical resolution in the documented reusable API; the generated molecular benchmark removes its scalar and does not exercise the large dynamic-range examples below. Main scientific acceptance and hosted CI remain pending independently of the findings.

## Findings

### ASTRAR1-001 — Blocking — P2: prevent unresolved Taylor cancellation from becoming a completed cubic score
<!-- gh-review-pr:finding=scalar-shift-cubic-cancellation -->

Locations: `analysis/sparse_response_ordering.py:214-220`, `:229-239`, `:291-315`, and `:350-365`; public identity-term support is documented in `analysis/RESPONSE_ORDERING.md:11-17` and `molecular_units` at `:128-159`.

The response calculation carries identity coefficients through the Hamiltonian powers and product Taylor recurrence, then subtracts the two large cubic columns. Adding a scalar to a unit must leave the cubic defect unchanged: both product and exact evolution acquire the same scalar exponential, and their first difference is degree three. The current arithmetic needlessly loses that invariance.

A concrete supported input is one qubit, reference `0`, with units `A = 2**40 I + X` and `B = Z`, in order `(A, B)`. Both units satisfy the internal-commutation contract. Independently,

`C = -[X,[X,Z]]/24 + [Z,[Z,X]]/12 = -Z/6 + X/3`,

so the off-diagonal cubic column is exactly `1/3`, irrespective of the scalar. Before subtraction, the exact product and Hamiltonian Taylor off-diagonal entries are respectively `2**79 + 2/3` and `2**79 + 1/3`. Binary64 cannot resolve either small remainder at that magnitude. Subtracting their floating representations cannot recover `1/3`; a zero or a large cancellation residual is reported instead. The resulting completed score can therefore have no meaningful relation to the stated cubic objective. The adjacent-update commutators also unnecessarily carry the same scalar through their subtractions.

Identity centering alone is insufficient for the full numerical-resolution problem. For the scalar-free units `A = M X`, `B = Z`, the exact defect is `-M**2 Z/6 + M X/3`. The product and Hamiltonian off-diagonal cubic terms are `M**3/6 + M/2` and `M**3/6 + M/6`. With `M=2**40`, their floating subtraction cannot resolve the required `M/3` component either. Taking evolution time `T=1/M` keeps `T*||H||` approximately one and the exact transported perpendicular amplitude approximately `1/3` (squared norm `1/9`), so this failure does not require an enormous dimensionless evolution time. This scalar-free example was raised as a coordinator hypothesis and independently checked here against the same nested-commutator identity; it is not an executed result.

These are algebraic/representation counterexamples established by static inspection, not executed reproductions. The supplied tests exercise moderate coefficient scales and only a small scalar coefficient; they do not protect scalar-shift invariance or this dynamic-range case.

Remedy: use a stable sparse defect evaluation, a checked fallback, or an explicit numerical-resolution refusal when the computed defect is not resolved; do not silently return it as a completed score. Direct sparse commutator/nested-commutator evaluation is one candidate design, not a verified replacement prescribed by this static review. Remove identity contributions from numerical response actions and diagonal-energy differences as part of preserving scalar invariance, while keeping original unit identities and declared rotation-cost accounting separate. Keep the independent conservative certificate's contract explicit. Add regressions for both examples and identity shifts across multiple units, checking cubic columns and coherent response against an independent analytic/commutator oracle, including the adjacent-update path. Validate the chosen fix and its resource behavior before claiming this finding resolved.

### ASTRAR1-002 — Question — P2: establish that the CH4 individual-Pauli control can reach a scored seed
<!-- gh-review-pr:finding=ch4-pauli-control-pair-cap -->

Locations: `numerical_audit/scalable_ordering/training.py:7`, `:122-147`, `:176-194`, `:238-241`; matched searched-arm claim in `PROTOCOL.md:22-29`.

The generated CH4 input has 1,722 Pauli words, hence 1,481,781 possible unordered pairs. `FirstResponse` retains at most 250,000 nonzero, non-HF-destination pair entries and builds that complete structure before evaluating any seed. Reaching the next qualifying pair produces `support_limit` and returns an initializer with no score. Increasing the constructor allowance from two to 32 CPU seconds cannot fix a storage refusal.

The retained two-second CH4 preflight has `ordering_status="cpu_limit"` and zero evaluations for the Pauli arm. It does not establish whether complete setup fits the pair bound. I did not execute a pair count or constructor, and do not claim that the actual qualifying count exceeds 250,000. That count/setup outcome is the material unresolved question. The 32-second calibration measured a cubic seed, not the Pauli pair structure.

Remedy: before interpreting the main panel as including a searched strong Pauli control, obtain a bounded, resource-only qualifying-pair count or setup-completion result on the frozen CH4 input, without reading physical errors. If the cap prevents scoring, provide a memory-bounded implementation that can score/search this control, or explicitly retain and report the arm as a refused initializer and narrow the strong-control conclusion. Any method/protocol change must precede the main order freeze. Record completed seeds, search candidates and refusal status separately in the later evidence review.

### ASTRAR1-003 — Nonblocking — P3: avoid reserving a full rescore inside each already-small seed slice
<!-- gh-review-pr:finding=incremental-search-slice-reserve -->

Locations: `analysis/sparse_response_ordering.py:435-445`, `:451`, `:545-559`.

After all seed scores complete, the remaining allowance is divided by the number of scored seeds, but every slice subtracts `1.25 * longest_seed` as a full-rescore reserve. With three seeds, the first seed cannot even enter its column calculation when one third of the remaining allowance is below that reserve, even if the total remaining allowance could support useful work plus a rescore. Only the last scored order has cached columns, so earlier seeds can additionally pay for another full column calculation. The largest-case retained resource probe measured about seven CPU seconds for the initial score, making this a relevant risk for the 32-second main allowance, although it does not establish the eventual main search counts.

This preserves correctness and the deadline but can make an insertion arm primarily select among initializers while leaving usable allowance unspent. The existing large-list regression confirms some incremental candidates on a tiny reachable subspace; it does not cover this timing regime.

Remedy: consider a shared rescore reserve, a minimum viable slice, or redistributing unused slices, with a deterministic-clock test showing that available time can fund a checked improvement. A documented deferral is reasonable. In all cases, report actual seed/candidate/rescore counts and avoid inferring successful insertion work from the configured budget alone. Do not retune using main physical-error outcomes.

### ASTRAR1-004 — Nonblocking — P3: make the standalone verification result cover the reported pair fields
<!-- gh-review-pr:finding=standalone-paired-report-verification -->

Locations: `numerical_audit/scalable_ordering/verify.py:99-121`, compared with the fields emitted at `study.py:159-185`.

The maintained verifier independently checks circuit resources and state errors, but its final `PAIRS.json` comparison checks only each winner/tie/censor decision. It does not check the saved pair error values/differences, CX/depth fields, construction-cost sums, or `SUMMARY.json`. Those derived fields can be inconsistent while this command still writes `status="passed"`.

The separately prepared `/tmp/qhat-scalable-render-report.py` already contains cross-file checks for these fields and the summary (inspected at lines 134-175 and 290-323). This narrows the concern: that planned reporting gate can protect this campaign if it is retained, source-bound and actually run. It is not part of the reviewed 23-path change or the documented standalone verification command.

Remedy: include the derived-data checks in the maintained verification/reporting path, or explicitly document the additional mandatory reporting gate and its narrower division of responsibility. A small test that changes one pair resource/cost field or one summary count should fail the final gate. No repeat propagation is needed merely to check this arithmetic and provenance.

## Coverage and mathematical assessment

All 23 changed paths were inspected: the workflow change; changed passages of the root and analysis READMEs; the complete new response-ordering documentation; `analysis/sparse_response_ordering.py`; `analysis/train_ordering.py`; all three new analysis test files; and all 14 files under `numerical_audit/scalable_ordering` (`PROTOCOL.md`, `README.md`, `REGISTRY.json`, `common.py`, `compilation.py`, `compiler.py`, `generate_inputs.py`, `identity.py`, `preflight.py`, `reference.py`, `study.py`, `test_pipeline.py`, `training.py`, and `verify.py`). Relevant unchanged definitions inspected include `analysis/strang_schedule.py`, `analysis/fermionic_coefficients.py`, the parent construction/commutation/decomposition definitions in `numerical_audit/followup_classes/exact_blocks.py`, package initialization, test import setup, packaging, and the earlier coherent cubic response definition.

The review traced the third-degree Taylor recurrence, the adjacent-swap identity `[(A+B)/8 + suffix/4,[A,B]]`, source-specific complex transport and projection, fused rotation affordability, and the independent exact-rational step bound. I found no algebraic discrepancy in those inspected definitions; ASTRAR1-001 concerns their floating evaluation. The certificate's loose bound and its state-metric conversion are separate from the response heuristic. No global ordering, finite-error improvement, chemistry-model accuracy, or hardware guarantee follows from the heuristic.

The compiler/replay review covered Y phases, least-significant occupation bits, routed CX/SWAP restoration, inverse-Clifford images, HF final frame up to global phase, rational angles, adjacent cleanup, primitive counts, dependency depths, and bounded state/trace inputs. Search and campaign inspection covered complete permutations, deadline/support fallbacks, completed-score retention, initialization/setup/cleanup charges, the three searched-arm starts and allowances, freeze-before-errors sequencing, model/source/environment identities, exact molecular partitioning, number-conservation checks, and prospective comparison/censor handling. The added CI step runs the small pipeline regression file; it does not itself launch the chemistry generation or the 96-circuit study.

## Evidence used and limits

This review performed **no project execution, tests, builds, package installation, network access, Git mutation or publication**, and spawned no helpers. Its sole write is this report. The assigned 900-second inspection limit and read-only execution restriction were respected. Static inspection and retained producer evidence are distinct; I have not independently rerun any numerical claim.

Evidence read included `TASK.json`, the supplied brief/inventory, the targeted test result/source-input records and `scalable-final-tests-a5.log`/resource receipt, the broader collection-failure record/log, the preflight freeze/completion and continuation records, the three paired initializer-profile trials, the cubic-seed resource profile, and the full-size stream-probe report. The targeted log records **157 passed in 6.50 seconds** at the reviewed source. The broader local attempt ended with **15 collection errors** for missing `tomlkit`, `qualtran` or `pyLIQTR`; it is not a broader test pass.

The six-arm altered-input preflight was frozen at `5a4a1916f0855844862f45a8e06572f186b8dcef`, used the earlier two-second training allowance, and completed through retained attempts plus a same-source continuation. It is not a fresh six-arm run at the reviewed head or at the 32-second allowance. The continuation reused four circuits, completed two, preserved the interrupted trace, and explicitly excluded earlier attempt costs from its own 116-second stage total; cumulative campaign cost must come from the retained external ledger. The later stream probe records byte-identical compressed output and state/DAG parity for one 1,840,120-gate trace at the new stream implementation. Those facts support bounded reuse of the mechanics evidence, not completion of the amended main comparison.

The initializer record retains three paired repetitions and equal returned orders; its median paired old/new CPU ratio is 8.3566. That is a measurement of `current()` on one altered CH4 input, not the full trainer or full campaign. The resource-only cubic probe measured 7.0312 CPU seconds for one first seed and about 52 MiB lifetime process RSS; it supports the preregistered shared allowance selection, not successful search or accuracy at that allowance.

Omissions: no complete rehash/replay of all retained evidence bytes; no independent re-generation of molecular tensors; no executed scalar-shift regression; no actual CH4 qualifying-pair count; no completed main-result or rendered-report review; no live hosted checks or repository merge-policy assessment. The prepared external report renderer was inspected only at the relevant authentication/derived-pair checks, not reviewed in full as a delivered artifact.

## Required acceptance gates still open

1. Resolve ASTRAR1-001 and the material control-coverage question, and verify the relevant changes before freezing the main orders. Nonblocking dispositions may be documented without requiring an unnecessary respin.
2. Produce and authenticate **all 24 main orders before any main validation errors**, then retain the declared **96 circuit outcomes, eight reference evolutions and 240 paired comparisons**, including losses, ties, censors, fallback orders and all attempted-work charges. None of that main panel is complete at this review.
3. Complete the independent main verification and the final derived-report/provenance checks. Review the resulting evidence and claims before publication; the altered preflight is not the larger-molecule study.
4. After publication, obtain the required full hosted CI evidence for the final head. The targeted local pass and collection-limited broader attempt do not replace it. The separately authorized merge phase and verification of actual merge CI also remain outside this review's completed work.

Preserve this report verbatim as the independent public review record; record integrator dispositions separately.
