# Independent targeted implementation review — round 2

Reviewed `bernalde/qhat` at `a9e8f5047ca3ee5be68fa658abe25779ebd9d7f4`, against merge-base `d8ff197fde4c140eae07dafa500f92b8f5a56d7a`, in `/tmp/qhat-scalable-ordering-2026-10-27/worktree`. The addressing delta is the single commit after the prior reviewed head `e5f6b1fe68c419cf4095a489beed2b5ea9b50c35`. Local HEAD, merge-base, clean worktree, the complete 24-path change inventory, and all 24 recorded SHA-256 hashes agree with `astra-r2/INVENTORY.json`.

**Disposition: the four prior findings are addressed; source and protocol are ready for the declared main run.** There are **0 open Blocking, 0 open Nonblocking and 0 open Question findings** from this targeted review. I found no new actionable issue in the addressing changes. This is a local, no-post assessment under the maintainer-review skill, not human approval, merge approval, or acceptance of unperformed scientific results. No PR exists.

## Finding-by-finding assessment

The classifications and severities below preserve the original finding identities. Each is closed by this review; they are not new open findings.

### ASTRAR1-001 — Blocking — P2 — Addressed
<!-- gh-review-pr:finding=scalar-shift-cubic-cancellation -->

Locations inspected: `analysis/sparse_response_ordering.py:248-272`, `:298-354`, `:372-428`, `:454-523`, `:595-604`, `:664-735`; `analysis/tests/test_response_resolution.py:71-205`; the numerical-refusal boundary in `numerical_audit/scalable_ordering/training.py:283-297`; and the revised response documentation.

The implementation separates the original units from response-only centered units. It removes explicit identity terms and subtracts each unit's diagonal expectation at one fixed occupation. Algebraically, replacing each unit by its scalar-centered version multiplies both the symmetric product and exact exponential by the same scalar exponential. Their coefficients agree below degree three, so the degree-three defect is unchanged. Scalar additions also vanish in the adjacent nested commutators. Public `apply()` and `hamiltonian_action()` continue using the original units, and rotation costs still use their original lengths, including supplied identities. The cache key distinguishes original from centered actions.

The transport kernel now sums supplied-coefficient energy differences directly with rational arithmetic, excluding identity terms, before converting the phase to floating point. This avoids the former subtraction of large total energies. Its LRU cache is explicitly bounded. Source-specific phases remain inside the coherent source sum before projection, so the correction does not replace coherent transport with independent source probabilities.

The two resolution checks operate on the full cubic-column norm and the final projected response norm relative to the pre-subtraction scales. They explicitly refuse the reviewed scalar-free stiff example rather than accepting its unresolved zero response. The commuting-Pauli fast path establishes zero defect algebraically. Training catches `NumericalResolution` separately and retains an earlier fully checked incumbent; incremental proposals still require a full checked score before replacing it. The new analytic tests derive cubic columns from nested commutators, exercise multiple large scalar shifts and coherent references, isolate transport with independently supplied columns, check original operator/cost semantics, and check both first-seed refusal and incumbent retention.

Consequence/remedy assessment: the original scalar-invariance defect and demonstrated unresolved-score failure have appropriate remedies. The documented `2^-24` screen is a conservative numerical refusal policy, **not a forward-error certificate or a proof that every possible cancellation is detected**. It can refuse small responses; that is an explicit supported outcome. Acceptance here does not strengthen its advertised guarantee. The separate literal-matrix objective audit for every final main order remains required before scientific acceptance. No further source change is requested for this finding.

### ASTRAR1-002 — Question — P2 — Addressed
<!-- gh-review-pr:finding=ch4-pauli-control-pair-cap -->

Locations inspected: `numerical_audit/scalable_ordering/training.py:7`, `:120-174`, `:177-253`; `PROTOCOL.md:45-52`; corresponding dictionary, bound, crossing and work-counter tests in `test_pipeline.py`; retained `PAULI_PAIR_COUNT.a1.json` and `pauli-setup-probe-a1/PROBE.json`, their producer scripts, and the probe resource receipt/log.

The independent count identifies **588,360** qualifying singleton-Pauli pairs among **1,481,781** unordered pairs. Its support predicate is appropriate: distinct nonzero singleton Paulis contribute a nonzero commutator column exactly when they anticommute, and nonzero XOR flip excludes the HF destination. Thus the former 250,000-entry limit was insufficient. The replacement stores one encoded-pair dictionary, with `lookup` aliasing that same object, and enforces a 750,000-entry limit before insertion. Decoding the key recovers the original pair orientation; the adjacent crossing update preserves the sign convention. Tests exercise the exact storage boundary and reconstruction parity across small orders.

The retained probe at the reviewed head actually completed the constructor on the frozen altered CH4 input: three seeds, 2,997 incremental candidates, one full move rescore and one accepted move, with all 588,360 entries. The terminal status is `evaluation_limit` with a completed score, not setup refusal. Recorded constructor CPU is 7.057328239 seconds; the service log records 8.005 seconds of guarded CPU, 157.9 MiB cgroup peak, and 179,892 KiB process lifetime RSS under the unchanged 512 MiB cap. These memory/time boundaries are different measurements and should retain their labels.

I read the count/probe scripts and independently matched their hashes, the model/freeze hashes and the returned-order hash to the retained records. The probe imports the study's numerical environment but does not allocate its reference object or evaluate physical errors. Its source map matches the reviewed implementation.

Consequence/remedy assessment: the setup-coverage uncertainty is resolved by measured completed scoring/search, not by increasing the allowance alone. This is evidence for the altered CH4 resource case. Actual main-arm completion, counts, refusals and costs must still be reported from the main constructors; the probe does not establish main accuracy or a completed six-arm study. No further source change is requested for this finding.

### ASTRAR1-003 — Nonblocking — P3 — Addressed
<!-- gh-review-pr:finding=incremental-search-slice-reserve -->

Locations inspected: `analysis/sparse_response_ordering.py:538-615`, `:667-695`; `analysis/tests/test_response_resolution.py:208-287`; and the new search-accounting documentation.

The reserve is now subtracted once from the global remaining search allowance, and full move rescoring may use the original global deadline. The best fully scored seed is visited first. Only its column snapshot is retained, with its combined entry count already checked against `max_entries`; the adjacent update path constructs new mappings and does not mutate that snapshot. Other seed columns are obtained as needed. Separate counters expose completed seeds, incremental candidates, and attempted/completed full move rescoring.

The deterministic-clock regression represents the relevant failure mode: three slices would each be too short if charged their own full-rescore reserve. It now requires a candidate from the best seed, avoids rebuilding that seed's columns, and incorporates an independently rescored improvement within the global deadline. The prior-core log records failure at the lack of any visited candidate. This is meaningful evidence of useful work and incumbent incorporation, not merely survival of a timeout.

Consequence/remedy assessment: the avoidable per-slice reservation and redundant best-seed rebuild are corrected. The shared reserve does not guarantee that every seed or improvement completes on a real molecule, and the documentation/counters make no such claim. No further change is requested.

### ASTRAR1-004 — Nonblocking — P3 — Addressed
<!-- gh-review-pr:finding=standalone-paired-report-verification -->

Locations inspected: `numerical_audit/scalable_ordering/verify.py:14-123`, `:145-161`, `:214-223`; the actual producer in `study.py:148-185`; and the derived-data fixture and negative tests at `test_pipeline.py:790-893`.

The maintained verifier now invokes a pure derived-data gate after the physical/resource checks and before writing its final passed result. It checks exact row/constructor/reference/pair inventories; constructor cost/status/evaluation provenance; all nine emitted pair metric fields; pair status and winner/tie/censor decisions; and summary counts, statuses, cubic comparison aggregates and metric ambiguity. Saved errors must equal their producer fields, while independently recomputed errors must agree within the stated tolerance and yield the same decision. Resources come from the independently audited DAGs. Censored pairs must omit completed metrics. Final verification records hashes of both derived inputs and the fields covered.

The tests use the real `study.summarize()` producer on synthetic rows, then reject changed error, difference, CX, depth, constructor-cost and summary fields without state propagation. They also check authenticated constructor costs, independent resource mismatches and censor handling. I traced the caller's authentication/physical-verification boundary; the pure helper is not being used as a replacement for those checks.

Consequence/remedy assessment: the documented standalone verification path now covers the arithmetic previously left to the external renderer. No additional propagation is needed to test these derived-data guards. Running the maintained verifier on the actual main panel, and checking the final rendered report against those authenticated outputs, remain pending acceptance work. No further source change is requested for this finding.

## Scope, evidence and limitations

This round verified the entire current changed-file identity and read the complete addressing delta across seven paths: the response implementation and documentation, new numerical-resolution tests, protocol amendment, pipeline tests, first-response training implementation, and maintained verifier. Relevant unchanged callers, producer arithmetic, source-binding helpers and CI configuration were also inspected. The prior review's full 23-path source assessment is carried forward only for unchanged portions. This round does not claim a new line-by-line review of every unchanged file.

Retained evidence distinguishes executed producer checks from this static review:

- The red log records eight prior-core regression failures. The identity-inclusive rotation fixture was subsequently corrected to nine steps, consistent with `3*r+2 <= 30`; its original failure occurred earlier at the incorrect cubic column.
- The response/CLI log records 95 passed. The combined five-file cohort records **188 passed in 10.49 seconds**; its seven before/after source hashes match the reviewed addressing files, and its result records source stability and exit code zero. I did not rerun these checks.
- The corrected fixed-seed cubic profile binds the current source and same altered model. It records 0.011937574 CPU seconds for context construction, 3.554213747 for the first score, and 49,084 KiB lifetime RSS. Applying the unchanged preregistered calibration rule still selects 32 CPU seconds for the three searched arms. This is a resource measurement, not a quality or full-search result.
- The earlier two-second six-arm altered preflight remains mechanics evidence at its own frozen source. It is not a new 32-second six-arm study at this head. Prior failed attempts and cumulative costs remain part of the campaign record.
- The broader local attempt is collection-limited: 15 errors from unavailable `tomlkit`, `qualtran` or `pyLIQTR`. It is not a full local pass and cannot replace hosted CI.

This review performed no project execution, tests, builds, installations, numerical/model calls, network access, Git mutations or publication, and spawned no helpers. Its sole write is this report. I did not independently regenerate chemistry, rerun a resource probe, replay a trace, reconstruct every historical artifact hash, evaluate main physical errors, or inspect live GitHub policy/checks. The round-one public report is preserved verbatim; integrator dispositions remain separate.

## Required acceptance gates still open

1. Freeze and authenticate **all 24 main orders before any main validation errors**. Retain actual completed work, fallback/refusal status and full construction charges for every arm. Source readiness and the altered-input probes do not establish completion of those orders.
2. Produce the declared **96 circuit outcomes, eight reference evolutions and 240 paired comparisons**, retaining losses, ties, numerical-floor ambiguity, censors and all attempted-work charges. This main panel is still **unperformed**; no main physical-error result or winner is inferred here.
3. Run independent objective, primitive/state/resource and derived-data verification on the actual main panel. Reconcile the rendered report and claims with authenticated outputs, and complete the final evidence-focused review. The response heuristic remains distinct from the separate conservative certificate.
4. After the deferred publication, obtain full hosted CI for the exact final head. A separately authorized merge phase and verification of actual merge CI remain mandatory. No local review or targeted pass completes those gates.

Preserve this report as the independent round-two public record. Any later dispositions or evidence should be recorded separately.
