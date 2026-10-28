# Exact supplied-response bounds (Julia, stdlib only)

This backend constructs independently checkable lower bounds for coherent response norms. It supports exact supplied frozen-window vectors and a separate full-order model using original rational Pauli coefficients plus a total supplied dyadic kernel. A certificate may correctly give zero. `REPORT.md` describes the retained restricted-window results and the full-order feasibility refusal. New driver regressions must be executed separately; they do not change those historical numerical results.

The core API is `src/MolecularResponseBounds.jl`. `ResponseBounds` provides frozen encoding, coordinate/separator bounds and complete comparison covers. `ExactResponse` parses original integers, freezes/loads a total kernel, computes the exact full incumbent UB, builds complete XOR-targeted raw atoms and handles original-word endpoint costs. `StudyAPI.molecular_record` assembles a full certificate; `IndependentCheck` reconstructs all load-bearing algebra separately. Read `PROTOCOL.md` for precise target, proof, caps, failure modes and proposed study gates.

Run only through root's owning resource supervisor, from this directory, with a fresh output directory and labels. Set `JULIA_EXE` to the existing Julia1.10-compatible executable. These are commands to execute later, not recorded results:

```sh
JULIA_EXE=julia
BOUNDS_RUN=/tmp/molecular-bounds-tests-01
python3 ../julia_cubic_feasibility/tools/resource_guard.py \
  --label molecular-bounds-tests-01 --logs "$BOUNDS_RUN/logs" \
  --cpu-seconds 120 --wall-seconds 600 --disk-mib 128 --output "$BOUNDS_RUN" -- \
  "$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=250M --project=. test/runtests.jl
```

After source/protocol/input freeze and a successful resource preflight, use a copy of `registry.example.toml` with paths relative to that copy's location. Options are finite and validated. An existing immutable kernel can be supplied with `kernel="PATH"`; otherwise `freeze_kernel=true` creates a new total table in the output directory. Source coefficients are never replaced by their Float64 conversions. `ideal_budget` or `budget_multiplier` and a declared full `order`/`order_field` select the incumbent. No optimizer is run.

The default Caps profile preserves512 units,4096 words and2,000,000 target probes. To explicitly admit the intended later14-mode panel, select `[caps] profile="larger14"` in the registry, or pass `--caps-profile larger14` to `study.jl`. That fixed profile has1024 units,8192 words and4,500,000 probes; n<=14,50m work,128MiB component storage and every other limit remain unchanged. Only optional `max_units`/`max_words` reductions within the selected profile are allowed in `[caps]`; memory, work, register, probe and other free overrides reject. An explicit CLI profile takes precedence over a registry profile; explicit registry count requests still apply. Both requests and all effective/hard limits are recorded.

```toml
[caps]
profile = "larger14"
max_units = 1024
max_words = 8192
```

Equivalent explicit CLI selection is `study.jl --caps-profile larger14 REGISTRY OUTPUT_DIRECTORY`. Count admission does not prove the complete algorithm fits: scope, support, bit, work, memory and time caps can still refuse. No terms are deleted to fit a profile.

```sh
python3 ../julia_cubic_feasibility/tools/resource_guard.py \
  --label molecular-bounds-dev-01 --logs "$BOUNDS_RUN/logs" \
  --cpu-seconds 360 --wall-seconds 1800 --disk-mib 128 --output "$BOUNDS_RUN" -- \
  "$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=250M --project=. study.jl registry.example.toml "$BOUNDS_RUN/results"
```

The driver preserves each case START, completed-phase records, final result/refusal and complete registry inventory when their declared output caps permit. A final result exceeding its output cap is replaced by a small `bounded-output-refusal-v1` incomplete receipt under the same cap, and subsequent cases continue. The receipt has `result_payload_persisted=false` and claims no verified certificate. If even that receipt does not fit, the exception propagates. Existing outputs and all other I/O errors still fail; no output is overwritten. It exits2 if any case did not complete. This is not a positive-bound gate: completed zero bounds remain successful exact NEGATIVE outcomes. Run independent checks in a separate guarded process (100CPU/500wall for each case; input/kernel/result/output are explicit):

```sh
"$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=250M --project=. independent_check.jl \
  INPUT.toml KERNEL.toml RESULT.toml NEW_CHECK.toml
```

The checker rejects changed executable source hashes. Both supplied input/kernel byte sizes are checked before either user artifact is hashed or read, inside the handled checker path; oversized artifacts return `input_limit` without numerical dispatch. Producer input and supplied-kernel admission also precede hashing. The result file has its own bounded admission before parsing/hashing.

A fully bound refused result requires the producer's available input and kernel digests to match the supplied bytes. Its status is `recorded_refusal`, `artifact_links_authenticated=true`, but `complete=false`, `scientific_verification_performed=false` and `positive_certificate_verified=false`: this authenticates a refusal, not its unperformed mathematics. A legitimate refusal before digest creation, or with a missing artifact, is instead `recorded_refusal_unbound`, with explicit per-artifact binding states and `artifact_links_authenticated=false`; CLI exit2 retains this narrower outcome. `-` may denote an unavailable artifact for such metadata replay. A changed available artifact digest is rejected, and a completed result requires both artifacts to be available and authenticated. The checker uses the same final-output receipt behavior as the producer and continues later quadruples when a receipt fits.

A resource-limited independent replay returns unknown/refused verification, never a verified factor. Its accounting is deliberately stricter than the producer's separate component budgets: source reconstruction, the still-live complete projection and the uncached rational Taylor replay share one monotone `Budget` in `independent_full_check`. The projection is not released before Taylor replay. Thus a legitimate near-cap producer certificate can refuse independent replay under the same declared Caps; there is no claim that every producer proof fits replay. Separate independently bounded replay components, with explicit live-storage accounting, are future performance work. No Budget is reset or relaxed here. Mathematical certification of an exact encoded target does not prove equality with a rounded global scorer, an exact-transcendental objective, or physical infidelity.

The independent CLI reconstructs bounded Caps from the full recorded `caps_declaration`, compares every effective/hard-limit field, and never passes arbitrary producer values into its Caps constructor. `independent_check.jl --caps-profile larger14 INPUT KERNEL RESULT OUTPUT` additionally asserts the expected profile; omitting that flag still validates the declaration rather than reverting larger results to defaults. Normal START/phase/final records carry the declaration, source/registry identity and effective profile; case records distinguish input admission from later refusal. A minimal output-limit receipt can omit that larger metadata and cannot stand in for a fully authenticated producer result. Completed zero bounds undergo the full scientific replay.

Driver-only regressions are included in `test/runtests.jl`. `test/prehash_probe.jl BACKEND_DIRECTORY` additionally exercises the actual checker path against either the immutable historical source snapshot (expected red) or the active source (expected green). Its sparse8MiB+1 file is never read: a method spy intercepts any attempted hash. New tests build fresh current-source records; historical producer records, source hashes and evidence are never rewritten to make them pass a changed executable-source check.

The optional `window_development.jl INPUT START NEW_OUTPUT` constructs one historical eight-label frozen target with unchanged ResponseCore/WindowSearch, checks the new bound independently, and never enumerates permutations. Starts157/168 are the existing CH4/CO development windows. Its J interface requires interior windows and reports `positive_full_bound=false` even if its restricted bound is positive. Do not substitute these old local results for the required full-molecule investigation.
