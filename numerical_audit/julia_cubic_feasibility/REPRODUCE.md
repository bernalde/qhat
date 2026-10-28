# Reproduction

The Julia package uses standard libraries only. Tested runtime: Julia 1.10.11. The study's numerical methods require no Python installation. The optional resource supervisor and retained legacy verification use the established Python environment.

Run from this directory. Choose a fresh output directory and fresh attempt labels; scientific files and receipts are never overwritten. Replays are additional work: retain their cumulative costs rather than resetting a partially consumed allocation.

Set `JULIA_EXE` to an existing compatible Julia binary (the recorded runtime was Julia 1.10.11). Set `REPLAY_DIR` to a new scratch output directory. These shell variables describe this replay only; they do not change saved Julia configuration. Original machine paths remain in the historical receipts.

```sh
JULIA_EXE=julia
REPLAY_DIR=/tmp/julia-ordering-replay-01
mkdir "$REPLAY_DIR"
```

The retained `tools/resource_guard.py` requires the Linux user systemd bus. It enforces quarter-core CPU, low priority, 512MiB memory, no swap, one numerical thread, finite wall/output limits, and memory-headroom monitoring. It does not call an external language model. The Julia CLI additionally sets one Julia/GC thread and a 300 MiB heap-size hint. Run the commands sequentially.

## Tests

```sh
python3 tools/resource_guard.py --label julia-dp-tests-replay-01 \
  --logs "$REPLAY_DIR/logs" --cpu-seconds 60 --wall-seconds 300 \
  --disk-mib 128 --output "$REPLAY_DIR" -- \
  "$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. test/runtests.jl
```

## Exact nine-unit comparison

```sh
python3 tools/resource_guard.py --label julia-dp-archive-replay-01 \
  --logs "$REPLAY_DIR/logs" --cpu-seconds 140 --wall-seconds 700 \
  --disk-mib 128 --output "$REPLAY_DIR" -- \
  "$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. archive_study.jl main "$REPLAY_DIR/ARCHIVE.toml" \
  evidence/inputs/archive/F2_r1p55.toml \
  evidence/inputs/archive/H2O_scale1p10.toml \
  evidence/inputs/archive/Li2_r2p95.toml
```

Expect every method result to be complete and its `validation.valid` to be true. Exhaustive runs must have 362,880 visits and complete Lehmer coverage. Minimizer counts are 13,440/44/264. Inspect per-method checkpoints if a later step fails; an empty final case list does not imply no work occurred. Shared-process first calls are not independent cold-process benchmarks. Timings will vary across machines.

## Cubic window and physical checks

The portable input copies are under `evidence/inputs/cubic/`. Under the same guard, invoke:

```sh
"$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. window_batch.jl window \
  evidence/inputs/cubic/CH4_r1p09.toml "$REPLAY_DIR/CH4_WINDOW.toml" \
  evidence/inputs/cubic/CO_r1p13.toml "$REPLAY_DIR/CO_WINDOW.toml"
```

Use a 20 CPU / 100 wall second allocation for this recorded-size batch. `structure` replaces `window` for the capped diagnostics; the diagnostic's own caps/statuses still apply. Large widths of partial potential graphs are not exact coherent-width certificates.

For physical replay, use `physical_study.jl archive evidence/inputs/archive_models "$REPLAY_DIR/ARCHIVE.toml" OUTPUT.toml`, or `physical_study.jl cubic evidence/inputs/cubic WINDOW_RESULT... OUTPUT.toml`. The archived model manifest binds model/domain bytes and parent orders. Cubic mode checks input and parent identity before simulating. To include a recorded legacy continuation control, put its unchanged `CASE.CONTINUATION.a1.toml` beside the corresponding window result. Preserve and inspect its `cpu_limit` status rather than treating it as an optimality certificate.

## Optional independent legacy verification

`tools/verify_legacy_oracles.py` uses packaged original model JSON and the repository's unchanged legacy exact/matrix/replay implementations. It defaults to this study directory and its enclosing QHAT repository; `--study-root` and `--repo-root` permit explicit paths. Use the existing pinned Python environment with NumPy/SciPy, not a new dependency installation. Run it inside the resource guard:

```sh
"$LEGACY_PYTHON_EXE" tools/verify_legacy_oracles.py \
  --repo-root "$QHAT_REPO_ROOT" "$REPLAY_DIR/LEGACY_CHECK.json"
```

Set `LEGACY_PYTHON_EXE` to the existing pinned Python environment with NumPy/SciPy and `QHAT_REPO_ROOT` to the full QHAT source checkout containing the retained legacy implementations. This explicit root is required when the study directory is copied outside QHAT. It checks the packaged recorded results; fresh legacy replay costs belong to its own bounded allocation. The copies in `evidence/checks/` are the original recording adapters, retained verbatim for provenance, and contain historical scratch paths. The portable tool does not depend on those paths. A standalone Julia-only bundle can reproduce the new solvers; optional legacy oracles additionally need the pinned QHAT source tree.

Verify the retained archive with `sha256sum --check MANIFEST.sha256`. The manifest excludes itself and records all current source, input, output, log, and receipt files. `evidence/source_at_measurements/ResponseCore.jl` retains the exact pre-hardening source used by the timing/window records; the later guard-only difference is described in the report.

## Cubic DP integration validation

This is the subsequent study in `CUBIC_DP_PROTOCOL.md` and `CUBIC_DP_REPORT.md`. It uses the same fixed eight-unit CH4 and CO windows as the original cubic pilot. Their full order, endpoints, input bytes, and r=15 remain fixed; the DP now optimizes the coherent cubic response square rather than the archived first-order objective. Its certificate is exact for the supplied binary64 response tables.

Run the following command through the resource guard with a fresh output directory and a 230 CPU / 1,150 wall second limit. The recorded main run includes one first call plus three warm calls per method and case. Each DP call includes factor construction and graph/decomposition, with separate phase fields. The main results contain source hashes and a bit-pattern digest of their common frozen response tables.

```sh
"$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. cubic_dp_study.jl main \
  evidence/inputs/cubic/CH4_r1p09.toml evidence/results/CH4_WINDOW.a2.toml "$REPLAY_DIR/CH4_MAIN.toml" \
  evidence/inputs/cubic/CO_r1p13.toml evidence/results/CO_WINDOW.a2.toml "$REPLAY_DIR/CO_MAIN.toml"
```

The representative preflight uses `preflight` instead of `main` and the CO triple only. It scales coefficients by 7/5 and obtains a fresh exhaustive oracle; it assumes no homogeneity of the transported cubic objective. The preflight records only first calls. Inspect all method completion statuses and per-method checkpoints; a refusal is not a partial optimum certificate.

Under a separate 50 CPU / 250 wall second guard, check fresh main results with:

```sh
"$JULIA_EXE" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. independent_cubic_dp_check.jl \
  evidence/inputs/cubic/CH4_r1p09.toml evidence/results/CH4_WINDOW.a2.toml "$REPLAY_DIR/CH4_MAIN.toml" "$REPLAY_DIR/CH4_CHECK.toml" \
  evidence/inputs/cubic/CO_r1p13.toml evidence/results/CO_WINDOW.a2.toml "$REPLAY_DIR/CO_MAIN.toml" "$REPLAY_DIR/CO_CHECK.toml"
```

This checker does not rerun DP or exhaustive optimization. It independently evaluates the scalar factors and exact dyadic response vectors on fixed nonoptimal permutations, checks rational normalization and the archived optimum/witness, and verifies full-order/source/input/physical-evidence bindings. Full Taylor parity is numerical and separate from the exact frozen-table certificate. The identical full orders reuse the old physical results; no new physical win follows from reproducing them.

`evidence/cubic_dp_validation/` contains the new logs, resource receipts, results, independent review, and pre-main freeze. Its `provenance/prior_MANIFEST.sha256` preserves the old study manifest, which describes the prior source state. Changed old files are retained under `source_at_prior_commit/`; the root `MANIFEST.sha256` describes the complete current bundle.
