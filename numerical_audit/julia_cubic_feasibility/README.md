# Julia ordering feasibility and exact-window pilot

The numerical implementation in this directory is Julia, with standard-library dependencies only. It provides an exact bag-order dynamic program, an adapter for the archived nine-unit first-order response problems, and a separate cubic Strang window pilot. It does not replace QHAT's existing optimizer or change its default behavior.

See `REPORT.md` for the original pilot, and `CUBIC_DP_REPORT.md` and `CUBIC_DP_PROTOCOL.md` for the subsequent cubic DP integration validation. `evidence/` retains input/output hashes, exact results, and resource receipts. The original baseline source is QHAT `9a19a118247aeca74f968084c74533c1fd529d85`, previously recorded as the same scientific tree as merge `aca1cc33e8e2bf8c61b99775f342800e7791b988`.

## What is implemented

- `StructureDP`: complete signed local relative-order tables with exact integer/rational arithmetic; checked tree decompositions; fixed endpoints; witness gluing; admission limits before allocating factorial state tables. Float scalar costs denote their exact dyadic values.
- `ArchivedWindows`: exact expansion of the archived first-order response square into factors on at most four unit labels. Full fixed-outside interference and the objective constant are retained. Independent DFS enumeration checks all Lehmer ranks, ties, and the returned objective.
- `ResponseCore`: real equal-flip Pauli units, independent cubic Taylor columns, HF diagonal transport/projection, exact ideal-rotation affordability, and conditioned contiguous-window response DELTA tables.
- `WindowSearch`: complete bounded window enumeration with a roundoff envelope for screening and exact dyadic evaluation of every overlapping contender; same-window insertion control. The certificate covers the frozen response tables, not the exact transcendental objective.
- `CubicDP`: exact integer expansion of the complex cubic response square into signed relative-order factors on at most six labels, followed by the checked bag-order DP. It retains full baseline/outside interference and returns the exact frozen-table score. It has separate factor-construction and DP admission limits.
- `Diagnostics`: capped streaming support/overlap diagnostics. Potential support graphs, partial scans, and heuristic widths are explicitly distinguished from exact coherent factor evidence.
- `PhysicalCheck`: original-word full-space product-formula circuits, exact central/boundary fusion, stable projective state error, and number/spin leakage. Only the Hamiltonian reference uses a checked fixed-number sector; the circuit is never projected.

The cubic path currently supports a single computational-basis HF reference and real equal-flip units. The pilot uses twelve modes; physical propagation is capped at twelve modes and a reference-sector dimension of 2,000. It is not a generalized coherent-reference implementation or a global approximation-ratio solver.

## Julia environment and tests

The recorded runtime is Julia 1.10.11, one Julia/BLAS/GC thread, no startup file, and a 300 MiB heap-size hint for the historical study; post-review tests use 250 MiB. The external guard enforces a quarter core, 512 MiB RAM, no swap, memory-headroom pauses, and finite CPU/wall/output budgets. No package installation is required.

From this directory, the test command inside that guard is:

```sh
julia --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. test/runtests.jl
```

The package can also be loaded with `using CubicOrderingFeasibility` under `--project=.`. Its submodules are exported explicitly; no Python package is loaded by the Julia numerical code.

## Replaying the archived comparison

Use a fresh output directory. The scripts refuse to overwrite scientific output files. `archive_study.jl main OUTPUT INPUT...` runs each exact domain's first invocations and three warm repetitions of both solvers, alternating method order. Every completed method gets its own checkpoint before later work starts. A study can finish while a DP is explicitly refused by an admission cap; inspect each solver status.

Run the following numerical command through the retained resource guard:

```sh
julia --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
  --heap-size-hint=300M --project=. archive_study.jl main /tmp/dp-replay/RESULT.toml \
  evidence/inputs/archive/F2_r1p55.toml \
  evidence/inputs/archive/H2O_scale1p10.toml \
  evidence/inputs/archive/Li2_r2p95.toml
```

First invocations share one process and JIT helpers; they are not independent cold-process comparisons. Warm timings include factor construction, graph/decomposition, solving, and witness validation; detailed phase fields are retained. Parsing and whole-process startup are separate. Historical Python timings alone are not a controlled speed comparison. The retained fresh Python oracle invokes the old solver unchanged.

The archived target is first-order `J_D=N/H`, with exact integer numerator and the original supplied binary64 weight ratios. Neither the denominator nor a Strang `r^4` factor is changed. Only the nine-unit contiguous window varies; every outside position stays fixed. The DP's exact restricted optimum does not certify the full molecular ordering or physical error.

## Cubic and physical entry points

`cubic_dp_study.jl preflight|main INPUT PRIOR_WINDOW_RESULT OUTPUT [INPUT PRIOR_WINDOW_RESULT OUTPUT...]` compares the cubic DP with complete enumeration of exactly the same prepared window. Main mode retains first calls and three warm repetitions. `independent_cubic_dp_check.jl` separately checks factorization at nonoptimal orders and authenticates the completed certificates; see `REPRODUCE.md` for commands. The original nine-unit benchmarks concern the first-order objective and must not be read as cubic DP speedups.

`window_study.jl preflight|window|structure INPUT OUTPUT` handles one case. `window_batch.jl` handles several input/output pairs in one Julia process. Cubic window inputs contain the authenticated units and original parent order. Source-only window selection never reads physical errors.

`physical_study.jl archive MODEL_DIR ARCHIVE_RESULT OUTPUT` replays certified archived solver witnesses. `physical_study.jl cubic MODEL_DIR WINDOW_RESULT... OUTPUT` compares the frozen cubic parent/window/insertion orders and any authenticated legacy continuation control. Its returned error is `1-|overlap|`; infidelity is also retained. Norm drift is `abs(squared_norm-1)`. Reference coupling/Hermiticity checks are numerical maximum-entry diagnostics, not interval or operator-norm certificates.

The TOML input copies are direct serializations of retained JSON data. Their source-byte hashes and parent/domain bindings are recorded. Python is used only for the unchanged legacy baselines, independent legacy oracles, and the optional external resource supervisor. No new Python optimization backend is introduced.

## Current validation and result contracts

[Post-review validation](evidence/postreview_validation_2026-10-28/README.md) records the current source period separately from historical measurements. Current cubic main completion requires every scheduled method/stage and full factorial coverage on every enumeration. A first DP admission refusal has a distinct outcome and no DP certificate. Library inclusion does not run any of the six entry points.

A capped structural scan returns `execution_status=complete`, `status=diagnostic_incomplete` and `scan_complete=false`; the batch still visits all cases. These fields distinguish execution from diagnostic coverage. Old retained structure records describe the earlier contract and are unchanged. Missing optional physical continuation controls are reported in `skips`. Resource wall limits include frozen time; unknown CPU measurements remain distinct from measured zero.

The next practical and mathematical work is tracked in [FOLLOW_UP.md](FOLLOW_UP.md).
