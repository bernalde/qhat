# Julia ordering feasibility and exact-window pilot

The numerical implementation in this directory is Julia, with standard-library dependencies only. It provides an exact bag-order dynamic program, an adapter for the archived nine-unit first-order response problems, and a separate cubic Strang window pilot. It does not replace QHAT's existing optimizer or change its default behavior.

See `REPORT.md` for measured outcomes, `PROTOCOL.md` for the frozen study, and `evidence/` for input/output hashes, exact results, and all retained resource receipts. The baseline source is QHAT `9a19a118247aeca74f968084c74533c1fd529d85`, previously recorded as the same scientific tree as merge `aca1cc33e8e2bf8c61b99775f342800e7791b988`.

## What is implemented

- `StructureDP`: complete signed local relative-order tables with exact integer/rational arithmetic; checked tree decompositions; fixed endpoints; witness gluing; admission limits before allocating factorial state tables. Float scalar costs denote their exact dyadic values.
- `ArchivedWindows`: exact expansion of the archived first-order response square into factors on at most four unit labels. Full fixed-outside interference and the objective constant are retained. Independent DFS enumeration checks all Lehmer ranks, ties, and the returned objective.
- `ResponseCore`: real equal-flip Pauli units, independent cubic Taylor columns, HF diagonal transport/projection, exact ideal-rotation affordability, and conditioned contiguous-window response DELTA tables.
- `WindowSearch`: complete bounded window enumeration with a roundoff envelope for screening and exact dyadic evaluation of every overlapping contender; same-window insertion control. The certificate covers the frozen response tables, not the exact transcendental objective.
- `Diagnostics`: capped streaming support/overlap diagnostics. Potential support graphs, partial scans, and heuristic widths are explicitly distinguished from exact coherent factor evidence.
- `PhysicalCheck`: original-word full-space product-formula circuits, exact central/boundary fusion, stable projective state error, and number/spin leakage. Only the Hamiltonian reference uses a checked fixed-number sector; the circuit is never projected.

The cubic path currently supports a single computational-basis HF reference and real equal-flip units. The pilot uses twelve modes; physical propagation is capped at twelve modes and a reference-sector dimension of 2,000. It is not a generalized coherent-reference implementation or a global approximation-ratio solver.

## Julia environment and tests

The recorded runtime is Julia 1.10.11, one Julia/BLAS/GC thread, no startup file, and a 300 MiB heap-size hint. The external guard enforces a quarter core, 512 MiB RAM, no swap, memory-headroom pauses, and finite CPU/wall/output budgets. No package installation is required.

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

`window_study.jl preflight|window|structure INPUT OUTPUT` handles one case. `window_batch.jl` handles several input/output pairs in one Julia process. Cubic window inputs contain the authenticated units and original parent order. Source-only window selection never reads physical errors.

`physical_study.jl archive MODEL_DIR ARCHIVE_RESULT OUTPUT` replays certified archived solver witnesses. `physical_study.jl cubic MODEL_DIR WINDOW_RESULT... OUTPUT` compares the frozen cubic parent/window/insertion orders and any authenticated legacy continuation control. Its returned error is `1-|overlap|`; infidelity is also retained. Norm drift is `abs(squared_norm-1)`. Reference coupling/Hermiticity checks are numerical maximum-entry diagnostics, not interval or operator-norm certificates.

The TOML input copies are direct serializations of retained JSON data. Their source-byte hashes and parent/domain bindings are recorded. Python is used only for the unchanged legacy baselines, independent legacy oracles, and the optional external resource supervisor. No new Python optimization backend is introduced.
