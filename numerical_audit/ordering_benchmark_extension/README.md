# Expanded conserving-ordering benchmark

This local extension completes the comparison panel on the same fourteen molecular targets and states used by the conserving-response research. It is separate from the earlier immutable `ordering_mechanism` artifact. Nothing is published remotely.

The completed design contains 4,312 rows: 4,186 saved evolved states and 126 infeasible budgets. Each case has 36 deterministic schedules and three families of twenty fixed random permutations; randoms are evaluated at budget 100U, while deterministic controls use their declared finite work grid. All schedules are frozen before their campaign's first newly measured error. Cases are correlated development/follow-up fixtures, not fourteen independent molecules.

Start with `RESULTS.md`, `TIMINGS.md` and `VALIDATION.md`. The central empirical distinction is between a conserving constrained search and a larger individual-Pauli search; matching examined positions changes that comparison. Symmetric controls and a common boundary-fusion sensitivity prevent a first-order improvement from becoming an unsupported universal cost claim. Prescriptions adapted from the literature are labeled as adaptations.

`extension_validation.zip` contains the exact source/input snapshot, frozen schedules, all rows and states, independent numerical checks, native V2 checkpoint histories and publication analysis. The 264 MB generated ZIP is retained locally outside Git history; its code,
record and verification evidence are tracked. A fresh clone therefore needs
the separate archive. Its `extension_validation.record.json` pins the SHA-256, size, inventory and hash preflight. `verify_archive.py` is valid with Python -O. It checks bytes, not physics.

`source/` supplies visible source and protocol copies for review; `source_manifest.json` binds them. For execution, extract the archive into a new directory retaining the package name `qhat`, and follow `REPRODUCE.md`. Do not run the visible source directory as though it contained the raw data, and do not import an unrelated live checkout. The archive's `qhat/` directory is the exact source/input subset used by this experiment, not an installation of every optional QHAT feature.

The fresh pinned-environment replay evidence is in `clean_replay/`. It checks archived import paths, regenerates all fourteen target references and conserving partitions, independently scores 56 final objectives, propagates five selected production circuits, and reconstructs the full analysis/tables. This complements the original audit of all saved metrics and 3,258 independently replayed circuits; it does not claim a second clean-environment propagation of all 4,186 states.

Runtime fields distinguish setup, search and propagation. Pair setup is shared, raw timer boundaries differ, and concurrent jobs prevent isolated wall-time speedup claims. The exact NH3 incremental-V2 search took 4,952 seconds; its positively scaled companion reuses the verified order. This reuse is recorded rather than counted as a fresh fast search.

The accompanying paper revision is in the adjacent paper research worktree, on `research/ordering-mechanism-2026-10-09`. Its `ORDERING_BENCHMARK_EXTENSION.md` explains the new mathematical exposition, prior-art audit and limited Lean bridge scope.
