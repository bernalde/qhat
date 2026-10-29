# Normal window backend v1

These installed assets are the product authority for `method='window'`. They
adapt the audited `ResponseCore` and insertion portion of `WindowSearch`; their
ancestry is pinned in `SOURCE.toml`. Research inputs, results, manifests and source
snapshots under `numerical_audit` remain frozen and are not imported at runtime.
This backend has no external Julia dependencies, DP or enumeration entry point.

The new target aggregates each pair's real column with its conditioned exterior
triple columns before transport. It retains the independent full-order Taylor
baseline and every coherent DELTA, including exterior interference. A nonzero
local unit action/atom/aggregate/DELTA below 2^-24 of its operand norm scale is
refused; computed zero atoms/DELTAs are admitted without claiming exactness.
Near-zero unit actions are additionally checked with their supplied exact
coefficients, admitting certified scalar/forbidden-action zeros. This is a
cancellation heuristic, not an
interval-error certificate. The resulting binary64 tables differ from historical
study tables and inherit no historical timing or optimum claims. Insertion uses
an outward roundoff envelope and exact dyadic vector rescoring when intervals
overlap. Every global proposal needs a resolved full Taylor decrease, followed by
the unchanged native scorer's independent acceptance.

The request/reply protocol is bounded TOML v1 (8 MiB request, 64 KiB reply). Exact input coefficient fractions
and the native float's exact dyadic time spelling are serialized as numerator /
denominator strings; each rational input is capped at 4096 bits. Exact centering
and diagonal phases use those coefficients with a 16384-bit intermediate cap;
exact supplied-table witnesses have an 8192-bit cap. Labels cross the Python
zero-based / Julia one-based boundary once. Full orders, reference occupation,
original word counts, time, pulse cap, parent score and work/storage limits are
bound by SHA-256 over the exact request bytes. Replies include original lengths,
a complete order, selected interior intervals, replayable before/proposal labels
and acceptance for each completed window, Q, J=Q/r^4 and r. The Python owner
replays that history from the authenticated parent, verifying each fixed
exterior and original global label transition.

The reader admits n<=20, m<=1024, total words<=8192; none of those limits promises
that setup/search fits a chosen CPU allowance. Additional bounded caps cover
50000 vector/destination entries, 100000 cached amplitudes, 250000 DELTA entries,
1000000 flat complex entries, cumulative action/atom work and total evaluations.
The reader boundary test parses 1024 labels and 8192 words on 20 bits without
running a large search. Every explicit numerical loop checks a separate scoped
budget callback. Startup/JIT and non-cooperative allocations are covered by the
parent's owned-PID watchdog. The supported production path is single-threaded,
single-request and has one worker for all windows. `with_budget` restores its
callback even on refusal; it is not a concurrent Julia library API.

Only atomic checkpoints of fully Taylor-checked global incumbents are published.
The Python owner authenticates and independently scores the final checkpoint; a
late interruption never promotes a partial candidate. A worker exit code of two
may leave a valid earlier checkpoint. No automatic download, environment
instantiation, package installation or saved Julia configuration change occurs.

Root-run checks (use the task's existing external CPU/memory guard):

```
julia --startup-file=no --history-file=no --threads=1 --gcthreads=1 --project=analysis/julia_window analysis/julia_window/tests.jl
python -m pytest analysis/tests/test_window_ordering.py analysis/tests/test_sparse_response_ordering.py analysis/tests/test_train_ordering_cli.py -q
```

Set `QHAT_JULIA_EXECUTABLE` to an existing Julia 1.10+ binary. Set
`QHAT_WINDOW_WHEEL` to a wheel built from the exact source to enable the installed
wheel subprocess/CLI check; the focused CI lane does both. The Julia tests use
only tiny coupled fixtures and a reader-only cap boundary. They make no molecular
performance or physical-accuracy claim.
