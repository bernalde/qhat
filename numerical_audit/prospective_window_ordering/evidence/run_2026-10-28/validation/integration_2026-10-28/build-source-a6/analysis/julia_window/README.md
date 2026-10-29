# Normal window backend: exact local target v2, protocol v1

These installed assets are the product authority for `method='window'`. They
adapt the audited `ResponseCore` and insertion portion of `WindowSearch`; their
ancestry is pinned in `SOURCE.toml`. Research inputs, results, manifests and source
snapshots under `numerical_audit` remain frozen and are not imported at runtime.
This backend has no external Julia dependencies, DP or enumeration entry point.

Before spawning, the Python owner boundedly parses installed `SOURCE.toml` and
verifies the fixed runtime files `worker.jl`, `ResponseCore.jl`, `WindowSearch.jl`
and `Project.toml` against its declared hashes. It never follows manifest paths.
The manifest is limited to8KiB/scalar metadata and each runtime file to1MiB;
this work is charged to serialization inside the same whole-call allocation.
Stats retain the manifest hash, verified runtime hashes and verification state.
Missing, malformed or mismatched assets refuse the optional stage while keeping
the checked native parent. This authenticates installed files before launch;
it is not a signature or protection against concurrent hostile file replacement.

The `exact_real_deltas_before_transport_v2` target computes local unit actions
and two/three-label atoms from supplied exact rational coefficients. It merges
coincident destinations, aggregates every pair's conditioned exterior triples,
and subtracts its baseline variant exactly BEFORE floating conversion/kernel.
HF intermediates remain present; only final factor output projects HF away.
The selector counts exact preprojection pair support, including HF, with the
same earliest-tie policy. Zero-flip/exterior units and all2/6 variants remain.
The public `atom_column` still returns Float64 columns; only its internal
selector/builder path retains rational columns through subtraction.

Exact zeros are authenticated. Each final nonzero real DELTA requires a proved
binary64 rounding cell using exact adjacent midpoints and ties-to-even.
Conversion proposes a candidate only; checking its two neighbors can repair
double rounding without changing global precision. Unsupported zero underflow,
nonfinite values, infinite-neighbor boundaries and bit exhaustion refuse.
Nonzero non-HF DELTAs erased to an all-zero complex product also refuse, even
if the approximate kernel itself is floating zero. There is no small-term
clipping or lowered local threshold. Explicit `guard=false` keeps the separate
legacy floating helper; a floating cache cannot service a guarded exact request.

The independent full-order Taylor/native baseline and all global response,
transport, score-underflow and final-acceptance guards remain mandatory. The
model is a declared floating baseline plus newly rounded exact local differences,
not an exact global-response or interval-error certificate. It differs from old
floating aggregation and historical study tables, inheriting no old optimum,
bound certificate or timing claim. Insertion uses
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

Exact centering retains the original rational target and original rotation-word
counts; its scalar conversion is unchanged by this repair. The repaired action
can differ from summing independently rounded native coefficients. It establishes
no identity with old frozen tables or the native floating aggregation. Full
Taylor/native parent compatibility and the native scorer's independent
acceptance remain required.

The existing `parent_fraction` (at most 0.60) allocates both parent CPU and
parent evaluations. For total evaluation ceiling N>=6, reserve R=2 native
scores (parent revalidation and a final candidate check) and Cmin=3 child
evaluations (child parent score, flat baseline and first trial). The parent cap
is `min(max(1,floor(f*(N-2))),N-5)`. After actual parent use P, the child can
receive `N-P-2`; an early parent stop releases unused evaluations. For N<6 the
native parent stays within N and the child is explicitly skipped. These are
allocations within the unchanged absolute ceiling, not additional evaluations.
The minimum child prefix does not promise a completed insertion sweep or an
improving checkpoint; exact tie witnesses and a checked changed-order update can
need additional evaluations. Stats label this three-evaluation prefix explicitly.
All hard refusals, checkpoint authentication and native
final acceptance still apply.

Stats retain the parent cap/use, native reserves and serialized child cap.
Reported evaluations combine completed native scores with the authenticated
child checkpoint prefix. Work after the last checkpoint can remain unknown;
its absence from that prefix is not a claim that it was free or completed.
Owned child startup/JIT/user/system CPU remains charged for the entire process.
The cooperative child CPU check uses that same process-total origin, while its
phase timings remain differences from phase starts. Child CPU/wall cutoffs that
preserve the native verification reserve retain those distinct diagnostics;
genuine whole-call exhaustion remains separate. The parent still kills/reaps its
owned child and obtains the exact wait4 receipt on every cutoff.
`stats.child_outcome` retains the child outcome even if later native verification
produces a different optional-stage status.

The reader admits n<=20, m<=1024, total words<=8192; none of those limits promises
that setup/search fits a chosen CPU allowance. Additional bounded caps cover
50000 vector/destination entries, 100000 combined local amplitude/kernel entries, 250000 DELTA entries,
1000000 flat complex entries, cumulative action/atom work and total evaluations.
Guarded local caches own only rational amplitudes; they never duplicate a float
amplitude store. They share a fixed64MiB accounted-byte allowance with reusable
complex kernel entries. Accounting uses actual GMP allocated limb capacity,
not current value bitlength, plus conservative object/dictionary/resize charges.
When either cache allowance fills, insertion stops and the value is recomputed;
inputs/terms are never thinned. Setup stats expose actual charged bytes, limit,
combined entries and skipped insertions; Python provenance records the declared
policy/caps. The unchanged checkpoint wire schema reports its usual counters.
Each exact product/addition bounds prospective intermediates to16384bits before
large integer arithmetic. Cache and arithmetic loops retain scoped deadlines
and the existing atom/action counter meanings. This is not a global RAM bound.
The reader boundary test parses 1024 labels and 8192 words on 20 bits without
running a large search. Every explicit numerical loop checks a separate scoped
budget callback. Startup/JIT and non-cooperative allocations are covered by the
parent's owned-PID watchdog. The normal worker's `--heap-size-hint=32M` guides Julia GC.
This guidance is separate from the unchanged64MiB exact-local cache allowance.
The 512MiB/no-swap owning cgroup is the hard memory bound, and the continued
prospective preflight must still satisfy its 400MiB admission gate. The hint does
not raise that gate or any numerical/resource limit; physical/independent drivers
retain their separate250M guidance. The supported production path is single-threaded,
single-request and has one worker for all windows. `with_budget` restores its
callback even on refusal; it is not a concurrent Julia library API.

Only atomic checkpoints of fully Taylor-checked global incumbents are published.
The Python owner authenticates and independently scores the final checkpoint; a
late interruption never promotes a partial candidate. A worker exit code of two
may leave a valid earlier checkpoint. No automatic download, environment
instantiation, package installation or saved Julia configuration change occurs.

Top-level `TrainingResult.status` describes the usable native incumbent under
this optional-stage contract; it is not an assertion that every requested search
phase finished. A completed native parent can retain `status='complete'` when
the optional backend fails. A completed optional stage may also finish after a
parent sub-search reached its allocated cap. Inspect `stats.parent_status` and
`stats.window_status` separately for those phase outcomes. Both remain explicit;
backend failure never promotes an unverified changed order.

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

For a separately guarded old/current scalar-action regression, the same test
file accepts `--exact-action-probe PATH_TO_RESPONSECORE`. It runs only three
tiny cancellation/zero checks against that read-only source and exits; old
conservative refusal is expected to fail the nonzero checks. The normal suite
also covers exact midpoint ties and either side beyond 256 bits, subnormal/range
refusals, centered HF/other-state actions, bit/budget interruption and unchanged
global projection/kernel/Taylor/score-underflow guards. V2 tests independently
construct literal rational Pauli matrices and dense Taylor polynomials for all
tiny orders, including coincident flips, large exterior cancellation before
conversion, suffix-exterior interference, HF intermediates, original costs,
allocated-capacity cache admission and local transport refusal. Test authorship alone is
not a recorded pass or a successful molecular preflight.
