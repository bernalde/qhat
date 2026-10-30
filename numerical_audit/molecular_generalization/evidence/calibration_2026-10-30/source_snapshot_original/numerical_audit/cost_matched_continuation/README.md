# Cost-matched continuation

This follow-up asks whether a fixed CPU allowance is better spent continuing
full-order cubic search or refining a supplied parent through Julia windows.
The previous prospective panel found some lower errors from sparse-parent
windows but charged substantially more construction CPU. This study compares
quality at equal total allowances; it does not establish a universal ordering.

`PROTOCOL.json` is the executable, prospective declaration. Two new geometries
from the established LiH/BeH2 families, two timing repetitions, and 32/64 CPU
second allowances produce four shared parents, sixteen continued constructors,
and eighty logical physical rows. All choices precede their physical errors.
Prospective results are not available yet. Energies are in Hartrees and times
in inverse Hartrees, with hbar=1; caps count ideal original-word rotations.

The same measured parent prefix is charged to each logical continuation. Both
methods charge complete worker startup, loading, search, scoring, serialization,
and teardown, including owned Julia startup/JIT. Deterministic restart planning
runs inside the measured shared parent process. The independent archive join,
physics, reference construction, replay, reporting, and reviews have separate
resource receipts. Actual study cost counts the shared parent once; attributed
comparison cost counts it for each method/budget. Failed attempts stay visible.

Two small optimizer repairs support the comparison: a stopped seed no longer
truncates later native restarts, and Julia avoids comparing an exact incumbent
norm to itself when an adjacent insertion is visited twice. Candidate counting
and the numerical target are unchanged. Native full-order search may change any
position; the window arm keeps the existing strictly interior width-eight
domain. Four deterministic full-order starts and ample search ceilings avoid
the old 3,000-evaluation control ending prematurely. Equal CPU does not imply
equal arithmetic or identical order domains. CPU limits are ceilings: neither
method is padded when it converges or reaches a declared support limit.

Run a complete old-model preflight before generating new inputs. The standalone
drivers use the pinned research Python environment and Julia 1.10.11; numerical
execution requires the owning quarter-core, 512 MiB/no-swap resource guard.
`generate_inputs.py --prepare --preflight PREVIOUS_DATA_ROOT --output NEW_DATA_ROOT`
first verifies and copies the representative preflight, then freezes all
producer and analysis sources/environment and audits geometry identity against
the merged baseline. Numerical source changes require a fresh preflight;
consumer-only analysis repairs can be replayed with both source versions kept.
`construct.py` separates algorithm execution from
the later OS-receipt join. Freeze every parent before any continuation, and
every continuation before assembling `panels.py` and measuring physical errors.
`physical.jl`, `objective.jl`, and `verify.jl` retain independent physical and
surrogate evidence. `analyze.py` consumes the archived source/protocol, exact
order grid, and scientific-cutoff resource ledger rather than current source
defaults. All output locations must be fresh.

Each canonical result is joined to its actual PLAN arguments, original output
path, guard allowance and receipt. Successful guarded input, panel, physical,
objective and independent jobs are also required. Execution admission is
separate from validity of a retained incumbent: genuine runtime/reply failures
remain visible but cannot support either advantage direction. Checked,
interrupted Julia checkpoints can remain eligible. A failed parent blocks its
dependent search; no new budget or substitute geometry is introduced.

The operator RUN_ROOT contains `RESOURCE_POOLS.json`, `jobs`, `logs`, a
`campaign` directory of frozen tools, and the scientific DATA_ROOT at `data`.
Disk receipts count cumulative RUN_ROOT storage, including DATA_ROOT and logs.
Terminal CPU requires process timing, a service footer or a final service
counter. A partial poll is retained as a diagnostic; unknown terminal CPU is
null and conservatively charged the full reservation. The preflight's measured
cgroup peak must fit the explicit 448 MiB admission estimate; the enforced
numerical profile stays at 512 MiB/no swap. Process RSS is retained separately.

The historical independent objective checker can abort conservatively when a
cubic response exceeds its support cap, retaining partial records. Such a
refusal does not establish objective verification or a completed comparison.
The admitted 14-mode register has 16,384 bitstrings, below that checker's
50,000-entry per-vector support cap; this generic refusal is outside this
study's admitted vector domain.
The audit ResponseCore and the budgeted production ResponseCore are separate,
intentional implementations; both source hashes are retained. Their numerical
agreement is distinct from a formal proof of either implementation.

The primary decision uses B=64 CPU seconds, T=1.35, and the 8U rotation cap.
Both case medians must favor the same arm by more than 1e-12, with no resolved
opposite primary repetition and all verification/cost gates passing. Otherwise
the conclusion is mixed, unresolved, or incomplete. The second time and CPU
allowance are descriptive. Timing repetitions are not independent molecular
samples. No result changes the default optimizer automatically.
