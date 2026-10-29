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
No new results are available yet.

The same measured parent prefix is charged to each logical continuation. Both
methods charge complete worker startup, loading, search, scoring, serialization,
and teardown, including owned Julia startup/JIT. The independent archive join,
physics, reference construction, replay, reporting, and reviews have separate
resource receipts. Actual study cost counts the shared parent once; attributed
comparison cost counts it for each method/budget. Failed attempts stay visible.

The production optimizer is unchanged. Native full-order search may change any
position; the window arm keeps the existing strictly interior width-eight
domain. Four deterministic full-order starts and ample search ceilings avoid
the old 3,000-evaluation control ending prematurely. Equal CPU does not imply
equal arithmetic or identical order domains. CPU limits are ceilings: neither
method is padded when it converges or reaches a declared support limit.

Run a complete old-model preflight before generating new inputs. The standalone
drivers use the pinned research Python environment and Julia 1.10.11; numerical
execution requires the owning quarter-core, 512 MiB/no-swap resource guard.
`generate_inputs.py` freezes sources/environment and audits geometry identity
against the merged baseline. `construct.py` separates algorithm execution from
the later OS-receipt join. Freeze every parent before any continuation, and
every continuation before assembling `panels.py` and measuring physical errors.
`physical.jl`, `objective.jl`, and `verify.jl` retain independent physical and
surrogate evidence. `analyze.py` consumes the archived source/protocol, exact
order grid, and scientific-cutoff resource ledger rather than current source
defaults. All output locations must be fresh.

The primary decision uses B=64 CPU seconds, T=1.35, and the 8U rotation cap.
Both case medians must favor the same arm by more than 1e-12, with no resolved
opposite primary repetition and all verification/cost gates passing. Otherwise
the conclusion is mixed, unresolved, or incomplete. The second time and CPU
allowance are descriptive. Timing repetitions are not independent molecular
samples. No result changes the default optimizer automatically.
