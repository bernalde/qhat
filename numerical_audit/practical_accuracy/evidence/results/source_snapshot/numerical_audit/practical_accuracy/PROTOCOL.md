# Practical target accuracy under total classical construction caps

This is a new controlled development study for paper issue #16. It reuses the
three authenticated molecular Hamiltonians from the completed studies, at T=2.11.
It makes no prospective molecular generalization claim; issue #17 owns that test.
The orbital basis, scalar removal, exact supplied coefficients and HF determinant
remain unchanged. No physical error is used to initialize or score an order.

## Frozen panel

Models: Li2_r2p95 (6 modes,2 electrons), F2_r1p55 (8 modes,6 electrons),
H2O_scale1p10 (8 modes,4 electrons). Symmetric formula only. Evaluate every r in
1,2,3,4,6,8,12,16,24,32,40,48,64,80,96,128. State-error targets are
0.0001,0.000001,0.00000001 in the inherited phase-insensitive state-distance
convention `1-|<exact|actual>|`, evaluated using a residual to avoid cancellation.
Absolute threshold ambiguity is 2e-12; ambiguous candidates do not qualify.
Errors at or below2e-12 receive an informational numerical-floor label; no such
case is removed or promoted to exact-zero evidence.

Four new methods use matched standalone total CPU caps of2 and8 seconds, three
repetitions each: conserving first-order JD, individual-Pauli first-order JD,
conserving symmetric cubic JD, and conserving cubic JD with a charged first-order
warm start. The first three use current placement, its reversal and canonical
order as three starts (Pauli uses the compatible flattened current placement).
The fourth spends45% of remaining nominal ordering time on first-order training,
then uses its result, its reversal and current placement as cubic starts. Every
method recomputes placement/model/method setup inside its own caller timer.

The nominal ordering allowance is half the total cap. Context setup, seed scores,
search, independently rescored accepted moves and cleanup all consume that half.
The established three-start insertion policy uses equal remaining CPU slices,
no unused carry, at most3000 candidate evaluations and5 moves per start. It keeps
completed scores and partial sweeps; a lack of seed/search time returns the valid
initializer and explicit incomplete status. A20ms nominal cleanup reservation
is not a claim that an operation cannot overrun. Actual elapsed work is decisive.
All scores use the inherited30U ideal-rotation normalization at this same T;
target errors never select the training score or a winning repetition.

Retain three frozen strong references per model: prior current-placement,
larger-search Pauli symmetric, and cubic-insertion symmetric orders at30U,T2.11,
authenticated against actual ordering merge8fd424f532cce04b595abc3445cc032ab03ab6c4.
Their original standalone ordering costs are inherited charges, explicitly outside
the matched-cap cohort. They are neither rerun nor advertised as matched CPU.
Main totals:72 matched curves plus9 references,1296 grid candidates,243 target rows
and918 paired decisions, including direct charged-warm-start comparisons and
comparisons against every inherited reference (explicitly outside matched CPU).
The method/cap/repetition dispatch is cyclically rotated by case and repetition;
one symmetric unmeasured ordering-plus-compilation warmup for every method
precedes measurements. Its compiler usesr3 and retains an exact trace, with no
validation error measured. The same rotated dispatch controls ordering and
measured compilation, rather than sorting method names for the latter.

## Construction and screening accounting

A standalone candidate charge is measured ordering CPU through its full return,
including input read/hash/validation, base/placement, method setup, initialization,
search, cleanup and layout, plus compilation CPU through its full return. The
latter includes a new model read/hash/validation, HF preparation, exact boundary
fusion, lowering/routing/cancellation, resource counting, gzip trace serialization,
hashing and cleanup. All CPU overruns are retained. A complete circuit qualifies
for a matched total cap only if this sum is at most its declared cap. Lowering
checks the remaining CPU and retains an exact partial trace on interruption;
counting/serialization can overrun and then make a completed trace ineligible.
All complete circuits, including ineligible ones, are physically screened.

Each measured order is reused across its16 grid points. Repeated standalone order
charges describe hypothetical one-circuit construction, not actual campaign CPU.
Report actual campaign work (including input/source authentication, warmups,
metadata, physical references, full-grid screening, verification and reporting)
separately. Screening is an offline comparison, not a deployable oracle for choosing
r. Correlated preparation is absent here because the state is HF; broader-state
preparation costs belong to #17. These caps do not include interpreter startup,
module import/JIT warmup or development; their actual costs remain in process logs.

## Compiler and physical checks

Reuse source-pinned line-pauli-exact-angle-0.1.0: line0..n-1, exact H/S/Sdg/Rz/CX,
nearest-neighbor SWAP-based routing restored after each remote CX, only adjacent
inverse cancellation and adjacent same-wire Rz combination, exact unit-boundary
fusion, inherited12M generated/2M stored gate guards. No hardware, synthesis-
optimality or separately depth-optimal claim. Count all HF preparation gates.

Freeze all orders and all complete/partial exact primitive traces before any
new validation error. Full-space extended-precision primitive replay preserves
leakage; the reference exponential uses the full target. Independently check
literal Kronecker Pauli actions, unfused ideal formulas, extended-precision
primitive replay, two-pass circuit DAG costs, state errors/leakage, all target
choices and pairwise decisions. Target choice minimizes CX, then unit-duration
ASAP depth, total gates, then r among complete, CPU-eligible, strictly qualifying
candidates. Preserve unmet, ambiguous, ineligible, censored, floor and all loss/tie
cases. Repetition-level comparisons precede any aggregation; no best repetition
is selected using validation errors. Raw errors and both classical/circuit costs
remain available for every complete candidate.

Source snapshots authenticate every executable dependency; current verifier
identity stays distinct from producer identity. Final extraction/replay and Git
blob checks cover every delivered file including NPZ binaries and paper omissions.
No old validation set is retuned. Any later source repair preserves the producing
snapshot and cumulative failed-attempt usage rather than overwriting history.

## Preflight and finite resource ceilings

Before main runs, use H2O with every coefficient scaled by7/5, renamed preflight,
both caps, one repetition, all four methods and all three frozen reference orders.
The references test pipeline mechanics on altered coefficients, not new inherited
cost validity. This full-path preflight has11 curves,176 candidates,33 target rows;
it includes102 paired decisions, independent checking, packaging, extraction
and full byte authentication against the independently replayed circuits. It may inform
implementation/resource feasibility, never tune the main panel on its errors.

Finite cumulative phase ceilings, frozen before execution (CPU seconds / wall
seconds / output MiB): tests120/900/64; preflight600/3600/256;
measurement3600/21600/512; independent verification3600/21600/512;
reporting300/1800/512. Failed attempts remain charged. Every local worker uses
quarter-core CPU,Nice19,one thread,MemoryHigh448MiB,MemoryMax512MiB,no swap;
start/resume above5.5GiB available,freeze only owned service below5GiB.
No numerical/reviewer overlap and no unapproved resource increase.

Prospective holds, sharper epsilon/finite-step proofs and actual manuscript
integration remain separate required follow-ups, not completed by this panel.
