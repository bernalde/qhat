# Sparse molecular ordering: frozen larger-case comparison

This study evaluates a bounded ordering algorithm; it does not select molecules,
reference states, times, budgets or repetitions using validation errors. The
registry is fixed before chemistry generation. All four inputs have twelve spin
modes and six active electrons: two N2 geometries, CO and tetrahedral CH4 in
STO-3G, with the stated contiguous occupied/virtual active windows. They are four
cases from three families; the two N2 geometries are dependent observations.
Geometry identity is checked against the numerical-audit JSON at the pinned base.
Unseen means absent from that recorded baseline, not absent from external data.

RHF uses one thread, 256 MiB internal memory, tolerance 1e-11 and at most 100
iterations. A failure is retained and stops the complete-panel claim; cases and
tolerances are not replaced. The producer retains orbitals, active tensors,
exact dyadic Pauli coefficients, the adjoint-average/whole-parent cutoff policy,
removed scalar, HF contraction and orthogonality checks. Complete equal-flip
units must commute internally and conserve particle number exactly.

## Orders and resource boundaries

Train on the HF occupation at time 1.0 with the 30U ideal-rotation budget. U is
the number of retained Pauli words. All six arms use the symmetric formula:
canonical units, coefficient-magnitude units, current-placement units,
first-response insertion, individual-Pauli insertion, and sparse cubic insertion.
The three searched arms receive equal 32-CPU-second deadlines, three declared starts,
up to 3000 candidate evaluations and five accepted moves per start. Actual setup,
search and return costs are recorded, including incomplete scores and refusals.
The cheap unsought baselines keep their actual lower charges. No validation error
or exact-H time propagation enters any constructor. Sorting, setup, seed scoring,
search and cleanup belong to each standalone ordering charge.

Training-resource amendment before any main orders or errors: the initial
two-second altered-input preflights could not finish one cubic seed on the
largest input. A separate, resource-only probe timed one fixed current-order
seed on that same 7/5-scaled CH4 Hamiltonian. Before the probe, the rule was fixed
as the smallest power of two at least
`4*(context CPU + first-score CPU + 0.7 initializer CPU)` seconds, minimum two
and maximum 32; a requirement above 32 would be a recorded resource refusal.
Measured context and score CPU were 0.005947819 and 7.031157889 seconds, selecting
32 seconds for all three searched arms. No physical-error fields, other seed
qualities or main outcomes informed this allowance. Actual incomplete work and
overruns remain charged. The standalone library's configurable two-second
default remains unchanged; this amendment belongs to this benchmark.

Pre-main review correction: the CH4 singleton-Pauli support was counted without
pair allocation or physical errors: 588,360 qualifying non-HF commutator pairs
among 1,481,781 unordered pairs. The constructor stores one bounded dictionary
with at most 750,000 such entries, without a duplicate pair list. Its actual
setup/search must complete a scored seed under the unchanged 512 MiB OS cap
before the main panel is described as containing a searched CH4 Pauli control.
Earlier 250,000-entry/two-second refusals remain retained; no accuracy outcome
informed the storage correction.

Freeze all 24 orders before evaluating any validation error. Reuse each order at
validation times 0.5 and 1.5 and ideal-rotation caps 10U and 30U. The affordable
step count is obtained from exact fused symmetric rotation accounting, not from
validation error. The declared panel has 96 circuits and eight reference
evolutions. Reference and validation work are offline campaign costs, never an
advertised deployable accuracy oracle. A separate conservative step selector may
provide a valid bound; the response score itself is not a finite-error bound.

## Compiled circuits and independent validation

Every candidate is actually lowered to exact rational H/S/Sdg/Rz/CX primitives
on a nearest-neighbor line, including HF preparation. Adjacent inverse/rotation
cleanup and equal-unit boundary fusion are the only optimizations. The compiler
is a separately versioned twelve-mode adaptation of the earlier eight-mode
prototype; frozen historical modules remain unchanged. Retain complete traces,
generated/stored gate counts, CX count, and ASAP/dependency depth. Censors and
overruns are retained. No native hardware or synthesis-optimality claim follows.

The independent validator constructs literal full-space Pauli actions and a CSR
Hamiltonian, and computes reference actions independently of the trainer. Its
4096-state exponential work is explicitly limited to twelve modes. It verifies
each compiled primitive trace using an exact inverse-Clifford Pauli frame,
numerically applying each actual Rz as its conjugated Pauli rotation. The final
frame must equal the expected HF bit flip up to global phase. Small tests also
compare this replay against literal per-gate matrices. Thus all primitive gates
are inspected without allocating a full state for every Clifford operation.

Independent ideal rotations, primitive replay, two-pass circuit DAG costs and
state metrics must agree within fixed numerical tolerances. The state metric is
1-|overlap|, evaluated through a residual; squared infidelity is recorded
separately. Keep norm drift and particle/spin-sector leakage. Comparisons are
paired within case/time/cap before aggregation. Report mixed outcomes, setup and
total costs, and numerical-floor rows. These small-space validation arrays are
not part of the sparse training algorithm or a proof of its large-case speed.

Every final frozen order also receives a separate literal-matrix objective audit:
nested commutators for the symmetric cubic defect and a commutator column for
the first-response baselines, with diagonal transport evaluated independently.
These checks never call the new trainer's score. A stored completed training
score must match this oracle within relative 1e-7 and absolute 1e-20; these are
numerical parity checks, not interval certificates. Fall-back orders with no
completed score retain that status and receive a separately labelled oracle.

Executable precision amendment before the main order freeze: primitive/ideal
circuit states and sparse-matrix input vectors use complex128. The separate
scaled Taylor reference still uses extended precision. Exact rational gate
angles, fixed parity/norm/leakage tolerances and the 1e-12 paired-comparison
ambiguity are unchanged. The earlier extended-precision preflight and its CPU
refusal are retained. This change responds to measured resource cost, without
using preflight accuracy to choose the main panel or ordering objective.

## Execution policy

Root schedules every heavy command under the task's shared phase allowances:
one worker, quarter-core CPU, 512 MiB RAM, no swap, and one native thread. Start
above 5.5 GiB available memory and freeze the owned worker below 5 GiB. Source,
environment, protocol, raw attempts and failures are retained. The task ledger,
not this document, owns actual finite CPU/wall allocations. Preflight may reveal
implementation/resource problems but cannot change the frozen main panel using
its accuracy outcomes. Any protocol amendment receives an explicit version and
preserves the original. No study-producer run is implied by source preparation.

Representative preflight selects the generated model with the largest retained
Pauli-word count, with a lexical case tie. It scales all supplied coefficients
by 7/5 under a separate preflight identity, freezes six orders, then checks six
circuits at time 1.5 and cap 30U against one reference. This rule uses resource
size only; main-case errors stay unexposed. Preflight raw accuracy values cannot
be used to select or retune the main scientific panel.
