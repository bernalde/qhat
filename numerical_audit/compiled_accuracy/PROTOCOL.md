# Compiled resources at fixed target accuracy

Issue #9 extends ideal Pauli rotations to a declared reproducible compiled-gate
model. This protocol is fixed before this study's errors. Implementation and
preflight must complete before any full-panel result is claimed. Numerical
checks do not prove interval error bounds, production-compiler optimality or
hardware performance. No device is contacted.

Use the three supplied rational models and HF determinants from issue #8,
with its fixed orbital/qubit basis and both times0.83,2.11 inverse Hartree.
Before observing any issue #8 validation error, this study precommitted all
eight of its methods at each time/model, selecting the fixed30U training cell.
No method, model, time or loss is selected using its observed accuracy. Inherit
the actual source/result hashes and construction charges, labelling development
dependencies explicitly until their actual merged source is verified.

Eight method/formula pairs per model/time are current placement(first/symmetric),
fresh same-time conserving J_D(first/symmetric), larger same-time individual-
Pauli J_D(first/symmetric), new cubic orientation(symmetric) and cubic insertion
(symmetric). These are orders trained at a30U ideal-rotation cap, not orders
optimized for compiled cost or newly optimized for every step count.

For each of48curves, predeclare r in1,2,4,8,16,32,64,128:384candidate circuits.
Freeze every feasible candidate's exact primitive trace, source/model/order/time,
step count, cost and SHA256 before any physical evaluation. No error-guided
binary search or presumed error monotonicity is used. Preserve censored compiler
candidates, unmet accuracy targets and all losses without replacement.

Compiler: repo-native prototype `line-pauli-exact-angle-0.1.0`, pinned by the
committed source hash before preflight and main. Native gates H,S,Sdg,Rz(angle),
CX, with undirected nearest-neighbor connectivity on the original register line
0..n-1. No ancillas, qubit remapping search, device calibration or noisy gates.

For each nonidentity Pauli word exp(-i theta P), sort support by qubit. X uses
H before/after. Y uses chronological Sdg,H before and H,S after. A support
parity ladder computes into the last support qubit, applies Rz(2theta), then
uncomputes parity and basis. Rz(phi)=diag(exp(-i phi/2),exp(i phi/2)).
Scalar terms must already be removed according to the supplied target policy;
the compiler rejects empty words rather than silently changing global phase.

A remote CX swaps its control along the line until adjacent to the target,
applies CX, then reverses every SWAP. Each SWAP is three adjacent CX gates.
The logical mapping is restored after each remote CX. A distance-d remote CX
therefore lowers to6(d-1)+1CX before cancellations. This is a fixed routing
policy, not an optimal route or production backend.

Apply exact center/boundary fusion of internally commuting complete units before
Pauli lowering, using the established circuit convention. Then remove only
globally adjacent inverse Clifford pairs or identical CX pairs, and combine
globally adjacent Rz on the same qubit. There are no commuting rewrites across
intervening gates and no approximate synthesis or error-based cancellation.
Angles use supplied rational coefficients, exact binary64 T.as_integer_ratio
and rational durations1/r or1/(2r); combine Rz angles as exact Fractions.
Convert to floating trigonometry only for numerical replay. This does not make
the real transcendental unitary an interval-certified object.

Include HF preparation from all-zero explicitly: H,S,S,H per occupied qubit is
exact X including global phase. The unoptimized common preparation has4Ne
single-qubit gates,zeroCX and4unit-duration layers; final reported costs include
any allowed adjacent cancellation at its circuit boundary. This computational
determinant preparation is not correlated-state preparation or a hardware cost.

Retain every candidate's entire deterministic gzip NDJSON primitive trace with
exact rational Rz numerators/denominators and target identity. Count total gates,
single-qubit gates,CX and ASAP depth with unit-duration gates. A second declared
CX-dependency depth treats single-qubit gates as zero-duration. Neither is a
calibrated execution time. Routing's CX gates are included. Ideal Pauli rotations
remain a separate column, not a substitute for compiled counts/depth.

The fixed real-target producer rejects odd-Y Hamiltonian terms. Exact unit
fusion requires internally commuting Pauli words: the producer checks integer
anticommutation parity and the independent reference checks literal matrix
commutators. A failed precondition rejects the target before applying fusion.

Declare per-candidate hard bounds before compilation: at most8qubits,128Pauli
terms/units,12million generated primitive insertions and2million stored final
gates. A candidate exceeding a bound is censored with its reason and retained
partial work; it cannot be selected. Record the largest feasible r on the fixed
grid for every method, without claiming a continuous-r resource optimum.

Accuracy targets are fixed projective state errors0.01 and0.001. Use the same
stable full-space normalized metric1-|overlap| as the earlier studies. Keep
actual norms,infidelities,number and joint(N_alpha,N_beta)leakage and all raw
grid errors. No sector projection or hardware noise model is applied.
For each target, select the feasible qualifying grid circuit with lowest CX,
then lowest total ASAP depth,then fewest total gates,then smallest r. A point
within2e-12 of the threshold is numerically ambiguous and is not certified as
attaining it. A curve without a qualifying point is explicitly unmet/censored.
No favorable case replacement or unperformed target comparison is allowed.

Independent verification checks literal Pauli/basis/routing unitaries, exact
angle and inverse rewrites, and the complete compiled trace's target binding.
Replay primitive circuits from all-zero and compare against separately computed
unfused ideal formulas within2e-10 in state and2e-12 in physical metrics. A separate
scaled-Taylor exact-H reference checks SciPy propagation. Recompute all resource
counts/depth, feasibility,censoring and target-selection decisions, retaining
every disagreement/ambiguous threshold. Complete unitary tests use small
registers and hand-derived phases, not only an implementation-mirroring counter.

Keep inherited ordering costs separate from newly measured circuit construction,
lowering/routing/cancellation,serialization/hashing and physical replay CPU/wall.
Measure complete caller returns and campaign process resources; do not omit
preparation/common work or sum repeated standalone charges as fresh campaign
work. A prototype result implies no practical classical speedup or optimal
hardware/compiler choice.

Finite cumulative ceilings before any tests/preflight:

| Phase | CPU seconds | Wall seconds | Output MiB |
| --- | ---: | ---: | ---: |
| Focused unitary/compiler tests | 120 | 900 | 64 |
| Largest altered-target complete preflight | 300 | 1800 | 128 |
| Full compile/freeze/physical panel | 1800 | 10800 | 256 |
| Independent trace/resource/physical verification | 1800 | 10800 | 256 |
| Report/archive/extraction/compact checks | 240 | 1800 | 128 |

Preflight uses the largest altered H2O7/5 input atT2.11, all8fixed comparator
orders from the authenticated issue8altered preflight's30Utraining cell, and
the complete8point r-grid including128. Thus64candidate traces, all compilation,
resource/depth/accuracy/selection checks,independent replay,report,archive/
extraction and strict delivery are exercised before main. Measure CPU/wall,
memory,retained/temporary/trace/archive disk and replay,then commit margin-bearing
main estimates before any main errors. Failed/interrupted attempts remain
charged to their original phase. No silent enlargement or easier replacement.

One local numerical/build/test worker and one BLAS/Numba thread, with no
overlapping reviewer. CPUQuota25%,Nice19,CPUWeight5,IOWeight10,MemoryHigh448M,
MemoryMax512M,MemorySwapMax0. Start/resume above5.5GiB available; below5GiB freeze
only this task's service. Existing pinned Python3.12/NumPy2.5.3/SciPy1.18.1,
no new dependency installs. Every heavy command uses guarded.py. Required
resource decisions preserve partial evidence and do not reset counters.

Conclusions describe this fixed prototype,connectivity,gate model,grid and small
related model family. They cannot establish optimal compilation,routing,
continuous-r target cost,asymptotic scaling,conservation as an accuracy mechanism,
large-molecule validity,noise tolerance or actual hardware execution.


## Retained preflight arithmetic correction

The first full64candidate preflight compiled and propagated all circuits, then
failed the fixed2e-12 independent physical-metric tolerance at the r64 cubic
insertion circuit. A retained diagnostic isolated raw norm drift: literal versus
unfused norm difference3.1532554345403696e-12, while projective-error differences
were below1e-20. This is a numerical replay issue, not an accuracy-tolerance or
compiler-circuit change. Both primitive replay implementations now accumulate
in NumPy complex long double, with long-double Hadamard normalization and
rational-angle trigonometry. The independent unfused formula uses long-double
coefficients/trigonometry too; the SciPy and scaled-Taylor exact-H references
remain distinct. No gatewise renormalization is used. Record actual available
long-double precision in execution evidence; this arithmetic is still numerical,
not interval-certified. Preserve the failed run and diagnostic, charge both to
the original preflight pool, and rerun the complete64candidate pipeline before
freezing main estimates. All grids, targets,tolerances,compiler/version,exact
circuits,routing,rewrites,selection rules and resource ceilings remain unchanged.


## Calibrated main estimates frozen after the complete preflight

`calibration/PREFLIGHT.json` records the passed64candidate/16decision altered
largest-target end-to-end run at producing source ae08e0a. All64circuits were
feasible; independent replay,counts/depth,target decisions,report,archive,
extraction andstrictchecks passed. Bothfailedpreflight anddiagnostic remain
charged to the original300CPU/1800wall/128MiB phase,with all resource records.

Measured preparation/compilation38.7952CPU seconds,physical propagation42.9614,
andindependent verification118.4774. Freeze main estimates1000CPU/4500wall/
80MiB output,and independent1500CPU/6500wall/110MiB. These round upward from
six largest case/time cells times a further2xCPU/wallmargin. Estimate256MiB
peak memory from the measured108.1MiB cgroup peak,below the512MiB hard cap.
Freeze reporting/archive/extraction estimate30CPU/240wall/112MiB; itsdiskbasis
is six times largest delivery+extracted delivery+archive+report with an extra
1.5xmargin. Mainresults remain in their original phase directory rather than
beingduplicated into reporting. Exactbasisand measurements are in
`calibration/MAIN_ESTIMATES.json`; original main1800/10800/256,independent
1800/10800/256,andreport240/1800/128 ceilings are unchanged. No maincompiled
candidate orphysical error has been generated at this calibration checkpoint.
