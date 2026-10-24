# New exact nine-unit response windows

Issue #6 asks how far a heuristic is from an optimum of a fixed surrogate.
This protocol is frozen before any new window enumeration or physical errors.
The completed four 8! windows are prior evidence and will not be rerun.

Inputs: the three authenticated prospective-reference molecular models
Li2_r2p95 (6 modes/2 electrons), F2_r1p55 (8/6), and H2O_scale1p10 (8/4),
with their supplied integer coefficient numerators and positive common
denominators. The fixed orbital basis, scalar-removal policy and complete
conserving units are inherited without regeneration. The initial occupation is
the HF determinant `(1<<electrons)-1`, T=1. No failed two-level correction is
retuned or admitted as an accepted method.

For each model, start from its authenticated prior HF `optimized_JD` order.
Consider every contiguous nine-unit window in that fixed order. Select the
window having the largest number of nonzero internal commutator-response pairs,
breaking ties at the earliest offset. This uses the Hamiltonian and determinant
only, with no new propagation or physical error. Freeze the selected nine labels,
the full order outside that contiguous window, all coefficients, all pair
amplitudes and supplied weights in `DOMAINS_FROZEN.json` before enumeration.
Exactly 9! = 362880 permutations are feasible for each of three windows. There
are no whole-order, noncontiguous-window, orbital or state optimizations.

Let `b_w(order)` be the integer coefficient of destination w in
`sum_{i<j in listed order}[G_i,G_j]|HF>`, using numerator units. Remove its HF
parallel component. The fixed objective is

`J_D(order) = sum_w rho_w b_w(order)^2 / (4 denominator^4)`.

The diagonal D consists of all diagonal Pauli terms. At T=1, supplied
`rho_w=float(np.sinc(gap_numerator/(2*pi*denominator))**2)`. Each binary64
weight is converted via `as_integer_ratio` to an exact dyadic rational, and a
common positive scale turns comparison into arbitrary-precision integer
arithmetic. Exact integer coefficient numerators are never multiplied through
floating complex values. The certificate concerns these *supplied* weights,
not exact transcendental sinc weights; independent math.sin checks are
floating sanity checks, not interval proofs. These supplied weights can differ
in the last bit from the prior heuristic's floating score implementation.
Its authenticated order is a feasible warm start; its old floating value does
not define this new fixed representation.

The solver enumerates all permutations with itertools, preserving outside pair
orientations. It retains the warm-start value, incumbent, valid lower bound,
integer and relative gaps, zero exact-integer tolerances, degeneracy count,
lexicographically smallest minimizing window, measured CPU/wall time and
termination. Complete enumeration yields incumbent=lower bound. An interrupted
enumeration retains its incumbent and valid lower bound zero, is labelled
incomplete, and prevents physical evaluation. Each window has a 120 CPU-second
solver limit inside the study budget; no continuation silently enlarges it.

The independent certificate checker reconstructs each Pauli action with integer
phase exponents, checks supplied coefficient/order/weight normalization, and
uses a different DFS construction. Each appended unit contributes its pairs
once. A Lehmer-rank bitset rejects duplicate leaves and verifies coverage of all
9! feasible orders. It reconstructs full warm-start/incumbent objectives,
outside-order consistency, integrality, lower bounds, gaps and tie breaking.
The independent per-window enumeration limit is 120 CPU seconds. Literal sparse
Pauli matrices separately check complete-unit conservation, internal commutation
and floating response norms, without certifying floating error bounds.

Only after all three complete certificates are written, freeze six physical
schedules per input: first-order inherited HF heuristic, new certified window,
current placement and larger-search individual-Pauli; symmetric reuse of the
inherited heuristic and certified window. No second-order-specific objective is
claimed. At both 10U and 30U ideal Pauli-rotation caps (U=retained Pauli count),
use maximal affordable steps with exact central/boundary fusion. This yields
36 full-space circuits. No particle-number or spin-sector projection is used.
Retain every exact schedule, output state, cost, error and leakage value.

The stable physical error is squared projective distance `1-|overlap|`,
calculated as orthogonal-residual infidelity divided by `1+|overlap|` on
normalized states. Compare every pair of the six schedules within case/budget;
an absolute difference <=1e-12 is a tie. An error <=1e-12 is floor-limited;
exclude ratios when either side is at/below this fixed floor. Retain all rows,
losses and exclusions. Independent replay uses literal sparse Pauli matrices,
unfused circuits and a scaled-Taylor exact-reference approximation, checking
states within 2e-10 and metrics within 2e-12, plus exact agreement of all floor
and W/T/L decisions. These numerical checks are distinct from kernel proofs.

Construction costs distinguish inherited heuristic/current/large-Pauli setup
and search records from new domain setup and exhaustive solve CPU/wall time.
The inherited records are not new matched-time trials. Independent verification,
physical propagation/reference and full service resource accounting remain
separate. No classical speedup, whole-order or physical-error optimality claim
is permitted from restricted surrogate certificates. No compiler or hardware
resources are measured.

Finite cumulative phase ceilings, declared before running:

| Phase | CPU seconds | Wall seconds | Output MiB |
| --- | ---: | ---: | ---: |
| Focused tests | 120 | 900 | 64 |
| Representative end-to-end preflight | 240 | 1800 | 128 |
| Main domain setup, enumeration and physical evaluation | 600 | 3600 | 128 |
| Independent enumeration and physical replay | 600 | 3600 | 128 |
| Reporting, archive, extraction and delivery checks | 120 | 900 | 128 |

The preflight uses the largest 8-mode/104-word H2O target with every coefficient
scaled by 7/5, one complete 9! solve, independent complete enumeration, all
12 circuits, reporting, gzip archive/extraction and strict checks. Its measured
time, memory, temporary/retained disk and replay costs determine margin-bearing
panel estimates before the panel begins. A failed/interrupted attempt consumes
the same phase allowance; it does not reset the budget. If the required study
cannot fit, retain partial evidence and identify the specific resource decision.

One numerical/build/test worker, one BLAS/Numba thread and no overlapping
reviewer client. `guarded.py` enforces CPUQuota=25%, Nice=19, CPUWeight=5,
IOWeight=10, MemoryHigh=448M, MemoryMax=512M, MemorySwapMax=0. Start/resume only
above 5.5 GiB available; below 5 GiB freeze only the owned task service. Maximum
register 8 modes, at most 128 units/Pauli words, sparse pair data and at most
256 amplitudes per propagated state. No additional dependencies are installed.

Expected scientific outcomes may be negative: a restricted J_D optimum can tie
the heuristic or yield larger physical errors. All such outcomes answer this
bounded question and remain explicit. The completed 8! pilot and this new 9!
study are distinct restricted certificates. Neither answers a whole-order
optimality question.


## Representative preflight and frozen panel estimates

Before main-panel enumeration or errors, the complete altered H2O preflight at
source5791a352a66a77da437bd11876e322e4eaac8d1d passed: all362880 producer
orders, independent DFS/Lehmer coverage,12 full-space circuits, report,
archive/extraction and strict delivery. NumPy2.5.3/SciPy1.18.1; n=8/P=104;
all coefficients scaled by the predetermined7/5. Algorithm phases used
prepare0.02734CPU seconds, certification1.49406, propagation0.20151,
independent verification1.96920, reporting0.00246, packaging0.01177,
archive/extraction/check0.04788; total3.75726CPU/15.01581wall. The service
used4.024CPU/18.06653controller wall, process peakRSS66920KiB. Retained
delivery396496bytes, gzip42411bytes, all preflight outputs1215721bytes.
These altered-target observations are not the prospective panel.

For the three-case panel, estimate main setup/enumeration/propagation60CPU
seconds/300wall and independent enumeration/replay60CPU/300wall. Treating
all three inputs as the largest shape, these CPU estimates exceed ten times
the corresponding preflight phases per case. This supplies margin for unit,
destination and internal-pair variation. Estimate reporting/packaging/extraction
8CPU/60wall, peakRSS160MiB, main retained/temporary disk16MiB and reporting/
extraction32MiB. No setup or phase is granted extra execution solely because
its estimate differs. The original cumulative600CPU/3600wall/128MiB main
and independent ceilings, and120CPU/900wall/128MiB reporting ceiling remain.

All failed attempts remain charged and retained. Preflight01 stopped before
enumeration because optional removed-scalar metadata was absent. The metadata
repair preserves absence and scales supplied values exactly. Preflight02
completed the altered solve/circuits but caught an independent-checker bug:
ordinary sin requires gap/(2den), whereas NumPy's normalized sinc already
includes pi. The producer's supplied-weight formula was correct and unchanged;
the independent check was corrected without loosening its5e-15 tolerance.
A hand-derived cosine-identity regression and the full13-test suite passed.
No main input, coefficient, window-selection rule or error threshold was tuned
using these failures; both partial attempts remain separate from preflight03.
