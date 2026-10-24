# Prospective symmetric-formula response ordering

Issue #8 asks for a rule that uses the symmetric formula's own leading defect,
rather than only applying first-order-trained orders to a symmetric circuit.
This protocol fixes the study before any new physical errors. The derivation
in [DERIVATION.md](DERIVATION.md) and all source are committed before preflight;
preflight measurements and margin-bearing estimates are appended before main.
No failed two-level reference correction is retuned or accepted.

Inputs are the three authenticated supplied prospective-reference models
Li2_r2p95 (6 modes/2 electrons), F2_r1p55 (8/6), H2O_scale1p10 (8/4), without
regeneration or orbital changes. Rational coefficient numerators, denominators,
scalar-removal policy and complete conserving units remain authenticated.
The determinant is HF, occupation (1<<electrons)-1. New times are 0.83 and 2.11
inverse Hartree. This small related family is descriptive evidence, not three
independent statistical replications or larger-molecule validation.

For chronological symmetric half/full/half steps, the producer evaluates the
actual third-degree Taylor polynomial on the determinant. With q3 the real
coefficient of (-i h)^3, define A3=q3-H^3|HF>/6. Verify both lower coefficients,
q1=H|HF> and q2=H²|HF>/2, for every evaluated order. Search does not use exact
Hamiltonian time evolution or a measured physical error. The diagonal reference
D is all diagonal Pauli terms. The new fixed surrogate is

`J_D2(pi,T)=sum_(w!=HF) T² sinc²(T(Ew-EHF)/2) A3_pi(w,HF)²`.

At fixed T and ideal rotation cap, optimize `J_D2/r_pi^4`, where r_pi is the
maximal affordable symmetric step count after exact center/boundary fusion.
The common positive T^4 factor omitted within each cell cannot change a ranking.
Supplied coefficients and binary64 weights are not an interval certificate;
diagonal transport approximates full-H transport. No exact physical-error,
global ordering or universal mechanism guarantee is claimed.

At each time, construct a fresh conserving first-order J_D control with the same
time-dependent weights, and a fresh individual-Pauli J_D control with a larger
nominal search allowance. Both use three seeds: fixed label order, decreasing
unit coefficient L1 order, and a compatible current-placement warm start. The
legacy current placement is separately named and retains its selection time T=1.
Conserving control search has 2 CPU seconds, 1000 attempted candidates total and
at most 5 accepted moves per start. Larger individual-Pauli control has 8 CPU
seconds, 20m(m-1) attempted candidates total and at most 10 moves per start.
First-order r does not depend on endpoint, so the same order is reused across
the two caps, rather than claiming a new search at each cap.

For each time/cap, cubic seeds are the freshly trained conserving first-order
order, its complete reverse, and current placement. An orientation-only control
chooses the lowest new score among these fixed three seeds, breaking ties by
the lexicographic order. Cubic insertion uses the same seeds, 2 nominal CPU
seconds, 1000 attempted positions total and at most 5 moves per start. Each
search divides CPU/candidate allowances among three starts without carrying
unused allowance. Retain attempted/completed candidates, incomplete sweeps,
interrupted candidate stages, accepted moves, independent final rescoring and
termination. A completed improvement from a partial sweep remains admissible.
The score and affordability of every seed and accepted order are checked.

Nominal allowances cover bounded insertion only, not complete constructors.
Measure model loading/base setup, legacy placement, method setup, initializer,
seed scores, cubic cache setup/cleanup, context cleanup, search and constructor
return separately. Standalone new-method charges include the upstream fresh
first-order initializer and all required common/cache work. Constructor charges
repeated across physical rows are not summed as actual campaign work. Shared
construction and full process resource ledgers remain separate. This is not a
matched-total-CPU comparison; all actual overruns and classical costs remain.

Eight methods per time/cap are fixed:

- Current placement, first-order and symmetric.
- Fresh conserving first-order J_D at the same T, first-order and symmetric.
- Larger-search individual-Pauli J_D at the same T, first-order and symmetric.
- New cubic orientation-only rule, symmetric.
- New cubic bounded insertion, symmetric.

Each of 3 models has 2 times, 2 caps and 8 methods: 96 full-space circuits.
Common ideal Pauli-rotation caps are 10U and 30U, U=retained Pauli count.
For unit lengths l, first-order cost is rU; symmetric cost is
`r*(2U-l_first-l_last)+l_first`. Use the largest integer r with cost<=cap.
These counts use actual complete internally commuting units with exact center
and boundary fusion. Formula-changing contrasts do not isolate ordering.
All inputs, times, controls, scores, orders, costs and complete Pauli circuits
are frozen in hashes before any main physical propagation. No target, seed,
threshold, loss or floor case is selected using its physical error.

Physical propagation is full 2^n, with no particle/spin-sector projection.
Retain exact-reference and propagation costs, actual steps/rotations, every
output state, stable projective error 1-|overlap|, infidelity, norm, particle-
number leakage and joint (N_alpha,N_beta) leakage. The fixed absolute W/T/L
tolerance and floor are 1e-12. Compare all 28 pairs per cell, retaining all
336 decisions. Exclude paired ratios when either side is at/below the floor;
retain the underlying rows and floor counts. Report paired ratios/differences
before descriptive summaries. Include every strong control and every loss.

Independent verification constructs literal sparse Pauli matrices. It checks
unit conservation and internal commutation, and reconstructs the full signed
cubic column with the separate nested-commutator formula
`sum_i(-[Gi,[Gi,Btail]]/24+[Btail,[Btail,Gi]]/12)`.
Checking the full column precedes projected norms, so a sign error cannot hide
behind a square. Ordinary sin(x)/x supplies an independent weight evaluation.
It checks input/seed/control identities, every accepted insertion, actual
work/cost closure, frozen schedules and affordable steps. Independent unfused
circuits and a scaled-Taylor exact-reference approximation check all states
within 2e-10 and metrics within 2e-12, plus exact agreement of every W/T/L and
floor decision. Score agreement tolerance is rtol 2e-9/atol 5e-18 and full cubic
column absolute tolerance is 2e-11. These are numerical checks, not Lean proofs.

Finite cumulative resource ceilings, fixed before tests or preflight:

| Phase | CPU seconds | Wall seconds | Output MiB |
| --- | ---: | ---: | ---: |
| Focused tests | 120 | 900 | 64 |
| Representative complete preflight | 240 | 1800 | 128 |
| Main construction and physical evaluation | 600 | 3600 | 128 |
| Independent objectives and full replay | 600 | 3600 | 128 |
| Reporting/archive/extraction/compact checks | 120 | 900 | 128 |

Preflight uses only largest H2O with all coefficients scaled by predetermined
7/5, T=2.11, all 8 arms and both caps: 16 circuits, all objectives/replay,
report/package/gzip archive/extraction/strict check. No main validation error
is observed. Measure compute, peak memory, retained/temporary/archive disk,
packaging/extraction and replay, and freeze margin-bearing main estimates
before executing main. Every failed/interrupted attempt consumes its original
phase ceiling and remains retained. If essential work does not fit, preserve
evidence and request the specific resource decision; no silent enlargement.

One local numerical/build/test worker, one BLAS/Numba thread and no overlapping
reviewer client. CPUQuota=25%, Nice=19, CPUWeight=5, IOWeight=10,
MemoryHigh=448M, MemoryMax=512M, MemorySwapMax=0. Start/resume only above 5.5 GiB
available; below 5 GiB freeze only this task's service. Bound registers at 8
modes, 128 units/Pauli terms, 256 state amplitudes and 8 MiB cubic action cache.
Use the existing verified Python3.12 NumPy2.5.3/SciPy1.18.1 environment with no
additional dependency installs. All worker phases use guarded.py.

Scientific conclusions may be negative. Improved cubic score can yield larger
physical error, and second-order controls may dominate first-order methods.
Retain those outcomes. Ideal rotations, future compiled resources and hardware
measurements remain distinct. No conservation-only, practical speedup,
asymptotic scaling or larger-molecule conclusion follows from this panel.
