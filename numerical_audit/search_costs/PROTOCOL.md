# Controlled classical search costs and declared scaling range

Issue #5. This protocol takes effect only after the bounded representative
preflight, a source commit and PROTOCOL_FROZEN/SCHEDULES_FROZEN records. No
physical error may be evaluated before all search schedules are frozen. The
current implementation is a draft until those steps and required reviews pass.

## Fixed inputs and factorial design

Reuse the authenticated Li2_r2p95, F2_r1p55 and H2O_scale1p10 retained molecular
models. Add exactly five synthetic extended-Hubbard targets with n=4,6,8,10,12
spin orbitals, half filling (first n/2 occupations), T=1 inverse Hartree and
hbar=1. The synthetic family is not new molecular active-space evidence.
The orbital basis remains fixed; complete conserving Hamiltonian units are
ordered. Preserve all cases and losses; do not replace a difficult input.

For synthetic inputs, epsilon_p=(2(floor(p/2)+1)+(p mod 2))/20; equal-parity
hopping t_pq=-3/[10(1+(q-p)/2)]; density U_pq=3/5 on the same spatial site,
otherwise 1/[10(1+q-p)]. Use all p<q density pairs. JW hopping uses complete
XX/YY partners and intervening Z strings. Diagonal words move individually.
Retain exact rational coefficients as integer numerators/common denominator,
and record the omitted scalar exactly. No coefficient or geometry is tuned.

Eight factorial optimizers: conserving/individual-Pauli partition x J_D/J_0
objective x cold/placement initialization. Cold starts are label and decreasing
coefficient-l1 orders. Placement starts are current-placement and the same l1
order. Flatten the conserving placement schedule into a compatible individual-
Pauli warm start, charging its conversion. Keep initializer and objective
effects separate in comparisons. The rejected two-level correction remains
rejected and is not retuned or included as an accepted initializer.

## CPU, setup, warmups and actual work

Total stand-alone construction CPU caps: 0.5 and 2 seconds. Process CPU time,
one thread. Separately charge common setup (stored input parse, target/unit
inventory and common diagonal/off ordering), initializer, per-method response
setup and search. Current placement is measured once per repetition; its full
construction cost is charged to each warm-start method. Pauli warm conversion
is additionally charged only where needed. Actual campaign execution is
reported separately without counting a cached initializer multiple times.
Original molecular integral generation is reused, not called new search work.

Three repetitions. One identical 0.5-CPU warmup per factorial arm/input before
measurement; no JIT is used. Predeclare cyclic arm order by repetition.
Every arm uses two starts, at most ten accepted moves/start and at most
20m(m-1) examined positions. Equal CPU caps do not imply equal arithmetic or
examined positions. Count every position, including zero commutators, and
record partial sweeps, setup visits, accepted moves and exact termination.

CPU charge for shared setup/placement is subtracted before local construction.
Deadline checks occur during pair setup and sparse indexing, and before search
candidates. A 0.1-CPU reservation is required for a full objective rescore;
candidate scans also retain a 0.002-CPU margin. Report measured rescore times
and any reservation violation. Setup may stop incomplete. A valid initializer
then remains the frozen order, with search score null, not falsely completed.
If an initializer/shared charge alone exceeds a cap, label the result outside
budget. Retain all actual work and overrun amounts; do not claim such a record
met its declared CPU cap or silently lend it more time. Record bookkeeping/GC
overhead so local and charged phase accounting closes.

Current placement is a separately timed baseline. A stronger individual-Pauli
J_D control gets 8 CPU seconds and the same starts/move/evaluation rules; label
its unequal classical allowance outside the matched-CPU factorial.

## Frozen circuits, verification and interpretation

There are 54 recorded orders per input (48 factorial, three larger-Pauli and
three current baselines), 432 total. Freeze every order before physical errors.
Evaluate first-order and symmetric reuse of each order at 10U ideal Pauli
rotations, 864 logical circuit comparisons. Symmetric execution uses exact
central/boundary fusion and the largest affordable integer step count. It
reuses a first-order objective; this does not complete issue #8. Retain full
2^n states, including leakage, absolute overlap error, infidelity, norm,
particle-number leakage and joint (N_alpha,N_beta)-sector leakage.

Independent sparse literal-Pauli scoring checks every frozen order and both
J_D/J_0 values, without crediting validation scores as in-budget search work.
Audit named arms, exact target permutations, all CPU/evaluation/sweep ledger
identities, objective values and unit boundaries. Replay circuits with an
unfused formula and a distinct sparse rotation implementation. Compare the
producer's SciPy exponential reference with an independently scaled Taylor
reference, retaining numerical tolerances and limitations. A numerical check
is not a Lean proof, and no optimized surrogate certifies physical accuracy.

Report each input size/structure, all repetitions and paired differences before
aggregation. Do not subtract separate medians to create a paired result.
Incomplete and out-of-budget records remain visible. Efficiency conclusions
are restricted to this implementation, budget policy and measured range;
larger synthetic mode count is not larger-molecule validation. Single-thread
CPU under a shared-machine quarter-core quota does not establish dedicated-
machine wall speedups. Ideal rotations are not compiled/hardware resources.

## Storage, resource preflight and finite campaign allowances

Refuse n>12/P>256 before exponential physical allocations. Search itself uses
sparse pair/destination responses, not a 2^n state array. Before pair allocation,
refuse more than 256 units or a conservative 1024 bytes per possible pair
estimate above 64 MiB. Record actual pair/destination counts, array payload,
the estimate and peak RSS across every declared size. Literal/circuit caches
are separately bounded at 64 MiB. The OS resident-memory cap remains binding;
Python container estimates are not a theorem about total RSS.

One local numerical/build/test worker at a time, one BLAS/Numba thread, one
build job. CPUQuota=25%, Nice=19, CPUWeight=5, IOWeight=10, MemoryHigh=448M,
MemoryMax=512M, MemorySwapMax=0. Start/resume only above 5.5 GiB available;
freeze only the owned service below 5 GiB. Serialize reviewers with heavy work.
Use the existing pinned Python 3.12 NumPy/SciPy environment, with pytest only
for tests; no unrelated installation or saved client-setting changes.

Finite preflight allowance: 120 CPU seconds/900 wall seconds/64 MiB outputs.
Use the largest declared synthetic structure with deterministic altered
coefficients; measure end-to-end setup, search, propagation, independent replay,
extraction and packaging, without inspecting panel ordering errors. The representative
preflight uses one repetition of all eight factorial arms at both CPU caps, the
larger-Pauli control and current placement (18 orders/36 circuits), with one
0.5-CPU warmup per factorial arm. Repetition/range expansion is estimated with margin. Before the
panel freeze, record estimates with at least 2x margin. Main construction plus
propagation: 1800 CPU seconds/10800 wall seconds/512 MiB outputs, shared across
phases and failed attempts. Independent verification: 1800/10800/512 MiB.
Extraction, reporting and replay: 300/1800/128 MiB separately. Preserve failed
attempts, partial work and cumulative use. Do not enlarge ceilings or quietly
substitute an easier scaling question if an essential step cannot fit.

The required independent Astra/Argon local Review/Address/Verify loop and
deferred publication remain distinct from the separate authorized merge phase.
No issue closure before complete criteria, merged provenance and post-merge CI.

## Frozen preflight estimates

Before prospective construction/errors, the altered largest12-mode/138-word
end-to-end preflight passed:18orders/36circuits,all eight factorial warmups,
construction14.556CPU seconds,propagation4.062CPU,independent verification5.300CPU,
extraction/packaging/replay0.148CPU;total24.067 algorithmCPU/96.295 wall seconds.
The guarded service measured24.674CPU/99.443controller wall seconds; full process
RSS is retained in the preflight resource ledger. All propagated-state, response
and metric comparisons passed. These are preflight observations, not prospective
ordering results. The altered target is scaled7/5 relative to the declared
synthetic12-mode target and is never included in the measured panel.

For the full8-input/3-repetition panel, treating every input as the largest
shape and multiplying every representative phase by24 then at least2 yields
conservative estimates:construction800CPU seconds,propagation250CPU seconds,
independent verification350CPU seconds,and extraction/reporting/packaging/replay
30CPU seconds. Estimated main-phase wall time with quota and margin is4800seconds,
verification1600seconds,extraction180seconds;peakRSS192MiB,main temporary disk
256MiB,retained panel/report/delivery80MiB,and extraction temporary disk128MiB.
Smaller inputs are retained rather than substituted; these bounds deliberately
overcount repeated warmups and largest-shape storage. The original finite
1800CPU/10800wall/512MiB main and independent-verification ceilings and
300CPU/1800wall/128MiB extraction ceiling remain binding and are not enlarged.
