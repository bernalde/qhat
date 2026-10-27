# Response to the independent numerical project review

The reviewed sources are QHAT research merge `c28aba45644a817bc200c946523807569c2956c2`
and paper research merge `71608ab0f1709372b9c2a994d1ace3588e4b9ade`.
This response covers all 21 numerical/core/Julia items in the completed initial
review sections A2, B and A5. The separate paper response covers the six region C
items. Whole-project review means the reported sampled scope, not every source
line or raw evidence row. Claude performed source inspection and selected scalar
and rendered-text inspection, with no project execution. Original reports and
their limitations remain separate from coordinator dispositions.
The client reported `claude-opus-5-5`; high effort was requested.

### CLD-A2-001 — Done

Unsupported mapping names now raise the intended `NotImplementedError` from
`HamiltonianConfiguration.fermion_to_qubit_name()`. This also prevents
`f2q_tag()` and filename construction from silently using `None`. The defensive
dispatch branch reads the configured mapping string without calling it.
Every previously supported JW/BK alias keeps its name and filename tag.

The focused regression reproduced eight failures and fourteen passing alias
controls on the unchanged reviewed baseline. It exercises both configuration
accessors, filename construction and the defensive dispatch error path without
a chemistry calculation. See
[`test_hamiltonian_mapping_validation.py`](../../analysis/tests/test_hamiltonian_mapping_validation.py).

### CLD-A2-002 — Already satisfied

The missing-test premise does not hold at the reviewed source. Existing tests in
`analysis/tests/test_time_evolution_audit.py` check the 4-qubit JW/BK determinant,
mapped number operator and evolution agreement. Its parameterized
`test_bk_reference_matches_every_mode_occupation_for_non_power_of_two_sizes`
also checks each mapped occupation at 3, 5 and 6 qubits and several electron counts.

The pinned OpenFermion implementation constructs `bravyi_kitaev_code(n).encoder`
from the same `_encoder_bk(n)` used by the generator, so the alleged two-encoder
discrepancy has no source basis here. The existing regressions are reused rather
than adding a duplicate test.

### CLD-A2-003 — Done

The root README now links the time-evolution entry point and numerical study tree,
including the practical-accuracy, transfer and molecular-certificate studies.
It also links this response so the independent-checking boundaries are discoverable.

### CLD-B-001 — Done, differently

The current `sufficient_epsilon_domain(a, b, c, d)` is a sufficient condition for
the **ordered claim that the first slope is smaller than the second**. In the
fixed report, `a` is current placement and `b` is diagonal prefix. Its certified
gap is `b.lo - a.hi`, and the positive-coupling threshold is that gap divided by
`c + d`, subject to the stated strictness and coupling-at-most-one restriction.
`no_certified_slope_gap` means this forward sufficient condition is unresolved;
it does not classify the reverse ordering. The separate interval `winner()`
function already checks both directions.

This explicitly documents the helper's orientation instead of broadening the
frozen two-schedule reporter into a general winner-selection API. In the retained
LiH data the forward gap is positive and the algebraic threshold is approximately
0.03319. The physical leading guard remains unresolved, and the eight sampled weak
guards and 66 independent finite-step comparisons retain their original labels.
The existing reporter, report and archive bytes remain unchanged.

### CLD-B-002 — Done, differently

The practical verifier's mathematical paths independently reconstruct objective
values, primitive versus ideal action, DAG resources and target selections.
Its input authentication, dispatch naming and layout validation also reuse
producer bookkeeping helpers (`study.bind`, `dispatch`, `curve_name` and
`validate_record`). Those shared checks alone are not an independent proof of
the bookkeeping contract. This boundary applies to the original replay as well
as the current verifier.

A separate standard-library structural audit now checks all 81 retained
construction layouts and the complete cyclic dispatch plan from literal
Hamiltonian terms and the declared protocol. It imports no producer bookkeeping
helpers. Negative cases require it to reject a missing rotation and a duplicated
unit. See
[`test_practical_retained_structure.py`](../../analysis/tests/test_practical_retained_structure.py).
The regular analysis test lane exercises this additional audit. The historical
verifier and measured replay remain intact; no historical run is relabeled as
having executed this new check.

### CLD-B-003 — Already satisfied

`transfer_holdouts.study.constructor()` measures CPU through the complete
constructor return and overwrites the standalone charge with
`reference_charge + caller_CPU_seconds`. Its existing
`test_outer_training_charge_includes_final_constructor_return` deliberately gives
the inner constructor a 1-second charge and the outer caller an 8-second charge,
then requires an 8.5-second total including the 0.5-second reference charge.
This is the requested overhead boundary and regression protection.

### CLD-A5-001 — Done

The matrix-free Pauli action now matches the dense left-to-right Kronecker order:
character zero acts on the most significant bit. Y sends zero to i times one and
one to minus i times zero. Real inputs and SciPy matrix-column calls retain
complex phases. Exhaustive one-, two- and three-qubit Pauli basis-column tests,
complex vectors and adjoints compare against independently built dense matrices.
The unchanged implementation failed all seven initial operator/coefficient
regression cases before these corrections. Old tests encoding the wrong basis
convention were deliberately migrated, rather than preserved as false oracles.

### CLD-A5-002 — Done

Automatic coefficient estimation and either circuit backend receive the same
once-reordered Pauli dictionary. The estimator no longer uses a separate loader
order. Call-boundary tests capture the estimator and builder inputs under
magnitude and lexicographical ordering.

### CLD-A5-003 — Done

The reference and fast estimators now use the complete ordered C21 index set,
including i=j and both distinct orientations, and uniform C22 pairs. Exact
enumeration, Monte Carlo sampling, tracking keys, completion counts and memory
reservations follow that domain. For N terms it contains N(N-1)(2N-1)/6 triples.
Uniform rejection sampling avoids allocating a cubic sampling population.

Independent dense nested-commutator oracles include [X,Z], which requires
(C1,C2)=(1,1/2), and [X,2Z] versus [2Z,X], which give C2=5/3 and 4/3.
Deterministic sampling-domain tests and asymmetric nested terms protect the
orientations. Existing 1/6 expectations were missing the diagonal contribution.
The formulas are the triangle-inequality expansion of
[Childs et al., Eq. 152](https://arxiv.org/html/1912.08854v3#S5.E152).
Monte Carlo outputs remain estimates, not certified upper bounds. Historical
speedup figures are explicitly not a new benchmark of the corrected domain.

### CLD-A5-004 — Done

Automatic fourth-order sizing now fails before estimation with an instruction
to supply `trotter_steps`. Manual fourth-order support remains. The removed
C2-based heuristic was dimensionally inconsistent; no C4 certificate is invented.

### CLD-A5-005 — Done

Driver and dispatcher consistently trim and lowercase the method name. The
dispatcher uses exact tuple membership, rejecting empty strings and substrings
such as `r`, `ramped` and `trotter`. Supported configuration factories already
provided the canonical name; the correction covers manually changed objects.

### CLD-A5-006 — Done

All evolution times must be finite and positive. QPE and phase-derived
eigenvalue-error analysis additionally enforce the centered physical branch cut
and configured phase margin. Ordinary or controlled time evolution may use a
longer override when phase-derived eigenvalue errors are disabled. Automatic
QPE precision uses ceil(log2(2π/(time*energy_error))) at the final time. The
automatic time is now the largest margin-safe time, without counterproductive
power-of-two energy padding; default times and resource estimates can change.
At bounds [-1,1], energy_error=.1 and factor1.01, the strict minimal precision is
5 qubits, while the discarded padded calculation chose6. Finite positive factors
below one still permit safe overrides and clamp automatic times below the
physical cut. Explicit phase registers retain their authority in the builder.

### CLD-A5-007 — Done

The analysis README and wrapper docstring now describe the `exp(-iHt)` convention
for the flattened/error-analysis route, bijective eigenvector-overlap matching,
principal phase ratios and their alias period, and evolution-operator/state
metrics. Dense-dispatch limitations and retained unused compatibility arguments
are explicit. The historical `original` backend's opposite time sign is called
out; historical reproduction semantics are not silently changed.

### CLD-A5-008 — Done

Selecting unsupported `LCPS` now reaches the supported-source error instead of
an undefined loader. The unfinished boson helper raises `NotImplementedError`
and points to the implemented binary mapping instead of returning None.

### CLD-A5-009 — Done

Dense text, sparse text, JSON and HamLib loaders accumulate repeated coefficients
before Hermiticity validation. Sparse spellings are sorted before accumulation.
Tests cover duplicate sums, cancellation, identity terms, reordered spellings,
imaginary contributions that cancel and a genuinely nonreal final coefficient.
No historical input data changes.

### CLD-A5-010 — Done

The optional coloring helper imports Gurobi lazily and reports a missing optional
dependency clearly. Ordinary reorder functions retain no Gurobi requirement.
The dependency-error test does not claim a licensed optimization was executed.

### CLD-A5-011 — Done

Julia operator errors use the largest singular value rather than spectral radius.
The iterative path scales tiny errors, checks convergence and two-sided residuals,
and handles zero operators. A four-dimensional difference of two unitary matrices
has norm 2 and spectral radius sqrt(2), independently exposing the original bug.
The iterative result remains a floating-point estimate, not a rigorous enclosure.

### CLD-A5-012 — Done

Julia documentation distinguishes global error at fixed total time from local
one-step error. Second-order global scaling is O(dt²), as used by the script.

### CLD-A5-013 — Done, differently

The omitted Chebyshev backend rescales the problem: its default Bessel arguments
are 1/3 for state propagation and 1/4 for operator propagation, with 12 and 16
compositions at total phase π. The reported unscaled order-10 tail estimate
therefore does not apply. There is no demonstrated reference-dominated panel.
The available finite-tail metadata is now retained by the caller and exposed
through an opt-in matrix-return API. Output explicitly calls the reference an
approximation and excludes arithmetic roundoff from the estimate. Small dense
exponential parity tests provide numerical checks, not an interval certificate.

### CLD-A5-014 — Done, differently

Both energy CLIs identify their full-Fock-space spectral convention and the
unverified sector provenance of supplied metadata. The repository metadata
generator also targets full-space extrema; a stored sector switch was not proven.
HF-start Arnoldi is explicitly a Krylov seed, with overlap/invariant-subspace
limitations, rather than a sector projection or global-minimum guarantee.
The JW-specific HF route rejects BK/unknown mapping metadata before reconstruction.
Projecting a leaking unitary onto a sector is deliberately not used as a repair.

### CLD-A5-015 — Already satisfied

The supported configuration factory always sets `phase_scale_factor`, defaulting
to 1.01; its existing regression checks that value. The fallback did not establish
a default aliasing bug. The driver's compatibility fallback is also aligned to
1.01 as part of the final-time validation correction.

## Validation and preserved scope

The first six items retain their earlier evidence: eight expected mapping failures
and fourteen passing alias controls on the unchanged baseline, followed by 32
passing focused checks. The new operator/coefficient regression reproduced seven
failures on the unchanged implementation. At the corrected source, 151 focused cases passed in Python 3.12 with Numba
JIT disabled and the task-owned tomlkit 0.13.3 parser. Eight new full-library
call-boundary cases were deselected because this lightweight environment lacks
Qualtran/pyLIQTR. The complete Python application CI also runs those cases, the
fourth-order backend tests and actual JIT execution. It remains a separate gate;
this local result is not full-environment CI parity.

No frozen study input, measured row, numerical report, archive, producer/verifier
snapshot, or circuit changes. Recent numerical study entry points do not invoke
these older analysis estimators or matrix-free operators. The separate paper PR
corrects status wording and rebuilds document artifacts, preserving every audited
Lean byte and earlier semantic assessment. No numerical panel or proof build is
relabeled as newly executed. Unknown quantum preparation cost, restricted domains,
unequal consumed CPU and software-versus-kernel boundaries remain explicit.

The bounded local Julia component suite passed using the exact repair source, Julia 1.11.9 and cached Arpack 0.5.4. The complete declared-project reference/CLI suite remains a distinct hosted CI gate.

## Original report identities

The original public reports are retained separately; private provider/account
metadata is excluded.

| Completed section | Original report SHA-256 |
|---|---|
| A | `99b099810473195fdbeeea4ca112e3756379de242eab8dc5cbe920ff0a8ed830` |
| A2 | `9bd58a973745eeda1aaf866fb2f574c0c86396d56bc4925505509565d04a62a2` |
| A5 | `d0aee63a0323db9bbd71283b310166ec33827c12a373a6e4ca0c60b7188367e9` |
| B | `3cae6a6b29d244041859561b1d9e809aa1e931469b0516900925c8fc2a1213c3` |


## Targeted review follow-up: core

The second review caught distinct default-path and configuration edges after the
initial regressions passed. These are corrections to the current repair, not
new measurements of the retained studies.

| Item | Disposition and evidence |
|---|---|
| R2-QC-001 | Done: remove counterproductive time padding, retain strict energy resolution; actual default-driver and independent minimal-register fixture. |
| R2-QC-002 | Done: factory validates and serializes mode, auto-exact and time limit; exact settings reach the estimator; unknown options fail clearly. |
| R2-QC-003 | Done: the independent81-layout audit binds model file SHA-256 and asserts canonical block union equals coefficients; corruption controls reject both gaps. |
| R2-QC-004 | Done, differently: phase alias guards apply to QPE/eigenvalue inference; finite-positive ordinary/controlled long-time overrides remain valid. |
| R2-QC-005 | Done, differently: preserve finite positive factors below one with physical safety checks and an accurate warning. |
| R2-QC-006 | Done: the analysis loader names supported JW/BK choices in a ValueError; its existing accepted strings and callable forms are preserved. |
| R2-QC-007 | Done: neither estimator may return zero for an empty sample phase; an expired budget raises clearly. Finite-time validation also rejects NaN. |
| R2-QC-008 | Done: matrix-free documentation states memory scaling and the remaining Python-loop wall-time limitation. |

Full-library driver/encoder integration and real JIT remain required hosted checks.
Earlier local counts describe their original revisions; follow-up execution is
recorded separately. No source-independent Monte Carlo certificate is claimed.


## Targeted review follow-up: Julia

| Item | Disposition and evidence |
|---|---|
| R2-QJ-01 | Already satisfied, with a new real-header regression: hamgen.py writes the exact metadata key and committed JW data uses jordan-wigner, which the guard accepts. |
| R2-QJ-02 | Done: a synthetic nondegenerate two-qubit fixture runs both actual single-file CLI solvers, both bond-parser solver paths and metadata-free fallback, plus the changed demonstration reference/operator call chain against dense oracles. Tiny-system Krylov dimensions are clamped to n and unconverged Arnoldi output fails clearly. Full demonstration main and molecular sweeps are outside these tests. |
| R2-QJ-03 | Done, differently: README explicitly carries the full-space/unverified-sector caveat into saved CSV/PNG artifacts, describes their legacy labels and metadata-versus-eigensolve source, and requires preserving metadata and warnings with exports. No sector switch is inferred. |
| R2-QJ-04 | Done: demonstration documentation uses approximate reference and estimated errors, compares observed ratios to expected asymptotic scaling, and removes the unimplemented recommendation claim. |
| R2-QJ-05 | Done, differently: pin checkout to the verified official v6 commit and bound inspected Arpack0.5/ArnoldiMethod0.4 APIs plus Julia1.11. No Manifest is committed; remaining exploratory packages resolve from the registry and full Pkg.status provenance remains required. These are compatibility checks, not a frozen full environment or local/hosted parity claim. |

These additional accepted-path tests have not run in the lightweight local
component environment. The full declared-project hosted lane must execute them.

Follow-up execution: 183 focused Python cases pass; 22 full-library cases remain
deselected locally. This includes the new pure phase/configuration, expired-sample,
model-binding and corruption checks. Numba JIT remains disabled locally. Paper
source/retained-evidence checks pass with exactly the intended pending-fourth-review
and missing-pair errors. The full Python/JIT and Julia hosted jobs remain required.
