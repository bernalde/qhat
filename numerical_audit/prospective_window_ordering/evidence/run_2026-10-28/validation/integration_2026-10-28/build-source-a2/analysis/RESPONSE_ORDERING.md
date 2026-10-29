# Sparse response ordering for molecular units

`sparse_response_ordering.py` constructs a permutation of supplied commuting
Hermitian Pauli units for a symmetric second-order product formula. It uses a
reference state, an evolution time and an ideal Pauli-rotation budget. Ordering
is a bounded classical optimization; it uses no exact evolution or validation
error. It is a heuristic for improving the response objective, not a global
optimizer or a certificate that one order has smaller finite physical error.

## Inputs and conventions

Words use OpenFermion-style tuples such as `((0, 'X'), (3, 'Z'))`; qubit zero is
the least significant occupation bit. Coefficients are real. Terms within each
unit must commute, allowing the unit exponential to be implemented as its Pauli
rotations. `molecular_units(terms, n_qubits)` groups compatible nonzero flip masks
and keeps diagonal terms separately. It rejects a noncommuting equal-flip group.
Identity terms are retained if supplied. A reference is one occupation integer
or a normalized mapping with at most sixteen nonzero complex amplitudes. This
classical representation supplies no quantum state-preparation algorithm.
Grouping alone does not certify particle- or spin-number conservation of an
arbitrary input; the molecular benchmark checks that property separately.

Time and coefficients use consistent units with hbar=1. Training uses floating
operator/state arithmetic and trigonometry. Diagonal gaps are summed from exact
supplied rational coefficients before conversion to a floating phase. Conservative certificates
instead interpret supplied rational values exactly, including the exact dyadics
represented by floats. They concern ideal unit exponentials; numerical rotation,
synthesis, SCF, integral, basis and hardware errors are separate.

## Response objective

For a permutation pi, the symmetric step is the forward half-step sequence,
the final full unit, and the reverse half-step sequence. The implementation
propagates only Taylor degrees zero through three on reachable bit strings.
Subtracting H^3/6 gives the cubic defect column C_pi. The common physical factor
i multiplying that column is omitted consistently.

Response actions first remove supplied identity terms and subtract each unit's
diagonal expectation at the smallest occupied reference basis index. These
scalar phases cancel from the cubic defect; the original units, public operator
actions and rotation lengths remain available without that centering. Diagonal
transport uses energy differences directly, avoiding subtraction of large total
energies. A bounded transport cache retains at most `max_entries` pairs.

Full columns and transported responses are refused when their norm is smaller
than 2^-24 times the corresponding pre-subtraction column scale. This screens
severe floating cancellation, including large scalar-free coefficient ratios;
it is not a forward-error certificate. Exact mutually commuting Pauli expansions
are handled as zero defect without cancellation. A `NumericalResolution`
refusal produces `status="numerical_limit"` during training and never a completed
zero score. The last earlier fully scored incumbent remains usable.

For each occupied reference component v, the destination u receives

`F_uv(T) = T sinc(T (D_u-D_v)/2) exp(i T (D_u-D_v)/2)`,

where sinc(x)=sin(x)/x and D is the diagonal of the supplied Hamiltonian. Columns
are weighted before coherent summation and projection perpendicular to the
reference state. The score is the resulting squared norm divided by r^4, where
r is the affordable fused symmetric step count. This matches the preceding
small dense response objective, while avoiding full Hilbert-space arrays.

For unit lengths l_i and U=sum l_i, r repeated steps cost
`r (2U-l_first-l_last) + l_first` ideal rotations. A single commuting unit needs
only U rotations and is exact. This is an ideal rotation model; actual line-routed
CX and depth are evaluated separately by the molecular benchmark.
The maintained `strang_schedule.plan_strang_schedule` supplies the affordable
step count and a lazy pulse iterator; the orderer reuses that endpoint.

Sparse support and cached unit-column entries have explicit bounds. A bounded
LRU avoids storing full unit matrices. Support may still grow substantially;
avoiding 2^n allocation does not establish inexpensive training for every molecule.
Timeouts and support refusals retain only fully completed scores. An initializer
can be returned with `score=None`; this must not be interpreted as zero error.
Training uses insertion moves from declared initial orders, retaining completed
candidates if a later sweep is interrupted. The default starts are magnitude,
reverse magnitude and canonical order; callers may supply a cheaper specialized
placement as an initializer. Context work, seed scores and search consume the
reported CPU allowance. Validation that cannot complete raises rather than
claiming a valid result. Existing frozen scientific implementations are unchanged.

For larger unit lists, insertion destinations are visited by adjacent swaps. If
A and B are adjacent and R is their current suffix sum, the new-minus-old cubic
column is `[(A+B)/8 + R/4, [A,B]]` applied to each reference component. This avoids
replaying the entire Taylor sequence for every trial. The source-specific
transport and affordable step count are recomputed before comparison. A full
Taylor rescore is required before an updated incumbent is accepted; interrupted
or numerically inconsistent proposals cannot replace the last checked incumbent.

Insertion search prioritizes the best fully scored seed and retains only that
seed's column snapshot, bounded by `max_entries`. One rescore reserve belongs to
the whole remaining allowance; the remaining search time is divided among seed
slices. Completed seed, incremental-candidate and attempted/completed full-rescore
counts are recorded separately, so a configured allowance does not imply that
an insertion improvement completed.

## A separate conservative step certificate

Let lambda be the sum of absolute Pauli coefficients and a=|T| lambda. Each
symmetric step agrees with exact evolution through degree two. At x=a/r<=1/2,
each degree-three-and-higher Taylor tail is bounded by
`x^3 / (6 (1-x/4))`. The product tail obeys the same multinomial norm bound.
Unitary telescoping therefore gives the deliberately loose bound

`delta = min(2, (2/3) a^3/r^2)`, for r>=max(1,ceil(2a)).

For normalized states, `1-|overlap| <= delta^2/2` and
`1-|overlap|^2 <= delta^2`. `conservative_steps` solves the corresponding rational
inequality by integer binary search. A step limit is an explicit refusal, never
a downward-clamped certificate. Zero time, zero coefficients and a single
commuting unit are exact. This bound is order independent and can be very loose;
it does not certify a cubic-ordering advantage or a useful compiled target range.

## Use

```python
from qhat.analysis.sparse_response_ordering import molecular_units, train_ordering

units = molecular_units(qubit_operator.terms, n_qubits)
result = train_ordering(units, n_qubits, hf_occupation, time=1.0,
                        rotation_budget=1000, cpu_seconds=2.0)
# Check result.status and result.score before using its estimated score.
ordered_units = [units[i] for i in result.order]
```

The larger molecular benchmark uses a shared 32-CPU-second construction
allowance, chosen from a separate resource-only seed timing before its main
order freeze. A two-second call can return an initializer with no completed
score on those inputs. Callers should choose and record a construction budget
appropriate to their input and inspect the returned status.

For an existing QHAT research-model JSON with integer coefficients and a common
denominator:

```sh
python -m qhat.analysis.train_ordering model.json --time 1 \
  --rotation-budget 1000 --cpu-seconds 2 --error-tolerance 1e-3 \
  --metric operator --output new-order.json
```

The command reports the response score and analytic certificate separately.
`fits_rotation_budget=false` means the certified target is not established under
that budget. It never uses exact-reference screening to select an order or step
count. The [independently validated larger-case comparison](../numerical_audit/scalable_ordering/evidence/README.md)
retains all 96 circuit outcomes. Cubic construction took 3.4–30 CPU seconds on
these four twelve-mode cases. The measured accuracy/cost tradeoffs are specific
to the panel; avoiding full-state arrays does not establish inexpensive training
on arbitrarily larger molecules.

## Opt-in Julia conditioned windows

`train_ordering(..., method='window', window_options=WindowOptions(...))` runs
bounded Julia insertion after a declared native sparse or coefficient-magnitude
parent. The default `method='sparse'` is unchanged. Import `WindowOptions` from
`qhat.analysis.window_ordering`. For example:

```python
from qhat.analysis.sparse_response_ordering import train_ordering
from qhat.analysis.window_ordering import WindowOptions
result = train_ordering(blocks, n_qubits, occupation, time=.37,
    rotation_budget=budget, cpu_seconds=12, method='window',
    window_options=WindowOptions(width=8, parent_method='sparse'))
```

The ordinary CLI exposes the same opt-in path:

```
python -m qhat.analysis.train_ordering model.json --time 1.6 \
  --rotation-budget 1000 --cpu-seconds 12 --method window \
  --window-width 8 --window-parent sparse --julia-executable /path/to/julia
```

Resolve Julia from that option, `QHAT_JULIA_EXECUTABLE`, or PATH. The backend uses
installed stdlib-only assets and never downloads Julia or installs packages.
Linux `/proc` and owned-PID `wait4` accounting, Python 3.11+ `tomllib`, and an
existing Julia 1.10+ binary are required for this optional stage. A missing or
failed backend retains the native parent and records `stats.window_status`.
`score=None` remains unscored. Invalid inputs or an infeasible minimum rotation
cap remain exceptions; admission interrupted before a valid initializer exists
also raises explicitly.

Initial support is one supplied basis occupation (arbitrary normalized global
phase allowed), real equal-flip internally commuting nonempty units, n<=20,
m<=1024, total original words<=8192, and bounded storage/work/rational bit sizes.
Empty units, zero coefficient words, mixed-flip/complex units and coherent
references retain a checked native parent with an explicit refusal. Terms are
never thinned or regrouped to gain admission. The unchanged current T=1/HF
control and supplied historical parents are not supported parent policies.
`starts` supplies seed permutations; each is validated and freshly scored, with
all cost charged to this invocation. It is not a free or trusted historical parent.

Widths are 2..9, default eight; a missing strictly interior window is reported as
`window_unavailable`. Windows keep global endpoints and original word counts, so
`r*(2U-l_first-l_last)+l_first <= rotation_budget` continues to use the actual cap.
Selection maximizes the number of nonzero local pair atoms, deterministically,
without inspecting physical errors. Each sequential window starts from the last
fully Taylor-checked global incumbent; the new table target aggregates real
conditioned columns before transport, preserving coherent exterior interference.
See [backend details](julia_window/README.md) for local cancellation refusals and
exact supplied-table screening. These are heuristic proposals, not optimality,
approximation-ratio or physical-error certificates.

`CPU_seconds` for the window method includes parent/setup, serialization, the
single owned child's startup/JIT and user+system CPU, final native verification
and cleanup. It does not use process-global `RUSAGE_CHILDREN`. Stats separate the
complete child receipt from informational phase timings (which are not added
again). Startup/JIT plus uninstrumented child work is reported as a residual,
not a separately measured JIT duration. Parent construction receives at most 60%
of the same CPU cap; a measured native full-score duration determines a reserved
final verification allowance. The remaining allocation includes Julia startup.
A parent-side CPU/wall watchdog kills and reaps the child. Default wall allowance
is max(1,4*CPU cap), configurable through `WindowOptions.wall_seconds` or
`--window-wall-seconds`; CPU and wall are distinct. Sampling and cleanup can have
a small observed overrun, which is reported without extending the declared cap.
There is no built-in hard memory containment; run under the documented external
memory guard when strict process-group limits are required.

A final authenticated checkpoint replaces the native parent only after the
existing full native scorer resolves an improvement greater than
`max(1e-30, abs(parent_J)*1e-12)` and agrees with Julia within its stated numerical
comparison tolerance. On a verification deadline/refusal or malformed checkpoint,
the native parent remains. An earlier atomic, fully checked child checkpoint may
still be independently verified after later worker failure. `status` describes
the usable native incumbent; `stats.window_status` describes the optional stage.
The result dataclass and CLI JSON fields are unchanged, and the optional
conservative certificate remains a separate order-independent calculation.
