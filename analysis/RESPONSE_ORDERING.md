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

Time and coefficients use consistent units with hbar=1. Input coefficients and
time are evaluated in floating arithmetic for training. Conservative certificates
instead interpret supplied rational values exactly, including the exact dyadics
represented by floats. They concern ideal unit exponentials; numerical rotation,
synthesis, SCF, integral, basis and hardware errors are separate.

## Response objective

For a permutation pi, the symmetric step is the forward half-step sequence,
the final full unit, and the reverse half-step sequence. The implementation
propagates only Taylor degrees zero through three on reachable bit strings.
Subtracting H^3/6 gives the cubic defect column C_pi. The common physical factor
i multiplying that column is omitted consistently.

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
count. The independently validated larger-case benchmark lives under
`numerical_audit/scalable_ordering`; measured outcomes and limits are reported
there after execution, not inferred from the implementation.
