# Symmetric product formulas under a rotation budget

`strang_schedule.py` is the maintained budget-to-pulses endpoint. It accepts
the sizes of a **nonempty** chronological list of internally commuting Pauli
groups. The caller supplies the groups, establishes internal commutation,
and flattens their terms in the same order. The planner does not inspect
operators or certify conservation. Counts are emitted Pauli rotations, not
compiled circuit depth.

```python
from qhat.analysis.strang_schedule import plan_strang_schedule

plan = plan_strang_schedule([3], budget=3)
if not plan.feasible:
    raise ValueError("No positive step count fits")
for term_index, numerator in plan.iter_pulses():
    duration = T * numerator / plan.duration_denominator
    # Apply exp(-1j * duration * terms[term_index]) to the current state.
```

For a single commuting group the entire target exponentiates exactly, so
one traversal at full time T suffices: `steps=1`, cost equal to the group's
number of terms. Each emitted numerator is 2 and the denominator is 2;
there is no second traversal or doubling of the evolution time.
An insufficient budget yields `feasible=False`, zero pulses and a `None`
denominator. Do not interpret that plan as a zero-step approximation.

For two or more groups, the first is outermost and the last central.
The largest feasible positive step count uses the fused cost
`r * (2*U - lengths[0] - lengths[-1]) + lengths[0]`.
The planner fuses the center and adjacent identical outer groups, preserving
within-group term order on both passes. Other possible commutation-based
fusions are outside this convention. Reversing the group list changes its
orientation. All durations are exact integer multiples of `T/(2*r)`.
Empty/scalar-only targets are outside this API; callers must handle their
identity evolution separately. Positive integer sizes and nonnegative
integer budgets are required. The pulse iterator is lazy.

The dependency-free CLI exposes the complete plan without propagating a state:

```bash
python analysis/strang_schedule.py --unit-lengths 3 --budget 3
python analysis/strang_schedule.py --unit-lengths 1 2 3 --budget 60
```

## Historical snapshots and the singleton correction

The frozen benchmark drivers in `numerical_audit/ordering_benchmark_extension/`
are historical replay sources. Their hashes and molecular results are
unchanged. They are **not redirected** to this new endpoint:

- `source/optimized/run.py` divides by zero when selecting singleton steps;
  its separate singleton pulse helper is already correct.
- `source/published/run.py` has the corresponding one-Pauli division by zero.
- `review_revision/source/strang_controls/run.py` has that selector problem
  and its unguarded pulse generator would accumulate twice the intended time.

New schedules should use the maintained planner above. Moving a historical
campaign to it requires a new source/design freeze and a separate output
location; never replace source hashes in an existing evidence manifest.
There is no claim that all older experiment frontends have been retrofitted.

## Tests

From the QHAT package's parent directory, using an environment with NumPy
and SciPy (the planner itself requires only the standard library):

```bash
python -m unittest qhat.analysis.tests.test_strang_schedule -v
```

Tests compare pulse matrices with independent dense matrix exponentials,
including a non-diagonal commuting singleton and noncommuting multi-group
palindromes, forward and reversed. Exhaustive small partitions use a literal
unmerged palindrome to check counts/maximality, plus budget boundary,
input-validation, duration and CLI regressions. These are numerical/software
checks, not additional Lean proofs or new molecular experiments.
