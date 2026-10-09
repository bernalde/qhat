# Initial-state sensitivity findings

The ordering ranking depends materially on the input state. At fixed `T=1`
and `r=100`, the fermionic signed order wins for the HF reference in
`CH2O_s1.25` and `HCN_s1.25`, but it loses for all five predeclared random
states in each of those two cases. This is a controlled state change: within
each comparison, Hamiltonian, state, time, Pauli set, and nominal exponential
count are identical across orderings.

Across the eight 10-qubit cases, the HF winners are JW signed in six cases and
fermionic signed in two. Across the 40 random-state comparisons, the winners
are JW signed in 38 and JW magnitude in two; fermionic signed wins none.
There are no winner ties at the predeclared absolute overlap-error tolerance
`1e-12`.

| Case | Fermionic / best JW error on HF | Range on five random states |
|---|---:|---:|
| CH2O_s1.25 | 0.529 | 1.434–1.911 |
| HCN_s1.25 | 0.910 | 1.228–1.775 |

The error is the stable full-state `1 - abs(overlap)`. These ratios do not
compare different normalizations or remove leakage. Fermionic signed retains
the fixed-spin sector to numerical precision (largest leakage `3.48e-36`),
while JW signed has leakage as large as `2.93e-8`; nevertheless JW signed often
has smaller total state error. Conservation alone therefore does not establish
a uniformly better approximation.

This is an exploratory sensitivity control, **not** a benchmark of molecular
ground states or a population-level performance claim. Five seeded complex
Haar states were sampled from each full fixed-`(N_alpha,N_beta)` sector. They
include components outside the HF-connected affine space and need not be
chemically representative low-energy states. The historical retained
Hamiltonians and preprocessing are held fixed; integral accuracy and the
scientific choice of active spaces are not assessed here.

## Reproducibility and numerical checks

Run directory: [results/initial-state-ergijt7d](results/initial-state-ergijt7d).
The directory contains the complete results, per-case designs, saved initial
states, immutable input/order/source hashes, and exact-reference crosschecks.
The five base seeds `2026100901` through `2026100905`, case-specific seed
streams, and state hashes were recorded before the tested evolution. The
reference Hamiltonian independently sums unit-Pauli sparse matrices with the
snapshot's hexadecimal coefficients, with no additional coefficient pruning.
Both reference and Trotter propagation use all 1024 amplitudes; no state is
postselected or normalized after sector projection.

All 144 ordering-state comparisons completed in 19.7 seconds on one thread.
The maximum norm defect was `3.87e-13`. Full-space `expm_multiply` agrees with
an independent sector eigendecomposition to at most `9.42e-14` in state norm.
The source hashes were checked again after the experiment.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /tmp/qhat-audit-venv/bin/python -m qhat.numerical_audit.initial_state_sensitivity --max-seconds 900
```

A separate regression extension validates the CLI's BK HF determinant for
register sizes 3, 5, and 6 by checking **every mapped mode occupation**, with
vacuum, partial, and full fillings. It does not introduce an arbitrary-excited-
determinant constructor; arbitrary encoded inputs remain supported through
state files. The focused standalone suite has 23 passing tests.
