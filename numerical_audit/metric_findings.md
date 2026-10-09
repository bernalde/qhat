# State-metric audit

The historical state-overlap objective is meaningful across orderings when the
Hamiltonian, initial state, physical time, and budget agree. The historical
`compare_states` implementation nevertheless had numerical robustness defects.
These defects alone do not establish that any reported molecular comparison is
incorrect. The static BCH cancellation diagnostic has a separate predictive
limitation; replacing its normalization does not make it a fixed-time fidelity
coefficient.

## Provenance

- Requested branch: `improve_ordering`, commit
  `5bc6ded98fe6ab52182c24f03139832b5f6349e6`.
- Updated upstream main: `d3eea1ea458908f75d317794a80400a07896c2a8`.
- Actual imported molecular study: `L-sweep`, commit
  `89c701652c4c64918825e279e84a7bea5924822a`. Original file hashes are in
  `source_provenance.json`.
- `analysis/error_analysis.py` was absent on the requested branch: its defects
  are inherited from updated main. `time_evolution.py` was byte-identical on
  requested branch and updated checkout before these edits.

`metrics_before_after.json` compares the original `compare_states` function,
extracted with Python AST from the untouched L-sweep worktree, against the
corrected local function. It is a synthetic regression witness, not molecular
benchmark data.

## Reproduced defects and changes

1. An outside-sector amplitude of `1e-10` previously produced zero error and
   zero leakage through subtraction from an overlap rounded to one. The
   corrected implementation retains `one_minus_overlap = 5e-21`, infidelity
   `1e-20`, and leakage `1e-20`.
2. A state with norm `1.000001` previously appeared perfect after clipping.
   Significant norm defects now raise; both norms and defects are recorded.
3. A normalized state with 25% leakage gives the same physical result before
   and after. The original sector overlap already counted leakage: the fix
   must not renormalize the approximate sector projection.
4. The shared metric helper evaluates the manuscript error as half the squared
   globally phase-aligned distance, then evaluates infidelity as `E*(2-E)`.
   It retains phase-sensitive distance as a different observable. The legacy
   molecular caller embeds the reference into the full space and computes
   leakage from outside-sector amplitudes.
5. Baseline and ablation consume the stable error directly. All four consumer
   CSV schemas carry `normalized_full_state_v2`, norms, and defects. Old-header
   CSV files cannot be resumed into this schema.
6. The standalone CLI now handles config return values and paths, normalizes
   its public fidelity input, transports reference occupation bits for BK,
   and rejects malformed states and invalid step counts. It prints overlap
   error, infidelity, and norm defects separately.
7. General `error_analysis` preserves its old phase-sensitive absolute and
   relative fields while adding normalized phase-independent metrics. Its
   energy diagnostic now matches eigenstates bijectively, labels the principal
   phase branch/energy alias period, and marks relative errors at zero energy
   undefined. Eigenvector-based energy matching remains ambiguous at degeneracy;
   it is not a universal ordering-quality score.

Phase-sensitive state distances can be appropriate, especially where controlled
unitaries make global phases observable. They are not the manuscript's
phase-insensitive state-overlap error. Relative energy errors also depend on
the chosen energy origin.

## Validation

All commands used the isolated environment and one BLAS/OpenMP/Numba thread.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /tmp/qhat-audit-venv/bin/python -m pytest analysis/tests/test_operators.py analysis/tests/test_error_analysis.py analysis/tests/test_state_metrics_audit.py analysis/tests/test_time_evolution_audit.py analysis/tests/test_legacy_state_comparison_audit.py -q
```

Result: **98 passed**. Synthetic four-qubit end-to-end baseline and ablation
runs completed 5/5 and 2/2 schedules, respectively. A BK CLI run selected the
correct two-electron determinant (basis index 8 instead of JW index 12).

The frozen audit's new reference machinery was independently checked using
literal Kronecker-product Pauli matrices and dense exponentials:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /tmp/qhat-audit-venv/bin/python -m pytest numerical_audit/test_reference_tools.py -q
```

Result: **6 passed**. Tests include odd-Y phases, complex initial states,
affine-space intertwining, retained leakage outside a conserved-number sector,
the chronological commutator sign, dense finite differences, defect columns,
an independent quadrature of the transported defect, nonzero global-phase
expectation subtraction, and convergence of the fixed-time prediction as the
step count increases. The dynamic prediction tests use invariant reference
sectors, as required by its restricted exact propagator. No live runner or
reference-helper source was changed during this validation.
