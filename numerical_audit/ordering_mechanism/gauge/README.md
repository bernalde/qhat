# Orbital-sign robustness diagnostic

The signed-first-owned ordering is **not robust to physically equivalent
orbital-sign choices** in this controlled H2O panel. The physical retained
Hamiltonians are exactly gauge conjugate, and matched fixed schedules give
identical final states after conjugation. Rebuilding the signed schedules
changes both their order and first ownership, and changes their error.

The frozen protocol is in `PROTOCOL.md`. There are 72 re-sorted evolutions and
60 explicit fixed-order controls, all HF, T=1, r=100, historical-policy tensors,
with six fixed spatial-sign patterns per case and the same sign for alpha/beta.
Cases were selected before this diagnostic to include previous successes and
a failure, plus a STO-6G control. They are not an independent test set for a
new algorithm. Patterns are fixed probes, not samples from a gauge distribution.

| Case | Owned E, original gauge | Owned E range over six gauges | JW-magnitude E, invariant | Owned vs magnitude wins / losses | Owned vs signed wins / losses |
|---|---:|---:|---:|---:|---:|
| STO-6G, 10q | 1.45060e-6 | 1.44926e-6–1.45809e-6 | 1.36637e-6 | 0 / 6 | 4 / 2 |
| HGBS-5, 12q | 5.79639e-10 | 3.32050e-10–8.10813e-10 | 4.44086e-10 | 1 / 5 | 3 / 3 |
| HGBS-5, 14q | 9.10971e-10 | 8.90318e-10–3.15319e-9 | 1.78798e-9 | 2 / 4 | 6 / 0 |
| HGBS-5, 16q | 2.10498e-9 | 2.10498e-9–2.39062e-7 | 3.00744e-7 | 6 / 0 | 1 / 5 |

There are no ties at absolute tolerance 1e-12. These counts are six gauge
representations of the **same physical problem** per row, not six independent
molecular benchmarks. At 16q all five nontrivial fixed patterns increase owned
error to 2.25055e-7–2.39062e-7. Meanwhile signed-JW improves from 2.28555e-7 to
1.26039e-8–1.95102e-8. Thus the 143-fold owned-versus-magnitude advantage in the
original gauge falls to only about 1.26–1.34-fold, and the owned-versus-signed
ranking reverses. At 14q the owned-versus-magnitude ranking reverses for four
of the five nontrivial gauges. At 12q a previously losing case can become a win.

This does not invalidate the earlier numbers or the symmetry/conservation
theorems. It shows that their accuracy advantage cannot be attributed only
to the physical molecule and active space: the sign convention is part of
the algorithm specification. Selecting the best gauge after seeing these
errors would be a different, error-informed algorithm and is not proposed.
A gauge-equivariant ordering criterion or a fixed canonical gauge warrants
study, followed by independent validation. Canonicalizing signs alone does
not imply improved accuracy.

The first-owner partition itself changes, not merely the order of fixed
blocks. In the 16q case, 216–272 Pauli words change owner under the five
nontrivial patterns. Complete flip-mask membership and all within-class
coefficient magnitudes remain exactly invariant in all 24 case/gauge pairs.
This is a positive structural control, not a propagation test for a new
complete-class ordering.

## Verification

- Every transformed one-/two-body tensor is checked against direct sign
  multiplication. Final retained Pauli keys and coefficients obey exact
  signed covariance, with no tolerance. The cutoff therefore changes no
  paired retained target in this experiment.
- Sparse `H_g = D H D` is checked exactly on the invariant affine closure.
  No number/spin projection is performed. Closure dimensions are 64, 256,
  1024, 4096 for the four cases.
- All 60 explicitly propagated fixed-order controls, and all magnitude-JW
  schedules, have **zero** state difference from the conjugated original
  saved trajectories. Magnitude-JW errors are bitwise invariant.
- The 10q reference is independently diagonalized in the full 1024-dimensional
  space using direct local Pauli actions: state difference 2.752e-15. All 18
  ten-qubit re-sorted trajectories are independently propagated with NumPy
  Pauli products, without QHAT's compiled-mask functions.
- `recheck_gauge.py` imports no QHAT code. Under `python -O` it checked all
  132 saved states and all 115 producer artifact hashes, exact tensor and
  coefficient covariance, complete-class magnitude controls, closure, norms,
  metrics, and leakage. Maximum metric difference is 2.711e-20. This is a
  saved-state audit, not a second independent large-space propagation.
- Eighteen deliberate wrong-sign coefficient controls are rejected; six
  case/gauge combinations induce no retained sign change and need no such
  control. There are no Python `assert` validation gates.
- Producer ran with at most two threads and an 8 GiB resident-memory cgroup
  cap. Peak RSS was 911312 KiB (about 890 MiB).

`results/` contains frozen producer source, designs, signs, schedules, blocks,
transformed tensors, full arrays on each recorded closure, results and hashes.
`independent_recheck.json` contains detailed ranges and independent check
outputs. Prior references were reused after checking the original saved
Hamiltonian, design, archive and array hashes; they were not newly propagated
for every gauge because exact unitary covariance provides the reference.

Reproduce with the pinned reviewed QHAT checkout at
`0be7d4857537668eec957dabaa8438ea5c0d33f2`, its pinned H2O environment, and
the preceding H2O data/reference archives:

```sh
systemd-run --user --scope -q -p MemoryMax=8G -p MemorySwapMax=0 \
  /tmp/qhat-h2o-venv/bin/python run_gauge.py \
  --qhat-root /path/to/pinned/qhat --output /tmp/new-gauge-results
/tmp/qhat-h2o-venv/bin/python -O recheck_gauge.py \
  --results /tmp/new-gauge-results --output /tmp/new-gauge-recheck.json
```

Output directories/files must be new. `--qhat-root` verifies the imported
checkout and its commit. The run cannot silently use a different live QHAT
checkout. No reviewed repository was edited by this experiment.
