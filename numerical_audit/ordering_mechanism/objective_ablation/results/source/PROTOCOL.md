# Objective and unit-partition ablation, fixed before new errors

2026-10-09. Exploratory algorithm-development control after the phase-2
diagonal-refinement candidate. It is not held-out validation.

Use exactly the same eight ten-qubit targets in /tmp/followup-class-comparison,
the same HF initial states, T=1, first-order r=100, complete invariant affine
closures, retained coefficients and scalar convention. Sixteen evolutions:

1. `class_refined_J0`: the phase-2 refined partition (individual diagonal
   Pauli words; complete nonzero flip-mask classes), but all rho_m=1. This
   is untransported HF BCH variance at T=1.
2. `class_pauli_JD`: individual Pauli words as units, sorted first by
   (flip mask, full Pauli label), using the original finite-time diagonal
   reference weights rho_m. D remains the sum of all diagonal units.
   This comparison does not conserve number after every unit: individual
   transverse Pauli words may leak. Keep that leakage, with no postselection.

Both use the identical phase-2 optimization: fixed-label and descending
unit-coefficient-l1 starts; at most ten best-improving insertion moves per
start, tolerance max(1e-30,1e-12*currentJ), deterministic scan/tie rules.
Select the lower final surrogate with the same near-tie rule. No hyperparameter
tuning, error-based selection, choice of best gauge, or exact-H response is
permitted. Save all sixteen orders before any new propagation. Existing
phase-2 JD orders and results stay frozen and are paired comparators only.

For individual Pauli units duplicate nonzero masks are permitted. Verify
exactly that every nonzero pair amplitude has destination different from HF;
real Pauli commutators with equal flip masks vanish, so the HF expectation
is zero. If this fails, stop rather than incorrectly using a squared norm as
variance. Every unit is real and Hermitian, each Pauli appears once, all target
coefficients are preserved, and all comparisons have the same rotations.

No new exact full-H transport diagnostics. References are reused after
target/state/hash checks. Save active schemas, arrays, source/input hashes
and all results including losses. Independently replay all sixteen circuits
using the no-QHAT saved-design checker. Run sequentially with one thread and
a 6 GiB resident-memory cap. The phase-2 search/helper source is copied
verbatim to phase2_reference.py and hashed before the run.
