# Exploratory refinement, declared after complete-class candidate results

This is a separate development stage, not a predeclared held-out validation.
The preceding diagonal-response stage lowered errors within complete classes
but lost to JW signed in six of eight cases. Hypothesis: forcing all diagonal
terms into one contiguous block prevents useful commutator cancellation.

Keep the same eight exact targets, HF,T1,r100 and HCN r50/200 anchors. Refine
only the zero-mask class into its individual nonzero diagonal Pauli words;
keep every nonzero flip-mask class intact. This preserves number conservation
and emits exactly the same number of Pauli rotations. All words are emitted
once; no fractional steps or hidden additional exponentials.

Use the same diagonal-reference objective, with reference D equal to the SUM
of all zero-mask units. Distinguish these units by their full Pauli word labels.
Start from fixed mask/word order and descending unit coefficient-l1, label ties.
For each start perform at most ten best-improving single-unit insertion moves.
Evaluate each insertion via its crossed adjacent-pair updates without dynamics.
Only strict decreases >max(1e-30,1e-12 currentJ) are accepted; deterministic
scan order breaks within-search ties. Choose lower final objective across the
two starts; near ties use lexicographic unit-index order. No circuit error or
exact transported response is used in candidate selection.

Twenty new evolutions: candidate and refined-unit magnitude, eight primary
cases plus HCN convergence anchors. Freeze every candidate before any new
errors. Validate exact algebraic pair amplitudes, composite insertion changes,
five orbital-sign covariance controls, reversal symmetry, all exact response
diagnostics, and independent NumPy circuit replay. Existing baseline errors
refer to the identical target/state/time/rotation budget. Report every case,
including losses. This run changes both the partition and optimizer relative
to the previous stage; do not attribute all gain to diagonal splitting alone.
A later factorial control is needed if gains warrant a causal claim.
