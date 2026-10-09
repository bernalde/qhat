# Exploratory signed BCH-l1 supplement

Declared after the original 106-run class comparison and its exact diagnosis
that the suffix coefficient-l1 bound is order independent for complete flip
classes. This supplement is explicitly exploratory; its 10 rows are separate
from the original campaign and are not included in its primary summaries.
No learning model or error-guided schedule search is used.

Reuse all eight exact-dyadic class decompositions and common targets already
saved by the original class comparison. Start from class_lex (ascending
integer flip mask). Let B(pi)=sum_(i before j) [H_i,H_j]. Its Pauli-l1 norm
retains signs and cancellation across different block pairs, unlike the
suffix-l1 bound. Carry out exactly two left-to-right passes of adjacent swaps.
For adjacent A,B the new sum is B(pi)-2[A,B]. Accept a swap iff the *exact
integer coefficient-l1 sum* decreases strictly. Ties are rejected. Freeze the
complete resulting order and both BCH sums before evaluating any state error.

This is a leading-defect/operator-norm surrogate. A lower value does not
certify a smaller finite-time overlap error, smaller transported variance,
or a better finite-r error bound. No local-optimality or global-optimality
claim is made after two passes. The initial order, pass count and objective
will not be tuned after observing errors.

Evaluate the resulting schedule for all eight cases at T=1,r=100 (eight
rows), and HCN scale1 at r=50 and200 (two anchors). Use the same HF initial
vectors, reference target and complete Pauli closure as before, with no
projection, same stable metric and rotation counts. Save exact order, initial
and final BCH coefficients/objectives, accepted swaps, source/parent-artifact
hashes and final vectors. Independently replay all ten using NumPy/SciPy
without QHAT, OpenFermion or Numba. Compare descriptively with class_lex and
JW signed, retaining the four-family/eight-related-fixture limitation.

The sector Frobenius suffix bound is recorded as a future theoretical idea
only: distinct commutator flip masks are Hilbert--Schmidt orthogonal, so a
sector Frobenius expression can retain square-root ordering information.
There is no Frobenius-order experiment in this supplement.
