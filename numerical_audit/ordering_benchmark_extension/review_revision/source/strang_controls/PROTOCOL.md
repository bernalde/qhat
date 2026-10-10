# Frozen Strang transfer and orientation controls

This review-driven control reuses the 14 exact Hamiltonians, HF initial states,
T=1, and previously selected JD orders in the immutable extension archive
SHA256 d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43.
No search is rerun and no ordering is selected using new propagation errors.

Six arms use the conserving JD, larger-search individual-Pauli JD, and
insertion-position-budget-matched Pauli JD orders in forward and reversed
unit order. The first chronological unit is outermost in each Strang
palindrome. Reversal reverses units, retaining the original word order inside
each conserving commuting unit. Each unit is emitted internally in that order
on both halves; internal order does not alter its exact commuting exponential.
The 84 forward conserving states are reused with byte-level provenance.
The other five arms give 420 new evolutions. Every one of 14 designs and all
420 schedules/budgets are frozen, together with runner/verifier source hashes,
before the first new evolution.

Budgets are b U for b=2,5,10,20,50,100, with U the unique nonzero consolidated
Pauli count. Each arm uses the largest integer r fitting the budget under
C(r)=r(2U-L_first-L_last)+L_first after exact declared-unit fusion. For
individual Pauli units this is 2r(U-1)+1. The same retained target, full affine
Pauli-invariant closure and reference state are used, with no postselection.
E=1-|<exact,approximate>| is evaluated by the normalized, phase-aligned
residual; absolute ties use 1e-12 and ratios require comparator E>1e-12.

Independent verification derives Pauli action directly in physical occupation
coordinates and replays an unmerged forward/reverse Pauli palindrome, with no
QHAT import. All new states are replayed. Existing baseline states are also
replayed when practical. Target sums, closure, reference vectors, norm,
particle/spin leakage and metrics are checked independently. Production
imports explicitly bind the qhat package search path to an extracted original
archive and check every imported qhat module against its manifest digest.

Threshold comparisons use epsilon=1e-4,1e-6,1e-8 and report lower/same/higher
first tested budget level, with both success/failure counts. Exact rotation
counts are also retained but within-level differences do not establish a
resolved minimum-work ranking. A uniform maximal commuting prefix/suffix
fusion cost sensitivity uses the existing states at their same r; it does not
create outcomes at extra steps that a different compiler might permit.

Forward and reverse orientations are compared separately. Any retrospective
minimum over them is an oracle, not the performance of a deployable selected
schedule. HF-state finite-T outcomes do not establish operator-norm, arbitrary-
state or compiled hardware performance. Search is not reoptimized for Strang.

All numerical processes cap address space at 6 GiB before numerical imports,
set numerical thread counts to one, and run serially. Timing intervals distinguish
propagation, independent replay and prior search; concurrent system load precludes
isolated speedup claims. Completion and all analysis require all planned rows.
