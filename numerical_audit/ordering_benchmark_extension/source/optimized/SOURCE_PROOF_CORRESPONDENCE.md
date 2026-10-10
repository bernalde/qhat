# Read-only source/proof correspondence review

Reviewed the two manuscript sources below against the recovered optimized benchmark implementation. No manuscript files were edited by this subtask.

- `paper/theory/conserving_ordering_method.tex`: `0195ef1df1b556875d34e32e482c186de223ab08c1bb780fb4c282e219dae731`
- `paper/theory/expanded_benchmark_methods.tex`: `03de926191fe923c38e70b1b8973f9e0314e29ebb0f6a960dc5d981c575cccd0`

The conserving refinement matches `common.normalized_case`: complete nonzero flip classes are retained and the zero-flip class is separated into individual diagonal words. Real same-flip Pauli words commute because their Y parities are even; the number-commutator matrix entries separate by flip destination. Thus completed units conserve number. The code does not assume that an individual intermediate off-diagonal Pauli rotation conserves number.

The response identity matches `optimize.pair_setup`, `amplitudes`, and `cost`. Integer coefficients have common denominator d, so pair amplitudes have denominator d² and the squared response has denominator d⁴. The specialized Pauli pair formula is checked exactly against the original general class algebra. The rho dictionary evaluates each distinct destination once while preserving first-encounter order and the original diagonal-gap arithmetic. All diagonal units contribute to D independent of their position in the order.

The code's chronological schedule is the reverse of written matrix order. This globally negates B and its destination amplitudes, preserving the variance objective, its insertion ranking, and its norm. The manuscript explicitly explains this reversal. The update subtracts twice the currently signed pair contribution; orientation is adjusted when the moved term crosses left rather than right. Exhaustive molecular-subproblem checks reconstruct the entire pair sum at every insertion position.

Both real-state variance and orbital-sign covariance arguments are sound in the stated scope. The conditional transport comparison follows from Duhamel and the commutator bound, then the reverse triangle inequality for the projected response. It is an asymptotic sufficient condition, and the manuscript correctly states that it did not certify the molecular rankings or a finite step threshold.

The Strang implementation is a product of exact commuting units. It fuses the full middle unit in a step and the full first unit across adjacent steps. For U words and endpoint unit lengths Lfirst, Llast, the count is r(2U-Lfirst-Llast)+Lfirst; all words still receive total duration1. Fused and unfused palindromic circuits are checked numerically. Requested budgets and actually emitted rotations are distinct fields.

The matched-search control receives the conserving JD search's count of examined insertion destinations, divided between the same two starting rules. It may finish a partial sweep. The text correctly distinguishes this from equal wall-clock work, excludes extra rescore/setup costs from the O(M²) scalar-update statement, and does not claim a formally verified floating-point bound.

During review, one wording issue was sent to the parent editor: the formalization-boundary paragraph initially referred to written proofs for the “insertion optimizer.” The actual displayed lemma proves insertion-update formulas, not implementation or optimizer correctness. The editor corrected this to “insertion-update formulas”; the manuscript hash recorded above includes that correction.

No mathematical error or additional implementation discrepancy was found in this correspondence review. This is a source-and-proof review, not a replacement for the independently maintained circuit replay and objective checks.
