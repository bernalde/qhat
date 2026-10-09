# Exploratory diagonal-reference ordering, frozen before candidate errors

Goal: test a gauge-covariant cheap ordering objective motivated by transported
BCH response, while diagnosing its failure to include off-diagonal dynamics.
This is algorithm-development data, not a held-out confirmatory study.

Inputs: all eight exact-dyadic ten-qubit targets frozen in the prior class
comparison. HF, T=1, r=100. HCN_s1.00 adds r50/200 anchors. Existing original
results serve as paired baselines only on the same saved target.

Units: complete flip-mask classes, ascending integer-mask labels. Identity
omitted. Define b_xy(v)=f_y(v)f_x(v xor y)-f_x(v)f_y(v xor x), with H_x acting
as f_x(v)|v xor x>. Compute these exactly as integers with common denominator.
Let A_m(order) sum b_xy with the order's pair signs for destination m=x xor y.
Objective J_D = 1/4 sum_m rho_m A_m² / denominator^4,
rho_m=T² sinc²(T*(D(v xor m)-D(v))/2), sinc(z)=sin(z)/z.
This is diagonal-reference transported variance, not full-H variance.

Candidate: start independently from ascending mask and descending class
coefficient-l1 (mask ties); make at most ten left-to-right adjacent passes,
accepting J decrease greater than max(1e-30, 1e-12*currentJ); stop on no moves.
Keep the lower objective final order; objective ties within that tolerance
use lexicographic mask order. Selection never reads Trotter errors or exact
full-H transported variance. Save all orders before propagating any states.

New arms: class_diagonal_response and class_magnitude; 8x2 + 2x2 anchors=20
new approximate evolutions. Same target, closure, r and rotation count.
Preserve number leakage, use stable normalized phase-insensitive overlap E.
Independently recompute every adjacent update from full pair sums; verify
orbital-sign covariance on five fixed spatial patterns and reversal symmetry.
Validate b_xy with independent matrix commutator actions.

Diagnostic after schedule freeze: compute exact transported response using
small-space eigendecomposition for class_lex, class_magnitude and candidate.
Compare J_D, initial BCH variance and exact J_H, and measured E. These exact
quantities are diagnostic only and are excluded from the algorithm.

No claim of practical success, novelty, universally correct rankings, or
finite-r guarantee is made from surrogate reduction alone. Preserve negative
results. Independent NumPy circuit replay and source/input hashes required.
