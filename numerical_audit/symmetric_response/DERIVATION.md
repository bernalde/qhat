# Symmetric cubic-defect response

This source derivation precedes prospective physical errors. The implementation
and tests remain unexecuted until the resource guard permits their worker.

For real Hermitian conserving units G_i, use the actual chronological symmetric
step with half durations for the outer units and a full duration for the center:

`S_pi(h) = exp(-i h G_pi(0)/2) ... exp(-i h G_pi(m-1)) ... exp(-i h G_pi(0)/2)`.

The palindrome reads the same in written matrix order and chronological order.
For its action on the HF determinant, write the third-degree Taylor truncation
as `sum_(k=0..3) (-i h)^k q_k`. Initially q_0 is the determinant and other q_k
are zero. Applying one unit of relative duration tau maps the *old* arrays to

`q'_k = sum_(l=0..k) tau^l/l! G^l q_(k-l)`.

The lower coefficients are H|v> and H²|v>/2 for every order. Define
`A3_pi|v> = q_3 - H³|v>/6`. Then

`(S_pi(h)-exp(-i h H))|v> = i h³ A3_pi|v> + O(h⁴)`.

This does not use an exact-H time propagator or any validation error in search.
The sign of the effective Hamiltonian correction is opposite to A3; a common
sign does not affect the squared projected response below.

A different independent construction uses nested commutators. For outer A/2,
inner B, outer A/2, direct multiplication of the third-degree series gives

`A3(A,B) = -[A,[A,B]]/24 + [B,[B,A]]/12`.

Indeed, after subtracting (A+B)³/6, the A²B, ABA, BA² coefficients are
-1/24, 1/12, -1/24; AB², BAB, B²A have 1/12, -1/6, 1/12. Nesting the remaining
symmetric tail contributes its own cubic defect unchanged at this order.
With B_i equal to the sum of all units after i in the order,

`A3_pi = sum_i (-[G_i,[G_i,B_i]]/24 + [B_i,[B_i,G_i]]/12)`.

The producer evaluates the actual step's polynomial arrays; the verifier uses
this separate nested-commutator expression and literal Pauli matrices. Compare
the full cubic column, including its sign and parallel component, before
comparing projected scores. A matching squared norm alone would miss a sign
error. These are numerical algorithm checks, not new Lean theorems.

For the same diagonal reference D and HF determinant v, use

`J_D_2(pi,T) = sum_(w!=v) T² sinc²(T(E_w-E_v)/2) |A3_pi(w,v)|²`.

Here sinc(x)=sin(x)/x with value one at zero. Real coefficient numerators and
their common denominator define the supplied target. The producer computes
exact integer diagonal-gap numerators before conversion to binary64; the
independent verifier reconstructs gaps from literal matrices and uses math.sin.
This score approximates full-H transport by D; it does not certify actual
physical error or a global ordering optimum.

At fixed time and common ideal Pauli-rotation cap B, let r_pi(B) be the largest
affordable step count after exact central/boundary fusion. A symmetric cubic
defect gives a leading global state change proportional to T²/r² times its
finite-time transported response. The order-dependent squared score is thus
`J_D_2(pi,T)/r_pi(B)^4`; the common T⁴ factor does not change its ranking.
Optimizing J_D_2 alone would miss the discrete endpoint-dependent step count.

The hand-derived two-unit regression makes the distinction concrete. For
A=(2/7)Z and B=(3/7)X, A3(A,B)=-(2/343)X+(6/343)Z. Reversing their roles gives
(4/343)X-(3/343)Z. At a Z-basis determinant the projected squared response
changes by a factor of four, although the first-order response norm is invariant
under reversal. The small X fixture tests algebra only; it is not a physical
conserving validation target. A separate complete XX/YY hopping fixture checks
conservation and endpoint-dependent affordability.

The study will retain a formula-specific orientation-only control (the best of
three fixed seed orders under this new objective) separately from bounded
insertion optimization and from first-order-trained symmetric reuse. All
physical comparisons, losses, floors, exact circuits and construction costs are
required before any accuracy/work conclusion. No novelty, universal mechanism,
finite-error guarantee, compiled resource or hardware claim follows here.
