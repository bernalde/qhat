# Mechanism-driven ordering: exploratory theory note

Prepared 2026-10-09. This note does not change the reviewed manuscript or
formalization. Its new statements have mathematical proofs below. A later
isolated Lean bridge verifies actual pair/BCH destination amplitudes and
orbital-sign covariance; see `lean/README.md` for its scope. The integral,
Frobenius, and approximation-guard results below remain prose proofs.
No claim of literature novelty is made for them.

## 1. Scope and notation

There are n spin orbitals and occupation states v in F_2^n. Addition of masks
means XOR. Let P_eta project onto occupations of Hamming weight eta. A complete
flip-mask class has the weighted-shift form

    H_x |v> = f_x(v) |v+x>.

The current molecular construction gives real symmetric H_x, distinct mask
labels x, and [H_x,N]=0. These imply f_x(v+x)=f_x(v), and f_x(v)=0 whenever
|v+x| != |v|. The scalar identity may be retained or removed consistently; it
changes neither commutators nor any energy gap below. Use a fixed indexed
decomposition: optimizing first ownership and optimizing an order of existing
blocks are different operations.

For a **written matrix-product** order pi define

    B_pi = sum_(a<b) [H_pi(a),H_pi(b)],
    K_pi = -i B_pi/2,
    Q_H,pi(T) = integral_0^T exp(i s H) K_pi exp(-i s H) ds,
    H = sum_x H_x.

These conventions match OrderingCovariance.lean. Chronological gate order is
the reverse written order. Reversing all factors changes B and Q by a minus
sign and therefore does not change their variance objectives. It does not
permit careless sign conventions for a single local swap.

The established endpoint in OrderingAnalysis.lean is

    E_pi(r) = T^2 Var_psi(Q_H,pi(T))/(2 r^2) + O(r^-3),

for normalized psi, Hermitian blocks, fixed T>0, and the paper's phase-insensitive
overlap error E=1-|<psi_exact,psi_approx>|. It is not an operator-norm error or
ground-energy error statement.

## 2. Exact covariance under orbital signs

Choose s in F_2^n, let chi_s(x)=(-1)^(s dot x), and let
U_s|v>=chi_s(v)|v>. In JW this is a product of Z gates. Direct evaluation gives

    U_s a_p U_s^* = (-1)^s_p a_p,
    U_s a_p^dagger U_s^* = (-1)^s_p a_p^dagger,
    U_s W(x,z) U_s^* = chi_s(x) W(x,z),
    U_s H_x U_s^* = chi_s(x) H_x,

where W(x,z)|v>=(-1)^(z dot v)|v+x>. For the last equality, the matrix element
acquires chi_s(v+x)chi_s(v)=chi_s(x). Occupied spectators contribute each orbital
twice and hence no extra sign. The same character multiplies every Pauli
coefficient of a canonical parent with active mask x.

Consequences, with exact or absolute-value-based cutoffs:

1. The zero pattern, flip labels, and Pauli supports do not change.
2. A fixed parent order has the same first owners before and after transformation.
   Every owned block transforms by U_s conjugation, including when its final
   coefficients include contributions from several parents.
3. Sorting parents by absolute coefficient, or complete classes by any absolute
   coefficient p-norm, gives the same indexed order if all ties are broken by
   fixed labels independent of signs. Consequently the constructed product
   formulas are conjugate and have identical corresponding errors.
4. Signed sorting need not give the same order. For parents sharing x, changing
   chi_s(x) from +1 to -1 reverses ascending signed priorities (away from ties).
   If their Pauli supports overlap this can change ownership itself. This is
   more than reordering a fixed set of owned matrices.
5. All conjugation-invariant commutator norms, sector Frobenius norms, and
   coefficient absolute norms are invariant. A full BCH coefficient-l1 objective
   for a fixed class order is invariant, even though it can depend on order.
6. U_s commutes with N. An HF determinant changes only by an irrelevant global
   sign. A multi-determinant initial state must be transformed by U_s; keeping its
   coefficient vector fixed would generally compare different physical states.

This is covariance under **orbital signs**, not arbitrary orbital rotations.
General rotations mix masks, change the chosen decomposition, and need not
preserve its product-formula error. For restricted spatial-orbital rephasing,
the corresponding alpha and beta spin orbitals must receive the same sign.

These observations are exact algebraic consequences, not assertions that a
signed heuristic is illegitimate. A representation-dependent heuristic needs
a fixed convention and a robustness assessment; one can also deliberately
optimize a gauge, but that would be a separate algorithm with a search cost.

## 3. Sector Frobenius suffix bound

Assume each H_x preserves the eta sector. For an operator A write
||A||_F,eta^2=Tr(P_eta A^* A P_eta), which equals the Frobenius norm squared of
its restriction when A preserves the sector. Set

    w_xy = ||[H_x,H_y]||_F,eta^2 >= 0.

### Proposition F1: exact suffix Pythagoras identity

For any finite set R of distinct partner masks,

    ||[H_x,sum_(y in R) H_y]||_F,eta^2 = sum_(y in R) w_xy.

Proof. Products H_x H_y connect only occupations differing by x+y, so their
commutator has that mask. For fixed x, distinct y give distinct x+y. Thus two
different commutators have disjoint matrix-entry supports, also after restricting
to the eta sector. Their Hilbert-Schmidt inner product is zero. Expand the norm
squared of the sum. This proof works for complex matrices and requires sector
preservation only to interpret restriction as a closed dynamics.

The existing suffix-commutator operator bound therefore gives

    Lambda_eta(pi) <= Phi_eta(pi)
      := sum_a sqrt(sum_(b>a) w_pi(a),pi(b)),

    ||S_pi(T/r)^r - exp(-iTH)||_op,eta <= T^2 Phi_eta(pi)/(2r).

This is genuinely order dependent, unlike the complete-mask suffix Pauli-l1
budget proved invariant in MaskCommutatorBudget.lean. It is an upper bound on
operator error, not a ranking theorem for actual state error. In particular,
Phi can distinguish an order from its reverse even though their leading overlap
variances are equal. A normalized Frobenius norm may be used to compare objective
values inside one fixed model, but one must restore sqrt(binomial(n,eta)) in an
operator-norm bound. It is wrong simply to substitute normalized Frobenius norm
for operator norm.

### Proposition F2: exact adjacent-swap rule

For adjacent x,y followed by R put

    a=sum_(z in R) w_xz, b=sum_(z in R) w_yz, c=w_xy.

Their contribution is sqrt(a+c)+sqrt(b), versus sqrt(b+c)+sqrt(a) after swapping.
If c>0, the original x,y order has smaller cost iff a>b. If c=0 both orders tie.
Proof: g_c(t)=sqrt(t+c)-sqrt(t)=c/(sqrt(t+c)+sqrt(t)) is strictly decreasing for
c>0 and t>=0. The original minus swapped cost is g_c(a)-g_c(b). Earlier suffixes
contain the same two matrices and have the same pair-weight sum; later ones are
unchanged.

Descending current weighted degree, deleting one vertex at a time and resolving
ties by mask label, gives an adjacent-local optimum. At the time x is selected,
its degree c+a is at least the prospective next vertex's c+b. Maintain degrees
in O(m^2) arithmetic after precomputing all pair weights. This is not a theorem
of global optimality. A four-vertex nonnegative-weight counterexample is saved
in frobenius_greedy_counterexample.json: greedy objective 6.6055512755 versus
optimal 5.9907047849. These abstract weights have not been asserted realizable
by a particular molecular Hamiltonian.

### Proposition F3: exact non-exponential pair-weight formula

Write a commutator of mask m in the signed XZ basis as

    C = sum_z d_z W(m,z).

For C preserving the eta sector,

    ||C||_F,eta^2
       = sum_(z,z') conjugate(d_z) d_z' K_eta(popcount(z+z')),
    K_eta(k) = sum_j (-1)^j binomial(k,j) binomial(n-k,eta-j).

Proof. C|v>=sum_z d_z(-1)^(z dot v)|v+m>. Square its amplitude and sum over all
|v|=eta. For z+z' of weight k, partition eta-occupations by the number j of their
ones inside that support. There are binomial(k,j)binomial(n-k,eta-j) such states,
each with sign (-1)^j. The expression includes cancellations between Pauli words;
they are not orthogonal after fixed-particle-number restriction.

No 2^n-sized state vector is required. Precompute K_eta(k) for k=0,...,n, then
evaluate a commutator with L stored signed-XZ coefficients in O(L^2) dot products.
One can combine equal z+z' values before evaluation. Use exact integers/rationals
where available; cancellation can make naive float evaluation numerically poor.
For conventional Hermitian Pauli coefficients, first convert their unit phases
to signed-XZ coefficients. Omitting that conversion gives an incorrect formula
for odd-Y, anti-Hermitian commutators.

For a molecular one-/two-body Hamiltonian the input Pauli expansion has polynomial
size in n, even though JW strings are long. Thus sparse Pauli commutator formation
and this trace formula are polynomial operations; the naive degree can still be
large and this is not a claim that the all-pairs implementation is practical at
arbitrary n. Fixed-spin sectors can use a product of two Krawtchouk factors, one
for each spin subset, if every block preserves both spin counts. For a general
selected affine closure, the simple Hamming-weight formula need not apply.

## 4. Recommended mechanism: diagonal-reference transported variance

The Frobenius objective is a useful baseline. The following criterion follows
the actual leading overlap coefficient more closely while remaining classically
computable without exact Hamiltonian evolution.

Let D=H_0 be the full occupation-diagonal class and V=H-D. Fix an HF determinant
|v>, physical T>=0, and complete-class order pi. Define

    Q_D,pi(T)=integral_0^T exp(i s D) K_pi exp(-i s D) ds,
    J_D(pi)=Var_v(Q_D,pi(T)),
    E_u=<u|D|u>, Delta_m=E_(v+m)-E_v,
    F_T(Delta)=integral_0^T exp(i s Delta) ds,
    rho_m=|F_T(Delta_m)|^2
          =T^2 sinc^2(T Delta_m/2),

where sinc(z)=sin(z)/z with sinc(0)=1. In NumPy, whose sinc argument is measured
in pi units, use T**2*np.sinc(T*Delta/(2*np.pi))**2. At exact Delta=0 use T^2;
do not divide by zero or silently remove a degenerate transition.

### Proposition D1: explicit determinant response

Set

    b_xy(v)=f_y(v) f_x(v+y)-f_x(v) f_y(v+x),
    A_m(pi)=sum_(a<b, pi(a)+pi(b)=m) b_pi(a),pi(b)(v).

Then

    [H_x,H_y]|v>=b_xy(v)|v+x+y>,
    B_pi|v>=sum_m A_m(pi)|v+m>,
    J_D(pi)=1/4 sum_(m!=0) rho_m |A_m(pi)|^2.

Proof. Apply the two weighted shifts in their written product order and subtract.
All pairs with the same XOR destination contribute to the same basis amplitude.
Conjugation by D multiplies a matrix entry (u,v) by exp(i s(E_u-E_v)). Integrate,
multiply by -i/2, and take the squared norm after removing the |v> component.
For distinct complete masks x!=y there is no m=0 pair. Real coefficients give
real b and A; the displayed absolute-square also accommodates complex values.

This is an exact formula for the **diagonal-reference** transported variance.
It is not an exact formula for Var_v(Q_H) unless an additional hypothesis makes
the two transports equal. Calling it the exact Trotter error would be wrong.

The interference is explicit: different pairs with identical XOR masks can
cancel before squaring, while different destinations remain orthogonal. The
gap filter retains the physical time T. It suppresses transitions undergoing
rapid phase cancellation under D, which a coefficient-l1 norm ignores.

### Proposition D2: O(1) adjacent updates and exact covariance

If adjacent written factors x,y are swapped, only m=x+y changes:

    A_m' = A_m-2b_xy,
    J_D' - J_D = rho_m (b_xy^2-A_m b_xy)       [real case].

For complex coefficients replace the right side by
rho_m (|b_xy|^2-Re(conjugate(A_m)b_xy)). This follows by expanding the single
changed absolute square, including its factor 1/4. Thus a cached candidate swap
can be scored and applied in O(1) arithmetic. A deterministic fixed-pass local
search need not call a propagator or use measured Trotter errors.

Under orbital signs, b_xy'=chi_s(x+y)b_xy and A_m'=chi_s(m)A_m, while D and rho_m
are unchanged. Therefore J_D and every adjacent delta are invariant. Starting
from a gauge-invariant mask order, using a sign-invariant improvement tolerance,
and resolving ties by fixed mask labels produces the same output order for every
orbital sign choice. This should be tested numerically; mathematical invariance
does not excuse inconsistent floating accumulation or sign-dependent tie rules.

Reversing the entire order sends A to -A and leaves J_D invariant, matching the
established reversal symmetry of the true leading coefficient.

### Proposition D3: explicit conservative ranking guard

Assume H,D,K_pi are Hermitian on a finite invariant sector and T>=0. Then

    ||Q_H,pi(T)-Q_D,pi(T)||
       <= integral_0^T (T-u)
            ||[V, exp(iuD) K_pi exp(-iuD)]|| du
       <= T^2 ||V|| ||K_pi|| =: epsilon_pi.

Proof. The superoperator Duhamel identity for the two conjugation groups writes
their difference on K as the integral of i[V,exp(iuD)Kexp(-iuD)], sandwiched by
unitaries from H. Unitary invariance of operator norm removes those sandwiches.
Integrate in s from 0 to T, interchange the triangular integration domain, and
use ||[V,A]||<=2||V||||A||. The integral of 2(T-u) is T^2.

For a normalized initial vector and the orthogonal projection onto its complement,
sqrt(Var_psi(Q))=||(I-|psi><psi|)Qpsi||. The reverse triangle inequality therefore
gives

    |sqrt(Var_v(Q_H,pi))-sqrt(J_D(pi))| <= epsilon_pi.

Consequently, if two orders satisfy

    sqrt(J_D(pi))+epsilon_pi
       < max(0, sqrt(J_D(sigma))-epsilon_sigma),

then their true leading variances have the corresponding strict ranking, and the
existing eventual-overlap theorem gives E_pi(r)<E_sigma(r) for all sufficiently
large r. No finite-r threshold follows without remainder constants.

All norms in this guard may be replaced by rigorously computed upper bounds,
such as Pauli coefficient-l1 bounds, at the price of looseness. A failed guard
is inconclusive, not evidence against an empirical improvement. The guard may
be far too loose at T=1 in chemistry and must not be presented as routinely
certifying the heuristic. A short-time alternative is

    ||Q_H,pi(T)-T K_pi|| <= T^2 ||[H,K_pi]||/2,

obtained by differentiating the exact conjugation. It gives another explicit
but potentially loose ranking guard for the untransported variance.

### Implementation and experiment implications

If class x stores k_x Pauli coefficients and M=sum_x k_x, evaluate f_x(v+y) for
all class pairs in O(m M) Pauli-character operations, plus diagonal energies for
the distinct XOR destinations. Cache f_x(v), b_xy, rho_m, and A_m. No state vector
of dimension 2^n or binomial(n,eta) is needed. The arithmetic and memory costs
of an explicit m-by-m b table are still material; stream or sparsify it if needed.

Recommended bounded first experiment:

- Compare a frozen gauge-invariant starting order against exactly specified
  fixed-pass adjacent descent of J_D, at the same Hamiltonian and T.
- Record construction time, objective before/after, accepted swaps, b/rho hashes,
  final mask order, actual overlap error, and particle-number leakage.
- Independently compare b_xy with direct small occupation matrices and compare
  J_D with a dense diagonal-transport calculation on tiny cases.
- Check zero gaps, zero commutators, reversal invariance, and all/specified orbital
  sign transformations. All priorities and tolerances must transform consistently.
- At a few small cases, compute actual Q_H independently to separate failure of
  the diagonal approximation from failure of the finite-r leading approximation.
- Use existing cases only for diagnosis. Freeze the rule before new hold-out
  molecular geometries/active spaces; do not claim eight familiar cases are a
  predictive validation set.
- Test successful large H2O cases on a non-HF state separately. An HF-specific
  objective does not promise improvement for correlated states or operator norm.

## 5. Formalization path and boundaries

Inspected read-only:

- MaskCommutatorBudget.lean already proves matrix mask multiplication,
  commutator support, actual coefficient disjointness, and suffix-l1 invariance.
- OrderingCovariance.lean proves the exact adjacent BCH/transport update and the
  state covariance criterion, then connects a strict variance gap to actual
  eventual overlap improvement.
- OrderingAnalysis.lean proves the reversal symmetry and actual overlap
  asymptotic comparisons.
- EncodingTransport.lean gives representation-transport results, but its existing
  permutation transport does not itself supply a diagonal sign-conjugation API.

Proposed minimal new theorem chain, if experiments justify development:

1. Diagonal occupation-sign unitary and weighted-shift covariance; deterministic
   indexed schedules then inherit conjugation covariance.
2. Sector Hilbert-Schmidt orthogonality from distinct matrix masks, followed by
   the suffix Frobenius formula and its adjacent objective identity.
3. Actual fixed-Hamming-weight Walsh trace/Krawtchouk identity, with unit-phase
   bridge for conventional Pauli words.
4. Weighted-shift commutator action on a determinant and diagonal-transport entry
   formula, giving J_D with its 1/4 and physical-time filter.
5. Adjacent J_D update/covariance; then the two-reference Duhamel norm bound and
   strict ranking corollary using the existing eventual-overlap endpoint.

Do not formalize an isolated scalar surrogate and then claim the physical
representation bridge is automatic. The most valuable endpoints are an actual
matrix identity and an explicit approximation guard. None of these would alone
establish useful molecular error reductions or publication-level novelty.

## 6. Limited prior-art check

This was a targeted check, not an exhaustive novelty search.

- Tranter et al., *Ordering of Trotterization: Impact on Errors in Quantum
  Simulation of Electronic Structure*, already study commutation-guided orderings
  and ordering through an error-operator norm. Therefore graph ordering or
  minimizing a BCH-derived surrogate is not by itself new.
  [Primary paper](https://arxiv.org/abs/1912.07555).
- Babbush et al., *Chemical Basis of Trotter-Suzuki Errors in Quantum Chemistry
  Simulation*, already motivate classically tractable approximate states for
  estimating chemistry-specific Trotter error. Therefore using HF information
  or a cheap state-dependent error estimate is not by itself new.
  [Primary paper](https://arxiv.org/abs/1410.8159).
- Kronenberger, Erakovic, and Reiher, *Trotter Error and Orbital Transformations
  in Quantum Phase Estimation*, examine orbital transformations and find that
  simple descriptors need not correlate well with Trotter error. Their general
  orbital rotations are broader than the exact sign covariance here. This is
  directly relevant caution for interpreting any cheap proposed objective.
  [Primary preprint](https://arxiv.org/abs/2602.18913);
  [Published paper](https://doi.org/10.1080/00268976.2026.2681062).

Potential contribution, if supported after a proper prior-art comparison: an
explicit complete-mask/determinant cancellation objective with physical-time
diagonal transport, exact O(1) adjacent scoring, representation covariance,
an approximation guard, and useful held-out computational performance. At
present this is a research direction, not an established novelty claim.

## 7. Independent prototype checks performed

`pair_frobenius.py` implements sparse signed-XZ commutator formation and the exact
fixed-weight character trace formula without importing QHAT. The exported weight
functions do not enumerate occupation states. `check_prototypes.py` independently
constructs small real number-preserving occupation matrices and uses enumeration
only as a validation oracle.

Executed with `/tmp/qhat-h2o-venv/bin/python check_prototypes.py`:

- 1,584 pair/sector identities on 2-, 3-, and 4-orbital matrices passed with exact
  rational arithmetic.
- An independent dense 4-orbital, 2-electron calculation checked the diagonal
  transport objective, every adjacent swap, and all 16 orbital-sign masks. The
  largest residual was 1.11e-16.
- At T=0.13, exact eigendecomposition gave a full-versus-diagonal transported
  operator difference 1.1070786850, below the conservative bound 1.5677330767.

These are representation and formula checks, not molecular benchmarks, formal
Lean proofs, or evidence that the surrogate predicts actual chemical errors.
Exact gap resonances can make a diagonal-reference contribution vanish; testing
nearby times is useful before interpreting an exceptionally small objective as
a robust physical mechanism. The summary is saved in `prototype_checks.json`.
