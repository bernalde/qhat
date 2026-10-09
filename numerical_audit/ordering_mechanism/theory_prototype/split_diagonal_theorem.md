# Split-diagonal response ordering: precise theorem and limits

Exploratory mathematical note, 2026-10-09. These full theorems have not been
formalized in Lean or checked for literature novelty. A subsequent isolated
matrix-level bridge verifies the actual pair/BCH destination amplitudes and
orbital-sign covariance; see `lean/README.md` for its narrower scope. They describe the
algorithm tested in `/tmp/ordering-diagonal-refinement/results`, rather than
claiming that its development-set success is a universal theorem.

## Theorem 1: conserving refinement and exact reference objective

Let H be a real symmetric n-orbital occupation matrix commuting with total
particle number N. Partition its actual Pauli expansion by flip mask:

    H=D+sum_(x!=0) H_x,      H_x|v>=f_x(v)|v+x>.

Write D=sum_alpha D_alpha, where D_alpha=c_alpha Z_zalpha are the individual
nonzero diagonal Pauli terms. Identity can be omitted consistently. Let the
indexed unit list G consist of those D_alpha and the complete nonzero H_x.
Indices of diagonal units remain distinct even though all have flip mask zero.

Then:

1. Every G_i commutes with N, and every ordered product of its exponentials
   conserves N. Each nonzero H_x has pairwise commuting real-Hermitian Pauli
   descendants, and each diagonal unit is already one Pauli term. Every retained
   Pauli word is therefore emitted exactly once per Trotter step, irrespective
   of the refined-unit order.
2. Let x_i be the unit's mask (possibly zero), set

       b_ij(v)=f_j(v)f_i(v+x_j)-f_i(v)f_j(v+x_i),
       A_m(pi)=sum_(a<b, x_pi(a)+x_pi(b)=m) b_pi(a),pi(b)(v),
       B_pi=sum_(a<b)[G_pi(a),G_pi(b)], K_pi=-i B_pi/2.

   For a normalized determinant |v>, physical T>=0, and

       Q_D,pi=integral_0^T exp(isD)K_pi exp(-isD) ds,
       rho_m=T^2 sinc^2(T(E_(v+m)-E_v)/2), E_u=<u|D|u>,

   the exact diagonal-reference variance is

       J_D(pi)=Var_v(Q_D,pi)=1/4 sum_(m!=0) rho_m A_m(pi)^2.

3. The diagonal/off-diagonal portion of A_x equals

       f_x(v) sum_alpha s_alpha,x Delta_alpha,x,
       Delta_alpha,x=E_alpha(v+x)-E_alpha(v),

   where s_alpha,x is +1 if D_alpha precedes H_x in written matrix order and
   -1 otherwise. For a diagonal Pauli word,

       Delta_alpha,x=c_alpha(-1)^(zalpha dot v)
                         [(-1)^(zalpha dot x)-1].

4. Orbital-sign conjugation U_s|v>=(-1)^(s dot v)|v> leaves D, rho, J_D and
   every objective difference invariant. A deterministic descent method with
   sign-independent initial priority and label-only ties returns the same unit
   order under every such sign convention.

Proof. Since N is diagonal, each complete flip-mask projection of [H,N]=0
vanishes separately. Each diagonal D_alpha commutes with N individually. Real
symmetry permits only even-Y words in one mask class; the binary symplectic
pairing of any two such words on the same flip mask is zero. Their exponentials
therefore factor exactly, and refinement partitions the same Pauli word list.

Applying G_iG_j and G_jG_i to |v> gives the displayed b_ij. Summing pairs by their
common destination gives B_pi|v>=sum_m A_m|v+m>. Conjugation by D multiplies each
entry (u,v) by exp(is(E_u-E_v)); integration supplies its scalar Fourier factor.
Its modulus squared is rho. The correction -i/2 gives the factor 1/4. Variance
removes the destination m=0. In the present real setting all diagonal entries
of B are zero; in particular diagonal/diagonal commutators vanish despite
repeated zero labels. Different nonzero destinations are orthogonal.

The diagonal placement identity follows by the two energies in D_alpha H_x
versus H_x D_alpha. The Pauli specialization follows from its two occupation
characters. Finally U_s multiplies f_i by (-1)^(s dot x_i), so each b_ij and
A_m acquire the common character (-1)^(s dot m). Squaring removes it; diagonal
energies are unchanged. Deterministic selection consequently makes the same
decisions. The signs s_alpha,x are constrained by one common order; this theorem
does not permit choosing an independent sign for every alpha,x.

## Theorem 2: exact insertion update

For adjacent units i,j in written order, put m=x_i+x_j and b=b_ij(v). The swap
changes only A_m, by A_m'=A_m-2b, and hence

    J_D(pi')-J_D(pi)=rho_m(b^2-A_m b).

If a unit i is moved right across an ordered set R of units, define

    delta_m=-2 sum_(j in R, x_i+x_j=m) b_ij(v).

Then the insertion changes the objective by

    J_D(pi')-J_D(pi)
      =1/4 sum_m rho_m [2 A_m delta_m+delta_m^2].

Moving i left instead changes the sign of delta. These formulas hold with
repeated zero masks. Proof: each crossed pair reverses its sole contribution
to the pair sum; all other pairs keep their relative order. Expand each changed
square. Multiple crossed terms reaching one destination must be combined before
using the displayed quadratic expression, or applied sequentially to the
updated amplitude. Treating them as independent deltas against a stale A would
miss cross terms.

After pair preparation, an adjacent score/update is constant time. Scanning all
possible destinations of one insertion takes O(m) scalar updates, and scanning
all moved units takes O(m^2) such updates per sweep. This describes arithmetic
scoring, not necessarily Python list-copy or metadata costs in the prototype.
An accepted reduction of J_D is exact as a reference-objective statement; it is
not by itself a guarantee on actual Trotter error.

## Theorem 3: a sufficient condition for a correct actual-error ranking

Restrict to a finite invariant particle-number sector. Let H=D+V, let all units
be Hermitian, and let psi=|v>. Define the true transported correction Q_H,pi by
replacing D with H in Theorem 1. For T>=0,

    ||Q_H,pi-Q_D,pi|| <= T^2 ||V||||K_pi||=:e_pi,
    |sqrt(Var_v(Q_H,pi))-sqrt(J_D(pi))| <= e_pi.

Therefore, if

    sqrt(J_D(pi))+e_pi < max(0,sqrt(J_D(sigma))-e_sigma),

then the true leading variances obey the same strict ranking. For T>0 the
existing product-formula asymptotic theorem implies E_pi(r)<E_sigma(r) for all
sufficiently large r. The supplied inequality is sufficient, not necessary.
Upper bounds on the operator norms are admissible, but may be loose.

Proof. A Duhamel identity for the two conjugation groups gives

    ||Q_H-Q_D|| <= integral_0^T (T-u)||[V,exp(iuD)Kexp(-iuD)]|| du.

Unitary norm invariance and ||[V,A]||<=2||V||||A|| give e_pi. Apply the reverse
triangle inequality to the response vectors after orthogonal projection away
from psi. The strict gap between the resulting intervals implies the strict
variance gap. The established expansion

    E_pi(r)=T^2 Var_v(Q_H,pi)/(2r^2)+O(r^-3)

then implies eventual improvement. No finite-r threshold is supplied unless
the cubic remainders are bounded explicitly.

An exact sufficient condition for Q_H=Q_D is
[V,exp(iuD)K_pi exp(-iuD)]=0 for every u in [0,T]; the simpler joint conditions
[V,D]=[V,K_pi]=0 imply it. These are restrictive conditions, not assumptions
made about the molecular experiments.

**No claim is made that this guard certified the eight development wins.**
Those wins were checked against actual propagations and exact transported
diagnostics. A reduction in J_D without the guard remains a heuristic choice.

## Theorem 4: controlled weak-off-diagonal regime

Fix D_alpha, complete nonzero classes V_x, an order pi, T>=0, and a determinant v.
Let H_epsilon=D+epsilon V, where V=sum_x V_x and every off-diagonal unit is
scaled by real epsilon. Write

    K_pi(epsilon)=epsilon K_1,pi+epsilon^2 K_2,pi,

where K_1 contains the diagonal/off-diagonal commutators and K_2 the
off-diagonal/off-diagonal ones. Let

    u_pi=(I-|v><v|) integral_0^T exp(isD)K_1,pi exp(-isD)|v> ds,
    a_pi=||u_pi||,
    C_pi=T||K_2,pi||+T^2||V||(||K_1,pi||+||K_2,pi||).

For |epsilon|<=1,

    |sqrt(Var_v(Q_H_epsilon,pi))-|epsilon|a_pi|
       <= C_pi epsilon^2.

Thus, if a_pi<a_sigma, the actual leading-error ranking is correct whenever

    0<|epsilon|<min(1,(a_sigma-a_pi)/(C_pi+C_sigma)),

with the zero-denominator case interpreted directly by the preceding inequality.
In particular, if insertion achieves a_pi=0 but a_sigma>0, the true leading
variance is O(epsilon^4) for pi versus a_sigma^2 epsilon^2+O(|epsilon|^3) for
sigma. This is a genuine scaling improvement in that restricted regime.

Proof. The pair commutator expansion has no diagonal/diagonal term, giving the
displayed decomposition of K. The difference between Q_D(epsilon) and epsilon
times its K_1 response has norm at most T epsilon^2||K_2||. Theorem 3 bounds the
additional full-transport difference by

    T^2 |epsilon| ||V|| ||K_pi(epsilon)||
       <=T^2 epsilon^2||V||(||K_1||+||K_2||).

Projection cannot increase norm, so the two estimates and the reverse triangle
inequality give the result. Comparing the intervals proves the strict ranking;
setting a_pi=0 gives the stated variance scaling.

This theorem treats fixed orders and fixed T while changing physical coupling
epsilon. It does not replace Trotter-step asymptotics, assert that molecular
epsilon=1 is sufficiently small, or promise a global optimum. If optimization
cancels the first-order reference response, both the remaining response and the
transport-approximation error can be O(epsilon^2). Their squared norms are then
both O(epsilon^4), so the relative accuracy of J_D may become poor even though
the order has improved the true error substantially. This is the approximation
floor exposed by the phase-2 HCN diagnostics.

## Publication and novelty scope

The useful result would combine: a conserving refinement, a physically motivated
finite-time state-response objective, an efficient exact insertion identity, an
explicit regime of reliable ranking, and independently confirmed predictive
performance. None of those goals licenses claiming that interleaving or greedy
insertion was first introduced here.

[Tranter et al.](https://arxiv.org/html/1912.07555), Sec. III, alternate Z and XY
sets for molecular H2; Secs. V.2–V.3 study error-operator insertion optimization.
[Hastings et al.](https://arxiv.org/abs/1403.1539), Sec. 4/Eq. (6), group
HF-canceling hopping contributions, and Sec. 5/Eqs. (7)–(8) use diagonal-gap
perturbation for coefficient corrections. Our explicit mask/determinant formula
and its guard should be compared carefully against these and subsequent work.
At present no exhaustive novelty conclusion has been established.
