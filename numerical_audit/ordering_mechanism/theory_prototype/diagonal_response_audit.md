# Independent review of the first diagonal-response experiment

Reviewed source: `/tmp/ordering-diagonal-response/run.py`, SHA-256
`9f42562b5b6b1bd9f25bf479fca05132b2f9a238e3be0cdedcf82eecec7152c2`.
Read-only review; no experiment sources or results were modified.

## Implementation assessment

No mathematical implementation defect was found.

- `fval` evaluates conventional even-Y Pauli words correctly: `i^number_of_Y`
  gives the factor `(-1)^(number_of_Y/2)`, multiplied by occupation Z parity.
- The coefficient convention numbers orbital q by the little bit `1<<q`;
  simulation occupation indices use `1<<(n-1-q)`. `reversebits` converts these
  conventions consistently. Independent direct-Pauli reconstruction confirms it.
- Pair amplitudes have the written commutator sign `H_x H_y-H_y H_x`.
  Emission uses chronological order, so its actual leading Hermitian correction
  is `+i B/2`. The diagnostic uses precisely this sign. The opposite convention
  in the theory note is its written-product convention and gives the same
  variance, because a global sign cancels in a variance.
- The diagonal gap filter includes the factor 1/4 outside the squared BCH
  amplitude and uses the correct NumPy sinc convention. T=1 is intentionally
  hard-coded throughout this experiment and matches the frozen protocol.
- Every accepted adjacent update is independently recomputed from all pair
  amplitudes; the objective and its tolerance are sign invariant.
- All eight stored mask arrays are strictly ascending. Therefore index-based
  lexicographic tie resolution is exactly the prescribed mask-label tie rule.
- Candidate orders are saved before new circuit errors or full-H response
  diagnostics are evaluated. This is a legitimate exploratory comparison, not
  a held-out evaluation.
- The exact-H transport uses the spectral gap `lambda_row-lambda_column`, as
  required by `exp(+isH) K exp(-isH)`. Removing the parallel-to-initial component
  computes variance, not merely the response norm.
- The leading prediction `J_H/20000` correctly equals `T^2 J_H/(2r^2)` at
  T=1,r=100.

## Independent checks actually executed

`audit_diagonal_response.py` imports neither QHAT nor the experiment's source.
It constructs the action of each conventional tensor-Pauli word directly on the
saved occupation indices, using local X/Y/Z rules. It then replays all 20 new
circuits with explicit cosine/sine Pauli rotations, rebuilds every exact reference
from the saved target, and recomputes the stable phase-insensitive error.

- Maximum approximate-state difference: 4.34e-16.
- Maximum metric difference: 6.46e-19.
- Maximum independently rebuilt exact-reference difference: 3.76e-15.

For CH2O_s1.00 and HCN_s1.25, it also constructs actual dense class matrices,
forms the whole BCH matrix by matrix multiplication against successive suffix
sums, and integrates transported response with adaptive vector quadrature.
These cases cover a small improvement and the largest diagonal-approximation
discrepancy. All three schedules per case were checked:

- Six independently integrated response diagnostics.
- Maximum J_H/J_D difference from the reported values: 6.94e-18.

The machine-readable record is `diagonal_response_independent_checks.json`.
These are independent small-space propagations and matrix/transport checks;
they are not Lean proofs or molecular integral-generation audits.

## What the mechanism diagnostics establish

The candidate decreases both J_D and exact J_H against both `class_lex` and
`class_magnitude` in all eight development cases. Candidate ratios J_H/J_D range
from 0.96245 to 1.10909. Measured candidate E divided by the leading prediction
ranges from 0.999602 to 1.000238. Thus the observed within-class improvement is
consistent with the intended leading-response cancellation, and neither a
wrong overlap normalization nor a pre-asymptotic ranking artifact explains it.

The diagonal approximation is not exact. For HCN_s1.25 the candidate J_D/lex
ratio is 0.39986 while J_H/lex is 0.48551: it overstates the improvement, although
it predicts its direction. The data support this surrogate on these development
comparisons; they do not prove that surrogate descent always improves J_H.

Measured candidate performance remains limited: eight wins against class-lex
with median error ratio 0.95137, seven wins against signed owned blocks, eight
wins against JW magnitude, but only two wins against JW signed. There is a
useful mechanism to investigate, not yet a strong general ordering result.

## Next justified refinement: distribute diagonal words

Let the zero-mask class be D=sum_alpha D_alpha, with each D_alpha diagonal,
and retain every nonzero complete-mask class H_x intact. For an initial
determinant v, define

    Delta_alpha,x = <v+x|D_alpha|v+x> - <v|D_alpha|v>.

Then direct action gives

    [D_alpha,H_x]|v> = Delta_alpha,x f_x(v)|v+x>.

If `s_alpha,x=+1` when D_alpha precedes H_x in written order and `-1` otherwise,
the total diagonal/off-diagonal contribution to the BCH amplitude at destination
x is

    f_x(v) sum_alpha s_alpha,x Delta_alpha,x.

A single unsplit D limits this contribution to `+/- Delta_x f_x(v)`. Splitting
D into individual Z words permits additional cancellations among those signed
gaps and with commutator paths between nonzero masks. For D_alpha=c_alpha Z_z,

    Delta_alpha,x = c_alpha (-1)^(z dot v) [(-1)^(z dot x)-1].

These formulas follow by evaluating the two diagonal energies before and after
the weighted shift. They do not need an approximate HF eigenstate assumption.
Different x cannot choose their signs independently: all signs must arise from
one common order, so each diagonal term's placement imposes a consistent cut
through the other factors.

Every diagonal Pauli word separately commutes with particle number, and all
nonzero-mask classes retain their conservation and exact commuting-Pauli
factorization. Each original Pauli rotation is still emitted once per step.
Thus this refinement preserves the conservation guarantee and primitive rotation
count, while enlarging the available orderings. Its single-mask pair-amplitude
formula remains valid with repeated zero masks; zero/zero commutators vanish.

Implementation requirements:

- Build the reference diagonal energy from **all** zero-mask pieces; selecting
  the first zero-mask block is no longer correct.
- Give diagonal pieces fixed distinct labels including their phase mask for
  tie-breaking. A flip-mask-only label is no longer injective.
- If one insertion crosses several blocks with the same destination, aggregate
  or sequentially update their amplitude changes; do not add independently
  evaluated quadratic deltas relative to the same stale original amplitude.
- Keep insertion budgets fixed before evaluating errors. Comparisons should
  separate the effects of decomposition refinement and more powerful search.

This is a justified next candidate but the basic ideas have precedents.
[Tranter et al.](https://arxiv.org/html/1912.07555), Sec. III, explicitly alternate
the Z and XY sets for H2; Secs. V.2–V.3 derive an error-operator insertion update
and an insertion heuristic. [Hastings et al.](https://arxiv.org/abs/1403.1539),
Sec. 4/Eq. (6), exploit HF cancellation by grouping bare and spectator-assisted
hopping, while putting the diagonal terms first. Their Sec. 5/Eqs. (7)–(8)
uses diagonal energy-gap perturbation to modify coefficients. That differs from
our fixed-H ordering objective. Interleaving, insertion, and energy-gap motivation
should not individually be advertised as new. The proposed combination still
needs a closer prior-art comparison and independent held-out performance.
