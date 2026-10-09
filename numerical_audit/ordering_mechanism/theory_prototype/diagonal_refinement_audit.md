# Independent review of split-diagonal insertion experiment

Read-only review of `/tmp/ordering-diagonal-refinement/results`.
Source `run.py` SHA-256:
`58263730c8c35b2dedbad845574c3f4f4286af5c0328940d001af85ddd18ecd7`.

## Checks and outcome

No defect was found in the tested circuit, objective, insertion, or exact-response
calculations. `audit_diagonal_refinement.py` imports neither QHAT nor the run's
source. It directly constructs conventional X/Y/Z action on the saved occupation
basis, rebuilds exact target evolution, and replays every new circuit.

- All 20 approximate evolutions replayed independently.
- Maximum saved-state discrepancy: 6.49e-16.
- Maximum stable overlap-error discrepancy: 6.60e-19.
- Maximum independently rebuilt exact-reference discrepancy: 3.76e-15.
- Every planned order is a permutation of the actual refined units.
- Every retained Pauli word appears exactly once across that partition.
- Every actual refined-unit matrix commutes with particle number on the saved
  Pauli-invariant closure, including states outside the target number sector.
- All 24 reported J_H/J_D diagnostics were checked by independent dense block
  multiplication and adaptive conjugation quadrature, rather than the run's
  weighted-shift BCH and analytic sinc transport. Maximum discrepancy: 1.35e-17.

Machine-readable evidence is `diagonal_refinement_independent_checks.json`.

The code correctly sums all zero-mask units when constructing reference D,
distinguishes diagonal pieces by their stable original word order, and performs
crossed-pair insertion updates sequentially on the current destination amplitude.
The stable magnitude sort uses fixed generation order to resolve repeated-zero
mask ties; that order comes from sorted Pauli-word labels and is gauge invariant.
The phase-1 representation, chronological sign, and sinc checks still apply.

## Mechanistic interpretation

The true full-H leading variance decreases against both refined-unit magnitude
and the contiguous diagonal/class-lex baseline in all eight development cases.
The measured errors are in the leading regime: candidate E divided by J_H/20000
ranges from 0.99599 to 1.00448 at r=100,T=1. Consequently the gains are consistent
with genuine leading-error cancellation rather than metric artifacts or a
coincidental finite-step crossover.

The approximation is less quantitatively accurate after optimization:

| Case | Candidate J_D / class-lex J_D | Candidate J_H / class-lex J_H | Candidate J_H / J_D |
|---|---:|---:|---:|
| CH2O_s1.00 | 0.02350 | 0.03067 | 1.262 |
| CH2O_s1.25 | 0.09464 | 0.10095 | 1.048 |
| H2S_s1.00 | 0.001413 | 0.004832 | 3.374 |
| H2S_s1.25 | 0.002086 | 0.002481 | 1.171 |
| HCN_s1.00 | 0.01273 | 0.09051 | 6.460 |
| HCN_s1.25 | 0.009444 | 0.16218 | 15.687 |
| SiH4_s1.00 | 0.004304 | 0.004876 | 1.115 |
| SiH4_s1.25 | 0.01400 | 0.01475 | 1.043 |

J_D therefore remains a successful choice criterion on these comparisons but
underestimates the residual true coefficient, strongly for HCN. Reporting it as
an accurate numerical error estimator after cancellation would be incorrect.
The sufficient transport guard has not been shown to certify these eight wins.

This experiment changes both the partition and the optimizer relative to phase
1. It does not isolate a causal effect of diagonal splitting. Its frozen rule is
eligible for a held-out test; the eight development wins are not such a test.

## Provenance cleanup needed before durable release

The refined `design.json` records inherit several inactive old orders/graphs
whose indices referred to the former partition. The planned run orders and the
24 actual response diagnostics use the new partition correctly, so this does
not invalidate the tested runs. Before release, preserve the original frozen
files, then provide explicit metadata identifying which fields are inherited
historical context and which orders actually index the refined units. A generic
consumer must not silently apply an inherited graph or inactive class schedule
to the new block array.

The final theorem/proof formulation and approximation caveats are in
`split_diagonal_theorem.md`. These remain mathematical notes, not new Lean
verification claims or an established literature-novelty assessment.
