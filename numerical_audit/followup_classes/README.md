# Follow-up class, parent and duplication experiments

This opt-in research module addresses discussion Q2/Q3/Q6/Q12/Q13/Q16. It
leaves the historical benchmark paths and previous frozen campaigns unchanged.
Use the locally pinned H2O environment (`numerical_audit/h2o_recovery/requirements.txt`).
The entry points require QHAT's parent directory on `PYTHONPATH` (the independent
rechecker does not). Run from the QHAT repository root:

```sh
PYTHONPATH=.. python numerical_audit/followup_classes/run_class_comparison.py --output /tmp/new-class-comparison
PYTHONPATH=.. python numerical_audit/followup_classes/count_duplication.py --output /tmp/new-duplication-counts
PYTHONPATH=.. python numerical_audit/followup_classes/slide23_example.py --output /tmp/new-slide23-example
python numerical_audit/followup_classes/recheck_saved.py --run /tmp/new-class-comparison --output /tmp/new-independent-check.json
PYTHONPATH=.. MPLCONFIGDIR=/tmp/qhat-class-mpl python -m unittest discover -s numerical_audit/followup_classes -p 'test*.py' -v
```

Each output path must be new. The saved designs contain exact dyadic
coefficients and schedules; input tensors, source snapshots and state arrays
are retained. Graphs are exact computational integer certificates for these
binary coefficients, not Lean certificates. The protocol defines a shared
whole-parent-truncated Hermitian target; this is a separate comparison from
historical final-Pauli-pruned runs. Maximum target changes are reported.

The completed comparison has 106 evolutions: 72 primary (eight fixtures,
nine methods), 18 HCN convergence-anchor rows, and 16 parent formulas with
matched rotation budgets. An independent NumPy/SciPy implementation imports
no QHAT/OpenFermion/Numba and reproduces every final vector within 1.88e-14.
All eight reference matrices also match independent occupation action.

The graph-coloring rules do not improve accuracy uniformly. Original-parent
color ordering loses to signed parent ordering on all eight primary fixtures.
Complete-class coloring beats class lexicographic ordering on two of eight.
These counts describe one deterministic heuristic, HF states and the stated
small molecular fixtures; they are not a general result against coloring.

A crucial negative finding concerns the sign-retaining *suffix Pauli-l1*
objective. It is identical for all complete-class orders tested, and the
adjacent-swap refinement accepts no swaps. For fixed flip mask x, commutators
[H_x,H_y] for distinct y have distinct x xor y masks. Consequently they have
disjoint Pauli support and the l1 norm of their sum is additive. Summing this
bound over the order merely counts each unordered pair once. Thus this
particular cheap bound cannot select an improved complete-mask ordering.
Its greedy schedule can change actual errors despite leaving the bound
unchanged. The true sector operator norm is not this coefficient-l1 norm.

The counting audit has 104 rows: 11 original fixtures plus 15 reconstructed
unpruned H2O inputs, each at four parent thresholds. Exact binary counts and
thresholded final-Pauli counts are kept separate. At threshold 1e-12 all 26
four-distinct ratios are 4; exact binary residuals prevent literal cancellation
in two of those inputs. Exact synthetic checks verify 24 permutations become
three eight-member classes and physical quartet ratios are 4 (generic free
weights can instead give 3). Finite whole-H ratios do not identify an
asymptotic limit. Already consolidated schedules have ratio exactly 1.

The slide-23 example uses actual hopping JW expansions, not invented Pauli
factors. Its grouped and contiguous-parent circuit matrices differ only by
4e-16, while the interleaved circuit is different. This is a pedagogical
factorization check, not evidence of a molecular accuracy improvement.

An explicitly **exploratory** ten-row supplement was declared after the flat
suffix-bound diagnosis. It starts from `class_lex` and makes two adjacent-swap
passes, accepting only exact strict decreases of the **global signed BCH-l1**
sum. Its eight primary errors improve relative to class_lex, but the median
ratio is 0.996982: a small improvement, not a general accuracy advantage or a
finite-time ranking theorem. The supplement remains separate from the 106
original rows, and does not change their predeclared treatments.

```sh
PYTHONPATH=.. python numerical_audit/followup_classes/analyze_class_bound.py --run /tmp/new-class-comparison --output /tmp/new-suffix-bound-diagnostic.json
PYTHONPATH=.. python numerical_audit/followup_classes/bch_order_supplement.py --parent-run /tmp/new-class-comparison --output /tmp/new-bch-supplement
python numerical_audit/followup_classes/recheck_saved.py --run /tmp/new-bch-supplement --output /tmp/new-bch-independent-check.json
```

`bch_supplement_protocol.md` records scope and candidate choices before the
supplement ran. Exact adjacent updates are checked by recomputing the final
full BCH sum. All ten vectors also pass the independent rechecker. The
Frobenius bound idea remains a future theoretical direction, with no extra
experiment. JW comparators sort their float64 simulation coefficients;
canonical parent priorities use the exact dyadic parent weights. This
interface detail and all resulting schedules are explicitly saved.
