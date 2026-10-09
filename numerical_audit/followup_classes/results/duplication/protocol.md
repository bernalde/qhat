# Parent, owned and complete-class controlled comparison

Version: exact-parent-common-target-class-comparison-v1.
Declared before evaluating the experiment errors. No learning model is used.

## Shared target

Use the eight existing 10-qubit molecular fixtures (CH2O, H2S, HCN, SiH4,
each scale 1 and 1.25), immutable tensor hashes from fixtures/manifest.json.
All molecular tensor coefficients are finite and exactly real. Treat their
binary values as exact dyadic numbers: exact CAR normalization, exact average
of each adjoint pair (the Hermitian part), discard entire canonical parents
with absolute weight <=1e-12, then exact JW summation. This explicitly defines
a new common target for this campaign. Report the Hermitian adjustment,
discarded parent l1 weight, and nonidentity Pauli l1 distance to raw and prior
final-Pauli-cutoff targets. Do not imply historical target identity.

There is no independent Pauli cutoff. Remove only exact zeros after exact
integer accumulation. Drop every scalar contribution from every schedule;
all reported state metrics are insensitive to this common global phase.
At the simulation interface, convert exact dyadics to float64 and report the
coefficient l1 discrepancy between each flattened schedule's float sum and
the reference target. No method receives a private target or repair.

## Frozen methods and graphs

1. Original-parent signed ascending canonical coefficient, canonical key ties.
2. Original-parent graph-color order.
3. First-owned blocks frozen once from method 1; signed owner order.
4. Same fixed owned blocks in graph-color order.
5. Complete 0/2/4 flip-mask classes, ascending integer flip mask.
6. Same classes in graph-color order.
7. Same classes using greedy sign-retaining suffix Pauli-l1 objective, followed
   by two left-to-right passes accepting only strict objective-decreasing
   adjacent swaps. The objective is sum_j ||coeff([G_j,sum_later G])||_1.
   It bounds the sector operator-norm objective but predicts no error ranking.
8. Individual JW signed ascending coefficients, lexicographic ties.
9. Individual JW magnitude descending, lexicographic ties.

Graphs use an edge exactly when the actual weighted block commutator is
nonzero in arbitrary-precision integer Pauli arithmetic. No numerical
commutation tolerance and no inherited-parent graph are used. Largest-degree
first greedy coloring breaks ties by initial block index. Colors are ordered
by their first-created integer label; inside colors use ascending block
index. Whole block factors remain intact. These computational zero checks
are not Lean certificates. Coloring alone changes no circuit at a matched
contiguous parent/class order and claims no accuracy advantage.

## Evolution and cost

Hartree-Fock input, physical time T=1, first-order formula. Primary r=100.
HCN scale 1 also uses r=50 and 200. Entries are chronological: later entries
multiply on the left. The common affine Pauli closure includes all parent
words, even those canceled from final H; no particle-sector projection is
applied. All methods target the same saved exact dyadic sum.

Consolidated methods use U rotations per step. Original parents use R and
retain duplicate words. Both parent methods additionally run at
max(1,floor(100 U/R)) steps, with actual counts and unused rotation budget
reported. This is a bounded integer budget comparison, not compiled gate
matching or exact budget equality. Eight primary cases yield 72 evolutions;
18 convergence-anchor and 16 matched-budget rows make 106 in total.

Metrics: normalized stable full-state E=1-|overlap|, infidelity, norm drift,
number leakage and spin-sector leakage, tie tolerance 1e-12. Independent
occupation-action matrix reconstruction and dense eigensolver references
check every small fixture. Save input tensors, exact coefficients, schedules,
graphs, source hashes/snapshots and all initial/exact/final state arrays.
No best-of-observed-errors schedule selection; the suffix rule uses only H.

## Counting and exact examples

A separate no-evolution audit measures all 11 original fixtures and all 15
reconstructed unpruned H2O tensors at parent thresholds 0,1e-14,1e-12,1e-10.
Save both exact-binary cancellation counts and documented thresholded Pauli
counts: floating tensors can violate analytic zero relations by roundoff.
The physical real-Coulomb quartet theorem predicts 4, not 3, in exact
arithmetic. Verify 24 permutations partition into three eight-member
Hermitian classes and each has eight actual JW words. Verify physical
three-parent and two-parent quartet ratios with exact coefficients.

A four-mode pedagogical demonstration uses h01,2h23,3h02,5h13 and actual
JW matrices. It distinguishes exact factorization within colors from outer
Trotter error; grouped and contiguous parent formulas must agree. It is no
molecular benchmark. No universal molecular asymptotic limit is inferred
from finite duplication ratios or qubit counts.
