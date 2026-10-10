# Expanded ordering benchmark results

Complete frozen comparison.

4312 planned rows; 4186 completed evolutions; 126 infeasible budgets. All 14 targets are fixed; four development molecular families at two geometries and three prospective families at two coefficient multipliers are not 14 independent molecular samples.

Candidate: complete nonzero flip classes with individually movable diagonal Paulis, determinant-response JD insertion. Comparisons below use budget 100U, T=1 and the same HF state and target. W/T/L refer to the candidate; absolute tie tolerance 1e-12. Ratios omit comparator denominators<=1e-12.

| Comparator | Dev W/T/L | Prospective W/T/L | Overall W/T/L | Median E candidate/comparator | Ratio cases |
|---|---:|---:|---:|---:|---:|
| BCH1 full-l1 insertion adaptation | 8/0/0 | 6/0/0 | 14/0/0 | 0.00488242 | 14 |
| Complete class BCH l1 | 8/0/0 | 6/0/0 | 14/0/0 | 0.0506589 | 14 |
| Complete class coloring | 8/0/0 | 6/0/0 | 14/0/0 | 0.0549145 | 14 |
| Complete class lexicographic | 8/0/0 | 6/0/0 | 14/0/0 | 0.0506051 | 14 |
| Complete class magnitude | 8/0/0 | 6/0/0 | 14/0/0 | 0.0510558 | 14 |
| Pauli response JD | 2/0/6 | 2/0/4 | 4/0/10 | 1.16995 | 14 |
| Pauli response JD, matched evaluations | 8/0/0 | 6/0/0 | 14/0/0 | 0.18284 | 14 |
| Conserving unweighted J0 | 2/6/0 | 2/4/0 | 4/10/0 | 1 | 14 |
| Complete class suffix l1 | 8/0/0 | 6/0/0 | 14/0/0 | 0.0510491 | 14 |
| Hastings consolidated, lexical completion | 8/0/0 | 6/0/0 | 14/0/0 | 0.0505891 | 14 |
| Hastings literal, lexical completion | 8/0/0 | 6/0/0 | 14/0/0 | 0.00453589 | 14 |
| Pauli coloring | 8/0/0 | 6/0/0 | 14/0/0 | 0.00625353 | 14 |
| Pauli lexicographic | 8/0/0 | 6/0/0 | 14/0/0 | 0.00463998 | 14 |
| Pauli magnitude | 8/0/0 | 6/0/0 | 14/0/0 | 0.0244594 | 14 |
| Pauli signed | 8/0/0 | 6/0/0 | 14/0/0 | 0.046486 | 14 |
| Owned coloring | 8/0/0 | 6/0/0 | 14/0/0 | 0.00355382 | 14 |
| Owned block l1 | 8/0/0 | 6/0/0 | 14/0/0 | 0.00700093 | 14 |
| Owned parent magnitude | 8/0/0 | 6/0/0 | 14/0/0 | 0.0418612 | 14 |
| Owned reverse signed | 8/0/0 | 6/0/0 | 14/0/0 | 0.0402481 | 14 |
| Owned signed | 8/0/0 | 6/0/0 | 14/0/0 | 0.0402438 | 14 |
| Parent coloring | 8/0/0 | 6/0/0 | 14/0/0 | 6.35947e-05 | 14 |
| Parent magnitude | 8/0/0 | 6/0/0 | 14/0/0 | 0.00103443 | 14 |
| Parent signed | 8/0/0 | 6/0/0 | 14/0/0 | 0.00142699 | 14 |
| Refined class magnitude | 8/0/0 | 6/0/0 | 14/0/0 | 0.051812 | 14 |
| Strang, complete class lexicographic | 0/0/8 | 0/0/6 | 0/0/14 | 4118.13 | 5 |
| Strang, conserving JD order | 0/0/8 | 0/0/6 | 0/0/14 | 9427.69 | 2 |
| Strang, Pauli magnitude | 0/0/8 | 0/2/4 | 0/2/12 | 1258.58 | 7 |
| Strang, Pauli signed | 0/0/8 | 4/0/2 | 4/0/10 | 184.783 | 9 |
| Tomesh-inspired beam32 | 8/0/0 | 6/0/0 | 14/0/0 | 0.00288443 | 14 |
| Tranter commutator | 8/0/0 | 6/0/0 | 14/0/0 | 0.00351977 | 14 |
| Tranter depleteGroups | 8/0/0 | 6/0/0 | 14/0/0 | 0.00290596 | 14 |
| Tranter equaliseGroups | 8/0/0 | 6/0/0 | 14/0/0 | 0.00197646 | 14 |
| Tranter reverseCommutator | 8/0/0 | 6/0/0 | 14/0/0 | 0.0027717 | 14 |
| Tranter V2-l1 order, first-order use | 8/0/0 | 6/0/0 | 14/0/0 | 0.00389944 | 14 |
| Strang, Tranter V2-l1 insertion | 0/0/8 | 2/0/4 | 2/0/12 | 783.18 | 8 |

## Symmetric use of the frozen conserving order

At budget 10U (less affected by the error floor), W/T/L now favor the symmetric conserving-order formula. This is an empirical transfer of a first-order surrogate, not a second-order objective certificate.

| Symmetric comparator | W/T/L | Median E conserving/comparator | Ratio cases |
|---|---:|---:|---:|
| Strang, complete class lexicographic | 14/0/0 | 0.0989086 | 14 |
| Strang, Pauli magnitude | 13/0/1 | 0.102486 | 14 |
| Strang, Pauli signed | 14/0/0 | 0.0630087 | 14 |
| Strang, Tranter V2-l1 insertion | 14/0/0 | 0.0662779 | 14 |

Fixed-random distributions are archived separately for every family/case; no order was selected using its measured error. Strang is a different formula order; its budgets account for exact commuting boundary fusion. Evaluation-matched Pauli search does not match classical wall time. Full methods, input hashes, threshold costs and leakage are in summary.json.
