# Response to the independent numerical project review

The reviewed sources are QHAT research merge `c28aba45644a817bc200c946523807569c2956c2`
and paper research merge `71608ab0f1709372b9c2a994d1ace3588e4b9ade`.
This response covers all six items in the completed numerical review sections.
The broader analysis/Julia and manuscript/proof review remains pending; this is
not a whole-project approval. The reviewer performed source inspection and
inspected the retained molecular scalar report; it performed no project execution.
The original reports are retained separately, with their scope and limitations.
The client reported `claude-opus-5-5`; high effort was requested.

### CLD-A2-001 — Done

Unsupported mapping names now raise the intended `NotImplementedError` from
`HamiltonianConfiguration.fermion_to_qubit_name()`. This also prevents
`f2q_tag()` and filename construction from silently using `None`. The defensive
dispatch branch reads the configured mapping string without calling it.
Every previously supported JW/BK alias keeps its name and filename tag.

The focused regression reproduced eight failures and fourteen passing alias
controls on the unchanged reviewed baseline. It exercises both configuration
accessors, filename construction and the defensive dispatch error path without
a chemistry calculation. See
[`test_hamiltonian_mapping_validation.py`](../../analysis/tests/test_hamiltonian_mapping_validation.py).

### CLD-A2-002 — Already satisfied

The missing-test premise does not hold at the reviewed source. Existing tests in
`analysis/tests/test_time_evolution_audit.py` check the 4-qubit JW/BK determinant,
mapped number operator and evolution agreement. Its parameterized
`test_bk_reference_matches_every_mode_occupation_for_non_power_of_two_sizes`
also checks each mapped occupation at 3, 5 and 6 qubits and several electron counts.

The pinned OpenFermion implementation constructs `bravyi_kitaev_code(n).encoder`
from the same `_encoder_bk(n)` used by the generator, so the alleged two-encoder
discrepancy has no source basis here. The existing regressions are reused rather
than adding a duplicate test.

### CLD-A2-003 — Done

The root README now links the time-evolution entry point and numerical study tree,
including the practical-accuracy, transfer and molecular-certificate studies.
It also links this response so the independent-checking boundaries are discoverable.

### CLD-B-001 — Done, differently

The current `sufficient_epsilon_domain(a, b, c, d)` is a sufficient condition for
the **ordered claim that the first slope is smaller than the second**. In the
fixed report, `a` is current placement and `b` is diagonal prefix. Its certified
gap is `b.lo - a.hi`, and the positive-coupling threshold is that gap divided by
`c + d`, subject to the stated strictness and coupling-at-most-one restriction.
`no_certified_slope_gap` means this forward sufficient condition is unresolved;
it does not classify the reverse ordering. The separate interval `winner()`
function already checks both directions.

This explicitly documents the helper's orientation instead of broadening the
frozen two-schedule reporter into a general winner-selection API. In the retained
LiH data the forward gap is positive and the algebraic threshold is approximately
0.03319. The physical leading guard remains unresolved, and the eight sampled weak
guards and 66 independent finite-step comparisons retain their original labels.
The existing reporter, report and archive bytes remain unchanged.

### CLD-B-002 — Done, differently

The practical verifier's mathematical paths independently reconstruct objective
values, primitive versus ideal action, DAG resources and target selections.
Its input authentication, dispatch naming and layout validation also reuse
producer bookkeeping helpers (`study.bind`, `dispatch`, `curve_name` and
`validate_record`). Those shared checks alone are not an independent proof of
the bookkeeping contract. This boundary applies to the original replay as well
as the current verifier.

A separate standard-library structural audit now checks all 81 retained
construction layouts and the complete cyclic dispatch plan from literal
Hamiltonian terms and the declared protocol. It imports no producer bookkeeping
helpers. Negative cases require it to reject a missing rotation and a duplicated
unit. See
[`test_practical_retained_structure.py`](../../analysis/tests/test_practical_retained_structure.py).
The regular analysis test lane exercises this additional audit. The historical
verifier and measured replay remain intact; no historical run is relabeled as
having executed this new check.

### CLD-B-003 — Already satisfied

`transfer_holdouts.study.constructor()` measures CPU through the complete
constructor return and overwrites the standalone charge with
`reference_charge + caller_CPU_seconds`. Its existing
`test_outer_training_charge_includes_final_constructor_return` deliberately gives
the inner constructor a 1-second charge and the outer caller an 8-second charge,
then requires an 8.5-second total including the 0.5-second reference charge.
This is the requested overhead boundary and regression protection.

## Validation and preserved scope

The unchanged baseline produced eight expected failures and fourteen passing
alias controls. After the fix, thirty mapping/BK/retained-structure cases passed,
along with the existing caller-boundary timing check and the strict-domain check
(thirty-two cases total). The targeted interpreter loaded the exact fix worktree.
Hosted CI and the broader independent review remain pending. No frozen study input, measured row, numerical report, archive,
producer/verifier snapshot, circuit, proof source, manuscript or PDF is changed
by these corrections. Scientific panel costs and earlier failed attempts retain
their original ledgers. Unknown quantum preparation cost, restricted certificate
domains, unequal consumed CPU and the software-versus-kernel distinction remain.

## Original report identities

These SHA-256 digests identify separately retained original reports; their bodies
and provider/account metadata are not included here.

| Completed section | Original report SHA-256 |
|---|---|
| A | `99b099810473195fdbeeea4ca112e3756379de242eab8dc5cbe920ff0a8ed830` |
| A2 | `9bd58a973745eeda1aaf866fb2f574c0c86396d56bc4925505509565d04a62a2` |
| B | `3cae6a6b29d244041859561b1d9e809aa1e931469b0516900925c8fc2a1213c3` |
