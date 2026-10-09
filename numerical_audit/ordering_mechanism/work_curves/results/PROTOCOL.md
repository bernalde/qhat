# Frozen-order finite-step work curves

Declared after the development, factorial, and objective-ablation outcomes.
This tests whether the selected conserving diagonal-response rule's gains
survive outside r=100. It does not redesign the rule or count as held-out data.

Reuse exactly the already frozen eight refined candidate orders. For each
target and the same HF state and T=1, evaluate r=1,2,5,10,20,50 for candidate,
complete-class lexicographic, owned signed, JW signed, and JW magnitude.
All 240 schedules are fixed before the first evolution. Reuse existing r=100
rows only with their exact target/design provenance. Do not choose a new order
for any step count or observed error. Save all outcomes including reversals.

Work is retained Pauli-word rotations times r. These methods emit identical
word multisets at a fixed r. Compiled gates, parallel depth, and cancellation
across rotation boundaries have not been priced. Report the least *tested*
r meeting each E threshold (1e-4,1e-6,1e-8); no untested global optimum or
monotonicity is inferred from the grid. Use stable normalized phase-insensitive
E=1-|overlap| and retain all number/spin leakage. Independent circuit replay
uses the prior no-QHAT checker; no new tensor generation is claimed.
