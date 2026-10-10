# Historical baseline word order

The fixed campaign reproduces the same ideal products as the historical baseline methods. There are 64 literal historical flat-sequence hash matches and 32 cases where words were permuted inside exactly commuting blocks.

The earlier `check_frozen.py` compared all 96 schedules with the sequences reconstructed from the serialized input documents. The historical development runner serialized each block in sorted order but executed its original dictionary order. Consequently, serialized-design equality alone did not establish literal equality to that historical execution. The corrected `archived_schedule_check.json` reports these two checks separately.

The 32 differences are `parent_signed`, `parent_color`, `owned_signed`, and `owned_color` for each of the eight development cases. `check_historical_internal_order.py` reconstructs the original dictionaries from the canonical fermionic components, reproduces all 32 historical flat hashes, and verifies that the new sequence has exactly the same ordered blocks and coefficients. It independently checks every internal Pauli pair by the binary symplectic commutation criterion.

For an internally commuting block, every permutation of the factors has the same product,

\[
\prod_{P\in B} e^{-ihc_P P}=e^{-ih\sum_{P\in B}c_P P}.
\]

Thus these internal permutations preserve the ideal formula for every time and step count, as well as the emitted Pauli rotation count. No frozen schedule or simulation was changed after discovering this historical serialization issue.

At the common 100U rotation budget and identical step count, comparisons against the historical saved states give:

- Maximum state-vector norm difference: 1.20e-14.
- Maximum absolute difference in the reported phase-invariant error E: 6.06e-18.
- Maximum absolute difference in particle-number leakage: 2.22e-36.

The remaining 64 baseline sequences match the literal historical hashes. Prospective owned-signed internal word order was separately aligned with its explicitly archived flat sequence before any fixed-campaign error was measured; that pre-execution correction is recorded in `baseline_alignment.json`. The distinction between 96 serialized-design matches, 64 literal historical matches, and 32 certified commuting-equivalent sequences should be retained in any review summary.
