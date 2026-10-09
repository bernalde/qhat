# Explicit coefficient-policy comparison

Compared the completed historical-retained and explicit-cutoff runs at T=1,
r=100 for all eleven fixtures and all three orderings (33 matched rows).
Within each run, all orderings use the same Hamiltonian, state and work count.
Across these two runs, the construction policy changes; this is a sensitivity
analysis, not an equal-Hamiltonian/cost comparison between methods.

No winning ordering changes. The largest relative change in the stable state
error E is **1.3753e-6** (about 1.38 parts per million), in NH3 fermionic signed.
This does not explain the much larger unresolved historical F2 discrepancy.
It also does not excuse the implicit-cutoff defect: a requested tolerance must
be respected, and small terms can matter for other inputs or accuracy targets.

NH3 has 5124 retained nonidentity Pauli terms instead of 5098 (28 added, two
removed). Its coefficient-l1 change is 3.62e-7. The other ten cases retain their
Pauli support with only roundoff-scale coefficient changes. Exact signed and
magnitude sorting can respond to tiny tie changes, so this comparison is not
an isolated perturbation of coefficients at a frozen schedule.

The corrected owned blocks have no unowned Pauli terms and no noncommuting
internal pairs on any fixture. The maximum reported owned-block number-
commutator coefficient-l1 residual is 2.56e-13 in NH3, consistent with a small
floating arithmetic/truncation residual rather than an exact symbolic proof.
The explicit runner rejects any imaginary retained Pauli coefficient instead
of silently projecting it to a real Hamiltonian.

Detailed paired rows, exponential counts and source hashes are in
[precision_comparison.json](precision_comparison.json).

| Case | Historical winner | Explicit-cutoff winner |
|---|---|---|
| CH2O_s1.00 | jw_signed | jw_signed |
| CH2O_s1.25 | fermionic_signed | fermionic_signed |
| H2S_s1.00 | jw_signed | jw_signed |
| H2S_s1.25 | jw_signed | jw_signed |
| HCN_s1.00 | jw_signed | jw_signed |
| HCN_s1.25 | fermionic_signed | fermionic_signed |
| SiH4_s1.00 | jw_signed | jw_signed |
| SiH4_s1.25 | jw_signed | jw_signed |
| f2 | jw_signed | jw_signed |
| li2 | fermionic_signed | fermionic_signed |
| nh3 | fermionic_signed | fermionic_signed |
