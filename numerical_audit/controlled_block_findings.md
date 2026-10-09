# What the fixed-block control establishes

The predeclared controlled campaign completed **728 evolutions** in 1,189.49 seconds, using eleven frozen Hamiltonians from seven molecular families. Its complete report and immutable evidence are in [results/controlled-blocks-b53wa_g3/REPORT.md](results/controlled-blocks-b53wa_g3/REPORT.md). The [protocol](controlled_block_protocol.md) was fixed before observing its outcomes.

The evidence supports a **conditional Hartree–Fock advantage over typical random block orders**, not consistent superiority over other deterministic priorities. Ownership was constructed once with signed parent priority and then frozen; none of the alternative schedules was allowed to change block membership or its internal Pauli sequence.

| Comparator, with fixed blocks | Signed wins / ties / losses at T=0.25 | At T=1 |
|---|---|---|
| Random-order median, twenty predeclared seeds | 9 / 0 / 2 | 9 / 0 / 2 |
| Descending absolute owner-parent coefficient | 5 / 2 / 4 | 6 / 0 / 5 |
| Descending retained block coefficient l1 norm | 6 / 2 / 3 | 7 / 0 / 4 |

Both primary times use r=100. The two cases where signed ordering loses to the random median are H2S_s1.00 and SiH4_s1.25. It beats all twenty sampled random orders for HCN_s1.25, Li2 and NH3 at both times; this is a comparison with the fixed sampled distribution, not a claim of global optimality. At T=1 the random median error is 222.62 times the signed error for Li2 and 160.88 times for NH3.

Reversing signed whole-block chronology changes E by at most 0.025% at T=1 and 0.002% at T=0.25. The signed-versus-reversed win/loss count should therefore not be interpreted as a substantial directional advantage. These nearly equal errors are consistent with the leading commutator changing sign under reversal while its leading variance is unchanged.

Across all 26 schedules in the three predeclared convergence anchors, the 156 ratios E(50)/E(100) and E(100)/E(200) lie between 3.997968 and 4.015979. This supports the expected r^-2 overlap-error scaling in the tested regime. All completed-block treatments retain particle and spin sector probabilities to a maximum observed leakage of 8.1424e-27. Maximum output/reference norm defect is 1.0045e-12. Finite coefficient-cutoff symmetry residuals remain recorded explicitly; no physical-sector projection is used.

The JW comparisons remain distinct external controls because they also change block contiguity. At T=1 signed fermionic ordering wins 8/11 against JW magnitude and 4/11 against JW signed. This controlled experiment therefore strengthens the case for a useful structured heuristic in selected HF settings; it does not establish that signed parent priority is uniformly best, that ownership itself is optimal, or that these benefits extend to arbitrary states. Geometries, times and random seeds are dependent repeats, not 728 independent molecules.

Independent artifact checks validate all 728 unique rows, thirteen pinned source/snapshot hashes, all eleven complete block designs, common Hamiltonian and budget identities, the E-to-infidelity conversion, and all 22 primary signed-rank calculations. Twenty targeted tests pass, including literal dense block-exponential oracles and synthetic log-plot edge cases. The scientific campaign used the frozen sources in `source_snapshot`. A separately documented, tested presentation-only patch preserves legends when an initial case is omitted, marks zero/below-floor errors, and makes the figures readable at manuscript width; see `presentation_provenance.json` and `post_presentation_artifact_validation.json`.

All code, inputs and outcomes remain local on the audit branch. No result or branch has been published.
