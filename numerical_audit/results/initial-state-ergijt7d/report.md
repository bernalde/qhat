# Initial-state sensitivity control

At T=1 and r=100, each ordering evolves the same initial state under the same retained Hamiltonian with the same nominal Pauli-exponential count. Exact reference evolution and Trotter propagation both use the full 1024-dimensional space; there is no sector postselection.

Five seeds were declared before evolution for each case. Complex Haar states cover the full fixed-(N_alpha,N_beta) sector, including components outside the HF-connected affine space. These are exploratory controls, not chemically representative low-energy states.

| Case | HF winner(s) | Haar wins: fermionic / JW signed / JW magnitude | Haar states with no HF winner best |
|---|---|---|---|
| CH2O_s1.00 | jw_signed | 0 / 5 / 0 | 0 / 5 |
| CH2O_s1.25 | fermionic_signed | 0 / 5 / 0 | 5 / 5 |
| H2S_s1.00 | jw_signed | 0 / 5 / 0 | 0 / 5 |
| H2S_s1.25 | jw_signed | 0 / 5 / 0 | 0 / 5 |
| HCN_s1.00 | jw_signed | 0 / 5 / 0 | 0 / 5 |
| HCN_s1.25 | fermionic_signed | 0 / 3 / 2 | 5 / 5 |
| SiH4_s1.00 | jw_signed | 0 / 5 / 0 | 0 / 5 |
| SiH4_s1.25 | jw_signed | 0 / 5 / 0 | 0 / 5 |

Completed 144 / 144 ordering-state comparisons in 19.7 seconds. Winners use a predeclared absolute overlap-error tie tolerance of 1e-12; ties count for every tied ordering.

Conservation and accuracy answer different questions: retaining a physical symmetry does not establish a uniformly better approximation. These data only test the specified states, times, budgets, and historical retained Hamiltonians. Input hashes, source hashes, state hashes, exact-reference crosschecks, norm defects, and leakage are stored in the adjacent JSON files. No ordering is selected separately for each metric after inspecting the outcome.
