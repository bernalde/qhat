# Corrected low-energy-state control

Eight frozen ten-qubit Hamiltonians; T=1, r=100; three initial states and three fixed chronological Pauli schedules. Ground and next refer to the retained active-space Hamiltonian in the (N_alpha,N_beta)=(2,2) sector, not a specified total-spin sector. All vectors are propagated in the full 1024-dimensional space without postselection.

The Hamiltonian and schedules come directly from controlled-blocks-b53wa_g3 with molecular-unit-weight-v1 coefficients. All initial vectors and designs were frozen before any tested Trotter evolution. Exact references use independent sparse full-space expm_multiply and sector diagonalization. The coherent state uses the predeclared canonical relative phase.

| Initial state | Proposed vs JW signed W/T/L | Proposed vs JW magnitude W/T/L | Best: proposed / JW signed / JW magnitude |
|---|---|---|---|
| sector_ground | 2 / 0 / 6 | 5 / 0 / 3 | 2 / 6 / 0 |
| sector_next | 2 / 0 / 6 | 4 / 0 / 4 | 1 / 6 / 1 |
| equal_coherent | 3 / 0 / 5 | 4 / 0 / 4 | 2 / 5 / 1 |
| all_states | 7 / 0 / 17 | 13 / 0 / 11 | 5 / 17 / 2 |

| Case | State | Proposed E | JW signed E | JW magnitude E | Winner set |
|---|---|---|---|---|---|
| CH2O_s1.00 | sector_ground | 1.87423412e-07 | 1.76828029e-07 | 5.62837044e-07 | jw_signed |
| CH2O_s1.00 | sector_next | 1.74805749e-08 | 4.67297307e-08 | 1.06828089e-07 | fermionic_signed |
| CH2O_s1.00 | equal_coherent | 1.02452027e-07 | 1.11778992e-07 | 3.3483311e-07 | fermionic_signed |
| CH2O_s1.25 | sector_ground | 1.27786851e-07 | 1.98623969e-07 | 1.5557577e-07 | fermionic_signed |
| CH2O_s1.25 | sector_next | 5.05728473e-08 | 7.9745378e-09 | 4.7617526e-08 | jw_signed |
| CH2O_s1.25 | equal_coherent | 9.10502169e-08 | 1.0301606e-07 | 1.01573388e-07 | fermionic_signed |
| H2S_s1.00 | sector_ground | 4.33882548e-08 | 4.41184194e-09 | 4.28415727e-08 | jw_signed |
| H2S_s1.00 | sector_next | 2.3874665e-08 | 2.06471888e-09 | 2.0910672e-08 | jw_signed |
| H2S_s1.00 | equal_coherent | 3.363146e-08 | 3.23828042e-09 | 3.18761225e-08 | jw_signed |
| H2S_s1.25 | sector_ground | 1.23881318e-07 | 1.69947016e-08 | 8.25970273e-08 | jw_signed |
| H2S_s1.25 | sector_next | 9.55007374e-08 | 2.03704234e-08 | 8.17846147e-08 | jw_signed |
| H2S_s1.25 | equal_coherent | 1.09691028e-07 | 1.86825625e-08 | 8.21908211e-08 | jw_signed |
| HCN_s1.00 | sector_ground | 3.37059642e-07 | 7.35335009e-08 | 2.82429462e-07 | jw_signed |
| HCN_s1.00 | sector_next | 1.40386995e-07 | 2.26232239e-08 | 1.83025591e-07 | jw_signed |
| HCN_s1.00 | equal_coherent | 2.46886328e-07 | 4.87637647e-08 | 2.32712429e-07 | jw_signed |
| HCN_s1.25 | sector_ground | 1.06573076e-07 | 3.86982371e-07 | 1.1984846e-07 | fermionic_signed |
| HCN_s1.25 | sector_next | 7.35875808e-08 | 9.58220325e-08 | 3.56817711e-08 | jw_magnitude |
| HCN_s1.25 | equal_coherent | 8.86404252e-08 | 2.41360372e-07 | 7.77312721e-08 | jw_magnitude |
| SiH4_s1.00 | sector_ground | 4.50910171e-08 | 3.29157117e-10 | 5.05931125e-08 | jw_signed |
| SiH4_s1.00 | sector_next | 3.08338882e-08 | 2.86369386e-10 | 3.11018424e-08 | jw_signed |
| SiH4_s1.00 | equal_coherent | 3.54790228e-08 | 3.19848083e-10 | 4.08466409e-08 | jw_signed |
| SiH4_s1.25 | sector_ground | 6.69482804e-08 | 5.22014038e-09 | 6.8387038e-08 | jw_signed |
| SiH4_s1.25 | sector_next | 2.58627451e-08 | 3.0298659e-09 | 2.58730622e-08 | jw_signed |
| SiH4_s1.25 | equal_coherent | 4.63173921e-08 | 4.03804243e-09 | 4.71435445e-08 | jw_signed |

Completed 72/72 Trotter evaluations in 9.87 seconds; tie tolerance and denominator floor 1e-12. All values and final vectors are retained. Eigenvalue gaps, cluster decisions, state hashes, HF overlaps, norms, residuals, and reference checks are saved per case.

Descriptive paired controls on eight frozen Hamiltonians from four molecular families. States and geometries are not independent molecular samples. Sector low-energy does not assert a total-spin quantum number or full-system chemical ground state.
