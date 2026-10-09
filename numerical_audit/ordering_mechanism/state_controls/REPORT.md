# H₂O larger-state control

This exploratory control addresses state dependence of two previously successful HF cases. It introduces no new ordering rule. Both are reconstructed HGBS-5 active spaces of one molecule, in the original frozen orbital gauge. They are not independent molecular samples.

The target coefficients and three chronological schedules are exactly those of QHAT revision `0be7d4857537668eec957dabaa8438ea5c0d33f2`. Every new evolution uses T = 1 and r = 100 in the full Hilbert space, with the same omitted scalar. No particle- or spin-sector projection is applied during propagation.

The state-preparation protocol was frozen before approximate errors: the first canonical vector of the complete ground-energy cluster, the first canonical vector of the next distinct-energy cluster in the full Nα = Nβ = 3 sector, and their equal coherent superposition. Initial arrays and hashes were written before approximate propagation. This differs deliberately from taking the second listed eigenvector when a ground cluster is degenerate.

| Qubits | State | Owned E | JW signed E | JW magnitude E | JW signed / owned | JW magnitude / owned |
|---:|---|---:|---:|---:|---:|---:|
| 14 | HF_archived | 9.109707e-10 | 6.708167e-09 | 1.787984e-09 | 7.364 | 1.963 |
| 14 | sector_ground | 9.163230e-10 | 6.780465e-09 | 1.811485e-09 | 7.4 | 1.977 |
| 14 | sector_next_distinct_energy | 1.686326e-08 | 3.263397e-08 | 4.554524e-08 | 1.935 | 2.701 |
| 14 | equal_coherent | 8.889792e-09 | 1.970722e-08 | 2.367837e-08 | 2.217 | 2.664 |
| 16 | HF_archived | 2.104981e-09 | 2.285549e-07 | 3.007445e-07 | 108.6 | 142.9 |
| 16 | sector_ground | 2.132458e-09 | 2.281335e-07 | 3.003417e-07 | 107 | 140.8 |
| 16 | sector_next_distinct_energy | 5.124720e-08 | 3.841909e-08 | 6.322615e-08 | 0.7497 | 1.234 |
| 16 | equal_coherent | 2.668983e-08 | 1.332763e-07 | 1.817839e-07 | 4.994 | 6.811 |

E = 1 − |〈normalized exact, normalized approximate〉|, evaluated through a stable phase-aligned distance. A ratio above one favors owned ordering. HF rows are rechecked archived states, not new propagations.

Both selected ground states are singlets and almost HF determinants; the next distinct-energy states are Ms = 0 triplets with negligible HF overlap. The coherent states have 〈S²〉 = 1 and Var(S²) = 1. They are spin mixtures, not excited singlets. Full-H eigenvector residuals and small off-sector target entries are measured and retained.

The independent checks include written-letter sparse Hamiltonian assembly; tensor-axis Hamiltonian action and a sampled twelve-rotation product checked against QHAT; dense/ARPACK eigenvalue agreement; full sparse exponential references compared with split-time action and eigenphase predictions. The standalone saved-state checker imports no QHAT and recomputes hashes, stable errors, and number/spin leakage. It is not a second independent full product-formula propagation.

Energy expectations are normalized and reported as drift from the exact evolved state, in Hartree. These are an additional observable, not interchangeable with overlap error. `coherent_diagnostics.json` reports relative phase and population outside the initial eigenpair; this is descriptive post-run analysis.

The 14-qubit case was run first. Its 33-second runtime justified extending the already frozen protocol to 16 qubits. That extension was authorized independently of observed errors. Each process had a 6 GiB resident-memory cap and two threads.

Run `recheck_saved.py --root <this directory> --output <new JSON path>` with the pinned environment; it supports Python `-O`. Scientific replay of `run_state_control.py` currently expects the original pinned QHAT checkout at `/home/bernalde/repos/qhat`; changing that checkout requires an isolated copy at the pinned revision or a separately documented path adaptation. Source snapshots and exact coefficient/schedule designs are included.

These controls do not resolve orbital-sign sensitivity, geometry transfer, comparative gate cost, or performance on unseen molecules. In particular, success in the original orbital gauge is not evidence of invariance of a signed coefficient ordering under physically equivalent orbital phases.
