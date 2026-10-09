# Prospective validation of the frozen diagonal-response ordering

The unchanged phase-2 candidate beat all four comparison schedules in each of six prospective cases. This is evidence of transfer to three additional molecular fixtures; it is not a universal guarantee or a gate-budget claim. The new panel used Li₂, F₂, and NH₃, each at fixed geometry with uniform Hamiltonian coefficient multipliers 1 and 5/4. These are three molecules with paired strength/time controls, not six independent molecules or six molecular geometries.

Candidate source SHA-256: `58263730c8c35b2dedbad845574c3f4f4286af5c0328940d001af85ddd18ecd7`. Its diagonal splitting, two starts, ten best-insertion moves per start, objective, acceptance threshold, and tie rules were unchanged. All six candidates were saved and hashed before any new approximate errors were measured. No exact dynamics enters their construction. Fixtures existed earlier, so this is prospective validation of a frozen algorithm, not blinded data collection.

All methods use the same screened exact dyadic Hermitian-parent target, HF initial state, T = 1, r = 100, scalar convention, and number of Pauli rotations. The coefficient multiplier is applied after the original whole-parent cutoff, preserving support. There is no final individual-Pauli cutoff. F₂ and NH₃ consequently retain 1184 and 5372 words, respectively, instead of the older final-cutoff counts 1176 and 5098. All baselines were rerun against the new common target.

| Case | Candidate E | Owned signed E | JW signed E | JW magnitude E | Complete-class lex E |
|---|---:|---:|---:|---:|---:|
| f2_coeff1p00 | 2.077437e-12 | 1.058153e-07 | 4.345745e-08 | 1.059283e-07 | 5.615985e-11 |
| f2_coeff1p25 | 4.056683e-12 | 2.350472e-07 | 1.552538e-07 | 2.356982e-07 | 7.562039e-11 |
| li2_coeff1p00 | 1.420604e-11 | 2.855298e-10 | 1.142612e-08 | 3.727978e-10 | 2.986668e-10 |
| li2_coeff1p25 | 5.333676e-11 | 7.774991e-10 | 1.239205e-08 | 9.130564e-10 | 8.235108e-10 |
| nh3_coeff1p00 | 2.353729e-10 | 9.631484e-10 | 1.039576e-09 | 9.948193e-10 | 9.168085e-10 |
| nh3_coeff1p25 | 5.169499e-10 | 2.249408e-09 | 2.376190e-09 | 2.345813e-09 | 2.009344e-09 |

E = 1 − |〈normalized exact, normalized approximate〉|, evaluated stably. Absolute tie tolerance and the denominator ratio floor are 1e-12. All comparisons are strict wins under that convention.

| Baseline | Candidate wins/ties/losses | Median candidate / baseline E |
|---|---:|---:|
| owned_signed | 6/0/0 | 0.059176844 |
| jw_signed | 6/0/0 | 0.0027737028 |
| jw_magnitude | 6/0/0 | 0.048261091 |
| class_lex | 6/0/0 | 0.059206444 |

F₂ also improves substantially with plain complete-class ordering: the candidate adds an 18.6–27.0-fold reduction beyond that baseline. Its larger comparison against JW therefore combines a useful partition with a useful ordering. NH₃ gains are more moderate, approximately fourfold. These distinctions matter when attributing the mechanism.

Additional classical ordering time, including pair-amplitude/objective construction and both optimization starts but excluding shared tensor-to-Pauli preprocessing, was about 0.167 s for Li₂, 1.25–1.27 s for F₂, and 20.0–20.2 s for NH₃. The refined partitions contain 121, 348, and 1056 units. These measurements are descriptive CPU timings, not asymptotic complexity bounds or end-to-end quantum resource savings.

Every full Hamiltonian in the retained affine closure was also assembled independently from occupation-basis CAR actions. The maximum Pauli/CAR matrix discrepancy was 1.54e-12 in Frobenius norm; exact-reference state differences were below 1.17e-13. The closure contains all intermediate Pauli flips, including canceled parent words, and retains number/spin leakage. No propagation postselection or sector projection is used.

Independent validation, without QHAT imports, recomputed all 30 saved-state metrics and all six diagonal-reference objectives. The latter were obtained from the second-order local defect of the actual chronological Pauli product, without using the candidate pair-amplitude functions. Maximum metric discrepancy: 1.323e-22. Maximum relative objective discrepancy: 1.031e-13.

An independently implemented in-place Pauli-pair kernel replayed all five Li₂ coefficient-1 schedules and the coefficient-1 F₂/NH₃ candidate schedules, including all 100 steps. Maximum full complex-state distance: 1.481e-15. Maximum overlap-error discrepancy: 2.134e-22. This is independent full propagation for seven selected runs; the other 23 have saved-state/target/reference checks. All checks also run under Python `-O`; a deliberately altered error was rejected.

Campaign runtime was 267.8 s, with peak RSS 1.066 GiB under a 6 GiB cap. The process used at most two threads; the imported frozen candidate retains its one-thread settings.

Limitations: the reported accuracy is HF-state overlap at one physical time and one step count. Small r = 100 errors alone do not establish lower cost at a useful target accuracy. The panel does not test arbitrary initial states, molecular geometries, large systems, gate synthesis, or the novelty of the criterion. A separate state-transfer study was deferred. Do not fold the coefficient-strength controls into the development geometry-scale counts.

The development labels s1.00/s1.25 refer to uniform molecular-coordinate scales: `codex_external_validation/PROTOCOL.md:12` and `external_validation.py:39` in Albert’s historical study explicitly specify scaled coordinates. The held-out labels deliberately use coeff1p00/coeff1p25 to distinguish their different meaning.

Replay: `recheck_validation.py --results results --output <new-json>` requires NumPy only. `independent_circuit_replay.py --results results --output <new-directory>` requires NumPy and Numba. Both import no QHAT. Scientific regeneration with `run_validation.py` currently expects the original pinned QHAT checkout at `/home/bernalde/repos/qhat`; source snapshots, exact target/schedule designs, tensors, and dependency pins are included for relocation. No reviewed repository was modified.
