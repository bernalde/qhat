# H2O reconstruction sweep and small state control

These are regenerated inputs, not a bitwise replay of unavailable historical tensors. Each within-case comparison fixes H, state, time, and Pauli exponential count. Active-space changes alter H, electron count, and frozen-core contributions.

126 HF evolutions completed: 15 cases × two tensor policies × three methods at T=1,r=100, plus three predeclared anchors at r=50/200. A separate prescribed two-case low-energy panel adds 36 evolutions. These repeated settings are not independent molecular samples.

## Primary outcomes

Wins/ties/losses for the proposed owned signed schedule; absolute tie tolerance 1e-12. Ratios suppressed when the comparator error is at or below 1e-12.

| Tensor policy | Comparator | Wins / ties / losses | Median owned/comparator ratio |
|---|---|---|---|
| historical_policy | jw_signed | 5 / 2 / 8 | 2.36224 |
| historical_policy | jw_magnitude | 5 / 4 / 6 | 1 |
| unpruned | jw_signed | 5 / 2 / 8 | 2.36224 |
| unpruned | jw_magnitude | 5 / 4 / 6 | 1 |

Both tensor policies use corrected downstream coefficient accumulation and the same final 1e-12 Pauli cutoff. Historical tensor policy preserves the generator’s old truncation, not its downstream metric/accumulation defects. The two policies are sensitivity arms and must not be pooled as independent cases.

## Absolute primary errors (historical tensor policy)

| Basis | Qubits / electrons | Owned signed E | JW signed E | JW magnitude E |
|---|---|---|---|---|
| hgbs-5 | 4 / 2 | 2.31550634e-12 | 2.31550634e-12 | 2.31550634e-12 |
| hgbs-5 | 6 / 2 | 2.37230008e-12 | 1.93378090e-12 | 2.37230008e-12 |
| hgbs-5 | 8 / 4 | 3.86368349e-11 | 2.61207593e-12 | 3.16253386e-12 |
| hgbs-5 | 10 / 4 | 5.73696558e-11 | 3.72932123e-12 | 4.13267137e-12 |
| hgbs-5 | 12 / 4 | 5.79638657e-10 | 2.00627024e-10 | 4.44086496e-10 |
| hgbs-5 | 14 / 6 | 9.10970735e-10 | 6.70816683e-09 | 1.78798396e-09 |
| hgbs-5 | 16 / 6 | 2.10498060e-09 | 2.28554930e-07 | 3.00744497e-07 |
| hgbs-5 | 18 / 8 | 1.28941295e-08 | 1.43104931e-06 | 1.53841968e-06 |
| hgbs-5 | 20 / 8 | 1.81274244e-08 | 2.55504209e-06 | 2.80266573e-06 |
| sto-6g | 4 / 2 | 2.30350684e-08 | 2.29724385e-10 | 2.30350684e-08 |
| sto-6g | 6 / 2 | 3.91706851e-08 | 1.65820326e-08 | 3.91706851e-08 |
| sto-6g | 8 / 4 | 2.31930374e-07 | 2.12462257e-08 | 2.36789270e-07 |
| sto-6g | 10 / 6 | 1.45060247e-06 | 1.84541704e-06 | 1.36636906e-06 |
| sto-6g | 12 / 8 | 1.07109274e-05 | 3.87356335e-06 | 3.71118548e-06 |
| sto-6g | 14 / 10 | 1.42830138e-05 | 3.91897999e-06 | 1.37815867e-05 |

## Low-energy panel

For both eight-qubit active models, the lowest fixed-(2,2) state is a singlet and the next is an Ms=0 triplet (S²=2, negligible variance). Their equal coherent superposition has <S²>=1 and variance1. All evolutions retain the full256-dimensional space.

| Tensor policy | Comparator | Owned wins / ties / losses across six case/state pairs |
|---|---|---|
| historical_policy | jw_signed | 2 / 0 / 4 |
| historical_policy | jw_magnitude | 5 / 0 / 1 |
| unpruned | jw_signed | 2 / 0 / 4 |
| unpruned | jw_magnitude | 5 / 0 / 1 |

## Checks and interpretation

The proposed schedule is uniquely best among the three in only four of fifteen HF cases. Against JW signed sorting its median error ratio is 2.36224. The selective benefits in larger HGBS-5 active spaces coexist with losses elsewhere; they do not establish a general accuracy advantage.

All 45 historical-policy primary errors agree with the archived table within 3.91993e-13 absolute. The largest relative discrepancy is 3.76% for a floor-adjacent error. This substantially supports the specific H2O table while not certifying identical missing historical tensor files.

Policy sensitivity: retained H and schedules are identical in 14 cases. Only the 20-qubit unpruned case adds four Pauli coefficients of about 1.15746e-12; its retained-H coefficient-l1 difference is 4.62985e-12. The maximum primary E change between policies is 9.74088e-21; outcomes are unchanged.

Validation maxima: {"max_affine_full_sample_distance": 0.0, "max_dense_pilot_distance": 3.7170592637564554e-14, "max_low_energy_dense_product_distance": 1.2631967801141352e-13, "max_low_energy_reference_dense_eigh_distance": 5.4508808602176836e-15, "max_norm_defect": 7.639444632445702e-13, "max_reference_off_spin_sector_frobenius": 7.792442607243671e-17, "max_reference_split_distance": 4.8254127770527304e-15, "peak_RSS_GiB": 4.0284271240234375}.

Small HF states reproduce independently dense-verified pilot states. All36 low-energy evolutions additionally agree with independently formed dense product matrices. The full sweep exact reference uses the retained Pauli sum on the complete Pauli-invariant closure, including round-off off-sector terms. A separate checker recomputes all162 final-state metrics through the projector-residual identity and verifies input/design/state hashes and Pauli closure. This checks the large saved results but is not a second independent large-state evolution.

Archived tally: {"tie": 4, "jw_magnitude_win": 6, "fermionic_win": 5}. Archived values are retained separately in archived_h2o_rows.csv and h2o_archived_panel.pdf. The new historical-convention classifications are {"historical_policy": {"jw_magnitude": {"loss": 6, "tie": 4, "win": 5}, "jw_signed": {"loss": 9, "tie": 1, "win": 5}}, "unpruned": {"jw_magnitude": {"loss": 6, "tie": 4, "win": 5}, "jw_signed": {"loss": 9, "tie": 1, "win": 5}}}.

See summary.json for all18 anchor/method entries (two step-convergence ratios each), per-case policy equality, archived-versus-reconstructed differences, raw leakage and norms. h2o_reconstruction.pdf displays paired absolute errors and ratios; it is not a qubit-only convergence law.
