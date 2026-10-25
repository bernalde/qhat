# Derived ordering on correlated/coherent state and time holdouts

At the primary 10U rotation cap, newly state/time-aware placement has 24 wins, 0 ties and 0 losses against transfer of the frozen HF/T=1 placement order, with median error ratio 0.325126 and 0 baseline-floor exclusions. This tests a new state/time-aware search, not transfer alone. The control comparisons below retain formula changes, unequal search work and all losses.

Against state_symmetric at the same 10U cap, the first-order placement candidate has 0 wins, 4 ties and 20 losses. This comparison changes the product formula; it does not isolate ordering alone.

Against hf_symmetric at the same 10U cap, the first-order placement candidate has 0 wins, 4 ties and 20 losses. This comparison changes the product formula; it does not isolate ordering alone.

This prospective panel compares transfer of frozen HF-trained schedules with newly state-aware diagonal-response placement and insertion. The orbital basis and retained target stay fixed within each family. All 432 circuits passed independent full-space replay. These are 24 related state/time inputs from three molecular families, not independent statistical replicates.

## State and objective scope

The inputs are sector ground and first-excited states and (g+e)/sqrt(2), (g+i e)/sqrt(2), at T=0.37 and 1.73 inverse Hartree. One-particle density spectra classify 12/12 selected states as correlated beyond a single Slater determinant under the declared 1e-8 idempotency-defect threshold. Participation ratios separately describe their expansion in the fixed orbital basis. These sector labels do not claim unrestricted global ground/excited states.

The one-particle density spectrum diagnoses departure from a single Slater determinant for pure fixed-number states; see [Gigena, Di Tullio and Rossignoli, Sec. II.A](https://arxiv.org/pdf/2001.03570). The declared numerical idempotency threshold and every occupation are retained in the state table.

State-aware J_D sums complex transported commutator amplitudes before the parallel-state projection and squared norm. It retains coherent interference. This is the generalized derived diagonal-response rule; earlier historical signed-parent low-energy tests do not validate it. The failed two-level reference correction is unchanged and excluded from new tuning.

## Coverage and all controls

| State | T | Informative families at 10U | Coverage |
| --- | --- | --- | --- |
| ground | 0.37 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| ground | 1.73 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| excited | 0.37 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| excited | 1.73 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| coherent_real | 0.37 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| coherent_real | 1.73 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| coherent_imag | 0.37 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |
| coherent_imag | 1.73 | Li2_r2p95, F2_r1p55, H2O_scale1p10 | covered |

Coverage requires at least two of three families with frozen-HF current-placement E>1e-12 in each state/time stratum. Failed strata remain inconclusive and are not replaced. All absolute errors and losses remain in the tables. The denominator floor is 1e-12; ratios exclude and count floor-limited baselines.

| Method / frozen-HF current | Budget | Wins / ties / losses (24 inputs) | Median error ratio | Floor exclusions |
| --- | --- | --- | --- | --- |
| state_current | 10U | 24/0/0 | 0.325126 | 0 |
| state_JD | 10U | 22/0/2 | 0.262639 | 0 |
| state_pauli_JD | 10U | 22/0/2 | 0.165579 | 0 |
| hf_pauli_large | 10U | 18/0/6 | 0.822882 | 0 |
| hf_symmetric | 10U | 24/0/0 | 0.000402254 | 0 |
| state_symmetric | 10U | 24/0/0 | 0.000257653 | 0 |
| jw_signed | 10U | 14/0/10 | 0.963666 | 0 |
| jw_magnitude | 10U | 8/0/16 | 1.43674 | 0 |
| state_current | 30U | 24/0/0 | 0.325089 | 0 |
| state_JD | 30U | 22/0/2 | 0.262248 | 0 |
| state_pauli_JD | 30U | 22/0/2 | 0.165186 | 0 |
| hf_pauli_large | 30U | 18/0/6 | 0.773636 | 0 |
| hf_symmetric | 30U | 24/0/0 | 4.46845e-05 | 0 |
| state_symmetric | 30U | 24/0/0 | 2.61038e-05 | 0 |
| jw_signed | 30U | 14/0/10 | 0.96359 | 0 |
| jw_magnitude | 30U | 8/0/16 | 1.43611 | 0 |

The state/time-aware placement candidate against every strong control at 10U:

| Denominator control | Candidate wins / ties / losses | Median candidate/control error ratio | Floor exclusions |
| --- | --- | --- | --- |
| hf_current | 24/0/0 | 0.325126 | 0 |
| state_JD | 8/4/12 | 1.20648 | 4 |
| state_pauli_JD | 8/4/12 | 1.81496 | 4 |
| hf_pauli_large | 15/0/9 | 0.475626 | 0 |
| hf_symmetric | 0/4/20 | 142.68 | 7 |
| state_symmetric | 0/4/20 | 305.379 | 7 |
| jw_signed | 20/4/0 | 0.306328 | 4 |
| jw_magnitude | 20/4/0 | 0.0610649 | 4 |

[Candidate/control contrasts](candidate_contrasts.csv) retain both budgets and every family separately. Losses remain in the win/tie/loss counts even when a denominator is below the descriptive ratio floor. State and time dependence are changed together in the new search; this panel does not isolate those two effects.

The two coherent phases tie in 108/108 paired physical-error comparisons at the declared 1e-12 threshold. [All phase pairs](coherent_phase_pairs.csv) retain errors, objective scores and schedule identity. These related phase inputs do not demonstrate general empirical phase sensitivity. Separate synthetic fixed-number/conserving-unit quadrature tests exercise nontrivial coherent interference; no panel case is replaced after this observation.

hf_current, hf_pauli_large, hf_symmetric and fixed JW controls transfer the previous HF/T=1 schedules unchanged. state_current runs generalized placement. state_JD and state_pauli_JD newly optimize the current state/time objective. Both symmetric arms reuse first-order schedules; neither is a second-order-aware optimization. State Pauli search has an explicitly larger position allowance, while the frozen HF Pauli control retains the earlier still-larger search; this is not a comparison of matched CPU budgets.

## Absolute primary errors, by family

| Family | State | T | hf_current | state_current | state_JD | state_pauli_JD | hf_pauli_large | hf_symmetric | state_symmetric | jw_signed | jw_magnitude |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| F2_r1p55 | coherent_imag | 0.37 | 4.405096e-07 | 3.307709e-07 | 2.470937e-07 | 1.738848e-09 | 5.963783e-07 | 1.090305e-10 | 1.035072e-10 | 5.543344e-07 | 5.374206e-07 |
| F2_r1p55 | coherent_imag | 1.73 | 1.925166e-04 | 1.023207e-04 | 9.524367e-05 | 8.773970e-07 | 2.579823e-04 | 9.245240e-07 | 8.774008e-07 | 2.125245e-04 | 2.063409e-04 |
| F2_r1p55 | coherent_real | 0.37 | 4.405096e-07 | 3.307709e-07 | 2.470937e-07 | 1.738848e-09 | 5.963783e-07 | 1.090305e-10 | 1.035072e-10 | 5.543344e-07 | 5.374206e-07 |
| F2_r1p55 | coherent_real | 1.73 | 1.925166e-04 | 1.023207e-04 | 9.524367e-05 | 8.773970e-07 | 2.579823e-04 | 9.245240e-07 | 8.774008e-07 | 2.125245e-04 | 2.063409e-04 |
| F2_r1p55 | excited | 0.37 | 8.658846e-07 | 5.737392e-31 | 5.475547e-31 | 8.079395e-31 | 1.184057e-06 | 2.615186e-30 | 2.615186e-30 | 6.831092e-31 | 8.079395e-31 |
| F2_r1p55 | excited | 1.73 | 3.788979e-04 | 3.188913e-30 | 2.144444e-30 | 2.119764e-30 | 5.086793e-04 | 2.077786e-30 | 2.036470e-30 | 3.235348e-30 | 2.119764e-30 |
| F2_r1p55 | ground | 0.37 | 1.512724e-08 | 8.541359e-09 | 5.131355e-09 | 4.933731e-09 | 8.474459e-09 | 2.113433e-10 | 2.796970e-11 | 1.108651e-06 | 1.074824e-06 |
| F2_r1p55 | ground | 1.73 | 6.089658e-06 | 4.244158e-06 | 2.302394e-06 | 4.164589e-06 | 5.224053e-06 | 1.778155e-06 | 2.348280e-07 | 4.249554e-04 | 4.125881e-04 |
| H2O_scale1p10 | coherent_imag | 0.37 | 9.348690e-08 | 1.171181e-08 | 4.553411e-08 | 3.755574e-08 | 7.692866e-08 | 5.207200e-11 | 3.274364e-11 | 4.562403e-08 | 5.715222e-07 |
| H2O_scale1p10 | coherent_imag | 1.73 | 3.344127e-05 | 3.982993e-06 | 1.747857e-05 | 1.299998e-05 | 3.049882e-05 | 3.467690e-07 | 1.828297e-07 | 1.436721e-05 | 1.734663e-04 |
| H2O_scale1p10 | coherent_real | 0.37 | 9.348690e-08 | 1.171181e-08 | 4.553411e-08 | 3.755574e-08 | 7.692866e-08 | 5.207200e-11 | 3.274364e-11 | 4.562403e-08 | 5.715222e-07 |
| H2O_scale1p10 | coherent_real | 1.73 | 3.344127e-05 | 3.982993e-06 | 1.747857e-05 | 1.299998e-05 | 3.049882e-05 | 3.467690e-07 | 1.828297e-07 | 1.436721e-05 | 1.734663e-04 |
| H2O_scale1p10 | excited | 0.37 | 1.693735e-07 | 6.862621e-09 | 2.717346e-08 | 2.653344e-08 | 1.526710e-07 | 3.668616e-11 | 3.312559e-12 | 1.099849e-08 | 3.032956e-07 |
| H2O_scale1p10 | excited | 1.73 | 6.113151e-05 | 3.104703e-06 | 1.137450e-05 | 9.833908e-06 | 6.052298e-05 | 2.539976e-07 | 2.847724e-08 | 3.909771e-06 | 1.010795e-04 |
| H2O_scale1p10 | ground | 0.37 | 1.759712e-08 | 9.624861e-09 | 2.983601e-08 | 4.221456e-08 | 1.186167e-09 | 6.745755e-11 | 5.685037e-11 | 8.024929e-08 | 8.397474e-07 |
| H2O_scale1p10 | ground | 1.73 | 5.722780e-06 | 4.043599e-06 | 8.933141e-06 | 1.536433e-05 | 4.739318e-07 | 4.395370e-07 | 3.682103e-07 | 2.482258e-05 | 2.458511e-04 |
| Li2_r2p95 | coherent_imag | 0.37 | 1.326737e-09 | 4.330251e-10 | 1.382693e-10 | 2.280635e-10 | 3.243244e-10 | 6.357701e-15 | 6.390331e-15 | 1.290962e-09 | 1.191974e-09 |
| Li2_r2p95 | coherent_imag | 1.73 | 6.322772e-07 | 2.047752e-07 | 6.530536e-08 | 1.076723e-07 | 1.341107e-07 | 6.595325e-11 | 6.588299e-11 | 6.033794e-07 | 5.617499e-07 |
| Li2_r2p95 | coherent_real | 0.37 | 1.326737e-09 | 4.330251e-10 | 1.382693e-10 | 2.280635e-10 | 3.243244e-10 | 6.357701e-15 | 6.390331e-15 | 1.290962e-09 | 1.191974e-09 |
| Li2_r2p95 | coherent_real | 1.73 | 6.322772e-07 | 2.047752e-07 | 6.530536e-08 | 1.076723e-07 | 1.341107e-07 | 6.595325e-11 | 6.588299e-11 | 6.033794e-07 | 5.617499e-07 |
| Li2_r2p95 | excited | 0.37 | 2.394399e-09 | 6.518992e-31 | 6.764548e-31 | 9.320257e-31 | 4.609227e-10 | 3.134169e-31 | 3.134169e-31 | 1.524092e-31 | 9.320257e-31 |
| Li2_r2p95 | excited | 1.73 | 1.142236e-06 | 1.857121e-31 | 1.471933e-31 | 2.535045e-31 | 2.185994e-07 | 1.857115e-31 | 1.857115e-31 | 1.903338e-31 | 2.535045e-31 |
| Li2_r2p95 | ground | 0.37 | 2.590721e-10 | 2.179269e-10 | 3.339169e-11 | 4.065181e-12 | 1.877258e-10 | 1.138273e-14 | 3.326208e-15 | 2.581923e-09 | 2.383945e-09 |
| Li2_r2p95 | ground | 1.73 | 1.222835e-07 | 1.027432e-07 | 1.586010e-08 | 1.254009e-09 | 4.962026e-08 | 1.179765e-10 | 3.428026e-11 | 1.206743e-06 | 1.123483e-06 |

[All rows](all_rows.csv) retain both ideal-rotation budgets, formula, affordable steps, actual rotations, absolute state error, infidelity, norm and number/joint spin-resolved-sector leakage. The historical spin_projection_leakage column denotes the joint (N_alpha,N_beta) sector, not pure S_z leakage. [All contrasts](contrasts.csv) retain each family and state/time stratum. Symmetric circuits use exact central/boundary fusion and their actual affordable integer steps. Individual-Pauli controls may leak, and propagation always retains the entire Hilbert space.

## Preparation and classical work

| Family | Preparation + correlation diagnostics CPU s | Wall s |
| --- | --- | --- |
| F2_r1p55 | 0.0392265 | 0.118325 |
| H2O_scale1p10 | 0.0555815 | 0.200162 |
| Li2_r2p95 | 0.0131488 | 0.0170055 |

Preparation includes retained-H construction, sector eigensystem, explicit state vectors and density-matrix diagnostics. It is charged once per family; these exact-diagonalization states do not supply an efficient or compiled hardware preparation routine. Original molecular integral generation is reused from authenticated prior evidence.

[State diagnostics](states.csv), [stand-alone arm construction costs](construction_costs.csv) and [actual campaign construction costs](campaign_construction_costs.csv) retain the detailed costs. Conserving feature setup is shared by placement/insertion; the stand-alone table assigns the full feature cost to each method, while the campaign table counts actual setup/search exactly once. state_symmetric inherits state_JD construction, not another search. Transferred HF schedules have zero new search, with historical construction costs listed separately. These single timing observations under a shared-machine quota do not establish controlled speedups or scaling.

## Verification and limitations

Independent literal-Pauli/BCH scoring and unfused full-space replay checked every row: maximum state-vector difference 2.49492e-15, physical metric difference 6.31006e-17, response-score difference 1.49186e-16. Synthetic tests include complex quadrature, coherent-versus-mixture separation, determinant reduction, updates, degenerate state selection, Slater correlation, pulse fusion and resource guards.

The frozen protocol and all state/schedule hashes precede physical errors. The earlier pre-error correlation-diagnostic amendment left state-file bytes and all 216 schedules unchanged. The numerical source and active inputs are authenticated in the result freeze. No schedule is chosen from physical errors, and no trial is replaced after viewing errors.

This bounded panel establishes only the reported state/time and family behavior. It does not isolate conservation as an accuracy mechanism, prove J_D ranks physical errors, demonstrate matched classical efficiency, provide Lean proofs, optimize a second-order objective, or measure compiled/hardware resources. Those remain separate issues #5–#9. The completed factorial, Tranter checks and initial 8! pilot remain prior completed evidence.
