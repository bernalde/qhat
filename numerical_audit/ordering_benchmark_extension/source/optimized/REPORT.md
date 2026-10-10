# Optimized and Strang controls: recovered campaign

Completed 728 propagations from 14 frozen molecular targets. All comparisons below use a requested budget of 100U Pauli rotations; Strang uses the largest feasible integer step count. The eight development geometries and six prospective coefficient-scaled targets are separate panels, not fourteen independent samples.

The candidate is the conserving refined JD first-order ordering. W/T/L compares its error against each comparator with an absolute tolerance of 1e-12. Ratios whose comparator error is at or below that floor are omitted. This is a comparison within this campaign, not a conclusion against all published/fixed controls.

| Comparator | Development W/T/L | Prospective W/T/L | All-case median candidate/comparator | Ratios |
|---|---:|---:|---:|---:|
| class_refined_J0 | 2/6/0 | 2/4/0 | 1 | 14 |
| class_pauli_JD | 2/0/6 | 2/0/4 | 1.16995 | 14 |
| class_pauli_JD_budget | 8/0/0 | 6/0/0 | 0.18284 | 14 |
| strang_class_refined_response | 0/0/8 | 0/0/6 | 9427.69 | 2 |
| strang_jw_magnitude | 0/0/8 | 0/2/4 | 1258.58 | 7 |
| strang_jw_signed | 0/0/8 | 4/0/2 | 184.783 | 9 |
| strang_class_lex | 0/0/8 | 0/0/6 | 4118.13 | 5 |

The matched Pauli arm controls insertion-destination evaluations, not wall-clock time or setup. Complete setup/search/evaluation counts, leakage, and least tested actual rotations at errors1e-4,1e-6,1e-8 are in `summary.json`. Strang changes formula order; its comparison is not evidence for or against ordering alone. The objective, pulse-fusion, checkpoint-provenance, and literal archived-schedule checks are separate artifacts. Independent replay/objective checks and complete cross-method rankings belong to the enclosing campaign.
