# Strang transfer and orientation controls

All 14 previously frozen targets and JD orders are retained. Five new arms give 420 evolutions; 84 forward conserving states are reused with exact provenance. No search or data-dependent selection is performed. The first chronological unit is outermost. Pauli arms are singleton products; conserving reverse reverses units while retaining each commuting unit's word order. All pairings use the same H, HF state, T=1 and nominal Pauli-rotation budget.

## Accuracy at 10U

Absolute ties use 1e-12; median ratios include comparator errors above 1e-12. A smaller ratio favors conserving Strang.

| Conserving orientation | Comparator | W/T/L | Median E ratio | Ratio cases |
|---|---|---:|---:|---:|
| class_refined_response__forward | class_pauli_JD__forward | 9/1/4 | 0.83079557 | 14 |
| class_refined_response__forward | class_pauli_JD__reverse | 12/0/2 | 0.63481937 | 14 |
| class_refined_response__forward | class_pauli_JD_budget__forward | 13/0/1 | 0.18525091 | 14 |
| class_refined_response__forward | class_pauli_JD_budget__reverse | 14/0/0 | 0.1188814 | 14 |
| class_refined_response__forward | class_refined_response__reverse | 11/0/3 | 0.56913442 | 14 |
| class_refined_response__reverse | class_pauli_JD__forward | 4/0/10 | 1.3410888 | 14 |
| class_refined_response__reverse | class_pauli_JD__reverse | 8/0/6 | 0.8059798 | 14 |
| class_refined_response__reverse | class_pauli_JD_budget__forward | 12/0/2 | 0.62920777 | 14 |
| class_refined_response__reverse | class_pauli_JD_budget__reverse | 12/0/2 | 0.34564411 | 14 |
| class_refined_response__reverse | class_refined_response__forward | 3/0/11 | 1.7649076 | 14 |

## Orientation robustness at 10U

| Pauli search | Conserving better for all four orientation pairs | Pauli better for all four | Orientation-dependent or tied |
|---|---:|---:|---:|
| class_pauli_JD | 3 | 2 | 9 |
| class_pauli_JD_budget | 11 | 0 | 3 |

## First successful tested budget level

Each entry is conserving lower / same / higher tested level. Same-level exact rotation differences are not evidence of a resolved minimum-work ranking. Success is checked separately; counts below use joint successes. No interpolation is used.

| Conserving orientation | E target | Comparator | Lower/same/higher | Candidate-only / comparator-only / neither success |
|---|---|---|---:|---:|
| class_refined_response__forward | 0.0001 | class_pauli_JD__forward | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 0.0001 | class_pauli_JD__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 0.0001 | class_pauli_JD_budget__forward | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 0.0001 | class_pauli_JD_budget__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 0.0001 | class_refined_response__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 1e-06 | class_pauli_JD__forward | 0/13/1 | 0/0/0 |
| class_refined_response__forward | 1e-06 | class_pauli_JD__reverse | 0/13/1 | 0/0/0 |
| class_refined_response__forward | 1e-06 | class_pauli_JD_budget__forward | 5/8/1 | 0/0/0 |
| class_refined_response__forward | 1e-06 | class_pauli_JD_budget__reverse | 5/9/0 | 0/0/0 |
| class_refined_response__forward | 1e-06 | class_refined_response__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__forward | 1e-08 | class_pauli_JD__forward | 3/10/1 | 0/0/0 |
| class_refined_response__forward | 1e-08 | class_pauli_JD__reverse | 4/10/0 | 0/0/0 |
| class_refined_response__forward | 1e-08 | class_pauli_JD_budget__forward | 6/8/0 | 0/0/0 |
| class_refined_response__forward | 1e-08 | class_pauli_JD_budget__reverse | 5/9/0 | 0/0/0 |
| class_refined_response__forward | 1e-08 | class_refined_response__reverse | 2/12/0 | 0/0/0 |
| class_refined_response__reverse | 0.0001 | class_pauli_JD__forward | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 0.0001 | class_pauli_JD__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 0.0001 | class_pauli_JD_budget__forward | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 0.0001 | class_pauli_JD_budget__reverse | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 0.0001 | class_refined_response__forward | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 1e-06 | class_pauli_JD__forward | 0/13/1 | 0/0/0 |
| class_refined_response__reverse | 1e-06 | class_pauli_JD__reverse | 0/13/1 | 0/0/0 |
| class_refined_response__reverse | 1e-06 | class_pauli_JD_budget__forward | 5/8/1 | 0/0/0 |
| class_refined_response__reverse | 1e-06 | class_pauli_JD_budget__reverse | 5/9/0 | 0/0/0 |
| class_refined_response__reverse | 1e-06 | class_refined_response__forward | 0/14/0 | 0/0/0 |
| class_refined_response__reverse | 1e-08 | class_pauli_JD__forward | 1/12/1 | 0/0/0 |
| class_refined_response__reverse | 1e-08 | class_pauli_JD__reverse | 2/12/0 | 0/0/0 |
| class_refined_response__reverse | 1e-08 | class_pauli_JD_budget__forward | 5/9/0 | 0/0/0 |
| class_refined_response__reverse | 1e-08 | class_pauli_JD_budget__reverse | 4/10/0 | 0/0/0 |
| class_refined_response__reverse | 1e-08 | class_refined_response__forward | 0/12/2 | 0/0/0 |

The uniform-fusion sensitivity retains each measured state and r, changing only exact word-rotation accounting via C_uniform(r)=r(2U-a-b)+a, where a,b are maximal pairwise-commuting prefix/suffix lengths. These are costs for the existing circuits, not simulations with additional steps. Detailed costs are in all_rows.csv and summary.json.

## Limits

Both orientations are reported individually. Reversal is equally good under first-order JD but its Strang error can differ; no error-informed orientation selector is presented. The matched-evaluation Pauli arm accepted only one or two moves per start in the original search, so it is a deliberately truncated best-improvement control, not evidence of general search-efficiency superiority. The larger-search arm is the substantive partition comparison. This controls search prescription and existing chosen orders, not separate reoptimization of a second-order objective. The panel comprises seven molecular families and correlated geometry/coefficient variants; no general arbitrary-state, operator-norm, or compiled hardware claim follows.
