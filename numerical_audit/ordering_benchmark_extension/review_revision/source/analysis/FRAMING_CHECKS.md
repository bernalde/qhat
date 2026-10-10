# Review framing checks

All statistics below are recomputed from the frozen production row files; no propagation or new order selection occurs. Original manifests and rows are unchanged.

Threshold levels (conserving lower / same / higher): {'same': 92, 'lower': 75, 'higher': 1}

| Target | Comparator | Lower / same / higher |
|---|---|---|
| 0.0001 | strang_class_lex | 0/14/0 |
| 0.0001 | strang_jw_magnitude | 0/14/0 |
| 0.0001 | strang_jw_signed | 4/10/0 |
| 0.0001 | tranter_v2_l1_strang | 3/11/0 |
| 1e-06 | strang_class_lex | 7/7/0 |
| 1e-06 | strang_jw_magnitude | 9/4/1 |
| 1e-06 | strang_jw_signed | 11/3/0 |
| 1e-06 | tranter_v2_l1_strang | 10/4/0 |
| 1e-08 | strang_class_lex | 6/8/0 |
| 1e-08 | strang_jw_magnitude | 7/7/0 |
| 1e-08 | strang_jw_signed | 9/5/0 |
| 1e-08 | tranter_v2_l1_strang | 9/5/0 |

The same-level fused rotation counts remain in threshold_levels.csv. They are exact costs of tested circuits, not resolved differences between minimum costs required to reach the target.

First-order candidate beats 840 of 840 fixed random outcomes. Random quartiles are in random_distributions.csv.

Per-case search depth, floor-aware formula-order comparisons and unweighted controls are in framing_checks.json.
