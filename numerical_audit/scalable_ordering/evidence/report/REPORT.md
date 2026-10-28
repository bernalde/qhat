# Sparse molecular ordering: derived completed-study report

Independent verification: **passed**; 24 frozen orders, 96 declared circuits, eight references and 240 paired comparisons. Circuit statuses: {"complete": 96}.

This panel evaluates four new recorded-baseline geometries from three families. The CH4 family already appeared in the baseline. Results below retain each setting and all cubic/control losses, ties and censors; they do not establish a universal winner or novelty.

Scientific source commit: `a21a449ca519b588fc22a9c3f66a64dddc0e1e54`. Every consumed input has an exact byte count and SHA-256 in `report_tables.json`.

The orders were originally constructed and frozen at `a9e8f5047ca3ee5be68fa658abe25779ebd9d7f4`. A row-packaging correction changed only study.py before successful execution. The renderer checks the original/effective freezes and amendment: all 24 constructor hashes are identical, with no retraining or retuning and no ordering-error result from the failed attempt. That attempt and its cost remain in the evidence archive.

## Frozen protocol and scope

Train at time 1.0 on HF occupations with a 30U ideal-rotation cap, where U is the retained Pauli-word count. Freeze all orders before validation at times 0.5 and 1.5 with caps 10U and 30U. The searched first-response, individual-Pauli and sparse-cubic arms have the same nominal 32 CPU-second allowance, up to 3,000 candidate evaluations and five accepted moves per declared start (current placement, its reverse and canonical). Actual costs and incomplete searches appear below. Canonical and magnitude are cheap baselines. All arms use the symmetric formula.

One worker, quarter-core CPU, one native thread, 512 MiB RAM, no swap; start above 5.5 GiB available and pause below 5 GiB; phase ledger owns total allowances.

Exact rational H/S/Sdg/Rz/CX primitives on a nearest-neighbor line; HF prefix included.

Reference propagation and independent full-space checking are offline validation. The response score is not a finite-time accuracy bound.

## Molecules

| Case | Family | Spin modes | Active electrons | Pauli words U | Conserving units | Chemistry CPU s |
| --- | --- | --- | --- | --- | --- | --- |
| N2_r1p10 | N2 | 12 | 6 | 390 | 145 | 9.84246 |
| N2_r2p00 | N2 | 12 | 6 | 690 | 175 | 0.338343 |
| CO_r1p13 | CO | 12 | 6 | 770 | 195 | 0.335999 |
| CH4_r1p09 | CH4 | 12 | 6 | 1722 | 351 | 0.904552 |

All cases use STO-3G and the registry's frozen-core/active windows. The two N2 geometries and repeated time/cap settings are dependent observations.

## Standalone order construction

CPU includes authentication, setup, initialization, search and returned layout. A null stored score means no completed training score; the separate final-order objective audit does not retroactively make that constructor a completed search. Completed seeds, candidate work and move rescoring below come from the recorded statistics. Cubic candidate work uses incremental_candidates; first/Pauli use candidate_evaluations. A dash means the counter is unknown, including unrecorded rescore attempts; it does not mean zero. Unsearched baselines have zero search counters, while their initializer work remains in CPU charges. Configured allowances alone do not establish completed search.

| Case | Arm | CPU s | Wall s | Status | Evaluations | Completed seeds | Candidate work | Rescores completed/attempted | Accepted moves | Stored score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N2_r1p10 | canonical | 0.00972318 | 0.00975425 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r1p10 | magnitude | 0.00847475 | 0.00843392 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r1p10 | current | 0.115244 | 0.488198 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r1p10 | first | 0.0894291 | 0.319444 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00193324 |
| N2_r1p10 | pauli | 0.380421 | 1.50301 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00145017 |
| N2_r1p10 | cubic | 3.37623 | 13.4903 | evaluation_limit | 3000 | 3 | 2996 incremental | 1/1 | 1 | 1.46139e-08 |
| N2_r2p00 | canonical | 0.107125 | 0.403902 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r2p00 | magnitude | 0.108997 | 0.485435 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r2p00 | current | 0.0894081 | 0.322939 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| N2_r2p00 | first | 0.170758 | 0.700748 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00150263 |
| N2_r2p00 | pauli | 1.21618 | 4.88583 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00596442 |
| N2_r2p00 | cubic | 5.83187 | 23.3046 | evaluation_limit | 3000 | 3 | 2996 incremental | 1/1 | 1 | 8.83905e-09 |
| CO_r1p13 | canonical | 0.110804 | 0.473408 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CO_r1p13 | magnitude | 0.0380118 | 0.186616 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CO_r1p13 | current | 0.110753 | 0.489175 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CO_r1p13 | first | 0.215553 | 0.816804 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.0104106 |
| CO_r1p13 | pauli | 1.42707 | 5.71017 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.0159292 |
| CO_r1p13 | cubic | 6.4049 | 25.6805 | evaluation_limit | 3000 | 3 | 2996 incremental | 1/1 | 1 | 1.9131e-08 |
| CH4_r1p09 | canonical | 0.0625922 | 0.210631 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CH4_r1p09 | magnitude | 0.176954 | 0.777714 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CH4_r1p09 | current | 0.425067 | 1.69994 | complete | 0 | 0 | 0 (unsearched) | 0/0 | 0 | — |
| CH4_r1p09 | first | 0.917392 | 3.61437 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00225795 |
| CH4_r1p09 | pauli | 7.03327 | 28.1099 | evaluation_limit | 3000 | 3 | 2997 candidates | 1/— | 1 | 0.00225191 |
| CH4_r1p09 | cubic | 29.9851 | 119.985 | cpu_limit | 2611 | 3 | 2605 incremental | 3/3 | 3 | 2.66672e-09 |

## Cubic versus controls: all 16 matched settings per control

A win means cubic has smaller state error by more than 1e-12; losses reverse that comparison. Differences within 1e-12 are ties. Ratios are cubic/control, formed within each setting before aggregation. Values below one favor cubic for the stated metric. Undefined zero-denominator ratios are excluded from medians and counted explicitly. Error ratios near the numerical floor are descriptive only.

| Control | Wins | Losses | Ties | Censors | Median error ratio (n) | Median CX ratio (n) | Median depth ratio (n) | Median construction CPU ratio (n) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | 16 | 0 | 0 | 0 | 0.269925 (16) | 1.04656 (16) | 1.03636 (16) | 4.19141 (16) |
| magnitude | 16 | 0 | 0 | 0 | 0.315652 (16) | 0.999621 (16) | 0.996612 (16) | 4.01517 (16) |
| current | 15 | 0 | 1 | 0 | 0.731055 (16) | 1.00013 (16) | 1.00031 (16) | 4.07133 (16) |
| first | 11 | 0 | 5 | 0 | 0.983491 (16) | 1.00013 (16) | 1.00034 (16) | 3.92443 (16) |
| pauli | 16 | 0 | 0 | 0 | 0.759643 (16) | 0.996461 (16) | 0.995643 (16) | 2.42383 (16) |

| Control | Complete floor pairs | Error ratios undefined | CX ratios undefined | Depth ratios undefined | Median error ratio above floor (n) |
| --- | --- | --- | --- | --- | --- |
| canonical | 0 | {} | {} | {} | 0.269925 (16) |
| magnitude | 0 | {} | {} | {} | 0.315652 (16) |
| current | 0 | {} | {} | {} | 0.731055 (16) |
| first | 0 | {} | {} | {} | 0.983491 (16) |
| pauli | 0 | {} | {} | {} | 0.759643 (16) |

The construction CPU ratio adds standalone ordering and that circuit's compilation. The same order is reused at four validation settings; the campaign totals charge it once.

## Raw six-arm settings

Errors are 1-|overlap|. CX and ASAP depth come from actual compiled primitive DAGs, including HF preparation. The listed lowest-error set includes every completed arm within 1e-12 of the setting's minimum; a censor prevents a complete six-arm ranking. Partial circuits keep their reason and partial costs in `report_tables.json`.

### N2_r1p10: time 0.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 3876 | 6.95579e-08 | 140880 | 135210 | 0.00972318 | 0.762111 |
| magnitude | complete | 5 | 3876 | 3.2083e-08 | 150840 | 142794 | 0.00847475 | 0.844063 |
| current | complete | 5 | 3876 | 1.79538e-08 | 150560 | 141569 | 0.115244 | 0.826745 |
| first | complete | 5 | 3876 | 1.46829e-08 | 150560 | 141529 | 0.0894291 | 0.884297 |
| pauli | complete | 5 | 3891 | 1.70315e-08 | 151430 | 142434 | 0.380421 | 0.9283 |
| cubic | complete | 5 | 3879 | 1.32413e-08 | 150792 | 141731 | 3.37623 | 0.834211 |

### N2_r1p10: time 0.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 11626 | 8.56328e-10 | 422640 | 405620 | 0.00972318 | 2.72227 |
| magnitude | complete | 15 | 11626 | 3.94974e-10 | 452520 | 428374 | 0.00847475 | 2.64044 |
| current | complete | 15 | 11626 | 2.21556e-10 | 451680 | 424699 | 0.115244 | 2.91251 |
| first | complete | 15 | 11626 | 1.81233e-10 | 451680 | 424579 | 0.0894291 | 2.75648 |
| pauli | complete | 15 | 11671 | 2.10191e-10 | 454290 | 427294 | 0.380421 | 2.80736 |
| cubic | complete | 15 | 11629 | 1.63549e-10 | 451912 | 424741 | 3.37623 | 2.90766 |

### N2_r1p10: time 1.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 3876 | 1.58628e-05 | 140880 | 135210 | 0.00972318 | 0.852628 |
| magnitude | complete | 5 | 3876 | 9.44481e-06 | 150840 | 142794 | 0.00847475 | 0.883296 |
| current | complete | 5 | 3876 | 5.76575e-06 | 150560 | 141569 | 0.115244 | 0.747002 |
| first | complete | 5 | 3876 | 4.57059e-06 | 150560 | 141529 | 0.0894291 | 0.888255 |
| pauli | complete | 5 | 3891 | 5.45832e-06 | 151430 | 142434 | 0.380421 | 0.893839 |
| cubic | complete | 5 | 3879 | 4.44879e-06 | 150792 | 141731 | 3.37623 | 0.838266 |

### N2_r1p10: time 1.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 11626 | 1.91674e-07 | 422640 | 405620 | 0.00972318 | 2.62973 |
| magnitude | complete | 15 | 11626 | 1.14598e-07 | 452520 | 428374 | 0.00847475 | 2.90664 |
| current | complete | 15 | 11626 | 7.0657e-08 | 451680 | 424699 | 0.115244 | 2.89451 |
| first | complete | 15 | 11626 | 5.6006e-08 | 451680 | 424579 | 0.0894291 | 2.69057 |
| pauli | complete | 15 | 11671 | 6.69511e-08 | 454290 | 427294 | 0.380421 | 2.84518 |
| cubic | complete | 15 | 11629 | 5.49153e-08 | 451912 | 424741 | 3.37623 | 2.93382 |

### N2_r2p00: time 0.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): first, cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 6876 | 8.95407e-09 | 190160 | 182345 | 0.107125 | 1.27386 |
| magnitude | complete | 5 | 6864 | 1.06806e-08 | 199096 | 188884 | 0.108997 | 1.17759 |
| current | complete | 5 | 6876 | 6.19017e-09 | 199400 | 188439 | 0.0894081 | 1.37286 |
| first | complete | 5 | 6876 | 4.48531e-09 | 199400 | 188424 | 0.170758 | 1.18592 |
| pauli | complete | 5 | 6891 | 6.04632e-09 | 200290 | 189504 | 1.21618 | 1.25795 |
| cubic | complete | 5 | 6876 | 4.48531e-09 | 199400 | 188424 | 5.83187 | 1.33758 |

### N2_r2p00: time 0.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): first, cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 20626 | 1.10383e-10 | 570480 | 547025 | 0.107125 | 4.18643 |
| magnitude | complete | 15 | 20584 | 1.31624e-10 | 596856 | 566204 | 0.108997 | 4.04156 |
| current | complete | 15 | 20626 | 7.63173e-11 | 597980 | 565089 | 0.0894081 | 4.15445 |
| first | complete | 15 | 20626 | 5.5283e-11 | 597980 | 565044 | 0.170758 | 4.23137 |
| pauli | complete | 15 | 20671 | 7.4542e-11 | 600650 | 568284 | 1.21618 | 4.40629 |
| cubic | complete | 15 | 20626 | 5.5283e-11 | 597980 | 565044 | 5.83187 | 4.16039 |

### N2_r2p00: time 1.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): first, cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 6876 | 2.95557e-06 | 190160 | 182345 | 0.107125 | 1.27101 |
| magnitude | complete | 5 | 6864 | 2.13514e-06 | 199096 | 188884 | 0.108997 | 1.36432 |
| current | complete | 5 | 6876 | 2.16555e-06 | 199400 | 188439 | 0.0894081 | 1.17076 |
| first | complete | 5 | 6876 | 9.41193e-07 | 199400 | 188424 | 0.170758 | 1.27727 |
| pauli | complete | 5 | 6891 | 2.10636e-06 | 200290 | 189504 | 1.21618 | 1.3275 |
| cubic | complete | 5 | 6876 | 9.41193e-07 | 199400 | 188424 | 5.83187 | 1.26204 |

### N2_r2p00: time 1.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): first, cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 20626 | 3.61408e-08 | 570480 | 547025 | 0.107125 | 4.10014 |
| magnitude | complete | 15 | 20584 | 2.60842e-08 | 596856 | 566204 | 0.108997 | 4.41388 |
| current | complete | 15 | 20626 | 2.65984e-08 | 597980 | 565089 | 0.0894081 | 4.15902 |
| first | complete | 15 | 20626 | 1.15091e-08 | 597980 | 565044 | 0.170758 | 4.06308 |
| pauli | complete | 15 | 20671 | 2.58743e-08 | 600650 | 568284 | 1.21618 | 4.24367 |
| cubic | complete | 15 | 20626 | 1.15091e-08 | 597980 | 565044 | 5.83187 | 4.18382 |

### CO_r1p13: time 0.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 7676 | 4.97519e-08 | 211440 | 202865 | 0.110804 | 1.56075 |
| magnitude | complete | 5 | 7676 | 6.10071e-08 | 221720 | 211694 | 0.0380118 | 1.71606 |
| current | complete | 5 | 7676 | 3.46458e-08 | 220772 | 210789 | 0.110753 | 1.51825 |
| first | complete | 5 | 7676 | 2.53311e-08 | 220772 | 210824 | 0.215553 | 1.49795 |
| pauli | complete | 5 | 7691 | 3.33999e-08 | 221462 | 211654 | 1.42707 | 1.5897 |
| cubic | complete | 5 | 7679 | 1.29054e-08 | 220938 | 210965 | 6.4049 | 1.52984 |

### CO_r1p13: time 0.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 23026 | 6.125e-10 | 634320 | 608585 | 0.110804 | 4.85736 |
| magnitude | complete | 15 | 23026 | 7.51017e-10 | 665160 | 635074 | 0.0380118 | 5.0524 |
| current | complete | 15 | 23026 | 4.2676e-10 | 662312 | 632359 | 0.110753 | 5.00632 |
| first | complete | 15 | 23026 | 3.12178e-10 | 662312 | 632464 | 0.215553 | 4.87576 |
| pauli | complete | 15 | 23071 | 4.11402e-10 | 664382 | 634954 | 1.42707 | 4.98447 |
| cubic | complete | 15 | 23029 | 1.58638e-10 | 662478 | 632555 | 6.4049 | 4.88339 |

### CO_r1p13: time 1.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 7676 | 1.20308e-05 | 211440 | 202865 | 0.110804 | 1.51516 |
| magnitude | complete | 5 | 7676 | 1.55659e-05 | 221720 | 211694 | 0.0380118 | 1.76665 |
| current | complete | 5 | 7676 | 9.43234e-06 | 220772 | 210789 | 0.110753 | 1.40687 |
| first | complete | 5 | 7676 | 7.07861e-06 | 220772 | 210824 | 0.215553 | 1.42942 |
| pauli | complete | 5 | 7691 | 9.26125e-06 | 221462 | 211654 | 1.42707 | 1.4521 |
| cubic | complete | 5 | 7679 | 3.40106e-06 | 220938 | 210965 | 6.4049 | 1.44145 |

### CO_r1p13: time 1.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 23026 | 1.45722e-07 | 634320 | 608585 | 0.110804 | 4.69245 |
| magnitude | complete | 15 | 23026 | 1.90216e-07 | 665160 | 635074 | 0.0380118 | 4.81778 |
| current | complete | 15 | 23026 | 1.15569e-07 | 662312 | 632359 | 0.110753 | 4.80344 |
| first | complete | 15 | 23026 | 8.68288e-08 | 662312 | 632464 | 0.215553 | 4.74301 |
| pauli | complete | 15 | 23071 | 1.13453e-07 | 664382 | 634954 | 1.42707 | 4.91276 |
| cubic | complete | 15 | 23029 | 4.1578e-08 | 662478 | 632555 | 6.4049 | 5.26922 |

### CH4_r1p09: time 0.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 17196 | 1.63391e-08 | 480640 | 460765 | 0.0625922 | 3.56175 |
| magnitude | complete | 5 | 17112 | 1.60772e-08 | 490168 | 468485 | 0.176954 | 3.40315 |
| current | complete | 5 | 17196 | 2.16615e-09 | 489952 | 466779 | 0.425067 | 3.26926 |
| first | complete | 5 | 17196 | 2.16615e-09 | 489952 | 466767 | 0.917392 | 3.69979 |
| pauli | complete | 5 | 17211 | 2.25685e-09 | 490122 | 467719 | 7.03327 | 3.2315 |
| cubic | complete | 5 | 17196 | 2.13684e-09 | 489952 | 466924 | 29.9851 | 3.05913 |

### CH4_r1p09: time 0.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): current, first, cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 51586 | 2.01353e-10 | 1.44192e+06 | 1.38228e+06 | 0.0625922 | 11.2338 |
| magnitude | complete | 15 | 51292 | 1.98123e-10 | 1.46981e+06 | 1.40484e+06 | 0.176954 | 10.9765 |
| current | complete | 15 | 51586 | 2.6706e-11 | 1.46985e+06 | 1.40033e+06 | 0.425067 | 10.7575 |
| first | complete | 15 | 51586 | 2.6706e-11 | 1.46985e+06 | 1.40029e+06 | 0.917392 | 11.0282 |
| pauli | complete | 15 | 51631 | 2.78307e-11 | 1.47036e+06 | 1.40315e+06 | 7.03327 | 11.4892 |
| cubic | complete | 15 | 51586 | 2.63443e-11 | 1.46985e+06 | 1.40076e+06 | 29.9851 | 11.0876 |

### CH4_r1p09: time 1.5, cap 10U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 5 | 17196 | 6.47444e-06 | 480640 | 460765 | 0.0625922 | 3.37474 |
| magnitude | complete | 5 | 17112 | 6.33366e-06 | 490168 | 468485 | 0.176954 | 3.75218 |
| current | complete | 5 | 17196 | 8.83068e-07 | 489952 | 466779 | 0.425067 | 3.28416 |
| first | complete | 5 | 17196 | 8.83068e-07 | 489952 | 466767 | 0.917392 | 3.4481 |
| pauli | complete | 5 | 17211 | 9.2192e-07 | 490122 | 467719 | 7.03327 | 4.0905 |
| cubic | complete | 5 | 17196 | 8.74128e-07 | 489952 | 466924 | 29.9851 | 3.56786 |

### CH4_r1p09: time 1.5, cap 30U

Lowest error among completed arms (1e-12 ambiguity): cubic. Censored arms: none.

| Arm | Status | Steps | Ideal rotations | Error | CX | ASAP depth | Ordering CPU s | Compile CPU s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | complete | 15 | 51586 | 7.86739e-08 | 1.44192e+06 | 1.38228e+06 | 0.0625922 | 11.8824 |
| magnitude | complete | 15 | 51292 | 7.69601e-08 | 1.46981e+06 | 1.40484e+06 | 0.176954 | 11.7846 |
| current | complete | 15 | 51586 | 1.08018e-08 | 1.46985e+06 | 1.40033e+06 | 0.425067 | 11.5887 |
| first | complete | 15 | 51586 | 1.08018e-08 | 1.46985e+06 | 1.40029e+06 | 0.917392 | 10.9877 |
| pauli | complete | 15 | 51631 | 1.13004e-08 | 1.47036e+06 | 1.40315e+06 | 7.03327 | 12.2166 |
| cubic | complete | 15 | 51586 | 1.06916e-08 | 1.46985e+06 | 1.40076e+06 | 29.9851 | 12.1371 |

## Every cubic/control outcome

| Case | Time | Cap U | Control | Outcome | Error ratio | CX ratio | Depth ratio | Floor pair |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N2_r1p10 | 0.5 | 10 | canonical | win | 0.190364 | 1.07036 | 1.04823 | no |
| N2_r1p10 | 0.5 | 30 | canonical | win | 0.190989 | 1.06926 | 1.04714 | no |
| N2_r1p10 | 1.5 | 10 | canonical | win | 0.280455 | 1.07036 | 1.04823 | no |
| N2_r1p10 | 1.5 | 30 | canonical | win | 0.286504 | 1.06926 | 1.04714 | no |
| N2_r2p00 | 0.5 | 10 | canonical | win | 0.500924 | 1.04859 | 1.03334 | no |
| N2_r2p00 | 0.5 | 30 | canonical | win | 0.500828 | 1.04821 | 1.03294 | no |
| N2_r2p00 | 1.5 | 10 | canonical | win | 0.318448 | 1.04859 | 1.03334 | no |
| N2_r2p00 | 1.5 | 30 | canonical | win | 0.318452 | 1.04821 | 1.03294 | no |
| CO_r1p13 | 0.5 | 10 | canonical | win | 0.259394 | 1.04492 | 1.03993 | no |
| CO_r1p13 | 0.5 | 30 | canonical | win | 0.259001 | 1.04439 | 1.03939 | no |
| CO_r1p13 | 1.5 | 10 | canonical | win | 0.282697 | 1.04492 | 1.03993 | no |
| CO_r1p13 | 1.5 | 30 | canonical | win | 0.285324 | 1.04439 | 1.03939 | no |
| CH4_r1p09 | 0.5 | 10 | canonical | win | 0.130781 | 1.01937 | 1.01337 | no |
| CH4_r1p09 | 0.5 | 30 | canonical | win | 0.130836 | 1.01937 | 1.01337 | no |
| CH4_r1p09 | 1.5 | 10 | canonical | win | 0.135012 | 1.01937 | 1.01337 | no |
| CH4_r1p09 | 1.5 | 30 | canonical | win | 0.135897 | 1.01937 | 1.01337 | no |
| N2_r1p10 | 0.5 | 10 | magnitude | win | 0.412721 | 0.999682 | 0.992556 | no |
| N2_r1p10 | 0.5 | 30 | magnitude | win | 0.414075 | 0.998656 | 0.991519 | no |
| N2_r1p10 | 1.5 | 10 | magnitude | win | 0.47103 | 0.999682 | 0.992556 | no |
| N2_r1p10 | 1.5 | 30 | magnitude | win | 0.479201 | 0.998656 | 0.991519 | no |
| N2_r2p00 | 0.5 | 10 | magnitude | win | 0.419948 | 1.00153 | 0.997565 | no |
| N2_r2p00 | 0.5 | 30 | magnitude | win | 0.420007 | 1.00188 | 0.997951 | no |
| N2_r2p00 | 1.5 | 10 | magnitude | win | 0.440811 | 1.00153 | 0.997565 | no |
| N2_r2p00 | 1.5 | 30 | magnitude | win | 0.441229 | 1.00188 | 0.997951 | no |
| CO_r1p13 | 0.5 | 10 | magnitude | win | 0.211539 | 0.996473 | 0.996556 | no |
| CO_r1p13 | 0.5 | 30 | magnitude | win | 0.211231 | 0.995968 | 0.996034 | no |
| CO_r1p13 | 1.5 | 10 | magnitude | win | 0.218495 | 0.996473 | 0.996556 | no |
| CO_r1p13 | 1.5 | 30 | magnitude | win | 0.218583 | 0.995968 | 0.996034 | no |
| CH4_r1p09 | 0.5 | 10 | magnitude | win | 0.132911 | 0.999559 | 0.996668 | no |
| CH4_r1p09 | 0.5 | 30 | magnitude | win | 0.13297 | 1.00003 | 0.997102 | no |
| CH4_r1p09 | 1.5 | 10 | magnitude | win | 0.138013 | 0.999559 | 0.996668 | no |
| CH4_r1p09 | 1.5 | 30 | magnitude | win | 0.138924 | 1.00003 | 0.997102 | no |
| N2_r1p10 | 0.5 | 10 | current | win | 0.737523 | 1.00154 | 1.00114 | no |
| N2_r1p10 | 0.5 | 30 | current | win | 0.738184 | 1.00051 | 1.0001 | no |
| N2_r1p10 | 1.5 | 10 | current | win | 0.771589 | 1.00154 | 1.00114 | no |
| N2_r1p10 | 1.5 | 30 | current | win | 0.777209 | 1.00051 | 1.0001 | no |
| N2_r2p00 | 0.5 | 10 | current | win | 0.724586 | 1 | 0.99992 | no |
| N2_r2p00 | 0.5 | 30 | current | win | 0.724384 | 1 | 0.99992 | no |
| N2_r2p00 | 1.5 | 10 | current | win | 0.43462 | 1 | 0.99992 | no |
| N2_r2p00 | 1.5 | 30 | current | win | 0.432699 | 1 | 0.99992 | no |
| CO_r1p13 | 0.5 | 10 | current | win | 0.372494 | 1.00075 | 1.00083 | no |
| CO_r1p13 | 0.5 | 30 | current | win | 0.371726 | 1.00025 | 1.00031 | no |
| CO_r1p13 | 1.5 | 10 | current | win | 0.360575 | 1.00075 | 1.00083 | no |
| CO_r1p13 | 1.5 | 30 | current | win | 0.359767 | 1.00025 | 1.00031 | no |
| CH4_r1p09 | 0.5 | 10 | current | win | 0.986469 | 1 | 1.00031 | no |
| CH4_r1p09 | 0.5 | 30 | current | tie | 0.986457 | 1 | 1.00031 | no |
| CH4_r1p09 | 1.5 | 10 | current | win | 0.989876 | 1 | 1.00031 | no |
| CH4_r1p09 | 1.5 | 30 | current | win | 0.989798 | 1 | 1.00031 | no |
| N2_r1p10 | 0.5 | 10 | first | win | 0.901822 | 1.00154 | 1.00143 | no |
| N2_r1p10 | 0.5 | 30 | first | win | 0.902422 | 1.00051 | 1.00038 | no |
| N2_r1p10 | 1.5 | 10 | first | win | 0.973352 | 1.00154 | 1.00143 | no |
| N2_r1p10 | 1.5 | 30 | first | win | 0.980524 | 1.00051 | 1.00038 | no |
| N2_r2p00 | 0.5 | 10 | first | tie | 1 | 1 | 1 | no |
| N2_r2p00 | 0.5 | 30 | first | tie | 1 | 1 | 1 | no |
| N2_r2p00 | 1.5 | 10 | first | tie | 1 | 1 | 1 | no |
| N2_r2p00 | 1.5 | 30 | first | tie | 1 | 1 | 1 | no |
| CO_r1p13 | 0.5 | 10 | first | win | 0.509468 | 1.00075 | 1.00067 | no |
| CO_r1p13 | 0.5 | 30 | first | win | 0.508166 | 1.00025 | 1.00014 | no |
| CO_r1p13 | 1.5 | 10 | first | win | 0.480471 | 1.00075 | 1.00067 | no |
| CO_r1p13 | 1.5 | 30 | first | win | 0.47885 | 1.00025 | 1.00014 | no |
| CH4_r1p09 | 0.5 | 10 | first | win | 0.986469 | 1 | 1.00034 | no |
| CH4_r1p09 | 0.5 | 30 | first | tie | 0.986457 | 1 | 1.00034 | no |
| CH4_r1p09 | 1.5 | 10 | first | win | 0.989876 | 1 | 1.00034 | no |
| CH4_r1p09 | 1.5 | 30 | first | win | 0.989798 | 1 | 1.00034 | no |
| N2_r1p10 | 0.5 | 10 | pauli | win | 0.777461 | 0.995787 | 0.995064 | no |
| N2_r1p10 | 0.5 | 30 | pauli | win | 0.778095 | 0.994765 | 0.994025 | no |
| N2_r1p10 | 1.5 | 10 | pauli | win | 0.815047 | 0.995787 | 0.995064 | no |
| N2_r1p10 | 1.5 | 30 | pauli | win | 0.82023 | 0.994765 | 0.994025 | no |
| N2_r2p00 | 0.5 | 10 | pauli | win | 0.741825 | 0.995556 | 0.994301 | no |
| N2_r2p00 | 0.5 | 30 | pauli | win | 0.741635 | 0.995555 | 0.994299 | no |
| N2_r2p00 | 1.5 | 10 | pauli | win | 0.446835 | 0.995556 | 0.994301 | no |
| N2_r2p00 | 1.5 | 30 | pauli | win | 0.444809 | 0.995555 | 0.994299 | no |
| CO_r1p13 | 0.5 | 10 | pauli | win | 0.386389 | 0.997634 | 0.996745 | no |
| CO_r1p13 | 0.5 | 30 | pauli | win | 0.385604 | 0.997134 | 0.996222 | no |
| CO_r1p13 | 1.5 | 10 | pauli | win | 0.367236 | 0.997634 | 0.996745 | no |
| CO_r1p13 | 1.5 | 30 | pauli | win | 0.366479 | 0.997134 | 0.996222 | no |
| CH4_r1p09 | 0.5 | 10 | pauli | win | 0.946824 | 0.999653 | 0.9983 | no |
| CH4_r1p09 | 0.5 | 30 | pauli | win | 0.946592 | 0.999653 | 0.9983 | no |
| CH4_r1p09 | 1.5 | 10 | pauli | win | 0.948161 | 0.999653 | 0.9983 | no |
| CH4_r1p09 | 1.5 | 30 | pauli | win | 0.946122 | 0.999653 | 0.9983 | no |

## Recorded CPU charges

| Recorded component | CPU s |
| --- | --- |
| chemistry_model_generation_CPU_seconds | 11.4214 |
| frozen_prepare_CPU_seconds | 58.9567 |
| standalone_constructor_CPU_seconds | 58.4113 |
| compilation_CPU_seconds | 365.688 |
| primitive_state_validation_CPU_seconds | 632.305 |
| reference_CPU_seconds | 0.464182 |
| independent_reference_CPU_seconds | 68.8962 |

Recorded successful-stage charges only; not total task CPU. Failed attempts, imports/startup, geometry scan, external review and independent replay costs beyond the recorded references require the separate task ledger. Each constructor was charged once in this total.

Frozen preparation CPU includes its constructor work and must not be added again to the standalone-constructor sum. No whole-task speed or memory claim is derived from these totals.

## Limits

- Four 12-mode, six-active-electron cases in three families. Two N2 geometries and repeated times/budgets are dependent observations; these counts are not independent sample sizes.
- The geometry audit concerns the pinned recorded baseline only. CH4 was already a baseline family; new geometries do not establish family generalization.
- All 24 orders were frozen before validation. Exact references and full 4096-state checking are offline validation costs; only training avoids Hilbert-space-sized arrays.
- Equal nominal search CPU limits do not imply equal consumed CPU, useful evaluations or construction cost. CPU-limited fallback orders, missing scores and censors remain reported.
- The cubic response is a heuristic, not a finite-error certificate. No globally best ordering or literature novelty is claimed.
- CX and depth are costs of this line compiler and its declared cleanup, including HF preparation; no native hardware, synthesis-optimality, correlated-state preparation, basis-error or SCF-error claim follows.
- Each ratio is computed within a matched case/time/cap before taking its median. Zero control values produce explicit undefined ratios and are excluded from that median. Numerical-floor pairs remain in outcome counts and are identified separately.
- Independent verification is numerical, not an interval or Lean certificate. This renderer checks the saved verification result and its inventories; it is not a formal Astra review.
