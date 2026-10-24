# Exact restricted response windows

Complete enumeration certified 3 new nine-unit windows, each over 9! = 362880 feasible orders. Outside units and their positions remain fixed. All incumbents equal the independently verified lower bounds, with zero exact-integer gap and zero tolerances.

The objective is J_D = sum_w rho_w b_w(order)^2/(4 denominator^4), with integer commutator amplitudes and the exact rational values of supplied binary64 rho_w = sinc^2(gap_w/(2 pi denominator)) at T=1. This does not certify exact transcendental sinc weights, the full ordering domain, a molecular ground-state problem, or finite-step physical error. Complete-unit coefficients remain the authenticated rational model coefficients.

Each window is the contiguous nine-unit subset of the previously frozen HF J_D heuristic order with the most nonzero internal response pairs (earliest offset on ties). Selection uses no physical errors. Prior heuristic/current/larger-Pauli construction timings are inherited records; setup and exhaustive-enumeration costs are newly measured. No matched-classical-cost advantage is claimed.

| Case | Offset | Initial J_D | Certified J_D | Minimizers | Enumeration CPU s |
| --- | ---: | ---: | ---: | ---: | ---: |
| F2_r1p55 | 30 | 0.00539793555 | 0.00539793555 | 13440 | 0.82175 |
| H2O_scale1p10 | 35 | 0.000700124605 | 0.000454990265 | 44 | 1.39173 |
| Li2_r2p95 | 12 | 2.26402651e-05 | 3.7157666e-07 | 264 | 0.82927 |

All six schedules per input were frozen after certification and before physical evaluation. The grid includes first-order heuristic/certified/current/larger-search individual-Pauli controls and symmetric reuse of both heuristic and certified orders, at common 10U and 30U ideal Pauli-rotation caps. Exact central and boundary fusion determines each affordable step count. Symmetric orders were not optimized for a second-order objective.

| Formula | Rotation cap/U | Certified vs heuristic W/T/L | Floor exclusions | Median paired error ratio |
| --- | ---: | ---: | ---: | ---: |
| 1 | 10 | 2/1/0 | 0 | 0.6654703246395094 |
| 1 | 30 | 2/1/0 | 0 | 0.6649621622401821 |
| 2 | 10 | 2/1/0 | 1 | 0.9372387020227405 |
| 2 | 30 | 1/2/0 | 1 | 0.9372685774998746 |

Physical error is the squared projective state distance 1-|<exact|actual>|, evaluated stably from the orthogonal residual in the full Hilbert space. Number and joint (N_alpha,N_beta) leakage are retained; no sector projection is applied. Independent verification uses literal Pauli matrices, unfused circuits and a scaled-Taylor reference. Every floor and pairwise win/tie/loss decision is checked independently; all losses and floor-limited rows remain in the tables. The fixed floor and absolute comparison tolerance are 1e-12. Ratios exclude pairs with either error at/below the floor. Three related molecular cases are descriptive observations, not independent statistical replications.

Independent replay: 36 rows, maximum state difference 1.86190061e-15, maximum metric difference 1.99840144e-15. Each integer certificate was independently re-enumerated by DFS with complete Lehmer-rank coverage, and its full-order objective was checked. Independent floating checks do not constitute an interval or real-transcendental certificate.

The earlier four eight-unit windows on CH2O_s1.00 and nh3_coeff1p00 remain completed prior restricted evidence. Their retained original result and source are included under prior/ with hashes and provenance; they were not rerun. This study extends their domain size, without converting either study into a whole-order or physical-error optimum.
