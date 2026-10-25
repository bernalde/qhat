# Compiled circuit resources at fixed accuracy

This is the complete prospectively frozen compiled-resource study.

Compiler `line-pauli-exact-angle-0.1.0` lowers to H/S/Sdg/Rz/CX on a nearest-neighbor line, with exact rational Rz angles, no ancillas, restored routing and only globally adjacent cancellation. This prototype supplies ideal compiled counts, not hardware execution or fault-tolerant synthesis.

All 48 order/formula curves and 384 candidates at steps [1, 2, 4, 8, 16, 32, 64, 128] were frozen before compiled-state errors. Orders are the eight fixed 30U-trained arms from numerical issue8; no order is retuned on compiled resources or these errors. There are 384 feasible circuits and 0 compiler-censored candidates.

The projective state-error targets are0.01 and0.001. Values within2e-12 of a target remain ambiguous and cannot establish attainment. Each curve selects the qualifying grid point with minimum CX, then unit-duration ASAP depth, total gates and step count; error need not be monotone. Unmet and censored cases remain explicit. Exact primitive traces include HF preparation and all allowed boundary cancellations.

| Target | Cubic insertion versus | Formula | Both met | CX W/T/L | Depth W/T/L | Cubic-only met | Control-only met | Both unmet |
| ---: | --- | ---: | ---: | --- | --- | ---: | ---: | ---: |
| 0.01 | current | 1 | 6 | 0/0/6 | 0/0/6 | 0 | 0 | 0 |
| 0.01 | first_JD_sameT | 1 | 6 | 0/0/6 | 0/0/6 | 0 | 0 | 0 |
| 0.01 | pauli_JD_large_sameT | 1 | 6 | 0/1/5 | 1/0/5 | 0 | 0 | 0 |
| 0.01 | current | 2 | 6 | 0/4/2 | 2/2/2 | 0 | 0 | 0 |
| 0.01 | first_JD_sameT | 2 | 6 | 0/4/2 | 2/2/2 | 0 | 0 | 0 |
| 0.01 | pauli_JD_large_sameT | 2 | 6 | 2/0/4 | 2/0/4 | 0 | 0 | 0 |
| 0.01 | cubic_orientation | 2 | 6 | 0/6/0 | 4/2/0 | 0 | 0 | 0 |
| 0.001 | current | 1 | 6 | 1/1/4 | 2/0/4 | 0 | 0 | 0 |
| 0.001 | first_JD_sameT | 1 | 6 | 1/1/4 | 2/0/4 | 0 | 0 | 0 |
| 0.001 | pauli_JD_large_sameT | 1 | 6 | 2/0/4 | 2/0/4 | 0 | 0 | 0 |
| 0.001 | current | 2 | 6 | 0/4/2 | 2/2/2 | 0 | 0 | 0 |
| 0.001 | first_JD_sameT | 2 | 6 | 0/4/2 | 2/2/2 | 0 | 0 | 0 |
| 0.001 | pauli_JD_large_sameT | 2 | 6 | 3/0/3 | 3/0/3 | 0 | 0 | 0 |
| 0.001 | cubic_orientation | 2 | 6 | 0/6/0 | 4/2/0 | 0 | 0 | 0 |

`all_candidates.csv` retains the entire error-versus-resource grid; `targets.csv` retains every attainment decision and selected one-circuit cost; `comparisons.csv` retains every control pair. `curves.csv` records each method’s largest feasible grid step. The reported depth belongs to its CX-minimizing selected circuit; it is not a separate depth optimum. A CX tie is distinct from an ASAP-depth tie. Both-met count/depth outcomes do not silently turn an unmet curve into a resource win.

`constructors.csv`, `curves.csv` and `COSTS.json` separate lowering, routing/cancellation, resource counting, exact-trace serialization/hashing, bookkeeping, source/upstream authentication, shared model setup and full-grid physical screening. Inherited ordering charges are standalone scenarios from issue8, not newly executed searches or a sum of campaign work. Per-row reference timings repeat the same shared solve. One chosen circuit cost is not the cost of screening the entire grid for target attainment.

Independent verification checks literal sparse primitive replay against unfused Pauli formulas and a separate scaled-Taylor target action, with state tolerance2e-10 and metric tolerance2e-12. Actual maximum differences: state3.65645654e-15, metric3.55271368e-15. An independent two-pass dependency DAG checks counts/depth. Target decisions, ambiguity and censoring agree. Exact canonical trace reproduction reuses the pinned producer compiler and is labelled separately from the independent replay/DAG oracle.

Results apply to these related finite molecular inputs, frozen times, exact compiler and supplied grid. They do not establish global compilation/search optimality, matched classical CPU, conservation as the mechanism, larger-molecule validity, hardware depth/timing, noise performance or a Lean proof.
