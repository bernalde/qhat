# Controlled search costs and declared synthetic scaling

This fixed panel retains 432 constructed orders and 864 full-space circuit comparisons: three retained molecular targets and five synthetic extended-Hubbard targets with 4–12 modes. The synthetic family extends the measured mode/term range; it supplies no larger-molecule validation. All repetitions and incomplete/out-of-budget records remain visible. 0 bounded-arm records exceeded their declared stand-alone CPU cap.

Conserving versus individual-Pauli first-order search has the following pooled wins/ties/losses: 0.5 CPU seconds: 37/0/59; 2.0 CPU seconds: 21/0/75. The complete table below includes both formulas/caps and every partition, initializer and objective comparison. These summaries pool the other declared factors and related repetitions; they do not establish a universal partition advantage or isolate conservation as a mechanism.

| Paired factor (left/right) | Formula | CPU cap | Wins/ties/losses | Median paired error ratio | Both within cap |
| --- | --- | --- | --- | --- | --- |
| conserving_vs_pauli | 1 | 0.5 | 37/0/59 | 1.19114 | 96/96 |
| conserving_vs_pauli | 1 | 2.0 | 21/0/75 | 1.56003 | 96/96 |
| conserving_vs_pauli | 2 | 0.5 | 51/6/39 | 0.968236 | 96/96 |
| conserving_vs_pauli | 2 | 2.0 | 39/6/51 | 1.08942 | 96/96 |
| placement_vs_cold | 1 | 0.5 | 66/18/12 | 0.565782 | 96/96 |
| placement_vs_cold | 1 | 2.0 | 60/18/18 | 0.81656 | 96/96 |
| placement_vs_cold | 2 | 0.5 | 56/18/22 | 0.864365 | 96/96 |
| placement_vs_cold | 2 | 2.0 | 57/18/21 | 0.851956 | 96/96 |
| JD_vs_J0 | 1 | 0.5 | 21/66/9 | 1 | 96/96 |
| JD_vs_J0 | 1 | 2.0 | 9/81/6 | 1 | 96/96 |
| JD_vs_J0 | 2 | 0.5 | 14/67/15 | 1 | 96/96 |
| JD_vs_J0 | 2 | 2.0 | 6/81/9 | 1 | 96/96 |

Eight factorial arms separate conserving/individual-Pauli partition, J_D/J_0 objective and cold/placement initialization. Caps are 0.5 and 2 process-CPU seconds, with three repetitions, cyclic arm order and symmetric warmups. [Search records](searches.csv) count common setup, charged warm initializer/conversion, method setup, search, positions, partial sweeps, rescoring, storage and actual overruns. [Shared costs](shared_setup.csv) record actual common/initializer execution once per input/repetition; charging an inherited initializer to several standalone methods does not mean the campaign recomputed it each time. [All warmups](warmups.csv) retain actual optimizer work and CPU failures. [Actual campaign construction](campaign_construction.csv) counts shared setup and initializer once for the warmup phase and each measured repetition, then adds actual local optimizer CPU; it excludes file serialization, import overhead and original input generation, which remain in the full process resource ledger. Original input generation remains separate from stored-input search setup.

[Paired factorial comparisons](paired_factorial.csv) distinguish partition, objective and initializer effects at each input/repetition/budget/formula. Aggregate counts are descriptive summaries of related pairs, not independent statistical replicates. Ratios are formed per pair before aggregation, excluding and counting only denominator errors at or below 1e-12. Losses remain in every count. The both-within-cap field identifies failures to satisfy a nominal matched budget; such rows are retained and never credited as successful in-budget comparisons. Incomplete sweeps remain eligible executed policies only under their explicitly retained stopping records.

[Strong controls](strong_controls.csv) include separately timed current placement and the 8-CPU-second individual-Pauli search, with both first-order and symmetric reuse. The larger Pauli control has an unequal classical budget. Comparisons changing formula do not isolate ordering. Symmetric controls reuse first-order objectives, with exact central/boundary fusion and actual affordable steps at 10U ideal rotations; they do not complete second-order-specific optimization.

Full leakage-preserving propagation, independent literal-Pauli surrogate scoring, an independently scaled Taylor reference and unfused replay validate every retained circuit. Preparation, search storage and circuit cache limits are declared before execution. The process peak RSS at each input end is cumulative for the worker, not an independently reset per-input peak. Sparse array payload estimates omit Python/container overhead; the OS memory ceiling and measured RSS remain separate controls.

Conclusions are limited to this implementation, target family, budget/stopping policy and declared scaling range. Quarter-core/shared-machine CPU observations do not establish dedicated-machine wall speedups, an asymptotic advantage, conservation as the accuracy mechanism, compiled resources or hardware results. The rejected two-level correction is unchanged; no failed holdout gate is retuned. Full absolute errors and sector leakage are in [all rows](all_rows.csv).
