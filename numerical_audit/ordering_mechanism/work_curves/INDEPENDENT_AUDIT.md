# Independent finite-step work-curve audit

Auditor: the `ordering_gauge_experiment` subagent, 2026-10-09. This is a read-only review of the frozen source and saved results. No producer files were changed and no new propagation was required.

**Assessment: the reported comparisons and threshold aggregation are correct for these eight development fixtures.** The conserving refined candidate beats each of complete-class lexicographic, owned signed, JW signed, and JW magnitude at every tested step count, with eight wins and no ties or losses for each comparison. The tested counts are 1, 2, 5, 10, 20, 50, and 100. A tie means an absolute difference in the overlap-error metric of at most `1e-12`.

## Reproduction and checked provenance

The independent checker imports NumPy, not QHAT, and uses explicit exceptions, so its checks remain active under `python -O`:

```sh
/tmp/qhat-h2o-venv/bin/python -O /tmp/ordering-work-curves/check_work_curves.py \
  --run /tmp/ordering-work-curves/results \
  --candidate /tmp/ordering-diagonal-refinement/results \
  --baseline /tmp/followup-class-comparison \
  --output /tmp/ordering-work-curves/independent_provenance_aggregation.json
```

The output is already present; use a new output filename for a replay. It reports 240 new rows, 40 reused `r=100` rows, and 280 distinct `(case, method, r)` records. Independently recomputing the stable normalized, phase-insensitive metric from all 280 saved states agreed to a maximum absolute difference of `2.168404344971009e-18`.

For every case, the checker establishes exact equality of the rationalized Pauli coefficients and denominator, electron and qubit counts, affine basis, target hash, initial-state hash, and tensor hash across the current run and both source runs. The initial, exact, and dense-exact reference arrays agree bitwise. Every flattened ordered word/coefficient schedule matches its actual source. In particular, the refined zero-mask singleton representation of `class_lex` flattens to exactly the same circuit as the original complete-class lexicographic representation. Every arm emits the same target Pauli multiset once per step, and the recorded rotation count equals its size times `r`.

All 40 reused rows agree exactly with their source rows apart from the added source-results hash; the recorded hash is correct. All active orders belong to the five tested methods, and the planned-run grid is correct. The inherited inactive graph, old schedule-hash, and suffix fields identified in the earlier refinement audit have been removed here. The frozen candidate helper is byte-identical to the phase-two source. All 24 threshold records were independently reconstructed from the saved errors.

The manifest SHA256 is `0e31a33e15644133038cf678158d92f8ec52f44c5dd60b08d0b82127dfea95f1`; the independent checker SHA256 is `063731f3b250a55241afadd22735b773d8d3a6a9a1a5990e808250adb005c35c`.

## Verified work claim

At overlap error `E <= 1e-6`, the least **tested** step counts and costs are:

| Case | Candidate r | JW signed r | Candidate Pauli rotations | JW signed Pauli rotations | JW/candidate work |
|---|---:|---:|---:|---:|---:|
| CH2O, scale 1.00 | 10 | 50 | 2,510 | 12,550 | 5 |
| CH2O, scale 1.25 | 20 | 50 | 5,020 | 12,550 | 2.5 |
| H2S, scale 1.00 | 2 | 10 | 502 | 2,510 | 5 |
| H2S, scale 1.25 | 5 | 20 | 1,375 | 5,500 | 4 |
| HCN, scale 1.00 | 20 | 50 | 8,860 | 22,150 | 2.5 |
| HCN, scale 1.25 | 20 | 100 | 8,860 | 44,300 | 5 |
| SiH4, scale 1.00 | 2 | 5 | 1,678 | 4,195 | 2.5 |
| SiH4, scale 1.25 | 5 | 10 | 4,135 | 8,270 | 2 |

Thus the claim of 2–5 times fewer retained Pauli-word rotations is supported on this grid. It is not a claim about optimal untested integer step counts, monotonicity between grid points, compiled two-qubit gates, circuit depth, or cancellations introduced by circuit compilation. The protocol states those limits correctly. Classical ordering-search cost is separate and is not included in this table.

## Independent circuit replay versus this audit

I inspected the separate `independent_recheck.json` produced by the existing no-QHAT replay checker. It records independent replay of all 240 newly propagated circuits with a maximum state distance of `1.3812056679245397e-15` and maximum metric difference of `7.155734338404329e-17`. The present audit independently rechecks provenance, frozen schedules, saved-state metrics, and aggregation; it does not claim to have repeated that propagation a second time. The reused `r=100` circuits have prior replay evidence and were matched exactly to those sources here.

## Remaining implementation hardening and interpretation

No scientific or active-schema error was found in the current artifacts. The runner's reuse path records a source-results hash but does not itself reject incompatible baseline targets, reference arrays, or schedules before appending old rows. Those conditions were independently checked and hold for this run. Before using the driver with other sources, add explicit compatibility checks; hashing a source file alone does not establish compatibility. The runner also relies on its execution environment for imported QHAT dependency identity, so final artifact packaging should retain the pinned dependency provenance already used by the research workflow.

These are eight development fixtures from four molecules at two geometry scales. The new step counts test finite-step persistence of a previously selected rule; they do not turn the fixtures into held-out data, establish performance for arbitrary states, or isolate the benefit of the diagonal-response weights. The separate objective and partition ablations remain necessary for those interpretations.
