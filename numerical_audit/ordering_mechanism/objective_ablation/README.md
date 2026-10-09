# Objective and partition ablation

These sixteen additional evolutions qualify the interpretation of the promising
phase-2 ordering result. They use the same eight exact-dyadic ten-qubit
development targets, HF states, T=1, r=100 and Pauli-rotation budgets. They
are paired development controls, not new independent validation cases.

**The major gains do not require the finite-time sinc weights.** Replacing
them by one produces exactly the same selected unit order in six of eight
cases, not merely an equivalent error. The untransported J0 method still
beats both signed-JW and magnitude-JW in all eight cases, with the same
median error ratio 0.07181 against signed-JW as the original refined JD rule.
JD improves the other two cases: H2S_s1.00 by about 7.61%, and SiH4_s1.25
by about 30.63%, relative to J0.

**Unconstrained Pauli ordering can be more accurate than the conserving
partition.** Applying JD and the same insertion-search budget to individual
Pauli units beats refined JD in six cases and loses in two. It also beats
both JW sorting baselines in all eight cases. Its median error ratio against
signed-JW is 0.05627; against conserving refined JD it is 0.81330.

| Case | Refined J0 / fixed refined JD error | Individual-Pauli JD / fixed refined JD error | Individual-Pauli number leakage |
|---|---:|---:|---:|
| CH2O_s1.00 | 1 (identical order) | 0.27203 | 2.524e-11 |
| CH2O_s1.25 | 1 (identical order) | 0.25550 | 2.975e-11 |
| H2S_s1.00 | 1.08237 | 0.58272 | 5.231e-12 |
| H2S_s1.25 | 1 (identical order) | 0.92992 | 4.824e-11 |
| HCN_s1.00 | 1 (identical order) | 0.96392 | 7.144e-11 |
| HCN_s1.25 | 1 (identical order) | 1.01078 | 3.773e-10 |
| SiH4_s1.00 | 1 (identical order) | 0.69668 | 1.379e-11 |
| SiH4_s1.25 | 1.44147 | 1.13768 | 3.387e-12 |

All reported non-ties remain non-ties at absolute error tolerance 1e-12.
The conserving-JD error divided by individual-Pauli-JD error ranges from
0.8790 to 3.9139, with median 1.2554. Hence conservation should not be described
as improving accuracy for free: here it constrains the search and can cost
accuracy, while ensuring a property the unconstrained circuits lack.

Individual-Pauli number leakage ranges from 3.3874e-12 to 3.7728e-10. Refined
J0 leakage is 2.3187e-38–1.7841e-37, consistent with roundoff after exactly
number-conserving completed units. No leakage was projected out. This is a
tradeoff between accuracy, structural preservation and classical search cost,
not evidence that the fermionic grouping is always the most accurate choice.

Both arms use two fixed starts and at most ten accepted best-insertion moves.
That matches **search rules and quantum work**, not classical runtime. There
are 89–181 refined units versus 251–839 individual Pauli units. Recorded
preparation totals are about 2.43 s for refined J0 and 44.20 s for Pauli JD
across the eight cases, including objective setup and both starts. This is a
small-instance timing observation, not an asymptotic speedup theorem. The
original phase-2 driver additionally repeated five gauge controls, so its
total preparation time is not a directly matched runtime comparator.

## Interpretation and prior art

The evidence supports pursuing **state-specific ordering optimization subject
to particle-number conservation**, with the finite-time weighting an optional
refinement. It does not support attributing the large gains chiefly to the
transport filter. The combination still requires held-out testing, mechanism
analysis, and comparison with established HF/error-operator estimators and
ordering methods. The focused primary-literature note is
`../ordering-diagonal-refinement/NOVELTY_SCOPE.md` in the staging layout; it
should travel alongside these results when archived.

The new method cannot be advertised as the first HF-informed ordering,
first BCH-error optimization, first insertion search, or first finite-time
spectral error kernel. Those ingredients have prior literature. No additional
tuning was performed after these controls, and the original JD candidate
orders remained fixed.

## Artifacts and verification

- `PROTOCOL.md` and `run.py` were frozen before new errors. The phase-2 helper
  source is copied verbatim as `phase2_reference.py` and its hash is checked.
- All sixteen orders were saved before any new propagation or loading of
  comparator error rows. `results/all_candidates_frozen.json` records that
  checkpoint. The independent check also verifies file ordering.
- Active designs contain only their actual scheduled order. Parent blocks
  are labeled as auxiliary data for reproducing the original common closure;
  stale inherited class orders, graphs and hashes were not copied.
- Exact checks ensure every retained Pauli appears exactly once with its
  original coefficient. Repeated nonzero masks in the individual-Pauli arm
  are allowed. All same-destination-zero commutator amplitudes vanish exactly,
  validating the variance formula without a missing expectation subtraction.
- `check_ablation.py` independently checks exact target/reference identity,
  active schemas, unchanged budgets, candidate restrictions, source/artifact
  hashes and the objective via direct local Pauli actions. It imports no QHAT.
- `recheck_saved.py` independently propagates the circuits and dense references
  with NumPy/SciPy. Both eight-row replays passed under `python -O`, with
  maximum state differences 4.30e-16 (J0) and 4.67e-16 (Pauli JD).
- Producer wall time was 47.61 s, one compute thread, 6 GiB resident-memory
  cgroup cap, observed peak RSS 911312 KiB (about 890 MiB).

Reproduce the producer using pinned QHAT commit
`0be7d4857537668eec957dabaa8438ea5c0d33f2`, the pinned H2O Python environment,
and the two preceding saved-design archives:

```sh
systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
  /tmp/qhat-h2o-venv/bin/python run.py \
  --parent-run /path/to/followup-class-comparison \
  --phase2-run /path/to/ordering-diagonal-refinement/results \
  --qhat-root /path/to/pinned/qhat --output /tmp/new-objective-ablation
```

Then invoke `recheck_saved.py --run` separately on each arm directory, using
fresh output paths. `results/artifact_manifest.json` hashes the frozen producer
outputs; the independent audit reports are deliberately outside that manifest.
