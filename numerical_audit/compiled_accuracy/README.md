# Compiled target-accuracy study preparation

The [protocol](PROTOCOL.md) freezes compiler/version, native gate set,
nearest-neighbor line, routing/cancellation, eight strong comparator arms,
three molecular inputs, both times, all eight step counts, both accuracy targets
and cumulative resource ceilings before execution. This complete source pipeline
is prepared but **UNEXECUTED**. No compiled-resource or accuracy result is claimed.

`compiler.py` retains every primitive and exact rational angle, including HF
preparation, restored routing and censored prefixes. `study.py` authenticates the
fixed 30U-trained issue8 orders, compiles every grid candidate before any new
physical errors, then propagates every feasible full-space circuit. Exact unit
fusion requires internal Pauli commutation, checked with integer Pauli parity.
It never projects a state into a symmetry sector or retunes an order on compiled
cost/error. Selection minimizes qualifying CX, then unit-duration ASAP depth,
total gates and steps, with all unmet, ambiguous and censored points retained.

`reference.py` independently uses literal sparse primitive/Pauli matrices,
unfused product formulas, scaled-Taylor target propagation and a two-pass circuit
DAG resource/depth audit. `verify.py` checks all states, resources, target decisions
and provenance. Deterministic canonical trace reproduction reuses the producer
compiler and is identified separately from the independent replay/DAG oracle.
The retained gzip bytes are integrity-pinned; canonical exact gate/angle records
are reproduced without depending on the compression-library encoding.

`reporting.py` retains the entire grid, every target decision and every comparator
pair, with CX and depth outcomes labelled separately. It separates inherited
standalone ordering charges, one-circuit construction phases and actual full-grid
compilation/physical-screening costs. Repeated standalone costs and shared
reference timings are never summed as fresh campaign work. `delivery.py`
authenticates the exact artifact inventory and regenerates reports.

`test_compiler.py` supplies full small-register literal-unitary checks, routing,
HF phase, exact rewrites and trace/resource fixtures. `test_reference.py` adds
closed-form target/primitive states, independent parallel depth, joint-sector
labels, nonmonotone target selection, threshold ambiguity, censoring and invalid
unit-fusion checks. `test_compiled_guard.py` covers budget stops. All tests remain
unexecuted at this preparation revision.

The representative preflight uses the largest altered H2O target (coefficients
scaled by 7/5), T=2.11, all eight fixed-order arms and all eight step counts:64
candidates and16target decisions. It runs compilation, full feasible propagation,
independent replay/DAG/selection, reporting, packaging, archive/extraction and
strict checking before calibrating and freezing main-panel estimates. A passing
pilot cannot substitute for the complete48curve/384candidate main study.

Use the verified Python3.12 numerical interpreter and pinned requirements.
Every local numerical/test phase is a payload of `guarded.py`, with an unused
label, original cumulative allowance and fresh owned output directory. The guard
uses one numerical thread, quarter-core CPU, low priorities,512MiB/no-swap and
memory headroom rules. A watchdog stop fails even when the supervisor returns
zero; require both termination metadata and the completed result artifact.
Do not run these payloads directly on the shared machine:

```
ABS_PYTHON -m pytest test_compiler.py test_compiled_guard.py test_reference.py -q
ABS_PYTHON preflight.py --output ABS_FRESH_PREFLIGHT --ordering-evidence ABS_I08_ALTERED_PREFLIGHT_DELIVERY --ordering-commit dcbb8667161946813854152a41453b401ad3e36b
ABS_PYTHON study.py prepare --output ABS_RESULTS --ordering-evidence ABS_MERGED_I08_DELIVERY --ordering-commit 8fd424f532cce04b595abc3445cc032ab03ab6c4
ABS_PYTHON study.py run --output ABS_RESULTS
ABS_PYTHON verify.py --output ABS_RESULTS
ABS_PYTHON reporting.py --output ABS_RESULTS --report ABS_FRESH_REPORT
ABS_PYTHON delivery.py pack --results ABS_RESULTS --report ABS_REPORT --delivery ABS_FRESH_DELIVERY
ABS_PYTHON delivery.py check --delivery ABS_DELIVERY
```

Keep all failed/interrupted attempts and cumulative CPU/wall/disk consumption.
Fresh output paths, source commits, workflow resumption or review packets do not
reset any allowance. Main execution, independent verification, evidence package,
paper companion and required independent reviews remain pending.
