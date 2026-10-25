# Controlled search costs and scaling

[PROTOCOL.md](PROTOCOL.md) declares a fixed eight-input molecular/synthetic panel, eight factorial optimizers, compatible charged warm starts, shared/method/search CPU accounting, symmetric warmups and repetitions, partial-sweep retention, sparse storage bounds and finite resource allowances. The retained study contains 432 constructed orders and 864 independently replayed full-space circuit comparisons. Read the [complete results](report/README.md): matched CPU caps do not show a uniform accuracy advantage from conserving grouping, and initializer effects remain separate from objective effects. All cases, losses, floors and incomplete sweeps remain visible.

Use the established pinned Python 3.12 numerical interpreter and the copied `guarded.py` wrapper. All numerical/test/build jobs on this shared machine run serially through that wrapper, with unused labels and absolute existing interpreter paths. Start/resume only above 5.5 GiB available; below 5 GiB freeze only the owned service. No local reviewer client overlaps numerical work. Subtract actual failed-attempt/phase use from each declared cumulative allowance.

Guard payloads, from this directory:

```
ABS_PYTHON -m pytest test_bounded.py -q
ABS_PYTHON preflight.py --output ABS_FRESH_PREFLIGHT
ABS_PYTHON study.py prepare --output ABS_FRESH_PANEL
ABS_PYTHON study.py run --output ABS_PANEL
ABS_PYTHON verify_costs.py --output ABS_PANEL
ABS_PYTHON reporting.py --results ABS_PANEL --output ABS_FRESH_REPORT
```

Preflight uses altered largest-shape coefficients, all factorial arms/caps and strong controls, both formulas, independent response/circuit/reference checks, compressed packaging and extraction. Its observations set conservative estimates before panel freeze. Preparation authenticates the unchanged prior inputs/helpers, commits scientific source identity and freezes every order before physical errors. Each data-producing phase preserves prior/partial attempts; use a new directory where required.

Reports retain every absolute error, related repetition, initializer/objective/partition pairing, nominal and actual CPU, incomplete sweep, overrun and strong-control comparison. No efficiency conclusion follows from equal examined positions or nominal caps when actual charges exceed them. Symmetric schedules reuse first-order objectives with exact pulse fusion; ideal rotations are separate from compiled or hardware resources. The two-level correction remains rejected and unchanged.

`deliver.py` copies verified results/report plus an explicit sanitized completed-phase execution ledger into a fresh delivery directory. `check_delivery.py` requires its full manifest. CI first authenticates and byte-regenerates the committed report, then replays physics in scratch and checks scientific decisions separately from platform-dependent diagnostic maxima. No account metadata, private model-review material or task control files belong in publication.


The corrected replacement panel measures optimizer calls through complete return;
the original internal-snapshot panel remains preserved at commit
8e03901e2063cfb8be15e360c762ac72d50caeca. Its omitted cleanup cost is not
retroactively treated as zero. `deliver.py --cumulative-ledger` retains the
sanitized continuing-campaign attempt ledger alongside raw process summaries.
It distinguishes completed phases from the earlier interrupted delivery checker,
which consumed61CPU seconds conservatively despite wrapper exit zero. Previous
and replacement attempts share the original finite phase ceilings. The ledger
covers the explicit labels completed before packaging; packaging/checking itself
is recorded separately in the owned task record, avoiding a self-hashing ledger.

Replacement outcomes can differ from the earlier internally timed panel because
CPU-limited search depends on measured execution. The observations do not isolate
the timing-boundary repair as the cause of every changed order or error. Every
replacement run is retained under the unchanged design; neither panel is selected
for favorable accuracy.


The delivery audit checks every warmup optimizer ledger, measured shared
initializer intervals and charges, and separately timed current-placement costs.
Full compact exports include warmup timing/cap fields and
[warmup shared charges](report/warmup_shared_setup.csv). The paper checker
reconciles campaign CPU totals and paired within-cap flags. Warmup shared-
constructor internal snapshots and common-setup internal snapshots were not
retained by the frozen producer; only their actual caller-side charges are
available. Missing diagnostics are not reconstructed or claimed as audited.
Current placement has its own timing schema and no bounded optimizer CPU cap.
These additional guards use fixed observations and producing source dea792a.
