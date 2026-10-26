# Practical accuracy and total construction cost

The complete frozen [issue #16](https://github.com/bernalde/Fermionic_Trotterization/issues/16)
panel passed independent objective, primitive/unfused circuit, DAG resource,
physical-metric, target-selection, archive/extraction and byte checks. All 81
orders, 1,296 candidates, 1,281 complete circuits, 15 censored prefixes, 243 target
rows and 918 paired comparisons remain retained. All 243 targets were met,
including all 216 matched-cap target rows. Seventy complete construction overruns
are retained and excluded from matched-cap selection. All 5,542 delivered files, including binary states and compressed traces, pass
staged Git blob checks. The original deterministic archive is retained unchanged.
Source review, paper integration, publication and merge CI are still pending.

At each target/cap, counts below compare cubic-trained symmetric ordering with
first-order-trained symmetric ordering. Each cell has nine paired observations:
three fixed molecular models times three construction repetitions. Values are
wins/ties/losses for the selected circuit under the declared CX-first rule.

| State-error target | Total construction cap (CPU s) | CX | Selected-circuit ASAP depth |
|---:|---:|---:|---:|
| 1e-4 | 2 | 6/0/3 | 6/0/3 |
| 1e-4 | 8 | 3/3/3 | 6/0/3 |
| 1e-6 | 2 | 6/0/3 | 6/0/3 |
| 1e-6 | 8 | 6/0/3 | 6/0/3 |
| 1e-8 | 2 | 6/0/3 | 6/0/3 |
| 1e-8 | 8 | 9/0/0 | 9/0/0 |

The tighter target and finer grid expose an advantage in a specific regime;
losses remain at the smaller cap and looser targets. Charged warm starts give
mixed additional changes. The complete tables also retain every matched-Pauli
and inherited strong-reference comparison. Inherited references are explicitly
outside matched CPU caps. These dependent observations are descriptive counts,
not nine independent molecular instances.

[PROTOCOL.md](PROTOCOL.md) defines the unchanged models, HF state, T=2.11,
step grid, all controls, exact supplied coefficients, compiler and accounting
boundaries. Compilation includes serialization/hash/cleanup; offline full-grid
physical screening and reference calculations are separate campaign costs.
Selected depth follows the CX-first selection rule and is not separately depth
optimal. This fixed-basis development panel establishes no new-family transfer,
hardware execution gain, global optimum or conservation-only mechanism.

The main measurement used 1,418.86 CPU seconds and 5,676.388 wall seconds.
The successful independent replay/delivery used 1,312.87 CPU seconds and
5,253.929 wall seconds. The preceding syntax-only wrapper failure remains
retained and conservatively charged separately; no measured circuit was rerun.
The original 3,600-CPU-second ceilings for those phases were unchanged.

The complete payload is under [evidence/](evidence/); its report includes every
candidate, target and pair. `delivery.tar.gz` has the identical inner payload
under its original `delivery/` root. `delivery_verification.json`,
`source_provenance.json`, `preflight_execution.json`, `execution_resources.json`
and `packaging_resources.json` record source identity, preflight feasibility,
the approved preflight extension and failed/successful attempt accounting.
The resource ledger states its cutoff explicitly; later work remains in the
campaign ledger. Original reviewer reports remain private workflow evidence.

With the pinned environment, check the exact committed delivery:

```sh
python numerical_audit/practical_accuracy/delivery.py check \
  --delivery numerical_audit/practical_accuracy/evidence --git-commit HEAD
```

The CI workflow additionally reruns the complete independent numerical checker
in a fresh output copy. Paper evidence depends on the actual numerical merge.

Independent replay requires the exact recorded Python release and NumPy/SciPy
versions. The producer's and checker's executable paths and full interpreter
build strings are separate provenance; they may differ across pinned machines.
Producer continuation still requires its full original environment identity.
The relocation regression uses the actual committed payload and rejects Python
release or NumPy version changes. Scientific data, tolerances and selection rules
are unchanged by this portability repair.

The original frozen `evidence/report/SUMMARY.json` groups only the 270 matched
cubic comparisons. All 648 comparisons against inherited strong controls remain
in `pairs.csv` and `COMPARISONS.json`. The additive
`all_comparison_groups.json` groups **all 918 comparisons**, with each control's
cost cohort explicit. Its CI check is mandatory. The original delivery manifest,
raw rows, binary states, traces and deterministic archive retain their bytes.

The repaired decision checker binds copied row identities, construction charges
and resources to authenticated candidates and independently audited DAG costs.
It derives paired outcomes from independent selections and costs, without using
the producer comparator. Focused tests reproduce and reject corrupted copied
rows and a comparison bug shared with the producer. The original full numerical
replay remains a separately identified execution record; post-review checks and
current-head full CI replay remain distinct.
