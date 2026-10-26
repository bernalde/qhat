# Practical accuracy and total construction cost

The complete frozen [issue #16](https://github.com/bernalde/Fermionic_Trotterization/issues/16)
panel passed independent objective, primitive/unfused circuit, DAG resource,
physical-metric, target-selection, archive/extraction and byte checks. All 81
orders, 1,296 candidates, 1,281 complete circuits, 15 censored prefixes, 243 target
rows and 918 paired comparisons remain retained. All 243 targets were met,
including all 216 matched-cap target rows. Seventy complete construction overruns
are retained and excluded from matched-cap selection. Source review, Git delivery,
paper integration, publication and merge CI are still pending.

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
