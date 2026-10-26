# Broader prospective molecular, state and time transfer

The frozen panel uses six authentic geometries in LiH, BeH2 and NH3, five
reference states, validation times 0.59, 1.41 and 2.73, eight arms and ideal-rotation
caps 10U and 30U. Each target has six spin modes, three active spatial orbitals,
two or four active electrons, and a nine-dimensional joint spin-number sector.
Full-space propagation retains leakage. Geometry screening covers the archived
prior validation records; it does not establish novelty across all chemistry.
The distinct LiH 1.841 preflight is excluded from the main results.

All 96 outer constructors, 1,332 unique complete circuits, 1,440 logical
observations, 90 reference evolutions and 5,040 arm comparisons are retained.
There are no censored main circuits. Both complete independent audits passed,
including the extracted archive replay with identical errors and decisions.
The 8,537-file payload, vectors, traces, costs, all losses/floors/ties and original
source snapshot are in `evidence/`. The original `delivery.tar.gz` has SHA256
`5873fbd34ffb369c7a4294724d213795f6ed358d3e5342a5fc2157cb9b6ac5e7`.

Producer and independent scientific source are commit
`a9992ef702b0d1d344778c5ea5ad08a048a20eef`; maintained documentation, CI and final
publication have separate pins. The practical prerequisite actually merged at
`da34971272a799b0dc6a6da6ccc028816feea138` and is integrated into this branch.

Training uses T=1.11, except the legacy current placement at T=1. Each frozen
schedule is evaluated at every validation time and cap. HF transfer and newly
state-aware first-order/cubic ordering are distinguished explicitly. The
first-trained symmetric arm is a separate formula control. The panel does not
isolate every initializer, objective and state effect as a complete factorial.
The LiH 1.713 HF certificate exposure at T=0.6 is recorded separately. The
transfer registry and exposure policy preceded those errors, and outcomes did
not select replacement inputs or methods.

Each accuracy comparison below contains 144 non-HF settings: six geometries,
four states, three times and two caps. Ties combine numerical-tolerance and
floor classes; the floor subset is shown. These dependent observations are
not 144 independent molecules.

| Method | Control | Wins/ties/losses | Floor ties |
|---|---|---|---:|
| State first-JD | HF first-JD | 138/6/0 | 0 |
| State cubic-JD symmetric | HF cubic-JD symmetric | 87/51/6 | 35 |
| State Pauli-JD symmetric | HF Pauli-JD symmetric | 84/31/29 | 15 |
| State cubic-JD symmetric | State first-trained symmetric | 77/66/1 | 54 |
| State cubic-JD symmetric | State Pauli-JD symmetric | 91/53/0 | 46 |
| State cubic-JD symmetric | HF Pauli-JD symmetric | 118/26/0 | 18 |

All controls and family/state/time/cap groups remain in `evidence/report/`.
Nominal search allowances differ; these are not matched-total CPU comparisons.
Reference preparation, input generation, initialization, search, returned-call
ordering/compilation charges and actual campaign work remain separate. Shared
logical arms are not multiplied into actual work. Non-HF compiled resources
cover evolution only; quantum preparation cost remains unknown. This finite
panel establishes no large-molecule, hardware, global-optimality or causal
conservation advantage. The floating audits are distinct from rigorous interval
certificates and Lean proofs.

`PROTOCOL.md`, `EXECUTION_PROTOCOL.md`, `REGISTRY.json` and the frozen evidence
record the predeclared design. `EXECUTION_CLOSEOUT.json` retains completed
resource attempts, including failures, through its stated cutoff. The main
producer used 550.14 CPU and 2,202.163 wall seconds. The two audits and packaging
used 455.29 CPU and 1,822.834 wall seconds, under the unchanged 512 MiB,
quarter-core profile. Paired reviews, publication/current-head CI, owning
handoff, separate eligible merges/actual CI and paper checks remain required.
