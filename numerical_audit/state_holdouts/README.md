# Correlated/coherent state and time holdouts

The [protocol](PROTOCOL.md) defines issue #4's fixed panel and the generalized
complex-state diagonal-response objective. This study distinguishes transfer
of frozen HF/T=1 orders from new state-aware placement and insertion searches.
It uses the three retained prospective molecular targets in a fixed orbital
basis. Exact diagonalization prepares the sector states; it is separately
charged and supplies no efficient state-preparation or hardware claim.

Use the existing Python 3.12 numerical interpreter. Only the pinned NumPy/SciPy
libraries (plus pytest for tests) are needed. No QHAT installation or optional
resource-estimation dependencies are needed; `response.py` loads the exact
neighboring determinant engine by file identity. The independent verifier uses
literal sparse Pauli matrices and unfused full-space circuits.

Run all local phases serially through `guarded.py`, with an unused label, the
absolute existing interpreter, and the protocol's finite per-phase allowances.
Never run a numerical phase outside the guard on this shared machine.

```
ABS_PYTHON -m pytest test_response.py -q
ABS_PYTHON study.py preflight --output ABS_FRESH_PREFLIGHT
ABS_PYTHON study.py prepare --output ABS_FRESH_PANEL
ABS_PYTHON study.py run --output ABS_PANEL
ABS_PYTHON verify.py --output ABS_PANEL
ABS_PYTHON study.py summary --output ABS_PANEL
```

The commands above are guard payloads. Example from the QHAT root:

```
python3 numerical_audit/state_holdouts/guarded.py --label i04-prepare --logs ABS_LOGS --cpu-seconds 1800 --wall-seconds 10800 --disk-mib 512 --output ABS_FRESH_PANEL -- ABS_PYTHON numerical_audit/state_holdouts/study.py prepare --output ABS_FRESH_PANEL
```

Preparation and propagation share one 1800 CPU second/10800 wall second panel
allowance; subtract preparation and earlier failed-attempt use before running
propagation. Independent verification has its own allowance. The wrapper's
initial memory-headroom wait does not start a worker. It freezes only its owned
service below the pause threshold and retains original logs/resource samples.
Each phase rejects existing outputs that it could overwrite; preserve partial
attempts. No destructive automatic resume is provided.

`prepare` authenticates the prior model/design hashes, requires committed active
source, fixes state vectors/phases and freezes every schedule before `run`.
`verify` independently checks the entire 432-row panel. `summary` authenticates
all verified artifacts and retains failed coverage strata, floors and losses.
The verifier binds state columns and named arm definitions independently of
producer metadata. It retains recomputed errors and checks floor/W/T/L
decisions explicitly; disagreement blocks a verified decision report and is
recorded as numerical ambiguity. The historical `spin_projection_leakage`
field measures the joint (N_alpha,N_beta) sector, not pure S_z leakage.
These are numerical checks, not Lean proofs or physical-error certificates.
