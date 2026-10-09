# Controlled H2O reconstruction follow-up

This directory records the discussion's Q7/Q15 follow-up on slide 26. It uses
all 15 reconstructed active spaces from `numerical_audit/h2o_recovery/data`.
Original tensors are still missing; these results support/reassess the reported
case comparisons, but matching numerical errors does not certify identical
historical input coefficients or orbital phases. The frozen input directory
and previous recovery/pilot results are not modified.

Read `PROTOCOL.md` and `LOW_ENERGY_PROTOCOL.md`. The main sweep has 126 HF
approximate evolutions (90 primary plus36 convergence anchors), with 36 additional
low-energy/coherent controls in two eight-qubit models. Tensor policies are
sensitivity arms, not independent molecules. The scientific output is in
`presentation/REPORT.md`, `presentation/summary.json` and the accompanying CSVs
and PDFs. `results/` and `low_energy/` retain design/source/state evidence.

## Reproduce

Use the recovery environment pinned in
`../h2o_recovery/environment-requirements.txt`; no SCF rerun is needed. Replace
the paths below with the selected QHAT checkout and fresh output directories.
The checkout's parent must be first on PYTHONPATH. The main driver checks its
imported QHAT source paths against `--qhat-root`; it records source bytes and
hashes, input hashes, settings, all schedules and states.

```sh
PYTHONPATH=/path/containing/qhat /path/to/venv/bin/python run_h2o_sweep.py \
  --qhat-root /path/containing/qhat \
  --data /path/containing/qhat/numerical_audit/h2o_recovery/data \
  --pilot /path/containing/qhat/numerical_audit/h2o_recovery/validation \
  --output /fresh/h2o-sweep --max-qubits 14

# Resume verifies all completed design/state/source hashes before larger cases.
PYTHONPATH=/path/containing/qhat /path/to/venv/bin/python run_h2o_sweep.py \
  --qhat-root /path/containing/qhat \
  --data /path/containing/qhat/numerical_audit/h2o_recovery/data \
  --pilot /path/containing/qhat/numerical_audit/h2o_recovery/validation \
  --output /fresh/h2o-sweep --max-qubits 20

PYTHONPATH=/path/containing/qhat /path/to/venv/bin/python run_h2o_low_energy.py \
  --qhat-root /path/containing/qhat \
  --sweep /fresh/h2o-sweep --output /fresh/h2o-low-energy

/path/to/venv/bin/python recheck_h2o_states.py \
  --sweep /fresh/h2o-sweep --low-energy /fresh/h2o-low-energy \
  --data /path/containing/qhat/numerical_audit/h2o_recovery/data \
  --output /fresh/h2o-independent-state-check.json

/path/to/venv/bin/python summarize_h2o.py \
  --sweep /fresh/h2o-sweep --low-energy /fresh/h2o-low-energy \
  --historical-csv presentation/archived_h2o_rows.csv \
  --output /fresh/h2o-presentation
```

To reproduce the recorded historical-table SHA, use the full pinned QHAT
89c701652c4c64918825e279e84a7bea5924822a file
`analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv`.
The filtered15-row CSV included here regenerates the same H2O values/figures but
naturally has a different whole-file SHA. Historical values are never mixed
with fresh simulations in a statistical sample.

`guard_checks.json` records behavior tests for wrong-source imports, changed
resume protocols and a tampered state file under `python -O`. No checks depend
on Python `assert`. Every check uses an explicit exception.
A separate checker imports no QHAT code and remeasures error by the projector
residual identity rather than the production phase-alignment formula. This
independently checks saved states/metrics; it does not independently re-evolve
large states. Small dense checks do independently re-evolve: existing4/6-qubit
CAR/dense pilot evidence and all 36 new8-qubit low-energy state evolutions.

## Interrupted execution

The initial combined14–20-qubit execution received SIGTERM (tool exit143)
during the first20-qubit evolution. No scientific exception or Python memory
error was produced, and the precise external termination cause was not exposed.
The completed<=18-qubit results were retained and verified on resume. The
partial 20-qubit design/reference-check files were preserved separately in the
local staging directory; they are excluded from the portable evidence payload
and do not contribute rows to any table. The interrupted run log remains
`run_large.log`. The subsequent20-qubit-only session used the
same frozen scientific protocol and source.

The figure is suitable for an appendix: paired active-space reconstruction with
absolute errors and method ratios. It does not establish an accuracy advantage
across general molecules, an asymptotic qubit scaling law, or a convergence
claim to the full molecular Hamiltonian. No learning method is included.
