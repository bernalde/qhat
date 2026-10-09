# Ordering-mechanism research checkpoint

Local experiments from 2026-10-09, isolated from the manuscript/QHAT revision
under external review. The paper worktree's
`theory/ordering_mechanism/RESEARCH_REPORT.md` explains the conclusions and
limitations; `ORDERING_RESEARCH_BRIEF.md` is the continuation handoff.

The frozen conserving candidate improves on four fixed-order baselines in
eight development and six prospective cases. It is not the most accurate
method in every comparison: optimizing individual Paulis wins against it in
six of eight development cases and leaks particle number. The finite-time
gap filter is not the main source of the gains. Read the ablations before
making a novelty or superiority claim.

## Contents

| Directory | Scope |
|---|---|
| `gauge/` | 132 old-rule orbital-sign/control evolutions; independent saved checks plus 18 small-space repropagations |
| `state_controls/` | 18 old-rule 14/16-qubit ground/triplet/coherent evolutions; reference and saved-state checks |
| `diagonal_response/` | First complete-class candidate, 20 development evolutions and exact response diagnostics |
| `diagonal_refinement/` | Frozen split-diagonal/insertion candidate, 20 evolutions; independent fairness, propagation and literature checks |
| `factorial/` | 16 missing arms of the partition/search 2x2 control, independently replayed |
| `objective_ablation/` | 16 unweighted-response/individual-Pauli controls, independently replayed |
| `heldout/` | Visible report/check summaries and complete 95-member ZIP of the prospective validation artifacts |
| `work_curves/` | 240 new fixed-order small-r evolutions, 40 reused r100 measurements, threshold records and independent provenance/replay checks |
| `theory_prototype/` | Written theory and independent matrix/prototype audits; numerical checks are not Lean proofs |
| `packaging/active_*_designs/` | Clean derived views of active candidate designs; frozen originals are preserved |

The standalone research Lean module and readable proofs are in the paper
worktree. Its 19 lemmas verify the matrix/amplitude/sign bridge, not the whole
new algorithm or analytic approximation guarantee.

## Reproduction and integrity

The original reviewed QHAT dependency commit is
`0be7d4857537668eec957dabaa8438ea5c0d33f2`. Producers use its pinned H2O Python
environment (`environment-requirements.txt`). Scripts state their input paths,
and some enforce this exact checkout. Merely launching a script from a new
directory does not select its QHAT import: explicitly control `PYTHONPATH`
and the runner's `--qhat-root` argument when supplied.

The portable integrity check, including all 94 embedded held-out files, is:

```sh
python -O verify_artifacts.py
```

This only checks bytes. It does not replace numerical replay or Lean checking.
All failure gates use explicit exceptions; optimization mode does not disable
them. The manifest does not authenticate itself: retain its separately recorded
SHA-256 in `artifact_manifest.sha256` or the paper checkpoint.

For independent NumPy/SciPy circuit reconstruction of either candidate stage,
factorial control, or work curves, use a fresh output file:

```sh
/path/to/python -O recheck_candidate_runs.py \
  --run diagonal_refinement/results --output /tmp/new-replay.json
```

This imports no QHAT or OpenFermion. It rebuilds exact references and actual
circuits from the saved design and arrays. Objective ablations have the same
checker inside their directory; invoke it separately for each arm. Gauge and
state controls have their own documented saved-state checkers.

The prospective ZIP is compressed because its exact-integer design JSON totals
about 104 MB expanded. It contains source, tensors, schedules, arrays, reports
and a per-file manifest. Extract to a fresh directory, then follow REPORT.md:

```sh
python -m zipfile -e heldout/heldout_validation.zip /tmp/heldout-replay
/path/to/python -O /tmp/heldout-replay/recheck_validation.py \
  --results /tmp/heldout-replay/results --output /tmp/new-heldout-check.json
```

`independent_circuit_replay.py` independently propagates seven specified
prospective runs with a separately derived in-place Pauli kernel. The other
23 prospective runs have target/reference/saved-metric checks, not independent
full propagation. No claim of reproducing unavailable original H2O tensors
is made by these new experiments.

## Frozen-schema and provenance notes

The two first candidate drivers inherited metadata from their parent designs.
After diagonal splitting, unused old class-color/suffix orders and graphs in
the refined raw design do not describe the new unit list. The active schedules
and all reported experiments were independently validated. Do not execute
inactive inherited orders. `packaging/extract_active_designs.py` generates
separate clean views with current schedule hashes, emitting every retained
Pauli exactly once. It never changes the frozen experiment.

The later factorial, ablation and work-curve designs retain only scheduled
orders; the prospective driver has its own explicit schema. The work-curve
runner's reused r100 compatibility was checked independently against targets,
states, schedules and stored results. Future generic reuse should enforce
these guards inside the runner as well.

`dependency_inventory.json` is a **post-run** source inventory from the
unchanged reviewed checkout, not a retroactive run-time import attestation.
Independent propagation checks validate the recorded circuits directly.
The mathematical targets are paired exactly, but this does not establish
physical integral-generation accuracy or broad molecular generalization.

Work comparisons price Pauli-word rotations. They do not include compilation,
two-qubit-gate cancellation, depth or hardware errors. Coefficient-scaled
prospective cases are not independent geometries, and all new candidate
performance uses HF initial states. The reviewed manuscript remains unchanged.
