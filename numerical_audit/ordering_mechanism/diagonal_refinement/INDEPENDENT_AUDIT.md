# Read-only independent audit of the diagonal-response candidates

Checked 2026-10-09 by a separate agent. No frozen source, design, result, or
reviewed-repository file was modified. The scripts and JSON reports alongside
this note are new independent audit artifacts.

**The active comparisons are fair, and the implemented phase-2 objective and
unit refinement are mathematically consistent.** One inactive-metadata issue
needs explicit documentation before packaging these designs for reuse.

## Directly checked

`independent_schedule_audit.py` imports no QHAT or OpenFermion. Under `python -O`
it checked the eight complete-class designs and the eight refined designs
against `/tmp/followup-class-comparison`, including the recorded parent
manifest/design hashes and frozen candidate-source hashes.

- Retained integer coefficients, denominator, tensor hash, electron count,
  qubit count, affine basis, target hash and initial-state hash agree exactly
  with the paired baselines. Saved initial/exact/dense-reference arrays are
  identical to their parent arrays.
- Every active order is a permutation of its actual unit list. Every retained
  Pauli coefficient occurs exactly once, with the original value. Each
  reported rotation count equals retained-word count times the declared r.
  Time is T=1 throughout. No split angles or additional rotations are hidden.
- Each refined zero-mask unit is one diagonal Pauli word; every nonzero mask
  retains its complete original class. Internal words commute pairwise.
- I evaluated every unit on **every full ten-qubit computational basis state**
  with exact integer arithmetic. A nonzero transition always preserves the
  occupation count. This checks number conservation for every completed unit,
  not for its intermediate individual transverse Pauli rotations.
- Complete stage: 520 units and 532480 exact basis actions. Refined stage:
  952 units and 974848 actions. All passed.
- Refining the diagonal class leaves the fixed lexicographic circuit exactly
  unchanged: flattening its singleton diagonal units reproduces the parent's
  complete-class word order. Thus the inherited lexicographic error is a valid
  baseline for the refined fixed order.
- The refined `pair_setup` builds D by summing **all** zero-mask units. It does
  not use only the first such unit. Unit indices distinguish repeated masks;
  nonzero masks remain distinct. Zero/zero commutators vanish exactly.
- A separate direct-local-Pauli implementation recomputed all pair-commutator
  actions and J_D values for all 20 schedules in each stage. Maximum absolute
  discrepancy is 2.1685e-19 for complete classes and 1.7348e-18 for refined
  units. All accepted moves strictly decrease the recorded objective.

`independent_circuit_recheck.json` records a separate no-QHAT NumPy/SciPy
replay of all 20 phase-2 circuits and dense exact references, using the prior
`recheck_saved.py` under `python -O`. This is actual independent propagation,
not merely recalculation of metrics from saved states. The root agent's
`/tmp/ordering-diagonal-response/independent_recheck.json` records the analogous
complete-class replay; I inspected that report rather than rerunning it.

## Formula and selection inspection

For unit i, the code uses its exact action
`H_i |v> = f_i(v) |v xor x_i>` with a common coefficient denominator. The
pair amplitude is the action of `[H_i,H_j]` on HF:

`f_j(v) f_i(v xor x_j) - f_i(v) f_j(v xor x_i)`.

Summing the pair signs for an order gives the leading BCH action on HF.
Transporting this perturbation with the occupation-diagonal reference D
gives the squared sinc energy-gap weights in J_D. All nonzero destinations
are orthogonal to HF; the diagonal commutator expectation vanishes. Therefore
the squared response norm equals its variance in this setting. The code's
factor 1/4 is consistent with the leading Hermitian BCH kernel i B / 2.
It is a **diagonal-reference** variance, not an identity for full-H evolution.

The insertion search follows crossed adjacent transpositions and updates the
affected destination amplitude by minus twice that pair's signed amplitude.
Full recomputation checks each accepted composite update. Stable iteration
order supplies deterministic ties for diagonal units with the same mask;
the original unit list is sorted by full Pauli word there.

All candidates are written before the code loads old error rows or propagates
any new circuits. Observed file timestamps corroborate this sequence
(`candidate_freeze_timestamps.json`). This is not a cryptographic preregistration,
but source control flow excludes circuit-error feedback within either stage.
Phase 2 was conceived after seeing phase-1 performance, so both remain
algorithm-development data. Exact transported J_H is computed only after
candidate freeze. Smaller J_D is not proved to imply smaller actual error.

## Inactive inherited schema caveat

The refined design replaces `blocks.class` and `class_flip_masks` but inherits
old class-color and suffix orders, class graphs, suffix-bound metadata and
old schedule hashes. Those inactive fields do **not** describe the refined
units. For example CH2O_s1.00 has 89 refined units but its inherited
`class_color` and `class_suffix_l1` arrays contain 35 old indices. Interpreting
them against the new unit list would omit or misidentify coefficients.

This does not affect any `planned_runs`, candidate circuits, comparisons or
independent replays checked here. Preserve the frozen outputs; document the
fields as inherited parent-only metadata, or provide a separately generated,
cleaned design view with provenance. Future drivers should separate parent
metadata from active-unit metadata explicitly and record new schedule hashes.

The drivers hash their own source and protocol, but unlike the gauge driver
they do not themselves enforce an import-root/commit check for every QHAT
dependency. Exact independent replay substantially checks numerical behavior;
the final research artifact should additionally record the actual imported
dependency paths/hashes and command/environment provenance. This is a
packaging recommendation, not evidence of a wrong import in these runs.

## Limits of this audit

I did not regenerate molecular tensors, independently prove a bound comparing
J_D and J_H, or establish novelty. I did not infer which change caused any
gain: phase 2 simultaneously refines the diagonal partition and replaces
adjacent passes with insertion search. A factorial control is needed for that
attribution. The primary-literature assessment is in `NOVELTY_SCOPE.md`.
