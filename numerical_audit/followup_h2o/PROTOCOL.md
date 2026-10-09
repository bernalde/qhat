# H2O slide-26 reconstruction: frozen protocol, version 1

Frozen before examining any new performance data. Original missing tensors are
not recoverable from their configuration labels. The 2026-10-09 regenerated
inputs are used, with fixed MO signs and input hashes; results are reconstruction,
not coefficient-identical historical replay.

- All 15 recovered cases, separately for historical generator truncation and
  unpruned spin tensors (30 Hamiltonians), at T=1 and r=100.
- Predetermined convergence anchors: STO-6G (4,4), HGBS-5 (4,8) and (6,10),
  additionally r=50 and r=200. Total: 126 approximate evolutions.
- Three fixed schedules: signed-ascending canonical-parent first-owner Pauli
  order; JW signed ascending; JW magnitude descending. Lexicographic ties as
  the audited QHAT implementation. Ownership constructed once per input.
- Corrected explicit coefficient accumulation with final Pauli cutoff 1e-12
  for every method. No intermediate Pauli cutoff. Historical *tensor* policy
  does not restore known bugs in downstream Pauli accumulation or error metric.
- Same retained Hamiltonian, scalar discarded for all methods, identical HF
  determinant, time and r within each paired comparison. Each schedule uses
  all retained Pauli words exactly once per step, so equal r means equal Pauli
  exponential count for these three methods. No gate/depth advantage claimed.
- E=1-|<normalized exact, normalized approximate>| computed by globally phase
  aligned state residual. Also infidelity, unaligned state distance, norms,
  particle-number and spin-sector leakage. No postselection or renormalization
  restricted to the physical sector. Save complete states in the common
  Pauli-invariant affine closure together with its embedding indices.
- Approximate dynamics uses the exact affine closure of every retained Pauli
  flip. The exact reference is exp(-i H) applied to the same HF vector using
  sparse or matrix-free action of the common retained Pauli sum. Reference
  action must preserve any round-off residual rather than silently project it
  away. Verify against independent dense/CAR pilot states for <=6 qubits;
  verify large sparse reference by tighter/split-time propagation and norm.
- New comparison convention: absolute tie tolerance 1e-12; suppress ratio when
  comparator E<=1e-12, but retain raw errors. Additionally classify using the
  archived relative 1e-6/absolute 1e-14 and floor exclusion conventions, clearly
  labeled. Reconstruct archived CSV values separately, never merge with new
  simulations as additional independent cases.
- Display paired absolute errors and ratios by basis, labeling active electrons
  and qubits. Changing active space changes H, electron count and frozen core;
  these curves are not a qubit-only asymptotic law or convergence to full H2O.
- Save inputs, implementation, schedule, initial/final state and environment
  hashes. Independent saved-state checking recomputes metrics and leakage.
- Run sequentially with at most two numerical threads and a 12 GiB address-
  space cap, target peak resident memory <10 GiB. Start with small cases and
  record costs before larger cases. Stop on scientific validation failure.
- This protocol contains no learning campaign, selection of best random
  schedule, or retrospective relabeling of the best JW method as a deployable
  baseline. Additional structural-method comparisons are a separate experiment.
