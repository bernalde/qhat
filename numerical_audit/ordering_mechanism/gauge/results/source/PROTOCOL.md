# Fixed orbital-sign sensitivity experiment

Frozen before performance evaluation, 2026-10-09. This is a representation
robustness diagnostic, not a best-gauge ordering algorithm or a random-sample
estimate of a gauge distribution.

Historical-policy reconstructed H2O inputs: STO-6G (6 occupied, 4 virtual)
spin orbitals (10 qubits), and HGBS-5 (4,8), (6,8), (6,10) (12,14,16 qubits).
The first is a control; the HGBS cases include one previous failure and two
previous successes. All use the HF determinant, physical time T=1, r=100,
first-order chronological Pauli rotations, and the prior final-Pauli cutoff
1e-12. No error-guided choice or optimization is performed.

For m spatial orbitals use these six signs in this fixed order:
1. All positive (baseline).
2. First spatial orbital negative, the rest positive.
3. Last spatial orbital negative, the rest positive.
4. Odd-index spatial orbitals negative (0-based).
5. First floor(m/2) spatial orbitals negative.
6. Indices congruent to 1 modulo 3 negative.
Alpha and beta signs are identical. These are distinct nontrivial spatial
patterns for every case, though a Hamiltonian may have symmetries making
some patterns equivalent on its reachable subspace. Such cases are retained.

Transform one- and two-body tensors by the corresponding products of signs,
rebuild signed-first-owned, signed-JW and magnitude-JW schedules from scratch,
and propagate all three. Check the transformed retained coefficients against
exact sign multiplication of the original retained coefficients; no tolerance
is accepted for this equality. Check common Pauli keys, scalar, flip masks,
affine closure and exact sparse H'=DHD. If cutoff or floating accumulation
breaks equality, stop and report it rather than asserting gauge equivalence.

The unchanged magnitude-JW schedule is a positive control. In addition evolve
each transformed coefficient list in the corresponding original fixed order,
and compare against D times the original final state. This tests conjugation
without conflating it with re-sorting. Baseline fixed and re-sorted trajectories
are identical by construction and are evaluated once. Total: 72 re-sorted
trajectories plus 60 explicit fixed-order controls = 132 propagations.

All Pauli rotations propagate on the same invariant affine closure as the
original experiment, retaining number/spin leakage. Exact references are the
prior saved reference transformed by D, after verifying original designs,
state hashes, and Hamiltonian coefficients. Initial states are transformed
by D, including any global phase. Record full final arrays, leakage, norms,
stable E=1-|<normalized exact, normalized approximate>| computed via half the
squared phase-aligned distance, schedule positions and owned-block changes.
No postselection. Tie tolerance/analysis floor is absolute 1e-12.

Independently check full-Hilbert-space Pauli action, dense eigendecomposition
reference, and NumPy Pauli-product propagation for the 10-qubit case using
direct per-site actions, not QHAT compiled masks. Include deliberate wrong-sign
negative controls that the exact covariance checks must reject.

Use at most two compute threads and an 8 GiB resident-memory cgroup cap.
Store source, input, protocol, and output hashes. Do not modify reviewed repos.
