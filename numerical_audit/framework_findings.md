# Local ordering and evolution audit, 2026-10-09

Source: QHAT `main` d3eea1ea458908f75d317794a80400a07896c2a8, incorporating
`improve_ordering` 5bc6ded98fe6ab52182c24f03139832b5f6349e6. Local branch:
`audit/fermionic-ordering-2026-10-09`. No commit or publication performed by this agent.

## Reproduced and repaired

1. `common/trotter_flattened.py`: collapsing a single-term repeat core multiplied
   its coefficient by r and then changed `num_steps` to one. Since the downstream
   time slice also becomes T rather than T/r, this evolves for rT. The repair
   leaves the coefficient unchanged when collapsing the repetitions. A Y-term
   witness at c=-0.61, T=0.23, hbar=1.7, r=3 showed max matrix-entry error
   0.1626306956 before this repair (after making circuit evaluation operational).

2. `common/pauli_string_evolution.py`: the Cirq decomposition previously used
   `PauliStringPhasorGate(P) ** (2 theta/pi)`, which fixes the positive eigenspace
   phase to one rather than exp(-i theta). For a Pauli with nonempty full support,
   this gives exp(i theta) exp(-i theta P). This phase matters for controlled
   evolution and absolute energy/spectral comparisons. Explicit exponents
   `exponent_neg=theta/pi`, `exponent_pos=-theta/pi` now preserve it. At Y,
   c=0.61,T=0.23,hbar=1.7, the old direct Cirq matrix differed in Frobenius norm by
   0.1166810932; its relative phase was 0.9965963806+0.0824357576i.

3. The same Cirq parity decomposition was given identity sites as active wires.
   In the installed Cirq 1.7.0 path, its XOR parity includes all supplied wires;
   XI can therefore rotate XZ rather than XI. The gate is now constructed only
   on nonidentity support and wired to those original positions. All-identity
   words use an explicit GlobalPhaseGate and leave physical wires idle.

4. Qualtran 0.7.0 uses `my_tensors`, not the old `add_my_tensors` hook. The
   implemented custom tensor protocol now uses the proper connection/bit-index
   pairs. CirqGateAsBloq uses a scalar QBit for one supported qubit, whereas the
   surrounding QHAT signature is always shaped; wiring now respects both cases.

These are framework correctness findings. We have not established that Albert's
historical paper figures used these paths. They must not be attributed to those
figures without tracing the actual historical entrypoint.

## Verified boundaries

The dense-Pauli reordering dispatcher preserves the input Hamiltonian and
coefficients. XYZ rejects mixed words instead of silently deleting them. A
flattened sequence of contiguous commuting blocks agrees with exact block Lie,
Strang, and fourth-order Suzuki products. Explicit chronological application is
first-to-last; the corresponding matrix factors appear right-to-left. Duplicate
Pauli labels at different nonadjacent positions retain both coefficients.

The older LCPSHamiltonian.commutes placeholder remains unconditional true. No
active benchmark/analysis caller was located outside its own examples; it was
not changed as part of this bounded production-path repair. The modern commuting
Pauli evolution validator uses the actual even-parity criterion.

## Validation

New file: `common/tests/test_ordering_matrix_oracles.py` (123 parameterized cases).
It uses independent SciPy exponentials on one to three qubits, checking total
physical time, identity phase, exact commuting evolution, noncommuting temporal
order, duplicate labels, block flattening, coefficient preservation, native
custom tensors versus circuit decomposition, and controlled-evolution phase.

Command, executed from QHAT root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/qhat-audit-matplotlib PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/home/bernalde/repos /tmp/qhat-audit-venv/bin/python -m pytest \
-p no:cacheprovider common/tests/test_ordering_matrix_oracles.py \
common/tests/test_pauli_string_evolution.py common/tests/test_trotter_flattened.py \
common/tests/test_commuting_pauli_string_evolution.py \
common/tests/test_alternation_reset.py analysis/tests/test_ordering.py -q
```

Result: 359 passed in 20.02 seconds; 59 upstream Qualtran static-cost deprecation
warnings, no failures. Output: `framework_test_results.txt`.
`git diff --check` passed for the changed evolution sources and new test file.

Focused environment: numpy 2.5.3, scipy 1.18.1, qualtran 0.7.0, cirq-core 1.7.0,
pytest 9.1.1. Installed pyLIQTR metadata is 0.3.3, below the repository's declared
>=1.4.2; the above tests do not constitute validation of a full declared QHAT
installation or its unrelated resource-estimation integrations.

No changes were needed to `analysis/ordering.py`. Production edits are limited
to the two evolution files described above. The exact gate decomposition has
changed, so any affected circuit resource counts should be regenerated; this
work verifies unitary semantics, not optimal gate synthesis.
