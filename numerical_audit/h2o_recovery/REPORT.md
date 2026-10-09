# H2O recovery: reconstructed inputs and integration checks

Completed locally on 9 October 2026. This is input reconstruction and a
bounded pipeline validation, not a replication of all historical slide-26
Trotter results. The original tensors remain unavailable.

## Recovered inputs

All 15 surviving H2O configurations were reconstructed from their original
geometry and active-space definitions. One RHF calculation per basis supplies
all its active spaces. The generator's relevant HF/active-space functions were
checked against the archived source by AST equality.

| Basis | Cases | Active qubits | RHF energy (hartree) | Largest HF tensor-contraction residual |
| --- | ---: | --- | ---: | ---: |
| STO-6G | 6 | 4, 6, 8, 10, 12, 14 | -75.67884634773328 | 1.43e-14 |
| HGBS-5 | 9 | 4, 6, 8, 10, 12, 14, 16, 18, 20 | -76.02442786572156 | 1.28e-13 |

These residuals measure consistency of tensors with their originating SCF
calculation, not uncertainty in the physical energy. RHF used the historical
default convergence tolerance 1e-9 and converged in 6 and 8 cycles respectively.

Each configuration has two separately labeled exports: the historical
OpenFermion spatial-to-spin truncation policy, and an unpruned spin expansion
of the same active spatial integrals. There are **30 tensor NPZ files**.
The largest coefficient-l1 difference between those policies is 6.07e-10.
This is distinct from the subsequently applied final Pauli cutoff of 1e-12.

The portable input payload (about 719 kB, excluding its top-level manifest)
contains all tensors, orbitals, exact resolved basis definitions, original and
relocated configurations, generation source snapshots, validation summaries
and logs. Optional full-HF HDF5 caches stay local and outside Git. The complete
61-package generation/validation environment is pinned in
`environment-requirements.txt`; the earlier audit environment was not changed.

Input manifest SHA-256:
`43e3ff2123b37978457a8a2473efbc3256fed8c1cb3ec8f0b00c31a785550284`.

## Independent validation

`validate_h2o_reconstruction.py` checked all 30 tensor files. It additionally
ran both coefficient policies at 4 and 6 qubits in both bases: four physical
configurations, eight configuration/policy combinations. Each used the same
Hamiltonian, HF state and T=1 across fermionic signed, JW signed and JW
magnitude schedules, at r=50, 100 and 200, for **72 evolutions**.

- Explicit occupation-basis CAR matrices agree with the retained Pauli sum
  plus scalar to at most 5.94e-13 in operator norm, within the final-cutoff
  allowance and numerical tolerance.
- Independent dense Pauli products agree with compiled matrix-free evolution
  to at most 8.66e-14 in state norm. No fixed-number projection is used.
- Maximum approximate-state norm drift is 9.09e-14.
- The dense owned-block number commutators vanish in these small cases;
  maximum measured owned-schedule number leakage is 2.35e-37.
- Initial, exact and approximate states, actual schedules, source identities,
  tensor hashes and stable full-state metrics are saved in `validation/`.

The small-space checks establish that the reconstructed inputs can feed the
audited numerical path. They do not establish an accuracy advantage, recover
the historical orbital gauge, or replace the full 15-case controlled sweep.
Repeated r values and coefficient policies are not additional independent
molecular examples. Canonical orbital signs/degenerate rotations can change
signed-coefficient sorting; saved MO coefficients are therefore part of the
input record.

## Integrity and scope

A portable copy without optional HF caches verifies successfully under
`python -O`. Appending bytes to one copied tensor is explicitly rejected under
`-O`; the original data were not modified. See `integrity_checks.json`.

The saved-state checker independently verifies all 168 saved arrays and all
72 distinct settings. Its alternative overlap formula agrees within
1.14e-21 in absolute E. Relocated evidence also verifies without reading
recorded original absolute tensor paths; invalid state paths are rejected.
See `independent_saved_state_review.json` and
`rechecker_portability_checks.json`.

The original manuscript numerical tables, reviewed campaigns, Lean files and
frozen review-5 archive were not changed. No remote branch was pushed and no
Overleaf synchronization was performed. A full reconstructed H2O comparison
can now be run with a declared choice of tensor policy and both fixed JW
comparators; it must retain the distinction from historical replication.
