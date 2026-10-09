# H2O tensor reconstruction for slide 26

This recovers **reconstructed** inputs for the 15 archived H2O cases. The original
tensor files are unavailable. Matching geometry, basis and active-space labels
does not establish coefficient-identical historical inputs or reproduce old
Trotter errors. In particular, canonical orbital signs can differ between SCF
runs and change a signed-coefficient sorting rule. The orbital coefficients,
energies and occupations are retained so this can be investigated.

The driver imports the selected QHAT checkout, verifies that its three relevant
HF/active-space functions have the same AST as Albert's archived functions, and
uses `hamgen.my_run_pyscf` once per basis. It applies each active space using
`hamgen.apply_active_space`. It does not compute 20-qubit Pauli decompositions,
allocate large state vectors or run Trotter simulations.

## Run

Use an isolated Python 3.12 environment with `requirements.txt`. This needs the
Hamiltonian-generation dependencies, not QHAT's Julia/resource-estimation stack.
All Python packages and source hashes are recorded. `basis-set-exchange` supplies
the HGBS-5 definition locally; the precise basis resolved by PySCF is frozen in
each output directory.

```sh
python3.12 -m venv /tmp/qhat-h2o-venv
/tmp/qhat-h2o-venv/bin/python -m pip install -r numerical_audit/h2o_recovery/requirements.txt

MPLCONFIGDIR=/tmp/qhat-h2o-mpl /tmp/qhat-h2o-venv/bin/python \
  numerical_audit/h2o_recovery/recover_h2o.py \
  --qhat-root /home/bernalde/repos/qhat \
  --historical-source /tmp/qhat-albert-study \
  --output /path/to/h2o_recovery_outputs \
  --threads 2 --max-memory-mb 4096 --address-space-mb 12288
```

`--plan-only` checks basis sizes and input/source provenance without running SCF
or writing results. `--bases sto-6g` generates the six smaller cases first; a
subsequent invocation with `--bases hgbs-5` adds the remaining nine. Selection of
bases is not part of the common provenance fingerprint, so these invocations
can share one output directory. On the current environment STO-6G has 7 spatial
orbitals and HGBS-5 has 65; the full latter four-index integral array is about
136 MiB. A size guard, PySCF memory budget, process address-space limit and thread
limits bound resource use. Each basis requires one SCF calculation and the full
spatial integral transformation; all active spaces then reuse those integrals.

Inputs come from the archived checkout's
`hamiltonian_generator/polyatomic_library/H2O/s-1.00/{sto-6g,hgbs-5}` directory.
The driver validates literal values without executing those archived configs,
and freezes their original bytes, including the historical absolute paths.
Relocated configs are supplementary records. Use this driver for isolated
generation: the ordinary QHAT HF routine otherwise lets OpenFermion write its
HDF5 cache under its installed default data directory. This driver supplies an
explicit local filename and runs temporary HDF5 writes inside the output tree.

## Cases and physical conventions

The geometry, in Angstrom, is O=(0,0,0), H=(0.757540797781,0,0.587079636589), and
H=(-0.757540797781,0,0.587079636589). The system is neutral H2O, ten electrons,
singlet RHF, with PySCF symmetry disabled as in the archived generator.

| Basis | Active occupied/vacant spin-orbital counts |
| --- | --- |
| STO-6G | (2,2), (2,4), (4,4), (6,4), (8,4), (10,4) |
| HGBS-5 | (2,2), (2,4), (4,4), (4,6), (4,8), (6,8), (6,10), (8,10), (8,12) |

The two numbers are **spin-orbital** counts. For a pair `(o,v)`, there are `o`
active electrons, `o+v` qubits, `(10-o)/2` frozen doubly occupied spatial orbitals,
and `(o+v)/2` active spatial orbitals immediately above the frozen core. Alpha
and beta spin orbitals are interleaved. The initial determinant fills the first
`o` spin orbitals, with QHAT's big-endian state-vector index
`((1 << o) - 1) << v`. Each active Hamiltonian includes nuclear repulsion and the
frozen-core scalar and one-body contributions. Different active sizes describe
different Hamiltonians; the sweep is not a pure qubit-scaling comparison.

## Two explicitly separated coefficient policies

Every case produces two `.tensors.npz` files (**30 files for 15 cases**), with keys `constant`, `one_body`,
and `two_body` in the native QHAT/OpenFermion `InteractionOperator` convention:

`H = constant + sum(h[p,q] a†p aq) + sum(v[p,q,r,s] a†p a†q ar as)`.

- `historical_policy/` calls the archived-equivalent QHAT active-space generator.
  OpenFermion 1.8.1 zeros spatial-to-spin integral entries smaller than `1e-8`.
  This happens **before** the two-body factor of `1/2`, so its effective stored
  two-body coefficient threshold is `5e-9`. A later requested Pauli cutoff of
  `1e-12` cannot restore these discarded integrals.
- `unpruned/` explicitly expands the same active spatial integrals into spin
  orbitals, retaining all numerical entries and including the two-body `1/2`.
  This is not an assertion that every tiny entry is physically nonzero. It
  exposes the generator's truncation separately from downstream pruning.

The historical tensors are compared bitwise with an independent spin expansion
with the historical cutoff applied. The norm of the coefficient difference is
reported. This confirms the reconstruction policy, **not** equality with missing
historical tensors. Downstream comparison arms must all use the same selected
tensor policy and final pruning policy.

## Validation and outputs

Each SCF run must converge and give the expected closed-shell occupations.
Each case checks Hermiticity, orbital orthonormality, the exact historical cutoff
policy and an independent occupied-determinant contraction of the spin tensors.
The latter must reproduce the full HF energy after frozen-core reduction; the
historical-cutoff energy change is bounded by the removed coefficient l1 norm.

`manifest.json` contains environment/source/configuration provenance and SHA-256
hashes for every completed output. Each basis directory contains:

- both sets of tensors and `summary.json` with validation residuals;
- `orbitals.npz` with MO coefficients/energies/occupations and AO overlap/core;
- `basis_definition.json` and generation log;
- original and relocated configs and relevant current/historical source files.

The full `hf_molecule.hdf5` is a potentially large **optional local cache**. Its
hash is recorded separately in `optional_cache_sha256`; if present, resuming
verifies it, but the portable artifact need not contain it. Keep this cache out
of Git/Overleaf. Tensor NPZ files and orbital metadata are the required portable
outputs.

Rerunning with the same output first checks provenance and every completed file.
Changed generation source/environment/configuration or files cause an explicit failure, including under Python
`-O`. Partial calculations are kept in `.<basis>.partial` for inspection and are
not silently reused. Use a fresh output directory after an intentional code or
environment change. A HEAD-only QHAT commit change is allowed when all captured
source bytes are unchanged; the original run's commit stays in its provenance.
The top-level manifest is the integrity trust anchor; retain
its hash with the resulting research artifact. Output hashes certify the saved
files; they do not promise byte-identical HDF5 metadata across new SCF runs.

No old result tables are overwritten, and this driver alone does not validate
the historical H2O Trotter errors. Small independent evolution checks and later
controlled reruns are separate steps.

## Completed recovery and independent pilot

See [REPORT.md](REPORT.md) for the completed 15-case reconstruction and
72 small-space integration checks. `environment-requirements.txt` pins the
complete installed environment, including the numerical pilot dependencies.
Install that file instead of the generator-only `requirements.txt` to rerun
both stages. The environment itself is deliberately outside this repository.

The archived inputs came from QHAT commit
`89c701652c4c64918825e279e84a7bea5924822a` (the historical `L-sweep` study).
If the scratch checkout is unavailable, obtain that commit in a separate
checkout/worktree and pass its root as `--historical-source`. The original
configs and relevant archived generator sources are also preserved in
`data/*/configs` and `data/*/source`; they are evidence, not executable
instructions to write to the original absolute paths.

Run the bounded independent pilot from the QHAT checkout, with its parent on
the Python import path. The example output directory must not already exist:

```sh
PYTHONPATH=/home/bernalde/repos /tmp/qhat-h2o-venv/bin/python \
  numerical_audit/h2o_recovery/validate_h2o_reconstruction.py \
  --qhat-root /home/bernalde/repos/qhat \
  --manifest numerical_audit/h2o_recovery/data/manifest.json \
  --output /tmp/h2o-validation-fresh
```

The pilot validates both coefficient policies in every active space, then
runs only the 4- and 6-qubit cases. It checks the Pauli construction against
independent occupation-basis CAR matrices and compiled evolution against
dense Pauli rotations. It saves full states and makes no particle-sector
postselection. Use the recorded `validation/` as immutable evidence; use a
fresh output directory for a new run.

The separate checker verifies the saved evidence and recomputes overlap error
using an independent projector-residual identity, without evolving states:

```sh
/tmp/qhat-h2o-venv/bin/python \
  numerical_audit/h2o_recovery/recheck_saved_states.py \
  --manifest numerical_audit/h2o_recovery/data/manifest.json \
  --validation numerical_audit/h2o_recovery/validation \
  --output /tmp/h2o-saved-state-check.json
```

Do not interpret these 72 integration checks as 72 independent molecular
examples or as the completed 15-case performance replication. No new result
has been inserted into the paper's numerical figures or tables.
