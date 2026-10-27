# QHAT Standalone: Self-Contained Trotterization Tools

This directory contains a self-contained version of the quantum simulation code with no external dependencies (outside of Julia packages).

## Contents

- **`src/`** - Core simulation modules (copied from parent repo)
- **`Li-Li_jw/`** - Test Hamiltonians for Li-Li molecule (4-16 qubits, Jordan-Wigner encoding)
- **`He-He/`** - Test Hamiltonians for He-He molecule
- **`Li-Li_bond/`** - Li-Li Hamiltonians across multiple bond lengths, bases, active spaces, and encodings
- **`trotter_expts.jl`** - Trotter error demonstration script
- **`trotter_ground_energy.jl`** - Read metadata reference energy and compute Trotter-derived energy estimates for one Hamiltonian
- **`bond_energy_curve.jl`** - Exact plus second-order Trotter energy-vs-bond-length curves
- **`Project.toml`** - Julia package dependencies

## Quick Start
Before the first run. Start julia REPL, hit `]` to get into package mode and run `activate .` and `instantiate` to have all the dependencies installed. Then outside the REP do the following: 

```bash
cd julia_trotter
julia --project=. trotter_expts.jl
```

This runs the demonstration on the smallest system (4 qubits). To use a different Hamiltonian:

```bash
julia --project=. trotter_expts.jl Li-Li_jw/Li-Li_2.90_hgbs-5_as-004-004_jw.dat  # 8 qubits
```

## Additional Scripts

### `trotter_ground_energy.jl`

Reads a single Hamiltonian, builds second-order Trotter product unitaries, and reports:
- metadata ground-state energy
- Trotter-derived ground-state energies for `nsteps = 1, 5, 10`

The default Hamiltonian is `He-He/He-He_2.40_hgbs-5_as-004-004_jw.dat`.
By default it uses Arpack on `U + U'`. Pass `--arnoldi` to use matrix-free Arnoldi with the Hartree-Fock start vector:

```bash
julia --project=. trotter_ground_energy.jl
julia --project=. trotter_ground_energy.jl --arnoldi
julia --project=. trotter_ground_energy.jl --arnoldi He-He/He-He_2.40_hgbs-5_as-004-004_jw.dat
```

### `bond_energy_curve.jl`

Builds a Li-Li potential energy curve across all bond lengths in `Li-Li_bond/` for a fixed:
- basis
- active occupied orbitals
- active vacant orbitals
- encoding

It prints rows as each bond length finishes:

```text
bond_length,exact_energy_hartree,trotter_1_step_hartree,trotter_5_step_hartree,exact_seconds,trotter_1_seconds,trotter_5_seconds,total_seconds,file
```

and saves a plot overlaying:
- exact ground-state energy
- second-order Trotter, 1 step
- second-order Trotter, 5 steps

Examples:

```bash
julia --project=. bond_energy_curve.jl
julia --project=. bond_energy_curve.jl --arnoldi
julia --project=. bond_energy_curve.jl --basis hgbs-5 --electrons 4 --vacant 8 --encoding jw
julia --project=. bond_energy_curve.jl --basis sto-6g --electrons 4 --vacant 8 --encoding jw
```

The default curve is `hgbs-5`, `as-004-004`, `jw`.

## Bond-Length Data Layout

The bond scan data is organized as:

```text
Li-Li_bond/<bond_length>/<basis>/Li-Li_<bond>_<basis>_as-EEE-VVV_<encoding>.dat
```

where:
- `EEE` = active occupied single-occupancy orbitals
- `VVV` = active vacant single-occupancy orbitals
- `encoding` = `jw` or `bk`

## What the Demonstration Shows

The `trotter_expts.jl` script provides a comprehensive analysis of Trotterization errors:

### 1. System Information
- Number of qubits and Hamiltonian terms
- Metadata energy as supplied, with sector provenance unverified
- L1 norm and normalization factor

### 2. State Vector Errors
Compares Trotter evolution against Chebyshev reference:
- **1st order Trotter**: Basic splitting, O(dt) global error at fixed total time
- **2nd order Trotter**: Strang splitting, O(dt²) global error at fixed total time
- **4th order Trotter**: Recursive composition, O(dt⁴) global error at fixed total time

The corresponding local, one-step errors are O(dt²), O(dt³), and O(dt⁵).

### 3. Operator Errors
Compares full unitary matrices (spectral norm):
- `‖U_trotter - U_reference‖` estimates for varying step counts
- Reports an observed ratio against the expected O(dt²) global scaling for
  2nd-order methods at fixed total time; a finite reference can mask that scaling

### 4. Commutator Error Bounds
Theoretical upper bounds from Childs et al. (PRX 2021):
- Exact nested commutator computation
- Comparison of bounds against estimated errors relative to the approximate reference
- Bound/estimate ratios, subject to reference accuracy

### 5. Error Scaling Analysis
- Reports the observed error ratio when doubling steps alongside the expected
  asymptotic ratio of four for global second-order error at fixed total time
- Reports the average bound/estimate ratio; no target-accuracy recommendations
  are computed

## Test Hamiltonians

Seven Li-Li molecule Hamiltonians with varying sizes generate by QHAT

The `Li-Li_bond/` directory contains bond-length scans for:
- bases `hgbs-5` and `sto-6g`
- active spaces `002-002`, `002-004`, `004-004`, `004-006`, `004-008`, `006-008`, `006-010`
- encodings `jw` and `bk`

## Implementation Details

- **Chebyshev order**: 10 (adjustable via `CHEBY_ORDER` constant)
- **Evolution time**: 1.0 (total time, split across Trotter steps)
- **Test step counts**: [10, 20, 50, 100]
- **Total runtime**: ~10-30 seconds for 4-qubit system

For the curve scripts:
- exact energies come from metadata if present, otherwise from sparse eigensolves of the Hamiltonian
- Trotter energies are extracted from the largest eigenvalues of `U + U'`
- the default Trotter eigensolver is Arpack on the full Trotter unitary
- passing `--arnoldi` switches to matrix-free Arnoldi actions for `U + U'`, which avoids constructing the full unitary


## Reference accuracy and energy scope

The demonstration uses an approximate Chebyshev reference, not an exact
exponential or an outward-rounded certificate. The backend subdivides the
propagation and uses Bessel arguments 1/3 (state) and 1/4 (operator) at the
default settings. Order 10 is therefore not an unscaled expansion at argument
π. The script reports the existing finite Bessel-tail estimates; these exclude
arithmetic roundoff. Errors near floating-point resolution need independent
reference validation and must not be interpreted as resolved Trotter errors.
The matrix backend retains its default matrix-return API; `return_info=true`
returns the matrix and the existing `(error, dt)` estimate metadata.

Operator errors use the largest singular value of the unitary difference,
including above the dense-size threshold. The iterative estimate checks
convergence and singular-vector residuals; it is not a rigorous norm bound.

Both energy CLIs warn that their operators act on full Fock space, without a
fixed electron-number or spin-sector restriction. Metadata energies are used
as supplied; their sector provenance is not checked. The Hamiltonian fallback
and default Arpack calculation target full-space extrema. These quantities
must not be labeled fixed-electron molecular ground energies without a separate
sector justification.

Saved artifacts inherit this scope caveat: the CSV and PNG are not standalone
claims about fixed-electron molecular ground energies. Their legacy `exact`
labels refer to a reference energy taken from a supplied `smallest eigenvalue`
metadata key when present, otherwise from a full-space Hamiltonian eigensolve;
the current row format does not identify that source separately. Retain the
source Hamiltonian metadata, these instructions and the stderr scope warnings
alongside exported tables/plots. The demonstration's printed ground-state energy
is likewise an unverified metadata reference, not an independently established
sector ground energy.

The optional Arnoldi method uses a Hartree-Fock start vector, which is only a
Krylov seed. It is not a sector projection or a guarantee of the global minimum;
invariant subspaces and overlap can affect the extrema found. Individual Pauli
rotations need not conserve particle number. Projecting a leaking unitary to a
sector would generally destroy unitarity and does not repair the existing
phase-to-energy interpretation. The implemented HF seed is JW-specific:
BK plus `--arnoldi` is rejected, and HF-seeded calls require explicit JW mapping
metadata. The full-space Arpack route remains available for BK inputs.

## Regression checks

Use the existing Julia project and installed dependencies; no installation is
performed by these commands:

```bash
julia --project=julia_trotter --startup-file=no julia_trotter/test/runtests.jl
julia --project=julia_trotter --startup-file=no julia_trotter/test/reference_and_cli_tests.jl
```

The first command covers small singular-value and scope/encoding fixtures.
The second loads the existing simulation/CLI dependencies and checks the
reference against a small dense exponential, backward-compatible return values,
actual CLI/backend rejection of unsupported HF encoding, and acceptance of a
committed JW header. A nondegenerate synthetic two-qubit fixture runs the
single-file CLI with both solvers, the bond parser with both solvers and its
metadata-free fallback, and the demonstration's changed reference/operator call
chain. It compares Arnoldi/Arpack energies and dense-oracle reference/norm values;
it does not execute the full demonstration or a molecular sweep.

The `Julia numerical regressions` workflow runs both commands on PRs and on pushes
to `main` or `research/ordering-mechanism-2026-10-09` that change this tree or its
workflow. Superseded PR runs are cancelled. It instantiates the declared project
on Julia 1.11 and records its resolved versions. Arpack 0.5 and ArnoldiMethod 0.4
compatibility bounds cover the eigensolver API families inspected for these regressions; Julia
1.11 is the hosted baseline. No full project Manifest is committed. The other
exploratory dependencies still resolve from the registry, so the recorded
`Pkg.status` output is essential provenance and these bounds are not a frozen
full dependency resolution. These checks use tiny matrices and no molecular
sweep; local component checks alone do not establish full-project CI parity.
