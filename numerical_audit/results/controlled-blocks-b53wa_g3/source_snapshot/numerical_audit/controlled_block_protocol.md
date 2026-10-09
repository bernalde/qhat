# Predeclared controlled whole-block experiment

Protocol `frozen-owned-block-permutations-v1`, declared before inspecting its simulation outcomes on 2026-10-09. The purpose is to isolate the benefit, if any, of signed-coefficient **sorting after ownership is fixed**, from the conservation benefit of keeping complete owned blocks together.

## Fixed construction and treatments

Use all eleven immutable tensor fixtures in `fixtures/manifest.json`: two geometries each for CH2O, H2S, HCN and SiH4, plus Li2, F2 and NH3. They represent seven molecular families, not eleven unrelated molecular systems. Construct each Hamiltonian with `molecular-unit-weight-v1`: no intermediate normal-order cutoff and a final consolidated Pauli coefficient cutoff of 1e-12. Retain its identity-free Pauli Hamiltonian for both exact and approximate evolution; report the scalar separately.

Call `build_molecular_ownership(..., ordering='signed_ascending')` exactly once per Hamiltonian. Freeze the resulting nonempty blocks, each block's owner, all final consolidated coefficients and each block's internal Pauli sequence. Empty parents are omitted. Never recompute ownership for a competing ordering. The owner's parent coefficient is the real coefficient of the lexicographically first canonical fermion component in that parent; it is a parent priority, not the sum of coefficients in its owned block.

Compare these 24 whole-block schedules:

1. Signed owner-parent coefficient ascending, with canonical component-key ties. This reproduces the frozen construction order.
2. Absolute owner-parent coefficient descending, with the same tie rule.
3. Retained owned-block coefficient l1 norm descending, with canonical component-key ties. This is a separately labeled diagnostic, not another definition of parent magnitude.
4. Reverse of the signed whole-block order, preserving the internal Pauli sequence of each block.
5. Twenty uniform random whole-block permutations. Base seeds are the consecutive integers 2026100901 through 2026100920 inclusive. Use NumPy PCG64 with `SeedSequence([base_seed, int(sha256(case_name)[:8],16)])`; save every resulting permutation explicitly. Reuse the same permutation at all times and step counts. These are draws with replacement: duplicate permutations, if any, are retained rather than rejected.

Include the final Pauli coefficient signed and magnitude orderings as two **external controls**. They change the contiguous block structure as well as chronology, so they do not isolate whole-block sorting. No treatment is chosen by its observed error.

## States, budgets and references

Use the same Hartree–Fock state within each comparison. Primary settings: T=0.25 and T=1, first-order formula, r=100 slices. Preselected convergence anchors: HCN_s1.00, Li2 and F2 at T=1 and r=50,200, alongside r=100. There are 26 schedules × (22 primary settings + 6 anchor settings) = **728 evolutions**.

The complete affine closure of all retained Pauli flips is an exact invariant subspace of **every individual Pauli**, so it retains all intermediate and final particle/spin leakage. Do not project into the physical particle sector. Validate representative whole-block and Pauli schedules against the full Hilbert-space implementation for one slice. Independently test the kernel, block exponentials, and ordering convention on literal dense Kronecker matrices. Check owned-block Pauli commutation and record full-space number and spin commutator coefficient bounds. Use sparse exact exponentiation of precisely the same retained Hamiltonian; crosscheck small references by dense eigendecomposition. Keep norm defects and number/spin leakage explicit.

All orderings receive exactly r times the retained Pauli count elementary Pauli exponentials. This is a matched algebraic budget, not a claim about compiled two-qubit gate count or wall time. Products are chronological: later listed factors multiply on the left.

## Metrics and planned summaries

Primary outcome: stable normalized full-state E = 1 - |<exact, approximate>|, computed from the phase-aligned residual. Also save infidelity, norm drift, and particle/spin leakage. No postselection or conditionally normalized sector metric is used to rank schedules.

For each Hamiltonian/time/budget report the deterministic errors, random median and empirical 25th/75th percentiles, and signed midrank among itself and the twenty random permutations. The signed midpercentile against random equals `(number of random errors smaller by >tol + half the tied count)/20`; smaller is better. An absolute tolerance of **1e-12** defines all ties. The interpretation floor is also **1e-12**: retain raw lower values, but omit ratios with a denominator at or below that floor. Report paired win/tie/loss counts separately for each primary time; retain case and molecular-family identifiers. State explicitly that repeated geometries, times, step counts and random permutations are dependent observations. Do not use 728 as an independent sample size, compute a naive population p-value, or claim the best random sample is a deployable algorithm.

The experiment is conditional on the signed-priority ownership partition and the HF state. It does not establish that this ownership assignment is optimal, test every state, or prove a universal accuracy advantage. The earlier random-state experiment remains relevant to that limitation.

## Reproducibility and stopping

Before tested propagation, save the protocol and source snapshots, package versions, fixture hashes, exact hexadecimal Hamiltonian coefficients, full fixed-block membership, owner weights, all schedule permutations/hashes, affine basis, and HF state hash. Save exact-reference state hashes and validation before comparisons. Each run writes a new directory; it never overwrites historical artifacts. Verify source hashes again at completion. A completion marker is written only after all predeclared rows exist. The default budget is one hour; a budget stop leaves an explicitly incomplete run and does not silently remove expensive cases. All computation uses one BLAS/OpenMP/Numba thread. The campaign remains local and unpublished.

Run with the isolated audit environment:

```sh
/tmp/qhat-audit-venv/bin/python -m qhat.numerical_audit.controlled_block_orders
```
