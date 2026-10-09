# Explicit coefficient precision and first-owner construction

This is an opt-in local correction. The five imported study files retain their
historical contents, pinned in `source_provenance.json`. Results from different
coefficient policies must be labeled separately, even when their Hamiltonians
agree to near machine precision: a change to a tiny coefficient or a sorting tie
can change first ownership and hence the product formula.

## Reproduced precision discrepancy

The historical scripts advertise a default coefficient tolerance of `1e-12` and
call `normal_ordered(operator)` followed by `compress(abs_tol=tolerance)`. In the
installed OpenFermion implementation, however, `SymbolicOperator.__iadd__` calls
`_issmall` with a default argument bound to `openfermion.config.EQ_TOLERANCE`,
which is `1e-8`. Normal ordering, Jordan--Wigner mapping, and operator sums can
therefore prune contributions before the final requested cutoff is applied.
Changing the config variable afterward would not reliably change default
arguments already bound by Python; this audit does not monkey-patch it.

Concrete observations:

```python
h = FermionOperator('0^ 0', 5e-10)
clean_fermion_operator(h, 1e-12).terms  # {}
jordan_wigner(h).terms                # {}
```

The coefficient is 500 times the requested cutoff, but the nonzero number
operator is lost. This is not equivalent to one clean final `1e-8` cutoff:
contributions can be discarded separately before several small contributions
combine to a larger final coefficient.

The same issue affects the legacy fermionic graph. For
`F = 1e-5 n_0` and `G = 1e-5 (a_0^dagger a_1 + a_1^dagger a_0)`, the true
commutator has operator norm `1e-10 > 1e-12`, but
`fermionic_terms_noncommute(F, G, 1e-12)` returns false. That graph is therefore
not an exact commutation certificate at its advertised tolerance. The opt-in
sorted first-owner construction does not use this graph. The historical graph
code is unchanged to retain a reproducible comparison, and its failure is a
strict expected-failure test.

## New coefficient API and bounded scope

`analysis/fermionic_coefficients.py` exports:

- `fermion_from_interaction(interaction)`: copy the scalar and one-/two-body
  tensors into raw fermionic coefficients, dropping only exact zeros. It follows
  OpenFermion's tensor convention without inserting another factor of one half.
- `normal_ordered_molecular(operator, *, atol)`: normalize unit-weight words,
  accumulate their actual coefficients directly, then apply the requested final
  cutoff.
- `jordan_wigner_molecular(operator, *, atol)`: map unit-weight words and combine
  their actual coefficients before applying the final Pauli cutoff.

Only input words of at most four ladder operators are supported. Unit-weight
normal-order expansions have integer coefficients. For a word of length at
most four, nonzero unit-weight JW coefficients are dyadic and at least `1/16`
in magnitude, well above OpenFermion's implicit arithmetic cutoff. Caching those
unit transforms is therefore safe from the particular small-coefficient pruning
identified here. Long words, nonfinite coefficients, and invalid cutoffs are
rejected rather than silently generalized.

Physical coefficients are combined with `math.fsum` separately on real and
imaginary parts. Results are inserted directly into operator dictionaries, not
through `+=`. This controls truncation placement; it does not turn floating
point arithmetic into an exact proof. Apply JW directly to the raw tensor
operator when no intermediate normal-order cutoff is intended.

## New ownership API

`analysis/fermionic_ownership.py` exports
`build_molecular_ownership(raw_operator, *, atol, ordering='signed_ascending')`.
Supported priorities are raw insertion, signed ascending, and magnitude
descending, with canonical component keys breaking ties. It returns:

- the untruncated normal-ordered Hamiltonian and canonical Hermitian parent
  pairs, preserving their supplied coefficients;
- the full JW scalar separately and the consolidated, final-cutoff nonidentity
  Pauli coefficients;
- parent priorities, actual owned Pauli blocks, the chronological flattened
  Pauli order, and any fallback Pauli keys;
- discarded coefficient weight, Hermitian/imaginary residuals, owner-block
  pairwise commutation, and particle-number commutator diagnostics.

Input canonical components must conserve number and be real up to the reported
absolute tolerance. Adjoint mismatches above that tolerance raise an error.
Smaller discrepancies are retained and reported; there is no Hermitian averaging
or real projection. Owner discovery uses nonzero parent contributions without
an intermediate threshold. Any subsequent real projection by an evolution
caller must be explicit and its imaginary residual inspected.

Independent Pauli truncation need not preserve exact particle-number symmetry.
The implementation therefore reports the consolidated Pauli coefficient l1 norm
of `[N,H]`, as well as the largest analogous residual among its actual owned
blocks. These are operator-norm upper bounds and numerical checks, not exact
Lean certificates. The total discarded Pauli coefficient l1 sum likewise bounds
the operator norm of the discarded nonidentity operator.

## Frozen-input measurements

`coefficient_precision.json` compares both policies on the eleven supplied
fixtures. NH3 changes from 5,098 to 5,124 final nonidentity Pauli terms (28 added,
2 removed), with coefficient l1 difference about `3.620837e-7`. The other ten
fixtures retain the same Pauli keys and differ by less than `5.3e-15` in
coefficient l1 sum. These measurements alone do not establish whether a revised
ordering improves the simulation error.

`ownership_precision_checks.json` records a separate direct test of the new
ownership API on all eleven raw tensor inputs, together with source hashes.
Every input builds successfully, with zero fallback terms and zero owner blocks
containing an anticommuting Pauli pair. Hermitian-parent l1 residual sums are at
most `3.53e-14`; they have not been averaged away.

For the explicit final `1e-12` cutoff, the largest owned-block number-commutator
l1 residual is `2.5592534e-13` (NH3), and the full NH3 Hamiltonian residual is
`1.0885511e-12`. Its discarded coefficient l1 sum is `1.1030927e-11`.
These nonzero values are disclosed rather than described as exact conservation.
Other owned-block residuals are at most `7.48e-16`. The previously inspected
legacy NH3 Hamiltonian was also numerically number-conserving; the hidden cutoff
finding does not establish a measurable symmetry failure in the old figure data.

## Independent oracles

`test_molecular_oracles.py` gives 35 passing tests and two strict expected
failures for the unmodified historical cutoff and graph issues. The tests cover
all five imported modules, including serial/parallel Numba Pauli action and
chronological first-, second-, and fourth-order evolution against independent
SciPy exponentials; contractions and scalars; actual owned-block conservation;
within-owner permutations; original-parent versus owned-block products; and
Hartree--Fock BCH norm evaluation against explicit dense commutators.

The new coefficient/ownership tests use occupation-basis ladder matrices as an
independent fermionic oracle, including coefficients below `1e-8`, cancellation
before the requested cutoff, tiny Hermitian discrepancies, and a four-mode
restriction of actual H2S integral data. Tiny-coefficient tests use absolute
matrix tolerances down to `1e-25`; the integral test separately permits
machine-relative summation roundoff from its approximately -400 Hartree scalar.

Execution uses `/tmp/qhat-audit-venv/bin/python`, with one BLAS/OMP/Numba thread
and `NUMBA_CACHE_DIR=/tmp/qhat-audit-numba`. Full output is retained in
`molecular_test_results.txt`. Framework fixes and their independent
359-test result are documented separately in `framework_findings.md`.
