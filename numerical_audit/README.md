# Fermionic ordering: local implementation and numerical audit

Audit date: 2026-10-09. Branch: `audit/fermionic-ordering-2026-10-09`.
All changes and results are local. No branch was pushed or PR opened.

**The overlap error is a meaningful comparison when the Hamiltonian, initial
state, time and work count are fixed. The supplied evidence does not establish
that signed fermionic ordering is generally more accurate.** The strongest
supported benefit is conservation by complete owned blocks under the exact
mathematical hypotheses. Accuracy depends
on the Hamiltonian, state, time and comparator.

## What was audited

The requested branch `improve_ordering` (`5bc6ded`) is already an ancestor of
current `main` (`d3eea1e`). The new branch was created from the former and then
fast-forwarded to the latter. The actual later molecular study is on `L-sweep`
(`89c7016`), so five of its benchmark modules and eleven frozen tensor fixtures
were imported with original hashes. See [source provenance](source_provenance.json)
and the [fixture manifest](fixtures/manifest.json).

The molecular benchmark uses matrix-free Pauli evolution. Separate bugs found
in the circuit framework must not be blamed for the paper's figures: they use
different execution paths.

## Main findings

1. **The principal error metric is valid, but its old evaluation is fragile.**
   For normalized full states, `E = 1 - |<exact, approximate>|` and infidelity
   `I = E(2-E)` rank methods identically. The replacement evaluates E from the
   phase-aligned residual, records norms and rejects significant norm drift.
   It retains leakage outside the physical sector. The old subtraction can
   return zero near its floating-point floor; clipping an unnormalized overlap
   can also conceal invalid propagation. See [metric evidence](metric_findings.md).
   In the fresh NH3 run at T=0.25, r=100, the old formula reverses the
   fermionic-versus-JW-magnitude winner even though both stable errors exceed
   the historical 1e-12 analysis floor.

2. **The BCH normalization is common across the compared Pauli permutations.**
   This was checked in the archived 20-case panel. It does not make those
   within-Hamiltonian comparisons unfair. However, the initial-state BCH norm
   is a short-time diagnostic, not the fixed-total-time error coefficient.
   F2 reverses its fermionic-versus-JW-signed ranking between T=0.25 and T=1.
   The transported BCH coefficient is independently verified against dense
   quadrature and finite-step convergence.

3. **One correlation has been overinterpreted.** Centering the two
   JW/fermionic log ratios within each case cancels the fermionic reference
   exactly. The resulting 40 plotted points are 20 dependent positive/negative
   JW contrasts, not 40 independent tests of the fermionic mechanism. The
   separate direct ratio correlation remains a descriptive association.
   See [provenance and statistics](provenance_and_statistics.md).

4. **The comparator matters.** Archived counts reproduce as 84 wins, 29 ties
   and 88 losses against fixed JW magnitude; against the retrospectively best
   of JW signed and magnitude they are 59, 14 and 128. Choosing the best using
   exact measured errors is an oracle baseline, not a free selection algorithm.
   These selected, related cases do not give a population win rate.

5. **Conservation does not imply a state-error advantage.** On eight frozen
   ten-qubit Hamiltonians at T=1 and r=100, fermionic signed wins two HF-state
   comparisons; JW signed wins six. For five predeclared random states per
   Hamiltonian in the same physical sector, JW signed wins 38 comparisons and
   JW magnitude two. Fermionic blocks still preserve the sector to roundoff.
   These are exploratory Haar-state controls, not chemically representative
   ground-state sampling. See [state sensitivity](results/initial-state-ergijt7d/report.md).

6. **The advertised coefficient tolerance was not the effective tolerance.**
   OpenFermion arithmetic can discard terms around its internal 1e-8 cutoff
   before the requested 1e-12 compression. The new opt-in molecular transform
   and ownership path applies a single final cutoff after accumulation.
   On frozen NH3, the nonidentity Pauli count changes from 5098 to 5124;
   the coefficient-l1 difference is 3.62e-7. The other ten Hamiltonians differ
   only at roundoff scale. Original data and corrected-input runs are kept
   separate. See [precision diagnostics](coefficient_precision.json) and
   [ownership diagnostics](ownership_precision_checks.json).
   At T=1, r=100, none of the eleven cases changes its winning method under
   this correction; the largest relative change in E is 1.38e-6. This defect
   does not explain the large earlier F2 discrepancy. See the
   [explicit-input comparison](reports/precision_comparison.md).

7. **Several framework defects are repaired.** These include repeated
   single-term evolution for the wrong total time, incorrect circuit phase and
   identity-wire treatment, obsolete tensor hooks, BK reference-state encoding,
   CLI configuration loading and non-bijective eigenstate matching. See
   [framework findings](framework_findings.md) and [metric findings](metric_findings.md).

## New calculations and validation

- Historical retained-Hamiltonian campaign: eleven frozen Hamiltonians,
  three orderings, T in {0.25, 1}, r in {25, 50, 100, 200}: **264 comparisons**.
  Hamiltonian and all saved ordering hashes must match before evolution;
  regenerated baseline orders must match as well. Exact references use the
  same retained Pauli sum. Affine reduction is closed under every Pauli factor
  and retains physical-sector leakage. See the [run](results/frozen-ntrgzhgz/)
  and [analysis/figures](reports/frozen-ntrgzhgz/REPORT.md).
- Explicit-cutoff construction: **33 additional comparisons**, eleven cases,
  three orderings, T=1, r=100. See [run](results/explicit-n0x4vv3g/).
  This changes coefficient arithmetic and can change exact sorting ties;
  it is not a rerun with an identical historical Hamiltonian.
- Initial-state sensitivity: **144 comparisons**, HF plus five fixed random
  states on each of eight ten-qubit cases. Full 1024-dimensional propagation
  is checked against independent sparse-Pauli and sector-eigensolver references.
- **501 tests passed; two strict expected failures** retain executable
  witnesses for the unmodified historical implicit-cutoff and graph-commutator
  bugs. The corrected opt-in molecular path passes its tests. Independent
  matrix oracles cover CAR contractions/scalars, JW mapping, consolidated
  ownership, number conservation, parent-versus-owned products, first/second/
  fourth-order formulas, complex states, leakage, and BCH chronology/projection.
  See [test output](test_results.txt).

The unchanged historical graph routine is not used by the corrected sorted
first-owner path. Do not use it as an exact commutation oracle at the requested
tolerance. The explicit construction reports floating Hermiticity and number-
commutator residuals; it does not turn floating-point truncation into a Lean
certificate of exact conservation.

## Reproduce

Use Python 3.12 and an editable checkout. The focused environment is recorded
in [environment.json](environment.json); the current interpreter is
`/tmp/qhat-audit-venv/bin/python`. Numerical/circuit audit dependencies are in
[requirements-audit.txt](requirements-audit.txt). This was not a full QHAT
installation validation: installed pyLIQTR 0.3.3 is below the project's declared
minimum 1.4.2. The frozen numerical path does not use that integration.

From the QHAT checkout:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_CACHE_DIR=/tmp/qhat-audit-numba
export MPLCONFIGDIR=/tmp/qhat-audit-mpl
python -m qhat.numerical_audit.run_frozen_audit
python -m qhat.numerical_audit.run_explicit_precision
python -m qhat.numerical_audit.initial_state_sensitivity
python numerical_audit/summarize_frozen_audit.py NEW_FROZEN_RUN --output NEW_REPORT_DIR
```

Each evolution runner writes a fresh directory, inputs/settings/source hashes,
incremental results and a completion/failure record or final summary. Completed aggregated
results are committed; redundant incremental per-row checkpoints may be ignored.
No SCF or molecular integrals are regenerated. Historical statistics additionally
need a checkout of `89c7016`; the scripts take that directory explicitly.

## Limits and consequences for the paper

The 201 historical valid comparisons have been rejoined and checked for recorded
settings, not all re-evolved. Their CSV keys lack Hamiltonian-content hashes.
The archive itself records earlier F2 error and NH3 term-count discrepancies;
the present reruns validate the later committed tensors, not the missing earlier
input lineage. The exact cause of those earlier differences remains unresolved.

The eleven frozen Hamiltonians represent seven molecular families. They support
reproducibility and counterexamples to universal ranking claims, not a broad
performance survey. Equal Pauli exponential counts do not guarantee equal
compiled gate cost. Transported BCH prediction is omitted for the largest NH3
sector to keep memory bounded. None of these floating-point experiments replaces
the exact mathematical hypotheses in the manuscript's Lean proofs.

The paper should separate structural conservation, state-error performance and
prediction. It should report both individual JW baselines, leakage and total
state error, fixed-time transported BCH tests, explicit coefficient policy and
input hashes. Interpret the centered correlation as a JW contrast; avoid claims
of universal signed-ordering superiority or of destructive accumulation from a
coherence factor that is merely smaller than another factor. The manuscript and
Lean sources were not edited during this numerical audit.
