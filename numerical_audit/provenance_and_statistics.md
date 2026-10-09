# Provenance, fair comparisons, and statistical interpretation

Audit date: 2026-10-09. This report reconstructs supplied results from archived
code/data and adds a small independent physics check. Reconstructing archived
statistics does **not** validate every historical molecular calculation.
The broader frozen-input reruns are reported separately.

## 1. The supplied branch is not the numerical-study revision

- User-supplied `origin/improve_ordering`:
  `5bc6ded98fe6ab52182c24f03139832b5f6349e6`, already an ancestor of
  `origin/main` (`d3eea1ea458908f75d317794a80400a07896c2a8`).
- Albert's later study is on `origin/L-sweep`:
  `89c701652c4c64918825e279e84a7bea5924822a`, committed 2026-10-06.
  This is the archived revision used below.
- The working implementation is Python, **not** the Julia experiments or the
  `common/trotter_flattened.py` circuit path. Bugs in those other paths should
  not be attributed to these molecular figures without additional evidence.

The historical execution chain is:

1. `analysis/benchmark_L_sweep_trotter.py` loads saved spin-orbital tensors,
   normal-orders fermionic terms and pairs each with its Hermitian adjoint.
2. `analysis/benchmark_b2_signed_coefficient_baseline.py` sorts signed or
   magnitude parent coefficients and builds the first-occurrence Pauli order.
   `analysis/benchmark_b2_coloring_robustness.py` supplies the mapping and
   initial-HF commutator evaluator.
3. `analysis/benchmark_b2_active_spaces_matrix_free.py` applies Pauli
   exponentials to a full state vector; its exact reference is computed in a
   particle-number or spin sector.
4. `analysis/benchmark_fermionic_structure_ablation.py` evaluates frozen
   first-owner blocks and ablations. The 20-case runner is
   `analysis/run_cancellation_hypothesis_validation.py`.
5. `analysis/analyze_jw_magnitude_baseline.py`,
   `analysis/analyze_fermionic_aware_performance.py`, and
   `analysis/analyze_cancellation_hypothesis_validation.py` produce the tables
   and figures discussed in the supplied manuscript.

Relevant history: `eb8c5c5` added the performance analysis (August 18),
`de139a8` the full cancellation panel (August 23), and `14f3d63` changed its
primary comparison to JW magnitude (August 24). The tip's
`FERMIONIC_ORDERING_CONDITIONS.md` retains earlier best-of-two-JW statistics;
its prose does not match all current CSVs/figures. That is document staleness,
not an inconsistency in the two differently defined baselines.

## 2. What the archived counts actually say

The audit rejoined all 675 selected schedule records from 29 raw CSV files.
There are 230 case/settings pairs, with fermionic results for 215 and 201
comparisons above the original `1e-12` error floor.

| Comparator | Fermionic wins | Ties | Fermionic losses |
|---|---:|---:|---:|
| Fixed JW descending magnitude | 84 | 29 | 88 |
| Minimum of JW signed and JW magnitude | 59 | 14 | 128 |

Both rows reproduce the archive. The first is a legitimate comparison with
one specified heuristic. The second is a stronger **retrospective oracle
baseline**: it uses the measured exact errors to choose one of two JW
heuristics separately for each case. It is useful context, but is not a
cost-free deployable selection rule. Report the individual JW methods as
well; do not silently substitute one comparator for another.

These are selected, related active-space calculations, not independent
molecular samples or a population win-rate estimate. Also, equal Pauli
exponential counts do not establish equal compiled two-qubit-gate cost.

All matched records agree on the recorded Pauli count, qubit count, nominal
exponential count and exact-sector dimension. However, the historical join
key is only `(tensor_path, steps, time, tolerance)`. It has no tensor-content,
Hamiltonian-coefficient, initial-state, or source-revision hash. Equal paths
and counts cannot establish that the input coefficients were identical
across separate runs. The audit found no count mismatch in these selected
rows; it does not infer coefficient mismatch merely from missing provenance.

## 3. The metrics are not all invalid

For normalized states and a shared Hamiltonian, initial state and total time,

`E = 1 - |<exact|approximate>|`

is a legitimate phase-insensitive state-error comparison. Standard
infidelity `I = 1 - |<exact|approximate>|^2` ranks two approximations identically
because `I = E(2-E)`. Their **ratios** need not be exactly identical; report
which quantity is used. A Hartree–Fock-state error is not an operator-norm
certificate or an all-initial-states comparison.
The historical sector reference does not by itself erase Trotter leakage:
its inner product with the unrenormalized sector slice equals the full-space
inner product. Normalizing that slice separately would instead introduce
postselection and change the comparison.

The historical implementation clips an **unnormalized** overlap to `[0,1]`
and subtracts it from one, without saving state norms. Near the numerical
floor, accumulated norm drift and cancellation can affect the result. The
later September mechanism scripts instead normalize both states and use a
stable orthogonal-residual calculation. The historical global error floor
alone does not certify that every close win/tie is resolved.

For the initial-state BCH diagnostic

`R_pi = ||B_pi psi_HF|| / sum_{i<j} ||[K_i,K_j] psi_HF||`,

the denominator is common to all permutations of the **same final Pauli
list and coefficients**. This was checked in all 60 rows of the 20-case
panel, together with common time, step count, term count and work count.
Every row has zero fallback Pauli terms. Thus the reported within-case
relative `R` is not invalid because one ordering received a different
normalization: that denominator cancels exactly in the ratio.

The actual limitation is that this is an **initial-state, short-time
commutator diagnostic**, not the finite-time leading fidelity coefficient.
At fixed total time, the latter involves the propagated/projected BCH
error. Absolute `R` values across different Hamiltonians or different
splittings also have different denominators and do not define a universal
accuracy scale. The original evaluator additionally prunes individual
commutator amplitudes at the coefficient tolerance; its output is a
numerical diagnostic, not an exact or certified bound.

## 4. One headline correlation needs a narrower interpretation

The reported full-panel case-centered correlation is reproduced:
`r = 0.8346046649216678`, using 40 plotted alternative/reference rows.
There are exactly **two alternatives per case**, JW signed and JW magnitude.
Let

`x_j = log(R_JW,j / R_ferm)` and `y_j = log(E_JW,j / E_ferm)`.

Subtracting the within-case mean gives the two points

`(+1/2 log(R_signed/R_magnitude), +1/2 log(E_signed/E_magnitude))`

and their negatives. **The fermionic reference cancels completely.** The
independent reconstruction verifies this identity to `3.6e-16`.
Consequently, this particular statistic measures agreement between the two
JW alternatives; it is not independent evidence that the fermionic schedule
creates a better cancellation mechanism.

The 40 centered points also form 20 deterministic positive/negative pairs.
The ordinary Pearson `p = 2.20e-11` treats 40 rows as independent; it should
not be presented as calibrated independent-case evidence. The original
report already includes a much wider case-bootstrap interval. A fresh
10,000-resample case bootstrap gives approximately `[0.042, 0.961]` for the
same centered statistic. An ordinary correlation on the 20 unrepeated JW
contrasts is `0.872`; it has a different centering/estimand and is not a
replacement p-value for the original statistic.

The **direct** 20-case JW-magnitude/fermionic ratio correlation is separate
and does involve the fermionic reference: `r = 0.9025090632664486`. Shared
references can induce dependence, but using ratios does not automatically
make this association spurious: numerator and denominator are different
physical measurements. The JSON reports the exact covariance decomposition
rather than asserting that the common reference invalidates it.

The 20 cases represent only ten molecular families and were selected using
historical favorable/control performance. An exploratory molecular-family
bootstrap gives about `[0.467, 0.987]` for the direct correlation. This
remains a strong descriptive association in the archived panel, with broad
uncertainty; neither this bootstrap nor the original binomial p-values
convert it into prospective evidence of universal superiority. Selection,
repeated families, and the frozen-input provenance gap remain.

## 5. The later coherence metric is descriptive, not an accuracy score

The September scripts decompose final infidelity as

`I_pi = D_pi + X_pi = D_pi C_pi`,

where `D_pi` sums squared projected, propagated local defects. For `r` steps,
Cauchy–Schwarz gives `0 <= C_pi <= r` when `D_pi > 0`. Hence:

- `C < 1` means net destructive cross terms; `C > 1` means net constructive
  accumulation. A smaller positive `C` need not mean destructive interference.
- At the same step count, comparing `C` can describe different accumulation,
  but `D` also depends on the ordering. Smaller `C` alone does not imply
  smaller final error. Comparing `C` across step counts requires its scaling
  to be addressed; `C/r` is a useful descriptive normalization.
- Since `C` is computed using the measured final error, the exact identity
  `I = DC` does not by itself provide a predictive or causal explanation.
  The separately computed transported BCH prediction is a different test.

The 129-row recorded mechanism campaign contains 63 decompositions and only
one with `C < 1`; its maximum `C/r` is `0.989974`. The prospective campaign
contains 128 decompositions among 272 rows, none with `C < 1`, and maximum
`C/r = 0.997717`. Both satisfy the recorded `I = DC` identity to floating-point
precision. These are internal checks of archived outputs, not new evolution
runs.

## 6. The archive already records important negative findings

`codex_first_experiment/.../SUMMARY.md` reports that reproducible fresh F2
and NH3 runs disagreed with earlier historical values. F2's reported error
changed from `5.86e-11` to `1.06e-7`; NH3's Pauli count changed from 5064 to
5098. These differences establish unresolved historical lineage, not their
cause. The first campaign's input directory was not committed, but the later
mechanism campaign includes the three frozen tensor files needed to rerun
F2, NH3 and Li2.

The later `codex_mechanism_experiment` records F2 ordering reversals with time
and a transported BCH predictor that explains them. The
`codex_external_validation` archive includes eight frozen 10-qubit molecular
Hamiltonians and 272 settings; its report says JW signed wins six of eight
anchor comparisons, fermionic signed wins two, and random whole-block
reordering sometimes improves on signed ordering. These are valuable
controls and negative results. They must be independently rerun before
being adopted as newly validated paper results; the broader audit handles
those reruns separately.

## 7. Independent check of the historical physics path

A fresh 10-qubit HCN calculation used the committed tensor
`09ede69f377f917bda17a6cb1a4a95ff855594cbaff4cf4e93d5214d0bdc6093`.
Independent OpenFermion sparse matrices were used rather than QHAT's Pauli
bit-action routine.

- Number-sector (dimension 210) and spin-sector (dimension 100) exact
  references, embedded in the full Hilbert space, agree with full-space
  sparse exponential evolution to at most `3.8e-16` in vector norm.
- The historical exact and Trotter references differ by scalar energy
  `alpha = -3.4650853304313642`. At `T = 0.37`, raw vector distance is
  `1.196`, while removing the scalar phase leaves about `2.0e-14`.
  This is a concrete demonstration that the old raw-vector-error column
  should not be used to measure Trotter error. The modulus-overlap metric
  is unaffected by this scalar phase.
- All three physical ordering lists agree exactly with independent
  chronological sparse Pauli-exponential action for one step. Reversing
  chronology changes these states by `5.6e-4` to `1.1e-3`; the convention
  must therefore be explicit in any BCH comparison.
- A separate arbitrary complex four-qubit state including odd-Y Pauli
  strings agrees to `9.4e-17`, checking phases not exercised by purely real
  molecular terms.

This supports the archived simulator's sign, chronology, bit convention and
sector embedding on the checked examples. It does not validate historical
input provenance, all numerical floors, or the performance claim.

## Reproduction

`provenance_statistics.py` reconstructs the historical counts, common
normalizations, correlation identities, and exploratory bootstraps;
`provenance_statistics.json` contains source hashes and numerical results.
`historical_pipeline_oracles.py` and `historical_pipeline_oracles.json`
contain the independent sparse-matrix checks and environment versions.

Both scripts take an archived `L-sweep` checkout path and a **new** JSON
output path. They refuse to overwrite the output. For example, from this
repository:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/qhat-audit-venv/bin/python \
  numerical_audit/provenance_statistics.py /tmp/qhat-albert-study \
  /tmp/qhat-statistics-new.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/qhat-audit-venv/bin/python \
  numerical_audit/historical_pipeline_oracles.py /tmp/qhat-albert-study \
  /tmp/qhat-oracles-new.json
```

The archived study checkout is read-only for these tests. Code and tensor
hashes in the outputs identify exactly what was checked. No historical CSV,
figure or manuscript was overwritten.
