# Correlated and coherent state/time holdouts (issue #4)

This protocol is fixed before panel ordering errors. It extends the derived
diagonal-response placement rule, not the historical signed-parent heuristic.
It does not revisit the failed two-level reference correction or the completed
factorial, Tranter, or 8! studies. Atomic units: coefficients in Hartree,
T in inverse Hartree, hbar=1. Fixed orbital basis and retained target per family.

## Inputs and states

Reuse all three authenticated prospective-reference models: Li2_r2p95 (6
qubits/2 electrons), F2_r1p55 (8/6), H2O_scale1p10 (8/4). No input replacement.
The source model, previous HF design and target hashes enter the freeze.
The number/spin sector has the same N_alpha,N_beta as the HF determinant.
Diagonalize its real retained H block for **state preparation only**. Select
the lowest two eigenvectors. Resolve a degenerate eigenspace deterministically
by projecting occupation basis vectors in increasing occupation order and
Gram-Schmidt (degeneracy threshold 1e-11 Hartree). Canonical phase: largest
absolute component, smallest occupation tie, positive real value. Record
eigenvalues, gap, residuals, participation ratio and maximum determinant weight.
Also record the one-particle density matrix occupations and Tr(gamma-gamma^2).
Fractional occupations (defect >1e-8) diagnose correlation beyond a single
Slater determinant; participation ratio alone diagnoses only multiple
determinants in the fixed basis. Keep these labels distinct.
Call these sector ground and first-excited states, not unrestricted global
ground states. Report determinant-like cases (participation ratio <=1+1e-8)
without replacing them. States are g, e, (g+e)/sqrt(2), (g+i e)/sqrt(2).
These last two are distinct fixed coherent phases, never an incoherent mixture.
Freeze explicit state vectors and their hashes before any circuit error.
Numerical eigensystem/state construction CPU and wall cost are separate;
this exact-diagonalization preparation is not an efficient hardware routine.

T = 0.37 and 1.73. Ideal Pauli rotation caps B=10U (primary) and 30U, where U
is the retained Pauli count. No interpolated accuracy or step counts.

## Objective, orders and controls

D is the unchanged full diagonal part. For any normalized complex psi, use
J_D(pi,psi,T)=||(I-psi psi*) Q_pi psi||^2, where
Q_pi=integral_0^T exp(itD) C_pi exp(-itD) dt and
C_pi=-i/2 sum_{a<b}[G_pi(a),G_pi(b)] (listed-order convention).
Q_wv=C_wv T exp(i T(D_w-D_v)/2) sinc(T(D_w-D_v)/(2 pi)).
Add complex amplitudes over all input occupations before projection/squaring.
Do not average determinant scores. Opposite chronological C sign has the
same variance. The generalized rule uses binary64 complex arithmetic; it
does not inherit the determinant implementation's exact integer comparison.

Nine arms, all implementing exactly the same retained nonidentity target:

1. hf_current: unchanged frozen HF current-placement schedule trained at T=1.
2. state_current: generalized state-aware placement, same l1 off-unit priority,
   best contiguous-D start, one diagonal sweep and earliest ties.
3. state_JD: conserving two-start insertion, label and l1 starts, position cap
   4m(m-1), at most 10 moves/start.
4. state_pauli_JD: individual-Pauli insertion, same starts/rule, cap 4P(P-1).
   This is an explicitly larger position allowance, not matched CPU.
5. hf_pauli_large: unchanged previous HF larger-search Pauli schedule.
6. hf_symmetric: unchanged previous HF conserving J_D order, symmetric execution.
7. state_symmetric: symmetric reuse of state_JD; no second-order optimization.
8. jw_signed: previous frozen HF individual-Pauli signed order.
9. jw_magnitude: previous frozen HF individual-Pauli magnitude order.

No state-aware search uses full-H dynamics or measured circuit errors. Its input
state may itself cost full-H state preparation, which is charged and disclosed.
Placement keeps complete internally commuting conserving units; Pauli controls
may leak. First-order steps=floor(B/U). Symmetric steps are the largest integer
with exact central and identical-boundary unit fusion, including orientation-
dependent first/last unit lengths. Retain full 2^n states for all propagation.
Freeze all 216 schedules (3 families x 4 states x 2 times x 9 arms) together
before the 432 circuit errors (two budgets) are evaluated.

## Coverage and reporting

Error E=1-|normalized overlap| via a stable residual; also report infidelity,
norm, number leakage and spin-projection leakage. Absolute tie and denominator
floor=1e-12. Descriptive ratios omit denominators <=floor and count exclusions;
all absolute errors, strict losses and floor-limited cases remain in tables.
Coverage for each state/time stratum requires hf_current E>floor at 10U in at
least 2 of the 3 families. A failed stratum is labelled inconclusive; no
replacement trial or claim of coverage for it. Report each family separately
and W/T/L for each arm versus hf_current. These are descriptive related inputs,
not independent statistical replicates or a universal transfer guarantee.
Report all setup/search/preparation costs and actual evaluations. This study
does not establish matched classical efficiency (issue #5) or hardware cost.

## Independent checks and finite resources

Synthetic complex-state tests compare response to dense quadrature and
incoherent-mixture counterexamples, determinant reduction, insertion updates,
degenerate-state selection, boundary fusion, invalid dimensions and budgets.
Independent panel verification builds literal Pauli operators, chronological
BCH commutators, dense D response and an unfused circuit propagation, checking
all schedules, all scores and all 432 state/error/leakage rows. Verify state
labels/eigen-residuals and source/input/state/schedule/result hashes.

One numerical/build/test worker at a time; one BLAS/Numba thread and build job.
CPUQuota=25%, Nice=19, CPUWeight=5, IOWeight=10; MemoryHigh=448M,
MemoryMax=512M, MemorySwapMax=0. Start/resume only above 5.5 GiB available;
below 5 GiB freeze only the owned service. Guard polls every 3s (at most 0.75
CPU second allowance overshoot). Do not alter unrelated processes.

Finite phase budgets: tests 120 CPU s/900 wall s/64 MiB outputs;
representative synthetic end-to-end preflight 120 CPU s/900 wall s/64 MiB;
panel construction+propagation 1800 CPU s/10800 wall s/512 MiB;
independent verification 1800 CPU s/10800 wall s/512 MiB. Extraction, reporting
and replay allowance 300 CPU s/1800 wall s/128 MiB separately. Preflight uses
the largest model's structure, deterministic altered coefficients and seeded
synthetic sector states, never panel error selection. Include extraction,
packaging and replay estimates with at least 2x time/memory/disk margin before
starting the panel. If essentials do not fit, retain partial artifacts and name
the gap; do not enlarge budgets or substitute a smaller scientific question.

Pin the existing Python 3.12 environment; no new runtime dependency. Preserve
failed attempts in fresh directories. Refuse n>8 or P>128 before exponential
allocation and refuse pair-vector storage above 64 MiB. New scientific-source
bytes require a new freeze or an explicit provenance-preserving correction;
never silently overwrite frozen evidence.
