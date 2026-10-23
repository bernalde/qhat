# Prospective physical-geometry decision experiment

Declared 23 October 2026 before generating these tensors or measuring their
Trotter errors. This is a bounded comparison, not a broad validation campaign.
The preceding mechanism studies are immutable and are not rerun here.

## Inputs and conventions

Generate neutral singlet RHF/STO-6G molecular integrals, with PySCF symmetry
disabled, convergence tolerance 1e-11, at these fixed geometries (Angstrom):

* Li2: Li=(0,0,0),(0,0,2.95); two active electrons in three spatial orbitals.
* F2: F=(0,0,0),(0,0,1.55); six active electrons in four spatial orbitals.
* H2O: O=(0,0,0), H=(+/-0.8332948775591,0,0.6457876002479);
  four active electrons in four spatial orbitals.

Freeze the lowest doubly occupied orbitals; take the active window immediately
above them. Retain orbital coefficients, occupations, energies, SCF convergence,
and AO/MO orthogonality. Explicitly expand spatial integrals into alternating
alpha/beta spin orbitals, including the InteractionOperator two-body factor
1/2, without OpenFermion's integral cutoff. Apply the existing exact real CAR,
signed-adjoint averaging, whole-parent cutoff 1e-12 and exact JW consolidation.
Omit the scalar identity consistently. Do not prune individual Pauli words.
Record the Hermitian adjustment and discarded-parent bound separately.

Mode q is bit 1<<q in this experiment's occupation-vector indexing. Independently
permute standard tensor-product matrices for verification. This does not change
the physical spin-orbital order or perform orbital optimization.

Each geometry contributes HF and the determinant of smallest Hamming distance
from HF in the same (N_alpha,N_beta) sector, breaking ties by occupation integer.
These six points are related within three molecular families, not six independent
molecular samples. They are determinants, not correlated-state validation.
Physical time T=1 in inverse Hartree (hbar=1). Rotation budgets are 5,10,20,50
times the retained nonidentity Pauli-word count U, without hardware gate claims.

## Candidate and fixed correction

Retain complete nonzero flip-mask units and separately movable diagonal words.
Fix the off-unit order by decreasing sum of absolute integer coefficients,
then flip mask. Visit diagonal terms by decreasing absolute coefficient, then
canonical Pauli word. The current rule uses the unchanged exact supplied-weight
one-pass code from the completed diagonal-rule experiment: start at the best
contiguous-D slot, visit each diagonal term once, scan all off-unit slots,
and choose earliest ties. Its source hash is retained.

The correction replaces reference D by A=D+h(|u><v|+|v><u|), where u is the
determinant connected to input v with largest absolute H_uv (smallest flip-mask
tie), and h=H_uv. If no coupling exists, use D. This reference is only a scoring
approximation; all physical circuits implement the same retained H. Choose u
before search, never from exact dynamics, response under H, or physical errors.

Evaluate Var_v(integral exp(itA) C_pi exp(-itA) dt) using just the two input
columns of C_pi, a 2x2 eigensystem, and diagonal energies at sparse destinations.
Do not diagonalize H to construct an ordering. C_pi=-i/2 sum_{a<b}[G_a,G_b]
in the listed-order convention; its opposite chronological sign has the same
variance. Sum coherent amplitudes before squaring, and remove the component
parallel to v. Starting slot, priority, one pass and earliest-slot policy match
the current rule. Corrected complex response scores use binary64 arithmetic;
the current rule retains its exact supplied-weight comparison guarantee.

## Comparators and budgets

Alongside current and corrected placement, run two-start best-improvement
insertion for conserving J_D, conserving J_0 and individual-Pauli J_D. Starts
are unit-label order and decreasing coefficient-l1 order. Cap their combined
examined positions at E=(k+1)(m+1), the current placement scan count, where k is
the number of diagonal terms and m is the number of off units. Include every
examined position, even a zero commutator, and allow ten accepted moves per
start. This matches positions, not wall time or arithmetic per position.

Also run an explicitly larger-budget individual-Pauli J_D control: the same
two starts and ten-move cap, at most 20 P(P-1) positions for P Pauli units.
This stronger control has unequal classical work and must be labelled as such.
Include fixed JW signed and magnitude orders (coefficient/word ties), plus
symmetric execution of the conserving J_D insertion order. The symmetric
control reuses a first-order order; it is not a second-order optimization.

All nine arms retain exactly the same target. First-order steps are floor(B/U).
Symmetric steps use the largest feasible count with exact fusion of identical
boundary units; charge actual Pauli rotations and record infeasible budgets.
The singleton-unit case is exact and uses one merged exponential per word.
Keep complete unit boundaries in conserving arms. Individual-Pauli controls
may leak; retain the full 2^n Hilbert space in every physical simulation.

Freeze all schedules for all six points, every active source/dependency hash,
input/state identity, score and search ledger before any new molecular error.
No schedule may be selected or changed using those errors.

## Measurements, verification and decision

Record stable normalized E=1-|overlap|, infidelity, norm drift and leakage
from number and spin projection. Absolute tie and denominator floor: 1e-12.
Record all results, including losses and floor exclusions, and distinguish
setup, search, generation/reference and propagation costs. One timing observation
under a shared-machine quota does not establish a wall-time speedup.

Independently check tensors/HF contraction, retained-target CAR versus literal
Pauli matrices, completed-unit conservation, every schedule sum and rotation
budget, response scores by dense reference integration on these small systems,
and all new physical circuits with a distinct propagation implementation.
Full-H exact dynamics is used only after schedule freeze for validation.
Synthetic tests cover the two-level formula, its zero-coupling reduction,
all small permutations, placement, insertion and symmetric pulse accounting.
Refuse n>8 and invalid work budgets before exponential allocations.

Primary comparison is at B=10U. A useful correction must win at least four
of six points against current placement, have no above-floor loss exceeding
2x, and median additional reference/setup/search CPU cost no greater than 2x
current construction. Report each condition independently and all other
budgets. This is a development decision rule, not a significance test or
universal guarantee. Compare both candidates against every strong control,
including the unequal-search control and symmetric formula.

Do not infer generalization to correlated states, larger spaces, different
encodings, hardware costs or second-order optimization. A failed correction
is retained rather than replaced after observing the holdout outcomes.

## Resources and review/publication

One numerical worker/thread, CPUQuota=25%, Nice=19, CPUWeight=5, IOWeight=10,
MemoryHigh=3G, MemoryMax=4G, MemorySwapMax=0. Start above 6 GiB available;
freeze only this service below 5 GiB and resume above 6 GiB. Leave other jobs
and the JW/BK project untouched. Fresh outputs per changed scientific source.

Review the local git trees with Astra and Argon. Retain original reports and
per-finding dispositions separately. At most three review rounds and three
verification passes. Publish PRs only after local convergence; merge through
the separately authorized merge workflow with live gates. Then create only
still-open future-work issues, traced to slide 57 and current evidence.
