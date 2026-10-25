# Quantitative molecular guard and finite-step certificate attempt

This is the numerical component of paper issue18. The generic sharper-epsilon
and finite-step Lean statements live in the separate paper proof worktree. They
remain unverified drafts until kernel checking; software interval certificates
must not be described as kernel-certified molecular proofs.

## Fixed input and orders, before any guard or Trotter-error calculation

Use the newly generated LiH_r1p713 input from the authenticated #17 pipeline:
STO-6G,Li at(0,0,0),H at(0,0,1.713)Angstrom,freeze spatial orbital0,active
orbitals1,2,3,six spin modes,two active electrons. The scalar-removed model JSON
SHA256 is ea38b9730eed2faf55fbd1fe22722d93b858b95fe4d804ecb72855d10176c134.
The chemistry producer is7c8c64acf0f13703df10b0335e0017da84e94f9b; its exact
source/environment/tensor/orbital/coefficient records are retained by #17.
This is a pinned development input until the actual prerequisite publication
and merge are verified. No input or order pair will be replaced after errors.

Compare the existing legacy current-placement order with the diagonal-prefix
order formed by the same established priorities: diagonal units by decreasing
absolute coefficient then word; nonzero complete flip classes by decreasing
coefficientL1 then mask. The current-placement algorithm retains its fixedT1
initializer rule. Freeze both actual permutations and integer unit actions before
computing any leading guard or physical error. Retain both comparison directions;
reference-response ordering does not choose a new pair or remove losses.

The written product is exp(-ihG1)...exp(-ihGm). State propagation therefore applies
its units in reverse written order. Use first-order product formulas only here.
The initial state is the HF determinant. Check every complete unit's exact
integer action on all64 occupation states before restricting to the invariant
Nalpha=Nbeta=1 sector of dimension9. The full-space invariance check is required;
this is not a leakage-removing projection of an approximate Pauli trajectory.

## Full predeclared panel

Physical time is the exact rational3/5 in atomic units(hbar1). Write
H_epsilon=D+epsilonV with D the unchanged diagonal units and V all complete
nonzero flip classes. Attempt the physical epsilon1 first. Retain the entire
scaled grid epsilon=2^-k for k=0,2,4,6,8,10,12,14,16,18,20, including every vacuous
or failed guard. A scaled result concerns that actual modified Hamiltonian and
cannot certify the physical epsilon1 molecule.

Use r=1,2,4,8,16,32 for both fixed schedules at every epsilon:132 finite circuits
and11 full-H reference evolutions. Report every interval and pair, not a selected
successful subset. No circuit error selects epsilon,T,r,an order,or a new input.

Derive actual K(epsilon)=epsilonK1+epsilon²K2 from the ordered commutator sum with
its physical -i/2 factor. Compute the diagonal transported K1 response coherently,
then project away the HF ray. Bound each a=norm(projected response) by directed
intervals. Use rigorously checked operator-norm upper bounds, permitting the
conservative sum of absolute matrix entries, in
C=T*normK2+T²*normV*(normK1+normK2). Record which direction,if any,has disjoint
intervals for the true response standard deviations. Zero gaps and unresolved
intervals are explicit outcomes, not successful certificates.

## Constructive finite-step bounds

Use96-term exponential Taylor polynomials, scaling any generator to norm at most
1/2. The generic draft bound is2*(1/2)^96. Propagate approximation errors through
finite products/powers with explicit norm bounds; no unspecified O-constant or
existential eventual-r threshold is a finite-step certificate.

All numeric enclosures use outward192-bit dyadic rectangular complex intervals.
Their arithmetic encloses polynomial-rounding error explicitly. The analytic
Taylor tail is added separately. Complete one-flip units may use their exact
signed-permutation/diagonal structure; the full-H reference may use a general
matrix polynomial with scaling and finite repetition. Independently recheck with
a generic matrix route and a separately implemented rational interval path at
at least160 bits, retaining its actual source,precision and limitations.

For a true overlap z and approximate overlap w with normz<=1 and norm(z-w)<=eta,
the squared-infidelity uncertainty is at most eta*(2+eta). Retain rigorous
finite intervals for1-normz². Strictly separated intervals prove the actual
finite-r ranking. Its ordering is the same as1-normz for unitary evolutions;
the existing Lean conversion theorem supplies that algebraic equivalence.
No normalization of an approximate state may hide a norm/rounding defect.

Exact integer/rational identities, directed interval arithmetic, generic Lean
bounds and an end-to-end kernel certificate are separate evidence layers. The
molecular example remains software-certified unless all concrete input facts
and interval witnesses are actually imported into and checked by Lean.

## Resources and completion

Before the main panel, exercise the full path on the same model with every
coefficient scaled7/5 and a distinct preflight label, at epsilon1 and1/256,r2.
This preflight cannot change the fixed main panel. Preserve failures and all
inconclusive enclosures. Initial cumulative numerical phase ceilings(CPU
seconds/wall seconds/outputMiB): tests120/900/64;preflight240/1800/256;
main1200/7200/256;independent1200/7200/256;reporting300/1800/256.
Numerical jobs retain quarter-coreCPU,Nice19,one thread,MemoryHigh448MiB,
MemoryMax512MiB,no swap,the existing headroom thresholds,and one worker.
Proof-only6GiB permission is a separate pending decision; it does not enlarge
these numerical limits. Freeze source/inputs before execution and retain full
attempt usage without resets.

Still required: implementation and tests of interval arithmetic/input identities,
full representative preflight,complete fixed panel,independent certificate audit,
actual Lean builds/axiom/semantic correspondence checks and endpoint inventory,
Astra/Argon convergence,publication,current-headCI,separate eligible merges and
actual-merge checks. A negative physical-epsilon guard attempt is a valid result;
missing computations or an unproved bridge remain gaps.
