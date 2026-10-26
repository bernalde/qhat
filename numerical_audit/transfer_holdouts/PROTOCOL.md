# Prospective molecular, state and time transfer

Issue17 extends the completed Li2/F2/H2O panel to six fixed geometries in three
other families: LiH,BeH2,NH3. REGISTRY.json fixes every coordinate, STO-6G basis,
charge/spin, three-spatial-orbital active window and frozen core before new input
or ordering-error generation. These are bounded active-space test problems, not
an equilibrium-geometry or chemical-accuracy claim. No failed/uninformative case
will be replaced after its errors. The distinct LiH1.841 chemistry preflight is
excluded from the main panel and cannot tune the main physical errors.

## Complete target panel and information rules

For each six-mode Hamiltonian prepare HF, the sector ground and first excited
states, and (g+e)/sqrt2 and (g+i e)/sqrt2. Fix Nalpha=Nbeta=Ne/2. Exact sector
preparation is charged as classical work and is not an efficient quantum-state
preparation protocol. Label correlation using the one-particle-density-matrix
idempotency defect: above1e-10 correlated, at/below1e-12 Slater-like, intermediate
indeterminate. Do not call every eigenstate correlated. Retain energy degeneracy,
small signal, floor and loss cases; report them by family/state/time. Eigenvector
phases use the largest component with a fixed lowest-index tie. Degenerate states
must retain their actual archived eigenbasis and an explicit degeneracy label;
no claim of a basis-independent coherent phase is made in such a cluster.

All searched response objectives use training time1.11 and training cap30U.
The legacy current-placement control/initializer retains its fixedT1 rule;
that exception is recorded and also disjoint from every validation time. Frozen validation
times0.59,1.41,2.73 are all disjoint from training. Search has no access to their
physical errors. Evaluate both10U and30U using the actual affordable fused
symmetric step count, or the exact first-order count as applicable. All five
states, all three times and all controls remain:6*5*3*8*2=1440 logical circuits.

The eight arms are current-placement first order; HF-trained first-order JD;
HF-trained larger-search individual-Pauli JD symmetric; HF-trained cubic JD2
symmetric; state-aware first-order JD; that same state-aware order symmetric;
state-aware larger-search individual-Pauli JD symmetric; state-aware cubic JD2
symmetric. HF-trained schedules are frozen before non-HF evaluation. State-aware
orders are separately trained on their stated reference state at1.11, then frozen
before any validation time errors. The larger-Pauli and strong symmetric controls
are retained. This study tests per-instance frozen algorithms on new molecular
inputs, not zero-input transfer of an order between incompatible unit lists.

HF and state-aware HF computations with identical inputs/methods are an explicitly
shared construction, not extra samples. Each first-order order used by two
formula arms is also built once. Deduplicate only exact identity of model,state,
training time,objective,partition,training cap,initializer and algorithm. Report
actual unique campaign work separately from repeated standalone charges.

## Generalized objective and numerical verification

For symmetric cubic defect C_pi, use the complete signed matrix columns required
by the initial state's support. Transport each entry by
F_uv=T*sinc(T*(E_u-E_v)/2)*exp(i*T*(E_u-E_v)/2), then sum
(Q_D psi)_u=sum_v F_uv*C_pi_uv*psi_v before projecting away psi and squaring.
No incoherent average of determinant scores is an allowed substitute. HF reduction
must agree with the completed cubic objective. Independently check the signed
operator with literal nested commutators and test real/imaginary coherent phases.
The analogous generalized first-order response retains coherent interference.

Use full64-dimensional propagation without projecting leakage away. Independently
verify source tensors/active-space HF contraction, exact supplied Pauli targets,
complete-unit conservation, state normalization/sector/correlation labels,
objectives and retained search moves, physical actions/errors and any compiled
resources. Bounds precede allocation:6 modes,at most256 Pauli terms/units, at most
9 occupied-spin-sector state-support entries. Refusals are retained without
substitution. The fixed8-mode/128-term legacy helper guards are not silently
relaxed; use validated adapters where the declared new domain requires them.

Use the final reviewed #16 accounting/compilation conventions, with any deliberate
study-policy differences explicitly recorded before execution. Keep input/SCF,
exact state preparation, model/method setup, initialization, search, cleanup,
compilation, screening and independent replay costs. Classical state-preparation
CPU is not a quantum preparation circuit cost. An end-to-end claim must retain
that unimplemented quantum preparation cost as unknown. Compiled counts are ideal
prototype policy measurements, not hardware execution time.

## Execution state

This is a prospective panel/implementation protocol. All six chemistry inputs have now been generated with the frozen registry.
All 30 reference states have independently verified sector, phase and RDM labels.
Reference-state preparation and the full circuit panel remain distinct phases. Full search-control allowances, full-path preflight
budgets, training/compilation source identities and main resource ceilings must
be frozen before their respective executions. The numerical memory profile stays
512MiB,no swap,one quarter-core worker,one thread,Nice19 and existing headroom gates.
No parallel numerical/build/review clients. Final publication depends on the
verified #16 conventions and requires complete independent Astra/Argon review,
current-head CI, separate merges/actual checks and all-criterion closeout.

Initial finite input/response-unit and chemistry phase ceilings, frozen before execution (CPU
seconds/wall seconds/outputMiB): input/response unit checks120/900/64;
chemistry preflight120/900/256; six-input generation240/1800/256.
All failed attempts remain charged. The chemistry preflight is not the required
full state/search/compilation preflight and does not complete the issue's panel.
The latter receives its own finite initial ceiling before it starts; these
chemistry phases may not be relabelled or restarted to reset a spent allowance.

Reference-state constructors independently reread/validate the input. Each
non-HF constructor recomputes its sector eigensystem and includes diagnostics
and cleanup through return; those repeated standalone preparations are actual
campaign work. HF construction performs no eigensystem. Correlation labels use
the previously frozen RDM thresholds; a ground/excited gap at or below1e-9 is
labelled near-degenerate, while canonical eigenbasis clustering uses1e-11.
All such cases remain in the panel.

Reference-only preflight and reference-panel preparation/independent-check
ceilings (each cumulative CPU seconds/wall seconds/outputMiB):120/900/256
and120/900/256 respectively. Their scope is input/operator/state preparation,
independent CAR/tensor/geometry/RDM checks and reference archive extraction;
no ordering-validation error is computed. They do not replace the later full
search/evolution/compilation preflight or reset any existing phase allowance.

Prepared training adapters retain three starts(current,reverse,canonical),
full caller-through-return context/setup/cleanup charges and coherent first-order
features. Conserving first-order search has a nominal2 CPU seconds/3000 candidate
cap/5 moves per start. The stronger individual-Pauli control retains8 CPU seconds,
20*m*(m-1) candidate evaluations(capped325120),10 moves per start. Cubic training
charges a complete fresh conserving first-order warm constructor, then a nominal
2 CPU seconds/3000 candidates/5 moves of generalized cubic search. These are
search allowances, not matched-total construction claims; setup/reference and
one-circuit compilation charges remain explicit. Reference preparation is charged
once at the outer cubic boundary. Real-panel training remains pending the
practical study's reviewed,frozen interface; no holdout-error comparison has run.

## Cross-study exposure record

The registry,state/time panel,both generalized response adapters and nominal
search controls were committed at86b29d161e3dbe7bedf43a8dc0c2457d5712ff2a
before issue18's molecular errors were calculated. That separate quantitative
study uses this panel's LiH1.713 input and HF state at exact time3/5,with two
fixed first-order schedules(current placement and diagonal prefix). Its
physical/scaled error and guard outcomes were produced at
5ed79f2a70093a2e3afeb6777e9ccf1c4439edcc. None of this panel's exact validation
times is3/5;the nearby0.59 is retained as originally declared.

This overlap is disclosed rather than treating the molecular inputs as
statistically independent of that quantitative example. The issue18 outcomes
do not select or retune this panel's cases,states,phases,times,budgets,methods
or initializers. The19pre-exposure response/state/training tests and original
source snapshots remain retained. Any later correctness repair must have an
independent mathematical/verification basis and preserve all declared losses
and controls;no favorable physical error selects a replacement input or arm.
Final source review must compare the executed adapters with that frozen source
and explain every change before making the prospective-transfer claim.

The [circuit and accounting addendum](CIRCUIT_POLICY.md) fixes HF preparation
versus external-state resource scopes before any transfer validation errors.
The new adapter and independent initial-state replay are implemented; their
additional toy checks now pass in the complete [78-check suite](CHECKS.md).
The original 19 pre-exposure checks remain historical evidence for the earlier
response/state/training source. The current suite also covers initializer-trace
retention and the new circuit/runner code.

The circuit addendum now fixes every logical identity and exact-sharing key:
all 1440 logical circuits remain, with 1332 unique complete-circuit candidates
and 96 outer training constructors before any censors. Ninety exact reference
evolutions and all 5040 within-setting arm pairs are planned. The separate
LiH1.841 full-path preflight retains all states/controls at T2.73 and both caps.
These counts describe the implemented plan; full-path execution, feasibility,
independent replay and the reviewed practical dependency remain pending.

The [full execution protocol](EXECUTION_PROTOCOL.md) prospectively fixes the
complete runner, initial full-study phase ceilings, paired-review interface gate,
warmups and stage-specific feasibility rules. These new full-study pools do not
reset any completed chemistry/reference/test allowance. The implementation is
unexecuted and requires its queued checks, reviewed dependency and complete
representative preflight before main validation.

The legacy initializer now retains its existing slot scores, tie choices and
move history instead of discarding them after obtaining the order. Its algorithm,
priorities, fixed T1 sinc weights and outputs are unchanged. This adds independent
checkability before transfer errors: a separate literal-matrix checker enumerates
every placement slot with exact integer BCH columns and the pinned binary64 weight
policy. First/cubic warm constructors retain those same initializer records.
The independent generalized objectives use sparse nested commutators and complex
column transport, without importing producer features or search/placement helpers.
The additional source checks passed in the complete suite; the real transfer
panel and its full representative preflight remain unexecuted.

Current implementation verification: all 78 checks passed at
`7b99593eb95e86ef3ea40a6c8621d6f34c08823f`.
This includes the unchanged legacy algorithm's complete retained decisions and
all new producer/plan/allocation/accounting/circuit/coherent-objective checks.
No full transfer preflight, new holdout error, independent panel audit or paired
source review is implied. Subsequent status-document changes do not alter the
executed Python sources; the exact tested commit remains retained in CHECKS.md.

The full independent checker now rederives the complete logical/sharing cohort,
all retained objective and insertion witnesses, exact legacy decisions, actual
and standalone charges, full-space reference/primitive/unfused states, circuit
DAG resources, every pair and the exact payload inventory. Producer plan/score/
fusion/comparison functions are not called. The independent input auditor shares
source-bound I/O/constants; its literal CAR/Pauli/state mathematics is separate. A partial trace can be counted only
in an explicit non-simulating resource audit; it cannot receive a physical error.
These additional integration checks are pending execution; the retained78-check
result predates this checker and its resource-only extension.

Complete report/delivery implementations now retain all unique and logical rows,
every pair and all costs/losses/floors. Raw resource snapshots reconcile every
completed attempt with its original pool and 512MiB numerical profile. Reports
require the complete independent audit and bind all scientific bytes. Packaging
uses deterministic archives, strict safe extraction, exact file inventories,
report regeneration and streaming Git-blob checks, including ignored binaries.
Their integration checks and full representative-path exercise remain pending;
none of this code authorizes bypassing the reviewed practical dependency.
