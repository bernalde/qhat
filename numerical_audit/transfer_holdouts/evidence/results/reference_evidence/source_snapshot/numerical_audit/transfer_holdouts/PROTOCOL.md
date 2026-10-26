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

All learned orders use training time1.11 and training cap30U. Frozen validation
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
