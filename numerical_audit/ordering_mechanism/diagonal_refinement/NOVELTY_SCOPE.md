# Primary-literature scope check for the diagonal-reference ordering objective

Checked 2026-10-09. This is a focused comparison, not an exhaustive novelty
search or a priority claim. Searches combined Trotter ordering, Hartree–Fock,
variance, diagonal reference, interaction-picture response, sinc, insertion,
and perturbative error. The statements below rely on the linked primary
papers, not search summaries or third-party paper reviews.

**The individual ingredients have strong antecedents.** A potentially useful
contribution is the particular efficient combination: a fixed molecular
target, conserving mask units with individually placeable diagonal terms,
exact determinant commutator amplitudes, and a finite-time diagonal-reference
response objective used to select an order without exact dynamics. The
present search does not establish that this combination is unprecedented.

| Primary source and exact location | Antecedent and distinction |
|---|---|
| Gorin, Prosen, Seligman, Žnidarič, *Dynamics of Loschmidt Echoes and Fidelity Decay* (2006), §2.1, Eqs. (15)–(16), [preprint PDF](https://arxiv.org/pdf/quant-ph/0607050) | Integrated interaction-picture perturbation variance already controls leading fidelity loss. Our factor for one minus absolute overlap and the product-formula specialization must be derived carefully, but the variance principle is established. Do not transfer their displayed remainder without checking its assumptions. |
| Hastings, Wecker, Bauer, Troyer, *Improving Quantum Algorithms for Quantum Chemistry* (2015), §4, Eq. (6), and §5, Eqs. (7)–(8), [v2 text](https://arxiv.org/html/1403.1539v2) | HF cancellation motivates molecular interleaving; a diagonal/off-diagonal reference also motivates gap-dependent coefficient corrections. Thus HF-informed cancellations, interleaving, and diagonal-reference reasoning are not new. Their correction changes coefficients and uses step-size-dependent gaps; our ordering retains the target and uses finite-physical-time response weights. |
| Tranter et al., *Ordering of Trotterization* (2019), **arXiv v1** §V.B–C, Eq. (11), Fig. 8, [PDF](https://arxiv.org/pdf/1912.07555) | They derive insertion-induced error-operator changes and greedily choose insertion positions using the error-operator norm. Neither insertion search nor incremental commutator updates should be claimed as a new general idea. Our present objective is determinant- and time-dependent, rather than their norm objective. These locators refer to the preprint; published numbering should be checked separately. |
| Babbush et al., *Chemical Basis of Trotter–Suzuki Errors in Quantum Chemistry Simulation* (2015), Eqs. (13)–(14), Figs. 9–10 and associated classical-ansatz discussion, [PDF](https://arxiv.org/pdf/1410.8159) | They use BCH error operators and approximate electronic states to estimate molecular Trotter energy errors. Cheap HF/CI-informed error estimation is established. Our task concerns first-order, finite-time state-overlap loss, not their second-order energy-shift estimator or energy correction. |
| Zhou, Zhao, Zhang, *Trotter Scars* (2026), main Eq. (6), supplemental spectral expansion Eq. (S21), [v2 text](https://arxiv.org/html/2603.29857v2) | The sine-over-gap transport kernel is explicit. Its squared sinc form follows when a single reference eigenstate is used. Their analysis uses the Hamiltonian spectrum and optimizes initial states; our candidate substitutes the occupation-diagonal reference and optimizes order. The kernel itself is not novel, and their result does not certify our substitution. |
| Mehendale et al., *Estimating Trotter Approximation Errors to Optimize Hamiltonian Partitioning for Lower Eigenvalue Errors* (2025), §II and Appendix A.4, [v3 text](https://arxiv.org/html/2312.13282v3), [published DOI](https://doi.org/10.1039/D5DD00185D) | Perturbative, task-relevant estimators already guide partition and resource comparisons. They assess second-order ground-energy errors and leave fragment orders as supplied by the partition methods. Our fixed-target first-order overlap ordering addresses a different optimization variable and metric. |
| Maxwell et al., *Practical Estimation of Trotter Error for Hamiltonian Simulation* (2026), Theorems 1 and 3 and BCH-estimation discussion, [v1 text](https://arxiv.org/html/2606.30738v1) | They analyze asymptotic spectral BCH information and efficient estimators. Their diagonal component is taken in the **full Hamiltonian eigenbasis**; our D is the **occupation-basis diagonal Hamiltonian**. Those are different objects. Their long-time analysis does not directly justify our fixed-time diagonal-reference ranking. |
| Kronenberger, Erakovic, Reiher, *Trotter Error and Orbital Transformations in Quantum Phase Estimation* (2026), orbital-basis versus ordering comparison surrounding Figs. 4–5, [v1 text](https://arxiv.org/html/2602.18913v1), [published DOI](https://doi.org/10.1080/00268976.2026.2681062) | They explicitly disentangle basis changes from the induced change in magnitude ordering by holding a reference ordering fixed. This is relevant precedent for representation-aware comparisons. Our diagnostic restricts to sign gauges, for which each corresponding Pauli factor is exactly conjugate; changes then isolate the sorting/ownership rule. Do not claim the general basis-ordering distinction is new. |
| Liu, Wei, Foss-Feig, Gorshkov, Childs, *Ordering-Aware Theory of Trotter Error*, submitted **2 October 2026**, §III–V, [v1 text](https://arxiv.org/html/2610.03672v1) | Very recent: computable ordering-sensitive local bounds and graph-based optimization for translation-invariant one-dimensional Hamiltonians. This rules out broad first claims for theoretically motivated efficient ordering optimization. The geometric-locality, higher-order bound and de Bruijn-cycle construction differ from our molecular determinant-response objective. |

The current manuscript bibliography already contains most of these papers;
the 2 October 2026 Liu et al. paper is a new item to consider. I have not
audited every theorem in these external sources. Their role here is to delimit
claims, not to import unchecked stronger hypotheses or results.

## What still needs distinguishing experimentally

1. Use the **same refined partition and optimizer** with all response weights
   set to one. This tests whether the finite-time gap weights add value beyond
   an untransported HF BCH-variance objective.
2. Compare the same objective/search on individual Pauli units, with the loss
   of per-unit number conservation disclosed. This tests the role of the
   conserving partition rather than assuming its accuracy benefit.
3. Keep the factorial partition/search controls and validate a frozen rule on
   cases not used in development before claiming a transferable algorithm.

These controls are scientifically relevant to novelty: a good empirical
combination can be publishable, but a familiar ingredient cannot be promoted
to a new theorem merely because it was newly formalized or implemented here.
The missing theoretical bridge remains a usable bound or criterion relating
the diagonal-reference response to the full-H response in a stated regime.

Suggested provisional description: “We investigate a gauge-equivariant,
determinant-based ordering heuristic obtained by evaluating the leading BCH
response under an occupation-diagonal reference. Its finite-time weights and
local search build on established perturbative and ordering ideas; the
molecular partition, efficient mask formulation, and empirical performance
are evaluated here.” This describes the method without asserting priority.
