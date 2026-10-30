# Molecular family generalization — report template

**Status: prospective draft. Main results do not exist. Every bracketed entry below is a placeholder.**

This study compares the frozen full-order control, width-eight window refinement and full-order annealing on the declared HF and H2S cases, with canonical and magnitude orders as context. The scope is the retained implementations, Hartree–Fock reference, two CPU allowances and two paired timing/seed repetitions. It does not establish a universal ordering method or an approximation guarantee. Older H2S first-order/Strang use is disclosed in `PRIOR_USE.json`.

## Evidence status and frozen scope

| Item | Recorded value / evidence |
|---|---|
| Analysis result status and `complete` | [pending] |
| Analysis JSON path and SHA-256 | [pending] |
| Protocol, registry and source-freeze hashes | [pending] |
| Numeric old-only annealing calibration and accepted preflight | [pending] |
| Cases, geometries, CAS(8e,6o), input hashes | [pending] |
| Inventory | [pending: main target 36 constructors, 144 logical physical rows, four exact references] |
| Objective, independent replay and resource gates | [pending] |

Both main cases use 12 modes/eight active electrons: circuits propagate all 4096 amplitudes; the checked exact-H reference uses the complete 495-state fixed-N sector. The reference is not restricted to the smaller fixed-spin block. No replacement or protocol amendment is justified by an unfavorable or floor-limited outcome.

## Window minus full — predeclared primary outcome

At B64, T = 1.35 and rotation cap 8U, use the raw within-repetition difference `error_window − error_full`.

| Case | Repetition 1 difference | Repetition 2 difference | Case median | Eligibility / gates |
|---|---|---|---|---|
| HF | [pending] | [pending] | [pending] | [pending] |
| H2S | [pending] | [pending] | [pending] | [pending] |

**Recorded decision: [pending].** State the exact analysis decision, retaining incomplete or mixed/unresolved outcomes. Candidate advantage requires both case medians below −1e−12, no resolved opposite repetition, and every required gate. The full-control rule is symmetric. Do not infer a confidence interval from two molecules.

## Annealing minus full — predeclared primary outcome

At the same primary condition, use `error_annealing − error_full`.

| Case | Repetition 1 difference | Repetition 2 difference | Case median | Eligibility / gates |
|---|---|---|---|---|
| HF | [pending] | [pending] | [pending] | [pending] |
| H2S | [pending] | [pending] | [pending] | [pending] |

**Recorded decision: [pending].** Apply the same frozen ambiguity and completeness rule independently. Do not select whichever named question, molecule, repetition or secondary setting looks most favorable.

## Quality and complete measured construction cost

[Insert one standalone SVG/PNG per case at T = 1.35 and 8U; figures pending.]

Plot every parent, canonical, magnitude, full, window and annealing record, including both repetitions and optimized B32/B64 points. Optimized x coordinates are the measured common parent prefix plus continuation charged CPU; parent and cheap points use their complete measured own CPU. The shared parent is counted once in actual study spending, while each optimized comparison includes its attributed prefix. The label B0 is an identity for parent/cheap controls, not a zero-cost claim.

Use raw projective error on y. A symlog scale with a linear region below 1e−12 retains exact zeros. Preserve floor and ineligible markers and list missing coordinates without fabricating a position. Do not clip or replace raw values, join independent budget runs into a claimed search trajectory, rank methods retrospectively, or add confidence bands. The figure remains diagnostic when the analysis is incomplete.

| Plot artifact / companion | Hash / status |
|---|---|
| Per-case SVG and PNG | [pending] |
| Unclipped `PLOT_DATA.json` and `PLOT_DATA.csv` | [pending] |
| Input result hash and `PLOT_PROVENANCE.json` | [pending] |
| Later plotting source hash, Python and Matplotlib versions | [pending] |

The private `plot_draft.py` is an unexecuted draft. When real analysis exists, run it in a fresh directory under the declared plotting/audit resource allowance. Review the rendered figures before publication. Its provenance must preserve the original producer map, analysis-consumer map and later plotting-consumer identity separately. A plot cannot complete a missing scientific gate.

## Cheap controls and the complete grid

Report canonical and magnitude errors, eligibility and complete measured costs separately for every case/repetition; include parent quality as context. Preserve identical orders and states as repeated logical observations without treating them as extra molecular evidence.

| Case | Rep | Arm | CPU identity | T | Cap/U | Raw error | Floor | Eligible | Constructor/state/check joins |
|---|---|---|---|---|---|---|---|---|---|
| [all 144 logical rows; no selection] | | | | | | | | | |

| Named question | Budget | T | Cap/U | Available pairs | Raw differences / floor flags | Gate status |
|---|---|---|---|---|---|---|
| [all 64 declared named pairs] | | | | | | |

Secondary settings are context. Record their raw losses, ties, gains and missing observations; they do not replace either primary question.

## Refusals, capacity stops and construction accounting

| Case/rep/arm/B | Attempt / receipt | Status / refusal | Scored fallback | Starts/proposals/native check | Eligibility effect |
|---|---|---|---|---|---|
| [every failed or capacity-limited attempt] | | | | | |

Retain unpromoted raw Julia proposed-best witnesses, native refusals, floor-censored objective comparisons, startup/serialization/teardown costs and interrupted work. A valid scored parent fallback is not a successful method execution. Record missing terminal CPU as unknown, never zero.

| Cost boundary | Raw CPU | Charged CPU | Wall | Unknown / failed attempts |
|---|---|---|---|---|
| Parent actual once per case/repetition | [pending] | [pending] | [pending] | [pending] |
| Optimized continuations, actual spending | [pending] | [pending] | [pending] | [pending] |
| Cheap controls, complete own spending | [pending] | [pending] | [pending] | [pending] |
| Actual study constructor total | [pending] | [pending] | [pending] | [pending] |
| Comparison attribution, including six prefixes per case/repetition | [pending] | [pending] | [pending] | [pending] |
| Calibration, input generation, physical/objective/replay, analysis and packaging | [pending] | [pending] | [pending] | [pending] |

Physical conditions never multiply constructor costs. Report cumulative phase limits, all retained attempts, measured memory/resource profile and any old-only amended preflight budget alongside observed costs.

## Provenance, replay and completion record

| Gate / artifact | Exact source or result identity | Outcome / limitation |
|---|---|---|
| Immutable producer commit, protocol and environment | [pending] | [pending] |
| Exact model JSON↔TOML, complete roster and parent/start graph | [pending] | [pending] |
| Raw constructor/worker/PLAN/receipt bindings | [pending] | [pending] |
| Every selected and raw proposed-best objective witness | [pending] | [pending] |
| Every logical physical row, distinct circuit and exact reference | [pending] | [pending] |
| Independent Python consumer and Julia replay sources | [pending] | [pending] |
| Package inventory, extraction and relocated replay | [pending] | [pending] |
| Later plotting consumer and figure/data hashes | [pending] | [pending] |
| Formal review passes and findings | [pending; do not infer from implementation helpers] | [pending] |
| PR URL, reviewed head and final CI | [pending] | [pending] |
| Merge commit and merge verification | [pending] | [pending] |
| Postmerge installed-path / exact-head checks | [pending] | [pending] |

Close with the measured scope and unresolved limitations. Main evidence remains pending until the frozen pipeline and all required gates finish. Missing evidence, capacity failures and adverse outcomes remain part of the report.
