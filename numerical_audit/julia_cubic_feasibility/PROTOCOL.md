# Julia DP and cubic-window feasibility pilot

Started 2026-10-27, America/New_York. User steering on 2026-10-28 prioritizes a direct comparison with the completed nine-unit enumeration cases. This protocol freezes the new comparison before measuring new solver outcomes. Prior studies and certificates remain unchanged.

## First decision: reproduce the archived nine-unit objectives

Use all three retained `exact_response_windows` domains: F2_r1p55, H2O_scale1p10, Li2_r2p95. Their nine labels, outside order, integer commutator amplitudes, supplied binary64 weight ratios encoded as integers, objective denominators, and old certificates are fixed. This is the earlier **first-order J_D**, not the cubic Strang objective J. No r^4 factor is introduced.

The new Julia factor construction expands the same coherent squared response into signed scalar factors on at most four labels, retaining the full outside response. Its checked bag-order DP must return the archived exact integer minimum. Its reconstructed full order must preserve the outside positions and attain that value under an independent amplitude evaluation. The old number of minimizers and canonical order are checked by fresh complete Julia enumeration; the DP need not count all minimizers or return the archive's particular tied optimum.

For performance, benchmark fresh Julia enumeration and Julia DP in the same worker/runtime on the same domain. Record parsing, factor construction, graph/decomposition, solve, and validation separately. Retain the whole cold-process resource receipt and each method's first invocation, then three warm repetitions, alternating method order by repetition. The first invocations share a process and JIT helpers; they are not independent cold-start comparisons. Report actual factorial bag states, factor-table entries, allocation estimates, and refusals. Historical Python timings are reference measurements, not a controlled language speed comparison. Do not infer scaling from three nine-unit cases.

Preflight on H2O with every integer weight multiplied by 7. Its exact minimum scales by 7; the set of minimizing orders stays the same. This checks the full adapter, both methods, reconstruction, serialization and replay without retuning the frozen domains. Main runs only after the preflight fits with margin.

DP allocation caps: 100,000 total bag permutations, 50,000 permutations per bag, 10,000,000 estimated factor/transition work, 128 MiB estimated state memory, and the module's coefficient-bit and factor-entry caps. Validate the decomposition and reject over-cap instances before allocating states. A refusal or a large heuristic width is not proof that every decomposition is infeasible.

## Second decision: cubic feasibility and local headroom

This remains the originally authorized pilot, distinct from the archived comparison. Use the retained CH4_r1p09 and CO_r1p13 development inputs and frozen cubic parent orders. Reference is the supplied HF occupation, T=1, ideal rotation budget 30U. These are known development cases, not new holdouts.

Structural audit: at most 150,000 examined atom scopes and 75,000 retained nonzero computed atoms. Stream destination overlaps and cap potential pair work at 1,000,000. Distinguish incomplete scans, floating support diagnostics, potential supergraphs, and exact factor evidence. Do not certify full-graph tractability from an incomplete scan or impossibility from a large heuristic width.

Select one interior eight-unit contiguous window per case by the largest number of computed nonzero internal cubic pair scopes, earliest start on ties. Selection uses operators/reference only, not physical errors. Freeze selection and parent before solving. Condition the full cubic response on the outside order and retain interference. The restricted optimum is for supplied frozen binary64 response DELTA tables. Float screening has an explicit roundoff envelope; all overlapping contenders are evaluated exactly as dyadic rationals. Report any numerical or resource refusal, and do not convert it to a zero score or certificate.

Compare with insertion within the same fixed window, and with a continuation of the existing optimizer from the same parent under equal additional CPU allowance if the representative preflight supports it. Include construction costs and distinguish restricted-window and full-order search domains. No global approximation ratio follows from a window optimum.

If improved changed orders are produced, freeze them before actual error evaluation. Physical validation uses the same T=0.5,1.5 and 10U,30U caps, full-space circuit propagation, number-leakage accounting, and projective state distance. Any sector restriction is allowed only for the exact number-conserving reference, with separate checks; circuit propagation retains the full space. Reuse unchanged-order physical results only when provenance/parameters match exactly. No physical advantage is claimed merely from a lower surrogate.

## Resources and evidence

Julia 1.10.11 already installed; standard libraries only. One worker, one Julia/BLAS/GC thread, CPUQuota25%, Nice19, MemoryHigh448M, MemoryMax512M, no swap. Start above5.5GiB available; freeze owned worker below5GiB. Record JIT/startup costs separately from warm solve costs. Julia is not presumed faster by language alone.

Finite phase ceilings, cumulative across retained failed/successful attempts:

| Phase | CPU seconds | Guarded wall seconds |
| --- | ---: | ---: |
| Runtime and tests | 180 | 900 |
| Representative preflight | 60 | 300 |
| Main structure | 160 | 800 |
| Main comparisons/windows | 160 | 800 |
| Independent verification | 180 | 900 |
| Conditional physical validation | 180 | 900 |

Total ceiling920CPU/4,600guarded wall seconds. Study output cap128MiB. Root dispatches all jobs through the existing resource guard; helpers do not execute numerical jobs. Every attempt gets a new label and retained logs/receipts; no counter resets. The full new source and input hashes must be recorded before main jobs. If essential work cannot fit, preserve evidence and report the precise incomplete decision.

The historical JSON inputs are serialized to TOML once, without numerical transformations. The new numerical implementation, exact DP, enumerator, window solver, and study entry points are Julia. Existing Python is used only as the retained baseline/independent legacy oracle or external resource supervisor. No new Python optimization backend is introduced.
