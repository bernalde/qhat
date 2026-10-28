# Cubic DP integration validation

2026-10-28, America/New_York. This extends the completed Julia feasibility pilot at5b9e57f5. Its original reports, certificates, numerical outcomes and resource receipts remain unchanged. This is a newly authorized cubic integration study, with its own finite accounting; it does not restart an unfinished prior phase.

## Frozen target

Use the existing CH4_r1p09 and CO_r1p13 eight-unit windows, their exact parent orders, input bytes, window positions157:164 and168:175, T=1, HF occupations, and budget30U. Do not reselect windows or train new orders. Global first/last units remain fixed, so r=15 in both domains.

The target is exactly the already-declared frozen-table objective:

`Q(pi) = || R_baseline + sum_S Delta_S(pi|S) ||^2`, `J=Q/r^4`.

All real/imaginary baseline and DELTA components denote their exact supplied binary64 dyadic values. Choose a common denominator D, scale vectors to integers, and optimize N=D^2 Q. Expand the full square into constant, per-scope linear/self, and cross factors involving at most six local labels. Keep the complete outside contribution inside the baseline; never replace the coherent square by a sum of independent norms. The original transcendental kernel and physical error remain distinct from this encoded target.

Factor tables and their aggregation must be exactly equal to the direct dyadic-vector evaluation for all tiny test permutations. Tests must include complex interference, negative cross terms, constant/zero cases, nonmonotone global labels, and a genuine six-label interaction. On both real cases, compare the DP result with the previously exhaustive40,320-order optimum, a fresh complete enumeration, and direct full Taylor scoring. A different tied optimum is allowed only after exact score and domain verification; both original domains recorded a unique optimum.

## Implementation and admission

All new numerical algorithms are Julia and use standard libraries. The existing generic bag-order DP validates its decomposition before states are allocated. Correct its conservative arithmetic-bit admission bound so integer denominators of1 add zero denominator-growth bits; retain a regression for real rational denominators. This changes admission accounting, not encoded costs or the optimizer's mathematical result.

Default limits:200,000 scalar-factor entries,1,000,000 exact-vector component entries,4,096 bits per exact component,128MiB estimated adapter/DP storage,100,000 total bag permutations and50,000 per bag. The new higher-arity adapter may use up to100,000,000 estimated DP work operations. Every numerical job also has the external512MiB cap. Allocation/work/CPU refusals remain explicit outcomes; no heuristic or partial result becomes an optimum certificate.

## Validation and timing

Preflight uses the CO model with coefficients scaled by7/5, the same fixed parent and same window positions. The scaling is declared before measuring its cubic objective and is not used to tune the main input/order/time. Exercise conversion, factor construction, checked graph/decomposition, DP, independent direct scoring, serialization and complete enumeration within the preflight pool. Inspect actual work and memory before main runs.

Main comparison uses both unchanged unscaled inputs. Retain first invocations and three warm repetitions of DP and complete enumeration, alternating method order. Each measured DP call includes its factor construction and graph/decomposition; report these separately from solve/validation. Parsing, common response-table construction, and JIT/process startup are separate. Also report an already-built-factor solve where useful, without treating its construction as free. Whole-method CPU cap60seconds; a cap/refusal is recorded and not retried with a weaker scientific question. The fast enumerator's floating screen plus exact contenders is a certified alternative for the same dyadic target, not an equal-arithmetic control.

Check archived integer/rational objective equality and returned full orders, input hashes and fixed endpoints. If the exact unique orders match the prior physical study, reuse that physical evidence with byte/parameter/order bindings; do not rerun unchanged physical errors or claim a new physical advantage. No global molecular optimum, approximation ratio, intrinsic width/hardness or scaling result follows from these two development cases.

## Resources

One root-owned quarter-core worker,512MiB/no swap, one numerical/GC/BLAS thread, existing headroom policy. Phase ceilings are cumulative across failures: tests120CPU/600wall seconds; representative preflight100/500; main validation240/1,200; independent checks100/500. Total560CPU/2,800guardedwall seconds and128MiB numerical output. Preserve unique attempt labels, original logs, receipts, checkpoints and partial results. Root alone integrates code and records local results; helpers do not run numerical commands or publish.
