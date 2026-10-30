# Prospective molecular-family generalization

This study tests full-order cubic insertion, width-eight window refinement, cheap fixed orders and a new full-order Julia annealing control on two families held out from current cubic development. HF and H2S are selected prospectively; older first-order/Strang H2S use is disclosed in [PRIOR_USE.json](PRIOR_USE.json). No universal never-studied claim is made.

Implementation and old-case preflight are in progress. Main execution is disabled until the complete protocol, calibration, source and environment are frozen. No new physical results are claimed.

The annealer is a classical optimizer written in Julia for the same cubic response J. It uses four fixed starts, unrestricted insertion proposals and geometrically cooling Metropolis acceptance. Temperature scale comes from fixed CO/CH4 development probes; no new-family physical error is used to choose settings. It supplies an empirical control, with no approximation-ratio or global-winner claim.

Current implementation checks pass 650 Julia mathematical assertions, 56 Python contract tests and the adapter admission/retention/objective checks. Old calibration is complete: all 66 native scores passed the Julia cross-checks, using 64 fixed insertion probes. The pooled median of 38 positive normalized deltas sets temperatures 0.04155802571112692 and 0.002013694406127226. The full old-case end-to-end preflight remains pending.

[Calibration evidence](evidence/calibration_2026-10-30/README.md) preserves the original producer, both incomplete attempts, the missing-row continuation and all costs: 661.83 CPU seconds and 4,832.43 wall seconds in total. The final reference check used a declared 768 MiB memory limit; constructor comparisons retain their 512 MiB limit. No new-family physical error entered the calibration.

See [REGISTRY.json](REGISTRY.json), [PROTOCOL.json](PROTOCOL.json) and [SCHEMA.md](SCHEMA.md) for the declared cases, methods, costs, row matrix, interfaces and admission limits. Existing study producers and production defaults are unchanged.
