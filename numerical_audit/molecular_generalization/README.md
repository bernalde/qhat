# Prospective molecular-family generalization

This study tests full-order cubic insertion, width-eight window refinement, cheap fixed orders and a new full-order Julia annealing control on two families held out from current cubic development. HF and H2S are selected prospectively; older first-order/Strang H2S use is disclosed in [PRIOR_USE.json](PRIOR_USE.json). No universal never-studied claim is made.

Implementation and old-case preflight are in progress. Main execution is disabled until the complete protocol, calibration, source and environment are frozen. No new physical results are claimed.

The annealer is a classical optimizer written in Julia for the same cubic response J. It uses four fixed starts, unrestricted insertion proposals and geometrically cooling Metropolis acceptance. Temperature scale comes from fixed CO/CH4 development probes; no new-family physical error is used to choose settings. It supplies an empirical control, with no approximation-ratio or global-winner claim.

Current implementation checks pass 650 Julia mathematical assertions, 39 Python contract tests and the adapter admission/retention/objective checks. Old calibration has verified all native scores and the CO cross-checks; CH4 Julia cross-checks and the full old-case end-to-end preflight remain pending. The bounded completion driver reuses authenticated completed work while retaining the original failed attempt and its costs. Original producer snapshots and later consumers are recorded separately.

See [REGISTRY.json](REGISTRY.json), [PROTOCOL.json](PROTOCOL.json) and [SCHEMA.md](SCHEMA.md) for the declared cases, methods, costs, row matrix, interfaces and admission limits. Existing study producers and production defaults are unchanged.
