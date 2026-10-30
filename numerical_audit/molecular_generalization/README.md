# Prospective molecular-family generalization

This study tests full-order cubic insertion, width-eight window refinement, cheap fixed orders and a new full-order Julia annealing control on two families held out from current cubic development. HF and H2S are selected prospectively; older first-order/Strang H2S use is disclosed in [PRIOR_USE.json](PRIOR_USE.json). No universal never-studied claim is made.

Implementation and old-case preflight are in progress. Main execution is disabled until the complete protocol, calibration, source and environment are frozen. No new physical results are claimed.

The annealer is a classical optimizer written in Julia for the same cubic response J. It uses four fixed starts, unrestricted insertion proposals and geometrically cooling Metropolis acceptance. Temperature scale comes from fixed CO/CH4 development probes; no new-family physical error is used to choose settings. It supplies an empirical control, with no approximation-ratio or global-winner claim.

Current implementation checks pass 650 Julia mathematical assertions and 36 Python contract tests. Calibration module loading also passes without executing calibration. Broader adapter admission tests and the full old-case preflight remain pending. Original producer snapshots and later consumers are recorded separately; all failed attempts and measured costs are retained.

See [REGISTRY.json](REGISTRY.json), [PROTOCOL.json](PROTOCOL.json) and [SCHEMA.md](SCHEMA.md) for the declared cases, methods, costs, row matrix, interfaces and admission limits. Existing study producers and production defaults are unchanged.
