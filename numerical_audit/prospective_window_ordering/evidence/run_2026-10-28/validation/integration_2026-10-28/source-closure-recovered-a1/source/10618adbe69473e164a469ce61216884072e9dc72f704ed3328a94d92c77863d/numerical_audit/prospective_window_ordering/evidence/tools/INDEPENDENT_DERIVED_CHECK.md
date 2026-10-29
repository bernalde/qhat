# Frozen independent-check input

This is a schema specification, not a populated or executed artifact. Final source head and final result paths are unknown. Root supplies actual paths/digests only after the final source fixes, required new preflight and whole main64-order freeze. Do not populate placeholders with invented outcomes.

Use `ROOT/INDEPENDENT_DERIVED_INDEX.main.a1.json`, schema `prospective-independent-derived-check-index-v1`. Every descriptor below is `{ "path": "ROOT-relative/path", "sha256": "full lowercase SHA-256", "bytes": optional_exact_integer }`. Relative paths cannot escape ROOT or follow symlinks outside it. Copy external source/dispatcher/guard files into inert evidence paths before referencing them; recorded historical absolute paths are never opened.

Required top-level fields:

- `expected_scientific_source_commit`: final full 40-character scientific-source commit, not an earlier reviewed head.
- `analysis_index`: descriptor of the final `ANALYSIS_INDEX.main.json`.
- `analyzer_summary`: descriptor of the final analyzer-generated JSON summary.
- `analyzer_source`: descriptor of the exact analyzer source whose digest appears in that summary. It must include the separately proposed primary-source inventory check and six-baseline comparison expansion; no earlier four-baseline summary can pass.
- `resource_pools`: descriptor of the unchanged original `RESOURCE_POOLS.json`.
- `original_campaign_freeze`: descriptor of the original campaign freeze containing `resource_pools_sha256`. Do not substitute a later reset pool.
- `guard_source`: descriptor of the frozen owning guard source, matching `resource_pools.guard.sha256`.
- `dispatcher_sources`: bounded list of descriptors of every root dispatcher source version named by retained job PLANs. Each PLAN's `dispatcher_sha256` must match one supplied source. These sources are hashed, never imported.
- `ledger_cutoff_UTC`: timezone-aware ISO timestamp frozen before launching this checker, after all jobs to be reconciled are terminal.
- `jobs`: complete list of all retained jobs started at or before that cutoff, not just successful jobs or main constructors. Each element has descriptors `plan`, `final`, `receipt`, `log`. Their paths must respectively be `jobs/LABEL.PLAN.json`, `jobs/LABEL.FINAL.json`, `logs/LABEL.resources.json`, `logs/LABEL.log`. The checker independently scans ROOT/jobs for omitted PLAN labels through the cutoff. A missing FINAL/receipt is a reconciliation failure, not zero work.
- Optional `input_reuse`: object with descriptors `sidecar`, `original_protocol`, `original_inventory`. Supply it when chemistry inputs were reused. The sidecar must bind the new protocol and execution commit, original protocol/inventory and chemistry-generation commit, exact original payload map, and `chemistry_rerun=false`. Preserve the original chemistry producer artifacts; do not invent a new SCF run. This validates the supplied reuse edge; a multi-edge history still needs its complete separately retained provenance chain.

The analysis index retains its normal source/input/plan/whole-order/per-case descriptors. Each case must provide inventory and PANEL, and may explicitly omit incomplete candidate/objective/primary/independent artifacts. Complete scientific success requires all four real graphs. Primary manifests must actually carry the new exact four-file `source_hashes`; the checker will not backfill old artifacts. Physical state files are hashed and their full-space byte sizes checked, but no state-vector metric is recalculated by this glue checker. Those numerical checks belong to the independently authored Julia checker already bound in the graph.

The final analyzer index's `phase_receipts` must use the actual dispatcher owner phase names for this independent contract. All such receipt descriptors must occur in `jobs`. Phase totals in the checker include every old failed and successful job through the cutoff exactly once. Analyzer `other_phase_totals` instead exclude current constructor receipts already charged in its construction totals; the checker compares that narrower raw-guard convention separately.

The check itself is launched after the cutoff, so its incomplete running PLAN must not be inserted into the frozen jobs list. The checker authenticates the old cumulative ledger, then its own cost and later jobs must be charged by root and included in subsequent final accounting. A cutoff is an explicit measurement boundary, not a campaign reset. Plans after the cutoff are not read as evidence for these totals.

Suggested future CLI (UNRUN):

```sh
python /path/to/frozen/independent_derived_check.py \
  /path/to/campaign/ROOT INDEPENDENT_DERIVED_INDEX.main.a1.json \
  /path/to/campaign/ROOT/INDEPENDENT_DERIVED_CHECK.main.a1.json
```

Root must dispatch this under its existing frozen independent/packaging pool, owned guard and single-job serialization. No unbounded retry or new pool is implied. The checker uses stdlib only and writes a new output with exclusive creation. Inputs: max32MiB per ordinary artifact, max512MiB total hashed bytes, max8192 artifacts, max2048 ledger jobs; full state files must be exactly262144 bytes. If the real final package exceeds a ceiling, record refusal and reassess the bounded specification explicitly; do not silently bypass it.

Exit0 requires complete derived agreement and an allowed unqualified under-ceiling claim. Exit2 preserves the diagnostic output for censored/incomplete inputs, source/schema/hash discrepancies, derived mismatch, missing raw accounting, or conservative accounting veto. `derived_agreement_with_accounting_veto` means the raw analyzer convention agrees but at least one call or phase crosses its conservative boundary. Its `unqualified_under_ceiling_claim_allowed=false` is authoritative for that claim; it must not be hidden in prose. Raw unknown CPU remains null alongside a separately identified reserved-allowance charge.
