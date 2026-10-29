# Evidence archive scope and historical limitations

The archive preserves the final scientific/product source, input payloads,
orders, complete physical and independent-check outputs, failed attempts, and
the original source declarations. Its copy manifest authenticates copied bytes;
it does not certify that every historical declaration has a byte provider.

The final 43-file scientific snapshot at G and 19-file product snapshot are
unchanged at H. H supplies the separately authenticated analysis and derived
checker. Their complete main audit covers 64 construction calls, 352 physical
rows, 384 pairs, 192 repetition medians, 12 summaries and 32 Pareto groups.

Three source versions declared by the earliest integration freezes could not be
recovered. Their original declarations remain, and the source-version binding
report deliberately retains `source_claim_complete: false`:

| Earliest declaration | Logical file | Missing SHA-256 |
| --- | --- | --- |
| PRODUCT_SOURCE_FREEZE.a1.json | .github/workflows/window-ordering.yml | 1fe9e733fcc8b9f4f3a2d48886fb79e75bb40ef4f13565dd09eff0748386455a |
| PROSPECTIVE_DRIVERS_SOURCE_FREEZE.a1.json | numerical_audit/prospective_window_ordering/panel.py | df31750f9b47ab5067f728eb00a7bcd657f917c6006b3ce98f9f086008e344e8 |
| PROSPECTIVE_DRIVERS_SOURCE_FREEZE.a1.json | numerical_audit/prospective_window_ordering/test/test_panel.py | 3b399f27b5fa71e202712098eca66f26cd17570657228f0e165228ce2fea167e |

Consequently those early integration attempts cannot be fully reproduced from
their exact declared source versions. Later source snapshots and checks do not
retroactively authenticate these missing versions. Seven other missing versions
were recovered with exact matching hashes; their retained bytes and concise
recovery record are provided separately.

Historical helpers also have mixed provenance: 22 versions have original retained
hash evidence and 16 have retrospective exact-file readbacks plus observed
filesystem continuity, as qualified in HISTORICAL_HELPER_READBACK.a1.json. The
latter are not original runtime-byte attestations.

The curation plan's earlier aspiration of complete historical source closure was
not achieved. This document records that limitation without rewriting the plan,
freezes, failed closure audit, measured outputs, or acceptance guards. Final
artifact review must assess this disclosed limitation; a successful copy or
relocated numerical audit does not remove it.

The main ledger cutoff is 2026-10-29T03:15:13.811928+00:00. Later indexing,
plotting, copying and audit costs are separate closeout records. Archive
selection has a later explicit cutoff and excludes its own running job. Files
generated after that cutoff must be accounted separately, without implying that
an earlier immutable manifest included future work.
