# Reproducibility follow-ups for future campaigns

The current archive preserves the historical records exactly. The analyzer and independent derived checker were replayed successfully from the delivered archive. Historical job-dispatch commands retain machine-specific paths and are provenance records; their portability is narrower than that demonstrated consumer replay.

The historical F/a4 terminal preflight decision retains its earlier partial-record timestamp and continuation text; that timestamp is not a final-write attestation. Some small memory peaks are absent from the early display-only parser output. F/a4 still fails its 400 MiB margin at 402.8 MiB; admitted G/a5 uses the later unit-aware parser. Some passing decision records omit a record timestamp entirely. Retained job receipts provide their own execution times without supplying missing declaration times.

Frozen script headers such as `DRAFT/UNRUN` describe their pre-execution authoring state. Execution claims depend on the hash-bound reports and receipts, not those unchanged headers. Historical API `scope` strings also require the actual `test_files`, counts and guard outcomes; the report links the full F cohort and narrower G checks separately.

Apply the following changes when preparing the next campaign template, before its first resource freeze. These follow-ups do not rewrite this study's frozen artifacts or retroactively authenticate missing information.

1. Record every resource extension with an explicit timezone-aware `declared_UTC`, `CPU_seconds_added`, `wall_seconds_added`, prior and new caps, and cumulative charged usage. Require the increments to agree with the cap differences before dependent work starts. The historical integration `RESOURCE_EXTENSION.a4.json` records prior/new caps and prior charged usage, but does not retain an explicit declaration timestamp or separate increment fields. Its exact declaration time is unavailable in that record.
2. Bind the next campaign's guard through a path relative to its declared root, or through an explicitly documented fallback to the archived guard. Always verify the frozen guard SHA-256 before execution. Demonstrate that the dispatcher can locate the correct guard when the original machine-specific path is absent, and that altered guard bytes are refused. The historical dispatcher resolves the absolute path in the frozen `RESOURCE_POOLS.json`; copying the archive alone does not make that dispatcher portable.

3. Give every future preflight decision an explicit record timestamp. When extending a partial decision, retain distinct `partial_recorded_UTC` and `final_recorded_UTC` fields and update the continuation text to its actual state. Do not reuse a partial timestamp as a terminal-write attestation. Continue recording memory peaks with explicit units and unknown values.

Keep these obligations distinct from current numerical claims, copied-byte integrity, hosted CI, and formal human review. A future campaign author should implement and validate the new template before freezing new results; the present report does not claim those future template changes have been implemented.
