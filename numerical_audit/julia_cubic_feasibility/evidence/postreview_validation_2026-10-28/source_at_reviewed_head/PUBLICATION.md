# Published study bundle

This branch publishes the first-order Julia DP pilot and its cubic Strang integration validation. The numerical sources, tests, input coefficients, result certificates, original review notes, resource logs and receipts are unchanged from the validated local studies. The full suite passed 1,815 assertions; both cubic windows reproduce their complete exhaustive optima. No new numerical experiment was run during publication.

The prior pilot is local commit `5b9e57f5873ec6b6948b70671a5b00d0a6cffb34`, retained as an ancestor of this publication. Cubic measurements were recorded at local commit `b82337d1eaf51cd9556f67e284a44304d721d8f1`; its measured numerical source hashes are checked in the retained independent results. That local commit also contained a private task-coordination record and is not the public branch head.

Publication omits only `evidence/cubic_dp_validation/provenance/STUDY_STATUS.json`, replacing its scientific status with `SCIENTIFIC_STATUS.json`. Local handoffs and agent/task coordination stay local. Original receipts keep their historical paths and bytes. `REPRODUCE.md` uses portable executable examples and makes the optional legacy oracle's full QHAT checkout explicit.

The original cubic-validation bundle inventory is retained as `evidence/cubic_dp_validation/provenance/measurement_bundle_MANIFEST.sha256`; it describes the measurement snapshot, including the omitted coordination record, and is not an inventory of this public tree. The root `MANIFEST.sha256` describes the complete current public bundle. The older pilot manifest and changed older source files remain under the same provenance directory.

Start with `CUBIC_DP_REPORT.md` for the current result and `REPORT.md` for the earlier first-order comparison. Cubic DP is optional: on these dense eight-label graphs it is about 31–34 times slower than certified enumeration. The studies establish neither a global molecular optimum nor a useful molecular approximation ratio, and do not replace the normal optimizer.
