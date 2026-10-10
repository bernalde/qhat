# Archive packaging and replay review

Reviewed 2026-10-09. Scope: `packaging/build_archive.py`, `verify_archive.py`, `replay_one.py`, `REPRODUCE.md`, and the binding fields consumed from independent-check records. No numerical propagation, expensive preparation, or changes to frozen scientific verifier/producer sources were performed.

## Findings and implemented resolutions

1. **Stale successful summaries did not bind the bytes being packaged.** The original builder accepted aggregate status/count fields only. A changed/deleted saved state, changed row, or changed frozen input after audit could therefore enter a new self-consistent archive. With authorization, `build_archive.py` now performs a hash-only preflight against existing evidence: exact case/row inventories; independent freeze/design/source/row/state bindings; analysis input hashes and counts; original QHAT inputs and recorded propagation-source hashes. It rejects stale evidence before writing the archive. This reuses prior scientific checks without claiming to rerun them.

2. **Output within an input tree could recursively bundle earlier archives on subsequent builds.** The builder now rejects an output path under either extension or QHAT, including resolved aliases. It also rejects symlinked packaged inputs rather than silently omitting required source paths, and verifies the temporary ZIP inventory and digests before the final rename.

3. **Independent-checker reuse is not fully bound to cached references.** `check_references.py` validates source hashes when resuming but not then-current input bytes/output-NPZ bytes; `check_campaign.py` checkpoint reuse omits a digest of the reference NPZ. Existing fresh audits remain valid. Editing these frozen checkers would invalidate their historical source hashes and require avoidable audit reruns. The parent maintainer updated the published recipe to require a newly created `fresh` output directory and explicitly delimit resume behavior. The packaging preflight validates original reference-input hashes and requires cached NPZ presence; its manifest/record explicitly says historical reference-cache output digests were never recorded. The archive manifest binds the packaged reference NPZ bytes from packaging onward, not retroactively.

## Verification evidence

- `packaging/check_preflight.py` uses tiny files rather than molecular states. All 15 probes passed under Python `-O`: valid fixture; changed/missing state; changed/missing row; changed design; changed producer/verifier source; changed freeze; stale summary counts; extra rows; duplicated audit checks; allowed outside output; and rejected output inside either input tree. Evidence: `packaging/preflight_checks.json`.
- Re-ran the seven existing archive-integrity acceptance/rejection probes under `-O`; all passed. Evidence: `packaging/review_verifier_checks.json`.
- Ran the new hash-only binding check against all completed fixed, optimized, and published-extra campaigns: 4,116 rows, including 4,004 saved completed evolutions, matched their records and original QHAT input bytes. Elapsed time was 11.76 seconds. Evidence: `packaging/preflight_completed_campaigns.json`.
- Final V2 and combined-analysis preflight were pending at review completion; the normal builder now requires them and cannot package a partial campaign. No claim is made that the final archive was already built or extracted.

## Other reviewed boundaries

- The verifier and new preflight use explicit exceptions, so their gates survive `-O`. The verifier checks canonical relative POSIX names, traversal, duplicate members, symlinks, inventory equality, sizes, and streamed member digests. It does not extract files. The external SHA-256 remains necessary to authenticate the complete archive against jointly altered members/manifest; the documentation states that distinction.
- Source selection includes the QHAT modules observed by the existing production-replay import record and preserves their package layout. `replay_one.py` explicitly resolves extension `common.py` and all loaded `qhat.*` module files against the selected roots. Fresh virtual-environment use, the fixed `qhat` directory name, and the prohibition on an unrelated editable checkout are stated. This is import-path validation, not an operating-system sandbox.
- Production replay regenerates the approximate state from the frozen schedule, uses an unmerged but equivalent symmetric palindrome, and compares to the saved state. It is accurately described as a production replay rather than an independent numerical implementation.
- Optional future hardening: validate row `T`, case/ordering/budget identity, and formula order explicitly inside `replay_one.py` instead of relying on the separately verified archive and independent campaign checks. Those bindings are already checked by the independent audit, so this is not a defect in the current measured results.
