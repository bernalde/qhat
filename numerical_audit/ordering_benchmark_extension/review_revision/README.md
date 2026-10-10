# Review revision: symmetric controls and reproduction corrections

Addresses the extension review at QHAT `2d67b236` alongside the paper's review6
carry-over corrections. No old molecular states or frozen schedule/production
source bytes were altered. The separate original archive remains required.

- `review_supplement.zip`: complete new control inputs, designs, 504 states
  (420 new evolutions and 84 unchanged baselines), independent replay records,
  corrected checker, 24 missing historical NPZs, archive/import guards, and
  recomputed framing analysis. `archive_record.json` is its external identity.
- `artifact_manifest.json`: complete supplement member inventory.
- `source/`: readable code, summaries and instructions; large arrays and
  per-row records are supplied by the ZIP. A Git clone alone does not contain
  the complete executable data artifact. `source_manifest.json` checks this copy.
- The original `extension_validation.zip` is unchanged and has SHA256
  `d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43`.

Verify the new ZIP using the parent `verify_archive.py` and its separately
retained archive_record.json values (expected archive SHA, manifest SHA and
member count). Extract it into a fresh directory; extract the original archive
separately. Follow `strang_controls/REPRODUCE.md` and
`archive_support/README.md` in the supplement. Run from a neutral directory,
using pinned Python3.12 and one numeric thread under a 6GiB cap. The documented
independent_verify.py corrects a checker-only HF bit-convention assertion;
the correction is explicit and did not change production states or schedules.

At 10U, forward conserving Strang versus the larger forward Pauli JD order
has 9/1/4 wins/ties/losses and median error ratio 0.831. Reversing only the
conserving order gives 4/0/10 and ratio 1.34. The method retains its conservation
and specified-search-cost benefits; this control does not isolate an accuracy
gain caused by the conserving partition. Read the complete orientation and
threshold tables rather than selecting a favorable orientation retrospectively.

The corrected H2O postprocessor is separate at
`../../followup_h2o/review_revision/README.md`; it reuses the original states.
Nothing in this directory constitutes public release or Lean verification of
Python numerical code.
