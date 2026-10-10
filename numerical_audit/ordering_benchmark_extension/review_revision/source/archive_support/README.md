# Archive-support supplement for extension review X-08–X-10

This supplement leaves `extension_validation.zip` and its frozen source and
manifests unchanged. It supplies the 24 development NPZ files omitted from that
archive, a portable 30-state comparison, and stricter external verification and
namespace guards. No new molecular propagation is claimed by these checks.

## Verify the original archive

Use a separately retained record for these identity values. The verifier now
requires `--expected-sha256`, or an explicit `--integrity-only` acknowledgement.
The latter detects internal inconsistency but cannot detect joint replacement
of a file and its manifest. `--expected-members` counts data files **excluding**
`artifact_manifest.json`; the original archive has 9,247 data files and 9,248 ZIP
entries including the manifest.

```sh
python -O archive_support/verify_archive.py /path/to/extension_validation.zip \
  --expected-sha256 d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43 \
  --expected-manifest-sha256 0f42cf3dcf15982ca1181bf0d8674d4cd7b15bc0ba7352b35ad01f7481f9bb13 \
  --expected-members 9247
```

Extract the original archive separately, preserving the sibling `extension/`
and `qhat/` directory names. The supplemental `historical/` directory is not
copied over any original file. Do not overwrite frozen evidence with new checks.

## Reconstruct the 30 historical-state comparisons

Use the pinned numerical environment described in the original reproduction
instructions. This checker requires only NumPy, uses one BLAS thread, applies
a 1 GiB address-space limit, and performs no state evolution.

```sh
python -O archive_support/check_historical_states.py \
  --archive-root /path/to/extracted-original-archive \
  --output /path/to/new-checks/historical_check.json
```

`historical_states_manifest.json` supplies every source path, NPZ key, file hash,
and corresponding extension output hash. It binds to the original archive's
`extension/optimized/archived_state_checks.json` rather than replacing it.
The 24 missing files total 170,487 bytes and originate in QHAT revision
`2d67b236e38834231552af3100801130e0cb5984`. The six remaining comparisons use the
already archived `qhat/numerical_audit/ordering_mechanism/heldout/heldout_validation.zip`.
The saved check verifies all 30 arrays bitwise, including dtype and shape;
maximum difference is zero. This compares saved outputs, not independent
regeneration of the historical calculations.

## Prevent namespace-package mixing

The original archived `qhat/` is a Python namespace package. A live checkout in
`PYTHONPATH`, an editable installation, or an unsafe working directory can merge
its namespace with the archive even when some modules resolve correctly. Use a
clean environment and a **neutral working directory**, such as `/tmp`, with no
live QHAT checkout on the import path. The wrapper below refuses any namespace
with more than one `qhat.__path__` entry and checks containment of every loaded
QHAT module before and after replay. It delegates to the unchanged frozen
`extension/packaging/replay_one.py` only after its initial guard succeeds.

```sh
cd /tmp
python -O /path/to/archive_support/replay_guarded.py \
  --archive-root /path/to/extracted-original-archive \
  --check-imports-only --output /path/to/new-checks/import_paths.json
python -O /path/to/archive_support/replay_guarded.py \
  --archive-root /path/to/extracted-original-archive \
  --campaign fixed --case CH2O_s1.00 --ordering jw_magnitude --budget 2 \
  --output /path/to/new-checks/production_first_order.json
```

The import-only check does not simulate. A production replay retains the frozen
runner's 6 GiB address-space cap. Run these sequentially; do not launch all cases
in parallel. Choose fresh output paths; the final output is published only after
the post-replay namespace guard passes. A failure can leave a `.guard-pending`
file for diagnosis, which must not be read as a successful final check.

## Exact original replay coverage

The original campaign independently checked saved-state metrics on all 4,186
states and independently propagated 3,258 selected circuits. For each of the
four 16-qubit F2/NH3 coefficient cases, 67 of 299 states were independently
propagated and **232 were not**. At budget `100U`, 92 of the 96 states per such
case were outside that independent replay subset. Those headline points rely
on the production propagator plus independent reference/metric/structural
checks; they do not have a blanket independent-propagation claim. The external
review's 38 additional replays are separate evidence and are not silently added
to the original 3,258 count. Subsequent review experiments have their own records.

## Adversarial checks

```sh
python -O archive_support/test_archive_support.py \
  --output /path/to/new-checks/support_tests.json
```

Eighteen saved tests pass, including rejection of missing identity mode, changed
members, joint member/manifest alteration against the original digest, wrong
manifest/count, missing/unlisted/duplicate/unsafe/symlink members, mixed QHAT
namespaces, escaped module files, and an escaped late import before final-output publication.
A mock runner exercises delegation without propagating a state. Every subprocess is invoked under `-O`;
checks use explicit exceptions, never `assert`. The test harness has a 1 GiB
address-space limit and each child a 30-second timeout. The positive namespace
check against the extracted original archive is also saved. These checks do not
execute numerical propagation or rebuild Lean.
