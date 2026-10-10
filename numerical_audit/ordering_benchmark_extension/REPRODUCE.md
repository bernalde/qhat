# Reproduce the expanded ordering comparisons

The archive contains `extension/` (new source, frozen schedules, rows, final
states and checks) and `qhat/` (the exact source and input subset needed by the
loaders). Earlier frozen artifacts remain separate and unchanged. The external
archive record supplies its SHA-256, byte count and inventory count. Run the
adjacent `verify_archive.py` with Python `-O` and the separately retained hash
using `--expected-sha256`. Optional `--expected-manifest-sha256` and
`--expected-members` also bind the inventory to the external record. The member
count excludes `artifact_manifest.json`. A deliberate `--integrity-only` mode
checks internal consistency without establishing release identity; omitting both
identity choices is an error. This is byte-integrity checking, not independent
simulation. The original ZIP and its frozen verifier remain unchanged; the
adjacent verifier incorporates the review's stricter interface.

```sh
python -O verify_archive.py extension_validation.zip \
  --expected-sha256 d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43 \
  --expected-manifest-sha256 0f42cf3dcf15982ca1181bf0d8674d4cd7b15bc0ba7352b35ad01f7481f9bb13 \
  --expected-members 9247
```

## Environment and paths

Extract to a new directory. Create a Python 3.12 environment and install
`qhat/numerical_audit/ordering_mechanism/environment-requirements.txt`. The
campaign used Python 3.12.12, NumPy 2.5.3, SciPy 1.18.1, Numba 0.68.0 and
OpenFermion 1.8.1. Do not install an unrelated live QHAT checkout as an editable
package. The archived `qhat/` has no `__init__.py` and is a namespace package:
a live checkout can merge into `qhat.__path__`. Use a neutral extraction or
working directory and a clean import path. The review supplement's
`archive_support/replay_guarded.py` requires exactly one archived namespace path
and checks module containment before and after the frozen replay; a merely
correct `common.py` path is not sufficient. Set `QHAT_ROOT` to the extracted `qhat` directory, `PYTHONPATH` to its
parent (retain the package directory name `qhat`), and writable `NUMBA_CACHE_DIR`/`MPLCONFIGDIR` locations. Numerical
threads are set to one by the scripts. The major runners impose a 6 GiB
address-space limit; sequential replay avoids multiplying that limit.

Below, run from the extraction directory, using the environment's Python.
Choose fresh output directories; do not overwrite the frozen evidence. The
following guard must succeed: it refuses an existing output root, so cached
independent checks cannot be mistaken for a new audit of changed inputs.

```sh
python -c "from pathlib import Path; Path('fresh').mkdir()"
python -O extension/independent/check_references.py --qhat-root qhat --output fresh/references
python -O extension/independent/check_partitions.py --qhat-root qhat --output fresh/partitions.json
python -O extension/independent/objective_check.py --campaign extension/optimized --output fresh/objectives
python -O extension/independent/check_campaign.py --campaign extension/fixed --references fresh/references --output fresh/fixed
python -O extension/independent/check_campaign.py --campaign extension/optimized --references fresh/references --output fresh/optimized
python -O extension/independent/check_campaign.py --campaign extension/published/extra --references fresh/references --output fresh/extra
python -O extension/independent/check_campaign.py --campaign extension/published/v2 --references fresh/references --output fresh/v2
```

These checks reconstruct target exponentials from occupation-basis Pauli
matrices; check every saved state's normalized errors, leakage, hashes and
rotation budgets; and independently replay every small-closure circuit plus
the predeclared large-case subset described in `check_campaign.py`. They do
not claim full independent propagation of every large-case row. Specifically,
each 16-qubit F2/NH3 coefficient case has 299 saved states, of which 67 were
independently propagated and 232 were not; 92 of its 96 states at `100U` were
outside that subset. The reviewer's 38 additional replays are separate evidence,
not part of the archived 3,258 count. No QHAT code
is imported by the independent scripts. Independent-checker resume validates
its recorded source/design/row/state fields but does not bind every cached
reference input or NPZ: use the new-directory guard above for a fresh audit.
Existing checkpoint evidence applies to the immutable verified archive only. The objective checker reconstructs
cumulative operator commutators instead of importing the optimizer pair cache.

For a fresh production replay, prefer the supplemental namespace guard below.
Run from `/tmp` or another neutral directory, using absolute paths and fresh
outputs. `--check-imports-only` checks isolation without starting a simulation.
The wrapper delegates to the original frozen runner and does not change it.

```sh
python -O /path/to/archive_support/replay_guarded.py --archive-root /path/to/extracted-original-archive --check-imports-only --output /path/to/fresh/import_paths.json
python -O /path/to/archive_support/replay_guarded.py --archive-root /path/to/extracted-original-archive --campaign fixed --case CH2O_s1.00 --ordering jw_magnitude --budget 2 --output /path/to/fresh/production_first_order.json
```

The original production replay commands below check imported module files but
lack the stricter namespace-cardinality gate. Use them only after verifying that
`qhat.__path__` contains exactly the extracted `qhat` directory:

```sh
python -O extension/packaging/replay_one.py --extension extension --qhat-root qhat --campaign fixed --case CH2O_s1.00 --ordering jw_magnitude --budget 2 --output fresh/production_first_order.json
python -O extension/packaging/replay_one.py --extension extension --qhat-root qhat --campaign optimized --case CH2O_s1.00 --ordering strang_class_refined_response --budget 2 --output fresh/production_strang.json
```

The symmetric replay uses an unmerged palindrome, equivalent to the archived
merged commuting-boundary circuit; the frozen row records the counted work.
These commands reuse the frozen *schedule*, not the saved approximate state.

Regenerate tables and publication figures without propagation:

```sh
python extension/analysis/summarize.py --root extension --output fresh/analysis
python extension/analysis/timings.py --root extension --output fresh/analysis
python extension/analysis/render.py --analysis fresh/analysis --output fresh/figures
```

## Producers and longer reconstruction

`fixed/run.py`, `optimized/run.py`, and `published/run.py` expose their actions
through `--help`. Every campaign freezes all fourteen schedules before its
first new error. Production checkpoint runners require their recorded source/input hash checks
before resume; this is distinct from the independent-checker cache limitation
described above. A
full regeneration should copy the source/protocol files into a new empty
campaign layout and point `QHAT_ROOT` at these archived inputs. It should not
silently remove frozen designs or change an existing manifest.

For optimized preparation use `optimized/checked_prepare.py`, which binds
checkpoints to input/source hashes; use `optimized/execute.py` for guarded
freeze and propagation. The sidecars document the post-crash checkpoint
adoption and the direct floating-score guard diagnosis. The original input
candidate's eight development orders are checked before propagation.

The expensive native incremental-V2 insertion has source in
`published/v2_exact.cpp`; compile with `g++ -O3 -std=c++20 -mpopcnt` and Boost
headers on x86-64 with POPCNT support. The compiled binary is not archived.
The exact compiler command/version is recorded in `build_record.json`.
Each incoming term is fsynced to a checkpoint; source and input hashes bind
resume. Exact positive uniform coefficient scaling permits reuse of the
complete NH3 order for its coefficient-multiplied case. This does not permit
reusing the evolved state or its physical error.

## Interpretation

Read `analysis/RESULTS.md`, the full `summary.json`, the per-campaign protocols,
and the manuscript. The fourteen cases are correlated molecular fixtures,
not fourteen independent molecular samples. Random controls are fixed sampled
permutations, not an accuracy-selected best order. First-order and Strang
formulas have a common rotation *budget*, with actual emitted work recorded.
Classical evaluation-count matching is not wall-time matching. Tomesh beam32
and the first-order use of Tranter's second-order insertion are explicit
adaptations. Raw bytes and numerical checks do not constitute Lean proofs.


## Historical optimized-state support added after review

The original 30-state bitwise comparison depends on 24 development NPZ files
that were omitted from the ZIP. The review supplement's `archive_support/`
contains those unchanged 170,487 bytes, their original source paths and hashes,
and `check_historical_states.py`. The other six historical vectors already
reside in the original nested heldout ZIP. Run the portable checker against the
separately extracted original archive to verify all 30 comparisons; no live QHAT
checkout or imports are needed. This supplement preserves the original archive
and manifests rather than silently changing their meaning. Its README records
the checks and the distinction between saved-vector comparison and new evolution.
