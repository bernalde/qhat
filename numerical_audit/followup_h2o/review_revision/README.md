# Review-6 reporting correction (no new simulation)

These reporter and guard revisions address R6-18 and R6-20. The original
`../summarize_h2o.py`, `../check_guards.py`, numerical states, designs and
source snapshots remain byte-for-byte unchanged as historical evidence.
Use the scripts in this directory for a new report. Output paths must be
explicit and fresh; they do not default to tracked output files.

The reporter calculates the uniquely-best count, median ratios, policy
agreement, extra-word counts and magnitudes, and spin diagnostics from the
supplied rows/designs. It does not print the expected findings as constants.
The shaded 4–10-qubit HGBS-5 region marks the near-floor regime. This is a
presentation update to the same measured points.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /path/to/python summarize_h2o.py \
  --sweep ../results --low-energy ../low_energy \
  --historical-csv ../presentation/archived_h2o_rows.csv \
  --output /fresh/h2o-review-report
PYTHONPATH=/parent/of/qhat /path/to/python check_guards.py \
  --qhat-root /parent/of/qhat --output /fresh/h2o-guard-checks.json
```

Run every numerical command under a 6 GiB memory cap and one numeric thread.
The guard checks import rejection and protocol/state corruption only; they do
not run a new molecular evolution. The reporting revision does not broaden
the independent-propagation coverage of the original experiment.

For reproducible tensors, use the frozen orbital arrays and orbital signs
under `numerical_audit/h2o_recovery/data`. A fresh HGBS-5 SCF solve can change
non-active orbital signs with thread count and shift active coefficients by
up to about 1.1e-10 in the reviewed environment. Regeneration from saved
orbitals and rerunning the SCF procedure are different reproducibility claims.

Integrity starts with the paper's committed
`theory/followup_evidence/provenance.json`, which hashes QHAT's original
`numerical_audit/followup_h2o/artifact_manifest.json`; the latter hashes the
saved evidence. Verify that chain before checking state metrics. A checker
alone cannot reject a coordinated change to states, hashes and metrics.
These historical manifests are not rewritten to describe this reporter:
this directory and its new outputs have their own revision evidence.
