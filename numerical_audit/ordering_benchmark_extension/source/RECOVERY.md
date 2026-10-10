# Durable recovery of the ordering comparison extension

The WSL reboot interrupted the previous extension campaign. The research
manuscript edits and all previously committed research artifacts survived;
source/preparation files and the new extension runs stored under `/tmp` did not.
The extension is being reconstructed here, with the pinned Python 3.12 package
environment in the sibling `benchmark-venv-py312` directory.

This directory is a working area, outside the reviewed paper and QHAT snapshots.
It will be packaged as a separate, manifested research artifact when complete.
Existing frozen manifests must not be rewritten. All work remains local.

* `fixed/`: exact target fixed orders and fixed random permutation controls.
* `optimized/`: conserving and individual-Pauli objective controls, Strang arms.
* `published/`: published algorithms and explicitly described adaptations.
* `independent/`: independently constructed Pauli actions, target exponential,
  circuit replay, saved-state metrics, and work-budget checks.

Completed schedules and individual states/rows are written durably, and native
long searches checkpoint progress. Output reuse requires matching input and
source/schedule hashes; existence alone is not validation. No numerical result
from the lost temporary campaign will be reported without regenerating it.
