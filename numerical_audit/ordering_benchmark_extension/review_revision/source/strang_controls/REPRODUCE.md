# Reproduction

Use Python 3.12.12 and the exact versions in environment-requirements.txt,
which is copied from the original extension archive. All filesystem paths are
arguments; the runner does not require a live QHAT checkout.

Extract the original extension_validation.zip to a fresh directory containing
qhat/, extension/, and artifact_manifest.json. Verify its SHA256 is
`d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43` using
verify_archive.py. Preserve the qhat directory name. Before any QHAT import,
the runner binds qhat.__path__ exclusively to that directory. Every imported
QHAT module is checked against the original archive manifest.

For a new independent run, create a NEW empty directory and copy run.py,
independent_core.py, PROTOCOL.md, independent_verify.py, analyze.py, and
CHECKER_CORRECTION.md there. With the pinned Python executable, run these
commands sequentially:

```sh
python run.py freeze --archive-root /absolute/path/to/extraction
python run.py run --archive-root /absolute/path/to/extraction
python independent_verify.py
python analyze.py
```

The first command freezes all 14 plans, inputs and relevant source digests
before propagation. Production and independent verification resume after
interruption and reject changed saved states or designs. Both support --cases
for diagnosing one fixture; final analysis requires all 504 rows.

Use independent_verify.py for the corrected physical-bit HF sanity check.
The original frozen run.py verify is retained for historical source identity
and should not be invoked; CHECKER_CORRECTION.md explains the correction.
Verification imports no QHAT. It recomputes all 14 exact reference states and
replays all 420 new states plus the 84 reused baseline states with independent
occupation-basis Pauli action and an unmerged palindrome.

The calculation uses each full recorded Pauli-invariant affine closure, not
the full 2^n array and not a fixed-particle-number sector. The dimensions are
64 for CH2O and H2S, 128 for HCN and Li2, 256 for SiH4, 4096 for F2, and 16384
for NH3. Both variants of each molecule have the same dimension. All sectors
reachable by the emitted Pauli factors are retained, including leakage.

For an independent check of the already saved control, freshly extract the
combined review artifact, enter its directory containing run.py and
ALL_PLANS_FROZEN.json, and run independent_verify.py with the pinned environment. The original extension archive and QHAT are not needed
for this check. Existing verification checkpoints are accepted only with
matching row and checker source hashes; remove the verification/
subdirectory, verification.json, and READY.json in that NEW extraction if a fresh replay of every circuit is
wanted. Keep original evidence unchanged.

The 6 GiB address-space cap is set before numerical imports, and numerical
thread variables are one. Run numerical commands sequentially. The cap does
not limit unrelated system processes.

Production seconds measures propagation after pulse-list construction;
reference construction, design loading, NPZ output and analysis are excluded.
The first call may include JIT compilation. Reused baseline timing metadata
retains its original, broader measurement scope. These are descriptive
measurements under concurrent system use, not isolated speedup claims.
