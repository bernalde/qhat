# Exact structural certificates and compilation illustration

Run `PYTHONPATH=/home/bernalde/repos /tmp/qhat-h2o-venv/bin/python structural_compilation.py --output /tmp/structural-compilation-replay` from this directory. Use the environment pinned in ../h2o_recovery/environment-requirements.txt. The script sets two-thread limits; run under a 6 GiB memory cap.

The committed results include executed source and SHA-256. Exhaustive four-mode certificates are independently checked by integer CAR matrices and JW matrices; the general theorem is in the paper's Lean development. Unknown means uncertified, not noncommuting. The 12-to-8 rotation compilation combines angles only within a certified contiguous region and preserves the unitary, hence the Trotter error. Hardware assumptions and gate counts are recorded in results/results.json. This is one illustration, not a molecular performance claim.
