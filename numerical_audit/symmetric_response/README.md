# Symmetric-formula-specific response ordering

Issue #8's formula-specific objective is derived in
[DERIVATION.md](DERIVATION.md); the prospective scientific and finite resource
limits are fixed in [PROTOCOL.md](PROTOCOL.md). `cubic.py` evaluates the actual
symmetric-step Taylor column and scores its diagonal-reference response divided
by the fourth power of the actual affordable fused step count. The independent
literal-Pauli verifier instead uses nested commutators, checks full signed columns,
and replays all circuits without fusion in the full Hilbert space.

The prepared controller retains fresh same-time first-order conserving and
larger individual-Pauli controls, their symmetric reuse, legacy placement,
cubic orientation-only and bounded cubic insertion. It freezes every circuit
before physical errors and records complete constructor charges and search work.

Use the existing pinned Python3.12 environment, one guarded worker/thread and
the quarter-core/512MiB/no-swap/headroom limits. Commit source before preflight
and append calibrated main estimates before main execution. Failed attempts
remain charged to the same phase. The completed main panel and independent verification are recorded below; the frozen derivation/protocol retain their prospective wording.

```
python -m pytest test_cubic.py test_symmetric_guard.py -q
python preflight.py --output FRESH_PREFLIGHT
python study.py prepare --output FRESH_RESULTS
python study.py run --output FRESH_RESULTS
python verify.py --output FRESH_RESULTS
python reporting.py --results FRESH_RESULTS --report FRESH_REPORT
python delivery.py pack --results FRESH_RESULTS --report FRESH_REPORT --delivery FRESH_DELIVERY
python delivery.py check --delivery FRESH_DELIVERY
```

These commands must run through guarded.py within the cumulative protocol
allowances. No accuracy, speedup, formal-proof, compiled-resource or hardware
conclusion follows from source preparation or a feasibility preflight.

The resource guard preserves the supervisor exit status separately and returns
failure after any enforced CPU/disk/wall stop. Require both no termination and
the completed validation artifact; a wrapper exit alone is not a scientific pass.


## Completed prospective panel

The complete fixed panel retains 96 full-space circuits on three supplied
molecular targets at both times and both rotation caps. All 336 pairwise and
floor decisions agree with independent replay. Maximum state/metric differences
are 7.692e-15/7.661e-15; the separate nested-commutator cubic-column difference
is 2.648e-15. See [all outcomes and costs](evidence/report/README.md).

Across the 12 descriptive case/time/cap pairs, cubic insertion versus symmetric
reuse of fresh first-order orders gives 7/5/0 wins/ties/losses; versus the
orientation-only cubic rule it also gives 7/5/0, and versus the larger-search
Pauli symmetric control 11/1/0. Floor-limited rows remain. Comparisons against
first-order formulas also change the formula and cannot isolate ordering.
The four cell-wise median constructor-cost ratios against the fresh conserving
initializer range from 8.36 to 10.19, with upstream/shared/cache costs included.
No matched-total-CPU advantage or practical speedup is established.

All exact schedules and 96 binary states are retained under evidence/;
delivery.py authenticates their manifest and byte-regenerates every report.
execution_resources.json records measured execution and fixed phase ceilings.
Full archive/extraction checks passed. This is numerical fixed-family evidence,
separate from Lean proofs, compiled gates/depth and hardware measurements.
