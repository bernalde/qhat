# Sparse molecular ordering: evidence and practical interpretation

**Independent numerical verification passed:** 24 frozen orders, 96 complete
circuit outcomes, eight reference evolutions and all 240 paired comparisons.
The largest difference between the two reference calculations is 3.3e-15.
The final scientific review and hosted CI are recorded separately in the PR
history; numerical verification alone does not supply those gates.

This experiment chooses a **per-instance order of commuting molecular units or
individual Pauli terms in fixed orbitals**. It does not optimize molecular
orbitals, spin-orbital indexing or the fermion-to-qubit mapping. Four STO-3G
active-space cases—N2 at 1.10 and 2.00 Å, CO at 1.13 Å and tetrahedral CH4 at
1.09 Å—each use 12 spin modes, six active electrons and an HF reference. These
are new recorded-baseline geometries in three families; CH4 was already a
baseline family. The two N2 geometries and repeated times/budgets are dependent
observations.

All six arms execute the symmetric second-order formula. The first-response and
individual-Pauli arms use a first-response **search objective**, not first-order
circuits. The cubic arm trains its sparse response at time 1 with a 30U rotation
cap; the frozen orders are then evaluated at times 0.5/1.5 and caps 10U/30U.
U counts retained Pauli words. Training reads no physical validation error.
Changing the individual-Pauli factorization also changes unit boundaries, so
these comparisons do not establish a causal effect of number conservation.

The **independently replayed** state-error metric, `1-|overlap|`, places cubic
lowest or tied in all 16 case/time/cap settings. Against first-response it has
11 wins and five ties; against cheap current placement it has 15 wins and one
tie, using the fixed absolute ambiguity of 1e-12. CO shows roughly twofold lower
error than first-response with nearly the same compiled CX/depth. Stretched N2
ties first-response. CH4 gains only about 1% over cheap current placement, with
one ambiguity tie, equal CX and slightly greater depth. This makes the extra
classical work consequential to the practical comparison:

| Case | Cubic construction CPU s | Completed seeds | Accepted moves |
| --- | ---: | ---: | ---: |
| N2, 1.10 Å | 3.3762 | 3 | 1 |
| N2, 2.00 Å | 5.8319 | 3 | 1 |
| CO, 1.13 Å | 6.4049 | 3 | 1 |
| CH4, 1.09 Å | 29.9851 | 3 | 3 |

These are complete standalone construction charges, including setup and local
cleanup, rather than elapsed wall times. Returned scores/incumbents are
completed; the configured 32-CPU-second allowance does not establish search
convergence. Exact propagation and 4096-state checking are offline audit costs.
The response remains a heuristic, with no literature-novelty, universal-ordering,
general finite-error guarantee, chemistry-model accuracy or hardware claim. Compiled
costs concern the declared nearest-neighbor-line primitive compiler, including
HF preparation.

Orders were trained/frozen at `a9e8f504`; execution uses `a21a449` after a
row-packaging correction. `SOURCE_AMENDMENT.json` and both freeze manifests
preserve all 24 identical constructor hashes, with no retraining or retuning.
The failed first execution produced zero error rows and remains archived with
its reference, trace and costs. The whole CPU limit stays 3,600 seconds. Before
independent replay, 220 unused measurement CPU seconds were prospectively
transferred to its allocation, now 1,200 seconds. All attempted-work costs,
including failures and continuations, remain in the cumulative ledger.

## Navigate the evidence

- `report/REPORT.md` and `report/report_tables.json`: raw six-arm settings,
  paired ratios/outcomes, actual search work and source/interpreter provenance.
- `INDEPENDENT_VERIFIED.json`: numerical replay, DAG/state/reference checks and
  the maintained pair/summary arithmetic gate.
- `main-study.tar.gz` and its inventory: unchanged models, tensors/orbitals,
  constructors, objective witnesses, traces/states, references and source lineage.
- `retained-attempts.tar.gz`, its inventory and `CAMPAIGN_COSTS.json`: earlier
  attempts, terminal guard receipts/logs, calibration and cumulative charges.
- `PACKAGE_MANIFEST.json` and `tools/`: payload hashes and exact renderer,
  packager, recovery/profiling scripts and guard. Archive verification compares
  the complete regular-file inventory and member bytes; raw data stay unchanged.
  The final review assessment and live CI are recorded in the PR history.

## Replay or rerender

Check archive/member hashes against the inventories before unpacking into a
fresh directory. The independent command requires scientific source bytes from
`a21a449` and the recorded Python 3.12/package/platform environment. Its strict
environment gate may reject a relocated or different host; preserve the original
records rather than rewriting their environment. Run heavy replay under the
retained guard with a separately recorded budget, one worker, quarter-core CPU,
512 MiB RAM, no swap and the recorded memory-headroom rules.

```sh
QHAT_EVIDENCE=/path/to/numerical_audit/scalable_ordering/evidence
QHAT_SOURCE=/path/to/qhat
QHAT_REPLAY=/fresh/replay-directory
QHAT_PYTHON=/path/to/matching-python-3.12
mkdir "$QHAT_REPLAY"
tar -xzf "$QHAT_EVIDENCE/main-study.tar.gz" -C "$QHAT_REPLAY"

# Heavy replay through the Linux/systemd user-service resource guard.
python3 "$QHAT_EVIDENCE/tools/claude-guarded.py" \
  --label scalable-retained-replay --logs "$QHAT_REPLAY/logs" \
  --cpu-seconds 1200 --wall-seconds 7200 --disk-mib 1536 \
  --output "$QHAT_REPLAY" -- \
  "$QHAT_PYTHON" "$QHAT_SOURCE/numerical_audit/scalable_ordering/verify.py" \
  --output "$QHAT_REPLAY/study-a2" \
  --destination "$QHAT_REPLAY/verification-fresh"

# Rerender saved, already-verified evidence without project/numerical imports.
python3 "$QHAT_EVIDENCE/tools/qhat-scalable-render-report.py" --study "$QHAT_REPLAY/study-a2" --verification "$QHAT_EVIDENCE" --destination "$QHAT_REPLAY/report-fresh" --project-root "$QHAT_SOURCE"
```

Rerendering uses only the standard library. Relocated paths and a different
interpreter change JSON provenance, so a relocated rerender is not claimed to be
whole-file byte-identical to the original report. Keep original artifacts and
all new outputs separate. A fresh training run is a separate reproduction; it
must not replace the retained orders or historical costs.

Historical packaging, recovery and profiling tools retain their recorded absolute
campaign paths. Inspect those paths before reuse. The renderer accepts paths on
the command line and can be used from a relocated evidence copy.

## Review record

The original [Astra round 1](reviews/astra-r1/REPORT.md) and
[targeted round 2](reviews/astra-r2/REPORT.md) reports are retained verbatim.
[Coordinator dispositions](reviews/astra-r1/DISPOSITIONS.md),
[verification pass 1](reviews/VERIFY_1.md) and
[verification pass 2](reviews/VERIFY_2.md) are separate historical records; their
then-pending main-study gates are satisfied by the evidence above. The final
evidence-focused review is published in the PR review history.

`CAMPAIGN_COSTS.json` covers terminal guarded attempts through report generation,
including failures. Packaging, later publication, coordinator/metadata overhead
and remote review/CI lie outside that snapshot; it is not whole-machine CPU.
