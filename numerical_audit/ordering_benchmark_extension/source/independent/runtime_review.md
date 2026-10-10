# Independent runtime accounting review

Reviewed the frozen optimized design metadata, all 728 optimized numerical rows, `optimized/runtime_report.py`, `optimized/runtime_report.json`, `optimized/RUNTIME_REPORT.md`, `optimized/timing_scope.json`, `optimized/run.py`, `optimized/optimize.py`, `analysis/timings.py`, published preparation/propagation timing boundaries, and the expanded manuscript's methods, comparison text and `findings.tex`. This was read-only, with lightweight JSON arithmetic and no molecular propagation or optimizer rerun.

**Assessment: accounting and substantive claims pass.** Shared setup is correctly identified, no setup is added twice in the existing runtime report, and the new `timings.py` labels its repeated per-arm setup values rather than adding them. One minor rounding correction is recommended below. The final V2-dependent timing output could not yet be checked: 12 of 14 V2 designs existed at review time; both NH3 designs were still pending.

## Independent recomputation

The 14 optimized designs contain four optimizer reports apiece. For every report, the recorded search seconds, examined insertion destinations and exact rescores equal the sum of the two per-start records. Counts in the two matched-evaluation arms agree case by case.

| Search arm | Examined insertion destinations | Sum of recorded search seconds | Median seconds/case |
|---|---:|---:|---:|
| Conserving JD | 52,446,880 | 24.213548924 | 0.1172471045 |
| Conserving J0 | 52,446,880 | 24.721806321 | 0.1242370010 |
| Individual-Pauli JD | 1,255,213,280 | 1401.660072732 | 2.8964916725 |
| Evaluation-matched individual-Pauli JD | 52,446,880 | 57.920052617 | 0.2045574880 |

All these values reproduce `runtime_report.json` and its rounded Markdown table. The printed 1.25 billion in `findings.tex` is a small truncation: 1,255,213,280 rounds to **1.26 billion** at three significant digits. “About 1.3 billion” or the exact count is also fine. The 52.4 million figure is correctly rounded.

There are **28 actual observed pair-cache setups**, two per case. Conserving JD/J0 have identical shared setup metadata, and both Pauli arms have identical shared setup metadata in all cases. Summing one representative from each group gives:

- Conserving cache: 2.007839957 seconds across 14 setups.
- Pauli cache: 11.050132440 seconds across 14 setups.
- Combined unique setup intervals: 13.057972397 seconds, **not 26.115944794** seconds.

These are the observed shared preparations. They do not mean a standalone isolated arm has zero setup cost, nor do they include target loading, initial-order construction, J0's alternative weighting map, or propagation. The source and timing-scope description agree about these exclusions.

For NH3 at coefficient scale one, independently read from its frozen design:

| Quantity | Exact recorded seconds | Manuscript rounded value |
|---|---:|---:|
| Conserving JD search | 10.785411105 | 10.8 |
| Individual-Pauli JD search | 658.849626775 | 659 |
| Matched Pauli search | 25.456986393 | 25.5 |
| Shared conserving setup | 0.847386713 | 0.85 |
| Shared Pauli setup | 5.008991997 | 5.01 |

The 728 saved optimized rows contain 98 states for each of four first-order arms and 84 states for each of four symmetric arms. Their recorded intervals sum to **689.475447047 seconds**, and their emitted rotations sum to **25,888,197**. Each method's state count, rotation total and time sum agrees with the runtime report. This is summed work across concurrent jobs, not elapsed campaign wall time.

## Timer boundary correspondence

The source confirms the following:

1. Optimizer search seconds are the sum of the two per-start elapsed intervals. They include initial amplitude reconstruction and checkpoint overhead between timed sweeps. Final checkpoint writing after its timestamp, pair setup, starting-order construction and circuit errors are excluded. Resumption accumulates checkpointed elapsed work, but does not recover unfinished work lost before a checkpoint. The scope file states this limitation.
2. Optimized propagation rows start after input and per-method term compilation. They include pulse construction for Strang, evolution, compressed state writing and the design hash. The timestamp precedes the state hash and physical metrics, which are excluded. Python's left-to-right dictionary evaluation makes the documented boundary accurate.
3. Fixed rows time evolution only. Published first-order rows likewise exclude metrics/output; published symmetric rows add merged-pulse assembly. These intervals cannot be directly ranked as uniform kernel or end-to-end runtime measurements.
4. Native Tranter V2 preparation records one order per case for its two formula arms. `timings.py` reads only the first-order arm's shared metadata and therefore does not double-count it. The 12 existing plans correctly use zero native search time for the two completed positive-scaling reuses. Native process timing includes process startup, search/checkpoint work and parsing the result, but excludes compilation and earlier payload preparation. Latest-invocation timing is not total operational recovery cost. The NH3 final native interval still needs checking after completion.
5. The new timing analyzer preserves raw setup values and group labels, summarizes search intervals separately, and explicitly says the repeated setup values must not be summed. Its `selected_start_accepted_moves` and `accepted_moves_both_starts` are different, correctly labeled quantities; they must not be added to infer another total.

## Presentation and scope

The manuscript states that the search measurements are descriptive, single-threaded observations under concurrent jobs and explicitly denies an isolated speedup interpretation. The method appendix and planned TIMINGS.md explain differing propagation boundaries, cache reuse, checkpointing and objective-evaluation versus wall-time matching. No unsupported isolated-runtime-speedup claim was found in the expanded comparison text, findings or figure caption.

The plotted comparisons concern accuracy versus Pauli-rotation work; the limited classical-work control matches examined insertion positions, not elapsed time or setup. Different classical work per position makes the differing recorded times compatible with identical evaluation counts. The figure's least-tested-work cost is quantum circuit work and should remain separate from classical search time.

No repository files were changed in this review. The remaining check is to run the final timing analyzer once all V2 designs are available, then verify its emitted tables against these recomputed totals and inspect the final NH3 native/reuse metadata. No repeat numerical simulation is needed for that final accounting check.

## Reviewed-file hashes

- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/analysis/timings.py`: `5ca54e6d59078925f4060f5e91eeb6428323b38fe74f68ddeccc8450f7e0b894`
- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/optimized/timing_scope.json`: `2a3ed83fcf2c05db583f24263f78bf379e224590f7ba380e7b6ea0957a8d4c14`
- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/optimized/runtime_report.json`: `20e5044e6ff03c91f0be826dd1ac1e8ebfdeca62e44c63621552f686521d08b6`
- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/optimized/runtime_report.py`: `81dcef67c3fb6457d540d51e744bb3d6f2d00302da68a840f6233a136411401b`
- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/optimized/run.py`: `3d3254a24d0918121403440546bed2b53886dda8949b1c861f75e4e024638fc3`
- `/home/bernalde/repos/fermionic-ordering-research/benchmark-extension/optimized/optimize.py`: `40db8427892702b645d58c77429e2a56c529393d2fe5c90e49eef30e56e65db5`
- `/home/bernalde/repos/fermionic-ordering-research/paper/theory/expanded_benchmark_evidence/findings.tex`: `4e0cda2de10aeacd5d9bcfd7146a67f6ad7ee86dfb02c4479a93b3d97578ec17`
- `/home/bernalde/repos/fermionic-ordering-research/paper/theory/expanded_ordering_benchmarks.tex`: `8c054ac47a7be8bfa111c828ff963df608b2808ae3e6ee85971054ea30fdb4f3`
- `/home/bernalde/repos/fermionic-ordering-research/paper/theory/expanded_benchmark_methods.tex`: `0dab94ca5a529691e4a58a2db11ce2000db98bafea09b84106fc02c3b930111e`
