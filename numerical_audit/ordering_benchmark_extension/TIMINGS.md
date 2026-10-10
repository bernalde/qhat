# Recorded runtimes and classical work

These are observed times from concurrent single-threaded research jobs, not an isolated speed comparison. See environment.json for hardware and software. All raw evolution rows retain their own seconds field.

## Optimizer setup and search

Pair-cache setup and search are separate. Conserving JD and J0 share one cache; the two Pauli arms share another. Repeated setup values in the CSV describe this shared cache and must not be summed as independent measurements. Search includes between-sweep checkpoint overhead. Evaluation matching limits examined insertion positions, not setup cost or wall time.

| Case | Conserving JD search (s) | Pauli JD search (s) | Evaluation-matched Pauli search (s) | Conserving / Pauli positions |
|---|---:|---:|---:|---:|
| CH2O_s1.00 | 0.091 | 0.842 | 0.099 | 156,640 / 1,255,000 |
| CH2O_s1.25 | 0.061 | 0.911 | 0.106 | 156,640 / 1,255,000 |
| H2S_s1.00 | 0.056 | 0.797 | 0.107 | 156,640 / 1,255,000 |
| H2S_s1.25 | 0.057 | 0.997 | 0.112 | 156,640 / 1,507,000 |
| HCN_s1.00 | 0.118 | 2.795 | 0.191 | 271,440 / 3,916,120 |
| HCN_s1.25 | 0.117 | 2.998 | 0.197 | 271,440 / 3,916,120 |
| SiH4_s1.00 | 0.289 | 11.569 | 0.561 | 651,600 / 14,061,640 |
| SiH4_s1.25 | 0.280 | 12.042 | 0.618 | 651,600 / 13,662,040 |
| f2_coeff1p00 | 0.523 | 22.544 | 1.906 | 2,415,120 / 28,013,440 |
| f2_coeff1p25 | 0.672 | 24.490 | 1.864 | 2,415,120 / 28,013,440 |
| li2_coeff1p00 | 0.101 | 1.388 | 0.186 | 290,400 / 2,119,000 |
| li2_coeff1p25 | 0.097 | 1.648 | 0.212 | 290,400 / 2,119,000 |
| nh3_coeff1p00 | 10.785 | 658.850 | 25.457 | 22,281,600 / 577,060,240 |
| nh3_coeff1p25 | 10.968 | 659.789 | 26.304 | 22,281,600 / 577,060,240 |

Full setup, both-start accepted moves, exact-positive rescores and objective counts are in optimizer_runtime.csv. Native Tranter V2 exact-search timings and reuse metadata are in timings.json. One V2 order supports two formula arms and its preparation is counted once. Exact positive-scaling reuse is not a fresh search.

## Propagation interval differences

- Fixed and published first-order row.seconds time evolution only, excluding metrics and output files. Compilation before evolution is excluded; first-call JIT inside evolution can be included.
- Published V2 Strang times merged-pulse assembly and evolution, excluding metrics and output files.
- Optimized row.seconds additionally includes pulse construction, compressed state output and the design hash; metrics and the state hash are excluded.
- These intervals must not be compared directly as algorithm timings. End-to-end case preparation has its own boundaries and may include verified structure reuse.

Raw seconds permit workload diagnosis. Accuracy-versus-rotation and accuracy-versus-evaluation comparisons remain the controlled comparisons in this campaign. A wall-time speedup claim would need fresh, repeated, isolated runs with matching measurement boundaries.
