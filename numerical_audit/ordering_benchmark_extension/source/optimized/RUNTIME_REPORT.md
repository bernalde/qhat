# Runtime and work of the optimized campaign

All numbers are observations on this machine, with one thread per process and a 6 GiB address-space cap. Multiple independent jobs ran concurrently. Counts are reproducible algorithmic work measures; wall times are descriptive.

## Search across all 14 cases

| Arm | Examined insertion destinations | Total search seconds | Median seconds/case | Largest case seconds |
|---|---:|---:|---:|---:|
| class_refined_response | 52,446,880 | 24.214 | 0.117 | nh3_coeff1p25: 10.968 |
| class_refined_J0 | 52,446,880 | 24.722 | 0.124 | nh3_coeff1p25: 11.432 |
| class_pauli_JD | 1,255,213,280 | 1401.660 | 2.896 | nh3_coeff1p25: 659.789 |
| class_pauli_JD_budget | 52,446,880 | 57.920 | 0.205 | nh3_coeff1p25: 26.304 |

The conserving JD and J0 arms share one pair cache per case; the two individual-Pauli arms share another. Counting each observed setup once:

- refined_conserving: 2.008s total over 14 setups; median 0.009s/case.
- individual_pauli: 11.050s total over 14 setups; median 0.022s/case.

Search timers include initial pair-amplitude reconstruction and between-sweep checkpoint overhead, but exclude pair setup and starting-order construction. The matched Pauli arm has exactly the conserving JD destination count on every case; it is not an equal-wall-time control.

## Circuit measurement intervals

| Arm | States | Pauli rotations | Sum of recorded seconds |
|---|---:|---:|---:|
| class_refined_response | 98 | 3,260,672 | 87.070 |
| class_refined_J0 | 98 | 3,260,672 | 86.606 |
| class_pauli_JD | 98 | 3,260,672 | 86.482 |
| class_pauli_JD_budget | 98 | 3,260,672 | 86.016 |
| strang_class_refined_response | 84 | 3,219,573 | 86.617 |
| strang_jw_magnitude | 84 | 3,223,464 | 86.390 |
| strang_jw_signed | 84 | 3,223,464 | 86.712 |
| strang_class_lex | 84 | 3,179,008 | 83.581 |

The 728 saved rows emit 25,888,197 rotations in total. Their recorded intervals sum to 689.475s; this is not wall-clock elapsed time because jobs ran concurrently.

These intervals include pulse construction, state compression/output, and a design hash as well as evolution. They exclude error metrics and row writing. The fixed campaign uses evolve-only timers, so the two columns must not be treated as identical runtime measurements. See `timing_scope.json` for exact boundaries. The deliberate NH3 process split reused completed rows; work in an interrupted unfinished row is not in the totals.
The largest recorded propagation-process peak is 639.9 MiB. This is a whole-process peak, possibly including earlier cases, and is not an optimizer-memory measurement.
