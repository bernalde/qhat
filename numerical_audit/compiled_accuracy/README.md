# Compiled target-accuracy study preparation

Issue #9's [protocol](PROTOCOL.md) declares the exact-angle prototype compiler,
nearest-neighbor line, routing/cancellation, all eight strong comparators, both
times, fixed step grid, target errors and cumulative resource limits before any
execution. `compiler.py` retains every primitive gate and exact rational angle,
including explicit HF preparation, and counts CX and declared gate-depth metrics.

`test_compiler.py` checks full literal unitaries, every three-qubit Pauli support,
restored nonlocal routing, preparation phase, exact adjacent rewrites, trace
roundtrips and resource/depth accounting. These tests are prepared and UNEXECUTED.
The full frozen-grid controller, independent replay, reporting, delivery,
representative preflight, CI and paper companion are still required. No compiled
resource, accuracy, speedup, interval-proof or hardware result is claimed.

The guard preserves the supervisor status separately and returns failure after
an enforced budget stop. Every execution must also check termination metadata
and the required completed evidence artifact. No budget or review reset follows
from an interrupted wrapper or an unperformed test.
