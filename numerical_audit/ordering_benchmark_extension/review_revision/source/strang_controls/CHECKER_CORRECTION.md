# Independent-checker correction

All production sources and plans were frozen before the 420 new evolutions.
Those sources, plans, inputs and resulting states remain byte-identical.

The first independent verification attempt stopped before replaying any circuit
because its extra HF sanity check compared the physical occupation-basis integer
with `(1 << electrons) - 1`. That low-bit mask is the logical orbital mask used
by some algebra helpers, not the physical q0-most-significant-bit encoding of
the recorded states. For CH2O_s1.00, n=10 and electrons=4, the saved HF basis
integer is 960 (1111000000), not 15 (0000001111). The saved design explicitly
states `Hartree-Fock, q0 most significant bit; no sector projection`.

The correct physical mask is
`((1 << electrons) - 1) << (n_qubits - electrons)`.

The frozen run.py has not been modified. `independent_verify.py` copies its
independent verifier with that corrected guard and adds its own source hash to
each verification checkpoint. It imports the unchanged independent Pauli-action
implementation and imports no QHAT modules. Reproduction uses this separate
checker. No production result, target, schedule, coefficient, state, budget,
reference, metric, threshold, or orientation is changed by the correction.
The failed attempt produced no passing verification records.
