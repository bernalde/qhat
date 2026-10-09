# Suggested manuscript interpretation

The reconstructed H2O study uses the 15 configurations underlying slide 26:
six STO-6G models (4–14 qubits) and nine HGBS-5 models (4–20 qubits), all at
one fixed geometry. Within each model all three methods use the identical
retained Hamiltonian, Hartree–Fock determinant, T=1 and r=100. All retained
Pauli words occur once per step for each method, so these comparisons also
match nominal Pauli-exponential work. Parent-factor formulas with duplicates
are a separate experiment and must not inherit this equal-work assertion.

For the historical tensor-truncation policy, the proposed owned schedule
has 5 wins,4 ties and 6 losses against JW magnitude sorting, recovering the
slide 26 tally with a stable metric. All 45 historical-policy primary errors
(15 cases times 3 schedules) agree with the archived table within 3.92e-13
absolute; relative differences can reach 3.76% for floor-adjacent 2.5e-12
errors. Against JW signed sorting it has 5 wins,
2 ties and 8 losses. It is the unique best of the three in only 4/15 cases,
all HGBS-5 models with 14–20 active qubits. Large improvements in those models
(e.g. E=2.105e-9 vs3.007e-7 for JW magnitude at16 qubits) coexist with losses
elsewhere. These facts support specific historical comparisons, not a
uniform advantage of signed parent priority.

The three prescribed convergence anchors (STO-6G 8q, HGBS-5 12q and16q)
have E(r)/E(2r) in [3.99969,4.00054] for r=50,100 and all three schedules.
This is consistent with the leading O(r^-2) overlap error; it does not establish
an error law in active-qubit count. Along the active-space sweep the target
Hamiltonian, active electron count and frozen-core contribution change.

The two prescribed8-qubit models also have a state control in the full
(N_alpha,N_beta)=(2,2) sector. Both ground states are singlets; their next
states are Ms=0 triplets (<S²>=2); the equal coherent mixture has <S²>=1 and
variance1. Across the six model/state combinations, the proposed schedule
beats/ties/loses to JW signed2/0/4 and to JW magnitude5/0/1. These are repeated
states of two models, not six independent molecules. Full-space propagation
retains particle and spin-projection leakage. All 36 policy/method/state
propagations pass independent dense product checks.

The historical tensor files are unavailable. New tensors were regenerated from
frozen configurations and orbital provenance, so agreement of errors does not
prove original coefficient-identical inputs. The unpruned-policy sensitivity
results are reported separately; consult presentation/summary.json for the
completed numerical difference and validation maxima. In the first 14 cases
both retained Hamiltonian and schedules are exactly identical across policies;
in 20q the unpruned arm retains 4 additional tiny Pauli words and requires a
larger Pauli-invariant closure. Neither policy arm is an additional independent
molecular sample.

Both policies have the same outcomes. The maximum primary overlap-error change
from the tensor policy is 9.74e-21; the retained-H difference is zero in the
first 14 cases and coefficient l1 = 4.62985e-12 in the 20-qubit case.
The completed independent checker passed all 162 final states and 276 saved
arrays under Python -O, with maximum independent E discrepancy 1.36e-19.
Peak resident memory was 4.03 GiB.
