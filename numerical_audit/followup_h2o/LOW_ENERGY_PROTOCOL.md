# H2O small low-energy supplement, frozen before state construction

Cases: STO-6G (4 occupied,4 virtual spin orbitals) and HGBS-5 (4,4), both
historical and unpruned tensor policies, using exactly the retained Hamiltonians
and three chronological schedules frozen in the main sweep design. Prepare
states in the full N_alpha=N_beta=2 sector, with no point-group restriction:
lowest and next canonical sector eigenvectors and (ground+next)/sqrt(2).
Use low_energy_state_sensitivity.py's existing 1e-10 Ha maximum-span degeneracy
clusters, coordinate-projector basis and phase rule (earliest coordinate within
1e-12 of largest magnitude positive real). Measure total S² and its variance;
fixing spin projection does not select a singlet. No changes to states based on
measured Trotter errors or total spin. Freeze state vectors and design before
approximate propagation.

T=1,r=100 for all three states and three fixed schedules: 36 approximate
evolutions, full 256-dimensional Hilbert space, all leakage retained. Exact
reference from full-space expm_multiply, checked against dense diagonalization.
Independent dense product matrices check each resulting approximate vector.
Metrics and tie/floor conventions are the main protocol's. This is a prescribed
small state control, not an independent molecular sample or an asymptotic test.
