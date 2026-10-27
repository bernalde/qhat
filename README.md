# QHAT: Quantum Hamiltonian Analysis Toolkit

QHAT ("Q-hat") supports research into Hamiltonian simulation on quantum computers.  The current
pieces are

- [`hamiltonian_generator`](/hamiltonian_generator): generate a Hamiltonian from a description of a molecule
- [`analysis`](/analysis): analysis of quantum algorithms involving Hamiltonians
- [`common`](/common): modules used by other parts of QHAT
- [`julia_trotter`](/julia_trotter): julia scripts to compute trotter error for time evolution and ground state energy.

- [`time_evolution.py`](time_evolution.py): compare exact, first-order and symmetric time evolution with explicit reference states.
- [`numerical_audit`](numerical_audit): ordering studies with frozen protocols, implementations and reproducibility evidence, including [practical accuracy](numerical_audit/practical_accuracy), [state transfer](numerical_audit/transfer_holdouts) and [molecular certificates](numerical_audit/molecular_guard).

See the README file in each directory for component-specific documentation. The [independent project review response](numerical_audit/reports/claude_project_review_2026-10-27.md) records current verification scope and finding dispositions.

If you need support or have a feature request, the best way to reach out to the team is to open an
issue in the GitHub repository.

## Release

This software has been approved for open source release (disclosure number O4626) by the Feynman
Center for Innovation at Los Alamos National Laboratory.  Please see the license file for licensing
details.
