# Reconstructed configuration. Coordinates are Angstrom.
# This native QHAT replay writes its standard historical-policy tensors.
general.file_stub = 'H2O_s-1.00_sto-6g'
general.output_directory = params['output_directory']
general.logfile = 'H2O_s-1.00_sto-6g_as-006-004.log'
general.print_verbose()
hamiltonian.add_atom('O', 0.0, 0.0, 0.0)
hamiltonian.add_atom('H', 0.757540797781, 0.0, 0.587079636589)
hamiltonian.add_atom('H', -0.757540797781, 0.0, 0.587079636589)
hamiltonian.basis = 'sto-6g'
hamiltonian.num_active_occupied = 6
hamiltonian.num_active_vacant = 4
hamiltonian.f2q_mapping = 'JW'
