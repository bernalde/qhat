# Reconstructed configuration. Coordinates are Angstrom.
# This native QHAT replay writes its standard historical-policy tensors.
general.file_stub = 'H2O_s-1.00_hgbs-5'
general.output_directory = params['output_directory']
general.logfile = 'H2O_s-1.00_hgbs-5_as-004-008.log'
general.print_verbose()
hamiltonian.add_atom('O', 0.0, 0.0, 0.0)
hamiltonian.add_atom('H', 0.757540797781, 0.0, 0.587079636589)
hamiltonian.add_atom('H', -0.757540797781, 0.0, 0.587079636589)
hamiltonian.basis = 'hgbs-5'
hamiltonian.num_active_occupied = 4
hamiltonian.num_active_vacant = 8
hamiltonian.f2q_mapping = 'JW'
