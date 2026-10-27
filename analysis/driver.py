# Disable qsharp's telemetry.  qsharp checks the environment variable
# QSHARP_PYTHON_TELEMETRY on import, so this needs to appear before
# that import occurs.
import os

os.environ["QSHARP_PYTHON_TELEMETRY"] = "none"

import logging
import math

from qhat.common.logging_utils import configure_logging
from qhat.analysis.algorithm import build_algorithm
from qhat.analysis.analysis import analyze_algorithm
from qhat.analysis.configuration import load_configuration
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.unitary import encode_as_unitary
from qhat.analysis.phase_parameters import automatic_phase_time, phase_qubits_for_time

logger = logging.getLogger("qhat.analysis.driver")

# =================================================================================================

def run():

    # Configuration _______________________________________________________________________________

    state, config_file_message, params_message = load_configuration()
    # Carry the same canonical spelling through encoding and downstream analysis.
    state.config_unitary.method = state.config_unitary.method.strip().lower()

    # Configure logging based on user settings
    logfile_path = state.config_general.get_output_path(state.config_general.logfile)
    configure_logging(level=state.config_general.loglevel, logfile=logfile_path)

    logger.info("=" * 99)
    logger.info("ANALYSIS DRIVER START")
    logger.info("=" * 99)
    logger.info(f"Logfile: {logfile_path}")
    logger.info(f"Git hash: {state.config_general.git_hash}")
    logger.info(config_file_message)
    if params_message:
        logger.info(params_message)

    # Hamiltonian _________________________________________________________________________________

    physical_hamiltonian = get_physical_hamiltonian(state.config_hamiltonian)

    # Compute Trotterization parameters ___________________________________________________________

    tevol_hbar = None
    P0 = None

    if state.config_unitary.method.strip().lower() == "ramped trotter":

        # first-pass computation of energy bounds
        Elo1, Ehi1 = physical_hamiltonian.compute_initial_energy_bounds(state.config_hamiltonian)

        # energy-shift Hamiltonian to center at zero
        # This maps eigenvalues from [Elo1, Ehi1] to [-(Ehi1-Elo1)/2, +(Ehi1-Elo1)/2]
        # With phase_scale_factor > 1, this enables phase angles in approximately [-π/s, +π/s]
        # where s is the scale factor. This ensures phases never hit exactly ±π, avoiding
        # aliasing ambiguity, while matching np.angle output range directly.
        E0 = (Elo1 + Ehi1) / 2
        physical_hamiltonian.energy_shift(-E0)
        Elo2 = Elo1 - E0
        Ehi2 = Ehi1 - E0
        logger.verbose(f"-- shifted bounds = [{Elo2}, {Ehi2})")
        # Apply phase scale factor to avoid ambiguity at ±π
        phase_scale = getattr(state.config_unitary, 'phase_scale_factor', 1.01)
        # Use the largest time within the physical phase margin. Padding the
        # energy interval to a power of two before sizing strict phase precision
        # can add a needless phase qubit at the default margin.
        tevol_hbar = automatic_phase_time(Elo2, Ehi2, phase_scale)
        logger.verbose(f"-- phase scale factor = {phase_scale}")
        logger.verbose(f"-- automatic evolution time = {tevol_hbar} * hbar")

        # check for user-defined timestep
        if getattr(state.config_unitary, 'timestep', None) is not None:
            tevol_hbar = state.config_unitary.timestep
            logger.verbose(f"-- user timestep override = {tevol_hbar} * hbar")

        # Every evolution needs a finite positive time. Phase branch cuts only
        # constrain QPE and analyses that infer energies from unitary phases.
        if not math.isfinite(tevol_hbar) or tevol_hbar <= 0:
            raise ValueError("timestep must be finite and positive")
        qpe_methods = ("qpe: qualtran textbook", "qpe: qualtran qubitization",
                       "qpe: pyliqtr qubitized")
        needs_unaliased_phase = (
            state.config_algorithm.method.strip().lower() in qpe_methods
            or getattr(state.config_analysis, 'enable_eigenvalue_errors', False)
            or any(request['source'] == 'approximate'
                   and request['operator_type'] == 'hamiltonian'
                   for request in getattr(state.config_analysis,
                                          '_operator_output_requests', [])))
        if needs_unaliased_phase:
            P0 = phase_qubits_for_time(tevol_hbar, Elo2, Ehi2, phase_scale,
                                      state.config_algorithm.energy_error)


    # Unitary _____________________________________________________________________________________

    unitary_hamiltonian = encode_as_unitary(
            state.config_unitary,
            physical_hamiltonian,
            tevol_hbar)

    # Algorithm ___________________________________________________________________________________

    algorithm = build_algorithm(
            state.config_algorithm,
            unitary_hamiltonian,
            P0)

    # Analysis ____________________________________________________________________________________

    state.store_results(analyze_algorithm(
        state.config_analysis,
        algorithm,
        state.config_unitary.method,
        approximate_time_evolution=unitary_hamiltonian,
        exact_hamiltonian=physical_hamiltonian,
        timestep=tevol_hbar,
        energy_shift=physical_hamiltonian.get_energy_shift(),
        config_general=state.config_general
    ))

    # Save Results ________________________________________________________________________________

    state.show_results()
    state.save_summary()

    logger.results("Analysis complete.")

# =================================================================================================

if __name__ == "__main__":
    run()
