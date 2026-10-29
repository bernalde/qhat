"""Phase-domain validation and precision at the actual chosen evolution time."""
import math


def phase_qubits_for_time(time, lower, upper, phase_scale, energy_error):
    if not math.isfinite(time) or time <= 0:
        raise ValueError("timestep must be finite and positive")
    if not math.isfinite(phase_scale) or phase_scale <= 0:
        raise ValueError("phase_scale_factor must be positive and finite")
    if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
        raise ValueError("Finite increasing physical energy bounds are required")
    radius = max(abs(lower), abs(upper))
    if time >= math.pi / radius:
        raise ValueError("timestep reaches the phase branch cut; choose a smaller timestep")
    if time > math.pi / (phase_scale * radius):
        raise ValueError("timestep exceeds the configured phase_scale_factor margin")
    if energy_error is None:
        return None
    if not math.isfinite(energy_error) or energy_error <= 0:
        raise ValueError("energy_error must be finite and positive")
    return max(0, math.ceil(math.log2(2*math.pi / (time*energy_error))))


def automatic_phase_time(lower, upper, phase_scale):
    """Largest automatic time inside the configured and physical phase limits.

    A positive factor below one does not provide a branch-cut margin by itself;
    clamp that case (and factor one) just below the physical branch cut.
    """
    if not math.isfinite(phase_scale) or phase_scale <= 0:
        raise ValueError("phase_scale_factor must be positive and finite")
    if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
        raise ValueError("Finite increasing physical energy bounds are required")
    radius = max(abs(lower), abs(upper))
    return min(math.pi / (phase_scale * radius),
               math.nextafter(math.pi / radius, 0.0))
