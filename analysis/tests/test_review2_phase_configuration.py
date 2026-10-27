"""Regression guards for reviewed phase-domain and estimator configuration behavior."""
import math
from types import SimpleNamespace

import pytest
from qhat.analysis.config_types import UnitaryConfiguration
from qhat.analysis.phase_parameters import automatic_phase_time, phase_qubits_for_time


def test_strict_default_precision_does_not_add_a_padding_qubit():
    time = automatic_phase_time(-1, 1, 1.01)
    precision = phase_qubits_for_time(time, -1, 1, 1.01, .1)
    assert time == pytest.approx(math.pi / 1.01)
    assert precision == 5
    assert 2 * math.pi / (time * 2**precision) <= .1
    assert 2 * math.pi / (time * 2**(precision - 1)) > .1


@pytest.mark.parametrize('scale', [.5, 1.0])
def test_small_margin_factor_keeps_safe_override_and_clamps_automatic_time(scale):
    automatic = automatic_phase_time(-1, 1, scale)
    assert 0 < automatic < math.pi
    assert phase_qubits_for_time(.5, -1, 1, scale, .1) == 7
    with pytest.raises(ValueError, match='branch cut'):
        phase_qubits_for_time(math.pi, -1, 1, scale, .1)


@pytest.mark.parametrize('scale', [0, -1, math.inf, math.nan])
def test_nonfinite_or_nonpositive_scale_rejected_by_factory_and_helpers(scale):
    for call in (
        lambda: UnitaryConfiguration().encode_ramped_trotter(phase_scale_factor=scale),
        lambda: automatic_phase_time(-1, 1, scale),
        lambda: phase_qubits_for_time(.5, -1, 1, scale, .1),
    ):
        with pytest.raises(ValueError, match='phase_scale_factor must be positive'):
            call()


@pytest.mark.parametrize('options,expected', [
    ({}, ('monte_carlo', False, 60)),
    ({'error_coeff_mode': 'exact', 'error_coeff_auto_exact': True,
      'error_coeff_time_limit': .25}, ('exact', True, .25)),
])
def test_coefficient_controls_survive_factory_and_toml(options, expected):
    cfg = UnitaryConfiguration()
    cfg.encode_ramped_trotter(energy_error=.1, **options)
    table = cfg._generate_TOML_table()
    names = ('error_coeff_mode', 'error_coeff_auto_exact', 'error_coeff_time_limit')
    assert tuple(getattr(cfg, name) for name in names) == expected
    assert tuple(table[name] for name in names) == expected


@pytest.mark.parametrize('options', [
    {'error_coeff_mode': 'typo'}, {'error_coeff_auto_exact': 'false'},
    {'error_coeff_time_limit': 0}, {'error_coeff_time_limit': -1},
    {'error_coeff_time_limit': math.nan}, {'error_coeff_time_limit': math.inf},
    {'error_coeff_time_limit': True},
])
def test_invalid_coefficient_controls_rejected_instead_of_silent_fallback(options):
    with pytest.raises(ValueError):
        UnitaryConfiguration().encode_ramped_trotter(energy_error=.1, **options)


def test_unknown_factory_keyword_fails_clearly():
    with pytest.raises(TypeError, match='error_coeff_mod'):
        UnitaryConfiguration().encode_ramped_trotter(energy_error=.1, error_coeff_mod='exact')


def _capture_driver(monkeypatch, *, method='QPE: qualtran textbook', timestep=None,
                    phase_scale=1.01, eigenvalue_errors=False):
    import qhat.analysis.driver as driver
    observed = {}
    general = SimpleNamespace(get_output_path=lambda name: name, logfile='unused.log',
                              loglevel='warning', git_hash='fixture')
    unitary = SimpleNamespace(method='ramped trotter', phase_scale_factor=phase_scale,
                              timestep=timestep)
    state = SimpleNamespace(config_general=general, config_unitary=unitary,
        config_algorithm=SimpleNamespace(method=method, energy_error=.1),
        config_hamiltonian=object(),
        config_analysis=SimpleNamespace(enable_eigenvalue_errors=eigenvalue_errors),
        store_results=lambda results: None, show_results=lambda: None, save_summary=lambda: None)
    h = SimpleNamespace(compute_initial_energy_bounds=lambda cfg: (-1, 1),
                        energy_shift=lambda shift: None, get_energy_shift=lambda: 0)
    monkeypatch.setattr(driver, 'load_configuration', lambda: (state, 'fixture', None))
    monkeypatch.setattr(driver, 'configure_logging', lambda **kwargs: None)
    monkeypatch.setattr(driver, 'get_physical_hamiltonian', lambda cfg: h)
    def encode(cfg, physical, time):
        observed['time'] = time
        return object()
    def build(cfg, encoded, p0):
        observed['precision'] = p0
        return object()
    monkeypatch.setattr(driver, 'encode_as_unitary', encode)
    monkeypatch.setattr(driver, 'build_algorithm', build)
    monkeypatch.setattr(driver, 'analyze_algorithm', lambda *args, **kwargs: {})
    driver.run()
    return observed


def test_driver_default_uses_strict_precision_at_unpadded_time(monkeypatch):
    result = _capture_driver(monkeypatch)
    assert result['time'] == pytest.approx(math.pi / 1.01)
    assert result['precision'] == 5


@pytest.mark.parametrize('method', ['time evolution', 'controlled time evolution'])
def test_non_qpe_long_override_is_valid_without_phase_energy_analysis(monkeypatch, method):
    assert _capture_driver(monkeypatch, method=method, timestep=10) == {
        'time': 10, 'precision': None}


@pytest.mark.parametrize('method,eigenvalue_errors', [
    ('QPE: qualtran textbook', False), ('QPE: qualtran qubitization', False),
    ('QPE: pyliqtr qubitized', False), ('time evolution', True),
    ('controlled time evolution', True),
])
def test_phase_consumers_reject_long_wrapping_override(monkeypatch, method, eigenvalue_errors):
    with pytest.raises(ValueError, match='branch cut'):
        _capture_driver(monkeypatch, method=method, timestep=10,
                        eigenvalue_errors=eigenvalue_errors)


def test_driver_accepts_safe_override_with_factor_below_one(monkeypatch):
    assert _capture_driver(monkeypatch, phase_scale=.5, timestep=.5) == {
        'time': .5, 'precision': 7}


@pytest.mark.parametrize('time', [0, -1, math.inf, math.nan])
def test_even_non_qpe_rejects_nonpositive_or_nonfinite_time(monkeypatch, time):
    with pytest.raises(ValueError, match='timestep must be finite and positive'):
        _capture_driver(monkeypatch, method='time evolution', timestep=time)


def test_exact_factory_controls_reach_estimator(monkeypatch):
    import qhat.analysis.unitary as unitary
    import qhat.analysis.trotter_coefficients_fast as fast
    import qhat.common.trotter_flattened as builder
    from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings
    cfg = UnitaryConfiguration()
    cfg.encode_ramped_trotter(energy_error=1, trotter_order='second order',
        error_coeff_mode='exact', error_coeff_auto_exact=True, error_coeff_time_limit=.25)
    h = Hamiltonian(LinearCombinationOfPauliStrings(dense={'Z': 2, 'X': 1}, num_qubits=1))
    seen = {}
    def estimate(terms, time_limit, **kwargs):
        seen.update(time_limit=time_limit, **kwargs)
        return 2, 5/3
    sentinel = object()
    monkeypatch.setattr(fast, 'trotter_error_estimator_fast', estimate)
    monkeypatch.setattr(builder, 'build_ramped_trotterized_unitary', lambda *args, **kwargs: sentinel)
    assert unitary.encode_ramped_trotter(cfg, h, 1) is sentinel
    assert seen == {'time_limit': .25, 'mode': 'exact', 'auto_exact': True}
