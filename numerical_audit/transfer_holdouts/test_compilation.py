"""Literal-matrix checks for transfer states, compiled traces and cost boundaries."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
from scipy.linalg import expm
import compilation as c
import circuit_reference as independent
import training
import transfer_common as k
from test_states import toy, literal


def target_record(model, formula, partition='conserving', reverse=False):
    size = len(model['blocks']) if partition == 'conserving' else len(model['coefficients'])
    order = list(range(size))
    if reverse:
        order.reverse()
    return dict(training.layout(model, order, partition), formula=formula)


def direct_formula(model, record, psi, T, steps):
    """Unfused dense exponentials, with no producer fusion or Pauli action."""
    offsets = np.cumsum([0] + record['unit_lengths'])
    units = []
    for i in range(len(record['unit_lengths'])):
        matrix = sum(int(coefficient) / int(model['denominator']) * literal(word, model['n'])
                     for word, coefficient in record['schedule'][offsets[i]:offsets[i + 1]])
        units.append(expm(-1j * T / steps / (2 if record['formula'] == 2 else 1) * matrix))
    state = psi.copy()
    for _ in range(steps):
        if record['formula'] == 1:
            for unit in units:
                state = unit @ state
        else:
            for unit in units[:-1]:
                state = unit @ state
            state = units[-1] @ (units[-1] @ state)
            for unit in reversed(units[:-1]):
                state = unit @ state
    return state


@pytest.mark.parametrize('formula', [1, 2])
@pytest.mark.parametrize('label', ['HF', 'coherent_real', 'coherent_imag'])
@pytest.mark.parametrize('reverse', [False, True])
def test_complete_trace_preserves_initial_state_and_matches_dense_and_sparse_replays(
        tmp_path, formula, label, reverse):
    model = toy()
    psi = np.zeros(16, complex)
    psi[3] = 1
    if label != 'HF':
        psi[3], psi[12] = 1 / np.sqrt(2), (1 if label == 'coherent_real' else 1j) / np.sqrt(2)
    path, trace = tmp_path / 'model.json', tmp_path / 'trace.ndjson.gz'
    k.write(path, model)
    record = target_record(model, formula, reverse=reverse)
    result = c.construct_candidate(path, k.sha(path), record, .37, 10, label, trace, {'case': 'toy'})
    assert result['status'] == 'complete'
    actual = c.Replay(model, psi).replay(trace, result['identity'])
    expected = direct_formula(model, record, psi, .37, result['steps'])
    reference = independent.Reference(model, psi)
    checked, resources = reference.audit_compiled_trace(trace, result['identity'], c.compiler.VERSION)
    assert np.max(abs(actual - expected)) < 2e-13
    assert np.max(abs(checked - expected)) < 2e-13
    assert np.max(abs(reference.ideal_formula(record, result['steps'], .37) - expected)) < 2e-13
    assert resources == result['compiled_resources']
    H = sum(int(coefficient) / 8 * literal(word, 4) for word, coefficient in model['coefficients'])
    exact = expm(-.37j * H) @ psi
    assert np.max(abs(reference.exact_action(.37) - exact)) < 2e-13
    assert c.metrics(exact, actual, model) == pytest.approx(independent.metrics(exact, checked, model), abs=2e-13)
    assert result['quantum_preparation_circuit_cost'] == ('HF prefix included' if label == 'HF' else 'unknown')


@pytest.mark.parametrize('formula', [1, 2])
@pytest.mark.parametrize('lengths', [[2], [2, 1], [1, 2, 3], [3, 1, 2]])
@pytest.mark.parametrize('cap', [10, 30])
def test_affordability_and_fusion_match_unfused_unit_sequence(formula, lengths, cap):
    steps = c.affordable_steps(lengths, cap, formula)
    expanded = []
    for _ in range(steps):
        if formula == 1:
            expanded.extend((i, 1 / steps) for i in range(len(lengths)))
        else:
            expanded.extend((i, 1 / (2 * steps)) for i in range(len(lengths) - 1))
            expanded.append((len(lengths) - 1, 1 / steps))
            expanded.extend((i, 1 / (2 * steps)) for i in reversed(range(len(lengths) - 1)))
    combined = []
    for unit, tau in expanded:
        if combined and formula == 2 and combined[-1][0] == unit:
            combined[-1] = (unit, combined[-1][1] + tau)
        else:
            combined.append((unit, tau))
    fused = list(c.fused_units(lengths, steps, formula))
    assert [i for i, _ in combined] == [i for i, _ in fused]
    assert [tau for _, tau in combined] == pytest.approx([float(tau) for _, tau in fused], abs=1e-14)
    counted = sum(lengths[i] for i, _ in combined)
    assert counted == c.rotation_count(lengths, steps, formula) <= cap * sum(lengths)
    if formula == 1 or len(lengths) > 1:
        # Independent expansion of one extra step must exceed the same cap.
        extra = ((steps + 1) * sum(lengths) if formula == 1 else
                 (steps + 1) * (2 * sum(lengths) - lengths[0] - lengths[-1]) + lengths[0])
        assert extra > cap * sum(lengths)


def test_leakage_from_split_pauli_formula_is_retained_in_full_space(tmp_path):
    # XX+YY annihilates an empty pair. Separating its two words by Z creates
    # real number leakage; projecting back to the reference sector would hide it.
    terms = [[[[2, 'Z']], '3'], [[[2, 'X'], [4, 'X']], '2'], [[[2, 'Y'], [4, 'Y']], '2']]
    model = dict(case='leakage_toy', n=6, electrons=2, denominator='5',
                 coefficients=terms, blocks=[[terms[0]], terms[1:]])
    blocks = [tuple((q, p) for q, p in word) for word, _ in sorted(
        [(tuple((q, p) for q, p in word), coefficient) for word, coefficient in terms])]
    order = [blocks.index(((2, 'X'), (4, 'X'))), blocks.index(((2, 'Z'),)),
             blocks.index(((2, 'Y'), (4, 'Y')))]
    record = dict(training.layout(model, order, 'pauli'), formula=1)
    path, trace = tmp_path / 'model.json', tmp_path / 'trace.ndjson.gz'
    k.write(path, model)
    result = c.construct_candidate(path, k.sha(path), record, .37, 10, 'HF', trace, {'case': 'leakage_toy'})
    psi = np.zeros(64, complex)
    psi[3] = 1
    state = c.Replay(model, psi).replay(trace, result['identity'])
    reference = independent.Reference(model, psi)
    checked, _ = reference.audit_compiled_trace(trace, result['identity'], c.compiler.VERSION)
    assert np.max(abs(state - checked)) < 2e-13
    values = c.metrics(reference.exact_action(.37), state, model)
    assert values['number_leakage'] > 1e-12
    assert values['joint_sector_leakage'] == pytest.approx(values['number_leakage'], abs=1e-13)
    assert len(state) == 64 and abs(np.linalg.norm(state) - 1) < 2e-13


def test_oversized_model_and_tampered_layout_fail_before_compiler_allocation(tmp_path, monkeypatch):
    model, path = toy(), tmp_path / 'model.json'
    record = target_record(model, 2)
    def forbidden(*args, **kwargs):
        pytest.fail('Compiler allocation must follow input validation')
    monkeypatch.setattr(c.compiler, 'Compiler', forbidden)
    record['schedule'][0][1] = '999'
    k.write(path, model)
    with pytest.raises(ValueError):
        c.construct_candidate(path, k.sha(path), record, .37, 10, 'HF', tmp_path / 'bad.gz', {})
    model['n'] = 7
    k.write(path, model)
    with pytest.raises(ValueError, match='Register bound'):
        c.construct_candidate(path, k.sha(path), record, .37, 10, 'HF', tmp_path / 'bad2.gz', {})


def test_reference_shape_checked_before_dtype_copy(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Dtype allocation must follow shape validation')
    monkeypatch.setattr(c.np, 'asarray', forbidden)
    with pytest.raises(ValueError, match='shape before allocation'):
        c.state_vector(toy(), np.zeros(32, complex))


def test_complete_compilation_charge_includes_trace_serialization_and_hash_return(tmp_path, monkeypatch):
    model, path, trace = toy(), tmp_path / 'model.json', tmp_path / 'trace.gz'
    k.write(path, model)
    elapsed = [0.]
    real = c.compiler.write_trace
    def delayed(*args, **kwargs):
        digest = real(*args, **kwargs)
        elapsed[0] += 7.
        return digest
    monkeypatch.setattr(c.compiler, 'write_trace', delayed)
    result = c.construct_candidate(path, k.sha(path), target_record(model, 2), .37, 10, 'HF', trace, {},
                                   clock=lambda: elapsed[0])
    assert result['compilation_CPU_seconds'] == 7.
    assert result['trace_sha256'] == k.sha(trace)


def test_censored_partial_trace_is_retained_and_cannot_be_replayed_as_complete(tmp_path, monkeypatch):
    model, path, trace = toy(), tmp_path / 'model.json', tmp_path / 'trace.gz'
    k.write(path, model)
    monkeypatch.setattr(c.compiler, 'MAX_GENERATED_GATES', 1)
    result = c.construct_candidate(path, k.sha(path), target_record(model, 2), .37, 10, 'HF', trace, {})
    assert result['status'] == 'censored' and result['compiled_resources'] is None
    assert result['partial_compiled_resources']['gates'] == 1 and trace.is_file()
    psi = np.zeros(16, complex)
    psi[3] = 1
    with pytest.raises(ValueError, match='Complete trace identity'):
        c.Replay(model, psi).replay(trace, result['identity'])
    with pytest.raises(ValueError, match='Censored trace'):
        independent.Reference(model, psi).audit_compiled_trace(trace, result['identity'], c.compiler.VERSION)


def test_trace_identity_and_HF_reference_cannot_be_silently_replaced(tmp_path):
    model, path, trace = toy(), tmp_path / 'model.json', tmp_path / 'trace.gz'
    k.write(path, model)
    result = c.construct_candidate(path, k.sha(path), target_record(model, 2), .37, 10, 'HF', trace, {'case': 'toy'})
    psi = np.zeros(16, complex)
    psi[12] = 1
    with pytest.raises(ValueError, match='HF prefix/reference'):
        c.Replay(model, psi).replay(trace, result['identity'])
    with pytest.raises(ValueError, match='HF reference'):
        independent.Reference(model, psi).audit_compiled_trace(trace, result['identity'], c.compiler.VERSION)
    psi[12], psi[3] = 0, 1
    with pytest.raises(ValueError, match='trace identity'):
        c.Replay(model, psi).replay(trace, dict(result['identity'], case='swapped'))
