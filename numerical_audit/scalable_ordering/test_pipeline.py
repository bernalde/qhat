"""Small independent checks for the new twelve-mode measurement pipeline."""
import importlib.util
import itertools
import json
import io
import gzip
from pathlib import Path
import sys
import zipfile
import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
NAMES = ('common', 'compiler', 'reference', 'training', 'compilation')
saved = {name: sys.modules.get(name) for name in NAMES}
loaded = {}
try:
    for name in NAMES:
        spec = importlib.util.spec_from_file_location('_scalable_pipeline_test_' + name, HERE / (name + '.py'))
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        sys.modules[name] = module
        spec.loader.exec_module(module)
        loaded[name] = module
finally:
    for name, old in saved.items():
        if old is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = old
k, compiler, reference, training, compilation = (loaded[name] for name in NAMES)


def model(n=3):
    raw = [([[0, 'Z']], '5'), ([[1, 'Z']], '-3'), ([[0, 'X'], [n - 1, 'X']], '2'),
           ([[0, 'Y'], [n - 1, 'Y']], '2')]
    return dict(n=n, electrons=2, denominator='10', coefficients=raw,
                blocks=[[raw[0]], [raw[1]], [raw[2], raw[3]]])


def record(data, order):
    blocks = [k.decode(b) for b in data['blocks']]
    return dict(partition='conserving', formula=2, order=list(order),
        unit_lengths=[len(blocks[i]) for i in order],
        schedule=[[[[q, p] for q, p in word], str(c)] for i in order for word, c in sorted(blocks[i].items())])


@pytest.mark.parametrize('order', list(itertools.permutations(range(3))))
def test_actual_primitive_frame_matches_literal_gates_and_ideal(order, tmp_path):
    data = model()
    layout = record(data, order)
    path = tmp_path / 'trace.gz'
    compiled = compilation.compile_candidate(data, layout, .5, 10, path, dict(test='frame'))
    ref = reference.Reference(data)
    actual, counts = ref.audit_trace(path, compiled['identity'], compiler.VERSION, literal=True)
    ideal = ref.ideal_formula(layout, compiled['steps'], .5)
    overlap = np.vdot(ideal, actual)
    assert np.max(abs(actual - ideal * overlap / abs(overlap))) < 2e-12
    counts.pop('actual_Rz_rotations')
    assert counts == compiled['compiled_resources']
    assert reference.metrics(ideal, actual, data)['error'] < 1e-20


def test_twelve_mode_register_and_long_line_routing(tmp_path):
    data = model(12)
    layout = record(data, (2, 0, 1))
    path = tmp_path / 'long-line.gz'
    compiled = compilation.compile_candidate(data, layout, .5, 10, path, dict(test='long-line'))
    ref = reference.Reference(data)
    actual, counts = ref.audit_trace(path, compiled['identity'], compiler.VERSION)
    ideal = ref.ideal_formula(layout, compiled['steps'], .5)
    assert abs(abs(np.vdot(ideal, actual)) - 1) < 2e-11
    assert counts['two_qubit_gates'] > 0
    assert ref.dimension == 4096


def test_reference_refuses_large_register_before_allocating(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('allocation must not be reached')
    monkeypatch.setattr(reference.np, 'arange', forbidden)
    with pytest.raises(ValueError, match='register'):
        reference.Reference(dict(n=13))


def test_taylor_reference_matches_exponential():
    ref = reference.Reference(model())
    assert np.max(abs(ref.series_action(1.5) - ref.exact_action(1.5))) < 2e-12


def literal(word, n):
    local = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
             'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    out = np.array([[1.]])
    labels = dict(word)
    for q in reversed(range(n)):
        out = np.kron(out, local[labels.get(q, 'I')])
    return out


def test_first_response_matches_independent_dense_commutators():
    data = model()
    raw = [k.decode(block) for block in data['blocks']]
    diagonal = {w: c for block in raw[:2] for w, c in block.items()}
    context = training.FirstResponse(raw, 10, 3, diagonal, float('inf'))
    matrices = [sum(c / 10 * literal(w, 3) for w, c in block.items()) for block in raw]
    D = (matrices[0] + matrices[1]).diagonal().real
    psi = np.eye(8)[:, 3]
    for order in itertools.permutations(range(3)):
        column = np.zeros(8, complex)
        for at, i in enumerate(order):
            for j in order[at + 1:]:
                column += -.5j * (matrices[i] @ matrices[j] - matrices[j] @ matrices[i]) @ psi
        gaps = D - D[3]
        weights = np.sinc(gaps / (2 * np.pi))**2
        weights[3] = 0
        expected = float(np.sum(weights * abs(column)**2))
        assert context.score(order, float('inf')) == pytest.approx(expected, rel=2e-12, abs=1e-18)


def test_current_baseline_retains_frozen_legacy_order():
    folder = HERE.parent / 'molecular_guard/evidence/results'
    data = json.loads((folder / 'model.json').read_text())
    expected = json.loads((folder / 'FROZEN.json').read_text())['orders']['current']
    assert training.current(data) == expected


def test_duplicate_terms_are_rejected():
    with pytest.raises(ValueError, match='Unique'):
        k.decode([[[[0, 'Z']], '1'], [[[0, 'Z']], '2']])


def test_pauli_timeout_keeps_correct_partition(monkeypatch, tmp_path):
    path = tmp_path / 'model.json'
    path.write_text(json.dumps(model()))
    def stopped(*args):
        raise training.Deadline('test deadline')
    monkeypatch.setattr(training, 'current', stopped)
    result = training.construct(path, k.sha(path), 'pauli')
    assert result['partition'] == 'pauli'
    assert result['order'] == list(range(len(model()['coefficients'])))
    assert result['status'] == 'cpu_limit'
    assert len(result['schedule']) == len(model()['coefficients'])


def test_npy_shape_is_rejected_before_numpy_allocation(monkeypatch, tmp_path):
    stream = io.BytesIO()
    np.lib.format.write_array_header_1_0(stream, dict(descr='<c16', fortran_order=False, shape=(1_000_000_000,)))
    path = tmp_path / 'malformed.npz'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('actual.npy', stream.getvalue())
    def forbidden(*args, **kwargs):
        raise AssertionError('np.load must not be called')
    monkeypatch.setattr(np, 'load', forbidden)
    with pytest.raises(ValueError, match='shape'):
        k.load_vector(path, 'actual', 3)


def test_npz_extra_member_is_rejected_before_numpy_load(monkeypatch, tmp_path):
    path = tmp_path / 'extra.npz'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('actual.npy', b'not an NPY')
        archive.writestr('extra.npy', b'not an NPY')
    def forbidden(*args, **kwargs):
        raise AssertionError('np.load must not be called')
    monkeypatch.setattr(np, 'load', forbidden)
    with pytest.raises(ValueError, match='member'):
        k.load_vector(path, 'actual', 3)


def test_bounded_npz_reads_actual_saved_state(tmp_path):
    path = tmp_path / 'valid.npz'
    state = np.arange(8, dtype=float).astype(np.clongdouble)
    np.savez_compressed(path, actual=state)
    assert np.array_equal(k.load_vector(path, 'actual', 3), state)


@pytest.mark.parametrize('where', ['header', 'gate'])
def test_gzip_line_bound_precedes_json_decode(where, tmp_path):
    path = tmp_path / 'oversized.gz'
    header = dict(n=3, compiler=compiler.VERSION, identity={})
    with gzip.open(path, 'wt') as stream:
        if where == 'gate':
            stream.write(json.dumps(header) + '\n')
        stream.write(' ' * 16385 + '\n')
    with pytest.raises(ValueError, match='line bound'):
        reference.Reference(model()).audit_trace(path, {}, compiler.VERSION)


def test_nested_commutator_oracle_matches_dense_cubic_taylor():
    data = model()
    matrices = [sum(int(c) / 10 * literal(w, 3) for w, c in block) for block in data['blocks']]
    H = sum(matrices)
    psi = np.eye(8)[:, 3]
    ref = reference.Reference(data)
    for order in itertools.permutations(range(3)):
        q = [psi.copy()] + [np.zeros(8, complex) for _ in range(3)]
        sequence = [(order[0], .5), (order[1], .5), (order[2], 1.), (order[1], .5), (order[0], .5)]
        for index, tau in sequence:
            out = [v.copy() for v in q]
            power = np.eye(8)
            factor = 1.
            for degree in range(1, 4):
                power = matrices[index] @ power
                factor *= tau / degree
                for kdegree in range(4 - degree):
                    out[kdegree + degree] += factor * power @ q[kdegree]
            q = out
        column = q[3] - H @ H @ H @ psi / 6
        gaps = (matrices[0] + matrices[1]).diagonal().real
        weights = np.sinc((gaps - gaps[3]) / (2 * np.pi))**2
        weights[3] = 0
        layout = record(data, order)
        expected = np.dot(weights, abs(column)**2) / compilation.affordable_steps(layout['unit_lengths'], 30)**4
        actual = ref.objective(layout)['value']
        assert actual == pytest.approx(expected, rel=1e-10, abs=1e-20)
