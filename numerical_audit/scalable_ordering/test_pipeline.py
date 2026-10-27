"""Small independent checks for the new twelve-mode measurement pipeline."""
import importlib.util
import copy
import hashlib
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
NAMES = ('common', 'compiler', 'reference', 'training', 'compilation', 'identity', 'study', 'verify')
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
k, compiler, reference, training, compilation, identity, study, verification = (loaded[name] for name in NAMES)


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


def test_circuit_precision_and_separate_extended_reference(tmp_path):
    data = model()
    ref = reference.Reference(data)
    layout = record(data, (2, 0, 1))
    trace = tmp_path / 'precision.gz'
    result = compilation.compile_candidate(data, layout, 1.5, 30, trace, {})
    actual, _ = ref.audit_trace(trace, result['identity'], compiler.VERSION)
    assert ref.initial.dtype == np.dtype(np.complex128)
    assert actual.dtype == np.dtype(np.complex128)
    assert ref.ideal_formula(layout, result['steps'], 1.5).dtype == np.dtype(np.complex128)
    extended = ref.series_action(1.5)
    assert extended.dtype == np.dtype(np.clongdouble)
    assert np.max(abs(extended - ref.exact_action(1.5))) < 2e-12


def test_full_compilation_charge_includes_gate_object_cleanup(monkeypatch):
    clock = [0.]
    class Cleanup:
        def __del__(self):
            clock[0] += 1.
    def inner(*args):
        owned = Cleanup()
        return dict(CPU_seconds=0., wall_seconds=0.)
    monkeypatch.setattr(compilation.time, 'process_time', lambda: clock[0])
    monkeypatch.setattr(compilation, '_compile_candidate', inner)
    result = compilation.compile_candidate(None, None, None, None, None, None)
    assert result['CPU_seconds'] == 1. and result['internal_CPU_seconds'] == 0.


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


def test_cubic_validation_deadline_retains_valid_placement(monkeypatch, tmp_path):
    path = tmp_path / 'model.json'
    path.write_text(json.dumps(model()))
    ordering = k.bootstrap()
    def interrupted(*args, **kwargs):
        raise ordering.SearchDeadline('input validation reached deadline')
    monkeypatch.setattr(ordering, 'train_ordering', interrupted)
    result = training.construct(path, k.sha(path), 'cubic')
    assert result['order'] == training.current(model())
    assert result['status'] == 'cpu_limit' and result['score'] is None


def test_cubic_numerical_refusal_preserves_placement_and_specific_status(monkeypatch, tmp_path):
    path = tmp_path / 'model.json'
    path.write_text(json.dumps(model()))
    ordering = k.bootstrap()
    def unresolved(*args, **kwargs):
        raise ordering.NumericalResolution('unresolved cubic subtraction')
    monkeypatch.setattr(ordering, 'train_ordering', unresolved)
    result = training.construct(path, k.sha(path), 'cubic')
    assert result['order'] == training.current(model())
    assert result['status'] == 'numerical_limit' and result['score'] is None


def test_bootstrap_cost_is_inside_constructor_deadline(monkeypatch, tmp_path):
    path = tmp_path / 'model.json'
    path.write_text(json.dumps(model()))
    clock = [0.]
    expired_cpu = k.TRAIN_CPU_SECONDS + 10.
    def expensive_bootstrap():
        clock[0] = expired_cpu
        return None  # Trainer access would fail if the expired cap were ignored.
    monkeypatch.setattr(training.time, 'process_time', lambda: clock[0])
    monkeypatch.setattr(k, 'bootstrap', expensive_bootstrap)
    result = training.construct(path, k.sha(path), 'cubic')
    assert result['status'] == 'cpu_limit' and result['CPU_seconds'] == expired_cpu


def test_complete_construction_charge_includes_return_cleanup(monkeypatch):
    clock = [0.]
    class Cleanup:
        def __del__(self):
            clock[0] += 1.
    def inner(*args):
        owned = Cleanup()
        return dict(CPU_seconds=0., wall_seconds=0.)
    monkeypatch.setattr(training.time, 'process_time', lambda: clock[0])
    monkeypatch.setattr(training, '_construct', inner)
    result = training.construct(None, None, None)
    assert result['CPU_seconds'] == 1. and result['internal_CPU_seconds'] == 0.


def test_first_response_bounds_stored_pairs_and_retains_fallback(monkeypatch):
    monkeypatch.setattr(training, 'MAX_STORED_PAIRS', 0)
    result = training.first_search(model(), 'pauli', [0, 1, 2], float('inf'))
    assert result['status'] == 'support_limit' and result['score'] is None
    assert sorted(result['order']) == list(range(4))
    assert result['stats']['seed_evaluations'] == result['stats']['candidate_evaluations'] == 0
    assert not result['stats']['pair_setup_completed'] and result['stats']['pair_entries'] is None


def test_many_commuting_singletons_do_not_reserve_quadratic_pair_storage():
    raw = [{((0, 'Z'),): 1}] * 1001
    context = training.FirstResponse(raw, 1, 0, raw[0], float('inf'))
    assert not context.pairs and context.score(range(1001), float('inf')) == 0.


def test_first_response_uses_one_encoded_pair_dictionary():
    raw = [{((0, 'Z'),): 2}, {((0, 'X'),): 3}, {((1, 'Z'),): 4}]
    context = training.FirstResponse(raw, 1, 0, {((0, 'Z'),): 2, ((1, 'Z'),): 4}, float('inf'))
    assert isinstance(context.pairs, dict) and context.lookup is context.pairs
    assert context.pairs == {1: (1, -6.)}
    assert context.amplitude((0, 1, 2), float('inf')) == {1: -6.}
    assert context.amplitude((1, 0, 2), float('inf')) == {1: 6.}


def test_first_response_dictionary_enforces_exact_small_storage_bound(monkeypatch):
    raw = [{((0, 'Z'),): 2}, {((0, 'X'),): 3}, {((0, 'Z'), (1, 'Z')): 4}]
    monkeypatch.setattr(training, 'MAX_STORED_PAIRS', 2)
    context = training.FirstResponse(raw, 1, 0, raw[0], float('inf'))
    assert len(context.pairs) == 2 and context.lookup is context.pairs
    monkeypatch.setattr(training, 'MAX_STORED_PAIRS', 1)
    with pytest.raises(training.SupportLimit, match='stored nonzero pair bound'):
        training.FirstResponse(raw, 1, 0, raw[0], float('inf'))


def test_first_response_crossings_match_full_reconstruction_for_every_small_order():
    data = model()
    raw = [k.decode(block) for block in data['blocks']]
    diagonal = {word: coefficient for block in raw[:2] for word, coefficient in block.items()}
    context = training.FirstResponse(raw, 10, 3, diagonal, float('inf'))
    for order in itertools.permutations(range(len(raw))):
        amplitude = context.amplitude(order, float('inf'))
        score = context.score(order, float('inf'))
        for at in range(len(order) - 1):
            trial, updated = list(order), dict(amplitude)
            incremental = context.crossing(updated, score, trial[at], trial[at + 1])
            trial[at], trial[at + 1] = trial[at + 1], trial[at]
            assert updated == pytest.approx(context.amplitude(trial, float('inf')), abs=1e-15)
            assert incremental == pytest.approx(context.score(trial, float('inf')), abs=1e-15)


def test_first_search_reports_completed_seed_candidate_and_rescore_work():
    result = training.first_search(model(), 'pauli', [0, 1, 2], float('inf'))
    stats = result['stats']
    assert result['score'] is not None and stats['completed_training_score']
    assert stats['pair_setup_completed'] and stats['pair_entries'] > 0
    assert stats['seed_evaluations'] == stats['unique_starts'] == 3
    assert stats['candidate_evaluations'] > 0
    assert result['evaluations'] == stats['seed_evaluations'] + stats['candidate_evaluations']
    assert stats['full_move_rescores'] >= result['accepted_moves']


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


class BatchInput(io.BytesIO):
    def __init__(self):
        super().__init__()
        self.requests = []

    def write(self, value):
        self.requests.append(value)
        return super().write(value)


class BatchOutput(io.BytesIO):
    def __init__(self, value):
        super().__init__(value)
        self.body_reads = []
        self.header_reads = []

    def read(self, size=-1):
        assert 0 <= size <= identity.MAX_JSON_BYTES
        self.body_reads.append(size)
        return super().read(size)

    def readline(self, size=-1):
        assert size == 256
        self.header_reads.append(size)
        return super().readline(size)


class BatchProcess:
    def __init__(self, payload):
        self.stdin, self.stdout = BatchInput(), BatchOutput(payload)
        self.returncode = None
        self.killed = False

    def poll(self):
        return self.returncode

    def wait(self, timeout=None):
        assert timeout == 5
        if self.returncode is None:
            self.returncode = 0
        return self.returncode

    def kill(self):
        self.killed, self.returncode = True, -9


def batch_blob(data):
    oid = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    return oid, oid.encode() + b' blob ' + str(len(data)).encode() + b'\n' + data + b'\n'


def test_geometry_inventory_accepts_authoritative_count_with_same_byte_caps():
    inventory = {f'{i:040x}': (1, f'{i}.json') for i in range(12936)}
    identity.validate_inventory(inventory)
    with pytest.raises(ValueError, match='storage bounds'):
        identity.validate_inventory({f'{i:040x}': (1, f'{i}.json') for i in range(16385)})
    with pytest.raises(ValueError, match='storage bounds'):
        identity.validate_inventory({'0' * 40: (identity.MAX_JSON_BYTES + 1, 'oversized.json')})


def test_geometry_scan_uses_one_batch_and_authenticates_every_blob():
    first, second = b'{"first":1}', b'{"second":2}\n'
    left, one = batch_blob(first)
    right, two = batch_blob(second)
    process = BatchProcess(one + two)
    invocations = []
    def factory(command, **kwargs):
        invocations.append((command, kwargs))
        return process
    result = list(identity.stream_blobs({left: (len(first), 'first.json'), right: (len(second), 'second.json')},
                                       popen=factory))
    assert result == [('first.json', first), ('second.json', second)]
    assert len(invocations) == 1 and invocations[0][0][-2:] == ['cat-file', '--batch']
    assert process.stdin.requests == [left.encode() + b'\n', right.encode() + b'\n']
    assert process.stdout.body_reads == [len(first), 1, len(second), 1]
    assert not process.killed and process.stdin.closed and process.stdout.closed


@pytest.mark.parametrize('header', ['wrong_oid', 'wrong_kind', 'wrong_size', 'huge_size', 'long', 'missing'])
def test_invalid_git_batch_header_precedes_blob_read_and_kills_owned_process(header):
    data = b'{}'
    oid, valid = batch_blob(data)
    choices = {
        'wrong_oid': b'0' * 40 + b' blob 2\n',
        'wrong_kind': oid.encode() + b' tree 2\n',
        'wrong_size': oid.encode() + b' blob 3\n',
        'huge_size': oid.encode() + b' blob 999999999999999999\n',
        'long': b'x' * 256 + b'\n',
        'missing': oid.encode() + b' missing\n',
    }
    process = BatchProcess(choices[header] + data + b'\n')
    with pytest.raises(ValueError, match='batch header'):
        list(identity.stream_blobs({oid: (2, 'fixture.json')}, popen=lambda *a, **k: process))
    assert process.stdout.body_reads == []
    assert process.killed and process.stdin.closed and process.stdout.closed


@pytest.mark.parametrize('failure', ['truncated', 'bad_delimiter', 'wrong_hash'])
def test_git_batch_body_authentication_and_owned_process_cleanup(failure):
    data = b'{"a":1}'
    oid, valid = batch_blob(data)
    payload = valid[:-3] if failure == 'truncated' else valid[:-1] + b'x' if failure == 'bad_delimiter' else valid.replace(b'"a"', b'"b"')
    process = BatchProcess(payload)
    with pytest.raises(ValueError, match='batch blob|authentication'):
        list(identity.stream_blobs({oid: (len(data), 'fixture.json')}, popen=lambda *a, **k: process))
    assert process.killed and process.stdin.closed and process.stdout.closed


def test_git_batch_real_committed_json_blob():
    relative = 'numerical_audit/scalable_ordering/REGISTRY.json'
    raw = identity.subprocess.check_output(['git', '-C', str(k.ROOT), 'show', 'HEAD:' + relative])
    oid, _ = batch_blob(raw)
    assert list(identity.stream_blobs({oid: (len(raw), relative)})) == [(relative, raw)]

# Append to numerical_audit/scalable_ordering/test_pipeline.py. This uses its
# existing literal(), training and pytest bindings; no Git/archive dependency.


def _exhaustive_current_order(data):
    """Small integer dense oracle: rebuild every candidate's entire objective."""
    raw = [{tuple((q, p) for q, p in word): int(c) for word, c in block}
           for block in data['blocks']]
    dimension, hf = 1 << data['n'], (1 << data['electrons']) - 1
    matrices = []
    masks = []
    for block in raw:
        matrix = [[0 for _ in range(dimension)] for _ in range(dimension)]
        block_masks = []
        for word, coefficient in block.items():
            dense = literal(word, data['n'])
            assert np.max(abs(dense.imag)) == 0  # Deliberate even-Y fixtures.
            for destination in range(dimension):
                for source in range(dimension):
                    matrix[destination][source] += coefficient * int(dense[destination, source].real)
            block_masks.append(sum(1 << q for q, p in word if p in ('X', 'Y')))
        assert len(set(block_masks)) == 1
        masks.append(block_masks[0])
        matrices.append(matrix)
    diagonal = sorted((i for i, mask in enumerate(masks) if mask == 0),
                      key=lambda i: (-abs(next(iter(raw[i].values()))), next(iter(raw[i]))))
    off = sorted((i for i, mask in enumerate(masks) if mask),
                 key=lambda i: (-sum(abs(c) for c in raw[i].values()), masks[i]))
    if not diagonal or not off:
        return diagonal + off

    base = {}
    for at, i in enumerate(off):
        for j in off[at + 1:]:
            # Literal [Ai,Aj]|HF>, rather than the producer's descriptors.
            for destination in range(dimension):
                if destination == hf:
                    continue
                coefficient = sum(matrices[i][destination][middle] * matrices[j][middle][hf]
                                  - matrices[j][destination][middle] * matrices[i][middle][hf]
                                  for middle in range(dimension))
                if coefficient:
                    base[destination] = base.get(destination, 0) + coefficient
    destinations = [hf ^ masks[i] for i in off]
    rows = [[matrices[j][destination][hf]
             * (matrices[i][destination][destination] - matrices[i][hf][hf])
             for j, destination in zip(off, destinations)] for i in diagonal]
    energies = [sum(matrices[i][occupation][occupation] for i in diagonal)
                for occupation in range(dimension)]
    ratios = {destination: float(np.sinc((energies[destination] - energies[hf]) /
                                       (2 * np.pi * int(data['denominator'])))**2).as_integer_ratio()
              for destination in set(base) | set(destinations)}
    common = max(denominator for numerator, denominator in ratios.values())
    weights = {destination: numerator * (common // denominator)
               for destination, (numerator, denominator) in ratios.items()}

    def exhaustive_objective(slots):
        amplitude = dict(base)
        for slot, row in zip(slots, rows):
            for position, (destination, coefficient) in enumerate(zip(destinations, row)):
                amplitude[destination] = amplitude.get(destination, 0) + (
                    coefficient if position >= slot else -coefficient)
        return sum(weights[destination] * coefficient**2
                   for destination, coefficient in amplitude.items())

    choices = range(len(off) + 1)
    start = min(choices, key=lambda slot: exhaustive_objective([slot] * len(diagonal)))
    slots = [start] * len(diagonal)
    for coordinate in range(len(diagonal)):
        candidates = []
        for slot in choices:
            proposed = list(slots)
            proposed[coordinate] = slot
            candidates.append(exhaustive_objective(proposed))
        slots[coordinate] = min(choices, key=lambda slot: candidates[slot])
    order = []
    for slot in choices:
        order.extend(i for i, position in zip(diagonal, slots) if position == slot)
        if slot < len(off):
            order.append(off[slot])
    return order


def _current_integer_fixture(seed):
    import random
    rng = random.Random(seed)
    scale = 1 << 68 if seed % 3 == 0 else 1
    draw = lambda: scale * rng.choice((-5, -3, -1, 1, 3, 5))
    blocks = [{((q, 'Z'),): draw()} for q in range(4)]
    blocks.append({((0, 'Z'), (2, 'Z')): draw()})
    for a, b, spectator in ((0, 2, None), (1, 2, None), (0, 1, 2), (2, 3, None)):
        coefficient = draw()
        extra = () if spectator is None else ((spectator, 'Z'),)
        blocks.append({tuple(sorted(((a, 'X'), (b, 'X')) + extra)): coefficient,
                       tuple(sorted(((a, 'Y'), (b, 'Y')) + extra)): coefficient})
    # Deliberately repeated mask and exact zeros test aggregate-destination
    # updates even outside the molecular producer's distinct-mask convention.
    if seed % 2 == 0:
        blocks.append({((0, 'X'), (2, 'X')): draw()})
    rng.shuffle(blocks)
    return dict(n=4, electrons=seed % 5, denominator=str(scale * (seed + 3)),
                blocks=[[[list(word), str(c)] for word, c in block.items()] for block in blocks])


@pytest.mark.parametrize('seed', range(12))
def test_current_exact_incremental_sweeps_match_exhaustive_integer_oracle(seed):
    data = _current_integer_fixture(seed)
    assert training.current(data) == _exhaustive_current_order(data)


@pytest.mark.parametrize('blocks', [
    [{((0, 'Z'),): 3}, {((1, 'Z'),): -3}],
    [{((0, 'X'), (1, 'X')): 2, ((0, 'Y'), (1, 'Y')): 2}],
    [{((3, 'Z'),): 3},
     {((0, 'X'), (1, 'X')): 2, ((0, 'Y'), (1, 'Y')): 2}],
])
def test_current_empty_side_and_all_zero_response_keep_exact_tie_order(blocks):
    data = dict(n=4, electrons=1, denominator='13',
                blocks=[[[list(word), str(c)] for word, c in block.items()] for block in blocks])
    assert training.current(data) == _exhaustive_current_order(data)
    assert training.current(data) == list(range(len(blocks)))


def test_current_optimized_sweeps_still_enforce_constructor_deadline():
    with pytest.raises(training.Deadline, match='deadline'):
        training.current(_current_integer_fixture(5), deadline=-float('inf'))

# Append these tests to numerical_audit/scalable_ordering/test_pipeline.py.
# They use that file's isolated `compiler`, `reference`, `k`, and `model` fixtures.
# Preparation only: no tests have been executed by the helper.
from collections import Counter as _StreamCounter
from fractions import Fraction as _StreamFraction


def _stream_unbatched_standard_emitter(path, n, gates, identity):
    """Independent old-format emitter: no compiler.encode/cache/batching calls."""
    header = dict(compiler=compiler.VERSION, n=n, identity=identity,
                  scope='Ideal prototype compiler; no hardware measurement')
    with path.open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            zipped.write((json.dumps(header, sort_keys=True, separators=(',', ':')) + '\n').encode())
            for gate in gates:
                if gate.kind == 'cx':
                    record = ['cx', gate.a, gate.b]
                elif gate.kind == 'rz':
                    record = ['rz', gate.a, str(gate.angle.numerator), str(gate.angle.denominator)]
                else:
                    record = [gate.kind, gate.a]
                zipped.write((json.dumps(record, separators=(',', ':')) + '\n').encode())


def _stream_compare_emitter(tmp_path, monkeypatch, gates):
    expected, actual = tmp_path / 'unbatched.gz', tmp_path / 'batched.gz'
    identity = dict(test='byte-identical', exact='fractional angles')
    _stream_unbatched_standard_emitter(expected, 12, gates, identity)
    writes, encodes = [], []
    original_gzip, original_encode = compiler.gzip.GzipFile, compiler.encode

    class RecordingGzip:
        def __init__(self, *args, **kwargs):
            self.stream = original_gzip(*args, **kwargs)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return self.stream.__exit__(*args)

        def write(self, data):
            writes.append(len(data))
            return self.stream.write(data)

    def counted_encode(gate):
        encodes.append(gate)
        return original_encode(gate)

    monkeypatch.setattr(compiler.gzip, 'GzipFile', RecordingGzip)
    monkeypatch.setattr(compiler, 'encode', counted_encode)
    digest = compiler.write_trace(actual, 12, iter(gates), identity)
    assert actual.read_bytes() == expected.read_bytes()
    assert digest == hashlib.sha256(expected.read_bytes()).hexdigest()
    return writes, encodes


@pytest.mark.parametrize('count', [0, 8191, 8192, 8193])
def test_stream_writer_preserves_gzip_bytes_at_exact_64k_boundaries(count, tmp_path, monkeypatch):
    # ['h',0] plus newline occupies exactly eight bytes in the declared format.
    gates = [compiler.Gate('h', 0)] * count
    writes, encodes = _stream_compare_emitter(tmp_path, monkeypatch, gates)
    expected_body = [] if not count else [8 * count] if count <= 8192 else [65536, 8]
    assert writes[1:] == expected_body
    assert len(encodes) == bool(count)
    assert all(size <= 65536 for size in writes[1:])


def test_stream_writer_finite_clifford_cache_never_caches_exact_angles(tmp_path, monkeypatch):
    primitives = [compiler.Gate(kind, q) for kind in ('h', 's', 'sdg') for q in range(12)]
    primitives += [compiler.Gate('cx', q, q + 1) for q in range(11)]
    primitives += [compiler.Gate('cx', q + 1, q) for q in range(11)]
    assert len(primitives) == 58
    gates = []
    for index in range(160):
        gates.extend(primitives)
        gates.append(compiler.Gate('rz', index % 12, angle=_StreamFraction(2 * index - 159, 203)))
    writes, encodes = _stream_compare_emitter(tmp_path, monkeypatch, gates)
    assert len(encodes) == 58 + 160
    assert sum(g.kind == 'rz' for g in encodes) == 160
    assert len(writes) > 2 and all(size <= 65536 for size in writes[1:])


def _stream_hf_lines():
    return [json.dumps([kind, q], separators=(',', ':')) + '\n'
            for q in range(2) for kind in ('h', 's', 's', 'h')]


def _stream_write_manual(path, lines):
    header = dict(n=3, compiler=compiler.VERSION, identity=dict(test='decoded-line-cache'))
    with path.open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            zipped.write((json.dumps(header) + '\n').encode())
            for line in lines:
                zipped.write(line.encode())


def _stream_spy_decoder(monkeypatch):
    original_loads, original_require = reference.json.loads, reference.k.require
    calls, checks, parsed = _StreamCounter(), _StreamCounter(), []

    def counted_loads(line):
        value = original_loads(line)
        if isinstance(value, list):
            calls[line] += 1
            parsed.append((value, list(value)))
        return value

    def counted_require(condition, message):
        checks[message] += 1
        return original_require(condition, message)

    monkeypatch.setattr(reference.json, 'loads', counted_loads)
    monkeypatch.setattr(reference.k, 'require', counted_require)
    return calls, checks, parsed


def test_stream_decoder_reuses_only_short_cliffords_and_keeps_every_gate_check(tmp_path, monkeypatch):
    short = '["h",2]\n'
    long = '["h",2]' + ' ' * 200 + '\n'
    positive, negative = '["rz",2,"1","7"]\n', '["rz",2,"-1","7"]\n'
    lines = _stream_hf_lines() + [short] * 400 + [long] * 4 + [positive, negative] * 2
    path = tmp_path / 'cached.gz'
    _stream_write_manual(path, lines)
    context = reference.Reference(model())
    calls, checks, parsed = _stream_spy_decoder(monkeypatch)
    actual, counts = context.audit_trace(path, dict(test='decoded-line-cache'), compiler.VERSION, literal=True)
    assert calls[short] == 1
    assert calls[long] == 4
    assert calls[positive] == calls[negative] == 2  # Even repeated short Rz is never cached.
    assert counts['gates'] == len(lines) and counts['actual_Rz_rotations'] == 4
    assert checks['Independent gate/line bound'] == len(lines)
    assert checks['Primitive gate support'] == len(lines)
    assert checks['H payload'] == 408
    assert checks['S payload'] == 4
    assert checks['Rz exact payload'] == checks['Canonical nonzero Rz angle'] == 4
    assert all(value == saved for value, saved in parsed)  # Cached JSON lists remain untouched.
    assert abs(np.vdot(context.initial, actual)) == pytest.approx(1., abs=2e-12)


def test_stream_decoder_cache_stops_at_128_short_raw_lines(tmp_path, monkeypatch):
    # Many whitespace spellings of valid Cliffords cannot grow the cache.
    lines = _stream_hf_lines()
    for q in (2, 1):
        for padding in range(110):
            line = f'["h",{q}]' + ' ' * padding + '\n'
            assert len(line) <= 128
            lines.extend([line, line])
    path = tmp_path / 'cache-cap.gz'
    _stream_write_manual(path, lines)
    context = reference.Reference(model())
    calls, checks, parsed = _stream_spy_decoder(monkeypatch)
    _, counts = context.audit_trace(path, dict(test='decoded-line-cache'), compiler.VERSION)
    assert len(calls) > 128
    assert sum(count == 1 for count in calls.values()) == 128
    assert all(count in (1, 2) for count in calls.values())
    assert counts['gates'] == len(lines)
    assert checks['Independent gate/line bound'] == checks['Primitive gate support'] == len(lines)
    assert checks['H payload'] == len(lines) - 4
    assert all(value == saved for value, saved in parsed)


@pytest.mark.parametrize('bad,expected', [
    ('["h",2,0]\n', 'H payload'),
    ('["s",2,0]\n', 'S payload'),
    ('["cx",0,2]\n', 'Primitive line CX'),
    ('["h",true]\n', 'Primitive gate support'),
    ('["rz",2,"2","4"]\n', 'Canonical nonzero Rz angle'),
    ('["rz",2,"0","1"]\n', 'Canonical nonzero Rz angle'),
    ('["rz",2,"1","-2"]\n', 'Rz exact payload'),
    ('["h",2,"' + 'x' * 200 + '"]\n', 'H payload'),
    ('["h",2]' + ' ' * 16385 + '\n', 'gate/line bound'),
])
def test_stream_decoder_cached_prefix_does_not_bypass_malformed_payloads(bad, expected, tmp_path, monkeypatch):
    lines = _stream_hf_lines() + ['["h",2]\n'] * 4 + [bad]
    path = tmp_path / 'malformed-after-cache.gz'
    _stream_write_manual(path, lines)
    context = reference.Reference(model())
    calls, checks, _ = _stream_spy_decoder(monkeypatch)
    with pytest.raises(ValueError, match=expected):
        context.audit_trace(path, dict(test='decoded-line-cache'), compiler.VERSION)
    assert calls['["h",2]\n'] == 1
    assert checks['Independent gate/line bound'] == len(lines)
    if len(bad) > 16384:
        assert calls[bad] == 0  # Line bound precedes decoder/cache access.


@pytest.fixture
def derived_report_data(tmp_path, monkeypatch):
    """Use the maintained producer on synthetic rows; no state work or real errors."""
    cases = ('synthetic-a', 'synthetic-b', 'synthetic-c', 'synthetic-d')
    frozen = dict(models={case: 'unused-model-hash' for case in cases})
    constructors = {(case, arm): dict(CPU_seconds=(index + 1) / 10, status='complete', evaluations=index)
                    for case in cases for index, arm in enumerate(k.ARMS)}
    rows, checks = [], []
    references = [dict(case=case, T=T) for case, T in itertools.product(cases, k.TIMES)]
    for case_index, case in enumerate(cases):
        for (ti, T), cap in itertools.product(enumerate(k.TIMES), k.CAPS):
            for index, arm in enumerate(k.ARMS):
                constructor = constructors[case, arm]
                label = f'{case}.{arm}.t{ti}.c{cap}'
                censored = case_index == 3 and arm == 'pauli' and ti == 1 and cap == 30
                row = dict(candidate=label, case=case, arm=arm, T=T, cap_multiple=cap,
                    status='censored' if censored else 'complete', CPU_seconds=.2 + index / 100,
                    ordering_CPU_seconds=constructor['CPU_seconds'], ordering_status=constructor['status'],
                    ordering_evaluations=constructor['evaluations'], compiled_resources=None)
                check = dict(candidate=label, status=row['status'])
                if not censored:
                    rank = index + 1 if (case_index + ti) % 2 else len(k.ARMS) - index
                    error = 0. if (case_index, ti, cap) == (0, 0, 10) else rank * (ti + 1) * 1e-5
                    resources = dict(gates=100 + 4 * index, single_qubit_gates=80 + 3 * index,
                        two_qubit_gates=20 + index, unit_duration_ASAP_depth=50 + index,
                        CX_dependency_depth=10 + index)
                    row.update(error=error, compiled_resources=resources)
                    check.update(error=error, resources=copy.deepcopy(resources))
                rows.append(row)
                checks.append(check)
                k.write(tmp_path / f'rows/{label}.json', row)
    k.write(tmp_path / 'ALL_ORDERS_FROZEN.json', frozen)
    k.write(tmp_path / 'RUN_COMPLETE.json', dict(rows=[r['candidate'] for r in rows],
        references={f"{r['case']}.{r['T']}": {} for r in references}))
    monkeypatch.setattr(study, 'bind', lambda output: None)
    summary = study.summarize(tmp_path)
    return dict(frozen=frozen, rows=rows, constructors=constructors, checks=checks,
                references=references, saved_pairs=k.read(tmp_path / 'PAIRS.json'), summary=summary)


def test_derived_gate_accepts_real_producer_arithmetic_without_propagation(derived_report_data, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Derived-data arithmetic must not propagate states')
    monkeypatch.setattr(verification, 'Reference', forbidden)
    monkeypatch.setattr(verification, 'metrics', forbidden)
    report, pairs = verification.validate_derived_data(**derived_report_data)
    assert report['status'] == 'passed'
    assert (report['circuits'], report['references'], report['all_pair_count']) == (96, 8, 240)
    assert len(pairs) == 240
    assert sum(r['independent_error_decision'] == 'censored' for r in pairs) == 5
    assert sum(r['independent_error_decision'] == 'tie' for r in pairs) >= 15
    assert set(report['checked_summary_fields']) >= {'statuses', 'cubic_error_comparisons', 'circuits'}


@pytest.mark.parametrize('field', [
    'error_left', 'error_right', 'error_difference', 'CX_left', 'CX_right',
    'depth_left', 'depth_right', 'construction_CPU_left', 'construction_CPU_right',
])
def test_derived_gate_rejects_tampered_saved_pair_fields(field, derived_report_data):
    pair = next(r for r in derived_report_data['saved_pairs'] if r['status'] == 'complete')
    pair[field] += 1
    with pytest.raises(ValueError, match='Independent saved pair field: ' + field):
        verification.validate_derived_data(**derived_report_data)


@pytest.mark.parametrize('field', ['circuits', 'references', 'all_pair_count'])
def test_derived_gate_rejects_tampered_summary_counts(field, derived_report_data):
    derived_report_data['summary'][field] += 1
    with pytest.raises(ValueError, match='Independent summary count: ' + field):
        verification.validate_derived_data(**derived_report_data)


def test_derived_gate_rejects_tampered_summary_status_counts(derived_report_data):
    derived_report_data['summary']['statuses']['complete'] += 1
    with pytest.raises(ValueError, match='Independent summary statuses'):
        verification.validate_derived_data(**derived_report_data)


def test_derived_gate_rejects_tampered_summary_cubic_aggregates(derived_report_data):
    group = derived_report_data['summary']['cubic_error_comparisons']['canonical']
    group[next(iter(group))] += 1
    with pytest.raises(ValueError, match='Independent summary cubic comparisons'):
        verification.validate_derived_data(**derived_report_data)


def test_derived_gate_checks_ordering_cost_against_authenticated_constructor(derived_report_data):
    key = next(iter(derived_report_data['constructors']))
    derived_report_data['constructors'][key]['CPU_seconds'] += 1
    with pytest.raises(ValueError, match='Independent row ordering cost/status'):
        verification.validate_derived_data(**derived_report_data)


def test_derived_gate_checks_pair_resources_against_independent_dag(derived_report_data):
    row = next(r for r in derived_report_data['rows'] if r['status'] == 'complete')
    row['compiled_resources']['two_qubit_gates'] += 1
    with pytest.raises(ValueError, match='Independent verified resource fields'):
        verification.validate_derived_data(**derived_report_data)


def test_derived_gate_rejects_metrics_attached_to_censored_pair(derived_report_data):
    pair = next(r for r in derived_report_data['saved_pairs'] if r['status'] == 'censored')
    pair['CX_left'] = 123
    with pytest.raises(ValueError, match='Censored pair must omit completed metrics'):
        verification.validate_derived_data(**derived_report_data)


def test_complete_small_study_runs_summarizes_and_independently_verifies(tmp_path, monkeypatch):
    """Exercise real phase composition; synthetic chemistry is not an SCF audit."""
    output, destination = tmp_path / 'study', tmp_path / 'independent'
    output.mkdir()
    models, constructors = {}, {}
    for index in range(4):
        case = f'synthetic-three-qubit-{index}'
        data = dict(model(), case=case)
        path = output / f'chemistry/inputs/{case}/model.json'
        k.write(path, data)
        models[case] = k.sha(path)
        for arm in k.ARMS:
            record = training.construct(path, models[case], arm)
            relative = f'constructors/{case}.{arm}.json'
            k.write(output / relative, record)
            constructors[relative] = k.sha(output / relative)
    frozen = dict(source_commit='synthetic integration fixture', sources=k.source_map(),
        environment=k.environment(), models=models, constructors=constructors,
        arms=list(k.ARMS), times=list(k.TIMES), caps=list(k.CAPS),
        train_time=k.TRAIN_TIME, train_CPU_seconds=k.TRAIN_CPU_SECONDS,
        expected_circuits=96, expected_references=8)
    k.write(output / 'ALL_ORDERS_FROZEN.json', frozen)
    # Input-generation authentication has separate coverage. Everything after
    # this boundary, including verify's source/model/constructor hashes, is real.
    monkeypatch.setattr(study, 'bind', lambda supplied_output: frozen)
    study.run(output)
    summary = study.summarize(output)
    independent = verification.verify(output, destination)
    complete = k.read(output / 'RUN_COMPLETE.json')
    assert len(constructors) == 24
    assert (len(complete['rows']), len(complete['references'])) == (96, 8)
    assert (summary['circuits'], summary['references'], summary['all_pair_count']) == (96, 8, 240)
    assert summary['statuses'] == {'complete': 96}
    assert independent['status'] == 'passed'
    assert (len(independent['candidates']), len(independent['references']), len(independent['pairs'])) == (96, 8, 240)
    assert independent['derived_data_validation']['status'] == 'passed'
    assert k.read(destination / 'INDEPENDENT_VERIFIED.json') == independent
