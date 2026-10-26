"""Construction accounting, exact fusion and cost-constrained selection guards."""
import copy
import gc
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pytest
import common
import construction as c
import study


def toy():
    z0 = [[[0, 'Z']], '3']
    z1 = [[[1, 'Z']], '-1']
    xx = [[[0, 'X'], [1, 'X']], '2']
    yy = [[[0, 'Y'], [1, 'Y']], '2']
    return dict(case='asymmetric_toy', n=2, electrons=1, denominator='7',
                blocks=[[z0], [z1], [xx, yy]], coefficients=[z0, z1, xx, yy])


def record(model):
    return dict(arm='first_JD', formula=2, T=common.TIME, partition='conserving',
                order=[1, 2, 0], unit_lengths=[1, 2, 1],
                schedule=[term for i in [1, 2, 0] for term in sorted(model['blocks'][i])])


@pytest.mark.parametrize('steps', [1, 3, 6, 24, 128])
def test_fine_grid_fusion_preserves_each_unit_duration(steps):
    fused = list(c.fused_units([2, 1, 3], steps))
    assert all(a[0] != b[0] for a, b in zip(fused, fused[1:]))
    for unit in range(3):
        assert sum((tau for i, tau in fused if i == unit), Fraction(0)) == 1
    # Independent count: 5 per symmetric step, one boundary fusion per join.
    assert len(fused) == 4 * steps + 1


def test_exact_compilation_and_unfused_literal_action_agree(tmp_path):
    model = toy()
    rec = record(model)
    path = tmp_path / 'model.json'
    common.write(path, model)
    trace = tmp_path / 'trace.gz'
    result = c.construct_candidate(path, common.sha(path), rec, 3, trace,
                                   {'test': 'noncoarse-grid'}, .1, 8.)
    assert result['eligible'] and result['status'] == 'complete'
    ref = common.reference.Reference(model)
    actual, costs = ref.audit_compiled_trace(trace, result['identity'], common.compiler.VERSION)
    ideal = ref.ideal_formula(rec, 3, common.TIME)
    assert np.max(abs(actual - ideal)) < 2e-14
    assert costs == result['compiled_resources']


def test_total_cap_counts_outer_constructor_return(monkeypatch, tmp_path):
    now = [0.]
    def compile_inner(*args):
        assert args[-2] == pytest.approx(1.8)
        now[0] += 1.81
        return {'status': 'complete'}
    monkeypatch.setattr(c, '_compile', compile_inner)
    result = c.construct_candidate(None, None, None, 3, None, None, .2, 2., clock=lambda: now[0])
    assert result['total_construction_CPU_seconds'] == pytest.approx(2.01)
    assert not result['eligible'] and result['CPU_overrun_seconds'] == pytest.approx(.01)


def test_context_cleanup_is_charged(monkeypatch):
    now = [0.]
    class Context:
        def __init__(self, *args):
            now[0] += .1
        def score(self, order, *args):
            now[0] += .01
            return dict(order=order, score=1.)
        def __del__(self):
            now[0] += .3
    def search(context, seeds, *args):
        now[0] += .2
        return dict(record=seeds[0])
    monkeypatch.setattr(c.cubic, 'FirstOrderContext', Context)
    monkeypatch.setattr(c.cubic, 'insertion', search)
    phases = {}
    result = c._search({'coefficients': [1]}, 'conserving', [[0], [0], [0]], 2., phases, 'test', lambda: now[0])
    gc.collect()
    assert result['best']['order'] == [0]
    assert phases['test_cleanup_CPU_seconds'] == pytest.approx(.3)
    assert sum(phases.values()) == pytest.approx(.63)


@pytest.mark.parametrize('method', common.METHODS)
def test_full_order_constructor_returns_exact_complete_partition(tmp_path, method):
    path = tmp_path / 'model.json'
    common.write(path, toy())
    result = c.construct_order(path, common.sha(path), method, 2.)
    c.validate_record(toy(), result['record'])
    assert result['ordering_CPU_seconds'] >= sum(result['phases'].values()) - 1e-12
    if method == 'cubic_warm_JD':
        assert result['warm_start'] is not None
        assert any(k.startswith('warm_') for k in result['phases'])
    assert result['record']['partition'] == ('pauli' if method == 'pauli_JD' else 'conserving')


def test_compiler_remainder_exhaustion_preserves_partial_trace(tmp_path):
    path = tmp_path / 'model.json'
    common.write(path, toy())
    trace = tmp_path / 'partial.gz'
    result = c.construct_candidate(path, common.sha(path), record(toy()), 3, trace, {}, 2., 2.)
    assert result['status'] == 'censored' and not result['eligible']
    assert trace.exists() and result['trace_sha256'] == common.sha(trace)
    assert result['partial_compiled_resources']['gates'] == 0


def candidate(name, steps, cx, depth, error, eligible=True):
    row = dict(candidate=name, steps=steps, compiled_resources=dict(two_qubit_gates=cx,
               unit_duration_ASAP_depth=depth, gates=100), error=error)
    return dict(candidate=name, status='complete', eligible=eligible), row


def test_selection_rejects_ineligible_and_ambiguous_cheap_candidates():
    pairs = [candidate('overrun', 1, 1, 1, 1e-10, False),
             candidate('ambiguous', 2, 2, 2, 1e-6 + 1e-13),
             candidate('deep', 3, 20, 80, 1e-7), candidate('shallow', 4, 20, 40, 2e-7)]
    panel = [c for c, row in pairs]
    rows = {row['candidate']: row for c, row in pairs}
    selected = study.select(panel, rows, 1e-6)
    assert selected['selected_candidate'] == 'shallow'
    assert selected['construction_ineligible_candidates'] == 1
    assert selected['ambiguous_candidates'] == ['ambiguous']
    assert selected['qualifying_candidates'] == 2


def test_unmet_ambiguous_and_censored_are_retained():
    first, row = candidate('ambiguous', 3, 10, 10, 1e-6)
    panel = [first, dict(candidate='partial', status='censored', eligible=False)]
    selected = study.select(panel, {'ambiguous': row}, 1e-6)
    assert selected['status'] == 'ambiguous' and selected['censored_candidates'] == 1
    assert study.select(panel, {'ambiguous': row}, 1e-8)['status'] == 'unmet'


def test_unequal_model_layout_cannot_silently_flatten_wrong_partition():
    model = toy()
    rec = record(model)
    broken = copy.deepcopy(rec)
    broken['unit_lengths'] = [2, 1, 1]
    with pytest.raises(ValueError, match='layout'):
        c.validate_record(model, broken)
    oversized = copy.deepcopy(model)
    oversized['n'] = 30
    with pytest.raises(ValueError, match='Register bound'):
        c.validate_record(oversized, rec)


def test_balanced_dispatch_keeps_every_repetition_and_budget():
    for case in range(3):
        cells = study.dispatch(case, 3)
        assert len(cells) == len(set(cells)) == 24
        assert set(cells) == {(m, b, r) for m in common.METHODS for b in common.BUDGETS for r in range(3)}
    assert study.dispatch(0, 3)[0][0] != study.dispatch(1, 3)[0][0]


def test_independent_target_decisions_reject_overruns(tmp_path):
    import verify
    curve = dict(case='toy', method='cubic_JD', cohort='matched', total_budget=2., repetition=0)
    curves = {'toy': curve}
    candidates, errors, costs = {}, {}, {}
    for step in common.GRID:
        label = 'toy.r' + format(step, '03d')
        candidates[label] = dict(candidate=label, status='complete', eligible=step != 1, steps=step)
        errors[label] = {'error': 1e-10 if step == 1 else 1e-6 if step == 2 else 1e-7}
        costs[label] = dict(two_qubit_gates=step, unit_duration_ASAP_depth=step * 2, gates=step * 3)
    result = verify.independent_choices(curves, candidates, errors, costs)
    assert result[1]['selected_candidate'] == 'toy.r003'
    assert result[1]['ambiguous_candidates'] == ['toy.r002']
    assert result[2]['status'] == 'unmet'


@pytest.mark.parametrize('change', ['missing_binary', 'changed_binary'])
def test_delivery_authentication_catches_omitted_and_changed_binaries(tmp_path, change):
    import delivery
    root = tmp_path / 'delivery'
    root.mkdir()
    binary = root / 'state.npz'
    binary.write_bytes(b'owned binary fixture')
    common.write(root / 'MANIFEST.json', dict(files={'state.npz': common.sha(binary)}))
    if change == 'missing_binary':
        binary.unlink()
    else:
        binary.write_bytes(b'altered fixture')
    with pytest.raises(ValueError, match='inventory|bytes'):
        delivery.check(root)


@pytest.mark.parametrize('method', common.METHODS)
def test_warmup_exercises_complete_order_and_compile_paths(tmp_path, method):
    path = tmp_path / 'model.json'
    common.write(path, toy())
    common.write(tmp_path / 'PROTOCOL_FROZEN.json', {'test': 'full constructor warmup'})
    warm = study.warmup_constructors(path, common.sha(path), method, tmp_path, 'toy')
    result = warm['warmup_compilation']
    assert result['status'] == 'complete'
    assert result['ordering_CPU_charge'] == warm['ordering_CPU_seconds']
    ref = common.reference.Reference(toy())
    state, cost = ref.audit_compiled_trace(tmp_path / result['trace_file'], result['identity'], common.compiler.VERSION)
    assert cost == result['compiled_resources']
    assert np.max(abs(state - ref.ideal_formula(warm['record'], 3, common.TIME))) < 2e-14


@pytest.mark.parametrize('change', ['location', 'python_release', 'numpy_version'])
def test_independent_replay_accepts_relocation_but_rejects_version_change(monkeypatch, change):
    # Real committed payload: preserve every frozen source/data check while
    # simulating a pinned CI interpreter in a different path/build location.
    output = common.HERE / 'evidence/results'
    observed = copy.deepcopy(common.read(output / 'PROTOCOL_FROZEN.json')['environment'])
    observed['executable'] = '/usr/local/bin/python3.12'
    observed['python'] = observed['python'].split()[0] + ' (CI build) [GCC test]'
    if change == 'python_release':
        observed['python'] = '3.12.99 (CI build) [GCC test]'
    elif change == 'numpy_version':
        observed['packages']['numpy'] = '0.0.0'
    monkeypatch.setattr(study, 'environment', lambda: observed)
    if change == 'location':
        protocol, frozen = study.bind(output, independent=True)
        assert protocol['source_commit'] == 'ae70c96d11572a990c1a591fd92f6c0e0c05c49d'
        assert frozen['complete'] == 1281
    else:
        with pytest.raises(ValueError, match='Numerical interpreter/package identity'):
            study.bind(output, independent=True)
