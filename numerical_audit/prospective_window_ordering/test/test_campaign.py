"""The complete order freeze must precede every per-case physical index."""
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import study_common as k
import construct_campaign as c


def prepared(root):
    data = root / 'data'
    data.mkdir()
    registry = k.read(HERE / 'REGISTRY.json')
    successes = {}
    for spec in registry['cases']:
        target = data / 'inputs' / spec['case']
        target.mkdir(parents=True)
        words = [[[[0, 'Z']], '1'], [[[13, 'Z']], '-3']]
        model = dict(spec, denominator='2', coefficients=words, blocks=[[words[0]], [words[1]]])
        k.write_new(target / 'model.json', model)
        k.write_model_toml(target / 'model.toml', model, source_json_sha256=k.sha(target / 'model.json'))
        successes[spec['case']] = {name: k.sha(target / name) for name in ('model.json', 'model.toml')}
    k.write_new(data / 'INPUT_PROTOCOL_FROZEN.json', dict(sources=k.source_map(), registry=registry,
                                                      configuration=k.configuration()))
    k.write_new(data / 'INPUTS_FROZEN.json', dict(status='complete', expected_cases=4,
                                               successes=successes, refusals={}))
    return data


def test_plan_preserves_venv_path_and_balances_fresh_repetitions(tmp_path):
    data = prepared(tmp_path)
    # Resolving this symlink would silently lose its environment selection.
    executable = tmp_path / 'venv-python'
    executable.symlink_to('/usr/bin/python3')
    output = tmp_path / 'PLAN.json'
    plan = c.plan(tmp_path, data, 'main', output, executable, Path('/usr/bin/julia'))
    assert len(plan['entries']) == plan['expected_calls'] == 64
    assert [row['arm'] for row in plan['entries'][:8]] == list(k.ARMS)
    assert [row['arm'] for row in plan['entries'][8:16]] == list(reversed(k.ARMS))
    assert str(executable) in plan['entries'][0]['command']
    assert str(executable.resolve()) not in plan['entries'][0]['command']
    assert len({row['label'] for row in plan['entries']}) == 64
    with pytest.raises(ValueError, match='unperformed'):
        c.freeze(tmp_path, output, tmp_path / 'ALL_ORDERS.json')
    assert not (tmp_path / 'ALL_ORDERS.json').exists()
    assert not list(tmp_path.glob('*.PANEL_INVENTORY.json'))


def test_only_whole_call_inventory_can_publish_per_case_indices(tmp_path):
    data = prepared(tmp_path)
    plan_path = tmp_path / 'PLAN.json'
    plan = c.plan(tmp_path, data, 'main', plan_path, Path(sys.executable), Path('/usr/bin/julia'))
    for row in plan['entries']:
        path = tmp_path / row['receipt']
        path.parent.mkdir(exist_ok=True)
        # A killed job has a receipt, but no invented order/record.
        k.write_new(path, dict(command=row['command'], exit_code=1))
    output = tmp_path / 'ALL_ORDERS.json'
    checkpoint = tmp_path / 'tmp' / 'interrupted-worker' / 'checked-order.json'
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_text('{"checked_order": [0, 1]}')
    c.freeze(tmp_path, plan_path, output)
    frozen = json.loads(output.read_text())
    assert frozen['inventory_complete'] and not frozen['all_calls_have_artifacts']
    assert len(frozen['artifacts']) == 64 and all(row['record'] is None for row in frozen['artifacts'])
    assert frozen['retained_temporary_and_checkpoint_artifacts'][0]['sha256'] == k.sha(checkpoint)
    indices = list(tmp_path.glob('*.PANEL_INVENTORY.json'))
    assert len(indices) == 4
    for path in indices:
        value = json.loads(path.read_text())
        assert value['all_orders_freeze']['sha256'] == k.sha(output)
        assert len(value['constructors']) == 16
        assert all(row['selected_attempt'] is None for row in value['constructors'])
        assert all('No constructor artifact' in row['censor_reason'] for row in value['constructors'])
