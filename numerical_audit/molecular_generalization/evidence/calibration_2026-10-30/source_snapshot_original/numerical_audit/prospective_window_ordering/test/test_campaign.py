"""The complete order freeze must precede every per-case physical index."""
import copy
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
    assert plan['attempt_tag'] == 'a1'
    assert len(plan['entries']) == plan['expected_calls'] == 64
    assert [row['arm'] for row in plan['entries'][:8]] == list(k.ARMS)
    assert [row['arm'] for row in plan['entries'][8:16]] == list(reversed(k.ARMS))
    assert str(executable) in plan['entries'][0]['command']
    assert str(executable.resolve()) not in plan['entries'][0]['command']
    assert len({row['label'] for row in plan['entries']}) == 64
    assert all(row['label'].endswith('-a1') and row['result'].endswith(row['arm'] + '.a1.json')
               for row in plan['entries'])
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
    assert frozen['attempt_tag'] == 'a1'
    assert frozen['inventory_complete'] and not frozen['all_calls_have_artifacts']
    assert len(frozen['artifacts']) == 64 and all(row['record'] is None for row in frozen['artifacts'])
    assert frozen['retained_temporary_and_checkpoint_artifacts'][0]['sha256'] == k.sha(checkpoint)
    indices = list(tmp_path.glob('*.PANEL_INVENTORY.json'))
    assert len(indices) == 4
    for path in indices:
        value = json.loads(path.read_text())
        assert value['attempt_tag'] == 'a1'
        assert path.name == 'main-' + value['case'] + '.PANEL_INVENTORY.json'
        assert value['all_orders_freeze']['sha256'] == k.sha(output)
        assert len(value['constructors']) == 16
        assert all(row['selected_attempt'] is None for row in value['constructors'])
        assert all('No constructor artifact' in row['censor_reason'] for row in value['constructors'])


def preflight_prepared(root):
    """Tiny exact source fixture; no SCF, numerical constructor or physics."""
    original = prepared(root)
    parent = next(iter(k.read(original / 'INPUTS_FROZEN.json')['successes']))
    source_model = k.read(original / 'inputs' / parent / 'model.json')
    case = parent + '_preflight_7over5'
    model = copy.deepcopy(source_model)
    model['case'] = case
    model['denominator'] = str(5 * int(source_model['denominator']))
    for word in model['coefficients']:
        word[1] = str(7 * int(word[1]))
    for block in model['blocks']:
        for word in block:
            word[1] = str(7 * int(word[1]))
    data = root / 'altered-data'
    target = data / 'inputs' / case
    target.mkdir(parents=True)
    k.write_new(target / 'model.json', model)
    k.write_model_toml(target / 'model.toml', model, source_json_sha256=k.sha(target / 'model.json'))
    successes = {case: {name: k.sha(target / name) for name in ('model.json', 'model.toml')}}
    k.write_new(data / 'INPUTS_FROZEN.json', dict(status='complete', expected_cases=1,
        preflight_only=True, source_protocol_path=str(original / 'INPUT_PROTOCOL_FROZEN.json'),
        successes=successes, refusals={}))
    return data


def retained_receipts(root, plan):
    for row in plan['entries']:
        path = root / row['receipt']
        path.parent.mkdir(exist_ok=True)
        k.write_new(path, dict(command=row['command'], exit_code=1))


@pytest.mark.parametrize('mode', ['main', 'preflight'])
def test_second_attempt_tag_separates_full_schedule_and_preserves_a1(tmp_path, mode):
    data = prepared(tmp_path) if mode == 'main' else preflight_prepared(tmp_path)
    first_path = tmp_path / 'PLAN.a1.json'
    first = c.plan(tmp_path, data, mode, first_path, Path(sys.executable), Path('/usr/bin/julia'))
    retained_receipts(tmp_path, first)
    old_result = tmp_path / first['entries'][0]['result']
    k.write_new(old_result, dict(retained_original='a1 fixture'))
    whole1 = tmp_path / 'ALL_ORDERS.a1.json'
    c.freeze(tmp_path, first_path, whole1)
    originals = {path: path.read_bytes() for path in [first_path, whole1, old_result,
        *tmp_path.glob('*.PANEL_INVENTORY.json')]}
    second_path = tmp_path / 'PLAN.a2.json'
    second = c.plan(tmp_path, data, mode, second_path, Path(sys.executable), Path('/usr/bin/julia'),
                    attempt_tag='a2')
    assert second['attempt_tag'] == 'a2'
    assert len(first['entries']) == len(second['entries']) == (64 if mode == 'main' else 8)
    assert second['source_snapshot'] == first['source_snapshot']
    assert second['source_freeze'] == first['source_freeze']
    assert second['input_freeze'] == first['input_freeze']
    for old, new in zip(first['entries'], second['entries']):
        for field in ('case', 'arm', 'repetition', 'cwd', 'attempt_CPU_seconds', 'attempt_wall_seconds'):
            assert new[field] == old[field]
        assert new['label'] == old['label'][:-2] + 'a2'
        assert new['result'] == old['result'][:-7] + 'a2.json'
        assert new['receipt'] == 'logs/' + new['label'] + '.resources.json'
        assert new['command'] == [str(tmp_path / new['result']) if arg == str(tmp_path / old['result']) else arg
                                  for arg in old['command']]
        assert not (tmp_path / new['result']).exists()
    retained_receipts(tmp_path, second)
    whole2 = tmp_path / 'ALL_ORDERS.a2.json'
    c.freeze(tmp_path, second_path, whole2, attempt_tag='a2')
    frozen = k.read(whole2)
    assert frozen['attempt_tag'] == 'a2' and frozen['plan_sha256'] == k.sha(second_path)
    for case in {row['case'] for row in second['entries']}:
        path = tmp_path / (mode + '-a2-' + case + '.PANEL_INVENTORY.json')
        value = k.read(path)
        assert value['attempt_tag'] == 'a2'
        assert value['all_orders_freeze'] == dict(path=whole2.name, sha256=k.sha(whole2))
        assert len(value['constructors']) == (16 if mode == 'main' else 8)
        assert all(row['attempts'][0]['receipt']['path'].endswith('-a2.resources.json')
                   for row in value['constructors'])
    assert all(path.read_bytes() == raw for path, raw in originals.items())


def test_attempt_tag_has_tight_ascii_format_and_cli_forwarding(tmp_path):
    data = preflight_prepared(tmp_path)
    for tag in ('a1', 'a2', 'a999'):
        assert c.checked_attempt_tag(tag) == tag
    for i, tag in enumerate(('', 'a0', 'a00', 'a01', 'a1000', 'A2', 'a-2', 'a2/x',
                             ' a2', 'a2\n', 'a２', None, 2, True)):
        output = tmp_path / ('invalid-' + str(i) + '.json')
        with pytest.raises(ValueError, match='Attempt tag'):
            c.plan(tmp_path, data, 'preflight', output, Path(sys.executable), Path('/usr/bin/julia'),
                   attempt_tag=tag)
        assert not output.exists()
    output = tmp_path / 'CLI.a2.json'
    assert c.main(['plan', '--root', str(tmp_path), '--inputs', str(data), '--mode', 'preflight',
        '--output', str(output), '--python', sys.executable, '--julia', '/usr/bin/julia',
        '--attempt-tag', 'a2']) == 0
    assert k.read(output)['attempt_tag'] == 'a2'
    assert all(row['label'].endswith('-a2') for row in k.read(output)['entries'])


@pytest.mark.parametrize('tag', ['a1', 'a2'])
@pytest.mark.parametrize('collision', ['result', 'inventory', 'partial_result'])
def test_plan_rejects_existing_target_namespace_without_overwrite(tmp_path, tag, collision):
    data = prepared(tmp_path)
    case = k.read(data / 'INPUT_PROTOCOL_FROZEN.json')['registry']['cases'][0]['case']
    if collision == 'inventory':
        prefix = 'main' if tag == 'a1' else 'main-' + tag
        retained = tmp_path / (prefix + '-' + case + '.PANEL_INVENTORY.json')
    else:
        retained = tmp_path / 'constructors-main' / case / 'r1' / ('canonical.' + tag + '.json')
        if collision == 'partial_result':
            retained = retained.with_name(retained.name + '.partial')
    retained.parent.mkdir(parents=True, exist_ok=True)
    retained.write_bytes(b'original retained artifact\n')
    output = tmp_path / 'NEW_PLAN.json'
    with pytest.raises(ValueError, match='collision'):
        c.plan(tmp_path, data, 'main', output, Path(sys.executable), Path('/usr/bin/julia'), attempt_tag=tag)
    assert retained.read_bytes() == b'original retained artifact\n'
    assert not output.exists() and not output.with_name(output.name + '.partial').exists()


def test_freeze_prechecks_all_index_collisions_before_whole_freeze_publication(tmp_path):
    data = prepared(tmp_path)
    plan_path = tmp_path / 'PLAN.a2.json'
    plan = c.plan(tmp_path, data, 'main', plan_path, Path(sys.executable), Path('/usr/bin/julia'), attempt_tag='a2')
    retained_receipts(tmp_path, plan)
    # Collide with the last case too: no earlier per-case index or whole freeze
    # may be published before the full index namespace is admitted.
    case = plan['entries'][-1]['case']
    retained = tmp_path / ('main-a2-' + case + '.PANEL_INVENTORY.json')
    retained.write_bytes(b'old index must survive\n')
    output = tmp_path / 'ALL_ORDERS.a2.json'
    with pytest.raises(ValueError, match='inventory.*collision'):
        c.freeze(tmp_path, plan_path, output)
    assert not output.exists()
    assert list(tmp_path.glob('*.PANEL_INVENTORY.json')) == [retained]
    assert retained.read_bytes() == b'old index must survive\n'


def test_freeze_tag_linkage_cli_assertion_and_legacy_a1_default(tmp_path):
    data = preflight_prepared(tmp_path)
    plan_path = tmp_path / 'PLAN.a2.json'
    plan = c.plan(tmp_path, data, 'preflight', plan_path, Path(sys.executable), Path('/usr/bin/julia'), attempt_tag='a2')
    retained_receipts(tmp_path, plan)
    output = tmp_path / 'ALL_ORDERS.a2.json'
    with pytest.raises(ValueError, match='Attempt tag differs'):
        c.freeze(tmp_path, plan_path, output, attempt_tag='a3')
    assert not output.exists()
    damaged = copy.deepcopy(plan)
    damaged['attempt_tag'] = 'a3'
    damaged_path = tmp_path / 'WRONG_TAG_PLAN.json'
    k.write_new(damaged_path, damaged)
    with pytest.raises(ValueError, match='tag/path linkage'):
        c.freeze(tmp_path, damaged_path, output)
    assert not output.exists() and not list(tmp_path.glob('*.PANEL_INVENTORY.json'))
    assert c.main(['freeze', '--root', str(tmp_path), '--plan', str(plan_path),
        '--output', str(output), '--attempt-tag', 'a2']) == 0
    assert k.read(output)['attempt_tag'] == 'a2'
    assert k.read(output)['plan_sha256'] == k.sha(plan_path)
    # Missing tag in an otherwise current-source plan keeps the old a1 schema
    # interpretation; it cannot relabel a2 entries as a1.
    legacy_path = tmp_path / 'LEGACY.a1.json'
    legacy = c.plan(tmp_path, data, 'preflight', legacy_path, Path(sys.executable), Path('/usr/bin/julia'))
    del legacy['attempt_tag']
    # This new fixture path avoids overwriting even the synthetic plan.
    missing_path = tmp_path / 'LEGACY_WITHOUT_TAG.json'
    k.write_new(missing_path, legacy)
    retained_receipts(tmp_path, legacy)
    legacy_output = tmp_path / 'ALL_ORDERS.legacy.json'
    c.freeze(tmp_path, missing_path, legacy_output)
    assert k.read(legacy_output)['attempt_tag'] == 'a1'
    case = legacy['entries'][0]['case']
    assert (tmp_path / ('preflight-' + case + '.PANEL_INVENTORY.json')).is_file()
