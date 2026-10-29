"""Prepare fresh-process job lists and freeze their entire order inventory.

This module never launches work. The root coordinator dispatches each planned
command through the cumulative resource guard, then calls freeze before physics.
Planning defaults to attempt tag a1; --attempt-tag a2 separates continued calls
without replacing old artifacts. Freeze derives that tag from its linked plan.
"""
from __future__ import annotations

import argparse
import re
import study_common as k


def descriptor(path, root):
    path = path.resolve()
    k.require(path.is_relative_to(root), 'Campaign artifact must remain within its root')
    return dict(path=str(path.relative_to(root)), sha256=k.sha(path))


def checked_attempt_tag(value):
    k.require(isinstance(value, str) and re.fullmatch(r'a[1-9][0-9]{0,2}', value) is not None,
              'Attempt tag must match a[1-9][0-9]{0,2}')
    return value


def inventory_path(root, mode, case, attempt_tag):
    prefix = mode if attempt_tag == 'a1' else mode + '-' + attempt_tag
    return root / (prefix + '-' + case + '.PANEL_INVENTORY.json')


def call_paths(root, mode, case, arm, repetition, attempt_tag):
    label = ('window-' + mode + '-' + case.replace('_', '-') + '-' + arm.replace('_', '-') +
             '-r' + str(repetition) + '-' + attempt_tag)
    result = root / ('constructors-' + mode) / case / ('r' + str(repetition)) / (arm + '.' + attempt_tag + '.json')
    return label, result


def require_fresh(path, name):
    partial = path.with_name(path.name + '.partial')
    k.require(not path.exists() and not path.is_symlink() and
              not partial.exists() and not partial.is_symlink(), name + ' collision; preserve existing artifact')


def plan(root, inputs, mode, output, python, julia, attempt_tag='a1'):
    root, inputs = root.resolve(), inputs.resolve()
    attempt_tag = checked_attempt_tag(attempt_tag)
    k.require(mode in ('main', 'preflight'), 'Known mode')
    require_fresh(output, 'Plan output')
    frozen_path = inputs / 'INPUTS_FROZEN.json'
    frozen = k.read(frozen_path)
    k.require(frozen['status'] == 'complete' and not frozen['refusals'], 'Complete declared input inventory')
    source = (inputs / 'INPUT_PROTOCOL_FROZEN.json' if mode == 'main' else
              k.Path(frozen['source_protocol_path']))
    protocol = k.read(source)
    k.require(protocol['sources'] == k.source_map(), 'Scientific source differs from source freeze')
    k.require(protocol.get('configuration') == k.configuration(), 'Frozen dispatch configuration differs')
    if mode == 'main':
        names = [row['case'] for row in protocol['registry']['cases']]
        k.require(len(names) == 4 and set(frozen['successes']) == set(names) and
                  frozen['expected_cases'] == 4, 'Original complete four-case roster')
    else:
        names = sorted(frozen['successes'])
        k.require(len(names) == 1 and frozen['expected_cases'] == 1 and frozen['preflight_only'],
                  'Distinct one-case altered preflight')
    entries = []
    for case in names:
        require_fresh(inventory_path(root, mode, case, attempt_tag), 'Per-case inventory')
    for case in names:
        model = inputs / 'inputs' / case / 'model.json'
        k.require(k.sha(model) == frozen['successes'][case]['model.json'], 'Frozen JSON differs')
        for repetition in range(1, (k.REPETITIONS if mode == 'main' else 1) + 1):
            arms = k.ARMS if repetition == 1 else tuple(reversed(k.ARMS))
            for arm in arms:
                label, result = call_paths(root, mode, case, arm, repetition, attempt_tag)
                require_fresh(result, 'Planned result')
                command = ['/usr/bin/env', 'TMPDIR=' + str(root / 'tmp'), 'JULIA_NUM_THREADS=1',
                    'JULIA_NUM_GC_THREADS=1', 'JULIA_PKG_OFFLINE=true', str(python.absolute()),
                    str(k.HERE / 'constructors.py'), str(model), k.sha(model), arm, str(result),
                    '--repetition', str(repetition), '--julia-executable', str(julia.absolute()),
                    '--source-freeze', str(source)]
                entries.append(dict(case=case, arm=arm, repetition=repetition, label=label,
                    command=command, cwd=str(k.ROOT), attempt_CPU_seconds=k.TRAIN_CPU_SECONDS,
                    attempt_wall_seconds=k.TRAIN_WALL_SECONDS,
                    result=str(result.relative_to(root)),
                    receipt='logs/' + label + '.resources.json'))
    expected = 64 if mode == 'main' else 8
    k.require(len(entries) == expected, 'Complete fresh-process call roster')
    value = dict(schema='prospective-window-constructor-plan-v1', mode=mode, attempt_tag=attempt_tag, entries=entries,
        root=str(root), inputs=str(inputs), source_freeze=descriptor(source, root),
        input_freeze=descriptor(frozen_path, root), expected_calls=expected,
        source_snapshot=protocol['sources'], physical_error_read=False)
    if mode == 'preflight':
        value['parent_input_freeze'] = descriptor(source.parent / 'INPUTS_FROZEN.json', root)
    # Admit the complete new namespace before creating directories/publishing a
    # plan; no retained result or per-case inventory can be replaced.
    for entry in entries:
        (root / entry['result']).parent.mkdir(parents=True, exist_ok=True)
    k.write_new(output, value)
    return value


def freeze(root, plan_path, output, *, attempt_tag=None):
    root = root.resolve()
    data = k.read(plan_path)
    tag = checked_attempt_tag(data.get('attempt_tag', 'a1'))
    if attempt_tag is not None:
        k.require(checked_attempt_tag(attempt_tag) == tag, 'Attempt tag differs from planned inventory')
    require_fresh(output, 'Whole-inventory output')
    k.require(data['source_snapshot'] == k.source_map(), 'Source changed after constructor planning')
    k.require(data['expected_calls'] == (64 if data['mode'] == 'main' else 8) and
              len(data['entries']) == data['expected_calls'], 'Complete planned roster required')
    for case in {entry['case'] for entry in data['entries']}:
        require_fresh(inventory_path(root, data['mode'], case, tag), 'Per-case inventory')
    records = []
    indices = {}
    for entry in data['entries']:
        label, planned_result = call_paths(root, data['mode'], entry['case'], entry['arm'], entry['repetition'], tag)
        k.require(entry['label'] == label and entry['result'] == str(planned_result.relative_to(root)) and
                  entry['receipt'] == 'logs/' + label + '.resources.json', 'Attempt tag/path linkage differs from plan')
        receipt = root / entry['receipt']
        k.require(receipt.is_file(), 'Attempt unperformed or unreconciled; no whole-inventory freeze: ' + entry['label'])
        parsed = k.read(receipt)
        k.require(parsed['command'] == entry['command'], 'Executed command differs from planned source/input/arm')
        record = root / entry['result']
        attempt = dict(receipt=descriptor(receipt, root))
        attempt['record'] = descriptor(record, root) if record.is_file() else None
        records.append(dict(case=entry['case'], arm=entry['arm'], repetition=entry['repetition'],
            receipt=attempt['receipt'], record=attempt['record'], exit_code=parsed['exit_code']))
        if entry['case'] not in indices:
            model = k.Path(data['inputs']) / 'inputs' / entry['case'] / 'model.json'
            index = dict(schema='prospective-window-panel-inventory-v1', mode=data['mode'], attempt_tag=tag,
                case=entry['case'], source_freeze=data['source_freeze'], input_freeze=data['input_freeze'],
                model_json=descriptor(model, root), model_toml=descriptor(model.with_suffix('.toml'), root),
                constructors=[])
            if data['mode'] == 'preflight':
                index['parent_input_freeze'] = data['parent_input_freeze']
            indices[entry['case']] = index
        selected = dict(arm=entry['arm'], repetition=entry['repetition'],
            selected_attempt=1 if record.is_file() else None, attempts=[attempt])
        if not record.is_file():
            selected['censor_reason'] = ('No constructor artifact after retained attempt; exit_code=' +
                str(parsed['exit_code']) + '; termination=' + str(parsed.get('termination')))
        indices[entry['case']]['constructors'].append(selected)
    temporary = root / 'tmp'
    auxiliary = sorted(path for path in temporary.rglob('*') if path.is_file()) if temporary.exists() else []
    k.require(len(auxiliary) <= 512 and sum(path.stat().st_size for path in auxiliary) <= 256 * 1024**2,
              'Bounded complete temporary/checkpoint inventory required')
    checkpoints = [dict(**descriptor(path, root), bytes=path.stat().st_size) for path in auxiliary]
    # First freeze ALL calls across ALL cases, before publishing per-case indices.
    k.write_new(output, dict(schema='prospective-window-all-orders-freeze-v1',
        mode=data['mode'], attempt_tag=tag, source_freeze=data['source_freeze'], plan_sha256=k.sha(plan_path),
        inventory_complete=True, expected_calls=data['expected_calls'], artifacts=records,
        retained_temporary_and_checkpoint_artifacts=checkpoints,
        all_calls_have_artifacts=all(row['record'] is not None for row in records),
        physical_error_read=False, no_selected_repetition=True))
    for case, index in indices.items():
        index['all_orders_freeze'] = descriptor(output, root)
        k.write_new(inventory_path(root, data['mode'], case, tag), index)
    print('ALL_ORDERS_FROZEN', data['mode'], len(records), flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('plan', 'freeze'))
    parser.add_argument('--root', type=k.Path, required=True)
    parser.add_argument('--output', type=k.Path, required=True)
    parser.add_argument('--inputs', type=k.Path)
    parser.add_argument('--mode', choices=('main', 'preflight'))
    parser.add_argument('--python', type=k.Path)
    parser.add_argument('--julia', type=k.Path)
    parser.add_argument('--plan', type=k.Path)
    parser.add_argument('--attempt-tag', help='Bounded planning tag (default a1); freeze derives its plan tag unless asserted here')
    args = parser.parse_args(argv)
    if args.action == 'plan':
        k.require(all(value is not None for value in (args.inputs, args.mode, args.python, args.julia)), 'Complete planning arguments')
        plan(args.root, args.inputs, args.mode, args.output, args.python, args.julia,
             attempt_tag='a1' if args.attempt_tag is None else args.attempt_tag)
    else:
        k.require(args.plan is not None, 'Planned inventory required')
        freeze(args.root, args.plan, args.output, attempt_tag=args.attempt_tag)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
