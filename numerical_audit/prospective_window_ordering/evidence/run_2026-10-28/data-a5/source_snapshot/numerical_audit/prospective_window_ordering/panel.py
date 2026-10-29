"""Freeze a physical candidate matrix from authenticated constructor artifacts.

Only bounded serialization/authentication lives here; no optimizer, propagation,
chemistry import or subprocess is invoked. Usage: panel.py INVENTORY.json NEW_DIR.
See assemble() and the source-only coordination handoff for inventory schema.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import tomllib

import study_common as k

SCHEMA = 'prospective-window-panel-inventory-v1'
MAX_JSON = 32 * 1024**2
HEX = re.compile(r'^[0-9a-f]{64}$')
STATUSES = {'complete', 'cpu_limit', 'evaluation_limit', 'support_limit',
            'numerical_limit', 'control_or_backend_refusal', 'unexpected_error'}
KINDS = {'retained_first_order_response', 'diagonal_reference_cubic_response'}
PROFILE = dict(CPUQuota='25%', CPUWeight='5', IOWeight='10', Nice='19',
               MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0',
               MemoryAccounting='yes', CPUAccounting='yes')
ROUNDING_CPU_SECONDS = .03  # GNU time's decimal receipt, never a cap extension.
ALLOCATION_TOLERANCE = 1e-9  # Remainder comparisons only; actual ceilings stay strict.
THREAD_KEYS = {'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'JULIA_NUM_THREADS',
               'JULIA_NUM_GC_THREADS', 'JULIA_IMAGE_THREADS'}
COMMAND_ENV_KEYS = THREAD_KEYS | {'TMPDIR', 'JULIA_PKG_OFFLINE'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value, name, *, positive=False):
    require(type(value) in (int, float), name + ': numeric value required')
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(name + ': numeric conversion overflows') from error
    require(math.isfinite(result) and (result > 0 if positive else result >= 0),
            name + ': finite nonnegative value required')
    return result


def integer(value, name, low=0, high=10**9):
    require(type(value) is int and low <= value <= high, name + ': integer domain')
    return value


def unique_object(items):
    result = {}
    for name, value in items:
        require(name not in result, 'Duplicate JSON object key: ' + name)
        result[name] = value
    return result


def decode_json(raw):
    def reject_constant(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    def finite_float(value):
        result=float(value)
        require(math.isfinite(result), 'Overflowing JSON float')
        return result
    result = json.loads(raw, object_pairs_hook=unique_object, parse_constant=reject_constant,
                        parse_float=finite_float)
    require(isinstance(result, dict), 'Top-level JSON artifact must be an object')
    return result


def raw_file(path, maximum=MAX_JSON):
    require(path.is_file() and path.stat().st_size <= maximum, 'Bounded artifact required: ' + str(path))
    with path.open('rb') as stream:
        raw = stream.read(maximum + 1)
    require(len(raw) <= maximum, 'Artifact grew beyond byte bound')
    return raw


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def relative_path(root, spelling):
    require(isinstance(spelling, str) and 0 < len(spelling) <= 1024,
            'Bounded relative artifact path required')
    relative = PurePosixPath(spelling)
    require(not relative.is_absolute() and '..' not in relative.parts and '\\' not in spelling,
            'Artifact path escape refused')
    target = root.joinpath(*relative.parts)
    require(target.resolve().is_relative_to(root.resolve()), 'Artifact symlink escape refused')
    return target


def artifact(root, descriptor, *, maximum=MAX_JSON, parse=True):
    require(isinstance(descriptor, dict) and set(descriptor) == {'path', 'sha256'},
            'Artifact descriptor must contain exactly path/sha256')
    require(isinstance(descriptor['sha256'], str) and HEX.fullmatch(descriptor['sha256']),
            'Artifact SHA-256 spelling')
    path = relative_path(root, descriptor['path'])
    raw = raw_file(path, maximum)
    require(digest(raw) == descriptor['sha256'], 'Artifact SHA-256 differs: ' + descriptor['path'])
    return path, decode_json(raw) if parse else raw


def permutation(order, size):
    require(isinstance(order, list) and len(order) == size and
            all(type(index) is int for index in order) and sorted(order) == list(range(size)),
            'Complete original zero-based order required')
    return order


def bounded_coefficient(value):
    require(type(value) in (str, int) and not isinstance(value, bool), 'Exact integer coefficient required')
    if isinstance(value, str):
        require(len(value) <= 1280 and re.fullmatch(r'[+-]?[0-9]+', value), 'Bounded integer spelling')
    result = int(value)
    require(result.bit_length() <= 4096, 'Coefficient integer bit bound')
    return result


def authenticate_model(model, registry, *, parent_case=None):
    require(isinstance(model, dict), 'Model object required')
    require(isinstance(registry, dict) and isinstance(registry.get('cases'), list) and
            all(isinstance(spec, dict) for spec in registry['cases']), 'Registry case objects required')
    case = model['case'] if parent_case is None else parent_case
    specs = [spec for spec in registry['cases'] if spec['case'] == case]
    require(len(specs) == 1, 'Unique frozen registry identity required')
    for key, value in specs[0].items():
        if key != 'case' or parent_case is None:
            require(model.get(key) == value, 'Model/registry identity differs: ' + key)
    for key in ('basis', 'unit', 'charge', 'spin'):
        require(model.get(key) == registry[key], 'Model/global registry identity differs: ' + key)
    require(type(model['n']) is int and model['n'] == 14 and
            type(model['electrons']) is int and model['electrons'] == 4,
            'Original fourteen-mode/four-electron target')
    denominator = bounded_coefficient(model['denominator'])
    require(denominator > 0, 'Positive original denominator required')
    units = k.blocks(model)  # Existing exact serialization contract; no kernels.
    return units


def authenticate_toml(raw, model, model_hash, units):
    require(len(raw) <= k.MAX_MODEL_BYTES, 'Model TOML byte bound')
    data = tomllib.loads(raw.decode('utf-8'))
    for key in ('case', 'n', 'electrons'):
        require(data.get(key) == model[key], 'JSON/TOML model identity differs: ' + key)
    require(data.get('partition') == 'conserving' and data.get('source_json_sha256') == model_hash,
            'TOML must bind exact source JSON and original conserving partition')
    require(isinstance(data.get('denominator'), str) and
            bounded_coefficient(data['denominator']) == bounded_coefficient(model['denominator']),
            'TOML original common denominator differs')
    raw_units = data['units']
    require(isinstance(raw_units, list) and len(raw_units) == len(units), 'Original unit-label cardinality differs')
    denominator = bounded_coefficient(model['denominator'])
    for raw_unit, unit in zip(raw_units, units):
        require(isinstance(raw_unit, dict), 'TOML unit object required')
        expected = sorted(unit.items())
        require(raw_unit.get('words') == [','.join(str(site) + letter for site, letter in word)
                                        for word, _ in expected],
                'TOML original words/labels differ from structured JSON ordering')
        numerators = raw_unit.get('numerators')
        require(isinstance(numerators, list) and len(numerators) == len(expected) and
                all(isinstance(value, str) for value in numerators), 'Exact TOML numerator strings required')
        require([bounded_coefficient(value) for value in numerators] ==
                [int(coefficient * denominator) for _, coefficient in expected],
                'TOML coefficient target differs from exact JSON integers')
    baseline = data.get('baseline_order', data.get('cubic_order'))
    require(isinstance(baseline, list) and all(type(index) is int for index in baseline) and
            sorted(baseline) == list(range(1, len(units) + 1)), 'TOML baseline original labels')


def source_binding(root, index):
    source_path, frozen = artifact(root, index['source_freeze'])
    require(isinstance(frozen, dict) and isinstance(frozen.get('sources'), dict), 'Source freeze map required')
    require(re.fullmatch(r'[0-9a-f]{40}', frozen.get('source_commit', '')), 'Source freeze commit identity')
    require(frozen['sources'] == k.source_map(), 'Scientific sources differ from exact source freeze')
    registry_raw = raw_file(k.HERE / 'REGISTRY.json', k.MAX_MODEL_BYTES)
    registry = decode_json(registry_raw)
    require(frozen.get('registry') == registry and len(registry['cases']) == 4 and
            len({spec['case'] for spec in registry['cases']}) == 4,
            'Original four-case source/registry freeze required')
    require(frozen['sources'].get('numerical_audit/prospective_window_ordering/REGISTRY.json') == digest(registry_raw),
            'Registry bytes missing from source freeze')
    return source_path, frozen, registry, digest(registry_raw)


def generation_protocol(root, index, frozen_input, source_freeze, registry):
    expected = frozen_input.get('protocol_sha256')
    require(isinstance(expected, str) and HEX.fullmatch(expected), 'Input freeze generation-protocol hash required')
    descriptor = index.get('input_protocol', index['source_freeze'])
    require(descriptor['sha256'] == expected, 'Explicit input_protocol descriptor required for distinct generation freeze')
    _, protocol = artifact(root, descriptor)
    require(protocol.get('sources') == source_freeze['sources'] and protocol.get('registry') == registry,
            'Input generation/source/registry identities differ')


def frozen_case(root, freeze_path, frozen, name, registry):
    require(isinstance(frozen.get('successes'), dict), 'Frozen input success map required')
    hashes = frozen['successes'][name]
    require(isinstance(hashes, dict) and all(key in hashes for key in ('model.json', 'model.toml')),
            'Frozen JSON/TOML input identities required')
    # This path convention is the one authored by generate_inputs.py.
    parent = freeze_path.parent
    model_path = relative_path(root, str((parent / 'inputs' / name / 'model.json').relative_to(root)))
    model_raw = raw_file(model_path, k.MAX_MODEL_BYTES)
    require(digest(model_raw) == hashes['model.json'], 'Frozen original model JSON hash differs')
    model = decode_json(model_raw)
    require(model.get('case') == name, 'Frozen inventory key/model case differs')
    units = authenticate_model(model, registry)
    toml_path = relative_path(root, str((parent / 'inputs' / name / 'model.toml').relative_to(root)))
    toml_raw = raw_file(toml_path, k.MAX_MODEL_BYTES)
    require(digest(toml_raw) == hashes['model.toml'], 'Frozen original model TOML hash differs')
    authenticate_toml(toml_raw, model, hashes['model.json'], units)
    for filename, field in (('tensors.npz', 'tensors_sha256'), ('orbitals.npz', 'orbitals_sha256')):
        require(hashes.get(filename) == model.get(field) and isinstance(hashes.get(filename), str) and
                HEX.fullmatch(hashes[filename]), 'Frozen tensor/orbital provenance differs')
        path = relative_path(root, str((parent / 'inputs' / name / filename).relative_to(root)))
        require(digest(raw_file(path)) == hashes[filename], 'Frozen tensor/orbital bytes differ')
    return model, hashes


def input_binding(root, index, source_freeze, registry):
    input_path, frozen = artifact(root, index['input_freeze'])
    require(isinstance(frozen.get('successes'), dict), 'Input success map required')
    require(frozen.get('status') == 'complete' and not frozen.get('refusals'), 'Complete pre-physical input freeze required')
    model_path, model = artifact(root, index['model_json'], maximum=k.MAX_MODEL_BYTES)
    toml_path, toml_raw = artifact(root, index['model_toml'], maximum=k.MAX_MODEL_BYTES, parse=False)
    case = model['case'];require(case == index['case'], 'Panel/model case identity')
    require(isinstance(case, str) and re.fullmatch(r'[A-Za-z0-9_.-]{1,200}', case), 'Portable bounded case identity')
    expected_cases = {spec['case'] for spec in registry['cases']}
    if index['mode'] == 'main':
        require(frozen.get('expected_cases') == 4 and set(frozen.get('successes', {})) == expected_cases and
                not frozen.get('preflight_only', False), 'Complete original four-case input freeze required')
        generation_protocol(root, index, frozen, source_freeze, registry)
        require(case in expected_cases, 'Main panel may not substitute a registry case')
        # Authenticate the complete input inventory, not just a selected success.
        original = {name: frozen_case(root, input_path, frozen, name, registry) for name in sorted(expected_cases)}
        units = authenticate_model(model, registry)
    else:
        require(frozen.get('expected_cases') == 1 and frozen.get('preflight_only') is True and
                set(frozen.get('successes', {})) == {case}, 'Distinct one-case altered-input preflight freeze required')
        parent_path, parent = artifact(root, index['parent_input_freeze'])
        require(isinstance(parent.get('successes'), dict), 'Parent input success map required')
        links = [parent_hash for field in ('parent_inventory_sha256', 'parent_input_freeze_sha256')
                 if (parent_hash := frozen.get(field)) is not None]
        require(links and all(value == index['parent_input_freeze']['sha256'] for value in links),
                'Preflight parent inventory hash differs')
        require(parent.get('status') == 'complete' and not parent.get('refusals') and
                parent.get('expected_cases') == 4 and set(parent.get('successes', {})) == expected_cases,
                'Preflight needs the complete original four-case inventory')
        generation_protocol(root, index, parent, source_freeze, registry)
        original = {name: frozen_case(root, parent_path, parent, name, registry) for name in sorted(expected_cases)}
        chosen = min(expected_cases, key=lambda name: (-len(original[name][0]['coefficients']), name))
        parent_model, parent_hashes = original[chosen]
        require(case == chosen + '_preflight_7over5' and model.get('preflight_parent_case') == chosen and
                model.get('preflight_parent_model_sha256') == parent_hashes['model.json'] and
                model.get('preflight_scale') == '7/5', 'Frozen largest-case/7-over-5 preflight identity differs')
        units = authenticate_model(model, registry, parent_case=chosen)
        require(model.get('tensors_sha256') == parent_hashes['tensors.npz'] and
                model.get('orbitals_sha256') == parent_hashes['orbitals.npz'],
                'Preflight tensor/orbital identities must match the authenticated parent')
        require(bounded_coefficient(model['denominator']) == 5 * bounded_coefficient(parent_model['denominator']),
                'Preflight denominator must be exactly parent denominator times five')
        for field in ('coefficients', 'blocks'):
            raw_parent = [parent_model[field]] if field == 'coefficients' else parent_model[field]
            raw_altered = [model[field]] if field == 'coefficients' else model[field]
            require(len(raw_parent) == len(raw_altered), 'Preflight original block cardinality differs')
            for before, after in zip(raw_parent, raw_altered):
                a, b = k.decode(before), k.decode(after)
                require([entry[0] for entry in before] == [entry[0] for entry in after],
                        'Preflight must preserve exact within-block word ordering')
                require(set(a) == set(b) and all(b[word] == 7 * value for word, value in a.items()),
                        'Preflight must scale every exact original numerator by seven without thinning/regrouping')
    hashes = frozen['successes'][case]
    require(isinstance(hashes, dict), 'Selected input hash object required')
    require(hashes.get('model.json') == index['model_json']['sha256'] and
            hashes.get('model.toml') == index['model_toml']['sha256'], 'Selected input artifacts differ from frozen inventory')
    authenticate_toml(toml_raw, model, index['model_json']['sha256'], units)
    models = ({name: dict(model=value[0], model_hash=value[1]['model.json'], units=k.blocks(value[0]),
                         model_path=relative_path(root, str((input_path.parent/'inputs'/name/'model.json').relative_to(root))))
               for name, value in original.items()} if index['mode'] == 'main' else {})
    models[case] = dict(model=model, model_hash=index['model_json']['sha256'], units=units, model_path=model_path)
    return model_path, toml_path, model, units, models


def execution_path(value, execution_root, root):
    require(isinstance(value, str) and len(value) <= 4096, 'Bounded argv artifact path')
    path = PurePosixPath(value)
    require(path.is_absolute() and '..' not in path.parts, 'Constructor argv artifact paths must be absolute')
    try:
        relative = path.relative_to(execution_root)
    except ValueError as error:
        raise ValueError('Constructor argv path is outside declared execution_root') from error
    return relative_path(root, str(relative))


def receipt_binding(root, receipt, *, arm, repetition, model_path, model_hash, source_path,
                    execution_root, execution_source_root, seen_receipts):
    require(isinstance(receipt, dict), 'Guard receipt object required')
    unit = receipt.get('unit')
    require(isinstance(unit, str) and unit and unit not in seen_receipts, 'Guard unit reused or missing')
    seen_receipts.add(unit)
    settings = receipt.get('settings', {})
    require(isinstance(settings, dict), 'Guard settings object required')
    require(all(settings.get(key) == value for key, value in PROFILE.items()), 'Mandatory guard profile differs')
    budget = receipt['budget']
    require(isinstance(budget, dict), 'Guard budget object required')
    cpu_cap = number(budget['cpu_seconds'], 'guard CPU cap', positive=True)
    wall_cap = number(budget['wall_seconds'], 'guard wall cap', positive=True)
    require(cpu_cap <= k.TRAIN_CPU_SECONDS and wall_cap <= k.TRAIN_WALL_SECONDS and
            number(float(settings['RuntimeMaxSec']), 'RuntimeMaxSec', positive=True) == wall_cap,
            'Guard attempt may not extend the declared construction ceiling')
    command = receipt['command']
    require(isinstance(command, list) and all(isinstance(token, str) and len(token) <= 4096 for token in command),
            'Structured bounded guard argv required')
    script = str(execution_source_root / 'numerical_audit/prospective_window_ordering/constructors.py')
    positions = [index for index, token in enumerate(command) if token == script]
    require(len(positions) == 1, 'Guard must directly execute the frozen constructors.py')
    prefix=command[:positions[0]]
    require(prefix and re.fullmatch(r'python(?:[0-9]+(?:\.[0-9]+)*)?', PurePosixPath(prefix[-1]).name),
            'Frozen constructor must be the direct Python script, not a wrapper/shell argument')
    overrides={}
    if len(prefix)>1:
        require(PurePosixPath(prefix[0]).name=='env', 'Only a literal env-assignment prefix is admitted')
        for setting in prefix[1:-1]:
            require(re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*=.*',setting), 'Literal environment assignment required')
            key,value=setting.split('=',1)
            require(key in COMMAND_ENV_KEYS and key not in overrides, 'Unknown/repeated command environment override')
            overrides[key]=value
    base_environment=receipt.get('environment', {})
    require(isinstance(base_environment, dict), 'Guard environment object required')
    environment=dict(base_environment);environment.update(overrides)
    require(all(environment.get(key) == '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                      'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS')),
            'Effective numerical thread profile differs')
    require(all(environment[key]=='1' for key in ('JULIA_NUM_THREADS','JULIA_NUM_GC_THREADS','JULIA_IMAGE_THREADS')
                if key in environment), 'Effective Julia thread profile differs')
    if 'JULIA_PKG_OFFLINE' in environment:
        require(environment['JULIA_PKG_OFFLINE'] in ('true','1'), 'Effective Julia offline profile differs')
    for key in ('TMPDIR','NUMBA_CACHE_DIR','MPLCONFIGDIR'):
        if key in environment:
            execution_path(environment[key],execution_root,root)
    tail = command[positions[0] + 1:]
    require(len(tail) >= 8, 'Constructor argv is incomplete')
    require(execution_path(tail[0], execution_root, root).resolve() == model_path.resolve() and
            tail[1] == model_hash and tail[2] == arm, 'Guard/model/arm identity differs')
    output = execution_path(tail[3], execution_root, root)
    options = tail[4:]
    require(len(options) % 2 == 0, 'Constructor options must use explicit value arguments')
    parsed = {}
    for name, value in zip(options[::2], options[1::2]):
        require(name in ('--repetition', '--source-freeze', '--julia-executable') and name not in parsed,
                'Unknown or repeated constructor option')
        parsed[name] = value
    require(parsed.get('--repetition') == str(repetition) and '--source-freeze' in parsed,
            'Guard repetition/source-freeze binding differs')
    require(execution_path(parsed['--source-freeze'], execution_root, root).resolve() == source_path.resolve(),
            'Guard source freeze PATH differs from authenticated source freeze')
    accounting = receipt['final_accounting']
    require(isinstance(accounting, dict), 'Guard final accounting object required')
    complete = (accounting.get('CPU_accounting_available') is True and
                accounting.get('process_CPU_timing_available') is True and
                receipt.get('guard_error') is None and not receipt.get('cleanup_errors') and
                type(receipt.get('exit_code')) is int and type(receipt.get('supervisor_exit_code')) is int)
    cpu = None if accounting.get('CPU_seconds') is None else number(accounting['CPU_seconds'], 'guard total CPU')
    wall = None if accounting.get('elapsed_seconds') is None else number(accounting['elapsed_seconds'], 'guard absolute wall')
    require(wall is None or abs(wall - number(receipt['elapsed'], 'receipt elapsed')) <= 1e-9,
            'Guard absolute wall fields differ')
    complete = complete and cpu is not None and wall is not None
    within = complete and cpu <= cpu_cap and wall <= wall_cap
    return output, dict(unit=unit, CPU_seconds=cpu, wall_seconds=wall,
        accounting_complete=bool(complete), within_attempt_ceiling=bool(within),
        CPU_ceiling=cpu_cap, wall_ceiling=wall_cap, exit_code=receipt.get('exit_code'),
        supervisor_exit_code=receipt.get('supervisor_exit_code'), termination=receipt.get('termination'),
        guard_error=receipt.get('guard_error'), cleanup_errors=receipt.get('cleanup_errors', []),
        frozen_seconds=(number(receipt['frozen_seconds'], 'guard frozen wall') if 'frozen_seconds' in receipt else None))


def constructor_record(record, *, model, units, arm, repetition, model_hash, source_hash, receipt_cpu, receipt_wall):
    require(isinstance(record, dict), 'Constructor record object required')
    partition = 'pauli' if arm == 'pauli' else 'conserving'
    physical = k.blocks(model, pauli=partition == 'pauli')
    require(record.get('schema') == 'prospective-window-constructor-v1' and record.get('artifact_complete') is True,
            'Complete constructor artifact schema required')
    require(record.get('arm') == arm and type(record.get('repetition')) is int and record['repetition'] == repetition and
            record.get('case') == model['case'] and record.get('partition') == partition and
            record.get('model_sha256') == model_hash and record.get('source_freeze_sha256') == source_hash,
            'Constructor arm/repetition/input/source identity differs')
    require(record.get('physical_error_read') is False and record.get('response_is_error_bound') is False,
            'Constructor error-oracle/certificate scope differs')
    require(number(record['training_time'], 'training time', positive=True) == k.TRAIN_TIME and
            integer(record['training_rotation_budget'], 'training budget') == k.TRAIN_CAP * len(model['coefficients']) and
            number(record['total_CPU_ceiling'], 'constructor CPU ceiling', positive=True) == k.TRAIN_CPU_SECONDS,
            'Frozen construction time/cap/allowance differs')
    require(integer(record['original_word_count'], 'original words', 1, k.MAX_TERMS) == len(model['coefficients']),
            'Constructor original word inventory differs')
    result = record['constructor'];require(isinstance(result, dict), 'Constructor result object required')
    require(result.get('partition') == partition, 'Returned partition differs')
    order = permutation(result['order'], len(physical))
    lengths = [len(physical[index]) for index in order]
    require(record.get('unit_lengths') == lengths and all(type(value) is int for value in record['unit_lengths']),
            'Original ordered word lengths differ')
    layout = [[[[site, letter] for site, letter in word], str(value)]
              for index in order for word, value in sorted(physical[index].items())]
    require(record.get('pulse_layout_sha256') == digest(json.dumps(layout, separators=(',', ':')).encode()),
            'Original pulse-layout hash differs')
    status = result.get('status');require(isinstance(status,str) and status in STATUSES, 'Unknown constructor status')
    kind = record.get('score_kind');require(kind is None or isinstance(kind,str) and kind in KINDS, 'Unknown score kind')
    expected_kind = ('retained_first_order_response' if arm in ('first', 'pauli') else
                     'diagonal_reference_cubic_response' if arm in ('cubic', 'window_sparse', 'window_magnitude') else None)
    score = result.get('score')
    if score is not None:
        score = number(score, 'constructor score')
        require(kind == expected_kind and kind in KINDS, 'Scored constructor/score-kind membership differs')
    elif kind is not None:
        require(kind == expected_kind, 'Unscored constructor/score-kind membership differs')
    require(expected_kind is None or status != 'complete' or score is not None,
            'Complete scored constructor requires a resolved score')
    components = result.get('score_components')
    require(kind != 'diagonal_reference_cubic_response' or score is None or components is not None,
            'Scored native constructor requires score components')
    if components is not None:
        require(isinstance(components, dict) and set(components) == {'value', 'cubic_response', 'steps'} and
                kind == 'diagonal_reference_cubic_response' and score is not None,
                'Native score components schema')
        require(number(components['value'], 'native score component') == score, 'Native score value differs')
        cubic=number(components['cubic_response'], 'native cubic response')
        steps=integer(components['steps'], 'native score steps', 1)
        expected_steps=(1 if len(lengths)==1 else
            (k.TRAIN_CAP*len(model['coefficients'])-lengths[0])//(2*sum(lengths)-lengths[0]-lengths[-1]))
        require(steps==expected_steps, 'Native score steps differ from original-word affordability')
        require((cubic == 0) == (score == 0), 'Native score/components unresolved zero differs')
        require(abs(score-cubic/steps**4)<=max(1e-30,abs(score)*1e-12), 'Native score normalization differs')
    parent_cpu = number(record['process_CPU_before_output'], 'parent CPU before output')
    parent_wall=number(record['wall_seconds_before_output'], 'constructor wall before output')
    require(receipt_wall is None or parent_wall<=receipt_wall, 'Constructor wall exceeds absolute guard wall')
    child_cpu = 0.
    stats = result.get('stats', {})
    require(isinstance(stats, dict), 'Constructor stats object')
    evaluations=integer(result['evaluations'], 'completed constructor evaluations',0,k.NATIVE_EVALUATIONS)
    moves=integer(result['accepted_moves'], 'accepted constructor moves')
    rescored=integer(stats.get('full_move_rescores',0), 'completed full move rescores')
    require(moves<=evaluations+rescored, 'Accepted moves exceed completed score work')
    if 'window_status' in stats:
        require(isinstance(stats['window_status'],str) and 0<len(stats['window_status'])<=256,
                'Window status must be a bounded phase string')
    if 'wall_seconds' in stats:
        require(receipt_wall is None or number(stats['wall_seconds'],'window whole wall')<=receipt_wall,
                'Window wall exceeds absolute guard wall')
    if 'CPU_seconds' in result:
        number(result['CPU_seconds'],'reported constructor CPU')
    if 'child_process_CPU_seconds' in stats:
        child_cpu = (None if stats['child_process_CPU_seconds'] is None else
                     number(stats['child_process_CPU_seconds'], 'owned child CPU'))
    require(receipt_cpu is None or parent_cpu + (child_cpu or 0.) <= receipt_cpu + ROUNDING_CPU_SECONDS,
            'Guard receipt is smaller than completed parent plus owned-child CPU')
    return dict(order=order, partition=partition, status=status, score=score, score_kind=kind,
        score_components=components, window_status=stats.get('window_status'),
        parent_CPU_before_output=parent_cpu, owned_child_CPU_seconds=child_cpu,
        reason=result.get('reason'), evaluations=evaluations, accepted_moves=moves)


def whole_orders_binding(root, index, registry, entries, models, source_path,
                         execution_root, execution_source_root, current):
    """Authenticate the whole attempted-call freeze, including other-case censors.

    Artifact presence does not establish successful construction, a complete
    physical panel, or whole-campaign comparison. Those claims remain separate.
    """
    descriptor=index['all_orders_freeze']
    _,frozen=artifact(root,descriptor)
    mode=index['mode'];repetitions=2 if mode=='main' else 1
    cases={spec['case'] for spec in registry['cases']} if mode=='main' else {index['case']}
    expected={(case,arm,repetition) for case in cases for arm in k.ARMS
              for repetition in range(1,repetitions+1)}
    require(frozen.get('schema')=='prospective-window-all-orders-freeze-v1' and
            frozen.get('mode')==mode and frozen.get('inventory_complete') is True and
            frozen.get('physical_error_read') is False and frozen.get('no_selected_repetition') is True,
            'Complete source-bound whole-order freeze schema/mode/flags required')
    require(type(frozen.get('expected_calls')) is int and frozen['expected_calls']==len(expected) and
            frozen.get('source_freeze')==index['source_freeze'],
            'Whole-order freeze source/call-count identity differs')
    require(isinstance(frozen.get('plan_sha256'),str) and HEX.fullmatch(frozen['plan_sha256']),
            'Whole-order freeze must bind the planned inventory hash')
    rows=frozen['artifacts']
    require(isinstance(rows,list) and len(rows)==len(expected), 'Whole-order freeze exact artifact roster')
    # Current final edges intentionally appear in both graphs. Only earlier
    # attempts seed collision sets; legitimate identical final references are
    # checked below and then registered once by the whole roster.
    finals={}
    for entry in entries:
        require(isinstance(entry,dict), 'Current constructor entry object required')
        attempts=entry['attempts']
        require(isinstance(attempts,list) and attempts, 'Unperformed call cannot satisfy Freeze C: actual receipt required')
        require(isinstance(attempts[-1],dict), 'Final attempt object required')
        finals[(index['case'],entry['arm'],entry['repetition'])]=attempts[-1]
    observed={};receipt_hashes=set();receipt_units=set();record_hashes=set();destinations=set();record_paths=set();censors=[]
    for retained in current:
        for attempt in retained['attempts'][:-1]:
            require('failure' not in attempt, 'Malformed earlier attempt cannot bind whole-order freeze')
            receipt_hashes.add(attempt['receipt']['sha256']);receipt_units.add(attempt['unit'])
            destinations.add(relative_path(root,attempt['constructor_output_path']).resolve())
            if attempt.get('constructor') is not None:
                record_hashes.add(attempt['constructor']['sha256'])
                record_paths.add(relative_path(root,attempt['constructor']['path']).resolve())
    for row in rows:
        require(isinstance(row,dict), 'Whole-order artifact object required')
        identity=(row['case'],row['arm'],integer(row['repetition'],'whole-freeze repetition',1,repetitions))
        require(identity in expected and identity not in observed, 'Duplicate or unknown whole-order identity')
        observed[identity]=row
        require(row.get('not_launched') is not True,
                'Unperformed whole-order call cannot satisfy Freeze C: actual receipt required')
        if identity in finals:
            require(finals[identity].get('receipt')==row.get('receipt') and
                    finals[identity].get('record')==row.get('record'),
                    'Current-case final descriptors differ from whole freeze')
        _,receipt=artifact(root,row['receipt'])
        require(row['receipt']['sha256'] not in receipt_hashes, 'Whole-order freeze reuses a guard receipt')
        receipt_hashes.add(row['receipt']['sha256'])
        target=models[identity[0]]
        destination,accounting=receipt_binding(root,receipt,arm=identity[1],repetition=identity[2],
            model_path=target['model_path'],model_hash=target['model_hash'],source_path=source_path,
            execution_root=execution_root,execution_source_root=execution_source_root,seen_receipts=receipt_units)
        require(destination.resolve() not in destinations, 'Whole-order freeze reuses a constructor destination')
        destinations.add(destination.resolve())
        require(accounting['accounting_complete'] and accounting['within_attempt_ceiling'],
                'Whole-order receipt accounting/attempt ceiling is incomplete')
        require(type(row.get('exit_code')) in (int,type(None)) and
                type(row.get('exit_code')) is type(receipt.get('exit_code')) and
                row.get('exit_code')==receipt.get('exit_code'), 'Whole-order receipt exit identity differs')
        if row.get('record') is None:
            require(not destination.exists(), 'Whole-order unindexed constructor output exists')
            censors.append(dict(case=identity[0],arm=identity[1],repetition=identity[2],
                receipt=row['receipt'],exit_code=row.get('exit_code'),
                reason=row.get('censor_reason','No final constructor artifact; attempted call retained')))
            continue
        record_path,record=artifact(root,row['record'])
        require(record_path.resolve()==destination.resolve(), 'Whole-order command/result path differs')
        require(record_path.resolve() not in record_paths, 'Whole-order freeze reuses a constructor record path')
        record_paths.add(record_path.resolve())
        require(row['record']['sha256'] not in record_hashes, 'Whole-order freeze reuses a constructor artifact')
        record_hashes.add(row['record']['sha256'])
        value=constructor_record(record,model=target['model'],units=target['units'],arm=identity[1],
            repetition=identity[2],model_hash=target['model_hash'],source_hash=index['source_freeze']['sha256'],
            receipt_cpu=accounting['CPU_seconds'],receipt_wall=accounting['wall_seconds'])
        if not clean_outer(accounting) or value['status']=='unexpected_error':
            censors.append(dict(case=identity[0],arm=identity[1],repetition=identity[2],
                receipt=row['receipt'],record=row['record'],exit_code=row.get('exit_code'),
                reason='Outer constructor did not complete cleanly; record retained but not admitted'))
    require(set(observed)==expected, 'Missing whole-order call identity')
    present=all(row.get('record') is not None for row in rows)
    require(type(frozen.get('all_calls_have_artifacts')) is bool and
            frozen['all_calls_have_artifacts']==present, 'Whole-order artifact presence flag differs')
    return dict(descriptor=descriptor,mode=mode,expected_calls=len(expected),
        inventory_complete=True,all_calls_attempted=True,all_calls_have_artifacts=present,
        all_calls_have_admitted_orders=not censors,censored_calls=censors,
        censored_call_count=len(censors),plan_sha256=frozen['plan_sha256'],
        current_case_final_descriptors_match=True,
        scope='Complete logical-call inventory with final artifacts; earlier attempts stay in per-case cost ledgers, and other-case censors do not establish a complete campaign comparison')


def clean_outer(accounting):
    return (accounting['exit_code']==0 and accounting['supervisor_exit_code']==0 and
            accounting['termination'] is None)


def allocations_preserve_budget(attempts):
    used_cpu=used_wall=0.
    for attempt in attempts:
        if (attempt.get('CPU_ceiling',float('inf'))>k.TRAIN_CPU_SECONDS-used_cpu+ALLOCATION_TOLERANCE or
                attempt.get('wall_ceiling',float('inf'))>k.TRAIN_WALL_SECONDS-used_wall+ALLOCATION_TOLERANCE):
            return False
        used_cpu+=attempt.get('CPU_seconds') or 0.
        used_wall+=attempt.get('wall_seconds') or 0.
    return True


def toml_scalar(value):
    if type(value) is bool:
        return 'true' if value else 'false'
    if type(value) is int:
        return str(value)
    if type(value) is float:
        require(math.isfinite(value), 'TOML output must be finite')
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    if isinstance(value, list):
        return '[' + ', '.join(toml_scalar(item) for item in value) + ']'
    raise ValueError('Unsupported bounded TOML scalar')


def candidate_bytes(header, rows):
    lines = [key + ' = ' + toml_scalar(value) for key, value in sorted(header.items())]
    for row in rows:
        lines.extend(['', '[[candidates]]'])
        lines.extend(key + ' = ' + toml_scalar(value) for key, value in sorted(row.items()))
    raw = ('\n'.join(lines) + '\n').encode('ascii')
    require(len(raw) <= k.MAX_MODEL_BYTES, 'Candidate manifest byte cap')
    return raw


def write_new_bytes(path, raw):
    temporary = path.with_name(path.name + '.partial')
    require(not path.exists() and not temporary.exists(), 'Refuse existing candidate or partial output')
    with temporary.open('xb') as stream:
        stream.write(raw)
    os.link(temporary, path)
    temporary.unlink()


def assemble(index_path, output):
    """Write PANEL.json last; CANDIDATES.toml exists only for a complete ledger.

    INVENTORY descriptors are relative to its parent directory. Absolute argv
    paths in receipts are mapped from explicit execution_root/source_root values
    (defaulting to the current artifact/source roots), never executed. Selected
    attempts are 1-based, must be the final attempted call, and all earlier costs
    remain charged against the same arm/repetition 64 CPU / 320 wall ceiling.
    """
    index_path = Path(index_path).resolve();output = Path(output)
    require(not output.exists(), 'Fresh panel output directory required')
    output.parent.mkdir(parents=True, exist_ok=True);output.mkdir()
    report = dict(schema='prospective-window-panel-audit-v1', status='incomplete',
        candidate_panel_complete=False, physical_error_read=False, constructors=[], failures=[],
        cost_scope='Each arm/repetition charged once across ALL retained attempts; repeated candidate provenance is not extra construction cost',
        guard_CPU_seconds=None, guard_wall_seconds=None, all_guard_accounting_complete=False)
    rows = []
    try:
        require(Path(k.__file__).resolve()==Path(__file__).resolve().with_name('study_common.py'),
                'Wrong study_common import identity')
        raw = raw_file(index_path);index = decode_json(raw)
        report['inventory_sha256'] = digest(raw)
        require(index.get('schema') == SCHEMA and index.get('mode') in ('main', 'preflight'), 'Panel inventory schema/mode')
        root = index_path.parent
        mode = index['mode'];repetitions = 2 if mode == 'main' else 1
        expected = {(arm, repetition) for arm in k.ARMS for repetition in range(1, repetitions + 1)}
        execution_root = PurePosixPath(index.get('execution_root', str(root)))
        execution_source_root = PurePosixPath(index.get('execution_source_root', str(k.ROOT)))
        require(execution_root.is_absolute() and execution_source_root.is_absolute(), 'Explicit absolute execution roots')
        source_path, frozen, registry, registry_hash = source_binding(root, index)
        model_path, toml_path, model, units, models = input_binding(root, index, frozen, registry)
        report.update(mode=mode, case=model['case'], n=model['n'], hf=str((1 << model['electrons']) - 1),
            source_commit=frozen['source_commit'], source_freeze=index['source_freeze'], input_freeze=index['input_freeze'],
            model_json=index['model_json'], model_toml=index['model_toml'], registry_sha256=registry_hash,
            execution_root=str(execution_root), execution_source_root=str(execution_source_root),
            expected_constructor_calls=len(expected), expected_logical_rows=44 * repetitions)
        entries = index['constructors']
        require(isinstance(entries, list) and len(entries) == len(expected), 'Exact frozen constructor roster required')
        require(all(isinstance(entry,dict) for entry in entries), 'Constructor entry objects required')
        identities = [(entry['arm'], integer(entry['repetition'], 'constructor repetition', 1, repetitions)) for entry in entries]
        require(len(set(identities)) == len(expected) and set(identities) == expected,
                'Duplicate/missing/unrecognized constructor identity')
        seen_receipts = set();seen_receipt_hashes = set();seen_destinations=set();seen_records=set();seen_record_paths=set();all_complete = True
        for entry in sorted(entries, key=lambda item: (item['repetition'], k.ARMS.index(item['arm']))):
            arm, repetition = entry['arm'], entry['repetition']
            retained = dict(arm=arm, repetition=repetition, status='censored', attempts=[],
                            guard_CPU_seconds=None, guard_wall_seconds=None, cost_complete=False)
            report['constructors'].append(retained)
            try:
                attempts = entry['attempts']
                require(isinstance(attempts, list) and len(attempts) <= 32, 'Bounded complete attempt ledger required')
                if not attempts:
                    require(entry.get('not_launched') is True and isinstance(entry.get('censor_reason'), str) and
                            entry['censor_reason'], 'Unlaunched censor requires explicit scheduler attestation')
                    retained.update(status='unperformed',reason=entry['censor_reason'],not_launched=True,
                                    cost_complete=False,cost_scope='unperformed; no completed-cost claim')
                    all_complete = False
                    continue
                selected = entry.get('selected_attempt')
                require(selected is None or integer(selected, 'selected attempt', 1, len(attempts)) == len(attempts),
                        'Selected attempt must be final; no outcome-based earlier selection')
                selected_value = None;selected_descriptor = None
                for position, attempt in enumerate(attempts, 1):
                    record = dict(attempt=position, status='incomplete')
                    retained['attempts'].append(record)
                    try:
                        require(isinstance(attempt,dict), 'Attempt object required')
                        record.update(receipt=attempt.get('receipt'),constructor=attempt.get('record'))
                        descriptor = attempt['receipt']
                        require(isinstance(descriptor,dict), 'Attempt receipt descriptor object required')
                        require(descriptor['sha256'] not in seen_receipt_hashes, 'Guard receipt charged more than once')
                        seen_receipt_hashes.add(descriptor['sha256'])
                        _, receipt = artifact(root, descriptor)
                        destination, accounting = receipt_binding(root, receipt, arm=arm, repetition=repetition,
                            model_path=model_path, model_hash=index['model_json']['sha256'], source_path=source_path,
                            execution_root=execution_root, execution_source_root=execution_source_root,
                            seen_receipts=seen_receipts)
                        record.update(accounting)
                        record['constructor_output_path'] = str(destination.relative_to(root))
                        require(destination.resolve() not in seen_destinations, 'Constructor output destination reused across attempts')
                        seen_destinations.add(destination.resolve())
                        if attempt.get('record') is not None:
                            record_path, original = artifact(root, attempt['record'])
                            require(record_path.resolve() == destination.resolve(), 'Guard output/constructor record path differs')
                            require(attempt['record']['sha256'] not in seen_records and record_path.resolve() not in seen_record_paths,
                                    'Constructor record reused across attempts')
                            seen_records.add(attempt['record']['sha256']);seen_record_paths.add(record_path.resolve())
                            value = constructor_record(original, model=model, units=units, arm=arm, repetition=repetition,
                                model_hash=index['model_json']['sha256'], source_hash=index['source_freeze']['sha256'],
                                receipt_cpu=accounting['CPU_seconds'],receipt_wall=accounting['wall_seconds'])
                            record.update(status='retained_record', constructor=attempt['record'],
                                          constructor_status=value['status'], scored=value['score'] is not None,
                                          outer_constructor_success=clean_outer(accounting))
                            if not clean_outer(accounting):
                                record['status']='failed_outer_record_retained'
                            if position == selected and clean_outer(accounting):
                                selected_value=value;selected_descriptor=attempt['record']
                        else:
                            require(not destination.exists(), 'Unindexed constructor output exists; retain/hash it explicitly')
                            record['status']='censored_no_constructor_record'
                            record['censor_reason']=attempt.get('censor_reason', 'No final constructor artifact; guard receipt retained')
                            partial=destination.with_name(destination.name+'.partial')
                            if partial.is_file():
                                record['retained_partial']=dict(path=str(partial.relative_to(root)),sha256=digest(raw_file(partial)))
                    except (ValueError, KeyError, TypeError, OSError, RecursionError) as error:
                        record['failure']=dict(type=type(error).__name__,message=str(error)[:500])
                costs = retained['attempts']
                complete_cost = all(attempt.get('accounting_complete') is True for attempt in costs)
                known_cpu=sum(attempt['CPU_seconds'] for attempt in costs if attempt.get('CPU_seconds') is not None)
                known_wall=sum(attempt['wall_seconds'] for attempt in costs if attempt.get('wall_seconds') is not None)
                retained.update(known_guard_CPU_seconds=known_cpu, known_guard_wall_seconds=known_wall,
                    guard_CPU_seconds=known_cpu if complete_cost else None,
                    guard_wall_seconds=known_wall if complete_cost else None,cost_complete=complete_cost)
                allocated_without_reset=allocations_preserve_budget(costs)
                within = complete_cost and known_cpu <= k.TRAIN_CPU_SECONDS and known_wall <= k.TRAIN_WALL_SECONDS and all(
                    attempt.get('within_attempt_ceiling') is True for attempt in costs) and allocated_without_reset
                retained['attempt_allocations_preserve_whole_call_ceiling']=allocated_without_reset
                retained['within_total_construction_ceiling']=bool(within)
                require(all('failure' not in attempt for attempt in costs), 'Malformed/unauthenticated retained attempt')
                require(complete_cost and within, 'Incomplete accounting or cumulative construction ceiling exceeded')
                if selected_value is None:
                    require(selected is None and isinstance(entry.get('censor_reason'), str) and entry['censor_reason'],
                            'Missing selected record requires explicit censor reason')
                    retained['reason']=entry['censor_reason'];all_complete=False
                    continue
                value=selected_value
                require(value['status']!='unexpected_error', 'Unexpected constructor error cannot supply an admitted candidate')
                retained.update(status='retained_order', selected_attempt=selected, selected_record=selected_descriptor,
                    constructor=value, scored=value['score'] is not None)
                base=dict(arm=arm,repetition=repetition,partition=value['partition'],order=[i+1 for i in value['order']],
                    constructor_arm=arm,constructor_repetition=repetition,constructor_status=value['status'],
                    constructor_scored=value['score'] is not None,constructor_record_sha256=selected_descriptor['sha256'],
                    constructor_source_freeze_sha256=index['source_freeze']['sha256'],
                    constructor_guard_CPU_seconds=known_cpu,constructor_guard_wall_seconds=known_wall,
                    constructor_training_time=k.TRAIN_TIME,constructor_training_rotation_budget=k.TRAIN_CAP*len(model['coefficients']))
                if value['score'] is not None:
                    base.update(score_kind=value['score_kind'],constructor_score=value['score'])
                else:
                    base['constructor_unscored_reason']=value['status']
                if isinstance(value['window_status'],str):
                    base['constructor_window_status']=value['window_status']
                for formula in k.FORMULAS[arm]:
                    for duration in k.TIMES:
                        for multiplier in k.CAPS:
                            rows.append(dict(base,formula=formula,time=duration,budget_multiplier=multiplier,
                                name=f'{arm}_r{repetition}_f{formula}_t{duration}_b{multiplier}'))
            except (ValueError, KeyError, TypeError, OSError, RecursionError) as error:
                retained['reason']=dict(type=type(error).__name__,message=str(error)[:500]);all_complete=False
        report['arm_totals']={}
        for arm in k.ARMS:
            selected=[entry for entry in report['constructors'] if entry['arm']==arm]
            complete=all(entry['cost_complete'] for entry in selected)
            report['arm_totals'][arm]=dict(expected_calls=repetitions,retained_calls=len(selected),accounting_complete=complete,
                CPU_seconds=sum(entry['guard_CPU_seconds'] for entry in selected) if complete else None,
                wall_seconds=sum(entry['guard_wall_seconds'] for entry in selected) if complete else None,
                known_CPU_seconds=sum(entry.get('known_guard_CPU_seconds',entry.get('guard_CPU_seconds') or 0.) for entry in selected),
                known_wall_seconds=sum(entry.get('known_guard_wall_seconds',entry.get('guard_wall_seconds') or 0.) for entry in selected))
        complete_cost=all(entry['cost_complete'] for entry in report['constructors'])
        report['all_guard_accounting_complete']=complete_cost
        if complete_cost:
            report['guard_CPU_seconds']=sum(entry['guard_CPU_seconds'] for entry in report['constructors'])
            report['guard_wall_seconds']=sum(entry['guard_wall_seconds'] for entry in report['constructors'])
        report['available_logical_rows']=len(rows)
        try:
            whole=whole_orders_binding(root,index,registry,entries,models,source_path,
                execution_root,execution_source_root,report['constructors'])
            report['all_orders_freeze']=whole
        except (ValueError,KeyError,TypeError,OSError,RecursionError) as error:
            report['failures'].append(dict(phase='all_orders_freeze',type=type(error).__name__,message=str(error)[:500]))
            all_complete=False
        if all_complete:
            require(len(rows)==44*repetitions and len({row['name'] for row in rows})==len(rows), 'Exact unique physical candidate roster')
            header=dict(schema='prospective-window-candidates-v1',status='complete',mode=mode,case=model['case'],
                input_sha256=index['model_toml']['sha256'],source_json_sha256=index['model_json']['sha256'],
                source_freeze_sha256=index['source_freeze']['sha256'],input_freeze_sha256=index['input_freeze']['sha256'],
                all_orders_freeze_sha256=index['all_orders_freeze']['sha256'],
                whole_constructor_artifacts_present=whole['all_calls_have_artifacts'],
                inventory_sha256=report['inventory_sha256'],source_commit=frozen['source_commit'],n=model['n'],
                hf=str((1<<model['electrons'])-1),expected_logical_rows=44*repetitions,
                expected_constructor_calls=len(expected),candidate_panel_complete=True,physical_error_read=False,
                construction_accounting_complete=True,construction_within_limits=True)
            encoded=candidate_bytes(header,rows)
            write_new_bytes(output/'CANDIDATES.toml',encoded)
            report.update(status='complete',candidate_panel_complete=True,
                candidates=dict(path='CANDIDATES.toml',sha256=digest(encoded)),
                scored_constructor_calls=sum(entry.get('scored',False) for entry in report['constructors']))
    except (ValueError, KeyError, TypeError, OSError, RecursionError) as error:
        report['failures'].append(dict(type=type(error).__name__,message=str(error)[:500]))
    # This append-only audit is the publication marker, including on refusal.
    # A physical dispatcher must require its complete status and candidate hash.
    k.write_new(output/'PANEL.json',report)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inventory',type=Path)
    parser.add_argument('output',type=Path,help='new append-only panel directory')
    args=parser.parse_args(argv)
    report=assemble(args.inventory,args.output)
    print(json.dumps(dict(status=report['status'],candidate_panel_complete=report['candidate_panel_complete'],
                          audit=str(args.output/'PANEL.json'))))
    return 0 if report['candidate_panel_complete'] else 2


if __name__=='__main__':
    raise SystemExit(main())
