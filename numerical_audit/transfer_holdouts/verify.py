"""Complete transfer audit using literal objectives and circuit reference actions.

Producer construction, plan, score, propagation, fusion and comparison functions
are not called. The separate CAR/tensor/reference auditor shares source-bound
I/O and constants but supplies independent literal input checks. Archived Python
sources are never run.
"""
import argparse
from collections import defaultdict
import hashlib
import importlib.metadata
import io
import itertools
import json
import math
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time
import zipfile
import numpy as np
import circuit_reference
import objective_reference
import verify_references

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATES = ('HF', 'ground', 'excited', 'coherent_real', 'coherent_imag')
ARMS = ('hf_current', 'hf_first_JD', 'hf_pauli_JD', 'hf_cubic_JD',
        'state_first_JD', 'state_first_symmetric', 'state_pauli_JD', 'state_cubic_JD')
METHODS = ('first_JD', 'pauli_JD', 'cubic_JD')
TIMES, CAPS, TRAIN_TIME = (.59, 1.41, 2.73), (10, 30), 1.11
VERSION = 'line-pauli-exact-angle-0.1.0'
FLOOR, PAIR_TOLERANCE, STATE_TOLERANCE, METRIC_TOLERANCE = 1e-12, 2e-12, 2e-10, 2e-12
require = circuit_reference.require
for _module in (circuit_reference, objective_reference, verify_references):
    require(Path(_module.__file__).resolve() == HERE / (_module.__name__ + '.py'), 'Wrong independent source binding')


def safe(name):
    path = PurePosixPath(name)
    require(not path.is_absolute() and '..' not in path.parts and str(path) == name, 'Unsafe artifact name')
    return name


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    path = Path(path)
    require(path.stat().st_size <= 8 * 1024**2, 'Audit JSON byte ceiling')
    return json.loads(path.read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def inventory(root):
    names, total, entries = set(), 0, 0
    for path in Path(root).rglob('*'):
        entries += 1
        require(entries <= 17024 and not path.is_symlink(), 'Audit directory/symlink ceiling')
        if path.is_file():
            total += path.stat().st_size
            require(len(names) < 16000 and total <= 512 * 1024**2, 'Audit payload storage ceiling')
            names.add(str(path.relative_to(root)))
    return names


def array(path, key, shape, digest):
    require(Path(path).stat().st_size <= 32768 and sha(path) == digest, 'State archive identity/size differs')
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        require(len(entries) == 1 and entries[0].filename == key + '.npy' and
                entries[0].file_size <= 16384, 'Bounded single-array state archive required')
        raw = archive.read(entries[0])
    stream = io.BytesIO(raw)
    version = np.lib.format.read_magic(stream)
    require(version in ((1, 0), (2, 0)), 'Supported NPY header required')
    parser = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
    actual_shape, _, dtype = parser(stream)
    require(actual_shape == shape and dtype.kind == 'c' and dtype.itemsize in (16, 32) and
            len(raw) - stream.tell() == math.prod(shape) * dtype.itemsize, 'State array shape/dtype before allocation')
    value = np.load(io.BytesIO(raw), allow_pickle=False)
    require(np.all(np.isfinite(value)), 'Nonfinite state archive')
    return value


def near(actual, expected, *, relative=2e-9, absolute=2e-18):
    require(type(actual) in (float, int) and math.isfinite(actual) and math.isfinite(float(expected)),
            'Nonfinite numeric witness')
    require(abs(actual - expected) <= absolute + relative * max(abs(actual), abs(expected)),
            'Independent numeric witness differs')


def verifier_sources():
    paths = {ROOT / 'numerical_audit/state_holdouts/response.py',
             ROOT / 'numerical_audit/prospective_reference/engine.py'}
    for module in list(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT.resolve()):
                paths.add(path)
    return {str(path.relative_to(ROOT.resolve())): sha(path) for path in sorted(paths)}


def positive_time(value):
    require(type(value) in (float, int) and math.isfinite(value) and value >= 0, 'Invalid measured CPU charge')
    return value


def expected_panel(cases, preflight):
    times = (2.73,) if preflight else TIMES
    roles = {
        'hf_current': ('HF', 'current', 1), 'hf_first_JD': ('HF', 'first_JD', 1),
        'hf_pauli_JD': ('HF', 'pauli_JD', 2), 'hf_cubic_JD': ('HF', 'cubic_JD', 2),
        'state_first_JD': (None, 'first_JD', 1), 'state_first_symmetric': (None, 'first_JD', 2),
        'state_pauli_JD': (None, 'pauli_JD', 2), 'state_cubic_JD': (None, 'cubic_JD', 2),
    }
    constructors, candidates, logical = [], {}, []
    for ci, case in enumerate(cases):
        constructors.append(dict(constructor=case + '.HF.current', case=case, state='HF', method='current',
                                 training_time=1., training_rotation_cap_units=None))
        for si, state in enumerate(STATES):
            offset = (ci + si) % 3
            for method in METHODS[offset:] + METHODS[:offset]:
                constructors.append(dict(constructor=case + '.' + state + '.' + method, case=case,
                                         state=state, method=method, training_time=TRAIN_TIME,
                                         training_rotation_cap_units=30))
            for ti, T in enumerate(times):
                for cap_index, cap in enumerate(CAPS):
                    offset = (ci + si + ti + cap_index) % 8
                    for arm in ARMS[offset:] + ARMS[:offset]:
                        training_state, method, formula = roles[arm]
                        training_state = training_state or state
                        name = case + '.' + training_state + '.' + method
                        candidate = name + '.eval.' + state + '.f' + str(formula) + '.t' + str(ti) + '.cap' + str(cap)
                        candidates.setdefault(candidate, dict(candidate=candidate, case=case, state=state,
                            time_index=ti, T=T, cap_multiple=cap, constructor=name, training_state=training_state,
                            method=method, partition='pauli' if method == 'pauli_JD' else 'conserving', formula=formula))
                        logical.append(dict(case=case, state=state, time_index=ti, T=T, cap_multiple=cap, arm=arm,
                            candidate=candidate, constructor=name,
                            reference_charge_states=[label for label in STATES if label in (training_state, state)],
                            compiled_resource_scope='HF preparation plus evolution' if state == 'HF' else
                                'evolution only; quantum preparation cost unknown'))
    require(len(constructors) == 16 * len(cases) and len(candidates) == 74 * len(cases) * len(times),
            'Independent fixed-panel cardinality differs')
    return constructors, list(candidates.values()), logical, times


def layout(model, record, partition):
    blocks = [{tuple((q, p) for q, p in word): int(c) for word, c in block} for block in model['blocks']]
    if partition == 'pauli':
        blocks = [{tuple((q, p) for q, p in word): int(c)} for word, c in sorted(
            (tuple((q, p) for q, p in word), int(c)) for word, c in model['coefficients'])]
    order = record['order']
    require(all(type(i) is int for i in order) and sorted(order) == list(range(len(blocks))), 'Incomplete order witness')
    expected = dict(order=order, partition=partition, unit_lengths=[len(blocks[i]) for i in order],
        schedule=[[[[q, p] for q, p in word], str(c)] for i in order for word, c in sorted(blocks[i].items())])
    require(record == expected, 'Frozen layout differs from actual complete target units')
    return expected


def steps_and_rotations(lengths, formula, cap):
    U = sum(lengths)
    steps = cap if formula == 1 else (1 if len(lengths) == 1 else
                                      max(0, (cap * U - lengths[0]) // (2 * U - lengths[0] - lengths[-1])))
    require(type(steps) is int and 1 <= steps <= 128, 'Independent affordable-step domain')
    sequence = []
    for _ in range(steps):
        if formula == 1:
            sequence.extend(range(len(lengths)))
        else:
            sequence.extend(list(range(len(lengths) - 1)) + [len(lengths) - 1] +
                            list(reversed(range(len(lengths) - 1))))
    previous, rotations = None, 0
    for index in sequence:
        if formula == 1 or index != previous:
            rotations += lengths[index]
        previous = index
    require(rotations <= cap * U, 'Independent rotation cap exceeded')
    return steps, rotations


def check_score(context, record, method, counters):
    order = record['order']
    checked = (context.cubic_score(order, TRAIN_TIME, 30 * len(context.model['coefficients'])) if method == 'cubic_JD'
               else context.first_score(order, TRAIN_TIME, 30 * len(context.model['coefficients']),
                                        'pauli' if method == 'pauli_JD' else 'conserving'))
    for field in ('JD2' if method == 'cubic_JD' else 'JD1', 'score'):
        near(record[field], checked[field])
    require(record['steps'] == checked['steps'], 'Training affordability witness differs')
    if method == 'cubic_JD':
        require(record['column_sources'] == checked['column_sources'], 'Coherent column support differs')
        columns = record['cubic_columns']
        require(isinstance(columns, list) and len(columns) == context.dimension and
                all(isinstance(row, list) and len(row) == len(context.support) for row in columns),
                'Cubic column shape before allocation')
        require(np.max(abs(np.asarray(columns) - np.asarray(checked['cubic_columns']))) <= 2e-10,
                'Independent signed cubic columns differ')
    counters['objective_records'] += 1


def check_search(context, record, method, starts, counters):
    require(len(record['seeds']) == 3 and [row['order'] for row in record['seeds']] == starts, 'Initializer differs')
    for seed in record['seeds']:
        check_score(context, seed, method, counters)
    result = record['result']
    pauli = method == 'pauli_JD'
    allowance, moves = (8., 10) if pauli else (2., 5)
    m = len(context.paulis) if pauli else len(context.blocks)
    cap = min(325120, 20 * m * (m - 1)) if pauli else 3000
    require(result['nominal_search_CPU_seconds'] == allowance and result['evaluation_cap'] == cap and
            result['max_moves_per_start'] == moves and len(result['starts']) == 3, 'Search controls changed')
    attempted = completed = accepted = 0
    for index, (seed, search) in enumerate(zip(record['seeds'], result['starts'])):
        require(search['initial_order'] == seed['order'] and search['initial_score'] == seed['score'], 'Start differs')
        current = seed
        require(len(search['accepted']) <= moves, 'Accepted move cap exceeded')
        for update in search['accepted']:
            source, destination = update['move']
            require(type(source) is int and type(destination) is int and 0 <= source < m and
                    0 <= destination < m and source != destination, 'Invalid insertion witness')
            trial = list(current['order'])
            unit = trial.pop(source)
            trial.insert(destination, unit)
            require(trial == update['record']['order'], 'Retained move does not reconstruct its order')
            check_score(context, update['record'], method, counters)
            require(update['record']['score'] < current['score'] - max(1e-30, 1e-12 * abs(current['score'])),
                    'Accepted move violates the frozen strict improvement rule')
            current = update['record']
            accepted += 1
        require(search['record'] == current, 'Retained start incumbent differs')
        check_score(context, search['record'], method, counters)
        limit = cap // 3 + (index < cap % 3)
        require(0 <= search['completed'] <= search['attempted'] <= limit, 'Search effort exceeds its slice')
        require(sum(row['attempted'] for row in search['sweeps']) == search['attempted'] and
                sum(row['completed'] for row in search['sweeps']) == search['completed'], 'Sweep work differs')
        attempted += search['attempted']
        completed += search['completed']
    best = min((row['record'] for row in result['starts']), key=lambda row: (row['score'], tuple(row['order'])))
    require(result['record'] == best and result['attempted'] == attempted and result['completed'] == completed and
            result['accepted_rescores'] >= accepted and result['attempted_rescores'] >= result['accepted_rescores'],
            'Winner or cumulative search accounting differs')
    check_score(context, result['record'], method, counters)
    require(record['record']['order'] == result['record']['order'], 'Frozen layout does not use the selected order')


def check_constructor(context, record, method, current, reference_cost, counters, *, nested=False):
    outer = positive_time(record.get('caller_CPU_seconds', record.get('actual_local_CPU_seconds')))
    counters['constructor_checks'] += 1
    if nested:
        counters['nested_warm_constructors'] += 1
    near(record['reference_preparation_CPU_charge'], reference_cost, absolute=1e-10, relative=0)
    near(record['standalone_reference_and_ordering_CPU_charge'], reference_cost + outer, absolute=1e-10, relative=0)
    phases = sum(positive_time(value) for value in record['phases'].values())
    require(outer + 1e-9 >= phases, 'Constructor phase costs exceed outer return charge')
    actual = (record['warm_start']['current_placement'] if method == 'cubic_JD' else record['current_placement'])
    for field in ('order', 'start_slot', 'slots', 'history', 'initial_slot_scores', 'exact_weight_denominator',
                  'exact_score_numerator', 'evaluations', 'diagonal_priority', 'off_priority'):
        require(actual[field] == current[field], 'Literal integer legacy placement differs: ' + field)
    near(actual['score'], current['score'])
    partition = 'pauli' if method == 'pauli_JD' else 'conserving'
    layout(context.model, record['record'], partition)
    if method == 'current':
        require(record['training_time'] == 1. and record['record']['order'] == current['order'], 'Current order differs')
        return
    require(record['training_time'] == TRAIN_TIME, 'Training time changed')
    if method == 'cubic_JD':
        warm = record['warm_start']
        check_constructor(context, warm, 'first_JD', current, 0., counters, nested=True)
        require(record['phases']['charged_first_order_warm_constructor_CPU_seconds'] + 1e-9 >=
                warm['actual_local_CPU_seconds'], 'Fresh first-order warm constructor not fully charged')
        start = warm['record']['order']
        starts = [start, list(reversed(start)), current['order']]
    else:
        require(record['current_order'] == current['order'], 'Saved current initializer differs')
        start = current['order']
        if partition == 'pauli':
            words = sorted(context.words)
            blocks = [sorted(tuple((q, p) for q, p in word) for word, value in block)
                      for block in context.model['blocks']]
            start = [words.index(word) for index in start for word in blocks[index]]
        starts = [start, list(reversed(start)), list(range(len(context.paulis) if partition == 'pauli' else len(context.blocks)))]
    check_search(context, record, method, starts, counters)


def decision(a, b):
    require(a['status'] in ('complete', 'censored') and b['status'] in ('complete', 'censored'),
            'Unknown independent circuit status')
    if a['status'] != 'complete' or b['status'] != 'complete':
        return dict(decision='unavailable', informative=False, a_status=a['status'], b_status=b['status'],
                    a_error=a.get('error'), b_error=b.get('error'))
    x, y = a['error'], b['error']
    require(all(type(value) in (float, int) and math.isfinite(value) and 0 <= value <= 1 + 1e-12 for value in (x, y)),
            'Independent finite comparison errors required')
    floor = max(x, y) <= FLOOR
    which = 'floor_tie' if floor else 'numerical_tie' if abs(x - y) <= PAIR_TOLERANCE else 'a' if x < y else 'b'
    return dict(decision=which, informative=not floor, a_status='complete', b_status='complete',
                a_error=x, b_error=y, a_minus_b=x - y, both_at_floor=floor)



def write_candidate_progress(receipt, name, physical):
    # Verification outcome and physical completion/censoring are separate.
    write(Path(receipt) / 'progress/candidates' / (name + '.json'),
          dict(status='passed', physical=dict(physical)))


def verify(output, receipt_dir):
    output, receipt = Path(output), Path(receipt_dir)
    require(not receipt.exists() and not receipt.resolve().is_relative_to(output.resolve()),
            'Fresh independent receipts outside the producer payload required;retain earlier attempts')
    receipt.mkdir(parents=True)
    cpu, wall = time.process_time(), time.monotonic()
    protocol, frozen, plan = (read(output / name) for name in
                              ('PROTOCOL_FROZEN.json', 'SCHEDULES_FROZEN.json', 'PLAN_FROZEN.json'))
    for name, digest in frozen['artifacts'].items():
        require(sha(output / safe(name)) == digest, 'Frozen transfer artifact changed')
    require(protocol['plan_sha256'] == sha(output / 'PLAN_FROZEN.json') and
            protocol['reviewed_practical_interface_sha256'] == sha(output / 'PRACTICAL_INTERFACE.json'), 'Protocol binding differs')
    for name, digest in protocol['sources'].items():
        raw = (output / 'source_snapshot' / safe(name)).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == digest and raw == subprocess.check_output(
            ['git', '-C', str(ROOT), 'show', protocol['source_commit'] + ':' + name]), 'Producer Git source differs')
    dependency = read(output / 'PRACTICAL_INTERFACE.json')
    require(dependency['state'] == 'reviewed_frozen_interface' and dependency['repository'] == 'bernalde/qhat' and
            dependency['issue'] == 16 and dependency['automated_handoff'] == 'READY FOR HUMAN REVIEW' and
            dependency['requested_reviewers'] == ['gpt-6-astra/high', 'gemini-4-argon-high/high'] and
            dependency['actionable_findings_remaining'] == 0 and dependency['controller_policy_version'] == 2 and
            dependency['cost_convention'] == 'caller_through_return_including_serialization_hash_cleanup',
            'Unreviewed or incompatible practical dependency')
    for name, digest in dependency['source_sha256'].items():
        raw = subprocess.check_output(['git', '-C', str(ROOT), 'show', dependency['head'] + ':' + safe(name)])
        require(hashlib.sha256(raw).hexdigest() == digest, 'Reviewed dependency Git source differs')
    require(sha(ROOT / 'numerical_audit/compiled_accuracy/compiler.py') ==
            dependency['source_sha256']['numerical_audit/compiled_accuracy/compiler.py'], 'Reviewed compiler differs')
    require({str(path.relative_to(output / 'source_snapshot')) for path in (output / 'source_snapshot').rglob('*') if path.is_file()} ==
            set(protocol['sources']), 'Producing source snapshot inventory differs')
    reference_root = output / 'reference_evidence'
    input_audit = verify_references.verify(reference_root, receipt_dir=receipt / 'input_reference_audit')
    registry = read(HERE / 'REGISTRY.json')
    cases = [registry['chemistry_preflight']['case']] if protocol['preflight'] else [row['case'] for row in registry['cases']]
    training, candidates, logical, times = expected_panel(cases, protocol['preflight'])
    require(plan['training'] == training and plan['circuits'] == candidates and plan['logical'] == logical and
            plan['states'] == list(STATES) and plan['arms'] == list(ARMS) and plan['validation_times'] == list(times) and
            plan['error_floor'] == FLOOR and plan['comparison_absolute_tolerance'] == PAIR_TOLERANCE,
            'Independent complete-panel design differs')
    require(plan['preflight'] == protocol['preflight'] and plan['training_time'] == TRAIN_TIME and
            plan['legacy_current_time'] == 1. and plan['cap_multiples'] == list(CAPS) and
            plan['pairs'] == [dict(a=a, b=b) for a, b in itertools.combinations(ARMS, 2)] and
            plan['expected_constructors'] == len(training) and plan['expected_unique_circuits'] == len(candidates) and
            plan['expected_logical_circuits'] == len(logical) and
            plan['expected_reference_evolutions'] == len(cases) * 5 * len(times) and
            plan['expected_comparisons'] == len(cases) * 5 * len(times) * 2 * 28 and
            plan['no_validation_errors_used_for_construction'] is True and
            frozen['unique_candidates'] == len(candidates) and frozen['logical_candidates'] == len(logical),
            'Declared panel counts,controls or sharing differ')
    counts = {'constructor_checks': 0, 'main_constructors': 0, 'warmup_constructors': 0,
              'nested_warm_constructors': 0, 'objective_records': 0, 'complete': 0, 'censored': 0}
    contexts, references, currents, records = {}, {}, {}, {}
    warm_ordering = warm_compilation = 0.
    warm_trace_files = set()
    for case in cases:
        model = read(reference_root / 'inputs_evidence/inputs' / case / 'model.json')
        require(type(model['n']) is int and 1 <= model['n'] <= 6, 'Independent register bound before allocation')
        dimension = 1 << model['n']
        archive = read(reference_root / 'references' / (case + '.json'))
        values = array(reference_root / safe(archive['vectors_file']), 'vectors', (dimension, 5), archive['vectors_sha256'])
        references[case] = archive
        for index, state in enumerate(STATES):
            contexts[case, state] = objective_reference.Reference(model, values[:, index])
        currents[case] = contexts[case, 'HF'].current_placement()
        warm = read(output / 'warmups' / (case + '.json'))
        require(set(warm) == {'current', 'first_JD', 'pauli_JD', 'cubic_JD'}, 'Warmup cohort differs')
        for method, record in warm.items():
            check_constructor(contexts[case, 'HF'], record, method, currents[case], 0., counts)
            counts['warmup_constructors'] += 1
            warm_ordering += positive_time(record['caller_CPU_seconds'])
            formulas = [1, 2] if method == 'first_JD' else [1] if method == 'current' else [2]
            require([row['identity']['formula'] for row in record['compiled']] == formulas, 'Compiler warmup paths differ')
            for compiled in record['compiled']:
                formula = compiled['identity']['formula']
                step_count, rotations = steps_and_rotations(record['record']['unit_lengths'], formula, 10)
                expected_identity = dict(purpose='compilation_warmup', case=case, method=method, state='HF', T=TRAIN_TIME,
                    formula=formula, steps=step_count, cap_multiple=10, preparation='HF_prefix', compilation_status='complete')
                expected_path = 'warmup_traces/' + case + '.' + method + '.f' + str(formula) + '.ndjson.gz'
                require(compiled['identity'] == expected_identity and compiled['trace_file'] == expected_path and
                        compiled['status'] == 'complete' and compiled['ideal_rotations'] == rotations and
                        compiled['trace_sha256'] == sha(output / expected_path), 'Warmup circuit identity differs')
                warm_trace_files.add(expected_path)
                warm_compilation += positive_time(compiled['outer_compilation_CPU_seconds'])
                require(compiled['outer_compilation_CPU_seconds'] + 1e-9 >= compiled['compilation_CPU_seconds'],
                        'Warmup final constructor return not charged')
                state, costs = contexts[case, 'HF'].audit_compiled_trace(output / expected_path, compiled['identity'], VERSION)
                expected = contexts[case, 'HF'].ideal_formula(dict(record['record'], formula=formula), step_count, TRAIN_TIME)
                require(np.max(abs(state - expected)) < STATE_TOLERANCE and costs == compiled['compiled_resources'], 'Warmup trace differs')
    for specification in training:
        name, case, state, method = (specification[key] for key in ('constructor', 'case', 'state', 'method'))
        record = read(output / 'constructors' / (name + '.json'))
        require(record['constructor_identity'] == specification and not record['warmup'], 'Constructor identity differs')
        check_constructor(contexts[case, state], record, method, currents[case], references[case]['states'][state]['CPU_seconds'], counts)
        records[name] = record
        counts['main_constructors'] += 1
        write(receipt / 'progress/constructors' / (name + '.json'), dict(status='passed', source_sha256=sha(output / 'constructors' / (name + '.json'))))
    errors, physical, exact_cache, state_files, reference_files = {}, {}, {}, set(), set()
    compilation_cost = 0.
    max_state = max_metric = 0.
    for specification in candidates:
        name, case, label = (specification[key] for key in ('candidate', 'case', 'state'))
        candidate_path = output / 'candidates' / (name + '.json')
        candidate = read(candidate_path)
        require(all(candidate[key] == value for key, value in specification.items()), 'Candidate identity differs')
        context = contexts[case, label]
        record = dict(records[specification['constructor']]['record'], formula=specification['formula'])
        steps, rotations = steps_and_rotations(record['unit_lengths'], record['formula'], specification['cap_multiple'])
        require(candidate['steps'] == steps and candidate['expected_ideal_rotations'] == rotations and
                candidate['trace_sha256'] == sha(output / safe(candidate['trace_file'])), 'Candidate circuit binding differs')
        constructor_path = 'constructors/' + specification['constructor'] + '.json'
        expected_identity = dict(candidate=name, case=case, constructor=specification['constructor'],
            constructor_sha256=sha(output / constructor_path), model_sha256=references[case]['model_sha256'],
            reference_vectors_sha256=references[case]['vectors_sha256'], state=label, T=specification['T'],
            formula=specification['formula'], steps=steps, cap_multiple=specification['cap_multiple'],
            preparation='HF_prefix' if label == 'HF' else 'external_state_not_compiled',
            compilation_status=candidate['status'])
        require(candidate['identity'] == expected_identity and candidate['trace_file'] == 'traces/' + name + '.ndjson.gz',
                'Complete circuit/reference/preparation identity differs')
        ordering_cost = records[specification['constructor']]['caller_CPU_seconds']
        refs = [state for state in STATES if state in (specification['training_state'], label)]
        reference_cost = sum(references[case]['states'][state]['CPU_seconds'] for state in refs)
        input_cost = positive_time(context.model['generation_CPU_seconds'])
        outer_compile = positive_time(candidate['outer_compilation_CPU_seconds'])
        require(outer_compile + 1e-9 >= positive_time(candidate['compilation_CPU_seconds']) and
                candidate['reference_charge_states'] == refs, 'Compilation return/reference sharing charge differs')
        for field, value in dict(ordering_CPU_charge=ordering_cost, standalone_reference_CPU_charge=reference_cost,
                input_generation_CPU_charge=input_cost, standalone_without_input_CPU_charge=reference_cost + ordering_cost + outer_compile,
                standalone_including_input_CPU_charge=input_cost + reference_cost + ordering_cost + outer_compile).items():
            near(candidate[field], value, absolute=1e-9, relative=0)
        require(candidate['quantum_preparation_circuit_cost'] == ('HF prefix included' if label == 'HF' else 'unknown'),
                'Quantum preparation cost is missing or overstated')
        compilation_cost += outer_compile
        stored = read(output / 'physical' / (name + '.json'))
        require(stored['candidate_sha256'] == sha(candidate_path), 'Physical row source differs')
        for key, value in candidate.items():
            require(stored[key] == value, 'Physical row changed a constructor/circuit field')
        if candidate['status'] == 'complete':
            require(candidate['ideal_rotations'] == rotations, 'Complete ideal rotation count differs')
            state, costs = context.audit_compiled_trace(output / candidate['trace_file'], candidate['identity'], VERSION)
            ideal = context.ideal_formula(record, steps, specification['T'])
            saved = array(output / safe(stored['state_file']), 'vector', (context.dimension,), stored['state_sha256'])
            difference = max(float(np.max(abs(state - saved))), float(np.max(abs(state - ideal))))
            max_state = max(max_state, difference)
            require(difference <= STATE_TOLERANCE and costs == candidate['compiled_resources'], 'Primitive/formula/DAG audit differs')
            key = (case, label, specification['time_index'])
            if key not in exact_cache:
                exact_cache[key] = context.exact_action(specification['T'])
                reference_path = 'exact_references/' + case + '.' + label + '.t' + str(specification['time_index']) + '.npz'
                archived = array(output / reference_path, 'vector', (context.dimension,), sha(output / reference_path))
                require(np.max(abs(archived - exact_cache[key])) <= STATE_TOLERANCE, 'Independent full-H reference differs')
                reference_files.add(reference_path)
            values = circuit_reference.metrics(exact_cache[key], state, context.model)
            state_files.add(stored['state_file'])
            for key, value in values.items():
                max_metric = max(max_metric, abs(value - stored[key]))
                require(abs(value - stored[key]) <= METRIC_TOLERANCE, 'Independent physical metric differs')
            errors[name] = dict(status='complete', **values)
            counts['complete'] += 1
        else:
            require(candidate['status'] == 'censored' and candidate['compiled_resources'] is None and
                    'state_file' not in stored and 'error' not in stored, 'Censored circuit promoted to physical evidence')
            _, costs = context.audit_compiled_trace(output / candidate['trace_file'], candidate['identity'], VERSION,
                                                    simulate=False, allow_partial=True)
            require(costs == candidate['partial_compiled_resources'], 'Partial trace resources differ')
            errors[name] = dict(status='censored')
            counts['censored'] += 1
        physical[name] = stored
        write_candidate_progress(receipt, name, errors[name])
    for case in cases:
        for label in STATES:
            for ti, T in enumerate(times):
                key = (case, label, ti)
                if key not in exact_cache:
                    context = contexts[case, label]
                    exact_cache[key] = context.exact_action(T)
                    name = 'exact_references/' + case + '.' + label + '.t' + str(ti) + '.npz'
                    archived = array(output / name, 'vector', (context.dimension,), sha(output / name))
                    require(np.max(abs(archived - exact_cache[key])) <= STATE_TOLERANCE, 'Censored-group full-H reference differs')
                    reference_files.add(name)
    groups = defaultdict(dict)
    for index, specification in enumerate(logical):
        row = read(output / 'rows' / ('logical-' + str(index).zfill(4) + '.json'))
        require(all(row[key] == value for key, value in specification.items()), 'Logical arm identity differs')
        for key, value in physical[specification['candidate']].items():
            if key not in specification:
                require(row[key] == value, 'Logical row differs from its shared physical circuit')
        require(row['correlation'] == references[specification['case']]['states'][specification['state']]['metadata'], 'Reference labels differ')
        groups[specification['case'], specification['state'], specification['time_index'], specification['cap_multiple']][specification['arm']] = errors[specification['candidate']]
    independent = []
    for (case, state, ti, cap), arms in groups.items():
        require(set(arms) == set(ARMS), 'Missing logical control')
        for a, b in itertools.combinations(ARMS, 2):
            independent.append(dict(case=case, state=state, T=times[ti], cap_multiple=cap, a=a, b=b, **decision(arms[a], arms[b])))
    actual_pairs = read(output / 'COMPARISONS.json')
    require(len(actual_pairs) == len(independent) == len(cases) * 5 * len(times) * 2 * 28, 'Comparison cohort incomplete')
    for actual, expected in zip(actual_pairs, independent):
        require(set(actual) == set(expected), 'Comparison fields differ')
        for key, value in expected.items():
            if key in ('a_error', 'b_error', 'a_minus_b') and value is not None:
                require(abs(actual[key] - value) <= METRIC_TOLERANCE, 'Comparison error witness differs')
            else:
                require(actual[key] == value, 'Independent pair classification differs')
    costs = read(output / 'PREPARE_COSTS.json')
    for field, value in dict(actual_outer_training_CPU_seconds=sum(row['caller_CPU_seconds'] for row in records.values()),
            actual_training_warmup_CPU_seconds=warm_ordering, actual_compilation_warmup_CPU_seconds=warm_compilation,
            actual_compilation_CPU_seconds=compilation_cost).items():
        near(costs[field], value, absolute=1e-8, relative=0)
    require(not costs['repeated_standalone_charges_are_campaign_work'] and
            costs['original_reference_preparation'] == read(reference_root / 'PREPARATION_COSTS.json'),
            'Actual campaign work double-counts standalone/reference charges')
    constructor_files = {'constructors/' + row['constructor'] + '.json' for row in training}
    expected_frozen = ({'PROTOCOL_FROZEN.json', 'PLAN_FROZEN.json', 'PRACTICAL_INTERFACE.json', 'ORDERS_FROZEN.json',
                        'PREPARE_COSTS.json'} | constructor_files | warm_trace_files |
                       {'warmups/' + case + '.json' for case in cases} |
                       {'reference_evidence/' + name for name in inventory(reference_root)} |
                       {'source_snapshot/' + name for name in protocol['sources']} |
                       {'candidates/' + row['candidate'] + '.json' for row in candidates} |
                       {'traces/' + row['candidate'] + '.ndjson.gz' for row in candidates})
    require(set(frozen['artifacts']) == expected_frozen, 'Frozen complete payload inventory differs')
    order_freeze = read(output / 'ORDERS_FROZEN.json')
    require(order_freeze['constructors'] == len(training) and order_freeze['artifacts'] ==
            {name: sha(output / name) for name in sorted(constructor_files)}, 'Order freeze differs')
    post_freeze = ({'RUN_STARTED.json', 'RUN_COSTS.json', 'RUN_COMPLETE.json', 'COMPARISONS.json'} |
                   {'physical/' + row['candidate'] + '.json' for row in candidates} |
                   {'rows/logical-' + str(index).zfill(4) + '.json' for index in range(len(logical))} |
                   state_files | reference_files)
    require(inventory(output) == expected_frozen | post_freeze | {'SCHEDULES_FROZEN.json'},
            'Unexpected, missing or omitted physical payload')
    require(read(output / 'RUN_STARTED.json')['schedule_freeze_sha256'] == sha(output / 'SCHEDULES_FROZEN.json'),
            'Evaluation did not bind the frozen schedules')
    completion = read(output / 'RUN_COMPLETE.json')
    require(completion['status'] == 'producer_complete_independent_audit_pending' and
            completion['unique_complete'] == counts['complete'] and completion['unique_censored'] == counts['censored'] and
            completion['logical_circuits'] == len(logical) and completion['comparisons'] == len(independent) and
            completion['schedule_freeze_sha256'] == sha(output / 'SCHEDULES_FROZEN.json'), 'Producer completion cohort differs')
    write(receipt / 'INDEPENDENT_ERRORS.json', errors)
    write(receipt / 'INDEPENDENT_COMPARISONS.json', independent)
    report = dict(status='passed', producer_sources=protocol['sources'], producer_source_commit=protocol['source_commit'],
        current_verifier_sources=verifier_sources(), environment=dict(python=sys.version,
            packages={name: importlib.metadata.version(name) for name in ('numpy', 'scipy')}),
        input_reference_audit_status=input_audit['status'], unique_circuits=len(candidates), logical_circuits=len(logical),
        comparisons=len(independent), reference_evolutions=len(exact_cache), max_state_difference=max_state, max_metric_difference=max_metric, **counts,
        schedule_freeze_sha256=sha(output / 'SCHEDULES_FROZEN.json'),
        independent_errors_sha256=sha(receipt / 'INDEPENDENT_ERRORS.json'),
        independent_comparisons_sha256=sha(receipt / 'INDEPENDENT_COMPARISONS.json'),
        physical_artifacts={name: sha(output / name) for name in sorted(post_freeze)},
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
        limitation='Floating numerical audit;no Lean proof,rigorous interval,clock certification or hardware measurement')
    write(receipt / 'VERIFIED.json', report)
    print('INDEPENDENT_TRANSFER_PASSED', len(logical), 'logical circuits;', len(independent), 'pairs', flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt-dir', type=Path, required=True)
    args = parser.parse_args()
    verify(args.output.resolve(), args.receipt_dir.resolve())
