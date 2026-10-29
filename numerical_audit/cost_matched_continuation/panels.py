"""Assemble the complete physical panel from authenticated frozen constructors.

Only bounded serialization and identity/cost checks run here. All constructors
must be present and scored before any panel is written; no physical error is read.
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

import common as k
from admission import execution_admission
from provenance import bind_worker, raw_cpu, relative

PROTOCOL_PATH = 'source_snapshot/numerical_audit/cost_matched_continuation/PROTOCOL.json'
SOURCE_PATH = 'INPUT_PROTOCOL_FROZEN.json'
KIND = 'diagonal_reference_cubic_response'
MAX_JSON = 32 * 1024**2
MAX_MODEL = 8 * 1024**2
COST_FIELDS = ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds')
HASH = re.compile(r'[0-9a-f]{64}')
UNINDEXED = object()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, 'Duplicate JSON key: '+key)
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError('Nonfinite JSON constant: '+value)


def digest(value):
    require(isinstance(value, str) and HASH.fullmatch(value), 'SHA-256 spelling')
    return value


def number(value):
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
            'Finite nonnegative cost/score required')
    return value


def label(value):
    return isinstance(value, str) and 1 <= len(value) <= 200 and all(' ' <= c <= '~' for c in value)


class Artifacts:
    def __init__(self, root):
        self.root = Path(root).resolve(strict=True)
        require(self.root.is_dir(), 'Retained run root must be a directory')
        self.digests = {}

    def path(self, spelling):
        require(isinstance(spelling, str) and 1 <= len(spelling) <= 1024 and '\\' not in spelling,
                'Bounded POSIX artifact path')
        relative = PurePosixPath(spelling)
        require(not relative.is_absolute() and all(part not in ('', '.', '..') for part in relative.parts)
                and str(relative) == spelling, 'Canonical relative artifact path')
        path = (self.root / spelling).resolve(strict=True)
        require(path.is_relative_to(self.root) and path.is_file(), 'Artifact escapes retained run or is not a file')
        return path

    def raw(self, spelling, expected=UNINDEXED, maximum=MAX_JSON):
        path = self.path(spelling)
        require(path.stat().st_size <= maximum, 'Artifact byte cap: '+spelling)
        with path.open('rb') as stream:
            raw = stream.read(maximum+1)
        require(len(raw) <= maximum, 'Growing artifact byte cap: '+spelling)
        actual = hashlib.sha256(raw).hexdigest()
        require(expected is UNINDEXED or digest(expected) == actual, 'Artifact hash differs: '+spelling)
        self.digests[spelling] = actual
        return raw

    def load(self, spelling, expected=UNINDEXED, maximum=MAX_JSON):
        return json.loads(self.raw(spelling, expected, maximum), object_pairs_hook=unique_object,
                          parse_constant=reject_constant)


def authenticate_model(artifacts, case, declared):
    require(isinstance(declared, dict) and {'model.json', 'model.toml'} <= set(declared) and
            set(declared) <= {'model.json', 'model.toml', 'tensors.npz', 'orbitals.npz'},
            'Declared original input file inventory')
    for name, expected in declared.items():
        artifacts.raw('inputs/'+case+'/'+name, expected, MAX_MODEL if name.startswith('model.') else MAX_JSON)
    json_path = 'inputs/'+case+'/model.json'
    model = artifacts.load(json_path, declared['model.json'], MAX_MODEL)
    require(model['case'] == case and type(model['n']) is int and model['n'] == 14 and
            type(model['electrons']) is int and model['electrons'] == 4, 'Original model identity')
    k.prior.blocks(model)
    toml_path = 'inputs/'+case+'/model.toml'
    target = tomllib.loads(artifacts.raw(toml_path, declared['model.toml'], MAX_MODEL).decode('utf-8'))
    require(target['case'] == case and type(target['n']) is int and target['n'] == 14 and
            type(target['electrons']) is int and target['electrons'] == 4 and
            target['partition'] == 'conserving' and target['denominator'] == model['denominator'] and
            target['source_json_sha256'] == declared['model.json'], 'TOML exact original JSON identity')
    units = target['units']
    require(isinstance(units, list) and len(units) == len(model['blocks']) and
            target['baseline_order'] == list(range(1, len(units)+1)) and
            all(type(value) is int for value in target['baseline_order']), 'Original TOML partition/order')
    for raw, unit in zip(model['blocks'], units):
        entries = sorted(k.prior.decode(raw).items())
        words = [','.join(str(site)+letter for site, letter in word) for word, _ in entries]
        require(unit['words'] == words and unit['numerators'] == [str(value) for _, value in entries],
                'TOML original word/coefficient partition differs')
    return model


def cost_evidence(artifacts, cost, record, parent, worker, config):
    name = cost['label']
    require(label(name) and '/' not in name, 'Bounded owning-job label')
    require(cost['final_path'] == 'campaign/jobs/'+name+'.FINAL.json' and
            cost['receipt_path'] == 'campaign/logs/'+name+'.resources.json', 'Canonical owning-job evidence paths')
    final = artifacts.load(cost['final_path'], cost['final_sha256'])
    receipt = artifacts.load(cost['receipt_path'], cost['receipt_sha256'])
    plan_path = 'campaign/jobs/'+name+'.PLAN.json'
    if 'plan_path' in cost or 'plan_sha256' in cost:
        require(cost.get('plan_path') == plan_path and 'plan_sha256' in cost, 'Canonical owning PLAN pin pair')
        plan = artifacts.load(plan_path, cost['plan_sha256'])
    else:
        # Original records omitted both PLAN fields; the independent ledger
        # still authenticates that PLAN and its complete command grammar.
        plan = artifacts.load(plan_path)
    require(final['phase'] == 'constructors' and final['label'] == name and
            final['receipt_sha256'] == cost['receipt_sha256'], 'Owning-job FINAL/receipt identity')
    require(cost['raw_CPU_seconds'] == raw_cpu(final, receipt) and
            cost['charged_CPU_seconds'] == final['charged_CPU_seconds'] and
            cost['wall_seconds'] == final['absolute_guarded_wall_seconds'], 'Canonical costs differ from OS evidence')
    require(final['exit_code'] == receipt['exit_code'] == record['guard_exit_code'] and
            cost['wall_seconds'] >= number(receipt['elapsed']), 'Guard exit/wall evidence differs')
    binding = bind_worker(plan, receipt, record, parent, worker, config)
    require(record.get('cost_binding', binding) == binding, 'Canonical cost binding differs from executed PLAN')
    if worker is not None:
        require(cost['charged_CPU_seconds'] >= number(worker['process_CPU_before_output']), 'CPU charge below worker observation')


def cost_check(artifacts, record, parent, config, worker):
    costs = record['costs']
    require(isinstance(costs, dict) and all(boundary in costs for boundary in ('prefix', 'continuation', 'total')),
            'Complete prefix/continuation/total cost record')
    for boundary in ('prefix', 'continuation', 'total'):
        for field in COST_FIELDS:
            if field == 'raw_CPU_seconds' and costs[boundary][field] is None:
                require(record['eligible'] is False, 'Unknown raw CPU cannot be eligible')
            else:
                number(costs[boundary][field])
        require(costs[boundary]['raw_CPU_seconds'] is None or
                costs[boundary]['charged_CPU_seconds'] >= costs[boundary]['raw_CPU_seconds'],
                'Conservative CPU below raw CPU')
    for field in COST_FIELDS:
        require(costs['prefix'][field] == parent['costs']['total'][field], 'Measured shared parent prefix differs')
        if costs['prefix'][field] is None or costs['continuation'][field] is None:
            require(field == 'raw_CPU_seconds' and costs['total'][field] is None, 'Unknown raw total differs')
        else:
            require(math.isclose(costs['total'][field], costs['prefix'][field]+costs['continuation'][field],
                                 rel_tol=0., abs_tol=1e-9), 'Prefix/continuation total differs')
    parent_arm = record['arm'] == 'parent'
    cap = config['parent']['CPU_seconds'] if parent_arm else record['total_CPU_budget']
    if record['eligible']:
        require(costs['total']['raw_CPU_seconds'] <= cap and costs['total']['charged_CPU_seconds'] <= cap,
                'Strict complete construction CPU ceiling')
    if parent_arm:
        require(all(costs['continuation'][field] == 0 for field in COST_FIELDS), 'Parent continuation must be zero')
        wall_cap = config['parent']['wall_seconds']
    else:
        require(costs['prefix'] == parent['costs']['prefix'], 'Frozen parent prefix receipt references differ')
        wall_cap = config['accounting']['total_wall_multiplier']*cap
    if record['eligible']:
        require(costs['total']['wall_seconds'] <= wall_cap, 'Strict complete construction wall ceiling')
    cost_evidence(artifacts, costs['prefix'] if parent_arm else costs['continuation'], record, parent, worker, config)


def parent_preparation(record):
    """Verify the measured witness structurally; never generate new starts."""
    plan = record.get('restarts')
    if plan is None:
        return
    material = (record['case']+':'+str(record['repetition'])+':cost-matched-v1').encode('utf-8')
    seed = int.from_bytes(hashlib.sha256(material).digest()[:8], 'big')
    require(isinstance(plan, dict) and type(plan.get('seed')) is int and plan['seed'] == seed and
            plan.get('interior_perturbation_counts') == [0, 1, 2, 4] and
            all(type(value) is int for value in plan['interior_perturbation_counts']), 'Frozen parent restart seed/counts')
    orders = plan.get('orders'); order = record['order']
    require(isinstance(orders, list) and len(orders) == 4, 'Four measured parent starts')
    for candidate in orders:
        k.validate_order(candidate, len(order))
        require(candidate[0] == order[0] and candidate[-1] == order[-1], 'Interior frozen parent starts')
    require(orders[0] == order and len({tuple(value) for value in orders}) == 4, 'Distinct frozen parent starts')


def worker_witness(artifacts, record, parent):
    spelling, expected = record['worker_output_path'], record['worker_output_sha256']
    relative(spelling)
    path = (artifacts.root/spelling).resolve()
    require(path.is_relative_to(artifacts.root), 'Worker path outside archive')
    fallback = record['controller_fallback']
    if expected is None:
        require(fallback and not path.exists() and record['result'] is None, 'Missing worker contradicts canonical result')
        worker = None
    else:
        worker = artifacts.load(spelling, expected)
        keys = ('case', 'repetition', 'arm', 'total_CPU_budget', 'protocol_sha256',
                'source_freeze_sha256', 'parent_sha256', 'model_json_sha256', 'model_toml_sha256')
        require(worker['schema'] == 'cost-matched-worker-v1' and
                all(worker[key] == record[key] for key in keys) and worker['result'] == record['result'], 'Worker identity/result differs')
        if not fallback:
            require(worker['result']['order'] == record['order'] and worker['result']['score'] == record['score'],
                    'Selected incumbent differs from worker')
            if record['arm'] == 'parent':
                require(worker.get('restarts') == record.get('restarts'), 'Measured restart preparation differs')
                parent_preparation(record)
            else:
                expected_starts = parent['restarts']['orders'] if record['arm'] == 'full' else [parent['order']]
                require(worker['starts'] == expected_starts, 'Continuation input starts differ from frozen parent')
    if fallback:
        require(record['eligible'] is False and (record['result'] is None or record['result'].get('score') is None),
                'Controller fallback contradicts a scored/eligible result')
        if record['arm'] == 'parent':
            require(record['score'] is None and record['order'] == list(range(len(record['order']))), 'Unscored parent fallback differs')
        else:
            require(record['order'] == parent['order'] and record['score'] == parent['score'], 'Controller fallback is not the frozen parent')
    method = execution_admission(record['arm'], None if worker is None else worker['result'],
        parent_order=None if record['arm'] == 'parent' else parent['order'],
        refusal=dict(type='missing_worker') if worker is None else worker.get('refusal'))
    if record['arm'] == 'parent' and record.get('restarts') is None:
        method['admitted'] = False
    if record['arm'] != 'parent' and not parent['eligible']:
        method['admitted'] = False
    require(record['execution_admission']['admitted'] is method['admitted'], 'Method execution admission differs')
    if record['eligible']:
        require(not fallback and record['guard_exit_code'] == 0 and worker is not None and method['admitted'],
                'Eligible constructor has an unverified execution/fallback')
    return worker


def constructor_identity(record):
    return record['case'], record['repetition'], record['arm'], record['total_CPU_budget']


def authenticate(artifacts):
    source = artifacts.load(SOURCE_PATH)
    source_sha = artifacts.digests[SOURCE_PATH]
    protocol = artifacts.load(PROTOCOL_PATH, source['configuration_sha256'])
    protocol_sha = artifacts.digests[PROTOCOL_PATH]
    k.validate_protocol(protocol)
    require(source['configuration'] == protocol and source['sources'] == k.source_map(),
            'Archived protocol or current producer/dependency sources differ from freeze')
    for name, expected in source['sources'].items():
        artifacts.raw('source_snapshot/'+name, expected)
    require(source['sources'].get('numerical_audit/cost_matched_continuation/PROTOCOL.json') == protocol_sha,
            'Protocol bytes absent from source freeze')
    inputs = artifacts.load('INPUTS_FROZEN.json')
    require(inputs['status'] == 'complete' and not inputs.get('refusals') and
            inputs['protocol_sha256'] == protocol_sha and inputs['source_freeze_sha256'] == source_sha,
            'Complete source/protocol-bound inputs required')
    parents = artifacts.load('PARENTS_FROZEN.json')
    constructors = artifacts.load('CONSTRUCTORS_FROZEN.json')
    require(parents['schema'] == 'cost-matched-parent-freeze-v1' and
            constructors['schema'] == 'cost-matched-constructor-freeze-v1', 'Exact constructor freeze schemas')
    mode = parents['mode']; require(mode in ('main', 'preflight'), 'Declared panel mode')
    repetitions = protocol['repetitions'] if mode == 'main' else protocol['preflight']['repetitions']
    require(type(repetitions) is int and repetitions == (2 if mode == 'main' else 1), 'Frozen repetition grid')
    cases = parents['cases']
    require(isinstance(cases, list) and all(label(case) and '/' not in case for case in cases) and
            len(cases) == len(set(cases)), 'Unique bounded case roster')
    expected_cases = {row['case'] for row in source['registry']['cases']} if mode == 'main' else {protocol['preflight']['case']}
    require(set(cases) == expected_cases and len(cases) == (2 if mode == 'main' else 1), 'Complete frozen case roster')
    require(set(inputs['successes']) == set(cases) and inputs['expected_cases'] == len(cases), 'Input/campaign case roster')
    for index in (parents, constructors):
        require(all(index[key] == value for key, value in dict(mode=mode, repetitions=repetitions, cases=cases,
            protocol_sha256=protocol_sha, source_freeze_sha256=source_sha,
            protocol_path=PROTOCOL_PATH, source_freeze_path=SOURCE_PATH).items()), 'Frozen constructor index identity')
        require(isinstance(index['entries'], list), 'Frozen constructor entries required')
    models = {case: authenticate_model(artifacts, case, inputs['successes'][case]) for case in cases}
    expected_parents = {(case, rep, 'parent', 0) for case in cases for rep in range(1, repetitions+1)}
    expected_children = {(case, rep, arm, budget) for case in cases for rep in range(1, repetitions+1)
                         for arm in protocol['arms'] for budget in protocol['total_CPU_budgets']}
    records, hashes = {}, {}
    for index, expected, parent_arm in ((parents, expected_parents, True), (constructors, expected_children, False)):
        require(len(index['entries']) == len(expected), 'Complete constructor count before assembly')
        observed = set()
        for entry in index['entries']:
            require(isinstance(entry, dict), 'Constructor descriptor required')
            case, rep = entry['case'], entry['repetition']
            arm = 'parent' if parent_arm else entry['arm']
            budget = 0 if parent_arm else entry['total_CPU_budget']
            require(type(rep) is int and type(budget) is int, 'Integer constructor identity')
            key = case, rep, arm, budget
            require(key in expected and key not in observed, 'Unknown/duplicate frozen constructor')
            observed.add(key)
            canonical = ('parents/'+case+'/rep'+str(rep)+'.json' if parent_arm else
                         'orders/'+case+'/rep'+str(rep)+'/'+arm+'_B'+str(budget)+'.json')
            require(entry['path'] == canonical, 'Canonical constructor path differs')
            record = artifacts.load(canonical, entry['sha256'])
            require(record['schema'] == 'cost-matched-constructor-v1' and constructor_identity(record) == key,
                    'Canonical constructor/index identity')
            require(record['protocol_sha256'] == protocol_sha and record['source_freeze_sha256'] == source_sha,
                    'Constructor source/protocol identity')
            require(record['model_json_sha256'] == inputs['successes'][case]['model.json'] and
                    record['model_toml_sha256'] == inputs['successes'][case]['model.toml'], 'Constructor input bytes differ')
            k.validate_order(record['order'], len(models[case]['blocks']))
            score = record['score']
            require(isinstance(score, dict), 'Explicit refusal: unscored constructor '+canonical)
            number(score['value']); number(score['cubic_response'])
            require(type(score['steps']) is int and score['steps'] > 0, 'Scored affordable steps')
            require(type(record['eligible']) is bool and type(record['controller_fallback']) is bool and
                    type(record['guard_exit_code']) is int, 'Explicit constructor eligibility/guard status')
            if parent_arm:
                parent = record
            else:
                parent_key = case, rep, 'parent', 0
                parent = records[parent_key]
                require(record['parent_sha256'] == hashes[parent_key], 'Continuation frozen-parent hash differs')
            worker = worker_witness(artifacts, record, parent)
            cost_check(artifacts, record, parent, protocol, worker)
            records[key] = record; hashes[key] = digest(entry['sha256'])
        require(observed == expected, 'Missing frozen constructor identity')
    return dict(protocol=protocol, protocol_sha256=protocol_sha, source_freeze_sha256=source_sha,
                mode=mode, repetitions=repetitions, cases=cases, inputs=inputs, records=records, hashes=hashes)


def toml_value(value):
    if type(value) is int:
        return str(value)
    if type(value) is float:
        require(math.isfinite(value), 'Finite TOML scalar')
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    if isinstance(value, list):
        return '['+','.join(toml_value(item) for item in value)+']'
    raise ValueError('Unsupported panel TOML scalar/array')


def panel_bytes(header, rows):
    lines = [key+'='+toml_value(value) for key, value in header.items()]
    for row in rows:
        lines.append('[[candidates]]')
        lines.extend(key+'='+toml_value(value) for key, value in row.items())
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= MAX_MODEL, 'Candidate TOML byte cap')
    return raw


def assemble(root):
    artifacts = Artifacts(root)
    admitted = authenticate(artifacts)
    protocol = admitted['protocol']; reps = admitted['repetitions']
    staged = []
    for case in admitted['cases']:
        rows = []
        for rep in range(1, reps+1):
            parent_sha = admitted['hashes'][case, rep, 'parent', 0]
            for arm, budget in [('parent', 0)]+[(arm, budget) for arm in protocol['arms'] for budget in protocol['total_CPU_budgets']]:
                key = case, rep, arm, budget
                record = admitted['records'][key]
                oid = case+'-rep'+str(rep)+'-'+arm+'-B'+str(budget)
                require(label(oid), 'Bounded constructor ID')
                for duration in protocol['physical']['times']:
                    for cap in protocol['physical']['rotation_cap_multipliers']:
                        name = oid+'-T'+format(duration, 'g')+'-C'+str(cap)
                        require(label(name), 'Bounded candidate name')
                        score = float(record['score']['value'])
                        rows.append(dict(name=name, order_id=oid, arm=arm, repetition=rep, cpu_budget=budget,
                            partition='conserving', formula=2, time=float(duration), budget_multiplier=cap,
                            order=[value+1 for value in record['order']], constructor_sha256=admitted['hashes'][key],
                            parent_sha256=parent_sha, claimed_score=score, constructor_score=score, score_kind=KIND))
        require(len(rows) == 20*reps and len({row['name'] for row in rows}) == len(rows), 'Complete unique panel rows')
        header = dict(version=1, input_sha256=admitted['inputs']['successes'][case]['model.toml'],
            source_json_sha256=admitted['inputs']['successes'][case]['model.json'],
            protocol_sha256=admitted['protocol_sha256'], source_freeze_sha256=admitted['source_freeze_sha256'],
            case=case, mode=admitted['mode'], repetitions=reps, expected_logical_rows=len(rows))
        staged.append((case, 'panels/'+case+'/candidates.toml', panel_bytes(header, rows)))
    panel_root = artifacts.root/'panels'; index_path = artifacts.root/'PANELS_FROZEN.json'
    require(not panel_root.exists() and not index_path.exists() and not index_path.with_name(index_path.name+'.partial').exists(),
            'Fresh panel outputs required; preserve all earlier attempts')
    # Every input, constructor and output byte sequence was checked before this first write.
    entries = []
    for case, spelling, raw in staged:
        path = artifacts.root/spelling; path.parent.mkdir(parents=True, exist_ok=False)
        temporary = path.with_name(path.name+'.partial')
        with temporary.open('xb') as stream:
            stream.write(raw)
        os.link(temporary, path); temporary.unlink()
        entries.append(dict(case=case, path=spelling, sha256=hashlib.sha256(raw).hexdigest(),
                            expected_logical_rows=20*reps))
    frozen = dict(schema='cost-matched-panels-frozen-v1', mode=admitted['mode'], repetitions=reps,
        cases=admitted['cases'], protocol_sha256=admitted['protocol_sha256'],
        source_freeze_sha256=admitted['source_freeze_sha256'], protocol_path=PROTOCOL_PATH,
        source_freeze_path=SOURCE_PATH, entries=entries,
        input_freeze_sha256=artifacts.digests['INPUTS_FROZEN.json'],
        parents_freeze_sha256=artifacts.digests['PARENTS_FROZEN.json'],
        constructors_freeze_sha256=artifacts.digests['CONSTRUCTORS_FROZEN.json'],
        physical_error_read=False)
    k.write_new(index_path, frozen)
    return frozen


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args(argv)
    assemble(args.root)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
