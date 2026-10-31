"""Freeze panels, then independently join emitted manifests to original bytes.

This consumer performs serialization, exact-integer model and resource metadata
checks only. It never computes a training score or a physical state.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import random
import re
import sys
import tomllib

if __package__:
    from . import common as k
else:
    import common as k

HASH = re.compile(r'[0-9a-f]{64}')
COST_FIELDS = ('raw_CPU_seconds', 'charged_CPU_seconds', 'wall_seconds')
IDENTITIES = (('parent', 0), ('canonical', 0), ('magnitude', 0), ('full', 32),
              ('full', 64), ('window', 32), ('window', 64), ('annealing', 32), ('annealing', 64))
HEADER_FIELDS = ('version', 'order_index_base', 'case', 'n', 'electrons', 'partition',
    'input_sha256', 'source_json_sha256', 'protocol_sha256', 'source_freeze_sha256', 'mode',
    'repetitions', 'expected_logical_rows', 'training_time', 'training_budget_multiplier',
    'times', 'rotation_cap_multipliers', 'protocol_path', 'source_freeze_path')
ROW_FIELDS = ('name', 'order_id', 'arm', 'repetition', 'cpu_budget', 'partition', 'formula',
    'time', 'budget_multiplier', 'order', 'constructor_sha256', 'claimed_score', 'score_kind', 'eligible')
require = k.require


def number(value):
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
            'Finite nonnegative numeric scalar required')
    return value


def digest(value):
    require(isinstance(value, str) and HASH.fullmatch(value) is not None, 'SHA-256 spelling')
    return value


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, 'Duplicate JSON key: '+key)
        result[key] = value
    return result


class Artifacts:
    def __init__(self, root):
        self.root = Path(root).resolve(strict=True)
        require(self.root.is_dir(), 'DATA must be a directory')
        self.hashes = {}

    def path(self, relative):
        require(isinstance(relative, str) and 0 < len(relative) <= 1024 and '\\' not in relative,
                'Bounded POSIX relative artifact path')
        name = PurePosixPath(relative)
        require(not name.is_absolute() and str(name) == relative and
                all(part not in ('', '.', '..') for part in relative.split('/')), 'Canonical relative artifact path')
        lexical = self.root/relative
        path = lexical.resolve(strict=True)
        require(path.is_relative_to(self.root) and path.is_file() and path == lexical,
                'Artifact must be a canonical contained real file')
        return path

    def raw(self, relative, expected=None, maximum=8*1024**2):
        path = self.path(relative)
        require(path.stat().st_size <= maximum, 'Artifact byte cap: '+relative)
        with path.open('rb') as stream:
            raw = stream.read(maximum+1)
        require(len(raw) <= maximum, 'Growing artifact byte cap: '+relative)
        actual = hashlib.sha256(raw).hexdigest()
        require(expected is None or digest(expected) == actual, 'Actual artifact digest differs: '+relative)
        self.hashes[relative] = actual
        return raw

    def json(self, relative, expected=None):
        def reject(value):
            raise ValueError('Nonfinite JSON constant: '+value)
        return json.loads(self.raw(relative, expected), object_pairs_hook=unique_object, parse_constant=reject)

    def toml(self, relative, expected=None):
        return tomllib.loads(self.raw(relative, expected).decode('utf-8'))

    def reference(self, value):
        require(isinstance(value, dict) and set(value) == {'path', 'sha256'}, 'Exact artifact reference fields')
        return self.json(value['path'], digest(value['sha256']))


def context(data, *, producer=False, preflight=False):
    artifacts = Artifacts(data)
    freeze = k.frozen(artifacts.root, producer=producer)
    require(freeze['schema'] == 'molecular-source-freeze-v1', 'Immutable source freeze schema')
    info, config = k.header(artifacts.root, producer=producer)
    require((info['mode'] == 'preflight') == preflight, 'Explicit --preflight must match immutable freeze')
    require(config['total_CPU_budgets'] == [32, 64] and config['repetitions'] == 2 and
            config['optimized_arms'] == ['full', 'window', 'annealing'] and
            config['cheap_arms'] == ['canonical', 'magnitude'] and
            config['training_time'] == 1.0 and config['training_rotation_multiplier'] == 30 and
            config['physical']['times'] == [.65, 1.35] and
            config['physical']['rotation_cap_multipliers'] == [8, 24] and
            config['physical']['formula'] == 2, 'Frozen scientific grid')
    if not preflight:
        require(config['main_execution_allowed'] is True and config['protocol_state'] == 'frozen_main',
                'Main execution disabled until old-only preflight and final freeze')
        gate = artifacts.json('source_snapshot/numerical_audit/molecular_generalization/PREFLIGHT_ACCEPTED.json')
        calibration = artifacts.json('source_snapshot/numerical_audit/molecular_generalization/CALIBRATION.json')
        require(gate['complete'] is True and gate['new_physical_errors_observed'] is False and
                calibration['complete'] is True and config['annealing']['calibration_sha256'] ==
                artifacts.hashes['source_snapshot/numerical_audit/molecular_generalization/CALIBRATION.json'],
                'Authenticated old-only prerequisites')
        require(len(info['cases']) == 2 and len(set(info['cases'])) == 2, 'Exactly two frozen primary cases')
    else:
        require(info['cases'] == ['LiH_r1p487_631g_14m4e'] and info['repetitions'] == 1,
                'Only the declared old-case preflight')
    inputs = artifacts.json('INPUTS_FROZEN.json')
    require(inputs['status'] == 'complete' and not inputs['refusals'] and
            set(inputs['successes']) == set(info['cases']), 'Complete exact input case inventory')
    require(all(inputs[key] == info[key] for key in info), 'Input freeze header joins')
    indexed_paths = set(); owned_plans = set()
    for filename, parent_only in (('PARENTS_FROZEN.json', True), ('CONSTRUCTORS_FROZEN.json', False)):
        index = artifacts.json(filename)
        require(all(index[key] == info[key] for key in info), 'Frozen constructor-index header')
        expected = {(case, rep, arm, budget) for case in info['cases'] for rep in range(1, info['repetitions']+1)
                    for arm, budget in IDENTITIES if (arm == 'parent') == parent_only}
        observed = set()
        for entry in index['entries']:
            identity = tuple(entry[key] for key in ('case', 'repetition', 'arm', 'cpu_budget'))
            require(identity in expected and identity not in observed, 'Frozen constructor-index exact identity')
            observed.add(identity)
            case, rep, arm, budget = identity
            require(entry['path'] == f'constructors/{case}.rep{rep}.{arm}.B{budget}.json', 'Canonical indexed constructor path')
            owner = artifacts.json(entry['path'], digest(entry['sha256']))
            indexed_paths.add(entry['path']); owned_plans.add(owner['own_plan']['path'])
        require(observed == expected, 'Complete independent constructor-index grid')
    require({str(path.relative_to(artifacts.root)) for path in (artifacts.root/'constructors').glob('*.json')} == indexed_paths,
            'No missing/extra unindexed constructor JSON')
    constructor_plans = {str(path.relative_to(artifacts.root)) for path in (artifacts.root/'jobs').glob('*.PLAN.json')
                         if artifacts.json(str(path.relative_to(artifacts.root))).get('phase') == 'constructors'}
    require(constructor_plans == owned_plans and len(owned_plans) == len(indexed_paths),
            'Every constructor attempt has one distinct retained owning PLAN; no omitted retry costs')
    return artifacts, freeze, info, config, inputs


def _conservation(unit, *, alpha=False):
    polynomial = {}
    for word, numerator in unit.items():
        for site, letter in word:
            if letter not in ('X', 'Y') or (alpha and site % 2):
                continue
            changed = tuple((q, ('Y' if p == 'X' else 'X') if q == site else p) for q, p in word)
            value = polynomial.get(changed, 0) + (-numerator if letter == 'X' else numerator)
            if value:
                polynomial[changed] = value
            else:
                polynomial.pop(changed, None)
    return polynomial


def authenticate_model(a, case, declared, spec):
    require(isinstance(declared, dict) and {'model.json', 'model.toml'} <= set(declared) <=
            {'model.json', 'model.toml', 'orbitals.npz', 'tensors.npz'}, 'Exact input file inventory')
    for name, sha in declared.items():
        a.raw(f'inputs/{case}/{name}', digest(sha), maximum=32*1024**2 if name.endswith('.npz') else 8*1024**2)
    model = a.json(f'inputs/{case}/model.json', declared['model.json'])
    n, electrons = model['n'], model['electrons']
    require(model['case'] == case and type(n) is int and type(electrons) is int and
            (n, electrons) == spec and 1 <= n <= 14 and 0 < electrons < n and math.comb(n, electrons) <= 2000,
            'Frozen case dimensions/full-register/complete fixed-N domain')
    den = model['denominator']
    require(isinstance(den, str) and den.isascii() and den.isdigit() and 0 < len(den) <= 1280 and
            0 < int(den) and int(den).bit_length() <= 4096, 'Exact decimal denominator')
    total = k.decode(model['coefficients'])
    units = [k.decode(block) for block in model['blocks']]
    require(0 < len(total) <= 8192 and 0 < len(units) <= 1024 and all(units) and
            sum(map(len, units)) == len(total) and {w: v for u in units for w, v in u.items()} == total,
            'Complete exact original coefficient partition')
    for unit in units:
        require(not _conservation(unit) and not _conservation(unit, alpha=True), 'Exact unit N/alpha conservation')
        masks = []
        for word in unit:
            require(all(0 <= q < n for q, _ in word) and sum(p == 'Y' for _, p in word) % 2 == 0,
                    'Real original Pauli words within register')
            x = sum(1 << q for q, p in word if p in ('X', 'Y'))
            z = sum(1 << q for q, p in word if p in ('Z', 'Y'))
            masks.append((x, z))
        require(len({x for x, _ in masks}) == 1 and all(((x & zz).bit_count()+(xx & z).bit_count()) % 2 == 0
                for x, z in masks for xx, zz in masks), 'Equal-flip internally commuting original unit')
    target = a.toml(f'inputs/{case}/model.toml', declared['model.toml'])
    require(all(target[key] == model[key] for key in ('case', 'n', 'electrons', 'denominator')) and
            target['partition'] == 'conserving' and target['source_json_sha256'] == declared['model.json'] and
            target['baseline_order'] == list(range(1, len(units)+1)) and
            all(type(v) is int for v in target['baseline_order']) and len(target['units']) == len(units),
            'Exact model JSON/TOML identity and original one-based partition')
    for unit, row in zip(units, target['units']):
        terms = sorted(unit.items())
        require(set(row) == {'words', 'numerators'} and
                row['words'] == [','.join(str(q)+p for q, p in w) for w, _ in terms] and
                row['numerators'] == [str(v) for _, v in terms], 'Exact original model word/numerator join')
    return model, units


def _sampled_cpu(receipt):
    values = []
    for sample in receipt.get('samples', []):
        service = sample.get('service', '')
        require(isinstance(service, str), 'Resource sample service text')
        for raw in re.findall(r'^CPUUsageNSec=([0-9]+)$', service, re.MULTILINE):
            values.append(int(raw)/1e9)
    return values


def measured_receipt(receipt):
    process = receipt.get('process_measurements', {})
    user, system = process.get('user_cpu_seconds'), process.get('system_cpu_seconds')
    raw = number(user)+number(system) if user is not None and system is not None else None
    charged = max([raw+.011]+_sampled_cpu(receipt)) if raw is not None else None
    return dict(raw_CPU_seconds=raw, charged_CPU_seconds=charged, wall_seconds=number(receipt['elapsed']))


def _worker_arguments(command, plan, record):
    require(isinstance(command, list) and 3 <= len(command) <= 40 and
            all(isinstance(v, str) and len(v) <= 4096 and '\x00' not in v for v in command), 'Bounded worker command')
    script = PurePosixPath(command[1]); producer = PurePosixPath(plan['producer_root'])
    require(producer.is_absolute() and str(producer) == plan['producer_root'] and '..' not in producer.parts and
            script == producer/'numerical_audit/molecular_generalization/construct.py', 'Exact original constructor source path')
    flags = {}; index = 2
    while index < len(command):
        key = command[index]
        require(key not in flags, 'Duplicate worker command flag')
        if key == '--worker':
            flags[key] = True; index += 1
        else:
            require(key in ('--data', '--root', '--case', '--repetition', '--arm', '--budget', '--raw',
                            '--cpu-cap', '--wall-cap', '--julia-executable') and index+1 < len(command),
                    'Unknown worker command flag')
            flags[key] = command[index+1]; index += 2
    require(flags.get('--worker') is True and
            (flags.get('--case'), int(flags.get('--repetition', '0')), flags.get('--arm'), int(flags.get('--budget', '0'))) ==
            (record['case'], record['repetition'], record['arm'], record['cpu_budget']), 'Executed constructor identity')
    original = flags.get('--data', flags.get('--root'))
    require(isinstance(original, str) and PurePosixPath(original).is_absolute() and '..' not in PurePosixPath(original).parts,
            'Original DATA absolute path')
    require(PurePosixPath(flags['--raw']) == PurePosixPath(original)/plan['expected_raw'] and
            float(flags['--cpu-cap']) == plan['cpu_cap'] and float(flags['--wall-cap']) == plan['wall_cap'],
            'Executed worker output/cap binding')
    return flags


def cost_check(a, record, parent, config):
    costs = record['costs']; arm = record['arm']; optimized = arm in k.ARMS
    require(isinstance(costs, dict) and set(costs) == {'prefix', 'continuation', 'total'}, 'Complete exact cost boundaries')
    for boundary in costs.values():
        require(isinstance(boundary, dict) and set(boundary) == set(COST_FIELDS), 'Complete measured cost fields')
        for field in COST_FIELDS:
            if boundary[field] is None:
                require(field != 'wall_seconds' and record['eligible'] is False, 'Unknown cost remains ineligible')
            else:
                number(boundary[field])
    for field in COST_FIELDS:
        if arm == 'parent':
            require(costs['continuation'][field] == 0, 'Parent has no continuation')
        else:
            prefix = parent['costs']['total'][field] if optimized else 0.
            require(costs['prefix'][field] == prefix, 'Optimized measured parent / cheap zero prefix')
        prefix, continuation = costs['prefix'][field], costs['continuation'][field]
        expected = None if prefix is None or continuation is None else prefix+continuation
        require(costs['total'][field] == expected, 'Exact prefix plus continuation total')
    receipt = a.reference(record['own_receipt']); plan = a.reference(record['own_plan'])
    identity = {key: record[key] for key in ('case', 'repetition', 'arm', 'cpu_budget')}
    require(plan['schema'] == 'molecular-job-plan-v1' and plan['phase'] == 'constructors' and
            plan['identity'] == identity and plan['source_freeze_sha256'] == record['source_freeze_sha256'] and
            plan['command'] == receipt['command'], 'Actual PLAN/receipt/constructor identity')
    require(record['own_plan']['path'] == 'jobs/'+plan['label']+'.PLAN.json' and
            record['own_receipt']['path'] == 'logs/'+plan['label']+'.resources.json', 'Canonical own plan/receipt label')
    digest(plan['guard_sha256'])
    # The guard itself is a frozen source, never an unauthenticated executable.
    guard = 'numerical_audit/molecular_generalization/tools/resource_guard.py'
    a.raw('source_snapshot/'+guard, plan['guard_sha256'])
    raw_path = plan['expected_raw']
    worker = a.json(raw_path, record.get('worker_output_sha256', record.get('raw_sha256')))
    require(record.get('worker_output_path', raw_path) == raw_path, 'Worker output reference identity')
    require(record.get('worker_output_sha256', record.get('raw_sha256')) is not None, 'Actual worker digest required')
    for key in ('case', 'repetition', 'arm', 'cpu_budget', 'model_json_sha256', 'model_toml_sha256',
                'protocol_sha256', 'source_freeze_sha256'):
        require(worker[key] == record[key], 'Worker/final constructor join: '+key)
    require(worker.get('result') == record['result'], 'Retain exact raw backend result')
    require(worker.get('parent_sha256') == record.get('parent_sha256'), 'Worker measured parent join')
    immutable = a.json('INPUT_PROTOCOL_FROZEN.json', record['source_freeze_sha256'])
    require(plan['producer_root'] == immutable['producer_root'] and
            plan['producer_data_root'] == immutable['producer_data_root'], 'Original producer/data roots bound to freeze')
    flags = _worker_arguments(plan['command'], plan, record)
    require(flags.get('--data', flags.get('--root')) == immutable['producer_data_root'], 'Executed original DATA root')
    own_cap = config['parent']['CPU_seconds'] if arm == 'parent' else (
        config['cheap']['CPU_ceiling'] if arm in k.CHEAP else record['cpu_budget']-costs['prefix']['charged_CPU_seconds'])
    wall_cap = config['parent']['wall_seconds'] if arm == 'parent' else (
        config['cheap']['wall_ceiling'] if arm in k.CHEAP else math.floor(5*record['cpu_budget']-costs['prefix']['wall_seconds']))
    require(plan['cpu_cap'] == own_cap and plan['wall_cap'] == wall_cap and
            receipt['budget']['cpu_seconds'] == own_cap and receipt['budget']['wall_seconds'] == wall_cap and
            receipt['budget']['disk_mib'] == plan['disk_mib'], 'Complete measured own-process cap')
    settings = receipt['settings']; environment = receipt['environment']
    require(all(str(settings[key]) == value for key, value in dict(CPUQuota='25%', CPUWeight='5', IOWeight='10',
            Nice='19', MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0').items()), 'Frozen numerical resource profile')
    require(all(environment[key] == '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'One numerical thread')
    own = measured_receipt(receipt)
    require(costs['prefix' if arm == 'parent' else 'continuation'] == own, 'Claimed own cost differs from terminal measured receipt')
    require(own['charged_CPU_seconds'] is None or own['charged_CPU_seconds'] >= number(worker['process_CPU_before_output']), 'OS charge below raw worker CPU')
    require(worker['guarded_CPU_cap'] == own_cap and worker['guarded_wall_cap'] == wall_cap and
            worker['physical_error_read'] is False and worker['response_is_error_bound'] is False, 'Raw worker bounds and construction scope')
    for boundary in costs.values():
        require(boundary['raw_CPU_seconds'] is None or boundary['charged_CPU_seconds'] is None or
                boundary['charged_CPU_seconds'] >= boundary['raw_CPU_seconds'], 'Charge below raw CPU')
    ceiling = config['parent']['CPU_seconds'] if arm == 'parent' else config['cheap']['CPU_ceiling'] if arm in k.CHEAP else record['cpu_budget']
    wall_ceiling = config['parent']['wall_seconds'] if arm == 'parent' else config['cheap']['wall_ceiling'] if arm in k.CHEAP else 5*record['cpu_budget']
    complete = (receipt['exit_code'] == 0 and receipt.get('termination') is None and receipt.get('guard_error') is None and
        not receipt.get('cleanup_errors') and own['raw_CPU_seconds'] is not None and
        costs['total']['charged_CPU_seconds'] is not None and costs['total']['raw_CPU_seconds'] <= ceiling and
        costs['total']['charged_CPU_seconds'] <= ceiling and costs['total']['wall_seconds'] <= wall_ceiling)
    require(record['eligible'] is False or complete, 'Eligible record requires complete measured costs')
    return dict(complete=complete, own=own, plan=plan, receipt=receipt, worker=worker, flags=flags)


def parent_starts(parent, case, rep, size):
    saved = parent.get('restarts')
    require(isinstance(saved, dict) and saved.get('interior_perturbation_counts') == [0, 1, 2, 4], 'Measured four-start parent preparation')
    material = (case+':'+str(rep)+':cost-matched-v1').encode()
    seed = int.from_bytes(hashlib.sha256(material).digest()[:8], 'big'); rng = random.Random(seed)
    orders = [list(parent['order'])]
    require(size >= 5, 'Four-start interior perturbation domain')
    for count in (1, 2, 4):
        for _ in range(100):
            trial = list(parent['order'])
            for _ in range(count):
                source, target = rng.sample(range(1, size-1), 2)
                trial.insert(target, trial.pop(source))
            if trial not in orders:
                orders.append(trial); break
        else:
            raise ValueError('Four-start preparation refused')
    require(saved['seed'] == seed and saved['orders'] == orders, 'Frozen exact deterministic restart plan')
    return orders


def execution_admitted(record, parent):
    result = record.get('result'); arm = record['arm']
    if not isinstance(result, dict) or result.get('score') is None:
        return False
    if record.get('refusal') is not None or result.get('refusal') is not None:
        return False
    if arm in ('parent', 'full', 'window'):
        policy = k.load_module('_molecular_retained_admission', k.ROOT/'numerical_audit/cost_matched_continuation/admission.py')
        return policy.execution_admission(arm, result, None if parent is None else parent['order'])['admitted']
    if arm in k.CHEAP:
        return result.get('status') == 'complete'
    stats = result.get('stats', {}); raw = result.get('raw_julia_candidate')
    witness = raw.get('witness', {}) if isinstance(raw, dict) else {}
    return (result.get('eligible') is True and result.get('status') in
        ('complete', 'checked_update', 'slice_complete', 'cpu_limit', 'wall_limit',
         'evaluation_limit', 'proposal_limit', 'action_limit') and
        stats.get('four_chain_coverage') is True and stats.get('native_parent_revalidated') is True and
        stats.get('native_verifications') == 1 and witness.get('start_scores_completed') == 4 and
        isinstance(witness.get('chain_proposals'), list) and len(witness['chain_proposals']) == 4 and
        all(type(v) is int and v > 0 for v in witness['chain_proposals']) and isinstance(raw.get('native_score'), dict))


def raw_claims(record, size, worker=None):
    found = []
    def walk(value):
        if isinstance(value, dict):
            for key, child in value.items():
                if key in ('raw_julia_candidate', 'raw_proposed_best') and child is not None:
                    found.extend(child if isinstance(child, list) else [child])
                else:
                    walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)
    walk(record.get('result'))
    if worker is not None:
        for claim in worker.get('observed_raw_witnesses', []):
            found.append(claim)
    claims = {}
    for raw in found:
        require(isinstance(raw, dict), 'Persisted raw proposed-best object')
        order = k.validate_order(raw['order'], size)
        score = k.finite_score(raw['score'])
        key = (tuple(order), score['value'])
        native = raw.get('native_score')
        if native is not None:
            native = k.finite_score(native)
        else:
            require(raw.get('native_refusal') is not None or record['eligible'] is False,
                    'Missing raw native score must remain explicitly refused/ineligible')
        item = dict(order=order, score=score, native_score=native,
                    native_refusal=raw.get('native_refusal', record.get('refusal') or 'Native score not persisted'))
        require(key not in claims or claims[key] == item, 'Conflicting duplicate raw claim roles')
        claims[key] = item
    stats = (record.get('result') or {}).get('stats', {})
    if record['arm'] == 'window' and record['eligible']:
        require(stats.get('raw_candidate_inventory_complete') is True, 'Window raw claim observation completed')
    if record['arm'] == 'window' and stats.get('child_checkpoint_status') is not None:
        require(bool(claims), 'Window checkpoint exists but raw proposed-best claim was not retained')
    if record['arm'] == 'annealing' and stats.get('native_verification_attempts', 0):
        require(bool(claims), 'Annealing proposal exists but raw claim missing')
    return list(claims.values())


def case_inventory(a, freeze, info, config, inputs, case):
    require(case in info['cases'], 'Declared case only')
    if info['mode'] == 'preflight':
        spec = (14, 4)
    else:
        specs = {v['case']: v for v in freeze['registry']['cases']}
        spec = (specs[case]['n'], specs[case]['electrons'])
        require(spec == (12, 8), 'Frozen main CAS(8e,6o) case')
    model, units = authenticate_model(a, case, inputs['successes'][case], spec)
    if info['mode'] == 'main':
        require(all(model[key] == value for key, value in specs[case].items()), 'Exact frozen geometry/active-space specification')
    size = len(units); records = {}; metadata = []; all_claims = []
    expected_paths = set()
    for rep in range(1, info['repetitions']+1):
        for arm, budget in IDENTITIES:
            relative = f'constructors/{case}.rep{rep}.{arm}.B{budget}.json'; expected_paths.add(relative)
            record = a.json(relative); sha = a.hashes[relative]
            for key, value in dict(case=case, arm=arm, repetition=rep, cpu_budget=budget, n=spec[0], electrons=spec[1],
                model_json_sha256=inputs['successes'][case]['model.json'], model_toml_sha256=inputs['successes'][case]['model.toml'],
                protocol_sha256=info['protocol_sha256'], source_freeze_sha256=info['source_freeze_sha256'], score_kind=k.SCORE_KIND).items():
                require(type(record[key]) is type(value) and record[key] == value, 'Canonical constructor identity: '+key)
            require(record['schema'] == 'molecular-constructor-v1' and type(record['eligible']) is bool and
                    record['physical_error_read'] is False and record['response_is_error_bound'] is False,
                    'Constructor schema, Boolean eligibility and scope')
            k.validate_order(record['order'], size)
            score = k.finite_score({key: record[field] for key, field in
                (('value', 'score'), ('cubic_response', 'cubic_response'), ('steps', 'steps'))})
            lengths = [len(units[i]) for i in record['order']]; total = sum(lengths)
            steps = 1 if size == 1 else (30*total-lengths[0])//(2*total-lengths[0]-lengths[-1])
            require(score['steps'] == steps, 'Native score original-word affordability')
            result = record['result']
            parent = records.get((rep, 'parent', 0)) if arm in k.ARMS else None
            fallback = result is None or result.get('score') is None
            require(type(record['controller_fallback']) is bool and record['controller_fallback'] == fallback,
                    'Explicit scored fallback identity')
            if fallback:
                require(parent is not None and record['eligible'] is False and record['order'] == parent['order'] and
                        score == {key: parent[field] for key, field in (('value', 'score'), ('cubic_response', 'cubic_response'), ('steps', 'steps'))},
                        'Only an explicitly ineligible checked-parent fallback can remain scored')
            else:
                require(result['order'] == record['order'] and result['score'] == score,
                        'Top-level selected/native backend score and order join')
            if arm in k.ARMS:
                require(parent is not None and record['parent_sha256'] == parent['artifact_sha256'], 'Common measured parent digest')
            else:
                require(record.get('parent_sha256') is None, 'Parent/cheap JSON has no external parent')
            if arm == 'parent':
                parent_starts(record, case, rep, size)
            if arm == 'canonical':
                require(record['order'] == list(range(size)), 'Original canonical order')
            if arm == 'magnitude':
                magnitude = sorted(range(size), key=lambda i: (-sum(abs(v) for v in units[i].values()), i))
                require(record['order'] == magnitude, 'Exact native magnitude original-label order')
            evidence = cost_check(a, record, parent, config)
            if arm in k.ARMS and not fallback:
                require(evidence['worker']['starts'] == (parent['restarts']['orders'] if arm in ('full', 'annealing') else [parent['order']]),
                        'Identical supplied parent/four starts')
            if arm == 'parent':
                require(evidence['worker']['restarts'] == record['restarts'], 'Measured worker parent preparation')
            admitted = execution_admitted(record, parent)
            require(record['eligible'] is False or (admitted and (parent is None or parent['eligible'])),
                    'Claimed eligible constructor has incomplete method/parent execution')
            record = dict(record, artifact_path=relative, artifact_sha256=sha,
                          analysis_eligible=record['eligible'] and admitted and evidence['complete'])
            records[(rep, arm, budget)] = record
            claims = raw_claims(record, size, evidence['worker'])
            all_claims.extend((record, claim) for claim in claims)
            metadata.append(dict(case=case, repetition=rep, arm=arm, cpu_budget=budget, path=relative, sha256=sha,
                eligible=record['eligible'], analysis_eligible=record['analysis_eligible'], costs=record['costs'],
                status=None if result is None else result.get('status'), raw_claims=len(claims), refusal=record.get('refusal') or (result or {}).get('refusal')))
    observed_paths = {str(p.relative_to(a.root)) for p in (a.root/'constructors').glob(case+'.rep*.json')}
    require(observed_paths == expected_paths, 'Missing/extra canonical constructor identities')
    return dict(case=case, model=model, units=units, records=records, constructors=metadata, raw_claims=all_claims)


def manifests(inventory, info):
    case = inventory['case']; model = inventory['model']
    header = dict(version=1, order_index_base=1, case=case, n=model['n'], electrons=model['electrons'], partition='conserving',
        input_sha256=next(iter(inventory['records'].values()))['model_toml_sha256'],
        source_json_sha256=next(iter(inventory['records'].values()))['model_json_sha256'],
        protocol_sha256=info['protocol_sha256'], source_freeze_sha256=info['source_freeze_sha256'], mode=info['mode'],
        repetitions=info['repetitions'], expected_logical_rows=36*info['repetitions'], training_time=1.0,
        training_budget_multiplier=30, times=[.65, 1.35], rotation_cap_multipliers=[8, 24],
        protocol_path=info['protocol_path'], source_freeze_path=info['source_freeze_path'])
    artifacts = []; candidates = []; by_identity = {}
    for (rep, arm, budget), record in inventory['records'].items():
        oid = f'rep{rep}.{arm}.B{budget}'; by_identity[(rep, arm, budget)] = oid
        artifacts.append(dict(order_id=oid, path=record['artifact_path'], sha256=record['artifact_sha256']))
        for duration in (.65, 1.35):
            for cap in (8, 24):
                row = dict(name=f'{oid}.T{duration}.R{cap}', order_id=oid, arm=arm, repetition=rep, cpu_budget=budget,
                    partition='conserving', formula=2, time=duration, budget_multiplier=cap,
                    order=[v+1 for v in record['order']], constructor_sha256=record['artifact_sha256'],
                    claimed_score=record['score'], score_kind=k.SCORE_KIND, eligible=record['eligible'])
                if arm not in k.CHEAP:
                    row['parent_sha256'] = record['artifact_sha256'] if arm == 'parent' else record['parent_sha256']
                if not record['eligible']:
                    row['constructor_failure'] = json.dumps(record.get('refusal') or (record['result'] or {}).get('refusal') or
                        dict(status=(record['result'] or {}).get('status'), reason='Constructor method/cost eligibility refused'), sort_keys=True)
                candidates.append(row)
    witnesses = []
    for record, claim in inventory['raw_claims']:
        oid = by_identity[(record['repetition'], record['arm'], record['cpu_budget'])]
        index = sum(row['order_id'] == oid for row in witnesses)+1
        row = dict(name=f'{oid}.raw{index}', order_id=oid, arm=record['arm'], repetition=record['repetition'],
            cpu_budget=record['cpu_budget'], constructor_sha256=record['artifact_sha256'], partition='conserving',
            order=[v+1 for v in claim['order']], claimed_score=claim['score']['value'], score_kind=k.SCORE_KIND,
            role='raw_proposed_best', selected=claim['order'] == record['order'])
        if record['arm'] in k.ARMS:
            row['parent_sha256'] = record['parent_sha256']
        if claim['native_score'] is None:
            row['native_refusal'] = json.dumps(claim['native_refusal'], sort_keys=True)
        else:
            row['native_score'] = claim['native_score']['value']
        witnesses.append(row)
    return dict(header, constructor_artifacts=artifacts, candidates=candidates), dict(header,
        constructor_artifacts=artifacts, expected_witnesses=len(witnesses), witness_inventory_complete=True, witnesses=witnesses)


def validate_case(a, freeze, info, config, inputs, case, *, emitted=True):
    inventory = case_inventory(a, freeze, info, config, inputs, case)
    panel, witnesses = manifests(inventory, info)
    if emitted:
        actual_panel = a.toml(f'panels/{case}.toml')
        actual_witnesses = a.toml(f'witnesses/{case}.toml')
        require(actual_panel == panel, 'Emitted panel differs from independent raw constructor/model/cost derivation')
        require(actual_witnesses == witnesses, 'Emitted raw-witness inventory drops/adds/changes a persisted claim')
    return inventory, panel, witnesses


def _public(inventory):
    return dict(case=inventory['case'], constructors=inventory['constructors'],
        constructor_count=len(inventory['records']), expected_logical_rows=4*len(inventory['records']),
        raw_claims=len(inventory['raw_claims']), all_constructor_eligible=all(v['analysis_eligible'] for v in inventory['records'].values()))



def write_manifest(path, manifest, table_names):
    """Bounded scalar/array TOML output without importing a numerical backend."""
    def scalar(value):
        if isinstance(value, str):
            return json.dumps(value, ensure_ascii=True)
        if type(value) is bool:
            return 'true' if value else 'false'
        if type(value) in (int, float):
            number(value)
            return repr(value)
        require(isinstance(value, list), 'Scientific TOML scalar/array schema')
        return '['+', '.join(scalar(item) for item in value)+']'
    lines = [name+' = '+scalar(value) for name, value in manifest.items()
             if name not in table_names or not value]
    for name in table_names:
        for row in manifest[name]:
            lines.append('\n[['+name+']]')
            lines.extend(key+' = '+scalar(value) for key, value in row.items())
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= 8*1024**2, 'Scientific TOML byte cap')
    path = Path(path); temporary = path.with_name(path.name+'.partial')
    require(not path.exists() and not temporary.exists(), 'Fresh manifest and partial required')
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('xb') as stream:
        stream.write(raw)
    os.link(temporary, path); temporary.unlink()


def produce(data, *, preflight=False):
    a, freeze, info, config, inputs = context(data, producer=True, preflight=preflight)
    report = dict(schema='molecular-panels-v1', **info, complete=False, cases_metadata=[], failures=[],
                  original_frozen_sources=freeze['sources'], producer_sources=k.source_map(), physical_error_read=False)
    prepared = []
    for case in info['cases']:
        try:
            inventory, panel, witnesses = validate_case(a, freeze, info, config, inputs, case, emitted=False)
            prepared.append((case, panel, witnesses)); report['cases_metadata'].append(_public(inventory))
        except Exception as error:
            report['failures'].append(dict(case=case, type=type(error).__name__, message=str(error)))
    if not report['failures']:
        for case, panel, witnesses in prepared:
            write_manifest(a.root/'panels'/f'{case}.toml', panel, ('constructor_artifacts', 'candidates'))
            write_manifest(a.root/'witnesses'/f'{case}.toml', witnesses, ('constructor_artifacts', 'witnesses'))
        report['complete'] = True
    report['artifact_hashes'] = a.hashes
    k.write_new(a.root/'PANELS.json', report)
    return report


def validate(data, *, preflight=False):
    a, freeze, info, config, inputs = context(data, preflight=preflight)
    inventories = {}
    for case in info['cases']:
        inventories[case] = validate_case(a, freeze, info, config, inputs, case)[0]
    expected = {f'constructors/{case}.rep{rep}.{arm}.B{budget}.json' for case in info['cases']
                for rep in range(1, info['repetitions']+1) for arm, budget in IDENTITIES}
    require({str(p.relative_to(a.root)) for p in (a.root/'constructors').glob('*.json')} == expected,
            'Global constructor inventory includes missing/extra case identities')
    return dict(schema='molecular-panel-consumer-v1', complete=True, **info,
        cases_metadata=[_public(value) for value in inventories.values()],
        original_frozen_sources=freeze['sources'], current_consumer_sources=k.source_map(), artifact_hashes=a.hashes,
        raw_inventory_independently_confirmed=True, _inventories=inventories, _artifacts=a, _freeze=freeze)


def public_validation(report):
    return {key: value for key, value in report.items() if not key.startswith('_')}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('produce', 'validate')); parser.add_argument('data', type=Path)
    parser.add_argument('--preflight', action='store_true'); parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    report = produce(args.data, preflight=args.preflight) if args.action == 'produce' else public_validation(validate(args.data, preflight=args.preflight))
    if args.output:
        k.write_new(args.output, report)
    return 0 if report['complete'] else 2


if __name__ == '__main__':
    sys.exit(main())
