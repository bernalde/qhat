"""One later native check of one frozen, historically censored HF order."""
import time
BOOT_CPU, BOOT_WALL = time.process_time(), time.monotonic()
import argparse
from dataclasses import asdict
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tomllib

sys.dont_write_bytecode = True


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for raw in iter(lambda: stream.read(1024**2), b''):
            digest.update(raw)
    return digest.hexdigest()


def read(path):
    require(Path(path).stat().st_size <= 4*1024**2, 'bounded JSON')
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def bad(_):
        raise ValueError('nonfinite JSON')
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=bad)


def relative(root, name):
    require(isinstance(name, str) and not Path(name).is_absolute(), 'relative artifact path')
    path = root/name
    require(path.resolve(strict=True).is_relative_to(root) and path.is_file(), 'contained artifact')
    return path


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('ascii')
    require(len(raw) <= 256*1024, 'bounded result')
    with path.open('xb') as stream:
        stream.write(raw)


def parity(actual, claimed):
    require(all(type(v) in (int, float) and math.isfinite(v) and v >= 0
                for v in (actual, claimed)), 'finite nonnegative scores')
    difference = abs(actual-claimed)
    tolerance = 1e-24+1e-8*max(actual, claimed)
    floor = max(actual, claimed) <= 1e-24
    return dict(verified=difference <= tolerance and not floor,
                within_tolerance=difference <= tolerance, numerical_floor=floor,
                absolute_difference=difference, tolerance=tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan_path = args.plan.resolve(strict=True)
    root = plan_path.parent
    plan = read(plan_path)
    require(plan['schema'] == 'molecular-hf-later-native-plan-v1', 'plan schema')
    output = args.output.resolve()
    require(str(output) == plan['output'] and output.parent == root and not output.exists(), 'fresh owned output')
    output.mkdir()
    result = dict(schema='molecular-hf-later-native-result-v1',
        role='later_supplementary_native_validation', plan_sha256=sha(plan_path),
        original_scientific_gate_complete=False, original_source_analysis_complete=False,
        original_primary_decisions=dict(window_minus_full='incomplete', annealing_minus_full='incomplete'),
        original_native_score_absent=True, original_native_refusal=plan['target']['native_refusal'],
        original_constructor_selected=False, original_censor_repaired=False,
        optimizer_called=False, chemistry_called=False, physical_propagation_called=False,
        score_attempts=0, completed_native_scores=0, complete=False, status='refused',
        actual_started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        target=plan['target'], input_hashes=plan['inputs'], frozen_sources=plan['sources'],
        runner_sha256=sha(__file__), phases={}, failures=[])
    native = token = None
    def check():
        require(time.process_time()-BOOT_CPU < 18., 'supplement CPU allowance exhausted')
        require(time.monotonic()-BOOT_WALL < 110., 'supplement wall allowance exhausted')
    def integrity():
        original = Path(plan['original_data']).resolve(strict=True)
        inventory = original/'EVIDENCE_INVENTORY.json'
        require(sha(inventory) == plan['original_inventory_sha256'], 'original inventory identity')
        data = read(inventory)
        require(len(data['files']) == 1111, 'original sealed roster')
        total = 0
        for name, ref in data['files'].items():
            check(); path = relative(original, name)
            require(path.stat().st_size == ref['bytes'] and sha(path) == ref['sha256'], 'original sealed bytes changed')
            total += ref['bytes']
        return dict(files=len(data['files']), file_bytes=total,
                    inventory_sha256=sha(inventory), all_original_bytes_match=True)
    try:
        require(sys.version == plan['runtime']['version'] and str(Path(sys.executable).resolve()) == plan['runtime']['resolved_executable']
                and sha(Path(sys.executable).resolve()) == plan['runtime']['executable_sha256'], 'pinned Python runtime')
        require(result['runner_sha256'] == plan['runner_sha256'], 'pinned runner')
        require(os.environ.get('PYTHONDONTWRITEBYTECODE') == '1' and sys.dont_write_bytecode, 'no bytecode writes')
        for name, ref in plan['inputs'].items():
            check(); path = relative(root, name)
            require(path.stat().st_size == ref['bytes'] and sha(path) == ref['sha256'], 'frozen supplementary input')
        source_root = (root/'source_snapshot').resolve(strict=True)
        freeze = read(root/'inputs/INPUT_PROTOCOL_FROZEN.json')
        for name, ref in plan['sources'].items():
            check(); path = relative(source_root, name)
            require(freeze['sources'][name] == ref['sha256'] and path.stat().st_size == ref['bytes']
                    and sha(path) == ref['sha256'], 'frozen original source')
        result['original_before'] = integrity()
        constructor = read(root/'inputs/CONSTRUCTOR.json')
        worker = read(root/'inputs/WORKER.json')
        model = read(root/'inputs/model.json')
        claim = constructor['result']['stats']['raw_julia_candidate']
        raw_order = claim['order']
        target = plan['target']
        order_sha = hashlib.sha256(json.dumps(raw_order, separators=(',', ':'), ensure_ascii=True,
                                              allow_nan=False).encode('ascii')).hexdigest()
        require(order_sha == target['order_sha256'] and len(raw_order) == 161 and
                all(type(i) is int for i in raw_order) and sorted(raw_order) == list(range(161)), 'fixed complete target order')
        require(worker['observed_raw_witnesses'] == [claim] and 'native_score' not in claim and
                claim['native_refusal'] == target['native_refusal'], 'original refused native attempt')
        require(claim['score'] == target['julia_claim'] and constructor['case'] == model['case'] == target['case']
                and model['n'] == 12 and model['electrons'] == 8, 'fixed model and claim')
        request = tomllib.loads((root/'inputs/REQUEST.toml').read_text())
        checkpoint = tomllib.loads((root/'inputs/CHECKPOINT.toml').read_text())
        objective = tomllib.loads((root/'inputs/OBJECTIVE.toml').read_text())
        require(request['time'] == '1/1' and request['budget'] == 18900 and request['reference'] == 255
                and request['n'] == 12 and request['max_entries'] == 50000 and request['max_cached_actions'] == 100000,
                'original native evaluation parameters')
        require(checkpoint['order'] == raw_order and checkpoint['score'] == target['julia_claim']['value']
                and checkpoint['cubic_response'] == target['julia_claim']['cubic_response'] and checkpoint['steps'] == 15,
                'original Julia checkpoint join')
        witness = next(row for row in objective['raw_witnesses'] if row['name'] == target['name'])
        archived = objective['checks'][witness['check_index']-1]
        require(witness['status'] == 'native_refusal_censored' and witness['selected'] is False
                and 'native_score' not in witness and archived['order'] == [i+1 for i in raw_order]
                and archived['status'] == 'replayed' and archived['steps'] == 15, 'original archived check join')
        expected_archived = dict(value=archived['replayed_score'], cubic_response=archived['details']['cubic_response'], steps=archived['steps'])
        require(expected_archived == target['archived_claim'], 'frozen archived comparison')
        phase = time.process_time()
        spec = importlib.util.spec_from_file_location('_hf_supplement_common',
            source_root/'numerical_audit/molecular_generalization/common.py')
        common = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = common
        spec.loader.exec_module(common)
        native = common.bootstrap()
        token = native._EXTERNAL_BUDGET_CHECK.set(check)
        blocks = common.blocks(model)
        require(len(blocks) == 161 and sum(map(len, blocks)) == 630, 'original units and words')
        loaded = {str(Path(module.__file__).resolve().relative_to(source_root)): sha(module.__file__)
                  for module in tuple(sys.modules.values()) if getattr(module, '__file__', None)
                  and Path(module.__file__).resolve().is_relative_to(source_root)}
        require(loaded == {name: ref['sha256'] for name, ref in plan['sources'].items()}, 'exact loaded original closure')
        result['loaded_sources'] = loaded
        result['phases']['imports_deserialization_CPU_seconds'] = time.process_time()-phase
        phase = time.process_time(); check()
        context = native.SparseResponse(blocks, 12, 255, time=1., max_entries=50000,
                                        max_cached_actions=100000, deadline=BOOT_CPU+18.)
        result['phases']['context_CPU_seconds'] = time.process_time()-phase
        phase = time.process_time(); check(); result['score_attempts'] += 1
        score = context.score(tuple(raw_order), 18900, BOOT_CPU+18.)
        result['completed_native_scores'] += 1
        result['phases']['native_score_CPU_seconds'] = time.process_time()-phase
        result['native_score'] = asdict(score)
        result['native_statistics'] = dict(context.stats)
        require(score.steps == 15, 'exact affordable step identity')
        comparisons = {}
        for role, other in (('original_julia', target['julia_claim']), ('archived_independent', expected_archived)):
            comparisons[role] = {field: parity(getattr(score, field), other[field]) for field in ('value', 'cubic_response')}
            comparisons[role]['steps_match'] = score.steps == other['steps']
        difference = abs(score.value-target['julia_claim']['value'])
        tolerance = max(1e-28, abs(score.value)*1e-8)
        result['original_window_comparison'] = dict(verified=difference <= tolerance,
            absolute_difference=difference, tolerance=tolerance)
        result['archived_policy_comparisons'] = comparisons
        require(all(row[field]['verified'] for row in comparisons.values() for field in ('value', 'cubic_response'))
                and all(row['steps_match'] for row in comparisons.values()) and difference <= tolerance,
                'supplement numerical comparison failed or censored')
        check(); result['original_after'] = integrity()
        require(result['original_after'] == result['original_before'], 'original graph unchanged')
        result.update(complete=True, status='verified_later_native_only')
    except Exception as error:
        result['failures'].append(dict(type=type(error).__name__, message=str(error)[:1024]))
    finally:
        if native is not None and token is not None:
            native._EXTERNAL_BUDGET_CHECK.reset(token)
        result['actual_finished_UTC'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        result['process_CPU_before_result_seconds'] = time.process_time()-BOOT_CPU
        result['wall_before_result_seconds'] = time.monotonic()-BOOT_WALL
        write_new(output/'RESULT.json', result)
        print(json.dumps(dict(status=result['status'], complete=result['complete'],
            score_attempts=result['score_attempts'], completed_native_scores=result['completed_native_scores'],
            native_score=result.get('native_score'), result_sha256=sha(output/'RESULT.json'))), flush=True)
    return 0 if result['complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
