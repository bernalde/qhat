"""Metadata-only independent readback; never imports the scorer or runner."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import re
import tomllib

ROOT = Path(__file__).resolve().parent
PLAN_SHA = '3853b2c8adfe3383fff493adf35b15b13ca068290c7f1628efd34bee9af2d581'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for raw in iter(lambda: stream.read(1024**2), b''):
            value.update(raw)
    return value.hexdigest()


def ref(path):
    return dict(path=str(path.relative_to(ROOT)), bytes=path.stat().st_size, sha256=sha(path))


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def compare(left, right):
    difference = abs(left-right)
    tolerance = 1e-24+1e-8*max(left, right)
    require(all(math.isfinite(x) and x >= 0 for x in (left, right)), 'finite score values')
    return dict(verified=difference <= tolerance and max(left, right) > 1e-24,
                within_tolerance=difference <= tolerance, numerical_floor=max(left, right) <= 1e-24,
                absolute_difference=difference, tolerance=tolerance)


def main():
    plan = read(ROOT/'PLAN.a1.json')
    result = read(ROOT/'work/RESULT.json')
    receipt_path = ROOT/'logs/mg-hf-native-supplement-a1.resources.json'
    receipt = read(receipt_path)
    log_path = ROOT/'logs/mg-hf-native-supplement-a1.log'
    log = log_path.read_text()
    require(sha(ROOT/'PLAN.a1.json') == PLAN_SHA == result['plan_sha256'], 'immutable executed plan')
    require(sha(ROOT/'native_validate.py') == plan['runner_sha256'] == result['runner_sha256'], 'executed runner')
    require(sha(ROOT/'tools/resource_guard.py') == plan['guard_sha256'], 'original guard')
    for name, identity in plan['inputs'].items():
        path = ROOT/name
        require(path.stat().st_size == identity['bytes'] and sha(path) == identity['sha256'], 'input byte readback')
    for name, identity in plan['sources'].items():
        path = ROOT/'source_snapshot'/name
        require(path.stat().st_size == identity['bytes'] and sha(path) == identity['sha256'], 'source byte readback')
    require(result['loaded_sources'] == {name: x['sha256'] for name, x in plan['sources'].items()}, 'six loaded source identities')
    require(receipt['command'] == plan['command'] and receipt['exit_code'] == receipt['supervisor_exit_code'] == 0
            and receipt['termination'] is None and receipt['guard_error'] is None and not receipt['cleanup_errors'], 'closed successful actual guard')
    require(receipt['budget'] == dict(cpu_seconds=20., wall_seconds=120, disk_mib=64), 'actual finite guard limits')
    require(all(receipt['settings'][key] == value for key, value in dict(CPUQuota='25%', CPUWeight='5', IOWeight='10',
            Nice='19', MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0').items()), 'actual OS profile')
    require(receipt['environment']['PYTHONDONTWRITEBYTECODE'] == '1' and '-B' in receipt['command'], 'bytecode disabled')
    require(result['complete'] is True and result['status'] == 'verified_later_native_only'
            and result['score_attempts'] == result['completed_native_scores'] == 1 and not result['failures'], 'one complete later score')
    constructor = read(ROOT/'inputs/CONSTRUCTOR.json')
    worker = read(ROOT/'inputs/WORKER.json')
    claim = constructor['result']['stats']['raw_julia_candidate']
    objective = tomllib.loads((ROOT/'inputs/OBJECTIVE.toml').read_text())
    raw = next(x for x in objective['raw_witnesses'] if x['name'] == 'rep2.window.B64.raw1')
    archived = objective['checks'][raw['check_index']-1]
    require(claim['native_refusal'] == plan['target']['native_refusal'] == result['original_native_refusal']
            and 'native_score' not in claim and worker['observed_raw_witnesses'] == [claim]
            and raw['selected'] is False and raw['status'] == 'native_refusal_censored'
            and 'native_score' not in raw, 'original missing native check preserved')
    order = claim['order']
    order_sha = hashlib.sha256(json.dumps(order, ensure_ascii=True, separators=(',', ':')).encode('ascii')).hexdigest()
    require(order_sha == plan['target']['order_sha256'] == result['target']['order_sha256'], 'fixed complete order identity')
    require(archived['order'] == [x+1 for x in order], 'archived one-based order identity')
    model = read(ROOT/'inputs/model.json')
    model_toml = tomllib.loads((ROOT/'inputs/model.toml').read_text())
    require(model_toml['source_json_sha256'] == sha(ROOT/'inputs/model.json') and model_toml['n'] == model['n'] == 12
            and model_toml['electrons'] == model['electrons'] == 8
            and int(model_toml['denominator']) == int(model['denominator']), 'JSON/TOML model identity')
    for left, right in zip(model['blocks'], model_toml['units'], strict=True):
        original = {tuple(tuple(x) for x in word): int(value) for word, value in left}
        toml_words = {}
        for spelling, numerator in zip(right['words'], right['numerators'], strict=True):
            word = tuple((int(x[:-1]), x[-1]) for x in spelling.split(',')) if spelling else ()
            toml_words[word] = int(numerator)
        require(original == toml_words, 'exact unit coefficients and original labels')
    lengths = [len(unit['words']) for unit in model_toml['units']]
    words = sum(lengths)
    steps = (18900-lengths[order[0]])//(2*words-lengths[order[0]]-lengths[order[-1]])
    native = result['native_score']
    require(words == 630 and len(lengths) == 161 and steps == native['steps'] == archived['steps'] == 15
            and native['value'] == native['cubic_response']/steps**4, 'unchanged boundary fusion and normalized score')
    independent = dict(value=archived['replayed_score'], cubic_response=archived['details']['cubic_response'], steps=archived['steps'])
    checked = {}
    for name, other in (('original_julia', claim['score']), ('archived_independent', independent)):
        checked[name] = {field: compare(native[field], other[field]) for field in ('value', 'cubic_response')}
        checked[name]['steps_match'] = native['steps'] == other['steps']
    require(checked == result['archived_policy_comparisons']
            and all(x[field]['verified'] for x in checked.values() for field in ('value', 'cubic_response')), 'independent tolerance readback')
    difference = abs(native['value']-claim['score']['value'])
    tolerance = max(1e-28, abs(native['value'])*1e-8)
    require(result['original_window_comparison'] == dict(verified=difference <= tolerance,
            absolute_difference=difference, tolerance=tolerance) and difference <= tolerance, 'original window tolerance readback')
    original = Path(plan['original_data'])
    inventory = read(original/'EVIDENCE_INVENTORY.json')
    require(sha(original/'EVIDENCE_INVENTORY.json') == plan['original_inventory_sha256'], 'original seal identity')
    for name, identity in inventory['files'].items():
        path = original/name
        require(path.stat().st_size == identity['bytes'] and sha(path) == identity['sha256'], 'all original sealed bytes unchanged')
    actual = {str(path.relative_to(original)) for path in original.rglob('*') if path.is_file()}
    require(actual == set(inventory['files'])|{'EVIDENCE_INVENTORY.json'}, 'no additional original DATA files')
    require(not list((ROOT/'source_snapshot').rglob('__pycache__')), 'no supplementary source bytecode')
    require(result['original_before'] == result['original_after'] and result['original_after']['files'] == 1111
            and result['original_after']['all_original_bytes_match'], 'original in-worker integrity checks')
    measurements = receipt['process_measurements']
    raw_cpu = measurements['user_cpu_seconds']+measurements['system_cpu_seconds']
    sampled = [int(m.group(1))/1e9 for row in receipt['samples']
               if (m := re.search(r'^CPUUsageNSec=(\d+)$', row['service'], re.M))]
    main_charge = max([raw_cpu+.011]+sampled)
    terminal = re.search(r'^CPU time consumed:\s*([0-9.]+)(ms|s)$', log, re.M)
    require(terminal is not None, 'terminal cgroup CPU retained')
    terminal_cpu = float(terminal.group(1))*(.001 if terminal.group(2) == 'ms' else 1.)
    quantum = .001 if terminal.group(2) == 'ms' else 10**(-len(terminal.group(1).partition('.')[2]))
    conservative = math.ceil(max(main_charge, terminal_cpu+quantum)*1000)/1000
    require(conservative <= 20 and receipt['elapsed'] <= 120 and receipt['final_accounting']['output_bytes'] <= 64*1024**2, 'actual costs inside pool')
    costs = dict(raw_GNU_CPU_seconds=raw_cpu, existing_main_methodology_charge_CPU_seconds=main_charge,
        terminal_cgroup_CPU_seconds=terminal_cpu, charged_CPU_seconds=conservative,
        charge_policy='One conservative charge: max(GNU user+sys+.011, sampled cgroup CPU, terminal cgroup display plus one display quantum), rounded up to milliseconds; alternative bases are not extra spending',
        wall_seconds=receipt['elapsed'], guarded_output_bytes=receipt['final_accounting']['output_bytes'],
        GNU_peak_RSS_KiB=measurements['max_rss_kib'], cgroup_sampled_memory_peak_bytes=max(
            int(m.group(1)) for row in receipt['samples'] if (m := re.search(r'^MemoryPeak=(\d+)$', row['service'], re.M))))
    proof = dict(schema='molecular-hf-later-native-readback-v1', metadata_only=True, numerical_calls=0,
        created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), readback_source=ref(Path(__file__)),
        plan=ref(ROOT/'PLAN.a1.json'), runner=ref(ROOT/'native_validate.py'), result=ref(ROOT/'work/RESULT.json'),
        receipt=ref(receipt_path), log=ref(log_path), complete=True, checks_passed=True, failures=[],
        native_score=native, independently_read_back_comparisons=checked, original_window_tolerance_verified=True,
        original_inventory_sha256=sha(original/'EVIDENCE_INVENTORY.json'), original_files_verified=len(inventory['files']),
        original_FILE_set_unchanged=True, all_original_bytes_unchanged=True, original_native_refusal_preserved=True,
        original_scientific_gate_complete=False, original_primary_decisions_both='incomplete', costs=costs,
        metadata_closeout_scope='Light coordinator readback; no scorer import, new numerical work, test or rerun; not a formal Review/Verify round')
    write(ROOT/'READBACK.json', proof)
    final = dict(schema='molecular-hf-later-native-final-v1', completed=True, status='verified_later_native_only',
        role='later_supplementary_native_validation', plan=ref(ROOT/'PLAN.a1.json'), result=ref(ROOT/'work/RESULT.json'),
        receipt=ref(receipt_path), log=ref(log_path), independent_metadata_readback=ref(ROOT/'READBACK.json'),
        native_score=native, costs=costs, original_scientific_gate_complete=False,
        original_source_analysis_complete=False, original_primary_decisions_both='incomplete',
        original_censor_repaired=False, original_constructor_promoted=False, original_32_64_costs_changed=False,
        original_DATA_unchanged=True, no_automatic_retry=True)
    write(ROOT/'FINAL.json', final)
    print(json.dumps(dict(readback_sha256=sha(ROOT/'READBACK.json'), final_sha256=sha(ROOT/'FINAL.json'),
                         status=final['status'], costs=costs), indent=2))


if __name__ == '__main__':
    main()
