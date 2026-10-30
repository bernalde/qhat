"""Owned-process Julia full-order annealing, with unchanged native final checks.

Call after the study bootstrap has rooted the qhat namespace in this checkout.
Settings accept the keys in DEFAULT_SETTINGS; calibrated=True, a 64-hex
calibration_sha256 and positive decreasing theta_start/theta_end are mandatory.
seed_material is normally ``case + ':' + str(repetition)``. All wire/result
orders use zero-based original labels. No validation-error input is accepted.

The caller owns the outer OS receipt and shared-parent attribution. This module
charges its complete Python call plus exactly one owned child, including child
startup/JIT, native revalidation, serialization, final verification and teardown.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time as _time
import tomllib
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_NAMES = (
    'analysis/sparse_response_ordering.py', 'analysis/strang_schedule.py',
    'analysis/julia_window/ResponseCore.jl', 'analysis/julia_window/Project.toml',
    'numerical_audit/molecular_generalization/src/FullOrderAnnealing.jl',
    'numerical_audit/molecular_generalization/anneal.jl',
    'numerical_audit/molecular_generalization/annealing.py')
DEFAULT_SETTINGS = dict(version=1, calibrated=False, calibration_sha256='',
    theta_start=0., theta_end=0., energy_floor=1e-30,
    max_evaluations=1_000_000, max_proposals=1_000_000, max_actions=100_000_000,
    cache_bytes=64*1024**2, native_reserve_factor=2., native_reserve_min=.05,
    cleanup_CPU_seconds=.05)
HEX = re.compile(r'[0-9a-f]{64}')


class Refusal(ValueError):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class AnnealingResult:
    order: tuple[int, ...] | None
    score: Any
    eligible: bool
    status: str
    refusal: dict | None
    raw_julia_candidate: dict | None
    native_promotion: str
    CPU_seconds: float
    stats: dict


def require(condition, message, code='invalid_input'):
    if not condition:
        raise Refusal(code, message)


def integer(value, name, low, high):
    require(type(value) is int and low <= value <= high, 'Integer domain: '+name)
    return value


def finite(value, name, *, positive=False):
    require(type(value) in (int, float), 'Numeric scalar required: '+name)
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise Refusal('invalid_input', 'Float range: '+name) from exc
    require(math.isfinite(result) and (not positive or result > 0), 'Finite scalar required: '+name)
    return result


def validate_settings(settings):
    require(isinstance(settings, Mapping) and len(settings) <= len(DEFAULT_SETTINGS)
            and set(settings) <= DEFAULT_SETTINGS.keys(), 'Annealing settings schema')
    out = dict(DEFAULT_SETTINGS, **settings)
    require(type(out['version']) is int and out['version'] == 1, 'Settings version')
    require(out['calibrated'] is True and isinstance(out['calibration_sha256'], str)
            and HEX.fullmatch(out['calibration_sha256']) is not None,
            'Frozen old-only calibration required', 'uncalibrated_settings')
    hi = finite(out['theta_start'], 'theta_start', positive=True)
    lo = finite(out['theta_end'], 'theta_end', positive=True)
    require(hi > lo, 'Temperature endpoints must decrease')
    require(type(out['energy_floor']) is float and out['energy_floor'] == 1e-30, 'Frozen energy floor')
    for key, low, high in (('max_evaluations', 4, 1_000_000), ('max_proposals', 1, 1_000_000),
                          ('max_actions', 1, 100_000_000), ('cache_bytes', 0, 64*1024**2)):
        integer(out[key], key, low, high)
    require(1 <= finite(out['native_reserve_factor'], 'native_reserve_factor') <= 10, 'Native reserve factor cap')
    for key in ('native_reserve_min', 'cleanup_CPU_seconds'):
        require(.01 <= finite(out[key], key) <= 10, 'Native/cleanup reserve cap')
    out['theta_start'], out['theta_end'] = hi, lo
    return out


def _native():
    native = importlib.import_module('qhat.analysis.sparse_response_ordering')
    schedule = importlib.import_module('qhat.analysis.strang_schedule')
    require(Path(native.__file__).resolve() == ROOT/'analysis/sparse_response_ordering.py'
            and Path(schedule.__file__).resolve() == ROOT/'analysis/strang_schedule.py'
            and native.plan_strang_schedule is schedule.plan_strang_schedule,
            'Native import is not from this exact checkout', 'source_identity_mismatch')
    return native


def source_hashes(check):
    out = {}
    for name in SOURCE_NAMES:
        check()
        path = ROOT/name
        require(path.is_file() and path.stat().st_size <= 1024**2,
                'Runtime source missing/byte cap: '+name, 'source_identity_mismatch')
        with path.open('rb') as stream:
            raw = stream.read(1024**2+1)
        require(len(raw) <= 1024**2, 'Growing source byte cap', 'source_identity_mismatch')
        out[name] = hashlib.sha256(raw).hexdigest()
    return out


def _toml(value):
    if type(value) is bool:
        return 'true' if value else 'false'
    if type(value) is int:
        return str(value)
    if type(value) is float and math.isfinite(value):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    if isinstance(value, (tuple, list)):
        return '['+','.join(_toml(item) for item in value)+']'
    raise Refusal('invalid_input', 'Unsupported TOML scalar')


def _fraction(value):
    q = Fraction(value)
    require(q.numerator.bit_length() <= 4096 and q.denominator.bit_length() <= 4096,
            'Coefficient exceeds rational bit cap')
    return str(q.numerator)+'/'+str(q.denominator)


def _seeds(material):
    require(isinstance(material, str) and material.isascii() and 0 < len(material) <= 256,
            'Bounded ASCII seed material required')
    return [str(int.from_bytes(hashlib.sha256((material+
        ':molecular-generalization-sa-v1:chain'+str(k)).encode('ascii')).digest()[:8], 'big'))
        for k in range(1, 5)]


def _score_dict(score):
    if isinstance(score, Mapping):
        require(set(score) == {'value', 'cubic_response', 'steps'}, 'Parent score fields')
        out = dict(score)
    else:
        out = {key: getattr(score, key) for key in ('value', 'cubic_response', 'steps')}
    for key in ('value', 'cubic_response'):
        out[key] = finite(out[key], key)
        require(out[key] >= 0, 'Negative score')
    integer(out['steps'], 'score steps', 1, (1 << 62)-1)
    require(abs(out['value']-out['cubic_response']/out['steps']**4)
            <= max(1e-30, abs(out['value'])*1e-12), 'Parent score normalization')
    require(not (out['cubic_response'] > 0 and out['value'] == 0), 'Parent score underflows')
    return out


def _close(actual, claimed):
    return abs(actual-claimed) <= max(1e-28, abs(actual)*1e-8)


def _checked_score(score):
    require(score is not None and math.isfinite(score.value) and math.isfinite(score.cubic_response)
            and min(score.value, score.cubic_response) >= 0 and type(score.steps) is int and score.steps > 0,
            'Finite complete native score required', 'native_score_refusal')
    require(not (score.cubic_response > 0 and score.value == 0), 'Native normalized score underflows',
            'native_score_refusal')
    return score


def _compiled_units(native, blocks, n, check):
    integer(n, 'register', 1, 14)
    require(isinstance(blocks, (tuple, list)) and 1 <= len(blocks) <= 1024, 'Original unit cap')
    require(all(isinstance(b, Mapping) and b for b in blocks), 'Original nonempty units required')
    require(sum(map(len, blocks)) <= 8192, 'Original word cap')
    compiled = native._compile(blocks, n)
    units, seen = [], set()
    for original, unit in zip(blocks, compiled):
        check()
        require(unit and len(original) == len(unit) and len({x for x, _, _, _ in unit}) == 1
                and all(not phase.imag for _, _, phase, _ in unit), 'Real equal-flip nonzero original units required')
        parities, signs, coefficients = [], [], []
        for flip, parity, phase, coefficient in unit:
            check()
            key = flip, parity
            require(key not in seen, 'Duplicate global original Pauli word')
            seen.add(key)
            parities.append(parity)
            signs.append(1 if (phase.real > 0) == (coefficient > 0) else -1)
            coefficients.append(_fraction(coefficient))
        units.append(dict(flip=unit[0][0], parities=parities, signs=signs, coefficients=coefficients))
    return compiled, units


def _coefficient_hash(units):
    raw = ''.join(str(u['flip'])+':'+ ';'.join(str(p)+','+str(s)+','+c for p, s, c in
        zip(u['parities'], u['signs'], u['coefficients']))+'\n' for u in units).encode('ascii')
    return hashlib.sha256(raw).hexdigest()


def _request(units, n, reference, parent, starts, score, budget, settings,
             seed_material, max_entries, max_cached_actions, sources, check):
    settings_sha = hashlib.sha256(json.dumps(settings, sort_keys=True, separators=(',', ':'),
        allow_nan=False).encode('ascii')).hexdigest()
    fields = dict(version=1, n=n, reference=reference, time='1/1', budget=budget,
        order=parent, starts=starts, word_lengths=[len(u['coefficients']) for u in units],
        parent_score=float(score.value), parent_cubic_response=float(score.cubic_response), parent_steps=score.steps,
        max_entries=max_entries, max_cached_actions=max_cached_actions, seed_material=seed_material,
        seeds=_seeds(seed_material), settings_sha256=settings_sha, coefficients_sha256=_coefficient_hash(units))
    lines = [key+'='+_toml(value) for key, value in fields.items()]
    for name, mapping in (('settings', settings), ('source_sha256', sources)):
        lines.append('['+name+']')
        lines.extend(_toml(key)+'='+_toml(value) for key, value in sorted(mapping.items()))
    for unit in units:
        check();lines.append('[[units]]')
        lines.extend(key+'='+_toml(value) for key, value in unit.items())
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= 8*1024**2, 'Request byte cap')
    return raw, fields


def _proc_cpu(pid):
    try:
        raw = Path('/proc/'+str(pid)+'/stat').read_text()
        fields = raw[raw.rfind(')')+2:].split()
        return (int(fields[11])+int(fields[12]))/os.sysconf('SC_CLK_TCK')
    except (OSError, ValueError, IndexError):
        return None


def _run_worker(executable, request, checkpoint, child_cap, wall_cap, check, stats):
    worker = HERE/'anneal.jl'
    argv = [executable, '--startup-file=no', '--history-file=no', '--threads=1', '--gcthreads=1',
        '--heap-size-hint=32M', '--project='+str(ROOT/'analysis/julia_window'), str(worker),
        str(request), str(checkpoint), repr(child_cap), repr(wall_cap)]
    stats['child_command'] = list(argv)
    env = os.environ.copy()
    env.update(JULIA_NUM_THREADS='1', JULIA_NUM_GC_THREADS='1', OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', JULIA_LOAD_PATH='@:@stdlib',
        JULIA_PKG_PRECOMPILE_AUTO='0', JULIA_IMAGE_THREADS='1')
    began = _time.monotonic()
    process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, env=env, start_new_session=True)
    reaped, termination = False, None
    try:
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                process.returncode = os.waitstatus_to_exitcode(status)
                reaped = True
                stats.update(child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,
                    child_peak_RSS_KiB=usage.ru_maxrss, child_returncode=process.returncode)
                break
            current = _proc_cpu(process.pid)
            if current is None:
                termination = 'owned_child_cpu_unavailable'
            elif current >= child_cap:
                termination = 'child_cpu_limit'
            elif _time.monotonic()-began >= wall_cap:
                termination = 'child_wall_limit'
            if termination is None:
                try:
                    check(current)
                except Refusal as exc:
                    if exc.code not in ('child_cpu_limit', 'child_wall_limit'):
                        raise
                    termination = exc.code
            if termination is not None:
                break
            _time.sleep(.005)
    finally:
        if not reaped:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            stats.update(child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,
                child_peak_RSS_KiB=usage.ru_maxrss, child_returncode=process.returncode)
        stats['child_wall_seconds'] = _time.monotonic()-began
        stats['child_reaped'] = True
    if termination:
        return termination
    if process.returncode:
        refusal = Path(str(checkpoint)+'.refusal')
        if refusal.is_file():
            with refusal.open('rb') as stream:
                reason = stream.read(129)
            if reason in {b'cpu_limit', b'wall_limit', b'action_limit', b'numerical_limit',
                          b'admission_or_support_limit', b'worker_failure'}:
                return 'child_'+reason.decode('ascii')
        return 'child_failure'
    return 'complete'


REPLY_FIELDS = {'version', 'request_sha256', 'settings_sha256', 'coefficients_sha256', 'source_sha256',
    'word_lengths', 'time', 'budget', 'julia_version', 'rng', 'score_scope', 'phases', 'status', 'order',
    'score', 'cubic_response', 'steps', 'start_scores_completed', 'evaluations', 'proposals',
    'accepted_moves', 'uphill_attempts', 'uphill_accepted', 'chain_proposals', 'chain_accepted',
    'chain_CPU_seconds', 'seeds', 'stop_reason', 'kernel_stats'}
REPLY_STATUSES = {'initializing', 'checked_update', 'slice_complete', 'complete', 'cpu_limit',
    'wall_limit', 'evaluation_limit', 'proposal_limit', 'action_limit', 'numerical_limit', 'support_or_admission_limit'}


def _bounded_reply(path):
    require(path.is_file() and path.stat().st_size <= 64*1024, 'Missing/oversized checkpoint', 'invalid_worker_reply')
    with path.open('rb') as stream:
        raw = stream.read(64*1024+1)
    require(len(raw) <= 64*1024 and raw.isascii() and len(raw.splitlines()) <= 512,
            'Checkpoint text cap', 'invalid_worker_reply')
    quote, escaped, nesting = 0, False, []
    for byte in raw:
        if quote:
            require(byte not in (10, 13), 'Multiline checkpoint string', 'invalid_worker_reply')
            if escaped:
                escaped = False
            elif quote == 34 and byte == 92:
                escaped = True
            elif byte == quote:
                quote = 0
        elif byte in (34, 39):
            quote = byte
        elif byte in (91, 123):
            nesting.append(93 if byte == 91 else 125)
            require(len(nesting) <= 3, 'Checkpoint nesting cap', 'invalid_worker_reply')
        elif byte in (93, 125):
            require(nesting and nesting.pop() == byte, 'Checkpoint brackets', 'invalid_worker_reply')
        elif byte == 35:
            raise Refusal('invalid_worker_reply', 'Checkpoint comments outside schema')
    require(not quote and not nesting, 'Unterminated checkpoint', 'invalid_worker_reply')
    try:
        return tomllib.loads(raw.decode('ascii'))
    except tomllib.TOMLDecodeError as exc:
        raise Refusal('invalid_worker_reply', 'Malformed checkpoint TOML') from exc


def _validated_reply(path, digest, fields, sources, settings, native, context, check):
    check();data = _bounded_reply(path);check()
    require(set(data) == REPLY_FIELDS and type(data['version']) is int and data['version'] == 1,
            'Checkpoint schema', 'invalid_worker_reply')
    require(data['request_sha256'] == digest and data['settings_sha256'] == fields['settings_sha256']
            and data['coefficients_sha256'] == fields['coefficients_sha256'] and data['source_sha256'] == sources,
            'Checkpoint request/source/settings/coefficients identity', 'invalid_worker_reply')
    order = native._permutation(data['order'], len(context.units))
    require(data['word_lengths'] == list(context.lengths) and all(type(i) is int for i in data['word_lengths'])
            and type(data['time']) is float and data['time'] == 1.
            and type(data['budget']) is int and data['budget'] == fields['budget']
            and data['seeds'] == fields['seeds'] and data['rng'] == 'Random.Xoshiro'
            and data['score_scope'] == 'full_order_transported_projected_cubic_J_not_error_bound',
            'Checkpoint model/seed/scope identity', 'invalid_worker_reply')
    require(isinstance(data['julia_version'], str) and 1 <= len(data['julia_version']) <= 64
            and data['status'] in REPLY_STATUSES and isinstance(data['stop_reason'], str)
            and len(data['stop_reason']) <= 1024, 'Checkpoint status/runtime fields', 'invalid_worker_reply')
    score, response = finite(data['score'], 'Julia J'), finite(data['cubic_response'], 'Julia response')
    require(min(score, response) >= 0 and not (response > 0 and score == 0), 'Julia score sign/underflow', 'invalid_worker_reply')
    steps = context.affordable_steps(order, fields['budget'])
    require(type(data['steps']) is int and data['steps'] == steps and
        abs(score-response/steps**4) <= max(1e-30, abs(score)*1e-12), 'Julia score normalization', 'invalid_worker_reply')
    integer(data['start_scores_completed'], 'completed starts', 1, 4)
    integer(data['evaluations'], 'completed scores', 1, settings['max_evaluations'])
    for key, cap in (('proposals', settings['max_proposals']), ('accepted_moves', settings['max_proposals']),
                     ('uphill_attempts', settings['max_proposals']), ('uphill_accepted', settings['max_proposals'])):
        integer(data[key], key, 0, cap)
    for key in ('chain_proposals', 'chain_accepted'):
        require(isinstance(data[key], list) and len(data[key]) == 4, 'Four chain counters', 'invalid_worker_reply')
        for value in data[key]:
            integer(value, key, 0, settings['max_proposals'])
    completed = sum(data['chain_proposals'])
    require(data['evaluations'] == data['start_scores_completed']+completed and
        completed <= data['proposals'] <= completed+4 and
        data['accepted_moves'] == sum(data['chain_accepted']) <= completed and
        all(a <= p for a, p in zip(data['chain_accepted'], data['chain_proposals'])) and
        data['uphill_accepted'] <= min(data['uphill_attempts'], data['accepted_moves']) and
        data['uphill_attempts'] <= completed, 'Inconsistent checkpoint counters', 'invalid_worker_reply')
    require(isinstance(data['chain_CPU_seconds'], list) and len(data['chain_CPU_seconds']) == 4,
            'Four chain CPU fields', 'invalid_worker_reply')
    for value in data['chain_CPU_seconds']:
        require(finite(value, 'chain CPU') >= 0, 'Negative chain CPU', 'invalid_worker_reply')
    require(isinstance(data['phases'], dict) and len(data['phases']) <= 8, 'Phase fields', 'invalid_worker_reply')
    for value in data['phases'].values():
        require(finite(value, 'phase CPU') >= 0, 'Negative phase CPU', 'invalid_worker_reply')
    expected = {'unit_actions', 'cache_hits', 'cache_misses', 'cache_skipped', 'cache_entries', 'cache_accounted_bytes', 'kernel_hits',
        'kernel_misses', 'peak_vector_entries', 'peak_working_entries', 'full_scores'}
    require(isinstance(data['kernel_stats'], dict) and set(data['kernel_stats']) == expected,
            'Kernel counter fields', 'invalid_worker_reply')
    for key, value in data['kernel_stats'].items():
        integer(value, key, 0, 1_000_000_000)
    require(data['kernel_stats']['peak_vector_entries'] <= context.max_entries and
        data['kernel_stats']['peak_working_entries'] <= 12*context.max_entries and
        data['kernel_stats']['cache_entries'] <= context.max_cached_actions and
        data['kernel_stats']['cache_accounted_bytes'] <= settings['cache_bytes'] and
        data['kernel_stats']['cache_accounted_bytes'] == 1024*data['kernel_stats']['cache_entries'],
        'Kernel support exceeds declared caps', 'invalid_worker_reply')
    check()
    return order, data


def _reply(*args):
    try:
        return _validated_reply(*args)
    except Refusal as exc:
        if exc.code != 'invalid_input':
            raise
        raise Refusal('invalid_worker_reply', str(exc)) from exc
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise Refusal('invalid_worker_reply', str(exc)[:1024]) from exc


def train_annealing(blocks, n_qubits, reference, *, parent_order, starts, time,
        rotation_budget, parent_score, settings, seed_material, cpu_seconds,
        wall_seconds, max_entries=50000, max_cached_actions=100000, checkpoint_path=None):
    began, wall_began = _time.process_time(), _time.monotonic()
    stats = dict(method='julia_full_order_annealing', response_is_error_bound=False,
        physical_error_read=False, child_process_CPU_seconds=0., child_reaped=False,
        native_parent_revalidated=False, native_verifications=0, native_verification_attempts=0,
        startup_JIT_charged=True, final_score_scope='supplied_parent_not_revalidated_in_this_call',
        fallback=True, child_outcome='not_started')
    order = score = candidate = None
    refusal = None;status = 'not_started';promotion = 'checked_parent_fallback';eligible = False
    token = context = native = None
    allowance = wall_allowance = 0.

    def check(live_child=0.):
        cpu = _time.process_time()-began+max(live_child, stats['child_process_CPU_seconds'])
        if cpu >= allowance or _time.monotonic()-wall_began >= wall_allowance:
            raise Refusal('whole_call_limit', 'Complete annealing call allowance exhausted')

    try:
        allowance = finite(cpu_seconds, 'CPU allowance', positive=True)
        wall_allowance = finite(wall_seconds, 'wall allowance', positive=True)
        stats.update(requested_CPU_seconds=allowance, requested_wall_seconds=wall_allowance)
        native = _native()
        integer(n_qubits, 'register', 1, 14)
        require(type(reference) is int and 0 <= reference < 1 << n_qubits,
                'Single HF occupation integer required; coherent reference unsupported')
        max_entries = integer(max_entries, 'max_entries', 1, 50000)
        max_cached_actions = integer(max_cached_actions, 'max_cached_actions', 0, 100000)
        require(finite(time, 'training time') == 1., 'Training time must be one')
        integer(rotation_budget, 'rotation budget', 0, (1 << 62)-1)
        require(isinstance(blocks, (tuple, list)) and 1 <= len(blocks) <= 1024,
                'Bounded materialized original units required')
        blocks = tuple(blocks)
        order = native._permutation(parent_order, len(blocks))
        supplied = _score_dict(parent_score)
        score = native.Score(**supplied)
        settings = validate_settings(settings)
        _seeds(seed_material)
        require(isinstance(starts, (list, tuple)) and len(starts) == 4, 'Four starts required')
        starts = tuple(native._permutation(p, len(blocks)) for p in starts)
        require(len(set(starts)) == 4 and starts[0] == order, 'Four distinct parent-based starts required')
        token = native._EXTERNAL_BUDGET_CHECK.set(check)
        compiled, units = _compiled_units(native, blocks, n_qubits, check)
        require(rotation_budget == 30*sum(map(len, compiled)), 'Training budget must be original 30U')
        context = native.SparseResponse(blocks, n_qubits, reference, time=1.,
            max_entries=max_entries, max_cached_actions=max_cached_actions)
        phase = _time.process_time()
        actual = _checked_score(context.score(order, rotation_budget))
        parent_score_CPU = _time.process_time()-phase
        stats['native_parent_score_CPU_seconds'] = parent_score_CPU
        require(actual.steps == supplied['steps'] and _close(actual.value, supplied['value']) and
                _close(actual.cubic_response, supplied['cubic_response']), 'Frozen native parent score differs',
                'parent_score_mismatch')
        score = actual
        stats.update(native_parent_revalidated=True, final_score_scope='native_parent_revalidated')
        require(hasattr(os, 'wait4') and Path('/proc/self/stat').is_file()
                and signal.getsignal(signal.SIGCHLD) == signal.SIG_DFL,
                'Owned-PID Linux accounting unavailable', 'owned_accounting_unavailable')
        executable = shutil.which(os.environ.get('QHAT_JULIA_EXECUTABLE', 'julia'))
        require(executable is not None, 'Julia executable unavailable', 'julia_unavailable')
        stats['julia_executable'] = str(Path(executable).resolve())
        stats['julia_project'] = str(ROOT/'analysis/julia_window')
        sources = source_hashes(check)
        raw, fields = _request(units, n_qubits, reference, order, starts, actual, rotation_budget,
            settings, seed_material, max_entries, max_cached_actions, sources, check)
        digest = hashlib.sha256(raw).hexdigest()
        stats.update(request_sha256=digest, settings_sha256=fields['settings_sha256'],
            coefficients_sha256=fields['coefficients_sha256'], source_sha256=sources,
            seeds=fields['seeds'], settings=settings)
        reserve = max(settings['native_reserve_min'], settings['native_reserve_factor']*parent_score_CPU)
        cleanup = settings['cleanup_CPU_seconds']
        stats['verification_reserve_CPU_seconds'] = reserve

        def execute(request, checkpoint):
            nonlocal candidate, order, score, eligible, status, promotion
            require(not request.exists() and not checkpoint.exists() and
                not Path(str(checkpoint)+'.tmp').exists() and not Path(str(checkpoint)+'.refusal').exists(),
                'Fresh request/checkpoint paths required', 'existing_output')
            request.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            with request.open('xb') as stream:
                stream.write(raw)
            check()
            child_cap = allowance-(_time.process_time()-began)-reserve-cleanup
            child_wall = wall_allowance-(_time.monotonic()-wall_began)-reserve-cleanup
            require(child_cap > cleanup and child_wall > cleanup,
                    'Startup/setup exhausted verification/search allowance', 'verification_reserve_unavailable')
            stats.update(child_allocated_CPU_seconds=child_cap, child_allocated_wall_seconds=child_wall)

            def worker_check(live_child):
                check(live_child)
                require(_time.process_time()-began+live_child+reserve+cleanup < allowance,
                    'Owned child exhausted reserved CPU allowance', 'child_cpu_limit')
                require(_time.monotonic()-wall_began+reserve+cleanup < wall_allowance,
                    'Owned child exhausted reserved wall allowance', 'child_wall_limit')

            outcome = _run_worker(executable, request, checkpoint, child_cap, child_wall, worker_check, stats)
            stats['child_outcome'] = outcome
            if not checkpoint.is_file():
                raise Refusal('missing_worker_checkpoint', 'No complete Julia score checkpoint; '+outcome)
            proposed, reply = _reply(checkpoint, digest, fields, sources, settings, native, context, check)
            stats['child_startup_JIT_and_unattributed_CPU_seconds'] = max(0.,
                stats['child_process_CPU_seconds']-sum(reply['phases'].values())-sum(reply['chain_CPU_seconds']))
            stats['interrupted_proposal_attempts'] = reply['proposals']-sum(reply['chain_proposals'])
            candidate = dict(order=list(proposed), score=dict(value=reply['score'],
                cubic_response=reply['cubic_response'], steps=reply['steps']), witness=reply,
                checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest())
            stats['native_verification_attempts'] += 1
            phase = _time.process_time()
            checked = _checked_score(context.score(proposed, rotation_budget))
            stats['native_verification_CPU_seconds'] = _time.process_time()-phase
            check()
            stats['native_verifications'] += 1
            require(checked.steps == reply['steps'] and _close(checked.value, reply['score']) and
                    _close(checked.cubic_response, reply['cubic_response']),
                    'Native/Julia candidate score mismatch', 'native_julia_score_mismatch')
            candidate['native_score'] = dict(value=checked.value, cubic_response=checked.cubic_response, steps=checked.steps)
            stats['julia_version'] = reply['julia_version']
            if checked.value < score.value-max(1e-30, abs(score.value)*1e-12):
                order, score = proposed, checked
                promotion = 'native_verified_improvement'
                stats['fallback'] = False
                stats['final_score_scope'] = 'native_candidate_verified'
            else:
                promotion = 'native_verified_parent_retained'
                stats['fallback'] = False
            valid_outcomes = {'complete', 'child_cpu_limit', 'child_wall_limit', 'child_action_limit'}
            valid_statuses = {'complete', 'checked_update', 'slice_complete', 'cpu_limit', 'wall_limit',
                              'evaluation_limit', 'proposal_limit', 'action_limit'}
            coverage = reply['start_scores_completed'] == 4 and all(p > 0 for p in reply['chain_proposals'])
            eligible = outcome in valid_outcomes and reply['status'] in valid_statuses and coverage
            stats['four_chain_coverage'] = coverage
            status = reply['status'] if eligible else 'ineligible_checked_incumbent'
            if not eligible:
                raise Refusal('ineligible_method_execution', 'Incomplete starts/chains or refused backend: '+outcome+'/'+reply['status'])

        if checkpoint_path is None:
            with tempfile.TemporaryDirectory(prefix='qhat-full-annealing-') as temporary:
                checkpoint = Path(temporary)/'checkpoint.toml'
                execute(Path(temporary)/'request.toml', checkpoint)
        else:
            checkpoint = Path(checkpoint_path).resolve()
            execute(Path(str(checkpoint)+'.request.toml'), checkpoint)
    except Exception as exc:
        eligible = False
        code = exc.code if isinstance(exc, Refusal) else (
            'numerical_limit' if native is not None and isinstance(exc, native.NumericalResolution) else
            'support_limit' if native is not None and isinstance(exc, native.ResourceLimit) else
            'whole_call_limit' if native is not None and isinstance(exc, native.SearchDeadline) else 'method_failure')
        status = code
        refusal = dict(type=type(exc).__name__, code=code, message=str(exc)[:1024])
    finally:
        if native is not None and token is not None:
            native._EXTERNAL_BUDGET_CHECK.reset(token)
        context = None
        stats['parent_process_CPU_seconds'] = _time.process_time()-began
        stats['observed_CPU_seconds'] = stats['parent_process_CPU_seconds']+stats['child_process_CPU_seconds']
        stats['wall_seconds'] = _time.monotonic()-wall_began
        if allowance and (stats['observed_CPU_seconds'] > allowance or stats['wall_seconds'] > wall_allowance):
            eligible = False;status = 'whole_call_overrun'
            refusal = dict(type='Refusal', code=status, message='Complete observed call exceeds requested allowance')
    return AnnealingResult(order, score, bool(eligible), status, refusal, candidate, promotion,
        stats['observed_CPU_seconds'], stats)
