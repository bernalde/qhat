"""Optional bounded Julia conditioned-window insertion for the ordinary trainer.

Only process control, input admission and serialization live here. The numerical
window construction/search is shipped as versioned Julia stdlib assets. A native
full score is authoritative before a proposal replaces the checked parent.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time as _time

from qhat.analysis import sparse_response_ordering as native


@dataclass(frozen=True)
class WindowOptions:
    width: int = 8
    max_windows: int = 1
    max_sweeps: int = 10
    parent_method: str = 'sparse'
    parent_fraction: float = .60
    julia_executable: str | None = None
    wall_seconds: float | None = None
    max_factor_entries: int = 250000
    max_flat_entries: int = 1000000
    max_actions: int = 2000000
    max_atoms: int = 200000

    def __post_init__(self):
        for name, low, high in (('width', 2, 9), ('max_windows', 1, 32),
                                ('max_sweeps', 1, 100), ('max_factor_entries', 1, 250000),
                                ('max_flat_entries', 1, 1000000), ('max_actions', 1, 10000000),
                                ('max_atoms', 1, 1000000)):
            native._integer(getattr(self, name), name, low, high)
        if self.parent_method not in ('sparse', 'magnitude'):
            raise ValueError('window parent_method must be sparse or magnitude')
        fraction = native._finite(self.parent_fraction, 'parent_fraction', positive=True)
        if fraction > .60:
            raise ValueError('window parent_fraction must not exceed .60')
        if self.wall_seconds is not None:
            native._finite(self.wall_seconds, 'wall_seconds', positive=True)
        if self.julia_executable is not None and (not isinstance(self.julia_executable, str)
                                                 or not self.julia_executable):
            raise ValueError('julia_executable must be a nonempty path/name')


def asset_path():
    """Installed-wheel location, independent of the working directory."""
    return Path(__file__).resolve().parent / 'julia_window' / 'worker.jl'


def _parent_evaluation_cap(total, fraction):
    """Reserve two native scores and the first three child prefix evaluations.

    Fewer than six evaluations fund only the parent. Otherwise apply the same
    parent fraction as the CPU allocation, while retaining at least one parent
    score and the child parent/flat-baseline/first-trial prefix.
    """
    if total < 6:
        return total
    return min(max(1, math.floor(fraction * (total - 2))), total - 2 - 3)


def _checked_score(score):
    if score is not None and score.cubic_response > 0 and score.value == 0:
        raise native.NumericalResolution('window normalized native score underflows')
    return score


def _fraction(value):
    value = Fraction(value)
    if value.numerator.bit_length() > 4096 or value.denominator.bit_length() > 4096:
        raise native.ResourceLimit('window rational exceeds 4096-bit admission cap')
    return f'{value.numerator}/{value.denominator}'


def _toml(value):
    # Schema-specific scalar/array writer. JSON string escaping is compatible
    # with these ASCII-only controlled protocol strings; no source interpolation.
    if type(value) is int:
        return str(value)
    if type(value) is float and math.isfinite(value):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    if isinstance(value, (tuple, list)):
        return '[' + ','.join(_toml(v) for v in value) + ']'
    raise ValueError('unsupported window TOML value')


def _request(compiled, context, order, score, budget, options, evaluations, check):
    fields = dict(version=1, n=context.n_qubits, reference=next(iter(context.reference)),
                  time=_fraction(context.time), budget=budget, order=order,
                  word_lengths=context.lengths, parent_score=float(score.value),
                  width=options.width, max_windows=options.max_windows, max_sweeps=options.max_sweeps,
                  max_entries=context.max_entries, max_cached_actions=context.max_cached_actions,
                  max_factor_entries=options.max_factor_entries, max_flat_entries=options.max_flat_entries,
                  max_actions=options.max_actions, max_atoms=options.max_atoms,
                  max_evaluations=evaluations)
    lines = []
    for key, value in fields.items():
        check()
        lines.append(key + '=' + _toml(value))
    for unit in compiled:
        check()
        lines.append('[[units]]')
        lines.append('flip=' + _toml(unit[0][0]))
        parities, signs, coefficients = [], [], []
        for _, z, phase, coefficient in unit:
            check()
            parities.append(z)
            # phase includes the signed coefficient; encode only i^Y.
            signs.append(1 if (phase.real > 0) == (coefficient > 0) else -1)
            coefficients.append(_fraction(coefficient))
        lines.append('parities=' + _toml(parities))
        lines.append('signs=' + _toml(signs))
        lines.append('coefficients=' + _toml(coefficients))
    raw = ('\n'.join(lines) + '\n').encode('ascii')
    if len(raw) > 8 * 1024**2:
        raise native.ResourceLimit('window request exceeds 8 MiB')
    return raw


def _domain(blocks, compiled, context, check):
    if not 1 <= context.n_qubits <= 20 or not 2 <= len(compiled) <= 1024:
        return 'unsupported_register_or_unit_cap'
    if sum(map(len, blocks)) > 8192:
        return 'unsupported_word_cap'
    if len(context.reference) != 1:
        return 'unsupported_coherent_reference'
    if context.max_entries > 50000 or context.max_cached_actions > 100000:
        return 'unsupported_storage_cap'
    for original, unit in zip(blocks, compiled):
        check()
        if not unit or len(original) != len(unit):
            return 'unsupported_empty_or_zero_word_unit'
        if len({x for x, _, _, _ in unit}) != 1 or any(phase.imag for _, _, phase, _ in unit):
            return 'unsupported_complex_or_mixed_flip_unit'
    return None


def _proc_cpu(pid):
    """Linux owned-PID live user+system CPU; never global RUSAGE_CHILDREN."""
    try:
        text = Path(f'/proc/{pid}/stat').read_text()
        fields = text[text.rfind(')') + 2:].split()
        return (int(fields[11]) + int(fields[12])) / os.sysconf('SC_CLK_TCK')
    except (OSError, ValueError, IndexError):
        return None


def _run_worker(executable, request, checkpoint, child_cap, wall_cap, check, stats):
    """Own exactly one child, watchdog startup/JIT, then obtain its wait4 receipt.

    Polling /proc enforces a CPU deadline during uninterruptible startup/JIT.
    wait4 is the authoritative final user+system receipt for this exact PID.
    The child is killed/reaped before any return, including a Python exception.
    """
    worker = asset_path()
    if not worker.is_file():
        raise native.ResourceLimit('installed Julia worker assets are missing')
    argv = [executable, '--startup-file=no', '--history-file=no', '--threads=1', '--gcthreads=1',
            '--heap-size-hint=250M', '--project=' + str(worker.parent),
            str(worker), str(request), str(checkpoint), repr(child_cap), repr(wall_cap)]
    env = os.environ.copy()
    env.update(JULIA_NUM_THREADS='1', JULIA_NUM_GC_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', JULIA_LOAD_PATH='@:@stdlib',
               JULIA_PKG_PRECOMPILE_AUTO='0', JULIA_IMAGE_THREADS='1')
    began = _time.monotonic()
    process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, env=env, start_new_session=True)
    reaped = False
    termination = None
    observed = 0.
    try:
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                process.returncode = os.waitstatus_to_exitcode(status)
                reaped = True
                stats['child_process_CPU_seconds'] = usage.ru_utime + usage.ru_stime
                stats['child_peak_RSS_KiB'] = usage.ru_maxrss
                stats['child_returncode'] = process.returncode
                break
            current = _proc_cpu(process.pid)
            if current is None:
                termination = 'owned_child_cpu_unavailable'
            else:
                observed = max(observed, current)
                if observed >= child_cap:
                    termination = 'child_cpu_limit'
            if _time.monotonic() - began >= wall_cap:
                termination = 'child_wall_limit'
            try:
                check(observed)
            except native.SearchDeadline:
                termination = 'whole_call_limit'
            if termination:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            _time.sleep(.005)
    finally:
        if not reaped:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            stats['child_process_CPU_seconds'] = usage.ru_utime + usage.ru_stime
            stats['child_peak_RSS_KiB'] = usage.ru_maxrss
            stats['child_returncode'] = process.returncode
        stats['child_wall_seconds'] = _time.monotonic() - began
    if not termination and process.returncode:
        refusal_path=Path(str(checkpoint)+'.refusal')
        if refusal_path.is_file():
            with refusal_path.open('rb') as handle:
                refusal=handle.read(128)
            if refusal in {b'cpu_limit',b'wall_limit',b'evaluation_limit',b'action_limit',b'atom_limit',
                           b'numerical_limit',b'admission_or_support_limit',b'worker_failure'}:
                stats['child_refusal']=refusal.decode('ascii')
                return 'child_'+stats['child_refusal']
    return termination or ('child_failure' if process.returncode else 'complete')


def _reply(path, digest, context, baseline, options, budget, check):
    try:
        import tomllib
    except ImportError as exc:
        raise native.ResourceLimit('window backend requires Python 3.11+ tomllib') from exc
    check()
    with path.open('rb') as handle:
        raw = handle.read(64*1024 + 1)
    if len(raw) > 64*1024:
        raise ValueError('window reply byte cap')
    # The reply schema needs only shallow, short ASCII tables/arrays. Refuse
    # deep nesting before TOML parsing so malformed workers cannot exhaust the
    # embedding Python stack or force an unbounded parser allocation.
    if (not raw.isascii() or len(raw.splitlines()) > 512
            or b'"""' in raw or b"'''" in raw):
        raise ValueError('window reply text cap')
    nesting = 0
    quote = 0
    escaped = False
    for at, byte in enumerate(raw):
        if at % 1024 == 0:
            check()
        if escaped:
            escaped = False
        elif quote == 34 and byte == 92:
            escaped = True
        elif quote:
            if byte == quote:
                quote = 0
        elif byte in (34,39):
            quote = byte
        elif byte == 91:
            nesting += 1
            if nesting > 3:
                raise ValueError('window reply nesting cap')
        elif byte == 93:
            nesting -= 1
            if nesting < 0:
                raise ValueError('window reply nesting differs')
    check()
    data = tomllib.loads(raw.decode('utf-8'))
    expected = {'version', 'request_sha256', 'order', 'word_lengths', 'score', 'cubic_response',
                'steps', 'status', 'checked_updates', 'windows', 'history', 'counters', 'phases'}
    if set(data) != expected or type(data['version']) is not int or data['version'] != 1:
        raise ValueError('invalid window reply schema')
    if data['request_sha256'] != digest:
        raise ValueError('window reply request hash differs')
    order = native._permutation(data['order'], len(baseline))
    if (not isinstance(data['word_lengths'], list)
            or any(type(i) is not int for i in data['word_lengths'])
            or data['word_lengths'] != list(context.lengths)):
        raise ValueError('window reply word lengths differ')
    windows = data['windows']
    if not isinstance(windows, list) or len(windows) > options.max_windows:
        raise ValueError('invalid reply windows')
    changed = set()
    for window in windows:
        check()
        if (not isinstance(window, list) or len(window) != 2
                or any(type(i) is not int for i in window)):
            raise ValueError('invalid reply window')
        lo, hi = window
        if not 1 <= lo <= hi < len(order) - 1 or hi - lo + 1 != options.width:
            raise ValueError('reply window is not the requested interior width')
        changed.update(range(lo, hi + 1))
    if any(order[i] != baseline[i] for i in range(len(order)) if i not in changed):
        raise ValueError('window reply changed its fixed exterior')
    history=data['history']
    if not isinstance(history,list) or len(history)>options.max_windows:
        raise ValueError('invalid window history')
    replay=list(baseline); replay_windows=[]
    for entry in history:
        check()
        if not isinstance(entry,dict) or set(entry)!={'positions','baseline_labels','proposal_labels','accepted'}:
            raise ValueError('invalid window history fields')
        interval=entry['positions']
        if (not isinstance(interval,list) or len(interval)!=2
                or any(type(i) is not int for i in interval)):
            raise ValueError('invalid history interval')
        lo,hi=interval
        if not 1<=lo<=hi<len(order)-1 or hi-lo+1!=options.width:
            raise ValueError('history interval is not the requested interior width')
        before,after=entry['baseline_labels'],entry['proposal_labels']
        if (not isinstance(before,list) or not isinstance(after,list)
                or any(type(i) is not int for i in before+after)
                or before!=replay[lo:hi+1] or sorted(before)!=sorted(after)
                or type(entry['accepted']) is not bool):
            raise ValueError('history labels differ from its checked baseline')
        if entry['accepted']:
            if before==after:
                raise ValueError('history accepted an unchanged window')
            replay[lo:hi+1]=after
            replay_windows.append(interval)
    if replay!=list(order) or replay_windows!=windows:
        raise ValueError('window reply is not its reported insertion history')
    if (type(data['checked_updates']) is not int or data['checked_updates'] != len(windows)
            or data['status'] not in ('parent_checked', 'checked_update', 'complete')):
        raise ValueError('invalid checkpoint status')
    score = native._finite(data['score'], 'reply score')
    cubic = native._finite(data['cubic_response'], 'reply cubic response')
    if score < 0 or cubic < 0:
        raise ValueError('negative reply score')
    steps = context.affordable_steps(order, budget)
    if type(data['steps']) is not int or data['steps'] != steps:
        raise ValueError('window reply pulse cost differs')
    if abs(score - cubic / steps**4) > max(1e-30, abs(score)*1e-12):
        raise ValueError('window reply normalization differs')
    if not isinstance(data['counters'], dict) or set(data['counters']) != {
            'evaluation', 'action', 'atom', 'full_score', 'flat_candidate', 'exact_witness', 'insertion_move'}:
        raise ValueError('invalid reply counters')
    for key, cap in (('evaluation', 1000000), ('action', options.max_actions), ('atom', options.max_atoms),
                     ('full_score', options.max_windows+1), ('flat_candidate', 1000000),
                     ('exact_witness', 1000000), ('insertion_move', options.max_sweeps*options.max_windows)):
        native._integer(data['counters'][key], 'reply ' + key, 0, cap)
    if sum(data['counters'][key] for key in ('full_score','flat_candidate','exact_witness'))>data['counters']['evaluation']:
        raise ValueError('reply scoring counters exceed evaluation count')
    if not isinstance(data['phases'], dict) or len(data['phases']) > 16:
        raise ValueError('invalid reply phases')
    for value in data['phases'].values():
        if native._finite(value, 'reply phase') < 0:
            raise ValueError('negative reply phase')
    check()
    return order, data


def train_window(*args, **kwargs):
    """Run an isolated call-scoped budget through the existing native scorer."""
    token = native._EXTERNAL_BUDGET_CHECK.set(None)
    try:
        return _train_window(*args, **kwargs)
    finally:
        native._EXTERNAL_BUDGET_CHECK.reset(token)


def _train_window(blocks, n_qubits, reference, *, time=1., rotation_budget,
                 cpu_seconds=2., max_evaluations=3000, max_moves=5, max_entries=50000,
                 max_cached_actions=50000, starts=None, window_options=None,
                 _began=None, _wall_began=None):
    began = _time.process_time() if _began is None else _began
    wall_began = _time.monotonic() if _wall_began is None else _wall_began
    allowance = native._finite(cpu_seconds, 'cpu_seconds', positive=True)
    native._integer(max_evaluations, 'max_evaluations', 0, 1000000)
    native._integer(max_moves, 'max_moves', 0, 1000)
    native._integer(rotation_budget, 'rotation_budget', 0, (1 << 62) - 1)
    options = WindowOptions() if window_options is None else window_options
    if not isinstance(options, WindowOptions):
        raise ValueError('window_options must be WindowOptions')
    child_evaluation_eligible = max_evaluations >= 6
    parent_evaluation_cap = _parent_evaluation_cap(max_evaluations, options.parent_fraction)
    wall_allowance = native._finite(options.wall_seconds if options.wall_seconds is not None
        else max(1.,4*allowance), 'window wall_seconds', positive=True)
    stats = dict(method='window', backend_version=1, accounting_version='owned-wait4-v1',
                 construction_scope='fresh_parent_total', parent_method=options.parent_method,
                 window_status='not_started', response_is_error_bound=False,
                 requested_CPU_seconds=allowance, requested_wall_seconds=wall_allowance,
                 child_process_CPU_seconds=0., verification_CPU_seconds=0.,
                 requested_evaluations=max_evaluations, parent_allocation_fraction=float(options.parent_fraction),
                 parent_allocated_evaluations=parent_evaluation_cap,
                 native_reserved_evaluations=2 if child_evaluation_eligible else 0,
                 native_parent_revalidation_evaluations=0,
                 native_final_reserved_evaluations=1 if child_evaluation_eligible else 0,
                 child_minimum_evaluations=3, child_allocated_evaluations=0,
                 child_completed_evaluations=None,
                 evaluation_accounting_scope='Completed native scores plus authenticated child checkpoint prefix; later uncheckpointed child work can remain unknown',
                 serialization_CPU_seconds=0., cleanup_CPU_seconds=0.,
                 native_verification_attempts=0, native_verifications=0,
                 selection='max_nonzero_pair', aggregation='real_columns_before_transport_v1')
    stats['storage_caps'] = dict(entries=max_entries, cached_actions=max_cached_actions,
        factor_entries=options.max_factor_entries, flat_entries=options.max_flat_entries,
        units=1024, words=8192, register_bits=20, coefficient_bits=4096,
        intermediate_rational_bits=16384, exact_witness_bits=8192)
    def check(live_child=0.):
        charged = _time.process_time()-began + max(live_child, stats['child_process_CPU_seconds'])
        if charged >= allowance or _time.monotonic()-wall_began >= wall_allowance:
            raise native.SearchDeadline('whole-call window allowance exhausted')
    def deadline():
        check()
        return min(began + allowance-stats['child_process_CPU_seconds'],
                   _time.process_time()+max(0.,wall_allowance-(_time.monotonic()-wall_began)))
    native._EXTERNAL_BUDGET_CHECK.set(check)
    # Native validation remains authoritative. Validation never silently thins a
    # term to fit the Julia domain; unsupported valid inputs keep a native parent.
    blocks = tuple(blocks)
    native._finite(time, 'time', positive=True)
    native._integer(max_entries, 'max_entries', 1, 1000000)
    native._integer(max_cached_actions, 'max_cached_actions', 0, 1000000)
    reference = native._reference(reference, native._integer(n_qubits, 'n_qubits', 1, 4096), max_entries)
    compiled = native._compile(blocks, n_qubits, deadline())
    canonical = tuple(range(len(compiled)))
    magnitude = tuple(sorted(canonical, key=lambda i: (-sum(abs(c) for _, _, c, _ in compiled[i]), i)))
    def affordable(order):
        lengths = [len(compiled[i]) for i in order if compiled[i]]
        return not lengths or native.plan_strang_schedule(lengths, rotation_budget).feasible
    seeds = [magnitude] if starts is None else [native._permutation(s, len(compiled)) for s in starts]
    if not seeds:
        raise ValueError('at least one seed is required')
    seeds = [s for s in seeds if affordable(s)]
    if not seeds:
        seeds = [tuple(sorted(canonical, key=lambda i: (len(compiled[i]), i)))]
        if not affordable(seeds[0]):
            raise ValueError('rotation budget cannot afford one symmetric step')
    parent = native.TrainingResult(seeds[0], None, 0., 0, 0, 'cpu_limit', {})
    best_order, best_score = parent.order, None
    evaluations = moves = 0
    context = None
    try:
        check()
        parent_allowance = allowance*options.parent_fraction-(_time.process_time()-began)
        stats['parent_allocated_CPU_seconds']=max(0.,parent_allowance)
        if parent_allowance<=0:
            raise native.SearchDeadline("parent setup exhausted its allocation")
        parent = native.train_ordering(blocks, n_qubits, reference, time=time,
            rotation_budget=rotation_budget, cpu_seconds=parent_allowance,
            max_evaluations=parent_evaluation_cap,
            max_moves=0 if options.parent_method=='magnitude' else max_moves,
            max_entries=max_entries,max_cached_actions=max_cached_actions,
            starts=seeds if options.parent_method=='magnitude' or starts is not None else None)
        best_order = parent.order
        evaluations, moves = parent.evaluations, parent.accepted_moves
        stats.update(parent_CPU_seconds=parent.CPU_seconds, parent_status=parent.status,
                     parent_evaluations=parent.evaluations, parent_stats=parent.stats,
                     parent_evaluation_cap_saturated=bool(parent_evaluation_cap and
                                                        parent.evaluations==parent_evaluation_cap))
        best_score = _checked_score(parent.score)
        if best_score is None:
            stats['window_status']='parent_unscored'
        elif not child_evaluation_eligible:
            stats['window_status']='insufficient_evaluation_budget'
        elif max_evaluations-evaluations < 2+stats['child_minimum_evaluations']:
            stats['window_status']='evaluation_limit'
        else:
            # Reuse the existing native kernel to establish a measured verification
            # reserve and keep its initialized caches for the final comparison.
            phase=_time.process_time()
            context=native.SparseResponse(blocks,n_qubits,reference,time=time,max_entries=max_entries,
                max_cached_actions=max_cached_actions,deadline=deadline())
            score_began=_time.process_time()
            verified=_checked_score(context.score(best_order,rotation_budget,deadline()))
            score_time=_time.process_time()-score_began
            evaluations+=1
            stats['native_parent_revalidation_evaluations']=1
            stats['parent_revalidation_CPU_seconds']=_time.process_time()-phase
            best_score=verified
            refused=_domain(blocks,compiled,context,check)
            if refused:
                stats['window_status']=refused
            elif len(compiled)<options.width+2:
                stats['window_status']='window_unavailable'
            elif sys.version_info<(3,11):
                stats['window_status']='unsupported_toml_reader'
            elif not hasattr(os,'wait4') or not Path('/proc/self/stat').is_file():
                stats['window_status']='unsupported_owned_cpu_accounting_platform'
            elif signal.getsignal(signal.SIGCHLD)!=signal.SIG_DFL:
                stats['window_status']='unsupported_owned_cpu_reaper'
            else:
                executable=options.julia_executable or os.environ.get('QHAT_JULIA_EXECUTABLE') or 'julia'
                executable=shutil.which(executable)
                if executable is None:
                    stats['window_status']='julia_unavailable'
                else:
                    reserve=max(.02,1.25*score_time)
                    cleanup_reserve=.02
                    remaining=allowance-(_time.process_time()-began)
                    stats['verification_reserve_CPU_seconds']=reserve
                    stats['verification_reserve_basis']='native_parent_revalidation_full_score'
                    if remaining<=reserve+cleanup_reserve:
                        stats['window_status']='verification_reserve_unavailable'
                    else:
                        phase=_time.process_time()
                        child_evaluations=max_evaluations-evaluations-1
                        raw=_request(compiled,context,best_order,best_score,rotation_budget,options,
                                     child_evaluations,check)
                        digest=hashlib.sha256(raw).hexdigest()
                        stats['request_sha256']=digest
                        with (asset_path().parent/'SOURCE.toml').open('rb') as handle:
                            manifest=handle.read(8193)
                        if len(manifest)>8192:
                            raise native.ResourceLimit('backend source manifest byte cap')
                        stats['backend_manifest_sha256']=hashlib.sha256(manifest).hexdigest()
                        check()
                        with tempfile.TemporaryDirectory(prefix='qhat-window-') as temporary:
                            request=Path(temporary)/'request.toml';checkpoint=Path(temporary)/'checkpoint.toml'
                            request.write_bytes(raw);check()
                            stats['serialization_CPU_seconds']=_time.process_time()-phase
                            child_cap=allowance-(_time.process_time()-began)-reserve-cleanup_reserve
                            child_wall=wall_allowance-(_time.monotonic()-wall_began)-reserve-cleanup_reserve
                            if child_cap<=0 or child_wall<=0:
                                stats['window_status']='verification_reserve_unavailable'
                            else:
                                stats['child_allocated_CPU_seconds']=child_cap
                                stats['child_allocated_evaluations']=child_evaluations
                                def worker_check(live_child):
                                    check(live_child)
                                    if (_time.process_time()-began+live_child+reserve+cleanup_reserve>=allowance
                                            or _time.monotonic()-wall_began+reserve+cleanup_reserve>=wall_allowance):
                                        raise native.SearchDeadline('preserve verification reserve')
                                outcome=_run_worker(executable,request,checkpoint,child_cap,child_wall,worker_check,stats)
                                stats['window_status']=outcome
                                if checkpoint.is_file():
                                    try:
                                        candidate,reply=_reply(checkpoint,digest,context,best_order,options,
                                                               rotation_budget,check)
                                        if reply['counters']['evaluation']>child_evaluations:
                                            raise ValueError('reply exceeds allocated evaluation cap')
                                        evaluations+=reply['counters']['evaluation']
                                        stats['child_completed_evaluations']=reply['counters']['evaluation']
                                        stats['child_evaluation_scope']='checkpoint_prefix'
                                        stats['child_counters']=reply['counters'];stats['child_phases']=reply['phases']
                                        stats['child_startup_JIT_and_unattributed_CPU_seconds']=max(0.,
                                            stats['child_process_CPU_seconds']-sum(reply['phases'].values()))
                                        stats['child_checkpoint_status']=reply['status']
                                        stats['child_checked_updates']=reply['checked_updates']
                                        stats['window_history']=reply['history']
                                        if candidate!=best_order:
                                            phase=_time.process_time()
                                            stats['native_verification_attempts']+=1
                                            try:
                                                actual=_checked_score(context.score(candidate,rotation_budget,deadline()))
                                                evaluations+=1
                                                check()
                                                stats['native_verifications']+=1
                                            finally:
                                                stats['verification_CPU_seconds']+=_time.process_time()-phase
                                            if (abs(actual.value-reply['score'])>max(1e-28,abs(actual.value)*1e-8)):
                                                stats['window_status']='native_julia_score_mismatch'
                                            elif actual.value<best_score.value-max(1e-30,abs(best_score.value)*1e-12):
                                                best_order,best_score=candidate,actual;moves+=1
                                                stats['window_status']='verified_improvement' if outcome=='complete' else 'verified_checkpoint_after_'+outcome
                                            else:
                                                stats['window_status']='native_no_resolved_improvement'
                                    except (ValueError,KeyError,TypeError,UnicodeError) as exc:
                                        stats['window_status']='invalid_worker_reply'
                                        stats['refusal_reason']=str(exc)[:240]
                                elif outcome=='complete':
                                    stats['window_status']='missing_worker_reply'
    except native.SearchDeadline:
        stats['window_status']='whole_call_limit'
    except native.NumericalResolution:
        stats['window_status']='numerical_limit'
    except native.ResourceLimit as exc:
        stats['window_status']='support_limit';stats['refusal_reason']=str(exc)[:240]
    except OSError as exc:
        stats['window_status']='child_launch_or_transport_failure';stats['refusal_reason']=str(exc)[:240]
    finally:
        cleanup=_time.process_time()
        context=None
        compiled=None
        blocks=()
        stats['cleanup_CPU_seconds']+=_time.process_time()-cleanup
    stats['wall_seconds']=_time.monotonic()-wall_began
    stats['parent_process_CPU_seconds']=_time.process_time()-began
    charged=stats['parent_process_CPU_seconds']+stats['child_process_CPU_seconds']
    stats['observed_CPU_seconds']=charged
    status = 'complete' if stats['window_status'] in ('complete','verified_improvement','native_no_resolved_improvement') else parent.status
    if best_score is None and stats['window_status']=='numerical_limit':
        status='numerical_limit'
    return native.TrainingResult(best_order,best_score,charged,evaluations,moves,status,stats)
