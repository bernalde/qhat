"""OLD-only, finite temperature calibration; no chemistry or physical errors.

CLI: calibration.py --output FRESH_DIRECTORY --cpu-cap SECONDS --wall-cap SECONDS
                   [--julia-executable PATH]

Only the two retained, byte-pinned CO/CH4 inputs below are accepted. Thirty-two
unrestricted insertion draws per case are declared before any scoring, including
repeated/nonpositive probes. Native, new Julia and archived Julia scores must
agree under the existing objective tolerance/floor rule. No favorable-order or
physical-error search chooses a temperature. The owning OS guard supplies the
authoritative startup-through-teardown cost; exact wait4 accounts for our child.
"""
from __future__ import annotations

import time as _time
BOOT_WALL = _time.monotonic()
import argparse
from fractions import Fraction
import hashlib
import importlib.machinery
import importlib.util
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROBES = 32
RTOL, ATOL, ENERGY_FLOOR = 1e-8, 1e-24, 1e-30
P_START, P_END = .8, .01
MIN_POSITIVE = 1
FINAL_CPU_RESERVE, FINAL_WALL_RESERVE = .5, 3.
OLD = 'numerical_audit/julia_cubic_feasibility/evidence/inputs/cubic/'
CASES = (
    dict(case='CO_r1p13', path=OLD+'CO_r1p13.toml', n=12, electrons=6, units=195, words=770,
         sha256='b78921326c04f331a0d95b8fadf1be600bb6bfd0caa65e464d24c7aa6de5fe52',
         legacy_source_model_sha256='8a9f89166544f72f9ecae892ad664a3f110f1a2c534fceef22359aaa1fcbdbe8'),
    dict(case='CH4_r1p09', path=OLD+'CH4_r1p09.toml', n=12, electrons=6, units=351, words=1722,
         sha256='00a283f3736e331b2c3bfc140213e1d5d1d3a62f119b34ab1450a1b091afbc97',
         legacy_source_model_sha256='32aeb49b2e541088efa0604e90dae046e8be58e81259e027032b65fe2d779e37'))
SOURCE_NAMES = (
    'analysis/sparse_response_ordering.py', 'analysis/strang_schedule.py',
    'analysis/julia_window/ResponseCore.jl', 'analysis/julia_window/Project.toml',
    'numerical_audit/molecular_generalization/src/FullOrderAnnealing.jl',
    'numerical_audit/molecular_generalization/anneal.jl',
    'numerical_audit/molecular_generalization/annealing.py',
    'numerical_audit/molecular_generalization/calibration.py',
    'numerical_audit/molecular_generalization/calibration.jl',
    'numerical_audit/prospective_window_ordering/objective_check.jl',
    'numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl',
    'numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl')


class Refusal(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Refusal(message)


def bounded_read(path, maximum=8*1024**2):
    require(path.is_file() and path.stat().st_size <= maximum, 'Bounded input required: '+str(path))
    with path.open('rb') as stream:
        raw = stream.read(maximum+1)
    require(len(raw) <= maximum, 'Growing input byte cap')
    return raw


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write_new(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.partial')
    require(not path.exists() and not temporary.exists(), 'Fresh artifact required: '+str(path))
    raw = (json.dumps(data, indent=2, allow_nan=False)+'\n').encode('ascii')
    require(len(raw) <= 32*1024**2, 'Calibration JSON byte cap')
    with temporary.open('xb') as stream:
        stream.write(raw)
    os.link(temporary, path)
    temporary.unlink()
    return digest(raw)


def adapter():
    """Standalone bootstrap, followed by the adapter's exact native identity check."""
    if 'qhat' not in sys.modules:
        spec = importlib.machinery.ModuleSpec('qhat', loader=None, is_package=True)
        spec.submodule_search_locations = [str(ROOT)]
        sys.modules['qhat'] = importlib.util.module_from_spec(spec)
    name = '_molecular_old_calibration_adapter'
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, HERE/'annealing.py')
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    module = sys.modules[name]
    require(Path(module.__file__).resolve() == HERE/'annealing.py', 'Calibration adapter source identity')
    return module, module._native()


def read_model(path, expected):
    """Decode original rational coefficients; retain the legacy header unchanged."""
    raw = bounded_read(path)
    require(digest(raw) == expected['sha256'], 'Retained old input bytes differ')
    data = tomllib.loads(raw.decode('ascii'))
    require(data['case'] == expected['case'] and data['n'] == expected['n'] and
            type(data['n']) is int and data['electrons'] == expected['electrons'] and
            type(data['electrons']) is int, 'Fixed old model identity')
    require(data.get('source_model_sha256') == expected['legacy_source_model_sha256'] and
            'source_json_sha256' not in data, 'Preserve legacy source_model_sha256 provenance')
    denominator = data['denominator']
    require(isinstance(denominator, str) and len(denominator) <= 1280 and denominator.isdigit(), 'Denominator spelling')
    denominator = int(denominator)
    require(0 < denominator and denominator.bit_length() <= 4096, 'Positive bounded denominator')
    units = data['units']
    require(isinstance(units, list) and len(units) == expected['units'] and len(units) <= 1024, 'Original unit inventory')
    blocks, seen = [], set()
    for unit in units:
        words, numerators = unit['words'], unit['numerators']
        require(isinstance(words, list) and isinstance(numerators, list) and
                1 <= len(words) == len(numerators) <= 8192, 'Original word arrays')
        block = {}
        for text, numerator in zip(words, numerators):
            require(isinstance(text, str) and len(text) <= 4*data['n'], 'Bounded Pauli spelling')
            require(isinstance(numerator, str) and len(numerator) <= 1280 and
                    re.fullmatch(r'-?[0-9]+', numerator), 'Numerator spelling')
            number = int(numerator)
            require(number != 0 and number.bit_length() <= 4096, 'Nonzero bounded coefficient')
            word = []
            for token in ([] if not text else text.split(',')):
                match = re.fullmatch(r'([0-9]+)([XYZ])', token)
                require(match is not None, 'Pauli token')
                site, letter = int(match.group(1)), match.group(2)
                require(0 <= site < data['n'], 'Pauli site outside register')
                word.append((site, letter))
            word = tuple(word)
            require(word == tuple(sorted(word)) and len({q for q, _ in word}) == len(word)
                    and word not in seen, 'Canonical distinct original Pauli labels')
            seen.add(word);block[word] = Fraction(number, denominator)
        blocks.append(block)
    require(len(seen) == expected['words'] and len(seen) <= 8192, 'Original word inventory')
    parent = data['cubic_order']
    require(isinstance(parent, list) and all(type(x) is int for x in parent) and
            sorted(parent) == list(range(1, len(blocks)+1)), 'Retained complete parent permutation')
    return dict(data=data, raw=raw, blocks=tuple(blocks), parent=tuple(x-1 for x in parent),
        n=data['n'], reference=(1 << data['electrons'])-1, words=len(seen))


def probe_plan(case, parent):
    require(case in {spec['case'] for spec in CASES}, 'Only predeclared old cases')
    require(isinstance(parent, (list, tuple)) and 2 <= len(parent) <= 1024 and
            all(type(x) is int for x in parent) and sorted(parent) == list(range(len(parent))), 'Parent permutation')
    material = case+':molecular-generalization-sa-v1:calibration-v1'
    seed = int.from_bytes(hashlib.sha256(material.encode('ascii')).digest()[:8], 'big')
    rng = random.Random(seed)
    pairs, orders = [], [list(parent)]
    for _ in range(PROBES):
        source = rng.randrange(len(parent));destination = rng.randrange(len(parent)-1)
        if destination >= source:
            destination += 1
        trial = list(parent);trial.insert(destination, trial.pop(source))
        pairs.append([source, destination]);orders.append(trial)
    return dict(seed_material=material, seed=str(seed), rng='Python random.Random',
        pair_policy='32 draws with replacement; preserve duplicates and nonpositive scores',
        pairs=pairs, orders=orders)


def comparison(actual, claimed):
    require(type(actual) in (int, float) and type(claimed) in (int, float)
            and math.isfinite(actual) and math.isfinite(claimed) and min(actual, claimed) >= 0, 'Finite nonnegative score parity')
    difference = abs(actual-claimed);tolerance = ATOL+RTOL*max(actual, claimed)
    within = difference <= tolerance;floor = max(actual, claimed) <= ATOL
    return dict(verified=within and not floor, within_tolerance=within, numerical_floor=floor,
        status='mismatch' if not within else 'floor_censored' if floor else 'verified',
        absolute_difference=difference, tolerance=tolerance)


def fit_temperature(probes, *, expected_count=2*PROBES):
    """Every fixed probe must be present/verified; only positive native deltas fit."""
    require(len(probes) == expected_count and all(row.get('verified') is True for row in probes),
            'Complete verified fixed calibration probes required')
    if expected_count == 2*PROBES:
        expected = {(spec['case'], index) for spec in CASES for index in range(1, PROBES+1)}
        observed = [(row.get('case'), row.get('index')) for row in probes]
        require(set(observed) == expected and len(set(observed)) == len(observed)
                and all(type(row.get('index')) is int for row in probes), 'Exact fixed probe roster required')
        for spec in CASES:
            require(len({row['parent_J'] for row in probes if row['case'] == spec['case']}) == 1,
                    'Calibration case changed its parent normalization')
    positive = []
    for row in probes:
        parent, trial = row['parent_J'], row['native_J']
        require(type(parent) in (int, float) and type(trial) in (int, float) and
                math.isfinite(parent) and math.isfinite(trial) and min(parent, trial) >= 0, 'Finite normalized calibration energies')
        delta = trial/max(parent, ENERGY_FLOOR)-parent/max(parent, ENERGY_FLOOR)
        require(math.isfinite(delta), 'Normalized calibration delta is not finite')
        if delta > 0:
            positive.append(delta)
    require(len(positive) >= MIN_POSITIVE, 'Insufficient positive verified normalized deltas')
    d = statistics.median(positive)
    start, end = -d/math.log(P_START), -d/math.log(P_END)
    require(math.isfinite(start) and math.isfinite(end) and start > end > 0, 'Temperature conversion outside resolution')
    return dict(positive_delta_count=len(positive), pooled_positive_median=d,
        theta_start=start, theta_end=end, initial_typical_uphill_acceptance=P_START,
        final_typical_uphill_acceptance=P_END, fit_scope='native training J only; no physical errors')


def source_hashes(check):
    result = {}
    for name in SOURCE_NAMES:
        check();result[name] = digest(bounded_read(ROOT/name, 1024**2))
    return result


def _request(path, model, spec, plan, native_rows, sources, a):
    fields = dict(version=1, case=spec['case'], model_sha256=spec['sha256'],
        legacy_source_model_sha256=spec['legacy_source_model_sha256'], n=model['n'],
        reference=model['reference'], words=model['words'], units=len(model['blocks']), time=1.,
        budget=30*model['words'], probe_count=PROBES, seed_material=plan['seed_material'],
        seed=plan['seed'], pairs=plan['pairs'], orders=plan['orders'])
    lines = [key+'='+a._toml(value) for key, value in fields.items()]
    lines.append('[source_sha256]')
    lines.extend(a._toml(key)+'='+a._toml(value) for key, value in sorted(sources.items()))
    for row in native_rows:
        lines.append('[[native_scores]]')
        lines.extend(key+'='+a._toml(value) for key, value in row['score'].items())
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= 8*1024**2 and not path.exists(), 'Fresh bounded calibration request')
    with path.open('xb') as stream:
        stream.write(raw)
    return digest(raw)


def run_child(executable, model, request, output, cpu_cap, wall_cap, check, a):
    require(not output.exists(), 'Fresh Julia calibration directory')
    argv = [executable, '--startup-file=no', '--history-file=no', '--threads=1', '--gcthreads=1',
        '--heap-size-hint=32M', '--project='+str(ROOT/'analysis/julia_window'), str(HERE/'calibration.jl'),
        str(model), str(request), str(output), repr(cpu_cap), repr(wall_cap)]
    env = os.environ.copy()
    env.update(JULIA_NUM_THREADS='1', JULIA_NUM_GC_THREADS='1', OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', JULIA_LOAD_PATH='@:@stdlib',
        JULIA_PKG_PRECOMPILE_AUTO='0', JULIA_IMAGE_THREADS='1')
    began = _time.monotonic();stats = dict(command=argv, reaped=False, raw_CPU_seconds=None)
    stderr_path=output.parent/(output.name+'.stderr.log')
    stderr_path.parent.mkdir(parents=True,exist_ok=True)
    require(not stderr_path.exists(),'Fresh Julia calibration diagnostic log')
    with stderr_path.open('xb') as diagnostic:
        process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=diagnostic, env=env, start_new_session=True)
    stats['stderr_path']=str(stderr_path)
    reaped, stop = False, None
    try:
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                reaped = True;break
            live = a._proc_cpu(process.pid)
            if live is None:
                stop = 'owned_child_cpu_unavailable'
            elif live >= cpu_cap:
                stop = 'child_cpu_limit'
            elif _time.monotonic()-began >= wall_cap:
                stop = 'child_wall_limit'
            else:
                try:
                    check(live)
                except Refusal as exc:
                    stop = str(exc)
            if stop:
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
        stats.update(reaped=True, returncode=process.returncode, raw_CPU_seconds=usage.ru_utime+usage.ru_stime,
            peak_RSS_KiB=usage.ru_maxrss, wall_seconds=_time.monotonic()-began,
            stop=stop or ('complete' if process.returncode == 0 else 'child_refusal'))
    return stats


def _admit_child_row(path, case, index, order, native_row, model_sha, request_sha, sources):
    raw = bounded_read(path, 64*1024);row = tomllib.loads(raw.decode('ascii'))
    legacy = next(spec['legacy_source_model_sha256'] for spec in CASES if spec['case'] == case)
    require(type(row.get('version')) is int and row['version'] == 1 and row['case'] == case and
            row.get('legacy_source_model_sha256') == legacy and type(row['index']) is int and row['index'] == index and
            row['order'] == order and all(type(x) is int for x in row['order']) and
            row['model_sha256'] == model_sha and row['request_sha256'] == request_sha and
            row['source_sha256'] == sources, 'Julia calibration row source/input/order join')
    require(row['native_claim'] == native_row['score'], 'Julia row changed native claim')
    require(all(type(row['native_claim'][key]) in (int, float) and math.isfinite(row['native_claim'][key])
            and row['native_claim'][key] >= 0 for key in ('value', 'cubic_response')), 'Finite numeric native claim')
    require(type(row['native_claim']['steps']) is int and row['native_claim']['steps'] == native_row['score']['steps'], 'Native step identity')
    require(isinstance(row.get('julia_version'), str) and 1 <= len(row['julia_version']) <= 64,
            'Julia runtime identity required')
    verified = row.get('status') == 'complete'
    comparisons = {}
    for name in ('full_order_score', 'archived_score'):
        claim = row.get(name)
        if not isinstance(claim, dict):
            verified = False;continue
        require(set(claim) == {'value', 'cubic_response', 'steps'} and type(claim['steps']) is int and claim['steps'] > 0,
                'Complete Julia score witness')
        parity = {key: comparison(claim[key], native_row['score'][key]) for key in ('value', 'cubic_response')}
        parity['steps_match'] = claim['steps'] == native_row['score']['steps']
        parity['normalized'] = abs(claim['value']-claim['cubic_response']/claim['steps']**4) <= max(1e-30, abs(claim['value'])*1e-12)
        verified = verified and all(parity[key]['verified'] for key in ('value', 'cubic_response')) and parity['steps_match'] and parity['normalized']
        comparisons[name] = parity
    if isinstance(row.get('full_order_score'), dict) and isinstance(row.get('archived_score'), dict):
        comparisons['full_vs_archived'] = {key: comparison(row['full_order_score'][key], row['archived_score'][key])
                                             for key in ('value', 'cubic_response')}
        verified = verified and all(value['verified'] for value in comparisons['full_vs_archived'].values())
    result = dict(row, verified=bool(verified), independent_python_parity=comparisons,
                  artifact_sha256=digest(raw))
    return result


def retain_probe_inventory(report, declarations):
    """Timeouts retain every declared identity and every completed native score."""
    present={(row['case'],row['index']) for row in report['probes']}
    summaries={row['case']:row for row in report['cases']}
    for declared in declarations:
        case=declared['case']
        native={row['index']:row for row in summaries.get(case,{}).get('native_rows',[])}
        parent=native.get(0,{}).get('score',{}).get('value')
        for index in range(1,PROBES+1):
            if (case,index) in present:
                continue
            trial=native.get(index,{}).get('score',{}).get('value')
            normalized=None if parent is None or trial is None else trial/max(parent,ENERGY_FLOOR)-parent/max(parent,ENERGY_FLOOR)
            report['probes'].append(dict(case=case,index=index,pair=declared['pairs'][index-1],
                order=declared['orders'][index],model_sha256=declared['sha256'],
                legacy_source_model_sha256=declared['legacy_source_model_sha256'],
                parent_J=parent,native_J=trial,normalized_delta=normalized,verified=False,
                status='unprocessed' if trial is None else 'native_scored_Julia_evidence_incomplete',
                positive=None if normalized is None else normalized>0,
                objective_at_or_below_floor=None if trial is None else trial<=ATOL,
                energy_at_or_below_floor=None if trial is None else trial<=ENERGY_FLOOR))
    case_order={spec['case']:i for i,spec in enumerate(CASES)}
    report['probes'].sort(key=lambda row:(case_order[row['case']],row['index']))


def run_calibration(output, cpu_cap, wall_cap, julia_executable=None):
    require(type(cpu_cap) in (int, float) and math.isfinite(cpu_cap) and 1 <= cpu_cap <= 3600 and
            type(wall_cap) in (int, float) and math.isfinite(wall_cap) and 5 <= wall_cap <= 18000,
            'Finite supplied calibration caps required')
    output = Path(output).resolve()
    require(not output.exists(), 'Fresh calibration output directory required')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='molecular-generalization-old-calibration-v1', complete=False, status='initializing',
        physical_errors_read=False, new_chemistry_generated=False, failures=[], cases=[], probes=[],
        Julia_jobs=[], limits=dict(CPU_seconds=cpu_cap, wall_seconds=wall_cap),
        objective=dict(time=1., budget_multiplier=30, energy_floor=ENERGY_FLOOR, rtol=RTOL, atol=ATOL),
        python=dict(executable=str(Path(sys.executable).resolve()), version=sys.version),
        external_OS_receipt_required_for_admission=True,
        accounting_scope='absolute CLI Python process CPU plus exact owned wait4 child CPU; external receipt includes final output/teardown')
    child_CPU = 0.;token = native = context = None;declarations=[]

    def check(live=0.):
        require(_time.process_time()+child_CPU+live < cpu_cap-FINAL_CPU_RESERVE, 'calibration_cpu_limit')
        require(_time.monotonic()-BOOT_WALL < wall_cap-FINAL_WALL_RESERVE, 'calibration_wall_limit')

    try:
        check();a,native = adapter();sources = source_hashes(check)
        prepared = []
        for spec in CASES:
            check();model = read_model(ROOT/spec['path'], spec);plan = probe_plan(spec['case'], model['parent'])
            target = output/'inputs'/(spec['case']+'.toml');target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as stream:
                stream.write(model['raw'])
            declaration = dict(spec, copied_input=str(target.relative_to(output)),
                legacy_header_field='source_model_sha256', seed_material=plan['seed_material'], seed=plan['seed'],
                rng=plan['rng'], pair_policy=plan['pair_policy'], pairs=plan['pairs'], orders=plan['orders'])
            declarations.append(declaration);prepared.append((spec,model,plan,target))
        report['declaration_sha256']=write_new(output/'DECLARATION.json',dict(schema=report['schema'], case_order=[s['case'] for s in CASES],
            cases=declarations, source_sha256=sources, physical_error_read=False,
            probe_count_per_case=PROBES, minimum_positive_deltas=MIN_POSITIVE,
            fit='pooled median of positive native normalized deltas after all rows verify',
            initial_acceptance=P_START, final_acceptance=P_END))
        report['source_sha256']=sources
        executable = shutil.which(julia_executable or os.environ.get('QHAT_JULIA_EXECUTABLE', 'julia'))
        require(executable is not None, 'Julia executable unavailable')
        require(hasattr(os, 'wait4') and Path('/proc/self/stat').is_file() and
                signal.getsignal(signal.SIGCHLD) == signal.SIG_DFL, 'Owned Linux child accounting unavailable')
        report['julia_executable']=str(Path(executable).resolve())
        token = native._EXTERNAL_BUDGET_CHECK.set(check)
        for spec,model,plan,target in prepared:
            case=spec['case'];summary=dict(case=case, model_sha256=spec['sha256'], native_rows=[], Julia_rows=[], complete=False)
            report['cases'].append(summary)
            native._compile(model['blocks'], model['n']) # Admission retains complete original partition.
            context = native.SparseResponse(model['blocks'], model['n'], model['reference'], time=1.,
                max_entries=50000, max_cached_actions=100000)
            native_rows=[]
            for index,order in enumerate(plan['orders']):
                check();score=a._checked_score(context.score(tuple(order),30*model['words']))
                row=dict(case=case,index=index,order=order,score=dict(value=score.value,
                    cubic_response=score.cubic_response,steps=score.steps), model_sha256=spec['sha256'],
                    legacy_source_model_sha256=spec['legacy_source_model_sha256'],
                    objective_at_or_below_floor=score.value<=ATOL, energy_at_or_below_floor=score.value<=ENERGY_FLOOR)
                native_rows.append(row);summary['native_rows'].append(row)
                write_new(output/'native'/case/('row'+str(index).zfill(3)+'.json'),row)
            context=None
            request=output/(case+'.REQUEST.toml')
            request_sha=_request(request,model,spec,plan,native_rows,sources,a)
            check();remaining=cpu_cap-_time.process_time()-child_CPU-FINAL_CPU_RESERVE
            remaining_wall=wall_cap-(_time.monotonic()-BOOT_WALL)-FINAL_WALL_RESERVE
            require(remaining>.2 and remaining_wall>1., 'Insufficient remaining calibration allowance')
            job_record=dict(case=case,raw_CPU_seconds=None,stop='dispatched',request_sha256=request_sha)
            report['Julia_jobs'].append(job_record)
            job=run_child(executable,target,request,output/'julia'/case,remaining,remaining_wall,check,a)
            job_record.update(job);child_CPU+=job['raw_CPU_seconds']
            all_rows=[]
            for index,order in enumerate(plan['orders']):
                path=output/'julia'/case/('row'+str(index).zfill(3)+'.toml')
                if not path.is_file():
                    summary['Julia_rows'].append(dict(index=index,order=order,verified=False,status='missing'))
                    continue
                row=_admit_child_row(path,case,index,order,native_rows[index],spec['sha256'],request_sha,sources)
                all_rows.append(row);summary['Julia_rows'].append(row)
            summary_path=output/'julia'/case/'SUMMARY.toml'
            child_summary={}
            if summary_path.is_file():
                summary_bytes=bounded_read(summary_path,64*1024)
                child_summary=tomllib.loads(summary_bytes.decode('ascii'))
                require(child_summary['case']==case and child_summary['model_sha256']==spec['sha256'] and
                    child_summary['request_sha256']==request_sha and child_summary['source_sha256']==sources,
                    'Julia calibration summary identity')
                summary['Julia_summary_sha256']=digest(summary_bytes)
            summary['complete']=len(all_rows)==PROBES+1 and all(row['verified'] for row in all_rows) and job['returncode']==0 and job['stop']=='complete' and child_summary.get('complete') is True and type(child_summary.get('processed_rows')) is int and child_summary['processed_rows']==PROBES+1
            parent_J=native_rows[0]['score']['value']
            byindex={row['index']:row for row in all_rows}
            for index in range(1,PROBES+1):
                trial_J=native_rows[index]['score']['value'];scale=max(parent_J,ENERGY_FLOOR)
                report['probes'].append(dict(case=case,index=index,pair=plan['pairs'][index-1],order=plan['orders'][index],
                    model_sha256=spec['sha256'],legacy_source_model_sha256=spec['legacy_source_model_sha256'],
                    parent_J=parent_J,native_J=trial_J,normalized_parent=parent_J/scale,normalized_trial=trial_J/scale,
                    normalized_delta=trial_J/scale-parent_J/scale,verified=summary['complete'] and byindex.get(index,{}).get('verified',False),
                    positive=trial_J/scale-parent_J/scale>0, objective_at_or_below_floor=trial_J<=ATOL,
                    energy_at_or_below_floor=trial_J<=ENERGY_FLOOR))
            if not summary['complete']:
                report['failures'].append(dict(type='CaseRefusal',message='Incomplete/refused/floor-censored or mismatched calibration case: '+case))
        require(len(report['cases'])==len(CASES) and all(row['complete'] for row in report['cases']), 'Both fixed old calibration cases must complete')
        report['fit']=fit_temperature(report['probes'])
        check();report.update(complete=True,status='complete')
    except Exception as exc:
        report['failures'].append(dict(type=type(exc).__name__,message=str(exc)[:2048]))
        report.update(complete=False,status='incomplete')
    finally:
        if native is not None and token is not None:
            native._EXTERNAL_BUDGET_CHECK.reset(token)
        context=None
        retain_probe_inventory(report,declarations)
        unknown=sum(job.get('raw_CPU_seconds') is None for job in report['Julia_jobs'])
        report['observed_cost_before_summary']=dict(Python_process_CPU_seconds=_time.process_time(),
            owned_Julia_CPU_seconds=None if unknown else child_CPU,
            known_owned_Julia_CPU_seconds=child_CPU,unknown_terminal_Julia_costs=unknown,
            total_CPU_seconds=None if unknown else _time.process_time()+child_CPU,
            wall_seconds=_time.monotonic()-BOOT_WALL)
    report_sha=write_new(output/'CALIBRATION.json',report)
    if report['complete']:
        frozen=dict(a.DEFAULT_SETTINGS,calibrated=True,calibration_sha256=report_sha,
                    theta_start=report['fit']['theta_start'],theta_end=report['fit']['theta_end'])
        a.validate_settings(frozen)
        write_new(output/'ANNEALING_SETTINGS.json',frozen)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cpu-cap',type=float,required=True)
    parser.add_argument('--wall-cap',type=float,required=True)
    parser.add_argument('--julia-executable')
    args=parser.parse_args(argv)
    report=run_calibration(args.output,args.cpu_cap,args.wall_cap,args.julia_executable)
    print('CALIBRATION',report['status'],args.output,flush=True)
    return 0 if report['complete'] else 2


if __name__=='__main__':
    raise SystemExit(main())
