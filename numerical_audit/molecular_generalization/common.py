"""Scientific identities and exact input serialization for this new study.

Producer snapshots are immutable. Later consumers record their own source map
and verify the producer snapshot, rather than pretending to be that producer.
"""
from __future__ import annotations

from contextlib import contextmanager, ExitStack
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARMS = ('full', 'window', 'annealing')
CHEAP = ('canonical', 'magnitude')
BUDGETS = (32, 64)
TIMES = (.65, 1.35)
CAPS = (8, 24)
MAX_TERMS, MAX_UNITS, MAX_MODEL_BYTES = 8192, 1024, 8*1024**2
SCORE_KIND = 'diagonal_reference_cubic_response'
MAX_CALIBRATION_FILES = 4096
MAX_CALIBRATION_FILE_BYTES = 32*1024**2
MAX_CALIBRATION_TOTAL_BYTES = 128*1024**2
MAX_CALIBRATION_MANIFEST_BYTES = 1024**2


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


prior = load_module('_molecular_prior_common', ROOT/'numerical_audit/prospective_window_ordering/study_common.py')
retained = load_module('_molecular_retained_common', ROOT/'numerical_audit/cost_matched_continuation/common.py')
require, sha, read, write_new = retained.require, retained.sha, retained.read, retained.write_new
decode, bootstrap = prior.decode, prior.bootstrap


def configuration():
    value = read(HERE/'PROTOCOL.json')
    require(value['schema'] == 'molecular-generalization-v1', 'Protocol schema')
    require(value['total_CPU_budgets'] == list(BUDGETS) and value['repetitions'] == 2,
            'Frozen CPU/repetition grid')
    require(value['optimized_arms'] == list(ARMS) and value['cheap_arms'] == list(CHEAP), 'Frozen arms')
    require(value['training_time'] == 1. and value['training_rotation_multiplier'] == 30,
            'Common cubic objective T1/30U')
    require(value['physical']['times'] == list(TIMES) and
            value['physical']['rotation_cap_multipliers'] == list(CAPS), 'Physical grid')
    old = retained.configuration()
    for group, keys in dict(parent=('CPU_seconds', 'wall_seconds', 'max_evaluations', 'max_moves'),
            native=('max_evaluations', 'max_moves', 'max_entries', 'max_cached_actions'),
            window=('width', 'max_windows', 'max_sweeps', 'parent_method', 'parent_fraction',
                    'max_entries', 'max_cached_actions', 'max_factor_entries', 'max_flat_entries',
                    'max_actions', 'max_atoms'),
            accounting=('CPU_cleanup_reserve', 'wall_cleanup_reserve', 'total_wall_multiplier')).items():
        require(all(value[group][key] == old[group][key] for key in keys), 'Retained control changed: '+group)
    require(value['native']['restarts']['count'] == 4 and
            value['native']['restarts']['insertions'] == [0, 1, 2, 4], 'Same four frozen restarts')
    return value


def registry():
    return read(HERE/'REGISTRY.json')


def validate_order(order, size):
    require(isinstance(order, (tuple, list)) and all(type(x) is int for x in order)
            and sorted(order) == list(range(size)), 'Complete original-label permutation')
    return list(order)


def blocks(model):
    """Admit dynamic active spaces without changing the archived 14m/4e loader."""
    n, e = model['n'], model['electrons']
    require(type(n) is int and 1 <= n <= 14 and type(e) is int and 0 < e < n
            and math.comb(n, e) <= 2000, 'Bounded fixed-N/full-register instance')
    require(0 < len(model['coefficients']) <= MAX_TERMS and 0 < len(model['blocks']) <= MAX_UNITS,
            'Original terms/units bound')
    denominator = int(model['denominator'])
    require(0 < denominator and denominator.bit_length() <= 4096, 'Exact denominator bound')
    total = decode(model['coefficients'])
    units = [decode(raw) for raw in model['blocks']]
    require(all(units) and sum(map(len, units)) == len(total) and
            {w: c for b in units for w, c in b.items()} == total, 'Exact original-target partition')
    require(all(site < n for word in total for site, _ in word), 'Word fits register')
    # Chemistry generation checks exact commutators once. The independent
    # consumer repeats them from authenticated input bytes. Constructors only
    # deserialize the already frozen partition; they do not import chemistry.
    from fractions import Fraction
    return tuple({word: Fraction(value, denominator) for word, value in unit.items()} for unit in units)


def write_model_toml(path, model, *, source_json_sha256):
    units = blocks(model)
    lines = dict(case=model['case'], n=model['n'], electrons=model['electrons'],
                 denominator=str(model['denominator']), partition='conserving',
                 source_json_sha256=source_json_sha256, baseline_order=list(range(1, len(units)+1)))
    records = []
    for unit in units:
        ordered = sorted(unit.items())
        records.append(dict(words=[','.join(str(q)+s for q, s in w) for w, _ in ordered],
            numerators=[str(int(c*int(model['denominator']))) for _, c in ordered]))
    write_toml(path, lines, {'units': records})


def write_toml(path, fields, tables):
    # Existing bounded scalar/array encoder; no new serialization dependency.
    bootstrap()
    from qhat.analysis.window_ordering import _toml
    lines = [key+' = '+_toml(value) for key, value in fields.items()]
    for name, rows in tables.items():
        for row in rows:
            lines.append('\n[['+name+']]')
            lines.extend(key+' = '+_toml(value) for key, value in row.items() if value is not None)
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= MAX_MODEL_BYTES, 'TOML byte bound')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('xb') as stream:
        stream.write(raw)


def contained(root, relative):
    root = Path(root).resolve()
    require(isinstance(relative, str) and relative and not Path(relative).is_absolute()
            and all(p not in ('.', '..', '') for p in relative.split('/')),
            'Canonical relative artifact path')
    path = root/relative
    require(path.resolve().is_relative_to(root), 'Artifact escaped data root')
    require(str(path.relative_to(root)) == relative, 'Canonical path spelling')
    return path


def _calibration_relative(value):
    require(isinstance(value, str) and 0 < len(value) <= 4096 and '\\' not in value
            and all(ord(c) >= 32 and ord(c) != 127 for c in value)
            and not value.startswith('/') and len(value.split('/')) <= 64
            and all(p not in ('', '.', '..') for p in value.split('/')),
            'Canonical calibration evidence path')
    return value.split('/')


@contextmanager
def _calibration_parent(root, relative, *, create=False):
    """Walk only real directories beneath a trusted root, including at open time."""
    parts = _calibration_relative(relative)
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    descriptor = None
    try:
        descriptor = os.open(root, flags)
        for part in parts[:-1]:
            if create:
                try:
                    os.mkdir(part, dir_fd=descriptor)
                except FileExistsError:
                    pass
            child = os.open(part, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        yield descriptor, parts[-1]
    except OSError as exc:
        raise ValueError('Calibration evidence I/O refused: '+relative) from exc
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _calibration_file(root, relative, *, declared=None, copy_to=None, capture=False,
                      limit=MAX_CALIBRATION_FILE_BYTES):
    """Stream a bounded regular file; a copy is published only after its checks."""
    with ExitStack() as stack:
        parent, name = stack.enter_context(_calibration_parent(root, relative))
        descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=parent)
        stream = stack.enter_context(os.fdopen(descriptor, 'rb'))
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and 0 <= info.st_size <= limit,
                'Bounded regular calibration evidence: '+relative)
        require(declared is None or info.st_size == declared['bytes'],
                'Calibration evidence size differs: '+relative)
        output = None
        if copy_to is not None:
            dest, leaf = stack.enter_context(_calibration_parent(copy_to, relative, create=True))
            descriptor = os.open(leaf+'.partial', os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                 0o644, dir_fd=dest)
            output = stack.enter_context(os.fdopen(descriptor, 'wb'))
        digest, count, chunks = hashlib.sha256(), 0, []
        while True:
            chunk = stream.read(min(1024**2, limit-count+1))
            if not chunk:
                break
            count += len(chunk)
            require(count <= limit, 'Growing calibration evidence byte cap: '+relative)
            digest.update(chunk)
            if output is not None:
                output.write(chunk)
            if capture:
                chunks.append(chunk)
        actual = dict(bytes=count, sha256=digest.hexdigest())
        require(count == info.st_size and (declared is None or
                all(actual[key] == declared[key] for key in ('bytes', 'sha256'))),
                'Calibration evidence hash/size differs: '+relative)
        if output is not None:
            output.flush()
            os.link(leaf+'.partial', leaf, src_dir_fd=dest, dst_dir_fd=dest, follow_symlinks=False)
            os.unlink(leaf+'.partial', dir_fd=dest)
        return actual, b''.join(chunks) if capture else None


def _calibration_json(raw):
    def unique(items):
        value = {}
        for key, item in items:
            require(key not in value, 'Duplicate calibration JSON key: '+key)
            value[key] = item
        return value
    def reject(value):
        raise ValueError('Nonfinite calibration JSON constant: '+value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=reject)


def calibration_evidence(study_root, expected_sha256, *, copy_to=None):
    """Authenticate/copy only manifest members; historical locators are metadata.

    Paths in the returned binding are relative to the archived study directory.
    Copies preserve every original byte, including all provenance and locators.
    """
    def digest(value):
        return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)
    manifest_info, raw = _calibration_file(study_root, 'CALIBRATION_EVIDENCE.json',
        capture=True, limit=MAX_CALIBRATION_MANIFEST_BYTES, copy_to=copy_to)
    manifest = _calibration_json(raw)
    require(isinstance(manifest, dict) and manifest.get('schema') == 'molecular-calibration-evidence-v1'
            and isinstance(manifest.get('files'), dict) and
            0 < len(manifest['files']) <= MAX_CALIBRATION_FILES, 'Bounded calibration manifest schema')
    prefix = manifest.get('prefix')
    require(len(_calibration_relative(prefix)) >= 2 and prefix.startswith('evidence/'),
            'Contained calibration evidence prefix')
    accepted = manifest.get('accepted_report')
    _calibration_relative(accepted)
    files = manifest['files']
    require(digest(expected_sha256) and manifest.get('accepted_sha256') == expected_sha256
            and accepted in files, 'Accepted calibration report declaration')
    total = 0
    for name, entry in files.items():
        _calibration_relative(prefix+'/'+name)
        require(isinstance(entry, dict) and set(entry) == {'bytes', 'sha256', 'original_locator'}
                and type(entry['bytes']) is int and 0 <= entry['bytes'] <= MAX_CALIBRATION_FILE_BYTES
                and digest(entry['sha256']) and isinstance(entry['original_locator'], str)
                and 0 < len(entry['original_locator']) <= 16384, 'Calibration evidence entry: '+name)
        total += entry['bytes']
    require(total <= MAX_CALIBRATION_TOTAL_BYTES, 'Calibration evidence total byte cap')
    require(files[accepted]['sha256'] == expected_sha256, 'Accepted report manifest digest')
    accepted_raw = None
    for name, entry in sorted(files.items()):
        _, contents = _calibration_file(study_root, prefix+'/'+name, declared=entry,
            copy_to=copy_to, capture=name == accepted)
        if name == accepted:
            accepted_raw = contents
    report_info, report_raw = _calibration_file(study_root, 'CALIBRATION.json', copy_to=copy_to, capture=True)
    require(report_info['sha256'] == expected_sha256 and report_raw == accepted_raw,
            'Top calibration must equal the exact accepted evidence report')
    report = _calibration_json(report_raw)
    require(isinstance(report, dict) and report.get('complete') is True
            and report.get('physical_errors_read') is False, 'Complete old-only calibration report')
    return dict(schema='molecular-frozen-calibration-evidence-v1',
        manifest='CALIBRATION_EVIDENCE.json', manifest_sha256=manifest_info['sha256'],
        manifest_bytes=manifest_info['bytes'], prefix=prefix, accepted_report=accepted,
        accepted_sha256=expected_sha256, report='CALIBRATION.json',
        declared_files=len(files), declared_bytes=total)


def source_paths():
    own = []
    for directory, directories, files in os.walk(HERE, followlinks=False):
        directories[:] = [name for name in directories if name not in ('test', 'evidence', '__pycache__')]
        own.extend(p for name in files if (p := Path(directory)/name).is_file()
                   and p.suffix in ('.py', '.jl', '.json', '.toml'))
    # Retain authenticated numerical kernels and all their existing study dependencies.
    calibration_inputs = [ROOT/'numerical_audit/julia_cubic_feasibility/evidence/inputs/cubic'/name
                          for name in ('CO_r1p13.toml', 'CH4_r1p09.toml')]
    return sorted(set(own + retained.source_paths() + calibration_inputs))


def source_map():
    return {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}


def committed_sources():
    head = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for name, digest in source_map().items():
        value = subprocess.check_output(['git', '-C', str(ROOT), 'show', head+':'+name])
        require(hashlib.sha256(value).hexdigest() == digest, 'Commit producer before execution: '+name)
    return head


def freeze(data, *, mode):
    data = Path(data)
    require(mode in ('main', 'preflight') and not data.exists(), 'Fresh named execution mode')
    config = configuration()
    if mode == 'main':
        require(config['main_execution_allowed'] is True and config['protocol_state'] == 'frozen_main',
                'Main disabled until authenticated old-only preflight')
        calibration = HERE/'CALIBRATION.json'
        require(calibration.is_file() and read(calibration)['complete'] is True and
                config['annealing']['calibration_sha256'] == sha(calibration), 'Frozen old-only calibration')
        gate = read(HERE/'PREFLIGHT_ACCEPTED.json')
        require(gate['complete'] is True and gate['new_physical_errors_observed'] is False,
                'Complete old-only end-to-end gate required')
    evidence = calibration_evidence(HERE, config['annealing']['calibration_sha256'])
    head, sources = committed_sources(), source_map()
    study = HERE.relative_to(ROOT)
    calibration_sources = {str(study/'CALIBRATION_EVIDENCE.json'): evidence['manifest_sha256'],
                           str(study/'CALIBRATION.json'): evidence['accepted_sha256']}
    require(all(sources.get(name) == digest for name, digest in calibration_sources.items()),
            'Calibration evidence joins committed producer sources')
    data.mkdir(parents=True)
    for name in sources:
        if name in calibration_sources:
            continue
        dest = data/'source_snapshot'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, dest)
    snapshot = data/'source_snapshot'/study
    snapshot.mkdir(parents=True, exist_ok=True)
    require(calibration_evidence(HERE, evidence['accepted_sha256'], copy_to=snapshot) == evidence,
            'Calibration evidence changed during freeze')
    write_new(data/'INPUT_PROTOCOL_FROZEN.json', dict(schema='molecular-source-freeze-v1',
        mode=mode, source_commit=head, producer_root=str(ROOT), producer_data_root=str(data.resolve()),
        sources=sources, configuration=config, calibration_evidence=evidence,
        registry=registry(), environment=prior.environment(), runtime=read(HERE/'ENVIRONMENT.json')))


def frozen(data, *, producer=False):
    data = Path(data)
    value = read(data/'INPUT_PROTOCOL_FROZEN.json')
    for name, digest in value['sources'].items():
        require(sha(contained(data/'source_snapshot', name)) == digest, 'Producer snapshot changed: '+name)
    require(value['configuration'] == read(data/'source_snapshot'/HERE.relative_to(ROOT)/'PROTOCOL.json'),
            'Snapshot protocol/configuration join')
    if str(HERE.relative_to(ROOT)/'CALIBRATION_EVIDENCE.json') in value['sources'] or 'calibration_evidence' in value:
        require('calibration_evidence' in value and value['calibration_evidence'] == calibration_evidence(
            data/'source_snapshot'/HERE.relative_to(ROOT),
            value['configuration']['annealing']['calibration_sha256']),
            'Frozen calibration evidence binding changed')
    if producer:
        require(value['sources'] == source_map() and value['configuration'] == configuration()
                and value['registry'] == registry(), 'Producer differs from immutable freeze')
        require(value['environment'] == prior.environment(), 'Producer environment changed')
    if value['mode'] == 'main':
        require(value['configuration']['main_execution_allowed'] is True and
                value['configuration']['protocol_state'] == 'frozen_main', 'Main phase authorization missing')
    else:
        require(value['mode'] == 'preflight', 'Known campaign mode')
    return value


def header(data, *, producer=False):
    data = Path(data)
    value = frozen(data, producer=producer)
    mode, config = value['mode'], value['configuration']
    cases = ([config['preflight']['case']] if mode == 'preflight' else
             [spec['case'] for spec in value['registry']['cases']])
    result = dict(mode=mode, cases=cases, repetitions=1 if mode == 'preflight' else 2,
        protocol_sha256=sha(data/'source_snapshot'/HERE.relative_to(ROOT)/'PROTOCOL.json'),
        source_freeze_sha256=sha(data/'INPUT_PROTOCOL_FROZEN.json'),
        protocol_path='source_snapshot/numerical_audit/molecular_generalization/PROTOCOL.json',
        source_freeze_path='INPUT_PROTOCOL_FROZEN.json')
    return result, config


def identities(repetitions):
    return [(rep, arm, budget) for rep in range(1, repetitions+1)
            for arm, budget in [('parent', 0), ('canonical', 0), ('magnitude', 0)]] \
            + [(rep, arm, b) for rep in range(1, repetitions+1) for arm in ARMS for b in BUDGETS]


def constructor_path(data, case, rep, arm, budget):
    require(isinstance(case, str) and case.isascii() and case.replace('_', '').replace('-', '').replace('.', '').isalnum(),
            'Bounded case label')
    return Path(data)/'constructors'/f'{case}.rep{rep}.{arm}.B{budget}.json'


def finite_score(value):
    require(isinstance(value, dict) and set(value) == {'value', 'cubic_response', 'steps'}, 'Complete score fields')
    require(all(type(value[key]) in (int, float) and math.isfinite(value[key]) and value[key] >= 0
                for key in ('value', 'cubic_response')) and type(value['steps']) is int and value['steps'] > 0,
            'Finite nonnegative score and positive steps')
    require(abs(value['value']-value['cubic_response']/value['steps']**4)
            <= max(1e-30, value['value']*1e-12), 'Cubic r^-4 normalization')
    return value
