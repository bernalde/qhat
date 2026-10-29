"""Bounded serialization, source identity and retained-control loading.

Optimization and propagation live in Julia or the existing named controls.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.machinery
import importlib.metadata
import importlib.util
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = '24b3bf5388dfd19196e4ef210ba73fa9c574a11e'
TRAIN_TIME = 1.0
TRAIN_CAP = 30
TRAIN_CPU_SECONDS = 64.0
TRAIN_WALL_SECONDS = 320.0
TIMES = (0.4, 1.6)
CAPS = (8, 24)
REPETITIONS = 2
ARMS = ('canonical', 'magnitude', 'current', 'first', 'pauli', 'cubic',
        'window_sparse', 'window_magnitude')
FORMULAS = {arm: ((1, 2) if arm in ('canonical', 'current', 'first') else (2,))
            for arm in ARMS}
MAX_TERMS = 8192
MAX_UNITS = 1024
MAX_MODEL_BYTES = 8 * 1024**2


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    path = Path(path)
    require(path.is_file() and path.stat().st_size <= 32 * 1024**2, 'Bounded JSON file')
    return json.loads(path.read_text())


def write_new(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.partial')
    require(not path.exists() and not temporary.exists(), 'Refuse existing result or partial output')
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    os.link(temporary, path)  # Atomic EEXIST rather than replacing a raced output.
    temporary.unlink()


def decode(raw):
    require(isinstance(raw, list) and len(raw) <= MAX_TERMS, 'Pauli term count')
    result = {}
    for word, value in raw:
        require(type(value) in (int, str), 'Integer coefficient spelling')
        if isinstance(value, str):
            require(len(value) <= 1280 and value.lstrip('+-').isdigit(), 'Bounded integer spelling')
        value = int(value)
        require(value.bit_length() <= 4096 and value != 0, 'Finite nonzero bounded numerator')
        require(isinstance(word, list) and len(word) <= 20, 'Word length')
        support = []
        for site, letter in word:
            require(type(site) is int and 0 <= site < 20 and letter in ('X', 'Y', 'Z'), 'Pauli site')
            support.append((site, letter))
        key = tuple(support)
        require(len({site for site, _ in key}) == len(key), 'Repeated word site')
        require(key == tuple(sorted(key)) and key not in result, 'Canonical unique word')
        result[key] = value
    return result


def blocks(model, *, pauli=False):
    require(type(model['n']) is int and model['n'] == 14 and
            type(model['electrons']) is int and model['electrons'] == 4,
            'Declared fourteen-mode/four-electron instance')
    require(0 < len(model['coefficients']) <= MAX_TERMS and
            0 < len(model['blocks']) <= MAX_UNITS, 'Original term/unit bound')
    denominator = int(model['denominator'])
    require(denominator > 0 and denominator.bit_length() <= 4096, 'Bounded positive denominator')
    total = decode(model['coefficients'])
    units = [decode(raw) for raw in model['blocks']]
    require(all(units) and sum(map(len, units)) == len(total) and
            {word: value for unit in units for word, value in unit.items()} == total,
            'Original units exactly partition the supplied Hamiltonian')
    require(all(site < model['n'] for word in total for site, _ in word), 'Word fits register')
    if pauli:
        units = [{word: value} for word, value in sorted(total.items())]
    return tuple({word: Fraction(value, denominator) for word, value in unit.items()} for unit in units)


def bootstrap():
    if 'qhat' not in sys.modules:
        spec = importlib.machinery.ModuleSpec('qhat', loader=None, is_package=True)
        spec.submodule_search_locations = [str(ROOT)]
        sys.modules['qhat'] = importlib.util.module_from_spec(spec)
    from qhat.analysis import sparse_response_ordering as ordering
    require(Path(ordering.__file__).resolve() == ROOT / 'analysis/sparse_response_ordering.py',
            'Wrong normal optimizer source')
    return ordering


def write_model_toml(path, model, *, source_json_sha256, order=None, pauli=False):
    """Serialize exact supplied coefficients, without a float round trip."""
    original = blocks(model, pauli=pauli)
    raw = [{word: int(coefficient * int(model['denominator']))
            for word, coefficient in unit.items()} for unit in original]
    order = list(range(len(raw))) if order is None else list(order)
    require(all(type(index) is int for index in order) and
            sorted(order) == list(range(len(raw))), 'Complete original-label order')
    quote = json.dumps
    lines = ['case = ' + quote(model['case']), 'n = ' + str(model['n']),
             'electrons = ' + str(model['electrons']),
             'denominator = ' + quote(str(model['denominator'])),
             'partition = ' + quote('pauli' if pauli else 'conserving'),
             'source_json_sha256 = ' + quote(source_json_sha256),
             'baseline_order = ' + repr([index + 1 for index in order])]
    for unit in raw:
        ordered = sorted(unit.items())
        lines.extend(['', '[[units]]', 'words = ' + json.dumps(
            [','.join(str(site) + letter for site, letter in word) for word, _ in ordered]),
            'numerators = ' + json.dumps([str(value) for _, value in ordered])])
    encoded = ('\n'.join(lines) + '\n').encode()
    require(len(encoded) <= MAX_MODEL_BYTES, 'TOML model byte bound')
    with Path(path).open('xb') as stream:
        stream.write(encoded)


def retained_training():
    """Load unchanged control source with its own unchanged common module.

    The scoped import alias is restored immediately; no numerical constant or
    source is monkeypatched. The old twelve-mode _construct wrapper is unused.
    Its current/first_search routines retain their own admission limits.
    """
    common_path = ROOT / 'numerical_audit/scalable_ordering/common.py'
    training_path = common_path.with_name('training.py')
    previous = sys.modules.get('common')
    try:
        specification = importlib.util.spec_from_file_location('common', common_path)
        common = importlib.util.module_from_spec(specification)
        sys.modules['common'] = common
        specification.loader.exec_module(common)
        specification = importlib.util.spec_from_file_location('_window_retained_training', training_path)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
    finally:
        if previous is None:
            sys.modules.pop('common', None)
        else:
            sys.modules['common'] = previous
    require(module.k is common and common.ROOT == ROOT, 'Retained control import identity')
    return module


def source_paths():
    own = [path for path in HERE.rglob('*') if path.is_file() and
           'evidence' not in path.relative_to(HERE).parts and
           (path.suffix in ('.py', '.jl', '.json', '.toml') or
            path.name in ('PROTOCOL.md', 'README.md'))]
    dependencies = [ROOT / name for name in (
        'analysis/sparse_response_ordering.py', 'analysis/strang_schedule.py',
        'analysis/train_ordering.py', 'analysis/window_ordering.py',
        'analysis/fermionic_coefficients.py', 'analysis/__init__.py',
        'numerical_audit/__init__.py', 'numerical_audit/followup_classes/exact_blocks.py',
        'numerical_audit/scalable_ordering/common.py', 'numerical_audit/scalable_ordering/training.py',
        'numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl')]
    assets = ROOT / 'analysis/julia_window'
    if assets.exists():
        dependencies.extend(path for path in assets.rglob('*') if path.is_file() and
                            path.suffix in ('.jl', '.toml', '.json', '.md'))
    return sorted(set(own + dependencies))


def source_map():
    return {str(path.relative_to(ROOT)): sha(path) for path in source_paths()}


def committed_sources():
    head = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for name, digest in source_map().items():
        recorded = subprocess.check_output(['git', '-C', str(ROOT), 'show', head + ':' + name])
        require(hashlib.sha256(recorded).hexdigest() == digest, 'Commit scientific source before execution: ' + name)
    return head


def environment():
    basis_relative = 'pyscf/gto/basis/pople-basis/6-31G.dat'
    basis = Path(importlib.metadata.distribution('pyscf').locate_file(basis_relative))
    require(basis.is_file(), 'Pinned basis file missing')
    return dict(python=sys.version, executable=str(Path(sys.executable).resolve()), platform=platform.platform(),
                basis=dict(path=basis_relative, sha256=sha(basis), bytes=basis.stat().st_size),
                packages={name: importlib.metadata.version(name) for name in
                          ('numpy', 'scipy', 'openfermion', 'openfermionpyscf', 'pyscf', 'h5py')})
