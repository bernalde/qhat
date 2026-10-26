"""Pinned source adapters and immutable practical-study constants."""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = HERE.parent
CASES = ('Li2_r2p95', 'F2_r1p55', 'H2O_scale1p10')
TIME = 2.11
GRID = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 40, 48, 64, 80, 96, 128)
TARGETS = (1e-4, 1e-6, 1e-8)
BUDGETS = (2., 8.)
REPETITIONS = 3
METHODS = ('first_JD', 'pauli_JD', 'cubic_JD', 'cubic_warm_JD')
REFERENCE_ARMS = ('current', 'pauli_JD_large_sameT', 'cubic_insertion')
ORDERING_MERGE = '8fd424f532cce04b595abc3445cc032ab03ab6c4'
AMBIGUITY = 2e-12
PHASE_LIMITS = {
    'tests': dict(CPU_seconds=120, wall_seconds=900, disk_MiB=64),
    'preflight': dict(CPU_seconds=600, wall_seconds=3600, disk_MiB=256),
    'measurement': dict(CPU_seconds=3600, wall_seconds=21600, disk_MiB=512),
    'independent': dict(CPU_seconds=3600, wall_seconds=21600, disk_MiB=512),
    'reporting': dict(CPU_seconds=300, wall_seconds=1800, disk_MiB=512),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(Path(module.__file__).resolve() == path.resolve(), 'Source import differs')
    return module


# Older scripts import their sibling by its original name. Authenticate that
# binding instead of silently accepting a same-named installed module.
_prior_path = list(sys.path)
sys.path.insert(0, str(AUDIT / 'compiled_accuracy'))
compiler = __import__('compiler')
require(Path(compiler.__file__).resolve() == AUDIT / 'compiled_accuracy/compiler.py',
        'Wrong compiler binding')
compiled = load('practical_prior_compiled_study', AUDIT / 'compiled_accuracy/study.py')
cubic = load('practical_prior_cubic', AUDIT / 'symmetric_response/cubic.py')
reference = load('practical_literal_reference', AUDIT / 'compiled_accuracy/reference.py')
sys.path.insert(0, str(AUDIT / 'prospective_reference'))
for _name in ('engine', 'one_pass'):
    _module = __import__(_name)
    require(Path(_module.__file__).resolve() == AUDIT / 'prospective_reference' / (_name + '.py'),
            'Wrong placement helper binding')
placement = load('practical_legacy_placement', AUDIT / 'prospective_reference/campaign.py')
objective_reference = load('practical_objective_reference', AUDIT / 'symmetric_response/independent.py')
sys.path[:] = _prior_path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def safe(name):
    path = PurePosixPath(name)
    require(not path.is_absolute() and '..' not in path.parts and str(path) == name,
            'Unsafe artifact name')
    return name


def source_paths():
    own = [p for p in HERE.iterdir() if p.suffix in ('.py', '.md')]
    dependencies = [
        AUDIT / 'compiled_accuracy' / name for name in
        ('compiler.py', 'study.py', 'reference.py')]
    dependencies += [AUDIT / 'symmetric_response' / name for name in ('cubic.py', 'independent.py')]
    dependencies += [AUDIT / 'prospective_reference' / name for name in
                     ('engine.py', 'one_pass.py', 'campaign.py')]
    return sorted(own + dependencies)


def sources():
    return {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}


def environment():
    return dict(python=sys.version, executable=str(Path(sys.executable).resolve()),
                packages={name: importlib.metadata.version(name) for name in ('numpy', 'scipy')})


def commit_sources():
    commit = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for path in source_paths():
        relative = str(path.relative_to(ROOT))
        blob = subprocess.check_output(['git', '-C', str(ROOT), 'show', commit + ':' + relative])
        require(blob == path.read_bytes(), 'Commit source before execution: ' + relative)
    return commit
