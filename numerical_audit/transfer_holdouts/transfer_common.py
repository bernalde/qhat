"""Source-bound helpers and the prospective transfer panel's fixed constants."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
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
TRAIN_TIME = 1.11
TIMES = (.59, 1.41, 2.73)
CAPS = (10, 30)
STATES = ('HF', 'ground', 'excited', 'coherent_real', 'coherent_imag')
ARMS = ('hf_current', 'hf_first_JD', 'hf_pauli_JD', 'hf_cubic_JD',
        'state_first_JD', 'state_first_symmetric', 'state_pauli_JD', 'state_cubic_JD')
FLOOR = 1e-12
CORRELATED_THRESHOLD = 1e-10
SLATER_THRESHOLD = 1e-12
DEGENERACY_LABEL = 1e-9


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(name, path):
    if name in sys.modules:
        module = sys.modules[name]
    else:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        except BaseException:
            del sys.modules[name]
            raise
    require(Path(module.__file__).resolve() == path.resolve(), 'Wrong imported source: ' + name)
    return module


response = load('transfer_original_response', AUDIT / 'state_holdouts/response.py')
engine = response.determinant


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
    require(not path.is_absolute() and '..' not in path.parts and str(path) == name, 'Unsafe artifact path')
    return name


def source_paths():
    own = [p for p in HERE.iterdir() if p.suffix in ('.py', '.md', '.json')]
    deps = [AUDIT / name for name in ('state_holdouts/response.py', 'prospective_reference/engine.py',
                                     'symmetric_response/cubic.py', 'prospective_reference/one_pass.py',
                                     'prospective_reference/campaign.py', 'compiled_accuracy/compiler.py')]
    return sorted(own + deps)


def sources():
    return {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}


def committed_sources():
    commit = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for relative in sources():
        blob = subprocess.check_output(['git', '-C', str(ROOT), 'show', commit + ':' + relative])
        require(blob == (ROOT / relative).read_bytes(), 'Commit scientific source before reference preparation')
    return commit


def environment():
    return dict(python=sys.version, executable=str(Path(sys.executable).resolve()),
                packages={name: importlib.metadata.version(name) for name in ('numpy', 'scipy')})
