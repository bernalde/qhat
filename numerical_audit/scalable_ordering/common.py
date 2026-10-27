"""Source-bound constants and helpers for the larger sparse-ordering study."""
import hashlib
import importlib.machinery
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TRAIN_TIME = 1.0
TIMES = (0.5, 1.5)
CAPS = (10, 30)
ARMS = ('canonical', 'magnitude', 'current', 'first', 'pauli', 'cubic')
TRAIN_CPU_SECONDS = 2.0
MAX_TERMS = 4096
MAX_UNITS = 1024
MAX_MODEL_BYTES = 8 * 1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    path = Path(path)
    require(path.stat().st_size <= 32 * 1024**2, 'JSON input byte bound')
    return json.loads(path.read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def load_vector(path, key, n):
    """Validate ZIP member and NPY header bounds before NumPy allocates a state."""
    import numpy as np
    path = Path(path)
    require(type(n) is int and 1 <= n <= 12 and key in ('actual', 'exact'), 'Bounded state identity')
    dimension = 1 << n
    require(path.stat().st_size <= 2 * 1024**2, 'Compressed state byte bound')
    with zipfile.ZipFile(path) as archive:
        members = archive.infolist()
        require(len(members) == 1 and members[0].filename == key + '.npy' and
                members[0].compress_type in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED) and
                not members[0].flag_bits & 1 and members[0].file_size <= 2 * 1024**2,
                'NPZ member identity/size bound')
        with archive.open(members[0]) as stream:
            require(np.lib.format.read_magic(stream) == (1, 0), 'Supported bounded NPY header version')
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream, max_header_size=4096)
            require(shape == (dimension,) and not fortran and not dtype.hasobject and
                    dtype.kind == 'c' and dtype.itemsize in (16, 32), 'State shape/dtype before allocation')
            expected = dimension * dtype.itemsize
            require(members[0].file_size == stream.tell() + expected, 'NPY payload byte count')
    with np.load(path, allow_pickle=False, max_header_size=4096) as payload:
        result = payload[key]
    require(result.shape == (dimension,) and np.all(np.isfinite(result)), 'Finite state vector')
    return result


def bootstrap():
    if 'qhat' not in sys.modules:
        spec = importlib.machinery.ModuleSpec('qhat', loader=None, is_package=True)
        spec.submodule_search_locations = [str(ROOT)]
        sys.modules['qhat'] = importlib.util.module_from_spec(spec)
    from qhat.analysis import sparse_response_ordering as ordering
    require(Path(ordering.__file__).resolve() == ROOT / 'analysis/sparse_response_ordering.py',
            'Wrong sparse ordering import')
    return ordering


def decode(raw):
    result = {}
    require(isinstance(raw, list) and len(raw) <= MAX_TERMS, 'Encoded Pauli term count')
    for word, value in raw:
        require(type(value) in (int, str) and not isinstance(value, bool), 'Integer coefficient numerator')
        require(isinstance(word, list) and all(type(q) is int and q >= 0 and p in ('X', 'Y', 'Z') for q, p in word),
                'Encoded Pauli support')
        key = tuple((q, p) for q, p in word)
        require(key not in result and len({q for q, p in key}) == len(key), 'Unique Pauli support/term')
        result[key] = int(value)
    return result


def blocks(model):
    require(type(model['n']) is int and 1 <= model['n'] <= 12, 'Study register bound before allocation')
    require(0 < len(model['coefficients']) <= MAX_TERMS and 0 < len(model['blocks']) <= MAX_UNITS,
            'Study term/unit count bound')
    denominator = int(model['denominator'])
    require(denominator > 0, 'Positive molecular denominator')
    result = [decode(raw) for raw in model['blocks']]
    total = decode(model['coefficients'])
    require(sum(map(len, result)) == len(total) and {w: c for b in result for w, c in b.items()} == total,
            'Complete units must partition the exact target')
    return [{word: value / denominator for word, value in block.items()} for block in result]


def source_paths():
    own = [p for p in HERE.iterdir() if p.suffix == '.py' and not p.name.startswith('test_')]
    own += [HERE / 'REGISTRY.json', HERE / 'PROTOCOL.md']
    dependencies = [ROOT / name for name in ('analysis/sparse_response_ordering.py', 'analysis/strang_schedule.py',
        'analysis/fermionic_coefficients.py', 'analysis/__init__.py', 'numerical_audit/__init__.py',
        'numerical_audit/followup_classes/exact_blocks.py')]
    return sorted(own + dependencies)


def source_map():
    return {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}


def committed_sources():
    head = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for relative in source_map():
        require(subprocess.check_output(['git', '-C', str(ROOT), 'show', head + ':' + relative]) ==
                (ROOT / relative).read_bytes(), 'Commit scientific source before execution: ' + relative)
    return head


def environment():
    return dict(python=sys.version, executable=str(Path(sys.executable).resolve()), platform=platform.platform(),
        packages={name: importlib.metadata.version(name) for name in
                  ('numpy', 'scipy', 'openfermion', 'openfermionpyscf', 'pyscf', 'h5py')})
