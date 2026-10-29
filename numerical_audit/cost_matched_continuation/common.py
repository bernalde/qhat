"""Shared declarations and exact-source loading for the continuation study."""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = '646c029d71fa7abd2c06ae64693843f4aeb5f1e8'
OLD = ROOT / 'numerical_audit/prospective_window_ordering'
TIMES = (.65, 1.35)
CAPS = (8, 24)
BUDGETS = (32, 64)
REPETITIONS = 2
ARMS = ('full', 'window')
PARENT_CPU = 8.
PARENT_WALL = 40
MAX_EVALUATIONS = 1_000_000
MAX_MOVES = 1000
CLEANUP_RESERVE = .25


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    path = Path(path)
    require(path.is_file() and path.stat().st_size <= 32*1024**2, 'Bounded JSON input')
    return json.loads(path.read_text())


def write_new(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.partial')
    require(not path.exists() and not temporary.exists(), 'Preserve existing output/partial')
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    import os
    os.link(temporary, path)
    temporary.unlink()


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    require(Path(module.__file__).resolve() == path.resolve(), 'Wrong helper source')
    return module


prior = _load('_continuation_prior_common', OLD/'study_common.py')


@contextmanager
def aliases(mapping):
    """Old standalone drivers use absolute imports; never alter their sources."""
    previous = {name: sys.modules.get(name) for name in mapping}
    try:
        sys.modules.update(mapping)
        yield
    finally:
        for name, module in previous.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def chemistry_helpers():
    with aliases({'study_common': prior}):
        identity = _load('_continuation_prior_identity', OLD/'identity.py')
        with aliases({'identity': identity}):
            producer = _load('_continuation_prior_chemistry', OLD/'generate_inputs.py')
    require(producer.k is prior and identity.k is prior and prior.ROOT == ROOT,
            'Chemistry helper identity mismatch')
    return producer, identity


def configuration():
    return read(HERE/'PROTOCOL.json')


def registry():
    return read(HERE/'REGISTRY.json')


def source_paths():
    own = [path for path in HERE.rglob('*') if path.is_file()
           and 'evidence' not in path.relative_to(HERE).parts
           and not {'__pycache__', 'test'} & set(path.relative_to(HERE).parts)
           and path.name != 'analyze.py'
           and (path.suffix in ('.py', '.jl', '.json', '.toml')
                or path.name in ('README.md', 'PROTOCOL.md'))]
    return sorted(set(own + prior.source_paths()))


def source_map():
    return {str(path.relative_to(ROOT)): sha(path) for path in source_paths()}


def committed_sources():
    head = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for name, digest in source_map().items():
        content = subprocess.check_output(['git', '-C', str(ROOT), 'show', head+':'+name])
        require(hashlib.sha256(content).hexdigest() == digest,
                'Commit all study/dependency source before execution: '+name)
    return head


def bootstrap():
    ordering = prior.bootstrap()
    require(Path(ordering.__file__).resolve() == ROOT/'analysis/sparse_response_ordering.py',
            'Wrong optimizer import')
    return ordering


def validate_protocol(value):
    require(value['schema'] == 'cost-matched-continuation-v1', 'Protocol schema')
    require(value['total_CPU_budgets'] == list(BUDGETS), 'CPU budget grid')
    require(value['physical']['times'] == list(TIMES) and
            value['physical']['rotation_cap_multipliers'] == list(CAPS), 'Physical grid')
    require(value['repetitions'] == REPETITIONS and value['arms'] == list(ARMS), 'Arm/repetition grid')
    require(value['training_time'] == 1. and value['training_rotation_multiplier'] == 30,
            'Objective replay requires T1/30U')
    require(value['native']['max_evaluations'] == MAX_EVALUATIONS and
            value['native']['max_moves'] == MAX_MOVES, 'Continuation evaluation/move ceiling')
    require(value['window']['width'] == 8 and value['window']['parent_method'] == 'magnitude',
            'Exact supplied-parent window continuation')
    return value


def frozen_protocol(output):
    value = read(Path(output)/'INPUT_PROTOCOL_FROZEN.json')
    require(value['sources'] == source_map(), 'Source changed after input freeze')
    require(value['registry'] == registry() and value['configuration'] == configuration(),
            'Registry/protocol changed after freeze')
    require(value['environment'] == prior.environment(), 'Numerical environment changed')
    return value


def validate_order(order, size):
    require(isinstance(order, list) and all(type(x) is int for x in order)
            and sorted(order) == list(range(size)), 'Complete zero-based original-label order')
    return order


def write_toml(path, fields, tables):
    """Write the bounded scalar/array scientific schema using the existing encoder."""
    ordering = bootstrap()
    from qhat.analysis.window_ordering import _toml
    require(Path(ordering.__file__).resolve().is_relative_to(ROOT), 'Wrong TOML boundary source')
    lines = [key+'='+_toml(value) for key, value in fields.items()]
    for name, records in tables.items():
        for row in records:
            lines.append('[['+name+']]')
            lines.extend(key+'='+_toml(value) for key, value in row.items())
    raw = ('\n'.join(lines)+'\n').encode('ascii')
    require(len(raw) <= 8*1024**2, 'TOML output byte bound')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('xb') as stream:
        stream.write(raw)
