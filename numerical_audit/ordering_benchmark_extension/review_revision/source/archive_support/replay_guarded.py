#!/usr/bin/env python3
"""Guard the frozen production replay against merged namespace-package imports.

Run from a neutral directory with a clean environment, not a live repository's
parent. This wrapper does not alter the frozen archive's replay source.
"""
import argparse
import importlib.util
import json
import os
import resource
import sys
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def check_qhat_imports(qhat_root):
    import qhat
    paths = [Path(p).resolve() for p in qhat.__path__]
    require(len(paths) == 1 and paths[0] == qhat_root,
            'QHAT namespace must have exactly one archived path; got ' + repr(paths))
    checked = {}
    for name, mod in list(sys.modules.items()):
        if name == 'qhat' or name.startswith('qhat.'):
            file = getattr(mod, '__file__', None)
            if file:
                resolved = Path(file).resolve()
                require(resolved.is_relative_to(qhat_root),
                        'QHAT module escaped archive: ' + name + ': ' + str(resolved))
                checked[name] = str(resolved.relative_to(qhat_root))
            for path in getattr(mod, '__path__', ()):
                require(Path(path).resolve().is_relative_to(qhat_root),
                        'QHAT namespace escaped archive: ' + name + ': ' + path)
    return {'qhat_namespace_paths': list(map(str, paths)), 'module_paths': checked}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive-root', type=Path, required=True)
    p.add_argument('--check-imports-only', action='store_true')
    p.add_argument('--campaign')
    p.add_argument('--case')
    p.add_argument('--ordering')
    p.add_argument('--budget', type=int)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(not a.output.exists(), 'Choose a new output file')
    root = a.archive_root.resolve()
    qhat = root / 'qhat'
    ext = root / 'extension'
    require(qhat.is_dir() and ext.is_dir(), 'Expected extracted qhat/ and extension/')
    require(not Path.cwd().resolve().is_relative_to(qhat),
            'Run from a neutral directory outside the QHAT package')
    os.environ['QHAT_ROOT'] = str(qhat)
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[key] = '1'
    sys.path.insert(0, str(root))
    before = check_qhat_imports(qhat)
    if a.check_imports_only:
        resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
        result = {'status': 'passed', 'scope': 'namespace/path check only; no numerical replay',
                  **before}
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result))
        return
    require(all(value is not None for value in (a.campaign, a.case, a.ordering, a.budget)),
            'Production replay requires campaign, case, ordering and budget')
    pending = a.output.with_name(a.output.name + '.guard-pending')
    require(not pending.exists(), 'Pending output exists; choose a new output file')
    frozen = ext / 'packaging/replay_one.py'
    spec = importlib.util.spec_from_file_location('_frozen_extension_replay', frozen)
    require(spec is not None and spec.loader is not None, 'Cannot load frozen replay')
    replay = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(replay)
    sys.argv = [str(frozen), '--extension', str(ext), '--qhat-root', str(qhat),
                '--campaign', a.campaign, '--case', a.case, '--ordering', a.ordering,
                '--budget', str(a.budget), '--output', str(pending.resolve())]
    replay.main()
    after = check_qhat_imports(qhat)
    result = json.loads(pending.read_text())
    result['namespace_guard'] = {'before': before, 'after': after}
    pending.write_text(json.dumps(result, indent=2) + '\n')
    pending.replace(a.output)


if __name__ == '__main__':
    main()
