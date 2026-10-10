#!/usr/bin/env python3
"""Compare all 30 historical optimized states with the frozen extension outputs.

Requires NumPy only. No QHAT imports or state propagation are performed.
"""
import argparse
import hashlib
import io
import json
import os
import resource
import zipfile
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
import numpy as np


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024**2), b''):
            h.update(b)
    return h.hexdigest()


def contained(base, relative):
    p = (base / relative).resolve()
    require(p.is_relative_to(base.resolve()), 'Path escapes selected root')
    return p


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive-root', type=Path, required=True)
    p.add_argument('--support', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(not a.output.exists(), 'Choose a new output file')
    support = a.support.resolve()
    archive = a.archive_root.resolve()
    manifest_path = support / 'historical_states_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    original_check_path = archive / 'extension/optimized/archived_state_checks.json'
    require(sha(original_check_path) == manifest['original_check_sha256'],
            'Original comparison record changed')
    original = json.loads(original_check_path.read_text())
    expected = {(r['case'], r['ordering']): r for r in original['checks']}
    require(len(manifest['comparisons']) == 30 and len(expected) == 30,
            'Expected all 30 historical comparisons')
    seen = set()
    archive_hashes = {}
    checks = []
    for r in manifest['comparisons']:
        pair = (r['case'], r['ordering'])
        require(pair not in seen and pair in expected, 'Unexpected or duplicate comparison')
        seen.add(pair)
        prior = expected[pair]
        require(prior['source_sha256'] == r['source_sha256']
                and prior['new_state_sha256'] == r['new_state_sha256'],
                'Supplement provenance does not match the frozen record')
        if r['storage'] == 'supplement':
            source = contained(support, r['supplement_file'])
            require(sha(source) == r['source_sha256'], 'Historical file digest mismatch')
        else:
            require(r['storage'] == 'original_archive_member', 'Unknown source kind')
            zpath = contained(archive, r['original_archive_file'])
            if str(zpath) not in archive_hashes:
                archive_hashes[str(zpath)] = sha(zpath)
            require(archive_hashes[str(zpath)] == r['source_sha256'],
                    'Historical nested archive digest mismatch')
            with zipfile.ZipFile(zpath) as z:
                blob = z.read(r['member'])
            require(hashlib.sha256(blob).hexdigest() == r['member_sha256'],
                    'Historical archive member digest mismatch')
            source = io.BytesIO(blob)
        current = contained(archive, r['extension_state_file'])
        require(sha(current) == r['new_state_sha256'], 'Extension state digest mismatch')
        with np.load(source, allow_pickle=False) as f:
            old = f[r['source_key']]
        with np.load(current, allow_pickle=False) as f:
            new = f['actual']
        require(old.shape == new.shape and old.dtype == new.dtype,
                'Array shape or dtype mismatch')
        bitwise = old.tobytes(order='C') == new.tobytes(order='C')
        require(bitwise, 'Historical state is not bitwise identical')
        checks.append({'case': pair[0], 'ordering': pair[1], 'bitwise_equal': bitwise,
                       'norm_difference': float(np.linalg.norm(old - new)),
                       'source_sha256': r['source_sha256'],
                       'new_state_sha256': r['new_state_sha256']})
    require(seen == set(expected), 'Incomplete historical coverage')
    result = {'status': 'passed', 'count': len(checks),
              'bitwise_equal_count': sum(r['bitwise_equal'] for r in checks),
              'maximum_norm_difference': max(r['norm_difference'] for r in checks),
              'checker_sha256': sha(__file__), 'manifest_sha256': sha(manifest_path),
              'scope': 'all 30 frozen historical-output comparisons; no new propagation',
              'checks': checks}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}))


if __name__ == '__main__':
    main()
