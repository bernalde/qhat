"""Compare prospective geometries with all recorded JSON geometries at a Git base."""
from collections import Counter
import hashlib
import json
import math
import subprocess
import transfer_common as k

BASE = 'b8994cac5e7e07583c293fa2c659a24eb2d51127'


def signature(geometry):
    k.require(isinstance(geometry, list) and 1 <= len(geometry) <= 32, 'Geometry shape bound')
    atoms = []
    for entry in geometry:
        k.require(isinstance(entry, (list, tuple)) and len(entry) == 2 and isinstance(entry[0], str), 'Invalid atom')
        symbol, position = entry
        k.require(len(position) == 3 and all(type(v) in (int, float) and math.isfinite(v) for v in position), 'Invalid coordinates')
        atoms.append((symbol, position))
    pairs = []
    for i, (a, x) in enumerate(atoms):
        for b, y in atoms[i + 1:]:
            pairs.append((*sorted((a, b)), math.sqrt(sum((u - v)**2 for u, v in zip(x, y)))))
    return sorted(Counter(a for a, x in atoms).items()), sorted(pairs)


def same(left, right):
    atoms, pairs = left
    other_atoms, other_pairs = right
    return atoms == other_atoms and len(pairs) == len(other_pairs) and all(
        p[:2] == q[:2] and abs(p[2] - q[2]) <= 1e-9 for p, q in zip(pairs, other_pairs))


def recorded_geometries(root=k.ROOT, base=BASE):
    listing = subprocess.check_output(['git', '-C', str(root), 'ls-tree', '-rlz', base, '--', 'numerical_audit/'])
    blobs = {}
    file_count = 0
    for entry in listing.split(b'\0'):
        if not entry:
            continue
        metadata, name = entry.split(b'\t', 1)
        if not name.endswith(b'.json'):
            continue
        mode, kind, oid, size = metadata.split()
        k.require(kind == b'blob' and mode in (b'100644', b'100755'), 'Unexpected baseline JSON tree entry')
        file_count += 1
        blobs.setdefault(oid.decode(), dict(size=int(size), path=name.decode()))
    k.require(0 < file_count <= 16000 and len(blobs) <= 8192 and
              sum(v['size'] for v in blobs.values()) <= 128 * 1024**2, 'Baseline JSON inventory bound')
    # communicate handles both pipes; duplicate blobs are decoded only once.
    ids = list(blobs)
    packed = subprocess.run(['git', '-C', str(root), 'cat-file', '--batch'],
                            input=('\n'.join(ids) + '\n').encode(), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=True).stdout
    found = []
    cursor = 0
    def visit(value, path):
        if isinstance(value, dict):
            if 'geometry' in value:
                try:
                    value_signature = signature(value['geometry'])
                except (ValueError, TypeError):
                    pass
                else:
                    found.append(dict(path=path, signature=value_signature))
            for child in value.values():
                if isinstance(child, (list, dict)):
                    visit(child, path)
        elif isinstance(value, list):
            for child in value:
                if isinstance(child, (dict, list)):
                    visit(child, path)
    for oid in ids:
        end = packed.index(b'\n', cursor)
        returned, kind, length = packed[cursor:end].split()
        size = int(length)
        k.require(returned.decode() == oid and kind == b'blob' and size == blobs[oid]['size'], 'Baseline Git blob differs')
        data = packed[end + 1:end + 1 + size]
        k.require(hashlib.sha1(b'blob ' + str(size).encode() + b'\0' + data).hexdigest() == oid, 'Baseline blob authentication failed')
        visit(json.loads(data), blobs[oid]['path'])
        cursor = end + 2 + size
    k.require(cursor == len(packed) and found, 'Incomplete or empty baseline geometry scan')
    return found, dict(base=base, JSON_files=file_count, unique_JSON_blobs=len(blobs),
                       recorded_geometry_occurrences=len(found), scanned_JSON_bytes=sum(v['size'] for v in blobs.values()))


def audit(registry):
    prior, metadata = recorded_geometries()
    checked = []
    for case in registry['cases']:
        candidate = signature(case['geometry'])
        matches = sorted({g['path'] for g in prior if same(candidate, g['signature'])})
        k.require(not matches, 'Prospective target matches a recorded baseline geometry: ' + case['case'])
        checked.append(dict(case=case['case'], baseline_matches=matches))
    return dict(status='passed', targets=checked, baseline=metadata, tolerance_Angstrom=1e-9,
                method='Atomiccompositionandlabelledpairdistances;rotation/translation/atom-permutationinvariant',
                scope='Everygeometryrecordedinthepublishedbaselinesnumerical-auditJSON;noassertionaboutunrecordedexternaldatasets')
