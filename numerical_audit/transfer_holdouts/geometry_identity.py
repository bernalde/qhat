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
              max(v['size'] for v in blobs.values()) <= 16 * 1024**2 and
              sum(v['size'] for v in blobs.values()) <= 512 * 1024**2, 'Baseline JSON streaming bound')
    # Read exactly one authenticated blob at a time. The total source inventory
    # exceeds128MiB, but no aggregate copy or parsed aggregate is retained.
    ids = list(blobs)
    found = []
    unparsed = []
    def visit(value, path):
        if isinstance(value, dict):
            if 'geometry' in value:
                try:
                    value_signature = signature(value['geometry'])
                except (ValueError, TypeError):
                    unparsed.append(dict(path=path, type=type(value['geometry']).__name__))
                else:
                    found.append(dict(path=path, signature=value_signature))
            for child in value.values():
                if isinstance(child, (list, dict)):
                    visit(child, path)
        elif isinstance(value, list):
            for child in value:
                if isinstance(child, (dict, list)):
                    visit(child, path)
    child = subprocess.Popen(['git', '-C', str(root), 'cat-file', '--batch'],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        for oid in ids:
            child.stdin.write((oid + '\n').encode())
            child.stdin.flush()
            header = child.stdout.readline()
            returned, kind, length = header.split()
            size = int(length)
            k.require(returned.decode() == oid and kind == b'blob' and size == blobs[oid]['size'],
                      'Baseline Git blob differs')
            data = child.stdout.read(size)
            k.require(len(data) == size and child.stdout.read(1) == b'\n', 'Truncated baseline Git blob')
            k.require(hashlib.sha1(b'blob ' + str(size).encode() + b'\0' + data).hexdigest() == oid,
                      'Baseline blob authentication failed')
            value = json.loads(data)
            visit(value, blobs[oid]['path'])
            del value, data
        child.stdin.close()
        k.require(child.wait() == 0 and found, 'Incomplete or empty baseline geometry scan')
    finally:
        if child.poll() is None:
            child.kill()
            child.wait()
    k.require(not unparsed, 'Unparsed baseline geometry fields require inspection: ' + repr(unparsed[:10]))
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
