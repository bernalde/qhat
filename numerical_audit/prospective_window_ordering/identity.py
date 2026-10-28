"""Bounded baseline geometry authentication, adapted from transfer_holdouts."""
from collections import Counter
import hashlib
import json
import math
import subprocess
import study_common as k

MAX_JSON_BLOBS = 16384
MAX_JSON_BYTES = 16 * 1024**2
MAX_TOTAL_JSON_BYTES = 512 * 1024**2


def validate_inventory(blobs):
    k.require(0 < len(blobs) <= MAX_JSON_BLOBS and
              all(type(size) is int and 0 <= size <= MAX_JSON_BYTES for size, _ in blobs.values()) and
              sum(size for size, _ in blobs.values()) <= MAX_TOTAL_JSON_BYTES,
              'Geometry scan storage bounds')
    k.require(all(isinstance(oid, str) and len(oid) == 40 and
                  all(c in '0123456789abcdef' for c in oid) for oid in blobs), 'Invalid baseline blob identity')


def stream_blobs(blobs, root=k.ROOT, *, popen=None):
    """Authenticate one bounded Git blob at a time through one owned process."""
    validate_inventory(blobs)
    factory = subprocess.Popen if popen is None else popen
    process = factory(['git', '-C', str(root), 'cat-file', '--batch'],
                      stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    try:
        for oid, (size, name) in blobs.items():
            process.stdin.write(oid.encode('ascii') + b'\n')
            process.stdin.flush()
            header = process.stdout.readline(256)
            k.require(0 < len(header) <= 255 and header.endswith(b'\n'), 'Bounded Git batch header')
            fields = header[:-1].split()
            k.require(len(fields) == 3 and fields[0] == oid.encode('ascii') and fields[1] == b'blob' and
                      fields[2].isdigit() and len(fields[2]) <= 8, 'Git batch header identity/type/size')
            announced = int(fields[2])
            k.require(announced == size and 0 <= announced <= MAX_JSON_BYTES,
                      'Git batch header size differs from inventory')
            raw = process.stdout.read(announced)
            k.require(len(raw) == announced and process.stdout.read(1) == b'\n', 'Truncated Git batch blob')
            k.require(hashlib.sha1(b'blob ' + str(announced).encode() + b'\0' + raw).hexdigest() == oid,
                      'Baseline blob authentication failed')
            yield name, raw
        process.stdin.close()
        k.require(process.wait(timeout=5) == 0, 'Git batch process failed')
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        for stream in (process.stdin, process.stdout):
            if stream is not None and not stream.closed:
                stream.close()


def signature(geometry):
    k.require(isinstance(geometry, list) and 1 <= len(geometry) <= 32, 'Geometry atom count')
    atoms = []
    for symbol, position in geometry:
        k.require(isinstance(symbol, str) and len(position) == 3 and
                  all(type(x) in (int, float) and math.isfinite(x) for x in position), 'Invalid geometry')
        atoms.append((symbol, position))
    pairs = []
    for i, (a, x) in enumerate(atoms):
        for b, y in atoms[i + 1:]:
            pairs.append((*sorted((a, b)), math.sqrt(sum((u - v)**2 for u, v in zip(x, y)))))
    return sorted(Counter(a for a, x in atoms).items()), sorted(pairs)


def same(a, b):
    return a[0] == b[0] and len(a[1]) == len(b[1]) and all(
        p[:2] == q[:2] and abs(p[2] - q[2]) <= 1e-9 for p, q in zip(a[1], b[1]))


def audit(registry):
    base = registry['base_commit']
    listing = subprocess.check_output(['git', '-C', str(k.ROOT), 'ls-tree', '-rlz', base, '--', 'numerical_audit/'])
    blobs = {}
    for entry in listing.split(b'\0'):
        if not entry:
            continue
        metadata, name = entry.split(b'\t', 1)
        if name.endswith(b'.json'):
            mode, kind, oid, size = metadata.split()
            k.require(kind == b'blob' and mode in (b'100644', b'100755'), 'Unexpected JSON tree entry')
            blobs.setdefault(oid.decode(), (int(size), name.decode()))
    validate_inventory(blobs)
    observed, malformed = [], []
    def visit(value, name):
        if isinstance(value, dict):
            if 'geometry' in value:
                try:
                    observed.append((name, signature(value['geometry'])))
                except (TypeError, ValueError):
                    malformed.append(name)
            for child in value.values():
                if isinstance(child, (dict, list)):
                    visit(child, name)
        elif isinstance(value, list):
            for child in value:
                if isinstance(child, (dict, list)):
                    visit(child, name)
    stream = stream_blobs(blobs)
    try:
        for name, raw in stream:
            visit(json.loads(raw), name)
    finally:
        stream.close()
    k.require(observed and not malformed, 'Incomplete baseline geometry scan')
    rows = []
    for case in registry['cases']:
        candidate = signature(case['geometry'])
        matches = sorted({name for name, prior in observed if same(candidate, prior)})
        k.require(not matches, 'Proposed case already occurs in recorded baseline: ' + case['case'])
        rows.append(dict(case=case['case'], matches=matches))
    return dict(base=base, cases=rows, unique_JSON_blobs=len(blobs), geometry_occurrences=len(observed),
                scanned_bytes=sum(size for size, _ in blobs.values()), tolerance_Angstrom=1e-9,
                scope='Recorded baseline numerical-audit JSON; no claim about external or unrecorded data')
