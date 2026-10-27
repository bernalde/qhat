"""Bounded baseline geometry authentication, adapted from transfer_holdouts."""
from collections import Counter
import hashlib
import json
import math
import subprocess
import common as k


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
    k.require(0 < len(blobs) <= 12000 and max(size for size, _ in blobs.values()) <= 16 * 1024**2 and
              sum(size for size, _ in blobs.values()) <= 512 * 1024**2, 'Geometry scan storage bounds')
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
    for oid, (size, name) in blobs.items():
        raw = subprocess.check_output(['git', '-C', str(k.ROOT), 'cat-file', 'blob', oid])
        k.require(len(raw) == size and hashlib.sha1(b'blob ' + str(size).encode() + b'\0' + raw).hexdigest() == oid,
                  'Baseline blob authentication failed')
        visit(json.loads(raw), name)
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
