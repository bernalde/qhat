#!/usr/bin/env python3
"""Stream-check an archive against external identity and its internal inventory.

All checks remain active under python -O. Integrity-only mode establishes internal
consistency, not identity against a separately retained trusted release record.
"""
import argparse
import hashlib
import json
import re
import stat
import zipfile
from pathlib import Path, PurePosixPath


def require(ok, message):
    if not ok:
        raise ValueError(message)


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def safe_name(name):
    p = PurePosixPath(name)
    require(bool(name) and not p.is_absolute() and '..' not in p.parts
            and '\\' not in name and str(p) == name and name not in ('.', '..'),
            'Unsafe member ' + name)


def unique_object(items):
    value = {}
    for key, item in items:
        require(key not in value, 'Duplicate manifest key ' + key)
        value[key] = item
    return value


def verify(path, expected=None, *, expected_manifest_sha256=None,
           expected_members=None, integrity_only=False):
    require((expected is not None) != integrity_only,
            'Supply an expected archive SHA-256 or explicitly select integrity-only mode')
    for value in (expected, expected_manifest_sha256):
        require(value is None or re.fullmatch('[0-9a-f]{64}', value),
                'Expected SHA-256 must have 64 lowercase hexadecimal characters')
    if expected_members is not None:
        require(type(expected_members) is int and expected_members >= 0,
                'Expected member count must be nonnegative')
    digest = file_hash(path)
    if expected is not None:
        require(digest == expected, 'Archive SHA-256 mismatch')
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        require(len(names) == len(set(names)), 'Duplicate member')
        for info in z.infolist():
            safe_name(info.filename)
            require(not stat.S_ISLNK(info.external_attr >> 16), 'Symlink member')
        manifest_bytes = z.read('artifact_manifest.json')
        manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()
        if expected_manifest_sha256 is not None:
            require(manifest_sha == expected_manifest_sha256,
                    'Manifest SHA-256 mismatch')
        manifest = json.loads(manifest_bytes, object_pairs_hook=unique_object)
        entries = manifest['files']
        require(isinstance(entries, dict), 'Invalid manifest file inventory')
        require('artifact_manifest.json' not in entries, 'Manifest cannot list itself')
        require(set(names) == set(entries) | {'artifact_manifest.json'},
                'Unlisted or missing members')
        if expected_members is not None:
            require(len(entries) == expected_members, 'Inventory member count mismatch')
        total = 0
        for name, entry in entries.items():
            safe_name(name)
            require(type(entry['bytes']) is int and entry['bytes'] >= 0,
                    'Invalid declared size ' + name)
            require(isinstance(entry['sha256'], str)
                    and re.fullmatch('[0-9a-f]{64}', entry['sha256']),
                    'Invalid declared digest ' + name)
            require(z.getinfo(name).file_size == entry['bytes'],
                    'Declared size mismatch ' + name)
            h = hashlib.sha256()
            size = 0
            with z.open(name) as f:
                for block in iter(lambda: f.read(4 * 1024**2), b''):
                    h.update(block)
                    size += len(block)
            require(size == entry['bytes'] and h.hexdigest() == entry['sha256'],
                    'Member digest/size mismatch ' + name)
            total += size
    return {'status': 'passed', 'archive_sha256': digest,
            'manifest_sha256': manifest_sha, 'members': len(entries),
            'zip_entries_including_manifest': len(names),
            'uncompressed_bytes': total,
            'identity_mode': 'integrity_only' if integrity_only else 'external_sha256',
            'scope': 'byte integrity; no simulation or proof checking'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive', type=Path)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--expected-sha256')
    g.add_argument('--integrity-only', action='store_true')
    p.add_argument('--expected-manifest-sha256')
    p.add_argument('--expected-members', type=int,
                   help='Number of inventoried data members, excluding artifact_manifest.json')
    a = p.parse_args()
    print(json.dumps(verify(a.archive, a.expected_sha256,
          expected_manifest_sha256=a.expected_manifest_sha256,
          expected_members=a.expected_members, integrity_only=a.integrity_only)))


if __name__ == '__main__':
    main()
