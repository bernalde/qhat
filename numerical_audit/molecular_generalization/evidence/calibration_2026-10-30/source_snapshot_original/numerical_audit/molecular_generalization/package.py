"""Self-contained evidence archive and bounded, byte-verified extraction.

Only the named DATA tree enters the archive. Local coordination/reviewer client
records and unrelated repositories are outside this input boundary.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import common as k

MAX_FILES = 20000
MAX_BYTES = 512*1024**2


def inventory(data):
    files, total = {}, 0
    for path in sorted(data.rglob('*')):
        k.require(not path.is_symlink(), 'Evidence tree must not contain symlinks')
        if not path.is_file() or path == data/'EVIDENCE_INVENTORY.json':
            continue
        name = str(path.relative_to(data))
        k.contained(data, name)
        size = path.stat().st_size
        total += size
        k.require(len(files) < MAX_FILES and total <= MAX_BYTES, 'Declared evidence file/byte admission')
        files[name] = dict(sha256=k.sha(path), bytes=size)
    return dict(schema='molecular-evidence-inventory-v1', files=files, uncompressed_bytes=total,
        source_freeze_sha256=k.sha(data/'INPUT_PROTOCOL_FROZEN.json'))


def pack(data, output):
    data, output = data.resolve(), output.resolve()
    k.require(not output.exists() and not output.is_relative_to(data), 'Fresh archive outside DATA')
    contents = inventory(data)
    index = data/'EVIDENCE_INVENTORY.json'
    if index.exists():
        k.require(k.read(index) == contents, 'Existing inventory differs; retain prior evidence')
    else:
        k.write_new(index, contents)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Deterministic metadata; bytes of scientific evidence remain unchanged.
    with output.open('xb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', filename='', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as archive:
            for name in sorted([*contents['files'], 'EVIDENCE_INVENTORY.json']):
                path = k.contained(data, name)
                info = tarfile.TarInfo('data/'+name)
                info.size, info.mode, info.mtime = path.stat().st_size, 0o644, 0
                with path.open('rb') as raw:
                    archive.addfile(info, raw)
    k.write_new(output.with_name(output.name+'.json'), dict(archive_sha256=k.sha(output), bytes=output.stat().st_size,
        inventory_sha256=k.sha(index), source_freeze_sha256=contents['source_freeze_sha256']))


def extract(archive, output):
    k.require(not output.exists(), 'Fresh relocated extraction')
    output.mkdir(parents=True)
    seen, total = set(), 0
    with tarfile.open(archive, mode='r|gz') as stream:
        for info in stream:
            k.require(info.isfile() and info.name.startswith('data/') and info.name not in seen
                      and len(seen) < MAX_FILES+1, 'Only unique regular evidence files')
            target = k.contained(output, info.name)
            total += info.size
            k.require(0 <= info.size <= MAX_BYTES and total <= MAX_BYTES+32*1024**2, 'Extraction byte admission')
            seen.add(info.name)
            target.parent.mkdir(parents=True, exist_ok=True)
            with stream.extractfile(info) as source, target.open('xb') as dest:
                remaining = info.size
                while remaining:
                    raw = source.read(min(1024**2, remaining))
                    k.require(raw, 'Truncated archive member')
                    dest.write(raw); remaining -= len(raw)
    data = output/'data'
    declared = k.read(data/'EVIDENCE_INVENTORY.json')
    k.require(inventory(data) == declared, 'Every relocated evidence file authenticated')
    k.write_new(output/'EXTRACTION_VERIFIED.json', dict(complete=True, archive_sha256=k.sha(archive),
        inventory_sha256=k.sha(data/'EVIDENCE_INVENTORY.json'), files=len(declared['files']), bytes=total))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('pack', 'extract'))
    parser.add_argument('source', type=Path); parser.add_argument('output', type=Path)
    args = parser.parse_args()
    (pack if args.mode == 'pack' else extract)(args.source, args.output)
