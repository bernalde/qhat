"""Verify and extract the frozen CMC evidence archives without executing them.

Usage: python extract_evidence.py --evidence EVIDENCE --output FRESH_DIRECTORY
Tar names are exactly manifest member paths; archive roots name output folders.
On refusal, retain any partial output and never write the success marker.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import tarfile
import zlib

MANIFEST = "COPY_MANIFEST.a1.json"
SCHEMA = "cost-matched-evidence-archive-v1"
REPORT_SCHEMA = "cost-matched-evidence-extraction-v1"
ROOTS = frozenset(("main-a1", "preflight-a1", "preflight-a2", "preflight-a3", "operators"))
MAX_MEMBERS = 10000
MAX_MEMBER_BYTES = 32 * 1024**2
MAX_TOTAL_BYTES = 128 * 1024**2
MAX_MANIFEST_BYTES = 32 * 1024**2
# Bound tar headers, long-name/PAX records and final padding as well as payloads.
MAX_TAR_OVERHEAD = MAX_MEMBERS * 4096 + 1024**2
MAX_ARCHIVE_BYTES = MAX_TOTAL_BYTES + MAX_TAR_OVERHEAD
BLOCK = 1024**2
HASH = re.compile(r"[0-9a-f]{64}")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    require(type(value) is str and 0 < len(value.encode("utf-8")) <= 1024 and
            "\\" not in value and all(ord(char) >= 32 and ord(char) != 127 for char in value),
            "Bounded POSIX path required")
    path = PurePosixPath(value)
    require(not path.is_absolute() and bool(path.parts) and str(path) == value and
            all(part not in (".", "..") for part in path.parts), "Canonical relative POSIX path required")
    return path


def integer(value, maximum):
    require(type(value) is int and 0 <= value <= maximum, "Bounded nonnegative integer required")
    return value


def digest(value):
    require(type(value) is str and HASH.fullmatch(value) is not None, "Lowercase SHA-256 required")
    return value


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result


def json_integer(value):
    require(len(value.lstrip("-")) <= 20, "JSON integer digit cap")
    return int(value)


def json_float(value):
    result = float(value)
    require(math.isfinite(result), "Nonfinite JSON number")
    return result


def reject_constant(value):
    raise ValueError("Nonfinite JSON constant: " + value)


def parse_manifest(raw):
    depth = 0
    quoted = escaped = False
    for byte in raw:
        if quoted:
            if escaped:
                escaped = False
            elif byte == 92:
                escaped = True
            elif byte == 34:
                quoted = False
        elif byte == 34:
            quoted = True
        elif byte in (91, 123):
            depth += 1
            require(depth <= 32, "JSON nesting cap")
        elif byte in (93, 125):
            depth -= 1
    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique_object,
                      parse_int=json_integer, parse_float=json_float, parse_constant=reject_constant)


def open_input(root, spelling):
    """Open relative input components without following any symlink."""
    parts = canonical(spelling).parts
    flags = os.O_RDONLY | os.O_NOFOLLOW
    directory = os.open(root, flags | os.O_DIRECTORY)
    try:
        for part in parts[:-1]:
            next_directory = os.open(part, flags | os.O_DIRECTORY, dir_fd=directory)
            os.close(directory)
            directory = next_directory
        descriptor = os.open(parts[-1], flags, dir_fd=directory)
    finally:
        os.close(directory)
    try:
        metadata = os.fstat(descriptor)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                "Input must be a regular file without hard links: " + spelling)
        return os.fdopen(descriptor, "rb")
    except BaseException:
        os.close(descriptor)
        raise


def hash_stream(stream, expected_size):
    hasher = hashlib.sha256()
    size = 0
    while True:
        block = stream.read(min(BLOCK, expected_size - size + 1))
        if not block:
            break
        size += len(block)
        require(size <= expected_size, "Input grew beyond declared bytes")
        hasher.update(block)
    require(size == expected_size, "Input byte count differs")
    return hasher.hexdigest()


def admit(manifest):
    require(type(manifest) is dict and manifest.get("schema") == SCHEMA, "Evidence manifest schema")
    archives = manifest.get("archives")
    require(type(archives) is list and len(archives) == len(ROOTS), "Complete bounded archive roster")
    roots, paths = set(), set()
    total_bytes = total_members = 0
    for archive in archives:
        require(type(archive) is dict and
                set(archive) == {"path", "root", "bytes", "sha256", "members", "total_uncompressed_bytes"},
                "Archive metadata fields")
        root = archive["root"]
        canonical(root)
        require(root in ROOTS and root not in roots, "Unknown/duplicate archive root")
        roots.add(root)
        name = archive["path"]
        canonical(name)
        require(name == "archives/" + root + ".tar.gz" and name not in paths, "Canonical unique archive path")
        paths.add(name)
        require(integer(archive["bytes"], MAX_ARCHIVE_BYTES) > 0, "Empty compressed archive")
        digest(archive["sha256"])
        members = archive["members"]
        require(type(members) is list and 0 < len(members) <= MAX_MEMBERS, "Bounded member roster")
        member_paths, member_bytes = set(), 0
        for member in members:
            require(type(member) is dict and set(member) == {"path", "bytes", "sha256", "source"},
                    "Member metadata fields")
            canonical(member["path"])
            require(member["path"] not in member_paths, "Duplicate declared member")
            member_paths.add(member["path"])
            member_bytes += integer(member["bytes"], MAX_MEMBER_BYTES)
            digest(member["sha256"])
            # Provenance is descriptive metadata, never a path to open/execute.
            source = member["source"]
            require(type(source) is str and 0 < len(source.encode("utf-8")) <= 4096 and
                    "\x00" not in source, "Bounded source provenance string")
        for path in member_paths:
            require(all(str(parent) not in member_paths for parent in PurePosixPath(path).parents
                        if str(parent) != "."), "Member file conflicts with a directory prefix")
        require(integer(archive["total_uncompressed_bytes"], MAX_TOTAL_BYTES) == member_bytes,
                "Archive uncompressed total differs")
        total_bytes += member_bytes
        total_members += len(members)
        require(total_bytes <= MAX_TOTAL_BYTES and total_members <= MAX_MEMBERS, "Global extraction limit")
    require(roots == ROOTS, "Missing archive roots")
    for name in ("total_member_bytes", "globaltotal_member_bytes", "total_uncompressed_bytes"):
        if name in manifest:
            require(integer(manifest[name], MAX_TOTAL_BYTES) == total_bytes, "Declared global byte total differs")
    for name in ("total_members", "member_count", "total_member_count"):
        if name in manifest:
            require(integer(manifest[name], MAX_MEMBERS) == total_members, "Declared global member count differs")
    if "source_commit" in manifest:
        require(type(manifest["source_commit"]) is str and
                re.fullmatch(r"[0-9a-f]{40}", manifest["source_commit"]) is not None, "Source commit spelling")
    return archives, total_members, total_bytes


class BoundedTarReader:
    """Bound all decompressed bytes, including tar metadata hidden by tarfile."""
    def __init__(self, stream, maximum):
        self.stream, self.maximum, self.used = stream, maximum, 0

    def read(self, size=-1):
        if size < 0:
            size = self.maximum - self.used + 1
        require(size <= self.maximum - self.used + 1, "Tar read exceeds remaining stream bound")
        data = self.stream.read(size)
        self.used += len(data)
        require(self.used <= self.maximum, "Decompressed tar stream byte cap")
        return data


class BoundedTarInfo(tarfile.TarInfo):
    @classmethod
    def frombuf(cls, buf, encoding, errors):
        info = super().frombuf(buf, encoding, errors)
        if info.type in (tarfile.XHDTYPE, tarfile.XGLTYPE, tarfile.GNUTYPE_LONGNAME):
            # Refuse a giant PAX/long-name body before tarfile reads it eagerly.
            integer(info.size, 16 * 1024)
        else:
            require(info.type in (tarfile.REGTYPE, tarfile.AREGTYPE), "Nonregular tar header")
            integer(info.size, MAX_MEMBER_BYTES)
        return info


def fresh_member(output, spelling):
    parts = canonical(spelling).parts
    directory = output
    for part in parts[:-1]:
        directory = directory / part
        try:
            directory.mkdir(mode=0o700)
        except FileExistsError:
            require(stat.S_ISDIR(directory.lstat().st_mode), "Output directory conflict/symlink")
    descriptor = os.open(directory / parts[-1], os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    return os.fdopen(descriptor, "wb")


def extract_archive(evidence, output, archive):
    expected = {member["path"]: member for member in archive["members"]}
    seen = set()
    destination = output / archive["root"]
    destination.mkdir(mode=0o700)
    with open_input(evidence, archive["path"]) as compressed:
        require(os.fstat(compressed.fileno()).st_size == archive["bytes"], "Archive size differs")
        require(hash_stream(compressed, archive["bytes"]) == archive["sha256"], "Archive SHA-256 differs")
        compressed.seek(0)
        with gzip.GzipFile(fileobj=compressed, mode="rb") as uncompressed:
            bounded = BoundedTarReader(uncompressed, archive["total_uncompressed_bytes"] + MAX_TAR_OVERHEAD)
            with tarfile.open(fileobj=bounded, mode="r|", tarinfo=BoundedTarInfo) as container:
                for info in container:
                    name = info.name
                    canonical(name)
                    require(name in expected and name not in seen, "Unindexed/duplicate archive member")
                    require(info.type in (tarfile.REGTYPE, tarfile.AREGTYPE) and info.sparse is None and
                            not info.linkname, "Only ordinary regular archive members are permitted")
                    member = expected[name]
                    require(integer(info.size, MAX_MEMBER_BYTES) == member["bytes"], "Member byte count differs")
                    seen.add(name)
                    hasher, size = hashlib.sha256(), 0
                    source = container.extractfile(info)
                    require(source is not None, "Unreadable regular archive member")
                    with source, fresh_member(destination, name) as target:
                        while size < member["bytes"]:
                            block = source.read(min(BLOCK, member["bytes"] - size))
                            require(bool(block), "Truncated archive member")
                            size += len(block)
                            hasher.update(block)
                            target.write(block)
                        require(source.read(1) == b"", "Member exceeds declared bytes")
                    require(hasher.hexdigest() == member["sha256"], "Member SHA-256 differs: " + name)
                # Include tarfile's prefetched bytes and validate gzip CRC/EOF.
                # A second tar or other hidden nonzero trailing payload refuses.
                while True:
                    block = container.fileobj.read(BLOCK)
                    if not block:
                        break
                    require(not any(block), "Nonzero trailing tar payload")
            require(seen == set(expected), "Missing indexed archive member")
        compressed.seek(0)
        require(hash_stream(compressed, archive["bytes"]) == archive["sha256"], "Archive changed during extraction")
    return {"path": archive["path"], "root": archive["root"], "sha256": archive["sha256"],
            "bytes": archive["bytes"], "members": len(seen),
            "total_uncompressed_bytes": archive["total_uncompressed_bytes"]}


def extract(evidence, output):
    source_path = Path(__file__)
    source_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    require(not evidence.is_symlink() and evidence.is_dir(), "Evidence directory required")
    evidence = evidence.resolve(strict=True)
    with open_input(evidence, MANIFEST) as stream:
        require(os.fstat(stream.fileno()).st_size <= MAX_MANIFEST_BYTES, "Manifest byte cap")
        raw = stream.read(MAX_MANIFEST_BYTES + 1)
    require(len(raw) <= MAX_MANIFEST_BYTES, "Manifest grew beyond byte cap")
    manifest = parse_manifest(raw)
    archives, count, payload_bytes = admit(manifest)
    archive_directory = evidence / "archives"
    require(not archive_directory.is_symlink() and archive_directory.is_dir(), "Regular archives directory required")
    require({path.name for path in archive_directory.iterdir()} ==
            {PurePosixPath(archive["path"]).name for archive in archives}, "Unindexed/missing archive file")
    # Reject every bad compressed input before creating any extraction output.
    for archive in archives:
        with open_input(evidence, archive["path"]) as stream:
            require(os.fstat(stream.fileno()).st_size == archive["bytes"] and
                    hash_stream(stream, archive["bytes"]) == archive["sha256"], "Archive bytes/hash differ")
    require(not os.path.lexists(output), "Fresh output directory required")
    output = output.parent.resolve(strict=True) / output.name
    output.mkdir(mode=0o700)
    reports = [extract_archive(evidence, output, archive) for archive in archives]
    with open_input(evidence, MANIFEST) as stream:
        require(hash_stream(stream, len(raw)) == hashlib.sha256(raw).hexdigest(), "Manifest changed during extraction")
    require(hashlib.sha256(source_path.read_bytes()).hexdigest() == source_sha256, "Extractor source changed during verification")
    report = {"schema": REPORT_SCHEMA, "status": "complete", "manifest_path": MANIFEST,
              "manifest_sha256": hashlib.sha256(raw).hexdigest(),
              "extractor_source_sha256": source_sha256,
              "source_commit": manifest.get("source_commit"), "archives": reports,
              "archive_count": len(reports), "member_count": count,
              "total_member_bytes": payload_bytes, "extracted_source_code_executed": False}
    with (output / "EXTRACTION_VERIFIED.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = extract(args.evidence, args.output)
    except (ValueError, TypeError, KeyError, OSError, EOFError, tarfile.TarError, zlib.error,
            OverflowError, RecursionError) as error:
        print(json.dumps({"schema": REPORT_SCHEMA, "status": "refused", "error": str(error)}, allow_nan=False))
        return 2
    print(json.dumps(report, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
