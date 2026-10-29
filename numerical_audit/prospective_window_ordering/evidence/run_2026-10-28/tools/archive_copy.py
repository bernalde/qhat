"""Copy an explicit frozen mapping; no discovery, rewriting, publication or resume."""
import argparse, hashlib, json, os, re, stat, sys
from pathlib import PurePosixPath
MAX_FILES, MAX_BYTES, MAX_MAP = 16384, 256 * 1024**2, 8 * 1024**2
CHUNK = 1024**2
DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
RESERVED = {"COPY_MANIFEST.json", "COPY_REFUSAL.json"}
def require(ok, message):
    if not ok: raise ValueError(message)
def relative(value):
    require(isinstance(value, str) and 0 < len(value) <= 1024, "bounded relative path")
    require("\0" not in value and "\\" not in value and all(p not in ("", ".", "..") for p in value.split("/")), "canonical relative path")
    return value.split("/")
def directory(path):
    require(path.startswith("/") and "\0" not in path and all(p not in ("", ".", "..") for p in path.split("/")[1:]), "absolute directory")
    fd = os.open("/", DIR_FLAGS)
    try:
        for part in path.split("/")[1:]:
            nxt = os.open(part, DIR_FLAGS, dir_fd=fd)
            os.close(fd); fd = nxt
        return fd
    except BaseException:
        os.close(fd)
        raise
def file_at(root_fd, name, write=False):
    parts, fd = relative(name), os.dup(root_fd)
    try:
        for part in parts[:-1]:
            if write:
                try:
                    os.mkdir(part, 0o700, dir_fd=fd)
                except FileExistsError:
                    pass
            nxt = os.open(part, DIR_FLAGS, dir_fd=fd)
            os.close(fd); fd = nxt
        flags = os.O_NOFOLLOW | os.O_CLOEXEC
        flags |= os.O_WRONLY | os.O_CREAT | os.O_EXCL if write else os.O_RDONLY | os.O_NONBLOCK
        result = os.open(parts[-1], flags, 0o600, dir_fd=fd)
        if not stat.S_ISREG(os.fstat(result).st_mode):
            os.close(result)
            raise ValueError("regular files only")
        return os.fdopen(result, "wb" if write else "rb")
    finally:
        os.close(fd)
def serialized(value):
    return (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()
def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate mapping JSON key")
        result[key] = value
    return result
def copy_one(src_fd, dst_fd, entry):
    digest, count = hashlib.sha256(), 0
    with file_at(src_fd, entry["source"]) as source:
        before = os.fstat(source.fileno())
        require(before.st_size == entry["bytes"], "source size differs from freeze")
        with file_at(dst_fd, entry["destination"], True) as target:
            while chunk := source.read(CHUNK):
                count += len(chunk)
                require(count <= entry["bytes"], "source grew beyond freeze")
                digest.update(chunk); target.write(chunk)
            target.flush(); os.fsync(target.fileno())
        after = os.fstat(source.fileno())
        require((before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                (after.st_size, after.st_mtime_ns, after.st_ctime_ns), "source changed during copy")
    require(count == entry["bytes"] and digest.hexdigest() == entry["sha256"], "source digest differs")
    digest, count = hashlib.sha256(), 0
    with file_at(dst_fd, entry["destination"]) as copied:
        while chunk := copied.read(CHUNK):
            count += len(chunk)
            require(count <= entry["bytes"], "copied size exceeds freeze")
            digest.update(chunk)
    require(count == entry["bytes"] and digest.hexdigest() == entry["sha256"], "copied bytes differ")
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root")
    parser.add_argument("mapping", help="canonical path relative to source_root")
    parser.add_argument("destination", help="fresh absolute directory; parent must exist")
    parser.add_argument("--mapping-sha256", required=True)
    parser.add_argument("--max-files", type=int, default=MAX_FILES)
    parser.add_argument("--max-bytes", type=int, default=MAX_BYTES)
    args = parser.parse_args()
    require(0 < args.max_files <= MAX_FILES and 0 < args.max_bytes <= MAX_BYTES, "bounded copy caps")
    require(re.fullmatch("[0-9a-f]{64}", args.mapping_sha256), "mapping SHA-256 spelling")
    src_fd, dst_fd, done = None, None, []
    try:
        src_fd = directory(args.source_root)
        with file_at(src_fd, args.mapping) as stream:
            require(os.fstat(stream.fileno()).st_size <= MAX_MAP, "mapping size cap")
            raw = stream.read(MAX_MAP + 1)
        require(len(raw) <= MAX_MAP and hashlib.sha256(raw).hexdigest() == args.mapping_sha256, "mapping bytes differ from freeze")
        plan = json.loads(raw, object_pairs_hook=unique_object)
        require(type(plan) is dict and set(plan) == {"schema", "files"} and plan["schema"] == "explicit-archive-copy-map-v1", "mapping schema")
        entries = plan["files"]
        require(type(entries) is list and 0 < len(entries) <= args.max_files - 2, "file count cap, including receipt reserve")
        destinations, total = set(), 0
        for entry in entries:
            require(type(entry) is dict and set(entry) == {"source", "destination", "bytes", "sha256", "origin_logical_path"}, "entry schema")
            for field in ("source", "destination", "origin_logical_path"): relative(entry[field])
            require(entry["destination"].split("/")[0] not in RESERVED and
                    entry["destination"] not in destinations, "reserved/duplicate destination")
            destinations.add(entry["destination"])
            require(type(entry["bytes"]) is int and 0 <= entry["bytes"] <= args.max_bytes and
                    isinstance(entry["sha256"], str) and re.fullmatch("[0-9a-f]{64}", entry["sha256"]), "file freeze")
            total += entry["bytes"]
        require(not any(str(p) in destinations for d in destinations for p in PurePosixPath(d).parents), "file/parent collision")
        manifest = dict(schema="explicit-archive-copy-manifest-v1", status="complete", mapping=args.mapping,
                        mapping_sha256=args.mapping_sha256, max_files=args.max_files, max_bytes=args.max_bytes, payload_bytes=total, files=entries)
        body = serialized(manifest)
        require(len(body) <= MAX_MAP and total + 2 * len(body) + 8192 <= args.max_bytes, "archive byte cap including receipt reserve")
        parent, leaf = args.destination.rsplit("/", 1)
        relative(leaf)
        parent_fd = directory(parent)
        try:
            os.mkdir(leaf, 0o700, dir_fd=parent_fd)
            dst_fd = os.open(leaf, DIR_FLAGS, dir_fd=parent_fd)
        finally:
            os.close(parent_fd)
        for entry in entries:
            copy_one(src_fd, dst_fd, entry)
            done.append(entry["destination"])
        with file_at(dst_fd, "COPY_MANIFEST.json", True) as stream:
            stream.write(body)
        print("complete", len(done), "files", total, "payload bytes")
        return 0
    except (OSError, ValueError, RecursionError) as error:
        if dst_fd is not None:
            receipt = dict(status="incomplete", mapping_sha256=args.mapping_sha256, copied_files=entries[:len(done)], reason=str(error)[:500])
            with file_at(dst_fd, "COPY_REFUSAL.json", True) as stream:
                stream.write(serialized(receipt))
        print("refused", type(error).__name__, file=sys.stderr)
        return 2
    finally:
        if src_fd is not None:
            os.close(src_fd)
        if dst_fd is not None:
            os.close(dst_fd)
if __name__ == "__main__":
    raise SystemExit(main())
