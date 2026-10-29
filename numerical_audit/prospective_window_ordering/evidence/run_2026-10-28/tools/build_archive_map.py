"""Root-only metadata/byte closure builder; never copies an archive or rewrites records."""
import argparse, csv, datetime as dt, hashlib, io, json, os, re, stat, subprocess
from pathlib import Path

ROOT = "runs/prospective_window_ordering_2026-10-28/a1"
INTEGRATION = "runs/integration_2026-10-28"
HEX = re.compile(r"[0-9a-f]{64}")
MAX_FILES, MAX_BYTES, MAX_JSON = 16382, 256 * 1024**2, 32 * 1024**2
MAX_HASH, MAX_RECOVERY, MAX_BLOB, MAX_GIT = 1024**3, 32 * 1024**2, 2 * 1024**2, 128
BAD = {"__pycache__", ".cache", "cache", "fontcache", ".pytest_cache", ".mypy_cache"}
REDS = ("source-review-independent-a1", "source-review-preflight-links-red-a1",
        "source-consumer-allocation-red-a1", "analysis-source-before-primary-bindings")

def require(ok, reason):
    if not ok: raise ValueError(reason)

def rel(value):
    require(isinstance(value, str) and 0 < len(value) <= 1024 and "\0" not in value and
            "\\" not in value and all(p not in ("", ".", "..") for p in value.split("/")), "unsafe relative path")
    return value

def instant(value):
    result = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    require(result.tzinfo is not None, "timezone-aware cutoff/timestamp required")
    return result.timestamp()

def encode(value):
    return (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()

def declarations(value):
    result = {}
    for key in ("source_hashes", "files", "test_source_hashes", "sources", "sha256", "changed_files"):
        items = value.get(key, {})
        if isinstance(items, dict):
            result.update({k: v if isinstance(v, str) else v.get("sha256") for k, v in items.items()
                           if isinstance(v, (str, dict))})
    if not result and all(isinstance(v, str) and HEX.fullmatch(v) for v in value.values()):
        result = value
    require(all(isinstance(v, str) and HEX.fullmatch(v) for v in result.values()), "invalid declared source digest")
    return result

class Builder:
    def __init__(self, args):
        self.a, self.base, self.cutoff = args, Path(args.r).absolute(), instant(args.cutoff)
        self.rows, self.digests, self.providers, self.gaps, self.late, self.uses = {}, {}, {}, [], [], set()
        self.hashed = self.recovered = self.git_calls = 0
        self.provenance = {}
        self.excluded_temporary = []
        self.record = self.base / rel(args.recovered)
        require(args.output_map.startswith("coordination/") and not self.safe(args.output_map).exists(), "fresh private mapping path")
        require(args.recovered.startswith(INTEGRATION + "/source-closure-recovered-") and
                "/" not in args.recovered[len(INTEGRATION) + 1:], "fresh owned recovery namespace required")
        self.safe(INTEGRATION)
        self.record.mkdir(mode=0o700)

    def safe(self, name, anchor=None):
        path = (self.base if anchor is None else anchor) / rel(name)
        for parent in (path, *path.parents):
            require(not parent.is_symlink(), "symlink source/ancestor refused")
        return path

    def digest(self, path):
        before = path.stat()
        require(stat.S_ISREG(before.st_mode), "regular source required")
        key = (str(path), before.st_size, before.st_mtime_ns, before.st_ctime_ns)
        if key not in self.digests:
            self.hashed += before.st_size
            require(self.hashed <= MAX_HASH, "total hash-byte cap")
            h, count = hashlib.sha256(), 0
            with path.open("rb") as stream:
                while chunk := stream.read(1024**2):
                    count += len(chunk)
                    require(count <= before.st_size, "source grew beyond reserved hash bytes")
                    h.update(chunk)
            after = path.stat()
            require(count == before.st_size and (before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                    (after.st_size, after.st_mtime_ns, after.st_ctime_ns), "source changed during hash")
            self.digests[key] = h.hexdigest()
        return self.digests[key], before.st_size

    def read(self, name):
        path = self.safe(name)
        require(path.stat().st_size <= MAX_JSON, "metadata size cap")
        with path.open("rb") as stream: raw = stream.read(MAX_JSON + 1)
        require(len(raw) <= MAX_JSON, "metadata grew beyond cap")
        return json.loads(raw)

    def add(self, name, destination, logical=None, expected=None, generated=False):
        rel(destination)
        path = self.safe(name)
        require(generated or path.stat().st_mtime <= self.cutoff, "selected source after cutoff")
        sha, size = self.digest(path)
        require(expected is None or sha == expected, "declared source digest mismatch: " + name)
        row = dict(source=name, destination=destination, bytes=size, sha256=sha, origin_logical_path=rel(logical or destination))
        require(destination not in self.rows or self.rows[destination] == row, "archive destination collision")
        self.rows[destination] = row
        self.providers.setdefault(sha, row)
        require(len(self.rows) <= MAX_FILES and sum(v["bytes"] for v in self.rows.values()) <= MAX_BYTES, "archive mapping cap")
        return row

    def public(self, name):
        return next((r["destination"] for r in self.rows.values() if r["source"] == name), None)

    def note(self, kind, **fields):
        self.gaps.append(dict(kind=kind, **fields))

    def write(self, leaf, data, destination=None, logical=None):
        self.recovered += len(data)
        require(self.recovered <= MAX_RECOVERY, "recovered source/metadata byte cap")
        path = self.record / rel(leaf)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream: stream.write(data)
        return self.add(str(path.relative_to(self.base)), destination or leaf, logical, hashlib.sha256(data).hexdigest(), True)

    def metadata_uses(self, value, owner):
        if isinstance(value, dict):
            for child in value.values(): self.metadata_uses(child, owner)
        elif isinstance(value, list):
            for child in value: self.metadata_uses(child, owner)
        elif isinstance(value, str):
            prefix = str(self.base) + "/coordination/"
            if value.startswith(prefix) and value.endswith((".py", ".jl")):
                self.uses.add((owner, str(Path(value).relative_to(self.base))))

    def collect_root(self):
        temporary = {}
        for freeze in sorted(self.safe(ROOT).glob("ALL_ORDERS_FROZEN.*.json")):
            if freeze.stat().st_mtime > self.cutoff: continue
            name = str(freeze.relative_to(self.base))
            for descriptor in self.read(name).get("retained_temporary_and_checkpoint_artifacts", []):
                relative = rel(descriptor["path"])
                require(relative.startswith("tmp/") and HEX.fullmatch(descriptor["sha256"]), "temporary descriptor")
                if any(p in BAD for p in Path(relative).parts) or relative.endswith(".lock"):
                    self.excluded_temporary.append(dict(freeze=name, path=relative, reason="explicit cache/lock exclusion")); continue
                previous = temporary.setdefault(relative, descriptor["sha256"])
                require(previous == descriptor["sha256"], "conflicting retained temporary versions")
        for parent, directories, files in os.walk(self.safe(ROOT), followlinks=False):
            directories[:] = sorted(d for d in directories if d not in BAD)
            for directory in directories: require(not (Path(parent) / directory).is_symlink(), "symlink directory")
            for leaf in sorted(files):
                path = Path(parent) / leaf
                local = str(path.relative_to(self.safe(ROOT)))
                if leaf.endswith(".lock") or any(p in BAD for p in Path(local).parts): continue
                if local.startswith("tmp/") and local not in temporary: continue
                if path.stat().st_mtime > self.cutoff:
                    self.late.append(str(path.relative_to(self.base))); continue
                self.add(str(path.relative_to(self.base)), local, expected=temporary.get(local))
        for local in temporary:
            if local not in self.rows: self.note("retained_temporary_not_archived", path=local)
        plans = [r for r in self.rows.values() if r["destination"].startswith("jobs/") and r["destination"].endswith(".PLAN.json")]
        closed = 0
        for row in plans:
            plan = self.read(row["source"]); self.metadata_uses(plan, row["source"])
            label = plan["label"]
            require(instant(plan["started_UTC"]) <= self.cutoff, "future-started selected job")
            names = ("jobs/" + label + ".FINAL.json", "logs/" + label + ".resources.json", "logs/" + label + ".log")
            if not all(name in self.rows for name in names):
                self.note("unclosed_job", label=label); continue
            final = self.read(ROOT + "/" + names[0])
            require(instant(final["completed_UTC"]) <= self.cutoff and final["label"] == label and
                    final["receipt_sha256"] == self.rows[names[1]]["sha256"], "terminal job/receipt binding")
            closed += 1
        require(len(plans) >= self.a.minimum_jobs, "historical job inventory below explicit minimum")
        return closed

    def collect_integration(self):
        root = self.safe(INTEGRATION)
        for path in sorted(root.iterdir()):
            if path.is_file() and path.suffix in (".json", ".log", ".pylog"):
                if path.stat().st_mtime > self.cutoff: self.late.append(str(path.relative_to(self.base))); continue
                name = str(path.relative_to(self.base))
                self.add(name, "validation/integration_2026-10-28/" + path.name)
                if path.suffix == ".json": self.metadata_uses(self.read(name), name)
        for path in sorted((root / "logs").iterdir()):
            if path.is_file() and path.suffix in (".json", ".log") and path.stat().st_mtime <= self.cutoff:
                name = str(path.relative_to(self.base))
                self.add(name, "validation/integration_2026-10-28/logs/" + path.name)
                if path.suffix == ".json": self.metadata_uses(self.read(name), name)
        for version in range(1, 7):
            manifest = INTEGRATION + "/BUILD_INPUT_MANIFEST.a" + str(version) + ".json"
            files = declarations(self.read(manifest))
            require(len(files) == 134, "build input manifest must retain exactly134 keys")
            prefix = "build-source-a" + str(version) + "/"
            for logical, sha in files.items():
                self.add(INTEGRATION + "/" + prefix + rel(logical), "validation/integration_2026-10-28/" + prefix + logical, logical, sha)
            wheel = "wheels-a" + str(version) + "/qhat-0.1.0-py3-none-any.whl"
            self.add(INTEGRATION + "/" + wheel, "validation/integration_2026-10-28/" + wheel)
        for directory in REDS:
            prefix = INTEGRATION + "/" + directory + "/"
            self.add(prefix + "SOURCE.json", "validation/integration_2026-10-28/" + directory + "/SOURCE.json")
            for logical, sha in declarations(self.read(prefix + "SOURCE.json")).items():
                self.add(prefix + rel(logical), "validation/integration_2026-10-28/" + directory + "/" + logical, logical, sha)

    def resolve(self, logical, sha, direct=None):
        if sha in self.providers: return self.providers[sha], self.provenance.get(sha)
        if direct:
            path = self.safe(direct)
            if path.is_file() and self.digest(path)[0] == sha:
                return self.add(direct, "tools/helpers/" + sha + "/" + Path(direct).name, logical, sha), None
        path = self.safe(logical, Path(self.a.worktree).absolute())
        if path.is_file() and path.stat().st_size <= MAX_BLOB and self.digest(path)[0] == sha:
            with path.open("rb") as stream: raw = stream.read(MAX_BLOB + 1)
            require(len(raw) <= MAX_BLOB and hashlib.sha256(raw).hexdigest() == sha, "worktree source changed after verification")
            self.provenance[sha] = dict(kind="exact_current_worktree_bytes", path=logical)
            return self.write("source/" + sha + "/" + logical, raw, "source_history/" + sha + "/" + logical, logical), self.provenance[sha]
        for commit in self.a.git_ref:
            if self.git_calls >= MAX_GIT: break
            self.git_calls += 1
            spec = commit + ":" + logical
            size = subprocess.run(["git", "-C", self.a.worktree, "cat-file", "-s", spec], capture_output=True, timeout=5)
            if size.returncode or not size.stdout.strip().isdigit() or int(size.stdout) > MAX_BLOB: continue
            blob = subprocess.run(["git", "-C", self.a.worktree, "cat-file", "blob", spec], capture_output=True, timeout=5)
            if blob.returncode or len(blob.stdout) != int(size.stdout) or hashlib.sha256(blob.stdout).hexdigest() != sha: continue
            oid = hashlib.sha1(b"blob " + str(len(blob.stdout)).encode() + b"\0" + blob.stdout).hexdigest()
            row = self.write("source/" + sha + "/" + logical, blob.stdout, "source_history/" + sha + "/" + logical, logical)
            self.provenance[sha] = dict(kind="read_only_git_blob", commit=commit, path=logical, git_sha1_blob=oid)
            return row, self.provenance[sha]
        return None, None

    def bind_sources(self):
        worklist = self.safe(self.a.worklist)
        require(worklist.stat().st_size <= MAX_JSON, "worklist size cap")
        with worklist.open("rb") as stream: raw = stream.read(MAX_JSON + 1)
        require(len(raw) <= MAX_JSON, "worklist grew beyond cap")
        records = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
        require(len(records) <= 4096 and records, "bounded source worklist required")
        manifests, bindings = {}, []
        for item in records:
            original, logical, sha = rel(item["manifest_R_relative"]), rel(item["logical_path"]), item["declared_sha256"]
            require(HEX.fullmatch(sha), "source worklist digest")
            if original not in manifests: manifests[original] = declarations(self.read(original))
            if manifests[original].get(logical) != sha:
                self.note("worklist_manifest_mismatch", manifest=original, logical=logical, sha256=sha); continue
            row, provenance = self.resolve(logical, sha)
            if row is None: self.note("missing_source_version", manifest=original, logical=logical, sha256=sha)
            bindings.append(dict(manifest_original_R=original, manifest_public_path=self.public(original), logical_path=logical,
                                 sha256=sha, public_archive_path=None if row is None else row["destination"],
                                 provider_R_path=None if row is None else row["source"], git_provenance=provenance))
        self.write("SOURCE_VERSION_WORKLIST.tsv", raw)
        return bindings

    def bind_helpers(self):
        covered, results = set(), []
        if self.a.helper_bindings:
            spec = self.read(self.a.helper_bindings)
            require(spec["schema"] == "archive-helper-evidence-v1" and len(spec["bindings"]) <= 512, "helper binding schema")
            for item in spec["bindings"]:
                source, logical, sha = rel(item["source"]), rel(item["logical_path"]), item["sha256"]
                require(source.startswith("coordination/") and HEX.fullmatch(sha), "explicit helper source")
                evidence = self.read(rel(item["evidence"])); value = evidence
                for token in item["pointer"].strip("/").split("/"):
                    value = value[int(token)] if isinstance(value, list) else value[token.replace("~1", "/").replace("~0", "~")]
                require(value == sha and self.public(item["evidence"]) is not None, "public helper digest evidence")
                row, provenance = self.resolve(logical, sha, source)
                uses = [(rel(owner), source) for owner in item["uses"]]
                require(uses and all(use in self.uses for use in uses), "helper use does not match retained metadata")
                if row is None: self.note("missing_helper_version", source=source, sha256=sha)
                else: covered.update(uses)
                results.append(dict(source=source, sha256=sha, evidence=item["evidence"], pointer=item["pointer"],
                                    uses=item["uses"], public_archive_path=None if row is None else row["destination"], git_provenance=provenance))
        for owner, source in sorted(self.uses - covered): self.note("unbound_helper_use", owner=owner, source=source)
        return results

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("r", "worktree", "cutoff", "recovered", "output_map", "worklist"): parser.add_argument("--" + name.replace("_", "-"), required=True)
    parser.add_argument("--helper-bindings")
    parser.add_argument("--extra-files")
    parser.add_argument("--git-ref", action="append", default=[])
    parser.add_argument("--minimum-jobs", type=int, default=184)
    args = parser.parse_args()
    require(0 <= args.minimum_jobs <= MAX_FILES and len(args.git_ref) <= 16 and
            all(re.fullmatch("[0-9a-f]{40}", ref) for ref in args.git_ref), "finite explicit job/Git scope")
    builder = Builder(args)
    jobs = builder.collect_root(); builder.collect_integration()
    if args.extra_files:
        extra = builder.read(args.extra_files)
        require(extra["schema"] == "archive-explicit-extra-files-v1" and len(extra["files"]) <= 1024, "extra files schema")
        for item in extra["files"]:
            row = builder.add(rel(item["source"]), rel(item["destination"]), rel(item["origin_logical_path"]), item["sha256"])
            require(row["bytes"] == item["bytes"], "extra file size mismatch")
    sources, helpers = builder.bind_sources(), builder.bind_helpers()
    report = dict(schema="archive-source-version-bindings-v1", source_claim_complete=not builder.gaps, cutoff=args.cutoff,
                  closed_jobs=jobs, after_cutoff_excluded=builder.late, excluded_retained_temporary=builder.excluded_temporary,
                  bindings=sources, helper_bindings=helpers, gaps=builder.gaps,
                  scope="Byte/source closure only; not artifact acceptance or final whole-campaign accounting.")
    builder.write("SOURCE_VERSION_BINDINGS.json", encode(report))
    output = builder.safe(rel(args.output_map))
    require(args.output_map.startswith("coordination/") and not output.exists(), "fresh private mapping path")
    with output.open("xb") as stream: stream.write(encode(dict(schema="explicit-archive-copy-map-v1", files=sorted(builder.rows.values(), key=lambda r: r["destination"]))))
    print("complete" if not builder.gaps else "incomplete", "files", len(builder.rows), "source bindings", len(sources), "gaps", len(builder.gaps))
    return 0 if not builder.gaps else 2

if __name__ == "__main__":
    raise SystemExit(main())
