"""UNRUN: real full-checker baseline and two narrowly targeted data mutations."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import shutil
import signal
import subprocess
import sys
import tempfile

CHECKER = "tools/independent_derived_check.py"
CHECKER_SHA = "4ca02076449b5a1f0e49c187e01c425b60330a5d0320895728fc68bab23fc41d"
INDEX = "INDEPENDENT_DERIVED_INDEX.main.H.a1.json"
REPORT = "INDEPENDENT_DERIVED_CHECK.main.H.a1.json"
MAX_FILE, MAX_BYTES, MAX_FILES = 32*1024**2, 512*1024**2, 8192


def need(ok, message):
    if not ok: raise ValueError(message)


def path(root, name):
    rel = PurePosixPath(name)
    need(isinstance(name, str) and 0 < len(name) <= 2048 and not rel.is_absolute() and
         ".." not in rel.parts and "\\" not in name, "Confined relative path required")
    result = root.joinpath(*rel.parts).resolve()
    need(result.is_relative_to(root), "Symlink escape")
    return result


def raw(file):
    need(file.is_file() and file.stat().st_size <= MAX_FILE, "Missing/oversize artifact")
    with file.open("rb") as stream: value = stream.read(MAX_FILE+1)
    need(len(value) <= MAX_FILE, "Artifact grew past bound")
    return value


def sha(file):
    return hashlib.sha256(raw(file)).hexdigest()


def read(file):
    return json.loads(raw(file))


def encoded(value):
    value = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()
    need(len(value) <= MAX_FILE, "JSON output bound")
    return value


def fresh(file, value):
    with file.open("xb") as stream: stream.write(value)


def disk_bytes(root):
    total = count = 0
    for item in root.rglob("*"):
        if item.is_file():
            count += 1; total += item.stat().st_size
            need(count <= 65536, "Campaign disk-scan file bound")
    return total  # Match guard: each hardlink pathname is charged again.


def replace_link(mirror, root, name, value):
    target = path(mirror, name); original = path(root, name)
    need(os.path.samefile(target, original), "Mutation anchor is not its cloned original")
    temporary = target.with_name(target.name+".exclusive-overlay")
    fresh(temporary, value); os.replace(temporary, target)  # Never open linked inode for writing.
    need(not os.path.samefile(target, original) and raw(target) == value, "Distinct overlay inode/readback required")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path); parser.add_argument("--report-sha256", required=True)
    args = parser.parse_args(); root = args.root.resolve()
    need(sha(path(root, REPORT)) == args.report_sha256 and sha(path(root, CHECKER)) == CHECKER_SHA, "Pinned report/checker differs")
    report = read(path(root, REPORT)); spec = read(path(root, INDEX))
    need(report["check_index_sha256"] == sha(path(root, INDEX)) and report["checker_sha256"] == CHECKER_SHA and
         report["status"] == "verified_complete_derived_panel" and report["derived_values_agree"] is True and
         report["complete_352_derived_claim_allowed"] is True and report["unqualified_under_ceiling_claim_allowed"] is True and
         not report["failures"] and not report["censored_cases"] and not report["boundary_exceedances"], "Actual complete baseline required")
    counts = dict(raw_rows_checked=352, paired_observations=384, paired_repetition_medians=192,
                  win_tie_loss_summaries=12, pareto_groups_checked=32)
    need(all(report[k] == v for k, v in counts.items()) and len(spec["jobs"]) == 184 and
         len(report["constructor_reconciliation"]) == 64, "Actual baseline inventory differs")
    bindings = copy.deepcopy(report["artifact_bindings"])
    for name in (INDEX, REPORT, CHECKER): bindings[name] = dict(sha256=sha(path(root, name)), bytes=path(root, name).stat().st_size)
    need(len(bindings) <= MAX_FILES and sum(v["bytes"] for v in bindings.values()) <= MAX_BYTES, "Mirror inventory bound")
    before = {name: sha(path(root, name)) for name in bindings}
    need(all(before[n] == v["sha256"] and path(root, n).stat().st_size == v["bytes"] for n, v in bindings.items()), "Original graph differs")
    analysis = read(path(root, spec["analysis_index"]["path"]))
    plan = read(path(root, analysis["plan"]["path"]))
    need(len(plan["entries"]) == 64 and plan["mode"] == "main", "Original main64 plan required")
    final_name = "jobs/"+plan["entries"][0]["label"]+".FINAL.json"
    final = read(path(root, final_name))
    need(type(final["charged_CPU_seconds"]) in (int, float) and final["charged_CPU_seconds"] > 0, "Nonzero real FINAL anchor required")
    summary_name = spec["analyzer_summary"]["path"]; summary = read(path(root, summary_name))
    first = summary["paired_observations"][0]; delta = first["error_delta"]
    need(type(delta) in (int, float) and math.isfinite(delta) and delta+1 != delta, "Numeric real paired anchor required")
    pair_key = tuple(first[k] for k in ("window", "baseline", "case", "repetition", "formula", "time", "budget_multiplier"))
    expected_pair = "paired_observations"+repr(pair_key)+": error_delta"
    out = path(root, "main-negative-guards.a1"); need(not out.exists(), "Fresh retained output directory required")
    tmp = path(root, "tmp"); tmp.mkdir(exist_ok=True)
    need(not any(tmp.glob("main-negative-guard-*")), "Prior owned fixture requires explicit cleanup")
    cloned = {name: value for name, value in bindings.items() if name != REPORT}
    logical_mirror = sum(v["bytes"] for v in cloned.values())
    # Bound retained baseline/mutant reports, overlays, indexes and per-mirror replacement copies.
    report_bound = min(MAX_FILE, len(encoded(report))+65536)
    reserve = 3*report_bound + 5*(len(encoded(spec))+65536) + 3*(len(encoded(summary))+len(encoded(final))+65536) + len(encoded(bindings)) + 1024**2
    ceiling = read(path(root, "RESOURCE_POOLS.json"))["output_MiB"]*1024**2
    need(ceiling == 256*1024**2 and disk_bytes(root)+logical_mirror+reserve <= ceiling, "Logical ROOT disk cap cannot admit mirror and retained evidence")
    out.mkdir(); results = []; completed = False
    fresh(out/"ORIGINAL_BINDINGS.json", encoded(bindings))
    try:
        for mode in ("baseline", "final-charge", "paired-delta"):
            need(disk_bytes(root)+logical_mirror+reserve <= ceiling, "Remaining logical disk budget insufficient")
            mirror = Path(tempfile.mkdtemp(prefix="main-negative-guard-", dir=tmp)); inode = mirror.stat().st_ino
            try:
                for name in cloned:
                    target = path(mirror, name); target.parent.mkdir(parents=True, exist_ok=True)
                    os.link(path(root, name), target)
                changed = copy.deepcopy(spec); overlay_name = None
                if mode == "final-charge":
                    value = copy.deepcopy(final); value["charged_CPU_seconds"] = 0
                    overlay_name = final_name; payload = encoded(value)
                    matches = [j for j in changed["jobs"] if j["final"]["path"] == final_name]
                    need(len(matches) == 1, "Exactly one canonical FINAL descriptor required")
                    descriptor = matches[0]["final"]
                elif mode == "paired-delta":
                    value = copy.deepcopy(summary); value["paired_observations"][0]["error_delta"] = delta+1
                    overlay_name = summary_name; payload = encoded(value); descriptor = changed["analyzer_summary"]
                if overlay_name is not None:
                    replace_link(mirror, root, overlay_name, payload)
                    fresh(out/(mode+".OVERLAY.json"), payload)
                    descriptor["sha256"] = hashlib.sha256(payload).hexdigest()
                    if "bytes" in descriptor: descriptor["bytes"] = len(payload)
                    index_bytes = encoded(changed); replace_link(mirror, root, INDEX, index_bytes)
                    fresh(out/(mode+".INDEX.json"), index_bytes)
                target_report = out/(mode+".CHECK.json")
                env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
                run = subprocess.run([sys.executable, str(path(mirror, CHECKER)), str(mirror), INDEX, str(target_report)],
                                     cwd=mirror, env=env, capture_output=True, timeout=25)
                need(target_report.is_file() and target_report.stat().st_size <= report_bound, "Bounded checker report required")
                checked = read(target_report)
                need(checked["checker_sha256"] == CHECKER_SHA and checked["check_index_sha256"] == sha(path(mirror, INDEX)), "Executed checker/index binding differs")
                if mode == "baseline":
                    need(run.returncode == 0 and checked["status"] == "verified_complete_derived_panel" and not checked["failures"] and
                         checked["derived_values_agree"] is True and checked["complete_352_derived_claim_allowed"] is True and
                         checked["unqualified_under_ceiling_claim_allowed"] is True and all(checked[k] == v for k, v in counts.items()), "Fresh mirror baseline failed")
                    anchor = "verified_complete_derived_panel"
                else:
                    anchor = "FINAL charge does not reconcile guard/footer/service" if mode == "final-charge" else expected_pair
                    need(run.returncode == 2 and checked["status"] == "incomplete_or_failed" and
                         checked["failures"] == [dict(type="ValueError", message=anchor)], "Mutation failed for a different reason")
                results.append(dict(mode=mode, exit_code=run.returncode, exact_outcome=anchor,
                    overlay_original_path=overlay_name, check_report_sha256=sha(target_report),
                    index_sha256=sha(path(mirror, INDEX)), overlay_sha256=None if overlay_name is None else sha(path(mirror, overlay_name)),
                    mirrored_files=len(cloned), baseline_counts=counts if mode == "baseline" else None))
            finally:
                need(mirror.parent == tmp and mirror.is_dir() and mirror.stat().st_ino == inode, "Owned temporary directory identity changed")
                shutil.rmtree(mirror)
        completed = True
    finally:
        unchanged = all(path(root, n).is_file() and sha(path(root, n)) == digest for n, digest in before.items())
        fresh(out/"OUTCOME.json", encoded(dict(status="passed" if completed and unchanged and len(results) == 3 else "incomplete_or_failed",
            original_bindings_unchanged=unchanged, results=results, checker_sha256=CHECKER_SHA,
            original_report_sha256=args.report_sha256, driver_sha256=sha(Path(__file__)), numerical_kernels_executed=False)))
        need(unchanged, "An original graph artifact changed")
    print(json.dumps(dict(status="main_guard_proofs_passed", baseline_passes=1, targeted_rejections=2, originals_unchanged=True)))


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(RuntimeError("Guard termination")))
    try: main()
    except Exception:
        print(json.dumps(dict(status="main_guard_proof_incomplete", detail="Retain bounded reports/outcome; no quality values printed")))
        raise SystemExit(2)
