"""UNRUN candidate: independent artifact/statistics/accounting check, stdlib only.

Does not import the analyzer, scientific packages, or optimization kernels.
ROOT contains the frozen campaign; CHECK_INDEX is a root-relative JSON file.
No study is launched and no state vector is numerically recomputed here.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import tomllib

CASES = ("LiH_r1p487_631g_14m4e", "LiH_r2p113_631g_14m4e",
         "BeH2_r1p369_631g_14m4e", "BeH2_r1p907_631g_14m4e")
ARMS = ("canonical", "magnitude", "current", "first", "pauli", "cubic",
        "window_sparse", "window_magnitude")
WINDOWS = ("window_sparse", "window_magnitude")
CONTROLS = ("canonical", "magnitude", "current", "first", "pauli", "cubic")
PATHS = tuple((a, f) for a in ARMS for f in ((1, 2) if a in ("canonical", "current", "first") else (2,)))
CONDITIONS = tuple(itertools.product((.4, 1.6), (8, 24)))
METRICS = ("projective_error", "infidelity", "norm_drift", "number_leakage", "spin_leakage")
EPS = 1e-12
STATE_EPS = 1e-10
SHA = re.compile(r"[0-9a-f]{64}\Z")
LIMIT_FILE = 32 * 1024**2
LIMIT_BYTES = 512 * 1024**2
LIMIT_FILES = 8192
PRIMARY_SOURCES = {
    "physical.jl": "numerical_audit/prospective_window_ordering/physical.jl",
    "src/WindowPhysical.jl": "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl",
    "src/WindowStudyIO.jl": "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl",
    "../julia_cubic_feasibility/src/ResponseCore.jl": "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl",
}
OBJECTIVE_SOURCES = {
    "objective_check.jl": "numerical_audit/prospective_window_ordering/objective_check.jl",
    "src/WindowStudyIO.jl": PRIMARY_SOURCES["src/WindowStudyIO.jl"],
    "../julia_cubic_feasibility/src/ResponseCore.jl": PRIMARY_SOURCES["../julia_cubic_feasibility/src/ResponseCore.jl"],
    "../scalable_ordering/training.py": "numerical_audit/scalable_ordering/training.py",
    "../../analysis/sparse_response_ordering.py": "analysis/sparse_response_ordering.py",
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def number(value, *, missing=False):
    if missing and value is None:
        return None
    need(type(value) in (int, float), "Expected numeric scalar")
    v = float(value)
    need(math.isfinite(v) and v >= 0, "Nonfinite/negative measured scalar")
    return v


def count(value, low=0, high=10**9):
    need(type(value) is int and low <= value <= high, "Integer domain")
    return value


def json_object(raw):
    def unique(pairs):
        out = {}
        for key, value in pairs:
            need(key not in out, "Duplicate JSON key: " + key)
            out[key] = value
        return out
    def floating(text):
        value = float(text)
        need(math.isfinite(value), "Nonfinite JSON number")
        return value
    value = json.loads(raw, object_pairs_hook=unique, parse_float=floating,
                       parse_constant=lambda x: (_ for _ in ()).throw(ValueError("Nonfinite JSON: " + x)))
    need(isinstance(value, dict), "Expected JSON object")
    return value


class Store:
    """All artifact accesses stay inside the current ROOT, including relocation."""
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.seen = {}
        self.bytes = 0

    def path(self, name, base=None):
        need(isinstance(name, str) and 0 < len(name) <= 2048, "Artifact path length")
        rel = PurePosixPath(name)
        need(not rel.is_absolute() and ".." not in rel.parts and "\\" not in name, "Artifact escape")
        p = ((self.root if base is None else base) / Path(*rel.parts)).resolve()
        need(p.is_relative_to(self.root), "Artifact symlink escape")
        return p

    def get(self, desc, kind="json", base=None, max_size=LIMIT_FILE):
        need(isinstance(desc, dict) and {"path", "sha256"} <= desc.keys(), "Hash-bound descriptor required")
        need(isinstance(desc["sha256"], str) and SHA.fullmatch(desc["sha256"]), "Digest spelling")
        p = self.path(desc["path"], base)
        name = str(p.relative_to(self.root))
        need(p.is_file() and p.stat().st_size <= max_size, "Missing/oversized artifact: " + name)
        if name not in self.seen:
            need(len(self.seen) < LIMIT_FILES, "Artifact count ceiling")
            h = hashlib.sha256()
            n = 0
            with p.open("rb") as stream:
                while block := stream.read(1024**2):
                    n += len(block)
                    self.bytes += len(block)
                    need(n <= max_size and self.bytes <= LIMIT_BYTES, "Hash byte ceiling")
                    h.update(block)
            self.seen[name] = {"sha256": h.hexdigest(), "bytes": n}
        observed = self.seen[name]
        need(observed["sha256"] == desc["sha256"], "Artifact digest differs: " + name)
        need("bytes" not in desc or desc["bytes"] == observed["bytes"], "Artifact byte count differs")
        if kind == "binary":
            return p, None
        with p.open("rb") as stream:
            raw = stream.read(max_size+1)
        need(len(raw) <= max_size and hashlib.sha256(raw).hexdigest() == desc["sha256"], "Artifact changed during read")
        value = json_object(raw) if kind == "json" else tomllib.loads(raw.decode()) if kind == "toml" else raw.decode()
        return p, value

    def child(self, parent, name, digest, kind="binary", max_size=LIMIT_FILE):
        return self.get({"path": name, "sha256": digest}, kind, parent.parent, max_size)


def row_id(row):
    need(row["arm"] in ARMS, "Logical arm")
    count(row["repetition"], 1, 2)
    count(row["formula"], 1, 2)
    count(row["budget_multiplier"], 8, 24)
    need(type(row["time"]) in (int, float) and row["time"] in (.4, 1.6), "Logical physical time")
    return row["arm"], row["repetition"], row["formula"], row["time"], row["budget_multiplier"]


def call_id(row):
    need(row["case"] in CASES and row["arm"] in ARMS, "Constructor identity")
    count(row["repetition"], 1, 2)
    return row["case"], row["arm"], row["repetition"]


def all_rows(case):
    return [(case, a, r, f, t, c) for r in (1, 2) for a, f in PATHS for t, c in CONDITIONS]


def indexed(rows, key, allowed=None):
    need(isinstance(rows, list) and len(rows) <= 4096, "Bounded row collection")
    result = {}
    for row in rows:
        identity = key(row)
        need(identity not in result, "Duplicate logical row: " + repr(identity))
        need(allowed is None or identity in allowed, "Unexpected logical row: " + repr(identity))
        result[identity] = row
    return result


def equal(actual, expected, where):
    """Reporting arithmetic gets only a few binary64 ULPs, never EPS slack."""
    if type(expected) is float and type(actual) in (int, float):
        need(math.isfinite(actual) and abs(actual - expected) <= 4 * max(math.ulp(expected), math.ulp(float(actual))), where)
    else:
        need(type(actual) is type(expected) and actual == expected, where)


def fields(actual, expected, where):
    for key, value in expected.items():
        need(key in actual, where + ": missing " + key)
        equal(actual[key], value, where + ": " + key)


def subtract(left, right):
    # An independent exact-dyadic subtraction followed by one final rounding.
    return float(Fraction(left) - Fraction(right))


def median(values):
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else float((Fraction(ordered[middle-1]) + Fraction(ordered[middle])) / 2)


def outcome(delta):
    return "win" if delta < -EPS else "loss" if delta > EPS else "tie"


def metrics(row):
    result = {key: number(row[key]) for key in METRICS}
    need(result["norm_drift"] <= 1e-9 and all(v <= 1 + 1e-10 for k, v in result.items() if k != "norm_drift"), "Metric domain")
    return result


def raw_measurement(receipt):
    f = receipt.get("final_accounting", {})
    cpu = number(f.get("CPU_seconds"), missing=True)
    wall = number(f.get("elapsed_seconds"), missing=True)
    unit = receipt.get("unit")
    need(isinstance(unit, str) and 0 < len(unit) <= 512, "Owned service identity required")
    required = dict(CPUQuota="25%", CPUWeight="5", IOWeight="10", Nice="19",
                    MemoryHigh="448M", MemoryMax="512M", MemorySwapMax="0",
                    CPUAccounting="yes", MemoryAccounting="yes")
    need(all(receipt.get("settings", {}).get(k) == v for k, v in required.items()), "Frozen owning guard settings")
    complete = (cpu is not None and wall is not None and f.get("CPU_accounting_available") is True
                and f.get("process_CPU_timing_available") is True and receipt.get("guard_error") is None
                and not receipt.get("cleanup_errors"))
    return {"raw_CPU_seconds": cpu, "raw_wall_seconds": wall, "raw_complete": complete, "unit": unit}


def footer_upper(log):
    """Interpret retained systemd CPU footers; no wall/CPU inference from prose."""
    scales = {"h": 3600, "min": 60, "s": 1, "ms": .001, "us": .000001, "µs": .000001, "ns": .000000001}
    upper = []
    for line in log.splitlines():
        if not line.startswith("CPU time consumed: "):
            continue
        text = line.split(": ", 1)[1]
        tokens = re.findall(r"([0-9]+(?:\.[0-9]+)?)(h|min|ms|us|µs|ns|s)", text)
        need(tokens and re.sub(r"([0-9]+(?:\.[0-9]+)?)(h|min|ms|us|µs|ns|s)", "", text).strip() in ("", "."), "Unrecognized CPU footer")
        upper.append(sum(float(v) * scales[u] for v, u in tokens) + .001)
    return max(upper) if upper else None


def final_charge(plan, receipt, final, log):
    measured = raw_measurement(receipt)
    cap = number(plan["attempt_CPU_seconds"])
    candidates = [cap if measured["raw_CPU_seconds"] is None else measured["raw_CPU_seconds"]]
    footer = footer_upper(log)
    if footer is not None:
        candidates.append(footer)
    service = final.get("final_service_snapshot", {}).get("CPUUsageNSec")
    service_cpu = int(service) / 1e9 if isinstance(service, str) and service.isdigit() else None
    if service_cpu is not None:
        candidates.append(service_cpu)
    charged = number(final["charged_CPU_seconds"])
    need(abs(charged - max(candidates)) <= 1e-9, "FINAL charge does not reconcile guard/footer/service")
    wall = number(final["absolute_guarded_wall_seconds"])
    need(abs(wall - number(receipt["elapsed"])) <= 1e-9, "FINAL absolute wall differs from receipt")
    if measured["raw_wall_seconds"] is not None:
        need(abs(wall - measured["raw_wall_seconds"]) <= 1e-9, "Raw/final wall boundary differs")
    return dict(measured, charged_CPU_seconds=charged, charged_wall_seconds=wall,
                charge_is_reserved_fallback=measured["raw_CPU_seconds"] is None,
                footer_CPU_upper=footer, final_service_CPU_seconds=service_cpu,
                charge_minus_raw_CPU=None if measured["raw_CPU_seconds"] is None else charged-measured["raw_CPU_seconds"],
                conservative_cap_exceeded=charged > cap or wall > number(plan["attempt_wall_seconds"]),
                raw_cap_exceeded=None if not measured["raw_complete"] else
                measured["raw_CPU_seconds"] > cap or measured["raw_wall_seconds"] > number(plan["attempt_wall_seconds"]))


def phase_accounting(store, spec):
    """Explicit frozen ledger cutoff excludes this check's not-yet-finished job."""
    _, pools = store.get(spec["resource_pools"])
    _, freeze = store.get(spec["original_campaign_freeze"])
    need(freeze["resource_pools_sha256"] == spec["resource_pools"]["sha256"], "Cumulative resource pool was reset")
    store.get(spec["guard_source"], "binary")
    need(pools["guard"]["sha256"] == spec["guard_source"]["sha256"], "Guard source binding")
    dispatchers = set()
    for d in spec["dispatcher_sources"]:
        store.get(d, "binary")
        dispatchers.add(d["sha256"])
    cutoff = datetime.fromisoformat(spec["ledger_cutoff_UTC"])
    need(cutoff.tzinfo is not None, "Timezone-aware ledger cutoff required")
    rows = spec["jobs"]
    need(isinstance(rows, list) and len(rows) <= 2048, "Bounded job ledger")
    by_label, by_receipt, units = {}, {}, set()
    for row in rows:
        pp, p = store.get(row["plan"])
        label = p["label"]
        need(isinstance(label, str) and re.fullmatch(r"[A-Za-z0-9-]+", label), "Job label")
        need(pp == store.path("jobs/" + label + ".PLAN.json") and label not in by_label, "Unique label/PLAN path")
        started = datetime.fromisoformat(p["started_UTC"])
        need(started.tzinfo is not None and started <= cutoff, "Job outside frozen ledger cutoff")
        need(p["pools_sha256"] == spec["resource_pools"]["sha256"] and p["dispatcher_sha256"] in dispatchers, "Job pool/dispatcher binding")
        need(p["phase"] in pools["phases"], "Unknown accounting phase")
        fp, f = store.get(row["final"])
        rp, r = store.get(row["receipt"])
        lp, log = store.get(row["log"], "text")
        need(fp == store.path("jobs/"+label+".FINAL.json") and rp == store.path("logs/"+label+".resources.json")
             and lp == store.path("logs/"+label+".log"), "Job label/artifact paths")
        need(f["label"] == label and f["phase"] == p["phase"] and f["receipt_sha256"] == row["receipt"]["sha256"], "PLAN/FINAL/receipt join")
        need(r["command"] == p["command"] and r["budget"]["cpu_seconds"] == p["attempt_CPU_seconds"]
             and r["budget"]["wall_seconds"] == p["attempt_wall_seconds"], "Executed allocation differs from PLAN")
        need(f["exit_code"] == r["exit_code"], "FINAL/receipt exit differs")
        charge = final_charge(p, r, f, log)
        need(charge["unit"] not in units and row["receipt"]["sha256"] not in by_receipt, "Duplicate owned service/receipt in phase ledger")
        units.add(charge["unit"])
        entry = dict(charge, label=label, phase=p["phase"], plan=p, final=f, receipt=r, descriptors=row)
        by_label[label] = entry
        by_receipt[row["receipt"]["sha256"]] = entry
    # Independently discover omitted historical PLANs up to cutoff, including failures.
    present = list(store.path("jobs").glob("*.PLAN.json"))
    need(len(present) <= 2048, "On-disk PLAN inventory bound")
    before = set()
    for path in present:
        need(path.resolve().is_relative_to(store.root) and path.stat().st_size <= LIMIT_FILE, "PLAN path/size")
        with path.open("rb") as stream:
            raw = stream.read(LIMIT_FILE+1)
        need(len(raw) <= LIMIT_FILE, "PLAN grew past byte bound")
        body = json_object(raw)
        started = datetime.fromisoformat(body["started_UTC"])
        need(started.tzinfo is not None, "PLAN timestamp timezone")
        if started <= cutoff:
            before.add(body["label"])
    need(before == set(by_label), "Ledger omits/adds retained jobs at cutoff")
    totals = {}
    for phase, ceiling in pools["phases"].items():
        group = sorted((e for e in by_label.values() if e["phase"] == phase), key=lambda e: datetime.fromisoformat(e["plan"]["started_UTC"]))
        cpu = wall = 0.
        for e in group:
            need(abs(number(e["plan"]["prior_phase_CPU_seconds"])-cpu) <= 1e-8 and
                 abs(number(e["plan"]["prior_phase_wall_seconds"])-wall) <= 1e-8, "Phase ledger reset/omitted attempt")
            need(cpu + number(e["plan"]["attempt_CPU_seconds"]) <= number(ceiling["CPU_seconds"]) + 1e-8 and
                 wall + number(e["plan"]["attempt_wall_seconds"]) <= number(ceiling["wall_seconds"]) + 1e-8, "Attempt not admitted within cumulative pool")
            cpu += e["charged_CPU_seconds"]
            wall += e["charged_wall_seconds"]
        complete = all(e["raw_complete"] for e in group)
        totals[phase] = dict(attempts=len(group), raw_accounting_complete=complete,
            raw_CPU_seconds=sum(e["raw_CPU_seconds"] for e in group) if complete else None,
            raw_wall_seconds=sum(e["raw_wall_seconds"] for e in group) if complete else None,
            known_raw_CPU_seconds=sum(e["raw_CPU_seconds"] for e in group if e["raw_CPU_seconds"] is not None),
            known_raw_wall_seconds=sum(e["raw_wall_seconds"] for e in group if e["raw_wall_seconds"] is not None),
            conservative_charged_CPU_seconds=cpu, conservative_charged_wall_seconds=wall,
            reserved_fallback_attempts=sum(e["charge_is_reserved_fallback"] for e in group),
            phase_cap_exceeded=cpu > ceiling["CPU_seconds"] or wall > ceiling["wall_seconds"],
            raw_unknown_labels=[e["label"] for e in group if not e["raw_complete"]])
    return by_label, by_receipt, totals


def source_inventory(store, index, spec):
    pp, protocol = store.get(index["input_protocol"])
    need(re.fullmatch(r"[0-9a-f]{40}", protocol.get("source_commit", "")), "Full scientific source commit")
    need(protocol["source_commit"] == spec["expected_scientific_source_commit"], "Wrong scientific source target")
    source = protocol["sources"]
    need(isinstance(source, dict) and 0 < len(source) <= 4096, "Source map size")
    for name, digest in source.items():
        store.get({"path": "source_snapshot/"+name, "sha256": digest}, "binary", pp.parent)
    registry = protocol["registry"]
    registry_name = "numerical_audit/prospective_window_ordering/REGISTRY.json"
    _, snapshot_registry = store.get({"path": "source_snapshot/"+registry_name,
                                      "sha256": source[registry_name]}, base=pp.parent)
    need(snapshot_registry == registry, "Protocol registry differs from bound source registry")
    need({s["case"] for s in registry["cases"]} == set(CASES) and len(registry["cases"]) == 4, "Original four-case roster")
    need(all(s["n"] == 14 and s["electrons"] == 4 for s in registry["cases"]), "Model register/electron contract")
    config = protocol["configuration"]
    need(config["training_time"] == 1 and config["training_cap_times_original_words"] == 30
         and config["fresh_repetitions"] == 2 and tuple(config["arms"]) == ARMS, "Construction configuration")
    p = config["physical"]
    need(p["times"] == [.4, 1.6] and p["cap_times_original_words"] == [8, 24]
         and p["win_tie_ambiguity"] == EPS, "Fixed physical conditions/floor")
    need(p["formulas"] == {a: [f for arm, f in PATHS if arm == a] for a in ARMS}, "Eleven physical paths")
    ip, inventory = store.get(index["input_freeze"])
    need(ip.parent == pp.parent and inventory["protocol_sha256"] == index["input_protocol"]["sha256"]
         and inventory["expected_cases"] == 4, "Input/protocol binding")
    need(set(inventory["successes"]) | set(inventory["refusals"]) == set(CASES)
         and not set(inventory["successes"]) & set(inventory["refusals"]), "Complete input status inventory")
    models = {}
    for case, hashes in inventory["successes"].items():
        need(set(hashes) == {"model.json", "model.toml", "tensors.npz", "orbitals.npz"}, "Four chemistry payloads required")
        base = ip.parent / "inputs" / case
        jp, j = store.get({"path": "model.json", "sha256": hashes["model.json"]}, base=base)
        tp, t = store.get({"path": "model.toml", "sha256": hashes["model.toml"]}, "toml", base)
        need(j["case"] == t["case"] == case and j["n"] == t["n"] == 14
             and j["electrons"] == t["electrons"] == 4, "Model identity")
        specification = next(s for s in registry["cases"] if s["case"] == case)
        need(all(j.get(k) == v for k, v in specification.items()) and
             all(j.get(k) == registry[k] for k in ("basis", "unit", "charge", "spin")), "Model/declared scientific instance differs")
        need(t["source_json_sha256"] == hashes["model.json"] and t["partition"] == "conserving", "Model serialization binding")
        for file, field in (("tensors.npz", "tensors_sha256"), ("orbitals.npz", "orbitals_sha256")):
            need(j[field] == hashes[file], "Chemistry payload link")
            store.get({"path": file, "sha256": hashes[file]}, "binary", base)
        need(0 < len(j["coefficients"]) <= 8192 and 0 < len(j["blocks"]) <= 1024, "Input size")
        lengths = [len(block) for block in j["blocks"]]
        need(sum(lengths) == len(j["coefficients"]) and len(t["units"]) == len(lengths), "Original word partition counts")
        for block, unit in zip(j["blocks"], t["units"]):
            pairs = sorted(block, key=lambda item: tuple((q, letter) for q, letter in item[0]))
            need(unit["words"] == [",".join(str(q)+letter for q, letter in word) for word, _ in pairs]
                 and [int(x) for x in unit["numerators"]] == [int(v) for _, v in pairs], "JSON/TOML word/coefficient mapping")
        need(int(j["denominator"]) == int(t["denominator"]) > 0, "Common exact denominator")
        models[case] = dict(hashes=hashes, lengths=lengths, words=sum(lengths), json_path=str(jp.relative_to(store.root)),
                            toml_path=str(tp.relative_to(store.root)))
    reuse = spec.get("input_reuse")
    if reuse is not None:
        _, side = store.get(reuse["sidecar"])
        _, old_protocol = store.get(reuse["original_protocol"])
        _, old_inventory = store.get(reuse["original_inventory"])
        need(side["chemistry_rerun"] is False and side["new_protocol_sha256"] == index["input_protocol"]["sha256"]
             and side["original_protocol_sha256"] == reuse["original_protocol"]["sha256"]
             and side["original_inventory_sha256"] == reuse["original_inventory"]["sha256"]
             and old_inventory["protocol_sha256"] == reuse["original_protocol"]["sha256"], "Chemistry reuse graph")
        need(side["original_generation_source_commit"] == old_protocol["source_commit"]
             and side["new_execution_source_commit"] == protocol["source_commit"]
             and side["original_payload_hashes"] == old_inventory["successes"] == inventory["successes"], "Producer versus execution identity/reused bytes")
    return protocol, inventory, models


def source_hashes(record, mapping, protocol):
    claimed = record.get("source_hashes")
    need(isinstance(claimed, dict) and set(claimed) == set(mapping), "Complete executing source inventory")
    for logical, repository in mapping.items():
        digest = claimed[logical]
        need(isinstance(digest, str) and SHA.fullmatch(digest) and protocol["sources"].get(repository) == digest,
             "Executing dependency differs: " + logical)


def authenticated_plan(store, index, protocol, models, jobs, receipts):
    _, plan = store.get(index["plan"])
    expected = set(itertools.product(CASES, ARMS, (1, 2)))
    need(plan["schema"] == "prospective-window-constructor-plan-v1" and plan["mode"] == "main" and plan["expected_calls"] == 64 and plan["source_snapshot"] == protocol["sources"]
         and plan["source_freeze"] == index["input_protocol"] and plan["input_freeze"] == index["input_freeze"]
         and plan["physical_error_read"] is False, "Main plan/source binding")
    entries = indexed(plan["entries"], call_id, expected)
    need(set(entries) == expected, "All 64 planned calls required")
    _, whole = store.get(index["all_orders_freeze"])
    need(whole["schema"] == "prospective-window-all-orders-freeze-v1" and whole["mode"] == "main" and whole["expected_calls"] == 64 and whole["inventory_complete"] is True
         and whole["physical_error_read"] is False and whole["no_selected_repetition"] is True
         and whole["plan_sha256"] == index["plan"]["sha256"] and whole["source_freeze"] == index["input_protocol"], "Whole-order freeze binding")
    frozen = indexed(whole["artifacts"], call_id, expected)
    need(set(frozen) == expected, "All 64 attempted-order entries required")
    records = {}
    used = set()
    for identity, entry in entries.items():
        need(entry["attempt_CPU_seconds"] == 64 and entry["attempt_wall_seconds"] == 320, "Original 64 CPU/320 wall ceilings changed")
        need(entry["label"] in jobs, "Unjoined constructor label")
        job = jobs[entry["label"]]
        need(job["phase"] == "constructors" and job["plan"]["command"] == entry["command"]
             and job["plan"]["attempt_CPU_seconds"] == 64 and job["plan"]["attempt_wall_seconds"] == 320,
             "Main planned/executed command or caps differ")
        artifact = frozen[identity]
        need(artifact["receipt"] == job["descriptors"]["receipt"] and artifact["receipt"]["path"] == entry["receipt"]
             and artifact["exit_code"] == job["receipt"]["exit_code"], "Main final receipt label/hash join")
        digest = artifact["receipt"]["sha256"]
        need(digest in receipts and digest not in used, "Constructor receipt duplication")
        used.add(digest)
        record = None
        if artifact.get("record") is not None:
            need(artifact["record"]["path"] == entry["result"], "Planned record path differs")
            _, record = store.get(artifact["record"])
            need(call_id(record) == identity and record["schema"] == "prospective-window-constructor-v1"
                 and record["artifact_complete"] is True and record["physical_error_read"] is False
                 and record["source_freeze_sha256"] == index["input_protocol"]["sha256"]
                 and identity[0] in models and record["model_sha256"] == models[identity[0]]["hashes"]["model.json"], "Constructor source/target binding")
        records[identity] = dict(job=job, artifact=artifact, record=record)
    need(whole["all_calls_have_artifacts"] is all(x["record"] is not None for x in records.values()), "Whole-order artifact presence claim")
    return records


def cost_inventory(store, inventory, panel, case, receipts, global_calls):
    expected = set(itertools.product(ARMS, (1, 2)))
    ledger = indexed(inventory["constructors"], lambda r: (r["arm"], r["repetition"]), expected)
    reported = indexed(panel["constructors"], lambda r: (r["arm"], r["repetition"]), expected)
    need(set(ledger) == set(reported) == expected, "Sixteen constructor ledgers per case")
    costs = {}
    for key in expected:
        attempts = ledger[key]["attempts"]
        p = reported[key]
        need(isinstance(attempts, list) and len(attempts) <= 32 and len(attempts) == len(p["attempts"]), "Attempt-ledger bound/link")
        jobs = []
        known_cpu = known_wall = 0.
        within = True
        for attempt, display in zip(attempts, p["attempts"]):
            need(attempt["receipt"] == display["receipt"] and attempt.get("record") == display.get("constructor"), "PANEL receipt/record link")
            digest = attempt["receipt"]["sha256"]
            need(digest in receipts, "Attempt omitted from cumulative ledger")
            job = receipts[digest]
            need(job["descriptors"]["receipt"] == attempt["receipt"], "Receipt descriptor alias")
            need(job["plan"]["attempt_CPU_seconds"] <= 64-known_cpu and job["plan"]["attempt_wall_seconds"] <= 320-known_wall, "Retry extends original whole-call caps")
            within = within and job["raw_complete"] and job["raw_cap_exceeded"] is False
            known_cpu += job["raw_CPU_seconds"] or 0.
            known_wall += job["raw_wall_seconds"] or 0.
            jobs.append(job)
        complete = bool(jobs) and all(j["raw_complete"] for j in jobs)
        cpu = known_cpu if complete else None
        wall = known_wall if complete else None
        need(p.get("guard_CPU_seconds") == cpu and p.get("guard_wall_seconds") == wall, "PANEL comparison measurement differs")
        if complete:
            need(p.get("cost_complete") is True and p["within_total_construction_ceiling"] is bool(within and cpu <= 64 and wall <= 320), "Raw guard ceiling classification differs")
        final = global_calls[(case, *key)]
        need(attempts and attempts[-1]["receipt"] == final["artifact"]["receipt"]
             and attempts[-1].get("record") == final["artifact"].get("record"), "Whole/case final attempt differs")
        if p.get("status") == "retained_order":
            need(p["selected_record"] == final["artifact"]["record"] and
                 p["selected_attempt"] == ledger[key]["selected_attempt"] == len(attempts), "Selected repetition/attempt differs")
            need(final["record"] is not None and p["constructor"]["order"] == final["record"]["constructor"]["order"]
                 and p["constructor"].get("score") == final["record"]["constructor"].get("score"), "Selected constructor body differs")
        charged_cpu = sum(j["charged_CPU_seconds"] for j in jobs)
        charged_wall = sum(j["charged_wall_seconds"] for j in jobs)
        costs[(case, *key)] = dict(cpu=cpu, wall=wall, known_cpu=known_cpu, known_wall=known_wall, complete=complete,
            within_raw=complete and within and cpu <= 64 and wall <= 320, reported=p, attempts=jobs,
            charged_cpu=charged_cpu, charged_wall=charged_wall,
            conservative_cap_exceeded=charged_cpu > 64 or charged_wall > 320 or any(j["conservative_cap_exceeded"] for j in jobs))
    return costs


NUMERIC_OBJECTIVE_STATUSES = frozenset(("verified", "mismatch", "floor_censored"))


def objective_tolerances(check):
    relative = number(check.get("relative_tolerance"))
    absolute = number(check.get("absolute_tolerance"))
    need(0 < relative <= 1e-8 and 0 < absolute <= 1e-24, "Objective tolerance domain")
    return relative, absolute


def objective_numeric_acceptance(record, relative, absolute):
    """Reconstruct compare_scores from serialized binary64 operands, not flags.

    Exact dyadic arithmetic rounds at each ordinary Julia binary64 operation.
    Serialized comparison fields must match exactly; no physical EPS/ULP slack
    can move an objective across its tolerance or Boolean floor boundary.
    A refusal supplies no numerical comparison and remains unaccepted.
    """
    relative, absolute = objective_tolerances(dict(relative_tolerance=relative,
                                                  absolute_tolerance=absolute))
    status = record.get("status")
    need(isinstance(status, str) and bool(status), "Missing objective status")
    if status not in NUMERIC_OBJECTIVE_STATUSES:
        for key in ("verified", "within_tolerance", "numerical_floor"):
            if key in record:
                need(record[key] is False, "Refused objective carries acceptance/floor flag")
        for key in ("constructor_score", "replayed_score", "absolute_difference",
                    "relative_difference", "tolerance"):
            if key in record:
                number(record[key])
        return False
    actual, claimed = number(record.get("replayed_score")), number(record.get("constructor_score"))
    maximum = max(actual, claimed)
    difference = abs(subtract(actual, claimed))
    scaled = float(Fraction(relative) * Fraction(maximum))
    tolerance = float(Fraction(absolute) + Fraction(scaled))
    relative_difference = float(Fraction(difference) / Fraction(maximum)) if maximum else 0.
    for key, expected in (("absolute_difference", difference), ("relative_difference", relative_difference),
                          ("tolerance", tolerance)):
        need(number(record.get(key)) == expected, "Objective comparison field differs: " + key)
    within, floor = difference <= tolerance, maximum <= absolute
    verified = within and not floor
    expected_status = "mismatch" if not within else "floor_censored" if floor else "verified"
    need(status == expected_status, "Objective status differs from numeric comparison")
    for key, expected in (("within_tolerance", within), ("numerical_floor", floor), ("verified", verified)):
        need(record.get(key) is expected, "Objective Boolean field differs: " + key)
    return verified


def objective_acceptance(check, candidates):
    relative, absolute = objective_tolerances(check)
    logical = indexed(check.get("logical_rows", []), row_id, set(candidates))
    numeric = check.get("checks", [])
    need(isinstance(numeric, list) and len(numeric) <= 256, "Objective check bound")
    # Validate even unused/deduplicated numerical records; none may conceal
    # forged flags or a malformed numeric status behind logical-row selection.
    accepted_checks = [objective_numeric_acceptance(q, relative, absolute) for q in numeric]
    positions = {key: i+1 for i, key in enumerate(candidates)}
    result = {}
    for key, row in logical.items():
        c = candidates[key]
        need(row["name"] == c["name"] and row["candidate_index"] == positions[key], "Objective logical position differs")
        if not c.get("constructor_scored", False):
            need("constructor_score" not in c and row["status"] == "unscored_censored" and "check_index" not in row, "Unscored objective must remain explicit")
            result[key] = (True, "unscored_no_numeric_claim")
        elif "check_index" not in row:
            need(row["status"] not in NUMERIC_OBJECTIVE_STATUSES | {"unscored_censored"}, "Missing scored objective check")
            result[key] = (False, row["status"])
        else:
            i = count(row["check_index"], 1, len(numeric))
            q = numeric[i-1]
            need(q["check_index"] == i and q["order"] == c["order"] and q["partition"] == c["partition"]
                 and q["score_kind"] == c["score_kind"] and q["constructor_score"] == c["constructor_score"]
                 and q["status"] == row["status"], "Objective record target differs")
            result[key] = (accepted_checks[i-1], q["status"])
    return result


def independent_metric_acceptance(primary, checked, hard):
    """Check S6 scalar witnesses; never reconstruct or import state evolution."""
    saved = metrics(checked["metrics"])
    independent = metrics(checked["independent_reference_metrics"])
    for field, compared in (("saved_reference_metric_discrepancies", saved),
                            ("independent_reference_metric_discrepancies", independent)):
        for metric in METRICS:
            difference = abs(subtract(primary[metric], compared[metric]))
            need(difference <= EPS, "Independent metric discrepancy: " + field + ": " + metric)
            need(number(checked[field].get(metric)) == difference, "Independent stored discrepancy differs: " + field + ": " + metric)
    need(checked["independent_unfused_replay"] is hard, "Independent replay scope")
    if hard:
        need(number(checked["unfused_state_norm_residual"]) <= STATE_EPS, "Unfused state residual")
        unfused_saved = metrics(checked["unfused_saved_reference_metrics"])
        unfused_independent = metrics(checked["unfused_independent_reference_metrics"])
        for metric in METRICS:
            need(abs(subtract(unfused_saved[metric], saved[metric])) <= EPS,
                 "Unfused saved-reference metric discrepancy: " + metric)
            difference = abs(subtract(primary[metric], unfused_independent[metric]))
            need(difference <= EPS, "Unfused independent-reference metric discrepancy: " + metric)
            need(number(checked["unfused_independent_reference_metric_discrepancies"].get(metric)) == difference,
                 "Unfused stored discrepancy differs: " + metric)


def physical_affordability(model, candidate, row):
    lengths = [1]*model["words"] if candidate["partition"] == "pauli" else model["lengths"]
    order = candidate["order"]
    need(all(type(i) is int for i in order) and sorted(order) == list(range(1, len(lengths)+1)), "Original complete label permutation")
    budget = row["budget_multiplier"] * model["words"]
    seq = [lengths[i-1] for i in order]
    total = sum(seq)
    if len(seq) == 1:
        steps, pulses = (1, total) if budget >= total else (0, 0)
    elif row["formula"] == 1:
        steps = budget // total
        pulses = steps*total
    else:
        steps = (budget-seq[0]) // (2*total-seq[0]-seq[-1])
        pulses = steps*(2*total-seq[0]-seq[-1])+seq[0]
    need(steps > 0 and row["rotation_budget"] == budget and row["steps"] == steps and row["ideal_rotations"] == pulses <= budget, "Physical step/pulse affordability")


def case_graph(store, entry, index, protocol, model, receipts, global_calls):
    """Authenticate raw provenance before exposing any primary metric to derivation."""
    case = entry["case"]
    ip, inventory = store.get(entry["inventory"])
    need(inventory["schema"] == "prospective-window-panel-inventory-v1" and ip.parent == store.root and inventory["case"] == case and inventory["mode"] == "main"
         and inventory["source_freeze"] == index["input_protocol"] and inventory["input_freeze"] == index["input_freeze"]
         and inventory["all_orders_freeze"] == index["all_orders_freeze"], "Case inventory/source graph")
    need(inventory["model_json"] == {"path": model["json_path"], "sha256": model["hashes"]["model.json"]}
         and inventory["model_toml"] == {"path": model["toml_path"], "sha256": model["hashes"]["model.toml"]}, "Case target descriptors")
    pp, panel = store.get(entry["panel"])
    need(panel["schema"] == "prospective-window-panel-audit-v1" and panel["case"] == case and panel["mode"] == "main" and panel["inventory_sha256"] == entry["inventory"]["sha256"]
         and panel["source_freeze"] == index["input_protocol"] and panel["source_commit"] == protocol["source_commit"]
         and panel["model_json"] == inventory["model_json"] and panel["model_toml"] == inventory["model_toml"], "Panel target/source graph")
    costs = cost_inventory(store, inventory, panel, case, receipts, global_calls)
    records = {}
    if entry.get("candidates") is None:
        need(panel.get("candidate_panel_complete") is not True, "Complete panel lacks candidates")
        return costs, records, False
    cp, c = store.get(entry["candidates"], "toml")
    need(store.path(panel["candidates"]["path"], pp.parent) == cp and panel["candidates"]["sha256"] == entry["candidates"]["sha256"], "Panel candidate descriptor")
    need(panel["status"] == "complete" and panel["candidate_panel_complete"] is True
         and c["schema"] == "prospective-window-candidates-v1" and c["status"] == "complete"
         and c["case"] == case and c["mode"] == "main" and c["n"] == 14 and c["hf"] == "15"
         and c["source_commit"] == protocol["source_commit"] and c["source_freeze_sha256"] == index["input_protocol"]["sha256"]
         and c["input_freeze_sha256"] == index["input_freeze"]["sha256"]
         and c["all_orders_freeze_sha256"] == index["all_orders_freeze"]["sha256"]
         and c["inventory_sha256"] == entry["inventory"]["sha256"] and c["expected_logical_rows"] == 88
         and c["input_sha256"] == model["hashes"]["model.toml"] and c["source_json_sha256"] == model["hashes"]["model.json"], "Candidate source/target graph")
    expected = {key[1:] for key in all_rows(case)}
    candidates = indexed(c["candidates"], row_id, expected)
    need(set(candidates) == expected, "Eighty-eight candidate identities required")
    for key, row in candidates.items():
        cost = costs[(case, key[0], key[1])]
        ctor = cost["reported"].get("constructor")
        need(cost["complete"] and cost["within_raw"] and cost["reported"]["status"] == "retained_order"
             and row["constructor_record_sha256"] == cost["reported"]["selected_record"]["sha256"]
             and row["order"] == [i+1 for i in ctor["order"]]
             and row["partition"] == ("pauli" if key[0] == "pauli" else "conserving"), "Candidate constructor/partition binding")
        need(row.get("constructor_scored", False) == (ctor.get("score") is not None), "Candidate scored status")
        if row.get("constructor_scored"):
            need(row["constructor_score"] == ctor["score"], "Candidate constructor score differs")
    objectives = {}
    if entry.get("objective_check") is not None:
        _, obj = store.get(entry["objective_check"], "toml")
        need(obj["input_sha256"] == model["hashes"]["model.toml"] and obj["candidates_sha256"] == entry["candidates"]["sha256"]
             and obj["case"] == case and obj["n"] == 14 and obj["hf"] == "15" and obj["training_time"] == 1
             and obj["training_budget_multiplier"] == 30, "Objective source/target")
        source_hashes(obj, OBJECTIVE_SOURCES, protocol)
        objectives = objective_acceptance(obj, candidates)
    if entry.get("physical") is None or entry.get("independent_check") is None:
        return costs, records, False
    primary_path, primary = store.get(entry["physical"], "toml")
    source_hashes(primary, PRIMARY_SOURCES, protocol)
    need(primary["case"] == case and primary["mode"] == "main" and primary["n"] == 14 and primary["hf"] == "15"
         and primary["expected_logical_rows"] == 88 and primary["input_sha256"] == model["hashes"]["model.toml"]
         and primary["source_json_sha256"] == model["hashes"]["model.json"] and primary["orders_sha256"] == entry["candidates"]["sha256"], "Primary source/target/order link")
    _, independent = store.get(entry["independent_check"], "toml")
    need(independent["case"] == case and independent["n"] == 14 and independent["hf"] == "15"
         and independent["manifest_sha256"] == entry["physical"]["sha256"]
         and independent["input_toml_sha256"] == model["hashes"]["model.toml"]
         and independent["source_json_sha256"] == model["hashes"]["model.json"]
         and independent["checker_source_sha256"] == protocol["sources"]["numerical_audit/prospective_window_ordering/independent.jl"]
         and number(independent["state_norm_tolerance"]) <= STATE_EPS and number(independent["absolute_metric_tolerance"]) <= EPS,
         "Independent source/target/precision link")
    primary_rows = indexed(primary["rows"], row_id, expected)
    checked_rows = indexed(independent.get("rows", []), row_id, set(primary_rows))
    references = set()
    for ref in independent.get("reference_checks", []):
        if ref["status"] == "complete":
            need(number(ref["state_norm_residual"]) <= STATE_EPS and number(ref["phase_aligned_residual"]) <= STATE_EPS, "Reference residual")
            references.add((ref["time"], ref["reference_file"], ref["reference_sha256"]))
    for key, p in primary_rows.items():
        candidate = candidates[key]
        need(p["name"] == candidate["name"] and p["order"] == candidate["order"] and p["partition"] == candidate["partition"], "Primary row/order identity")
        physical_affordability(model, candidate, p)
        values = metrics(p)
        for prefix in ("actual", "reference"):
            state, _ = store.child(primary_path, p[prefix+"_file"], p[prefix+"_sha256"], max_size=16*16384)
            need(store.seen[str(state.relative_to(store.root))]["bytes"] == 16*16384, "Full-space state byte length")
        q = checked_rows.get(key)
        if q is None or q.get("status") != "complete":
            continue
        need(q["name"] == p["name"] and q["saved_actual_sha256"] == p["actual_sha256"]
             and (key[3], p["reference_file"], p["reference_sha256"]) in references, "Independent state/reference binding")
        hard = key[3:] == (1.6, 24)
        independent_metric_acceptance(values, q, hard)
        need(number(p["reference_norm_drift"]) <= 1e-9 and abs(number(q["metrics"]["reference_norm_drift"])-p["reference_norm_drift"]) <= EPS, "Reference norm drift")
        accepted, status = objectives.get(key, (False, None))
        records[(case, *key)] = dict(values=values, validated=accepted, objective_status=status,
                                      candidate=candidate, primary=p, independent=q)
    complete = (set(primary_rows) == set(checked_rows) == expected and len(records) == 88
                and primary.get("status") == "complete" and primary.get("physical_panel_complete") is True
                and independent.get("status") == "complete" and not independent.get("failures")
                and independent.get("checked_logical_rows") == 88 and set(objectives) == expected
                and all(v[0] for v in objectives.values()) and
                (obj.get("complete") is True if any(c.get("constructor_scored") for c in candidates.values())
                 else obj.get("status") == "no_scored_orders_supplied"))
    return costs, records, complete


def paired_record(window, baseline, cpu, wall):
    """One comparison, with raw metrics independent of objective acceptance."""
    out = dict(classification="censored", error_delta=None, infidelity_delta=None,
               CPU_delta=cpu, wall_delta=wall, baseline_over_window_error_ratio=None,
               fully_validated=False)
    if window is None or baseline is None:
        out["ratio_suppression_reason"] = "Unknown/incomplete physical metric"
        return out
    w, b = window["values"], baseline["values"]
    delta = subtract(w["projective_error"], b["projective_error"])
    valid = window["validated"] and baseline["validated"]
    out.update(classification=outcome(delta), error_delta=delta,
               infidelity_delta=subtract(w["infidelity"], b["infidelity"]), fully_validated=valid,
               window_error=w["projective_error"], baseline_error=b["projective_error"],
               window_floor=w["projective_error"] <= EPS, baseline_floor=b["projective_error"] <= EPS)
    if min(w["projective_error"], b["projective_error"]) <= EPS:
        reason = "At or below absolute physical floor/ambiguity"
    elif not valid:
        reason = "Objective check censored or unverified"
    else:
        ratio = b["projective_error"]/w["projective_error"]
        reason = None if math.isfinite(ratio) else "Ratio overflows"
        if reason is None:
            out["baseline_over_window_error_ratio"] = ratio
    out["ratio_suppression_reason"] = reason
    return out


def derivations(records, costs):
    pairs, medians, summaries = {}, {}, {}
    for window, baseline in itertools.product(WINDOWS, CONTROLS):
        observations = []
        for case, rep, condition in itertools.product(CASES, (1, 2), CONDITIONS):
            t, cap = condition
            wkey, bkey = (case, window, rep, 2, t, cap), (case, baseline, rep, 2, t, cap)
            wc, bc = costs.get(wkey[:3]), costs.get(bkey[:3])
            deltas = [None if wc is None or bc is None or wc[axis] is None or bc[axis] is None
                      else subtract(wc[axis], bc[axis]) for axis in ("cpu", "wall")]
            record = paired_record(records.get(wkey), records.get(bkey), *deltas)
            key = window, baseline, case, rep, 2, t, cap
            pairs[key] = record
            observations.append((key, record))
        complete_groups = 0
        for case, (t, cap) in itertools.product(CASES, CONDITIONS):
            group = [(rep, pairs[(window, baseline, case, rep, 2, t, cap)]) for rep in (1, 2)]
            available = [(rep, p) for rep, p in group if p["error_delta"] is not None]
            validated = [p for _, p in group if p["fully_validated"]]
            full = len(available) == 2
            complete_groups += full
            observed_median = median([p["error_delta"] for _, p in available])
            value = dict(expected_repetitions=2, available_repetitions=[r for r, _ in available],
                complete_repetitions=full, median_paired_error_delta=observed_median if full else None,
                observed_only_median_paired_error_delta=observed_median,
                fully_validated_median_paired_error_delta=median([p["error_delta"] for p in validated]) if len(validated) == 2 else None)
            for axis in ("CPU_delta", "wall_delta"):
                values = [p[axis] for _, p in group if p[axis] is not None]
                value["median_paired_"+axis] = median(values) if len(values) == 2 else None
            medians[(window, baseline, case, 2, t, cap)] = value
        categories = Counter(p["classification"] for _, p in observations)
        known = [p for _, p in observations if p["error_delta"] is not None]
        summaries[(window, baseline, 2)] = dict(expected_pairs=32, physically_checked_pairs=len(known),
            fully_validated_pairs=sum(p["fully_validated"] for _, p in observations),
            **{k: categories[k] for k in ("win", "tie", "loss", "censored")},
            complete_repetition_groups=complete_groups, expected_repetition_groups=16,
            complete_panel_median_of_paired_deltas=median([p["error_delta"] for p in known]) if len(known) == 32 else None,
            ratios_retained=sum(p["baseline_over_window_error_ratio"] is not None for _, p in observations))
    return pairs, medians, summaries


def compare_derived(summary, records, costs):
    pairs, medians, comparisons = derivations(records, costs)
    pair_key = lambda r: (r["window"], r["baseline"], r["case"], r["repetition"], r["formula"], r["time"], r["budget_multiplier"])
    median_key = lambda r: (r["window"], r["baseline"], r["case"], r["formula"], r["time"], r["budget_multiplier"])
    summary_key = lambda r: (r["window"], r["baseline"], r["formula"])
    for name, expected, key in (("paired_observations", pairs, pair_key),
                               ("paired_repetition_medians", medians, median_key), ("comparisons", comparisons, summary_key)):
        provided = indexed(summary[name], key, set(expected))
        need(set(provided) == set(expected), "Complete derived inventory: " + name)
        for identity, fields_expected in expected.items():
            fields(provided[identity], fields_expected, name+repr(identity))
    return dict(paired_observations=len(pairs), paired_repetition_medians=len(medians), win_tie_loss_summaries=len(comparisons))


def compare_raw(summary, records, costs):
    expected = {key for case in CASES for key in all_rows(case)}
    supplied = indexed(summary["raw_rows"], lambda r: (r["case"], *row_id(r)), expected)
    need(set(supplied) == expected, "All 352 raw logical rows required")
    for key in expected:
        reported = supplied[key]
        actual = records.get(key)
        cost = costs.get(key[:3])
        fields(reported, {"physical_checked": actual is not None,
                         "fully_validated": bool(actual and actual["validated"]),
                         "constructor_CPU_seconds": None if cost is None else cost["cpu"],
                         "constructor_wall_seconds": None if cost is None else cost["wall"]}, "Raw row "+repr(key))
        if actual is None:
            need(reported["status"] == "censored" and not reported.get("metrics"), "Unknown row invented metrics/status")
        else:
            fields(reported, {"status": "physical_checked", "physical_floor": actual["values"]["projective_error"] <= EPS,
                              "objective_accepted": actual["validated"]}, "Raw verified row "+repr(key))
            fields(reported["metrics"], actual["values"], "Raw metrics "+repr(key))
            # Exact embedded source rows prevent output relabeling or dropped fields.
            need(reported["physical"] == actual["primary"] and reported["independent"] == actual["independent"]
                 and reported["candidate"] == actual["candidate"], "Summary embeds a different raw source row")
            if actual["objective_status"] is None:
                need(reported.get("objective") is None, "Missing objective check was invented")
            else:
                need(reported["objective"]["status"] == actual["objective_status"], "Objective status differs")
    return len(expected)


def compare_pareto(summary, records, costs):
    keys = set(itertools.product(CASES, (1, 2), (.4, 1.6), (8, 24)))
    provided = indexed(summary["pareto_memberships"], lambda r: (r["case"], r["repetition"], r["time"], r["budget_multiplier"]), keys)
    need(set(provided) == keys, "Thirty-two Pareto condition groups")
    for case, rep, t, cap in keys:
        target = provided[(case, rep, t, cap)]
        need(target["expected_paths"] == 11 and set(target["sets"]) == {"CPU_seconds", "wall_seconds"}, "Full eleven-path/two-axis Pareto contract")
        for axis, column in (("CPU_seconds", "cpu"), ("wall_seconds", "wall")):
            points, excluded = {}, set()
            for arm, formula in PATHS:
                record = records.get((case, arm, rep, formula, t, cap))
                cost = costs.get((case, arm, rep))
                name = arm+"_f"+str(formula)
                if record is None or cost is None or cost[column] is None:
                    excluded.add(name)
                else:
                    points[name] = dict(path=name, cost=cost[column], projective_error=record["values"]["projective_error"],
                        physical_floor=record["values"]["projective_error"] <= EPS,
                        fully_validated=record["validated"], objective_status=record["objective_status"])
            # Independently form the directed dominance graph for the full path set.
            edges = {p: set() for p in points}
            for left, right in itertools.permutations(points, 2):
                a, b = points[left], points[right]
                if dominates(a["cost"], a["projective_error"], b["cost"], b["projective_error"]):
                    edges[right].add(left)
            supplied = target["sets"][axis]
            ps = indexed(supplied["points"], lambda p: p["path"], set(points))
            need(set(ps) == set(points), "Pareto eligible points")
            for path, value in points.items():
                fields(ps[path], value, "Pareto point "+path)
            need(set(supplied["dominated_by"]) == set(edges), "Pareto dominance domain")
            need(all(set(supplied["dominated_by"][p]) == edges[p] and len(supplied["dominated_by"][p]) == len(edges[p]) for p in edges), "Pareto dominance edges")
            frontier = {p for p in edges if not edges[p]}
            need(set(supplied["pareto_paths"]) == frontier and len(supplied["pareto_paths"]) == len(frontier), "Pareto frontier")
            ex = indexed(supplied["excluded"], lambda p: p["path"], excluded)
            need(set(ex) == excluded, "Pareto exclusion inventory")
    return 32


def dominates(left_cost, left_error, right_cost, right_error):
    return (left_cost <= right_cost and left_error <= right_error + EPS and
            (left_cost < right_cost or left_error < right_error - EPS))


def compare_costs(summary, costs, receipts, index):
    supplied = indexed(summary["constructor_costs"], call_id, set(itertools.product(CASES, ARMS, (1, 2))))
    need(set(supplied) == set(costs), "Constructor cost inventory differs")
    used_receipts, used_units = set(), set()
    for key, cost in costs.items():
        fields(supplied[key], dict(CPU_seconds=cost["cpu"], wall_seconds=cost["wall"], accounting_complete=cost["complete"],
                                  known_CPU_seconds=cost["known_cpu"], known_wall_seconds=cost["known_wall"]), "Constructor cost "+repr(key))
        expected_receipts = {j["descriptors"]["receipt"]["sha256"] for j in cost["attempts"]}
        actual_receipts = {a["descriptor"]["sha256"] for a in supplied[key]["attempts"]}
        need(expected_receipts == actual_receipts and len(actual_receipts) == len(supplied[key]["attempts"]), "Summary cost attempts")
        for job in cost["attempts"]:
            digest, unit = job["descriptors"]["receipt"]["sha256"], job["unit"]
            need(digest not in used_receipts and unit not in used_units, "A constructor service charged more than once")
            used_receipts.add(digest)
            used_units.add(unit)
    def totals(group, expected):
        known = [c for c in group if c["complete"]]
        return dict(expected_calls=expected, known_calls=len(known),
            CPU_seconds=sum(c["cpu"] for c in known) if len(known) == expected else None,
            wall_seconds=sum(c["wall"] for c in known) if len(known) == expected else None,
            known_CPU_seconds=sum(c["known_cpu"] for c in group), known_wall_seconds=sum(c["known_wall"] for c in group))
    fields(summary["construction_totals"], totals(list(costs.values()), 64), "Construction totals")
    need(set(summary["construction_by_arm"]) == set(ARMS), "By-arm cost inventory")
    for arm in ARMS:
        group = [c for key, c in costs.items() if key[1] == arm]
        fields(summary["construction_by_arm"][arm], dict(totals(group, 8), indexed_calls=len(group)), "Arm costs "+arm)
    # Analyzer comparison phases exclude receipts already counted in construction;
    # the full cumulative phase ledger above includes them once in their owner phase.
    phase = defaultdict(list)
    hashes, units = set(), {}
    for entry in index.get("phase_receipts", []):
        digest = entry["receipt"]["sha256"]
        need(digest in receipts and digest not in hashes, "Unknown/duplicate summary phase receipt")
        hashes.add(digest)
        job = receipts[digest]
        need(job["unit"] not in units, "Repeated summary phase service")
        units[job["unit"]] = digest
        need(entry["phase"] == job["phase"], "Phase label differs from dispatch owner")
        phase[entry["phase"]].append(job)
    need(set(summary["other_phase_totals"]) == set(phase), "Summary phase inventory")
    for name, rows in phase.items():
        unique = [j for j in rows if j["descriptors"]["receipt"]["sha256"] not in used_receipts]
        complete = all(j["raw_complete"] for j in unique)
        fields(summary["other_phase_totals"][name], dict(indexed_attempts=len(unique), quarantined_attempts=0,
            accounting_complete=complete, CPU_seconds=sum(j["raw_CPU_seconds"] for j in unique) if complete else None,
            wall_seconds=sum(j["raw_wall_seconds"] for j in unique) if complete else None,
            known_CPU_seconds=sum(j["raw_CPU_seconds"] for j in unique if j["raw_CPU_seconds"] is not None),
            known_wall_seconds=sum(j["raw_wall_seconds"] for j in unique if j["raw_wall_seconds"] is not None)), "Raw comparison phase "+name)
    return used_receipts


def run(root, check_index):
    store = Store(root)
    path = store.path(check_index)
    need(path.is_file() and path.stat().st_size <= LIMIT_FILE, "Check-index size")
    with path.open("rb") as stream:
        raw = stream.read(LIMIT_FILE+1)
    need(len(raw) <= LIMIT_FILE, "Check index grew past bound")
    spec = json_object(raw)
    need(spec.get("schema") == "prospective-independent-derived-check-index-v1", "Independent index schema")
    result = dict(schema="prospective-independent-derived-check-v1", status="incomplete_or_failed",
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        check_index_sha256=hashlib.sha256(raw).hexdigest(), executed_scientific_kernels=False,
        ambiguity=EPS, final_source_head=spec.get("expected_scientific_source_commit"),
        failures=[], censored_cases=[], source_coverage_gaps=[
            "Does not rerun chemistry, objective kernels, propagation, Krylov or state-vector metrics; consumes separately bound independent checks.",
            "Hash/source identity is not a Git server signature or live CI/approval check.",
            "Ledger totals stop at the explicit pre-check cutoff; this check and later jobs need subsequent accounting."],
        complete_352_derived_claim_allowed=False, unqualified_under_ceiling_claim_allowed=False)
    try:
        _, index = store.get(spec["analysis_index"])
        _, summary = store.get(spec["analyzer_summary"])
        need(index["schema"] == "prospective-window-analysis-index-v1" and index["mode"] == "main", "Main analysis index required")
        need(summary["schema"] == "prospective-window-report-data-v1" and summary["mode"] == "main"
             and summary["analysis_index_sha256"] == spec["analysis_index"]["sha256"]
             and summary["ambiguity"] == summary["ratio_floor"] == EPS, "Analyzer summary input/precision binding")
        store.get(spec["analyzer_source"], "binary")
        need(summary["analyzer"]["sha256"] == spec["analyzer_source"]["sha256"], "Compared analyzer executable identity")
        jobs, receipts, phase_totals = phase_accounting(store, spec)
        result["phase_totals"] = phase_totals
        result["ledger_cutoff_UTC"] = spec["ledger_cutoff_UTC"]
        result["accounted_jobs"] = [{k: v for k, v in e.items() if k not in ("plan", "final", "receipt")} for e in jobs.values()]
        protocol, inventory, models = source_inventory(store, index, spec)
        global_calls = authenticated_plan(store, index, protocol, models, jobs, receipts)
        entries = indexed(index.get("cases", []), lambda r: r["case"], set(CASES))
        records, costs, completeness = {}, {}, {}
        for case in CASES:
            if case not in entries or case not in models:
                result["censored_cases"].append(dict(case=case, reason="Missing case graph or authenticated input"))
                completeness[case] = False
                continue
            try:
                cc, rr, full = case_graph(store, entries[case], index, protocol, models[case], receipts, global_calls)
                costs.update(cc)
                records.update(rr)
                completeness[case] = full
                if not full:
                    result["censored_cases"].append(dict(case=case, reason="Incomplete objective/primary/independent graph", checked_rows=len(rr)))
            except (ValueError, KeyError, TypeError, OSError, OverflowError) as error:
                completeness[case] = False
                result["censored_cases"].append(dict(case=case, reason=str(error), invalid_graph=True))
        # Recover declared partial costs without turning unknown total work into zero.
        for key, frozen in global_calls.items():
            if key in costs:
                continue
            matched = [frozen["job"]]
            for item in index.get("partial_attempts", []):
                if call_id(item) != key:
                    continue
                job = receipts.get(item["receipt"]["sha256"])
                need(job is not None, "Unaccounted partial attempt")
                if job not in matched:
                    matched.append(job)
            costs[key] = dict(cpu=None, wall=None, complete=False, within_raw=False,
                known_cpu=sum(j["raw_CPU_seconds"] or 0. for j in matched), known_wall=sum(j["raw_wall_seconds"] or 0. for j in matched),
                attempts=matched, charged_cpu=sum(j["charged_CPU_seconds"] for j in matched),
                charged_wall=sum(j["charged_wall_seconds"] for j in matched),
                conservative_cap_exceeded=any(j["conservative_cap_exceeded"] for j in matched), reported={})
        result["raw_rows_checked"] = compare_raw(summary, records, costs)
        result.update(compare_derived(summary, records, costs))
        result["pareto_groups_checked"] = compare_pareto(summary, records, costs)
        compare_costs(summary, costs, receipts, index)
        result["constructor_reconciliation"] = [dict(case=k[0], arm=k[1], repetition=k[2],
            raw_guard_CPU_seconds=c["cpu"], raw_guard_wall_seconds=c["wall"],
            known_raw_guard_CPU_seconds=c["known_cpu"], known_raw_guard_wall_seconds=c["known_wall"],
            final_conservative_CPU_seconds=c["charged_cpu"], final_conservative_wall_seconds=c["charged_wall"],
            raw_guard_accounting_complete=c["complete"], final_charge_crosses_call_ceiling=c["conservative_cap_exceeded"],
            labels=[j["label"] for j in c["attempts"]]) for k, c in sorted(costs.items())]
        boundary = [row for row in result["constructor_reconciliation"] if row["final_charge_crosses_call_ceiling"]]
        result["boundary_exceedances"] = boundary
        complete = (inventory["status"] == "complete" and len(models) == 4 and all(completeness.values())
                    and len(records) == 352 and all(r["validated"] for r in records.values())
                    and len(costs) == 64 and all(c["complete"] for c in costs.values()))
        fields(summary["completion"], dict(expected_logical_rows=352, expected_constructor_calls=64,
            authenticated_independent_physical_rows=len(records), fully_validated_rows=sum(r["validated"] for r in records.values()),
            missing_or_censored_rows=352-len(records)), "Summary completion")
        if summary["completion"].get("complete_352_main_claim_allowed") is True:
            need(complete and not summary["failures"], "Analyzer overclaims complete main")
        result["derived_values_agree"] = True
        result["complete_352_derived_claim_allowed"] = complete and not summary["failures"]
        result["unqualified_under_ceiling_claim_allowed"] = (complete and not boundary and
            all(not t["phase_cap_exceeded"] and t["raw_accounting_complete"] for t in phase_totals.values()))
        result["status"] = "verified_complete_derived_panel" if result["complete_352_derived_claim_allowed"] else "verified_censored_derivations"
        if boundary or any(t["phase_cap_exceeded"] for t in phase_totals.values()):
            result["status"] = "derived_agreement_with_accounting_veto"
    except (ValueError, KeyError, TypeError, OSError, OverflowError, RecursionError) as error:
        result["failures"].append(dict(type=type(error).__name__, message=str(error)))
    result["artifact_bindings"] = store.seen
    result["bytes_hashed"] = store.bytes
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("check_index", help="ROOT-relative immutable check-index JSON")
    parser.add_argument("output", type=Path, help="New output; never replaces an existing result")
    args = parser.parse_args(argv)
    result = run(args.root, args.check_index)
    raw = (json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()
    need(len(raw) <= 32*1024**2, "Output bound")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as stream:
        stream.write(raw)
    return 0 if result["status"] == "verified_complete_derived_panel" and result["unqualified_under_ceiling_claim_allowed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
