"""UNRUN source-only candidate: immutable main analysis index, then check index.

No scientific imports or execution. Hashes artifacts without interpreting
physical metrics; the separate analyzer/checker determine outcomes.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path, PurePosixPath
import re

G = "722013aa6d9a0a48776e26de5f356b5a441c3f61"
ANALYZER_SHA = "a18db4dbf2adf85ccc1bee1a9f17aafd87f04bc337d8b5e06b5eea11ef442bc9"
CHECKER_SHA = "4ca02076449b5a1f0e49c187e01c425b60330a5d0320895728fc68bab23fc41d"
CASES = ("LiH_r1p487_631g_14m4e", "LiH_r2p113_631g_14m4e",
         "BeH2_r1p369_631g_14m4e", "BeH2_r1p907_631g_14m4e")
ARMS = ("canonical", "magnitude", "current", "first", "pauli", "cubic",
        "window_sparse", "window_magnitude")
LIMIT = 32 * 1024**2


def need(ok, message):
    if not ok:
        raise ValueError(message)


def timestamp(value):
    result = datetime.fromisoformat(value)
    need(result.tzinfo is not None, "Timezone-aware cutoff/job timestamp required")
    return result


def decode(raw):
    def pairs(items):
        out = {}
        for key, value in items:
            need(key not in out, "Duplicate JSON key")
            out[key] = value
        return out
    def constant(value):
        raise ValueError("Nonfinite JSON constant")
    result = json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)
    need(isinstance(result, dict), "Metadata object required")
    return result


class Files:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.seen = {}
        self.total = 0

    def path(self, name):
        need(isinstance(name, str) and 0 < len(name) <= 2048, "Bounded relative path required")
        rel = PurePosixPath(name)
        need(not rel.is_absolute() and ".." not in rel.parts and "\\" not in name, "Path escape")
        path = self.root.joinpath(*rel.parts).resolve()
        need(path.is_relative_to(self.root), "Symlink escape")
        return path

    def desc(self, name, optional=False):
        path = self.path(name)
        if optional and not path.exists():
            return None
        need(path.is_file() and path.stat().st_size <= LIMIT, "Missing/oversize artifact")
        if name not in self.seen:
            need(len(self.seen) < 8192, "Artifact-count ceiling")
            h = hashlib.sha256(); size = 0
            with path.open("rb") as stream:
                while block := stream.read(1024**2):
                    size += len(block); self.total += len(block)
                    need(size <= LIMIT and self.total <= 512*1024**2, "Hash byte ceiling")
                    h.update(block)
            self.seen[name] = h.hexdigest()
        return dict(path=name, sha256=self.seen[name])

    def bound(self, descriptor):
        need(isinstance(descriptor, dict) and set(descriptor) == {"path", "sha256"}, "Exact path/hash descriptor required")
        need(self.desc(descriptor["path"]) == descriptor, "Artifact hash differs")
        return descriptor

    def read(self, name):
        descriptor = self.desc(name)
        with self.path(name).open("rb") as stream:
            raw = stream.read(LIMIT+1)
        need(len(raw) <= LIMIT and hashlib.sha256(raw).hexdigest() == descriptor["sha256"], "Metadata changed during read")
        return decode(raw)

    def fresh(self, name):
        path = self.path(name)
        need(not path.exists() and not path.with_name(path.name+".partial").exists(), "Append-only output required")
        return path

    def write(self, name, value):
        path = self.fresh(name)
        raw = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()
        need(len(raw) <= LIMIT, "Output byte ceiling")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(raw)


def collect_ledger(files, cutoff, pools_desc, dispatchers):
    plans = sorted(files.path("jobs").glob("*.PLAN.json"))
    need(len(plans) <= 2048, "Ledger job ceiling")
    jobs, phases, labels, units, receipts = [], [], set(), set(), set()
    pools = files.read(pools_desc["path"])
    for path in plans:
        name = "jobs/" + path.name
        plan = files.read(name)
        started = timestamp(plan["started_UTC"])
        if started > cutoff:
            continue
        label = plan["label"]
        need(isinstance(label, str) and re.fullmatch(r"[A-Za-z0-9-]+", label) and
             name == "jobs/"+label+".PLAN.json" and label not in labels, "Unique PLAN label/path required")
        need(plan["pools_sha256"] == pools_desc["sha256"] and plan["dispatcher_sha256"] in dispatchers and
             plan["phase"] in pools["phases"], "Historical pool/dispatcher/phase binding differs")
        row = dict(plan=files.desc(name), final=files.desc("jobs/"+label+".FINAL.json"),
                   receipt=files.desc("logs/"+label+".resources.json"), log=files.desc("logs/"+label+".log"))
        final = files.read(row["final"]["path"]); receipt = files.read(row["receipt"]["path"])
        need(started <= timestamp(final["completed_UTC"]) <= cutoff, "Job must be terminal by frozen cutoff")
        need(final["label"] == label and final["phase"] == plan["phase"] and
             final["receipt_sha256"] == row["receipt"]["sha256"] and
             final["exit_code"] == receipt["exit_code"], "PLAN/FINAL/receipt join differs")
        need(receipt["command"] == plan["command"] and receipt["budget"]["cpu_seconds"] == plan["attempt_CPU_seconds"] and
             receipt["budget"]["wall_seconds"] == plan["attempt_wall_seconds"], "Executed allocation/command differs")
        unit = receipt.get("unit")
        need(isinstance(unit, str) and unit and unit not in units and row["receipt"]["sha256"] not in receipts,
             "Missing/duplicate owned unit or receipt")
        labels.add(label); units.add(unit); receipts.add(row["receipt"]["sha256"])
        jobs.append(row); phases.append(dict(phase=plan["phase"], receipt=row["receipt"]))
    return jobs, phases


def prepare(args):
    files = Files(args.root)
    files.fresh(args.analysis_index); files.fresh(args.seed)
    need(args.analysis_index != args.seed, "Distinct append-only outputs required")
    cutoff = timestamp(args.cutoff)
    need(cutoff <= datetime.now(timezone.utc), "Future cutoff prohibited")
    protocol_desc = files.desc("data-a6/INPUT_PROTOCOL_FROZEN.json")
    protocol = files.read(protocol_desc["path"])
    need(protocol["source_commit"] == G and len(protocol["sources"]) == 43, "Frozen G source required")
    for name, digest in protocol["sources"].items():
        need(files.desc("data-a6/source_snapshot/"+name)["sha256"] == digest, "Frozen scientific source differs")
    input_desc = files.desc("data-a6/INPUTS_FROZEN.json"); inputs = files.read(input_desc["path"])
    need(inputs["protocol_sha256"] == protocol_desc["sha256"] and
         set(inputs["successes"]) | set(inputs["refusals"]) == set(CASES), "Input inventory binding differs")
    for case, payloads in inputs["successes"].items():
        need(set(payloads) == {"model.json", "model.toml", "tensors.npz", "orbitals.npz"}, "Four input payloads required")
        for name, digest in payloads.items():
            need(files.desc("data-a6/inputs/"+case+"/"+name)["sha256"] == digest, "Input payload differs")
    reuse = dict(sidecar=files.desc("data-a6/INPUT_REUSE_PROVENANCE.json"),
                 original_protocol=files.desc("data-a2/INPUT_PROTOCOL_FROZEN.json"),
                 original_inventory=files.desc("data-a2/INPUTS_FROZEN.json"))
    side = files.read(reuse["sidecar"]["path"]); old = files.read(reuse["original_inventory"]["path"])
    old_protocol = files.read(reuse["original_protocol"]["path"])
    need(side["chemistry_rerun"] is False and side["new_execution_source_commit"] == G and
         side["new_protocol_sha256"] == protocol_desc["sha256"] and
         side["original_protocol_sha256"] == old["protocol_sha256"] == reuse["original_protocol"]["sha256"] and
         side["original_inventory_sha256"] == reuse["original_inventory"]["sha256"] and
         side["original_generation_source_commit"] == old_protocol["source_commit"] and
         side["original_payload_hashes"] == old["successes"] == inputs["successes"], "Original producer/reuse edge differs")
    pools = files.desc("RESOURCE_POOLS.json"); freeze = files.desc("CAMPAIGN_FREEZE_A.json")
    need(files.read(freeze["path"])["resource_pools_sha256"] == pools["sha256"], "Original pool reset prohibited")
    guard = files.desc(args.guard_source)
    need(files.read(pools["path"])["guard"]["sha256"] == guard["sha256"], "Archived guard differs")
    dispatcher_descs = [files.desc(name) for name in args.dispatcher_source]
    need(0 < len(dispatcher_descs) <= 64 and len({d["sha256"] for d in dispatcher_descs}) == len(dispatcher_descs),
         "Bounded unique dispatcher source versions required")
    dispatcher_hashes = {d["sha256"] for d in dispatcher_descs}
    jobs, phases = collect_ledger(files, cutoff, pools, dispatcher_hashes)
    analyzer = files.desc(args.analyzer_source)
    need(analyzer["sha256"] == ANALYZER_SHA, "Archived current analyzer differs")
    checker = files.desc(args.checker_source)
    need(checker["sha256"] == CHECKER_SHA, "Archived independent checker differs")
    plan_desc = files.desc("main-plan.a1.json"); plan = files.read(plan_desc["path"])
    whole_desc = files.desc("ALL_ORDERS_FROZEN.main.a1.json"); whole = files.read(whole_desc["path"])
    identity = lambda row: (row["case"], row["arm"], row["repetition"])
    expected = set(itertools.product(CASES, ARMS, (1, 2)))
    need(plan["schema"] == "prospective-window-constructor-plan-v1" and plan["mode"] == "main" and
         plan["attempt_tag"] == "a1" and plan["expected_calls"] == 64 and
         plan["source_freeze"] == protocol_desc and plan["input_freeze"] == input_desc and
         plan["source_snapshot"] == protocol["sources"] and plan["physical_error_read"] is False and
         len(plan["entries"]) == 64 and {identity(e) for e in plan["entries"]} == expected,
         "All 64 G main planned calls required")
    need(whole["schema"] == "prospective-window-all-orders-freeze-v1" and whole["mode"] == "main" and
         whole["inventory_complete"] is True and whole["expected_calls"] == 64 and
         whole["source_freeze"] == protocol_desc and whole["plan_sha256"] == plan_desc["sha256"] and
         whole["physical_error_read"] is False and whole["no_selected_repetition"] is True and
         len(whole["artifacts"]) == 64 and {identity(e) for e in whole["artifacts"]} == expected,
         "All 64 attempted-order freeze entries required; successes not assumed")
    by_receipt = {r["receipt"]["sha256"]: r for r in jobs}
    frozen = {identity(e): e for e in whole["artifacts"]}
    for entry in plan["entries"]:
        artifact = frozen[identity(entry)]; job = by_receipt.get(artifact["receipt"]["sha256"])
        need(job is not None and artifact["receipt"] == job["receipt"] and
             job["plan"]["path"] == "jobs/"+entry["label"]+".PLAN.json" and
             entry["attempt_CPU_seconds"] == 64 and entry["attempt_wall_seconds"] == 320,
             "Main call missing at cutoff or original caps differ")
        executed = files.read(job["plan"]["path"])
        need(executed["phase"] == "constructors" and executed["command"] == entry["command"] and
             executed["attempt_CPU_seconds"] == 64 and executed["attempt_wall_seconds"] == 320 and
             artifact["receipt"]["path"] == entry["receipt"], "Main executed source/allocation differs")
        if artifact.get("record") is not None:
            files.bound(artifact["record"])
    cases, partial, missing_panels = [], [], []
    for case in CASES:
        inventory_desc = files.desc("main-"+case+".PANEL_INVENTORY.json")
        inventory = files.read(inventory_desc["path"])
        need(inventory["case"] == case and inventory["mode"] == "main" and
             inventory["source_freeze"] == protocol_desc and inventory["input_freeze"] == input_desc and
             inventory["all_orders_freeze"] == whole_desc, "Case inventory/freeze link differs")
        need(len(inventory["constructors"]) == 16 and
             {(r["arm"], r["repetition"]) for r in inventory["constructors"]} == set(itertools.product(ARMS, (1, 2))),
             "All 16 case call ledgers required")
        for call in inventory["constructors"]:
            key = (case, call["arm"], call["repetition"])
            need(key in expected, "Foreign call in case inventory")
            for attempt in call["attempts"]:
                receipt = files.bound(attempt["receipt"])
                need(receipt["sha256"] in by_receipt and by_receipt[receipt["sha256"]]["receipt"] == receipt,
                     "Case attempt absent from cutoff ledger")
                if attempt.get("record") is not None:
                    files.bound(attempt["record"])
                if receipt != frozen[key]["receipt"]:
                    partial.append(dict(case=case, arm=key[1], repetition=key[2],
                                        receipt=receipt, record=attempt.get("record")))
        panel = files.desc("main-panel/"+case+"/PANEL.json", optional=True)
        if panel is None:
            missing_panels.append(case)
            continue  # Readers retain the roster and recover frozen/partial costs.
        cases.append(dict(case=case, inventory=inventory_desc, panel=panel,
            candidates=files.desc("main-panel/"+case+"/CANDIDATES.toml", optional=True),
            objective_check=files.desc("main-objectives/"+case+".toml", optional=True),
            physical=files.desc("main-physical/"+case+"/PHYSICAL.toml", optional=True),
            independent_check=files.desc("main-independent/"+case+".toml", optional=True)))
    index = dict(schema="prospective-window-analysis-index-v1", mode="main", expected_source_commit=G,
        input_protocol=protocol_desc, input_freeze=input_desc, plan=plan_desc, all_orders_freeze=whole_desc,
        cases=cases, partial_attempts=partial, phase_receipts=phases, historical_artifacts=[],
        phase_receipt_scope="Every retained job started and terminal through explicit cutoff; all failures retained",
        ledger_cutoff_UTC=args.cutoff)
    spec = dict(schema="prospective-independent-derived-check-index-v1", expected_scientific_source_commit=G,
        analyzer_source=analyzer, resource_pools=pools, original_campaign_freeze=freeze,
        intended_independent_checker_source=checker,
        guard_source=guard, dispatcher_sources=dispatcher_descs, ledger_cutoff_UTC=args.cutoff,
        jobs=jobs, input_reuse=reuse)
    files.write(args.analysis_index, index)
    spec["analysis_index"] = files.desc(args.analysis_index)
    files.write(args.seed, dict(schema="prospective-main-index-seed-v1", spec_without_summary=spec,
        missing_case_panels=missing_panels, builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        accounting_boundary="Builder/analyzer/checker and all jobs after cutoff require later cumulative reconciliation",
        source_bindings=files.seen))
    print(json.dumps(dict(status="indices_prepared_not_analyzed", ledger_jobs=len(jobs), indexed_cases=len(cases),
                          analysis_index_sha256=spec["analysis_index"]["sha256"], seed_sha256=files.desc(args.seed)["sha256"])))


def finish(args):
    files = Files(args.root); files.fresh(args.output)
    seed = files.read(args.seed)
    need(seed["schema"] == "prospective-main-index-seed-v1", "Index seed schema differs")
    for name, digest in seed["source_bindings"].items():
        need(files.desc(name)["sha256"] == digest, "Prepared metadata/source changed")
    spec = seed["spec_without_summary"]
    need(spec["expected_scientific_source_commit"] == G and "analyzer_summary" not in spec, "Unfinished G seed required")
    jobs, phases = collect_ledger(files, timestamp(spec["ledger_cutoff_UTC"]), spec["resource_pools"],
                                 {d["sha256"] for d in spec["dispatcher_sources"]})
    need(jobs == spec["jobs"], "Cutoff ledger changed after preparation")
    index = files.read(files.bound(spec["analysis_index"])["path"])
    need(index["phase_receipts"] == phases, "Analyzer phase inventory differs from cutoff ledger")
    summary_desc = files.desc(args.summary); summary = files.read(args.summary)
    need(summary["schema"] == "prospective-window-report-data-v1" and summary["mode"] == "main" and
         summary["analysis_index_sha256"] == spec["analysis_index"]["sha256"] and
         summary["analyzer"]["sha256"] == spec["analyzer_source"]["sha256"] == ANALYZER_SHA,
         "Actual summary source/index binding differs")
    # No success requirement: censored, refused and failed results remain evidence.
    spec["analyzer_summary"] = summary_desc
    spec["index_builder_seed"] = files.desc(args.seed)
    files.write(args.output, spec)
    print(json.dumps(dict(status="derived_index_bound_not_checked", ledger_jobs=len(jobs),
                          check_index_sha256=files.desc(args.output)["sha256"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    actions = parser.add_subparsers(dest="action", required=True)
    prepare_parser = actions.add_parser("prepare")
    for flag in ("cutoff", "guard-source", "analyzer-source", "checker-source", "analysis-index", "seed"):
        prepare_parser.add_argument("--"+flag, required=True)
    prepare_parser.add_argument("--dispatcher-source", action="append", required=True)
    finish_parser = actions.add_parser("finish")
    for flag in ("seed", "summary", "output"):
        finish_parser.add_argument("--"+flag, required=True)
    args = parser.parse_args()
    (prepare if args.action == "prepare" else finish)(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(json.dumps(dict(status="index_preparation_refused", failure_type=type(error).__name__,
                              detail="Metadata/source/cutoff/path mismatch; no physical values printed")))
        raise SystemExit(2)
