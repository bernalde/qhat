"""DRAFT/UNRUN: positive G/a5 graph, only after gate PASS; never print quality."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tomllib

R = Path("/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28")
ROOT = R / "runs/prospective_window_ordering_2026-10-28/a1"
WT = Path("/tmp/qhat-normal-window-2026-10-28/worktree")
TOOLS = WT / "numerical_audit/prospective_window_ordering/evidence/tools"
G = "722013aa6d9a0a48776e26de5f356b5a441c3f61"
ANALYZER_SHA = "c6b93ee687ca05c24c898e3204ef83f53f21b4d243a95a86f089c45e1eec8b7e"
DERIVED_SHA = "15cda988729a3fb852aa77aeedf1f152ce2965b0bae929a712961817dc161ee9"
GATE_CHECKER_SHA = "5dd601c640375d3c60982eca0baa65687cb10e945c73f58d27e6bdbd051d705a"
CASE = "LiH_r1p487_631g_14m4e_preflight_7over5"
LIMIT = 32 * 1024**2


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def raw(path):
    with path.open("rb") as stream:
        value = stream.read(LIMIT + 1)
    need(len(value) <= LIMIT, "Metadata byte limit")
    return value


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024**2):
            h.update(block)
    return h.hexdigest()


def desc(name):
    path = ROOT / name
    need(path.resolve().is_relative_to(ROOT.resolve()) and path.is_file(), "Missing/confined original artifact")
    return dict(path=name, sha256=sha(path))


def new_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def analyze(index_name, summary_name):
    output = ROOT / summary_name
    need(not output.exists() and not Path(str(output) + ".partial").exists(), "Fresh analyzer output required")
    run = subprocess.run([sys.executable, str(TOOLS / "analyze_study.py"), str(ROOT), str(output),
                          "--mode", "preflight", "--index", index_name], capture_output=True)
    need(output.is_file(), "Analyzer did not retain a summary")
    return run.returncode, json.loads(raw(output))


def main():
    need(sha(TOOLS / "analyze_study.py") == ANALYZER_SHA and
         sha(TOOLS / "independent_derived_check.py") == DERIVED_SHA,
         "Evidence-tool source differs from pinned G preparation")
    need(sha(R / "coordination/check_preflight_gate_a5.py") == GATE_CHECKER_SHA,
         "Resource gate checker source differs")
    critical = set()
    def bound(name):
        critical.add(ROOT / name)
        return desc(name)
    decision_path = ROOT / bound("PREFLIGHT_DECISION.a5.json")["path"]
    gate = json.loads(raw(decision_path))
    gate_descriptor = bound("PREFLIGHT_GATE_FROZEN.a5.json")
    frozen_gate = json.loads(raw(ROOT / gate_descriptor["path"]))
    need(gate["schema"] == "prospective-window-preflight-decision-v2" and
         gate["decision"] == "passed_resource_and_mechanics_gate" and gate["source_commit"] == G and
         gate["gate_sha256"] == gate_descriptor["sha256"] and gate["checker_sha256"] == GATE_CHECKER_SHA and
         gate["quality_values_inspected_or_used"] is False and gate["windows_completed"] is True and
         gate["objective_replay_complete"] is True and gate["physical_rows_checked"] == 44 and
         gate["independent_reference_times_checked"] == 2 and gate["all_prior_failures_retained"] is True and
         gate["all_resource_pools_unchanged"] is True and
         gate["main64_calls_352_rows"] == "ADMITTED; not yet performed", "G resource/mechanics gate PASS required")
    need(frozen_gate["source_commit"] == G and frozen_gate["attempt_tag"] == "a5" and
         frozen_gate["preflight_gate_checker_sha256"] == GATE_CHECKER_SHA, "Frozen G gate binding differs")
    index = dict(schema="prospective-window-analysis-index-v1", mode="preflight", expected_source_commit=G,
        input_protocol=bound("data-a6/INPUT_PROTOCOL_FROZEN.json"), input_freeze=bound("data-a6/INPUTS_FROZEN.json"),
        preflight_input_freeze=bound("preflight-data-a5/INPUTS_FROZEN.json"),
        plan=bound("preflight-plan.a5.json"), all_orders_freeze=bound("ALL_ORDERS_FROZEN.preflight.a5.json"),
        partial_attempts=[], historical_artifacts=[], phase_receipts=[],
        phase_receipt_scope="Exactly G preflight16 plus heap32 input-identity2; not setup or the complete historical pool ledger")
    entry = dict(case=CASE, inventory=bound("preflight-a5-" + CASE + ".PANEL_INVENTORY.json"),
        panel=bound("preflight-panel-a5/PANEL.json"), candidates=bound("preflight-panel-a5/CANDIDATES.toml"),
        objective_check=bound("preflight-objectives-a5.toml"), physical=bound("preflight-physical-a5/PHYSICAL.toml"),
        independent_check=bound("preflight-independent-a5.toml"))
    index["cases"] = [entry]
    bound("RESOURCE_POOLS.json")
    labels = ["window-preflight-LiH-r1p487-631g-14m4e-preflight-7over5-" + arm + "-r1-a5"
              for arm in ("canonical", "magnitude", "current", "first", "pauli", "cubic",
                          "window-sparse", "window-magnitude")]
    labels += ["window-preflight-" + stage + "-a5" for stage in
               ("alter", "plan", "freeze", "panel", "objectives", "physical", "independent", "gate")]
    labels += ["window-heap32-input-" + stage + "-a6" for stage in ("reuse", "audit")]
    need(set(gate["measurements"]) == set(labels[:15]), "Gate measurement job set differs")
    for number, label in enumerate(labels):
        plan = json.loads(raw(ROOT / bound("jobs/" + label + ".PLAN.json")["path"]))
        final = json.loads(raw(ROOT / bound("jobs/" + label + ".FINAL.json")["path"]))
        receipt = bound("logs/" + label + ".resources.json"); bound("logs/" + label + ".log")
        phase = "preflight" if number < 16 else "identity_inputs"
        need(plan["phase"] == final["phase"] == phase and final["receipt_sha256"] == receipt["sha256"],
             "G phase/terminal receipt binding differs")
        if label in gate["measurements"]:
            measurement = gate["measurements"][label]
            need(measurement["label"] == label and measurement["receipt_sha256"] == receipt["sha256"] and
                 measurement["final_accounting_sha256"] == sha(ROOT / ("jobs/" + label + ".FINAL.json")) and
                 0 <= measurement["owned_cgroup_peak_MiB_upper"] <= 400,
                 "PASS measurement binding/margin differs")
        index["phase_receipts"].append(dict(phase=phase, receipt=receipt))
    physical_path = ROOT / entry["physical"]["path"]; independent_path = ROOT / entry["independent_check"]["path"]
    primary = tomllib.loads(raw(physical_path).decode()); independent = tomllib.loads(raw(independent_path).decode())
    for row in primary["rows"]:
        for field in ("actual_file", "reference_file"):
            path = (physical_path.parent / row[field]).resolve()
            need(path.is_relative_to(ROOT.resolve()), "State path escapes root")
            critical.add(path)
    originals = {path: sha(path) for path in critical}
    try:
        index_name = "ANALYSIS_INDEX.preflight.a5.G_graph_a1.json"
        new_json(ROOT / index_name, index)
        code, positive = analyze(index_name, "ANALYSIS_SUMMARY.preflight.a5.G_graph_a1.json")
        need(code == 0 and positive["status"] == "complete_declared_panel" and len(positive["raw_rows"]) == 44 and
             positive["completion"]["authenticated_independent_physical_rows"] == 44 and
             positive["scientific_source"]["commit"] == G and positive["analyzer"]["sha256"] == ANALYZER_SHA and
             positive["analysis_index_sha256"] == sha(ROOT / index_name), "Positive G graph must check44 rows")
        for name, value in positive["source_artifact_SHA256_map"].items():
            path = ROOT / name
            if path != ROOT / index_name:
                originals.setdefault(path, value["sha256"])
        spec = importlib.util.spec_from_file_location("G_scalar_surface_check", TOOLS / "independent_derived_check.py")
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        key = lambda row: tuple(row[k] for k in ("arm", "repetition", "formula", "time", "budget_multiplier"))
        bykey = {key(row): row for row in primary["rows"]}
        checked_keys = {key(row) for row in independent["rows"]}
        need(len(bykey) == len(primary["rows"]) == len(checked_keys) == len(independent["rows"]) == 44 and
             checked_keys == set(bykey), "Distinct paired original scalar rows required")
        for row in independent["rows"]:
            need(row["status"] == "complete", "Original independent row incomplete")
            module.independent_metric_acceptance(bykey[key(row)], row, row["time"] == 1.6 and row["budget_multiplier"] == 24)
        result = dict(status="G_positive_graph_passed", source_commit=G, physical_rows_checked=44,
            independent_scalar_surface_rows_checked=44, phase_receipts=18, B1_tamper_repeated=False,
            resource_gate="passed_resource_and_mechanics_gate", main_graph_checked=False,
            gate_decision_sha256=sha(decision_path),
            analysis_index_sha256=sha(ROOT / index_name),
            analysis_summary_sha256=sha(ROOT / "ANALYSIS_SUMMARY.preflight.a5.G_graph_a1.json"))
    finally:
        need(all(path.is_file() and sha(path) == before for path, before in originals.items()),
             "An original critical artifact changed")
    print(json.dumps(result))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(json.dumps(dict(status="failed", failure_type=type(error).__name__,
            failure=str(error) if isinstance(error, AssertionError) else "Metadata/IO/schema failure; no quality values printed")))
        raise SystemExit(2)
