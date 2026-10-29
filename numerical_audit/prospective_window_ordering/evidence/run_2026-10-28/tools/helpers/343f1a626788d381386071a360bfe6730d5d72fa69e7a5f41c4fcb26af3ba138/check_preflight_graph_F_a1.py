"""DRAFT/UNRUN: F preflight graph plus exact B1 tamper; never print quality."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tomllib

R = Path(__file__).resolve().parents[1]
ROOT = R / "runs/prospective_window_ordering_2026-10-28/a1"
WT = Path("/tmp/qhat-normal-window-2026-10-28/worktree")
TOOLS = WT / "numerical_audit/prospective_window_ordering/evidence/tools"
F = "a45e836458274a4ff2f37477951c2e92eabb218d"
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


def edit_toml(source, target, section, key, replacement):
    before = raw(source); parsed = tomllib.loads(before.decode())
    expected = copy.deepcopy(parsed)
    (expected if section is None else expected[section])[key] = replacement
    active = section is None; hits = 0; output = []
    spelling = re.escape(key)
    pattern = re.compile(r"^(\s*(?:" + spelling + r'|"' + spelling + r'"|' +
                         "'" + spelling + r"')\s*=\s*)([\"'])([0-9a-f]{64})\2(\s*(?:#.*)?(?:\r?\n)?)$")
    for line in before.decode().splitlines(keepends=True):
        if line.lstrip().startswith("["):
            active = line.strip() == "[" + section + "]" if section is not None else False
        match = pattern.fullmatch(line) if active else None
        if match:
            hits += 1
            line = match[1] + match[2] + replacement + match[2] + match[4]
        output.append(line)
    encoded = "".join(output).encode()
    need(hits == 1 and encoded != before, "Precise single TOML field edit unavailable")
    need(tomllib.loads(encoded.decode()) == expected, "TOML edit changed another semantic field")
    with target.open("xb") as stream:
        stream.write(encoded)


def analyze(index_name, summary_name):
    output = ROOT / summary_name
    need(not output.exists() and not Path(str(output) + ".partial").exists(), "Fresh analyzer output required")
    run = subprocess.run([sys.executable, str(TOOLS / "analyze_study.py"), str(ROOT), str(output),
                          "--mode", "preflight", "--index", index_name], capture_output=True)
    need(output.is_file(), "Analyzer did not retain a summary")
    return run.returncode, json.loads(raw(output))


def messages(summary):
    return [str(item.get("message", "unlabelled failure")) for item in
            summary.get("failures", []) + [f for c in summary.get("cases", []) for f in c.get("failures", [])]]


def main():
    critical = set()
    def bound(name):
        critical.add(ROOT / name)
        return desc(name)
    index = dict(schema="prospective-window-analysis-index-v1", mode="preflight", expected_source_commit=F,
        input_protocol=bound("data-a5/INPUT_PROTOCOL_FROZEN.json"), input_freeze=bound("data-a5/INPUTS_FROZEN.json"),
        preflight_input_freeze=bound("preflight-data-a4/INPUTS_FROZEN.json"),
        plan=bound("preflight-plan.a4.json"), all_orders_freeze=bound("ALL_ORDERS_FROZEN.preflight.a4.json"),
        partial_attempts=[], historical_artifacts=[], phase_receipts=[],
        phase_receipt_scope="Exactly F preflight16 plus input-identity2; not the complete historical pool ledger")
    entry = dict(case=CASE, inventory=bound("preflight-a4-" + CASE + ".PANEL_INVENTORY.json"),
        panel=bound("preflight-panel-a4/PANEL.json"), candidates=bound("preflight-panel-a4/CANDIDATES.toml"),
        objective_check=bound("preflight-objectives-a4.toml"), physical=bound("preflight-physical-a4/PHYSICAL.toml"),
        independent_check=bound("preflight-independent-a4.toml"))
    index["cases"] = [entry]
    bound("RESOURCE_POOLS.json"); bound("PREFLIGHT_GATE_FROZEN.a4.json")
    labels = ["window-preflight-LiH-r1p487-631g-14m4e-preflight-7over5-" + arm + "-r1-a4"
              for arm in ("canonical", "magnitude", "current", "first", "pauli", "cubic",
                          "window-sparse", "window-magnitude")]
    labels += ["window-preflight-" + stage + "-a4" for stage in
               ("alter", "plan", "freeze", "panel", "objectives", "physical", "independent", "gate")]
    labels += ["window-source-review-input-" + stage + "-a5" for stage in ("reuse", "audit")]
    for number, label in enumerate(labels):
        plan = json.loads(raw(ROOT / bound("jobs/" + label + ".PLAN.json")["path"]))
        final = json.loads(raw(ROOT / bound("jobs/" + label + ".FINAL.json")["path"]))
        receipt = bound("logs/" + label + ".resources.json"); bound("logs/" + label + ".log")
        phase = "preflight" if number < 16 else "identity_inputs"
        need(plan["phase"] == final["phase"] == phase and final["receipt_sha256"] == receipt["sha256"],
             "F phase/terminal receipt binding differs")
        index["phase_receipts"].append(dict(phase=phase, receipt=receipt))
    decision_path = ROOT / "PREFLIGHT_DECISION.a4.json"; critical.add(decision_path)
    gate = json.loads(raw(decision_path))
    need(gate["decision"] == "failed_memory_margin" and gate["main64_352_admitted"] is False and
         gate["criterion_MiB"] == 400 and gate["actual_max_cgroup_display_upper_MiB"] == 402.8 and
         gate["source_commit"] == F and gate["physical_quality_inspected"] is False, "Retained failed gate differs")
    physical_path = ROOT / entry["physical"]["path"]; independent_path = ROOT / entry["independent_check"]["path"]
    primary = tomllib.loads(raw(physical_path).decode()); independent = tomllib.loads(raw(independent_path).decode())
    for row in primary["rows"]:
        for field in ("actual_file", "reference_file"):
            path = (physical_path.parent / row[field]).resolve()
            need(path.is_relative_to(ROOT.resolve()), "State path escapes root")
            critical.add(path)
    originals = {path: sha(path) for path in critical}
    try:
        index_name = "ANALYSIS_INDEX.preflight.a4.json"
        new_json(ROOT / index_name, index)
        code, positive = analyze(index_name, "ANALYSIS_SUMMARY.preflight.a4.F_graph_a1.json")
        if code != 0 or positive["status"] != "complete_declared_panel":
            print(json.dumps(dict(status="positive_graph_failed", failures=messages(positive))))
        need(code == 0 and positive["status"] == "complete_declared_panel" and len(positive["raw_rows"]) == 44 and
             positive["completion"]["authenticated_independent_physical_rows"] == 44, "Positive F graph must check44 rows")
        for name, value in positive["source_artifact_SHA256_map"].items():
            path = ROOT / name
            if path != ROOT / index_name:
                originals.setdefault(path, value["sha256"])
        spec = importlib.util.spec_from_file_location("F_scalar_surface_check", TOOLS / "independent_derived_check.py")
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        key = lambda row: tuple(row[k] for k in ("arm", "repetition", "formula", "time", "budget_multiplier"))
        bykey = {key(row): row for row in primary["rows"]}
        need(len(bykey) == len(independent["rows"]) == 44, "Distinct original scalar rows required")
        for row in independent["rows"]:
            need(row["status"] == "complete", "Original independent row incomplete")
            module.independent_metric_acceptance(bykey[key(row)], row, row["time"] == 1.6 and row["budget_multiplier"] == 24)
        tamper_primary = physical_path.with_name("PHYSICAL.B1_TAMPER_F_a1.toml")
        edit_toml(physical_path, tamper_primary, "source_hashes", "physical.jl", "0" * 64)
        tamper_independent = ROOT / "preflight-independent-a4.B1_TAMPER_F_a1.toml"
        edit_toml(independent_path, tamper_independent, None, "manifest_sha256", sha(tamper_primary))
        variant = copy.deepcopy(index)
        variant["cases"][0]["physical"] = desc(str(tamper_primary.relative_to(ROOT)))
        variant["cases"][0]["independent_check"] = desc(tamper_independent.name)
        variant_name = "ANALYSIS_INDEX.preflight.a4.B1_TAMPER_F_a1.json"; new_json(ROOT / variant_name, variant)
        rejected_code, rejected = analyze(variant_name, "ANALYSIS_SUMMARY.preflight.a4.B1_TAMPER_F_a1.json")
        need(rejected_code == 2 and rejected["status"] != "complete_declared_panel" and
             any("Primary physical scientific source differs" in m for m in messages(rejected)),
             "B1 variant must fail specifically on primary source binding")
        print(json.dumps(dict(status="F_graph_and_B1_tamper_passed", physical_rows_checked=44,
            independent_scalar_surface_rows_checked=44, phase_receipts=18, B1_source_rejection=True,
            resource_gate="failed_memory_margin", main_admitted=False)))
    finally:
        need(all(path.is_file() and sha(path) == before for path, before in originals.items()),
             "An original critical artifact changed")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(json.dumps(dict(status="failed", failure_type=type(error).__name__,
            failure=str(error) if isinstance(error, AssertionError) else "Metadata/IO/schema failure; no quality values printed")))
        raise SystemExit(2)
