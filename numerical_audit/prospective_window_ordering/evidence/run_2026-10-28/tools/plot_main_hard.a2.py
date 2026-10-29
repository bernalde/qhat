"""UNRUN display helper; fixed T=1.6, cap=24U. No scientific recomputation."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path, PurePosixPath

FLOOR = 1e-12
CASES = ("LiH_r1p487_631g_14m4e", "LiH_r2p113_631g_14m4e",
         "BeH2_r1p369_631g_14m4e", "BeH2_r1p907_631g_14m4e")
PATHS = ("canonical_f1", "canonical_f2", "magnitude_f2", "current_f1", "current_f2",
         "first_f1", "first_f2", "pauli_f2", "cubic_f2", "window_sparse_f2", "window_magnitude_f2")
COLORS = ("#1f77b4", "#aec7e8", "#ff7f0e", "#2ca02c", "#98df8a", "#d62728",
          "#ff9896", "#9467bd", "#8c564b", "#e377c2", "#17becf")
METRICS = ("projective_error", "infidelity", "norm_drift", "number_leakage", "spin_leakage")
MARKERS = {1: "o", 2: "s"}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def confined(root, name):
    need(isinstance(name, str) and 0 < len(name) <= 2048, "Relative path required")
    rel = PurePosixPath(name)
    need(not rel.is_absolute() and ".." not in rel.parts and "\\" not in name, "Path escape")
    path = root.joinpath(*rel.parts).resolve()
    need(path.is_relative_to(root), "Symlink escape")
    return path


def load(root, name, expected=None, binary=False):
    path = confined(root, name)
    need(path.is_file() and path.stat().st_size <= 32*1024**2, "Missing/oversize input")
    with path.open("rb") as stream:
        raw = stream.read(32*1024**2+1)
    need(len(raw) <= 32*1024**2, "Input byte ceiling")
    digest = hashlib.sha256(raw).hexdigest()
    need(expected is None or digest == expected, "Input hash differs")
    return (None if binary else json.loads(raw)), dict(path=name, sha256=digest)


def numeric(value):
    need(type(value) in (int, float) and math.isfinite(value) and value >= 0, "Invalid measured scalar")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path); parser.add_argument("check_index")
    parser.add_argument("check_report"); parser.add_argument("output_directory")
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--checker-source", default="tools/independent_derived_check.py")
    args = parser.parse_args(); root = args.root.resolve()
    spec, index_desc = load(root, args.check_index)
    report, report_desc = load(root, args.check_report, args.report_sha256)
    need(spec["schema"] == "prospective-independent-derived-check-index-v1" and
         report["schema"] == "prospective-independent-derived-check-v1" and
         report["check_index_sha256"] == index_desc["sha256"], "Derived index/report binding")
    allowed = {"verified_complete_derived_panel", "verified_censored_derivations", "derived_agreement_with_accounting_veto"}
    need(report["status"] in allowed and report.get("derived_values_agree") is True and
         not report["failures"] and report["raw_rows_checked"] == 352, "Actual checked derived agreement required")
    _, checker = load(root, args.checker_source, report["checker_sha256"], binary=True)
    if "intended_independent_checker_source" in spec:
        need(checker == spec["intended_independent_checker_source"], "Intended checker source differs")
    summary, summary_desc = load(root, spec["analyzer_summary"]["path"], spec["analyzer_summary"]["sha256"])
    analysis, analysis_desc = load(root, spec["analysis_index"]["path"], spec["analysis_index"]["sha256"])
    _, analyzer = load(root, spec["analyzer_source"]["path"], spec["analyzer_source"]["sha256"], binary=True)
    need(analysis["schema"] == "prospective-window-analysis-index-v1" and analysis["mode"] == summary["mode"] == "main" and
         summary["schema"] == "prospective-window-report-data-v1" and summary["analysis_index_sha256"] == analysis_desc["sha256"] and
         summary["analyzer"]["sha256"] == analyzer["sha256"], "Analyzer graph/source binding")
    need(report["final_source_head"] == spec["expected_scientific_source_commit"] == analysis["expected_source_commit"] ==
         summary["scientific_source"]["commit"] and summary["ambiguity"] == summary["ratio_floor"] == report["ambiguity"] == FLOOR,
         "Scientific source/precision binding")
    for descriptor in (summary_desc, analysis_desc, analyzer):
        need(report["artifact_bindings"][descriptor["path"]]["sha256"] == descriptor["sha256"], "Checked artifact differs")
    expected = set(itertools.product(CASES, PATHS, (1, 2), (.4, 1.6), (8, 24)))
    keyed = {(r["case"], r["path"], r["repetition"], r["time"], r["budget_multiplier"]): r for r in summary["raw_rows"]}
    need(len(summary["raw_rows"]) == len(keyed) == 352 and set(keyed) == expected, "Complete logical roster required")
    rows, source_rows = [], []
    for case, path, rep in itertools.product(CASES, PATHS, (1, 2)):
        r = keyed[(case, path, rep, 1.6, 24)]; source_rows.append(r)
        need(path == r["arm"]+"_f"+str(r["formula"]), "Path/formula identity differs")
        need(type(r["physical_checked"]) is bool and type(r["objective_accepted"]) is bool, "Explicit check flags required")
        values = r.get("metrics") or {}
        need(r["physical_checked"] or not values, "Unchecked row cannot invent physical metrics")
        point = dict(case=case, path=path, repetition=rep, time=1.6, budget_multiplier=24,
            status=r["status"], physical_checked=r["physical_checked"], objective_accepted=r["objective_accepted"],
            objective_status=(r.get("objective") or {}).get("status"), censor_reason=r.get("censor_reason"))
        point.update({m: numeric(values[m]) if r["physical_checked"] else None for m in METRICS})
        for axis in ("CPU", "wall"):
            cost = r["constructor_"+axis+"_seconds"]
            point[axis+"_seconds"] = None if cost is None else numeric(cost)
            point[axis+"_plotted"] = r["physical_checked"] and cost is not None
        error = point["projective_error"]
        point.update(display_floor=error is not None and error <= FLOOR,
                     displayed_error=None if error is None else max(error, FLOOR))
        rows.append(point)
    need(any(p["physical_checked"] for p in rows), "No checked hard-condition physical data to display")
    output = confined(root, args.output_directory); need(not output.exists(), "Fresh output directory required")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    output.mkdir(parents=True, exist_ok=False)
    bindings = dict(check_index=index_desc, check_report=report_desc, analyzer_summary=summary_desc,
                    analysis_index=analysis_desc, analyzer_source=analyzer, checker_source=checker,
                    plotter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    payload = dict(schema="normal-hard-condition-plot-data-v1", expected_points=88, rows=rows, source_rows=source_rows,
        bindings=bindings, scientific_source=report["final_source_head"], derived_status=report["status"],
        ledger_cutoff_UTC=report["ledger_cutoff_UTC"], complete_352_derived_claim_allowed=report["complete_352_derived_claim_allowed"],
        unqualified_under_ceiling_claim_allowed=report["unqualified_under_ceiling_claim_allowed"],
        display_policy=dict(time=1.6, cap_times_original_words=24, error_floor=FLOOR, palette=dict(zip(PATHS, COLORS)),
            repetitions=MARKERS, x="raw guard CPU / absolute whole-call wall; not conservative service charge",
            y="projective_error; floor only for display", missing="null; CSV empty field; never zero-imputed",
            overlap="No jitter; coincident observations remain separate table rows", intervals="none"))
    with (output/"plotted_data.json").open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False); stream.write("\n")
    csv_rows = [dict(p, summary_sha256=summary_desc["sha256"], check_report_sha256=report_desc["sha256"],
                     check_index_sha256=index_desc["sha256"], scientific_source=report["final_source_head"]) for p in rows]
    with (output/"plotted_data.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0])); writer.writeheader(); writer.writerows(csv_rows)
    legend = [Line2D([], [], color=c, label=p.replace("_", " ")) for p, c in zip(PATHS, COLORS)]
    legend += [Line2D([], [], color="black", marker=m, linestyle="none", label="repetition "+str(r)) for r, m in MARKERS.items()]
    legend += [Line2D([], [], color="black", marker="o", markerfacecolor="none", linestyle="none", label="objective unaccepted"),
               Line2D([], [], color="black", marker="v", linestyle="none", label="at/below display floor")]
    ymax = max([FLOOR*10]+[p["projective_error"]*1.5 for p in rows if p["physical_checked"]])
    for axis, label in (("CPU", "Raw guard whole-call CPU (s)"), ("wall", "Absolute whole-call wall (s)")):
        fig, panels = plt.subplots(2, 2, figsize=(13, 9), sharey=True)
        for panel, case in zip(panels.flat, CASES):
            group = [p for p in rows if p["case"] == case]
            for p in group:
                if not p[axis+"_plotted"]: continue
                color = COLORS[PATHS.index(p["path"])]; x, y = p[axis+"_seconds"], p["displayed_error"]
                panel.scatter(x, y, s=48, marker=MARKERS[p["repetition"]], edgecolors=color,
                              facecolors=color if p["objective_accepted"] else "none", linewidths=1.4)
                if p["display_floor"]: panel.scatter(x, y, s=16, marker="v", color="black", zorder=4)
            panel.set(title=case.replace("_631g_14m4e", ""), xlabel=label, ylabel="Projective physical error", xscale="linear", yscale="log")
            panel.set_ylim(FLOOR/2, ymax); panel.grid(alpha=.2)
            counts = (sum(p[axis+"_plotted"] for p in group), sum(not p["physical_checked"] for p in group),
                      sum(p[axis+"_seconds"] is None for p in group), sum(p["physical_checked"] and not p["objective_accepted"] for p in group))
            panel.text(.98, .98, "22 expected; %d plotted\n%d physical missing/censored; %d cost unknown\n%d objective unaccepted" % counts,
                       transform=panel.transAxes, ha="right", va="top", fontsize=8)
        fig.suptitle("Fixed hard condition: T=1.6, budget=24U — both fresh repetitions")
        fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(.5, .055), ncol=5, fontsize=8)
        fig.text(.5, .018, "Errors <=1e-12 displayed at 1e-12; table values unchanged. Coincident points may overlap; no jitter or intervals.\n"
                 "Raw costs differ from conservative charges. Derived status: "+report["status"], ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .19, 1, .95))
        for suffix in ("png", "svg"): fig.savefig(output/(axis.lower()+"_vs_error."+suffix), dpi=180, bbox_inches="tight")
        plt.close(fig)
    manifest = dict(bindings=bindings, outputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir())})
    with (output/"MANIFEST.json").open("x") as stream: json.dump(manifest, stream, indent=2); stream.write("\n")
    print(json.dumps(dict(status="display_artifacts_written", expected_points=88, figures=2, derived_status=report["status"])))


if __name__ == "__main__":
    main()
