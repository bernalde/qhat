#!/usr/bin/env python3
"""Standalone consumer for molecular-analysis-v1; no scientific gate is changed.

Run only after real analysis exists: plot.py ANALYSIS.json FRESH_DIRECTORY
Use --preflight explicitly for the old-case diagnostic. No training or physics
module is imported. Raw values, missing records and provenance accompany figures.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

GRID = (("parent", 0), ("canonical", 0), ("magnitude", 0), ("full", 32),
        ("full", 64), ("window", 32), ("window", 64), ("annealing", 32), ("annealing", 64))
COLORS = dict(parent="#666666", canonical="#111111", magnitude="#009E73",
              full="#0072B2", window="#D55E00", annealing="#CC79A7")
MARKERS = dict(parent="P", canonical="s", magnitude="D", full="o", window="^", annealing="v")
FLOOR = 1e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, "Duplicate JSON key: "+key)
        result[key] = value
    return result


def reject_nonfinite(value):
    raise ValueError("Nonfinite JSON constant: "+value)


def finite_nonnegative(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def identity(row):
    return tuple(row[key] for key in ("case", "repetition", "arm", "cpu_budget"))


def collect_points(result):
    require(result["schema"] == "molecular-analysis-v1", "Expected the actual analysis schema")
    cases, reps = result["cases"], result["repetitions"]
    require(isinstance(cases, list) and cases and len(cases) == len(set(cases)) and
            all(isinstance(case, str) and re.fullmatch(r"[A-Za-z0-9_.-]{1,160}", case) for case in cases),
            "Declared unique case names required")
    require(type(reps) is int and reps == (1 if result["mode"] == "preflight" else 2), "Frozen repetitions")
    require(result["mode"] in ("main", "preflight") and len(cases) == (1 if reps == 1 else 2), "Frozen case/mode grid")
    physical = [row for row in result.get("retained_physical_rows", [])
                if row.get("time") == 1.35 and row.get("budget_multiplier") == 8]
    costs = result.get("constructor_costs", {}).get("constructors", [])
    by_physical, by_cost = {}, {}
    for row in physical:
        by_physical.setdefault(identity(row), []).append(row)
    for row in costs:
        by_cost.setdefault(identity(row), []).append(row)
    case_status = {row["case"]: row.get("complete") for row in result.get("cases_metadata", [])}
    points = []
    for case in cases:
        for rep in range(1, reps+1):
            for arm, budget in GRID:
                key = (case, rep, arm, budget)
                observed, receipts = by_physical.get(key, []), by_cost.get(key, [])
                row = observed[0] if len(observed) == 1 else {}
                cost = receipts[0] if len(receipts) == 1 else {}
                own = cost.get("actual_own", {})
                attributed = cost.get("comparison_attributed", {})
                # Complete optimized comparison includes the measured parent.
                complete_cost = attributed if arm in ("full", "window", "annealing") else own
                x = complete_cost.get("charged_CPU_seconds")
                y = row.get("projective_error")
                flags = []
                if len(observed) != 1:
                    flags.append("missing error" if not observed else "duplicate errors")
                if len(receipts) != 1:
                    flags.append("missing cost" if not receipts else "duplicate costs")
                if not finite_nonnegative(x):
                    flags.append("CPU unavailable")
                if not finite_nonnegative(y):
                    flags.append("error unavailable")
                if row.get("eligible") is False or cost.get("eligible") is False:
                    flags.append("ineligible")
                elif row.get("eligible") is not True or cost.get("eligible") is not True:
                    flags.append("eligibility unavailable")
                floor = row.get("numerical_floor") is True or (finite_nonnegative(y) and y <= FLOOR)
                if floor:
                    flags.append("floor")
                if finite_nonnegative(y) and row.get("numerical_floor") is not (y <= FLOOR):
                    flags.append("floor flag missing/different")
                points.append(dict(case=case, repetition=rep, arm=arm, cpu_budget=budget,
                    complete_construction_charged_CPU_seconds=x, raw_projective_error=y,
                    actual_own=own, comparison_attributed=attributed,
                    constructor_eligible=cost.get("eligible"), physical_eligible=row.get("eligible"),
                    case_evidence_complete=case_status.get(case), analysis_complete=result.get("complete"),
                    floor=floor, flags=flags, raw_physical_matches=observed, raw_cost_matches=receipts))
    expected = {identity(point) for point in points}
    require(set(by_physical) <= expected and set(by_cost) <= expected, "Unexpected primary/cost identities; refuse silent omission")
    return points, physical, costs


def write_json(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def write_csv(path, points):
    fields = ("case", "repetition", "arm", "cpu_budget", "complete_construction_charged_CPU_seconds",
              "raw_projective_error", "constructor_eligible", "physical_eligible", "floor", "flags",
              "actual_own", "comparison_attributed", "raw_physical_matches", "raw_cost_matches")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for point in points:
            writer.writerow({key: json.dumps(point[key], allow_nan=False) if isinstance(point[key], (dict, list))
                             else point[key] for key in fields})


def plot_case(plt, Line2D, case, points, source, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 10,
                         "axes.titlesize": 11, "svg.fonttype": "none", "savefig.dpi": 300}):
        fig, (ax, table) = plt.subplots(1, 2, figsize=(11, 9), gridspec_kw={"width_ratios": [2.2, 1]})
        fig.subplots_adjust(left=.085, right=.985, top=.84, bottom=.30, wspace=.15)
        fig.suptitle(case+" — T = 1.35, rotation cap = 8U", fontsize=12, y=.965)
        fig.text(.085, .90, "Source analysis status: "+str(source.get("status", "unavailable"))+
                 "; complete = "+str(source.get("complete", False)), fontsize=9)
        fig.text(.085, .865, "One original unpromoted HF window native-validation timeout keeps both primary decisions incomplete; all plotted selected rows were independently replayed.", fontsize=8)
        for index, point in enumerate(points, 1):
            x, y, arm = point["complete_construction_charged_CPU_seconds"], point["raw_projective_error"], point["arm"]
            label = arm+(f" B{point['cpu_budget']}" if arm in ("full", "window", "annealing") else "")+f" / r{point['repetition']}"
            status = ", ".join(point["flags"]) or "recorded"
            table.text(0, .98-(index-1)*(.98/len(points)), f"{index:2}. {label}\n     {status}", fontsize=8,
                       va="top", linespacing=1.05)
            if not (finite_nonnegative(x) and finite_nonnegative(y)):
                continue  # No fabricated coordinate for a missing value.
            size = {0: 42, 32: 55, 64: 90}[point["cpu_budget"]]
            ax.scatter(x, y, s=size, marker=MARKERS[arm], edgecolors=COLORS[arm],
                       facecolors=COLORS[arm] if point["repetition"] == 1 else "none", linewidths=1.1, zorder=3, clip_on=False)
            if "ineligible" in point["flags"] or "eligibility unavailable" in point["flags"]:
                ax.scatter(x, y, s=size+20, marker="x", color="black", linewidths=1.0, zorder=5)
            if point["floor"]:
                ax.scatter(x, y, s=25, marker="v", facecolors="none", edgecolors="black", linewidths=.8, zorder=6)
            # Offset text only: coincident raw coordinates remain unchanged.
            level = (index-1) % 3
            label_y = 6+12*level if point["repetition"] == 1 else -12-12*level
            ax.annotate(str(index), (x, y), xytext=(4, label_y),
                        textcoords="offset points", fontsize=7)
        ax.set_yscale("linear")
        top = max((p["raw_projective_error"] for p in points if finite_nonnegative(p["raw_projective_error"])), default=0.)
        right = max((p["complete_construction_charged_CPU_seconds"] for p in points if finite_nonnegative(p["complete_construction_charged_CPU_seconds"])), default=0.)
        ax.set_ylim(0., 1.5*max(top, FLOOR)); ax.set_xlim(0., max(1., 1.1*right))
        ax.axhline(FLOOR, color="#777777", linestyle=":", linewidth=.8)
        ax.text(.015, FLOOR, " 1e−12 floor flag threshold", transform=ax.get_yaxis_transform(),
                fontsize=7, va="bottom", color="#555555")
        ax.set_xlabel("Measured complete construction charged CPU (s)")
        ax.set_ylabel("Raw projective error (linear; values unchanged)")
        ax.grid(True, which="major", color="#dddddd", linewidth=.5)
        ax.spines[["top", "right"]].set_visible(False); table.axis("off")
        handles = [Line2D([], [], color=COLORS[arm], marker=MARKERS[arm], linestyle="none", label=arm)
                   for arm in COLORS]
        fig.legend(handles=handles, ncol=6, loc="center", bbox_to_anchor=(.52, .19), frameon=False, fontsize=8)
        fig.text(.085, .13, "Filled/open: repetition 1/2. Small/large optimized markers: B32/B64. ×: ineligible or eligibility unavailable. ▽: floor.", fontsize=8)
        fig.text(.085, .09, "Optimized x includes the measured common parent prefix plus continuation. Parent and cheap controls use complete own costs.", fontsize=8)
        fig.text(.085, .05, "Shared parent spending is counted once in study totals. Exact raw coordinates are retained; missing values have no invented location. No confidence bands.", fontsize=8)
        for extension in ("svg", "png"):
            target = output/(case+".primary_quality_cost."+extension)
            require(not target.exists(), "Fresh figure output required")
            fig.savefig(target, dpi=300, metadata={"Creator": "Molecular generalization plotting consumer"} if extension == "svg" else None)
        plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args(argv)
    require(args.result.is_file(), "Existing analysis result required")
    observed = args.result.stat().st_size
    require(observed <= 64*1024**2, "Bounded analysis result required")
    with args.result.open("rb") as stream:
        raw = stream.read(observed+1)
        require(len(raw) == observed and not stream.read(1), "Analysis result size changed during read")
    result = json.loads(raw, object_pairs_hook=unique_object, parse_constant=reject_nonfinite)
    require((result["mode"] == "preflight") == args.preflight, "Explicit preflight mode must match the result")
    points, physical, costs = collect_points(result)
    args.output.mkdir(parents=True, exist_ok=False)
    script = Path(__file__).read_bytes()
    with (args.output/"plot_consumer.py").open("xb") as stream:
        stream.write(script)
    write_json(args.output/"PLOT_DATA.json", dict(schema="molecular-primary-plot-data-v1", result_sha256=sha(raw),
        points=points, raw_primary_rows=physical, raw_cost_records=costs,
        actual_shared_study_costs=result.get("constructor_costs"), primary_outcomes_unchanged=result.get("primary")))
    write_csv(args.output/"PLOT_DATA.csv", points)
    os.environ["MPLCONFIGDIR"] = str(args.output/".matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    for case in result["cases"]:
        plot_case(plt, Line2D, case, [point for point in points if point["case"] == case], result, args.output)
    outputs = {path.name: sha(path.read_bytes()) for path in args.output.iterdir() if path.is_file()}
    write_json(args.output/"PLOT_PROVENANCE.json", dict(schema="molecular-plot-consumer-v1",
        created_UTC=datetime.now(timezone.utc).isoformat(), input_result_path=str(args.result.resolve()),
        input_result_sha256=sha(raw), input_analysis_complete=result.get("complete"),
        protocol_sha256=result.get("protocol_sha256"), mode=result["mode"], repetitions=result["repetitions"],
        original_producer_sources=result.get("original_frozen_sources"),
        original_source_freeze_sha256=result.get("source_freeze_sha256"),
        analysis_consumer_sources=result.get("current_consumer_sources"),
        plotting_consumer=dict(source_sha256=sha(script), python=sys.version,
                               matplotlib=matplotlib.__version__, backend=matplotlib.get_backend()),
        outputs=outputs, scientific_status="Presentation only; preserves input completeness and never supplies a missing scientific gate"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
