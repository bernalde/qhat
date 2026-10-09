#!/usr/bin/env python3
"""Summarize completed frozen-input audits without changing their artifacts.

Usage:
  python numerical_audit/summarize_frozen_audit.py RUN [RUN ...] --output NEW_DIR

Counts describe the supplied fixed cases, not a molecular population. Repeated
runs are checked separately and never counted as independent physics samples.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qhat-frozen-summary-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ORDERS = ("fermionic_signed", "jw_signed", "jw_magnitude")
LABELS = {"fermionic_signed": "Fermionic signed", "jw_signed": "JW signed",
          "jw_magnitude": "JW magnitude"}
COLORS = {"fermionic_signed": "#0072B2", "jw_signed": "#D55E00", "jw_magnitude": "#009E73"}
MARKERS = {"fermionic_signed": "o", "jw_signed": "s", "jw_magnitude": "^"}
CASE_LABELS = {"f2": r"F$_2$ (16 qubits)", "li2": r"Li$_2$ (12 qubits)",
               "nh3": r"NH$_3$ (16 qubits)"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def dump(path, obj):
    with Path(path).open("x") as stream:
        json.dump(obj, stream, indent=2, allow_nan=False)
        stream.write("\n")


def key(row):
    return row["case"], float(row["T"]), int(row["steps"]), row["ordering"]


def classify(left, right, args, *, floor=True):
    if floor and min(left, right) <= args.error_floor:
        return "below_floor"
    if math.isclose(left, right, rel_tol=args.tie_rtol, abs_tol=args.tie_atol):
        return "tie"
    return "left_better" if left < right else "right_better"


def distribution(values):
    values = [float(x) for x in values if x is not None and math.isfinite(x)]
    if not values:
        return {"count": 0, "min": None, "median": None, "max": None}
    return {"count": len(values), "min": min(values), "median": float(np.median(values)), "max": max(values)}


def load_runs(paths):
    """Require completed data and equal physical hashes before pooling repeats."""
    rows, checks, provenance, repeats = {}, {}, [], []
    for supplied in paths:
        run = supplied.resolve()
        manifest = read(run / "manifest.json")
        completion = read(run / "completion.json")
        if completion.get("status") != "complete":
            raise ValueError(f"Run is not complete: {run}")
        data = read(run / "results.json")
        expected = {(c, float(t), int(r), o) for c in manifest["cases"]
                    for t in manifest["times"] for r in manifest["steps"] for o in ORDERS}
        if len(data) != len(expected) or {key(row) for row in data} != expected:
            raise ValueError(f"Missing or duplicate grid rows: {run}")
        if completion["rows"] != len(data):
            raise ValueError(f"Completion count disagrees: {run}")
        for case in manifest["cases"]:
            check = read(run / f"{case}.checks.json")
            if not (check["hamiltonian_and_all_schedule_hashes_match"]
                    and check["regenerated_baseline_orders_match"]):
                raise ValueError(f"Failed physical hash validation: {run}/{case}")
            if case in checks:
                old = checks[case]
                if old["input_sha256"] != check["input_sha256"] or old["fingerprints"] != check["fingerprints"]:
                    raise ValueError(f"Different physical input/order under reused case label: {case}")
            else:
                checks[case] = check
        for row in data:
            if row["ordering"] not in ORDERS:
                raise ValueError(f"Unexpected ordering: {row['ordering']}")
            for field in ("one_minus_overlap", "infidelity", "exact_norm", "approx_norm"):
                if not math.isfinite(row[field]):
                    raise ValueError(f"Non-finite metric {field}: {key(row)}")
            if not (0 <= row["one_minus_overlap"] <= 1 and 0 <= row["infidelity"] <= 1):
                raise ValueError(f"Invalid physical error: {key(row)}")
            if key(row) in rows:
                old = rows[key(row)]
                delta = row["one_minus_overlap"] - old["one_minus_overlap"]
                repeats.append({"case": row["case"], "T": row["T"], "steps": row["steps"],
                                "ordering": row["ordering"], "run": str(run),
                                "absolute_difference": delta,
                                "relative_difference": delta / old["one_minus_overlap"] if old["one_minus_overlap"] else None,
                                "within_repeat_tolerance": math.isclose(row["one_minus_overlap"], old["one_minus_overlap"], rel_tol=1e-8, abs_tol=1e-20)})
            else:
                rows[key(row)] = row
        provenance.append({"path": str(run), "manifest_sha256": sha(run / "manifest.json"),
                           "results_sha256": sha(run / "results.json"),
                           "completion_sha256": sha(run / "completion.json"),
                           "checks_sha256": {c: sha(run / f"{c}.checks.json") for c in manifest["cases"]},
                           "manifest": manifest, "completion": completion})
    return sorted(rows.values(), key=key), checks, provenance, repeats


def summarize(rows, checks, args):
    grouped = defaultdict(dict)
    for row in rows:
        grouped[(row["case"], row["T"], row["steps"])][row["ordering"]] = row
    comparisons, prediction_pairs, anchor, time_outcomes = [], [], [], defaultdict(Counter)
    for (case, T, steps), group in sorted(grouped.items()):
        if set(group) != set(ORDERS):
            raise ValueError(f"Incomplete ordering comparison: {case, T, steps}")
        counts = {r["nominal_exponential_count"] for r in group.values()}
        if counts != {steps * checks[case]["pauli_terms"]}:
            raise ValueError(f"Unmatched Pauli work counts: {case, T, steps}")
        for other in ORDERS[1:]:
            left, right = group["fermionic_signed"], group[other]
            a, b = left["one_minus_overlap"], right["one_minus_overlap"]
            outcome = classify(a, b, args)
            naive = classify(left["naive_unnormalized_one_minus_overlap"], right["naive_unnormalized_one_minus_overlap"], args)
            comparison = {"case": case, "T": T, "steps": steps, "baseline": other,
                          "outcome_for_fermionic": {"left_better": "win", "right_better": "loss"}.get(outcome, outcome),
                          "baseline_E_over_fermionic_E": b / a if a else None,
                          "fermionic_E": a, "baseline_E": b,
                          "legacy_naive_fermionic_E":left["naive_unnormalized_one_minus_overlap"],
                          "legacy_naive_baseline_E":right["naive_unnormalized_one_minus_overlap"],
                          "old_naive_classification_disagrees": naive != outcome,
                          "old_naive_outcome": naive}
            comparisons.append(comparison)
            if steps == args.anchor_steps:
                time_outcomes[(T, other)][comparison["outcome_for_fermionic"]] += 1
        if steps == args.anchor_steps:
            errors = {o: group[o]["one_minus_overlap"] for o in ORDERS}
            best = min(errors, key=errors.get)
            winners = [o for o in ORDERS if classify(errors[o], errors[best], args, floor=False) == "tie"]
            anchor.append({"case": case, "T": T, "steps": steps, "errors": errors,
                           "best_orderings_with_tolerance": winners,
                           "winner_below_floor": errors[best] <= args.error_floor})
        for first, second in itertools.combinations(ORDERS, 2):
            a, b = group[first], group[second]
            if min(a["one_minus_overlap"], b["one_minus_overlap"]) <= args.error_floor:
                continue
            if abs(a["infidelity"] - b["infidelity"]) <= args.rank_gap * max(a["infidelity"], b["infidelity"]):
                continue
            actual_sign = np.sign(a["infidelity"] - b["infidelity"])
            norms = checks[case]["initial_BCH"]
            static_a, static_b = norms[first]["norm"] ** 2, norms[second]["norm"] ** 2
            pair = {"case": case, "T": T, "steps": steps, "first": first, "second": second,
                    "actual_winner": first if actual_sign < 0 else second,
                    "static_HF_BCH_correct": bool(np.sign(static_a - static_b) == actual_sign),
                    "dynamic_correct": None}
            if "predicted_infidelity" in a and "predicted_infidelity" in b:
                pair["dynamic_correct"] = bool(np.sign(a["predicted_infidelity"] - b["predicted_infidelity"]) == actual_sign)
            prediction_pairs.append(pair)
    scaling = []
    series = defaultdict(list)
    for row in rows:
        series[(row["case"], row["T"], row["ordering"])].append(row)
    for (case, T, ordering), group in sorted(series.items()):
        eligible = sorted((r for r in group if r["one_minus_overlap"] > args.error_floor), key=lambda r:r["steps"])
        pairs = []
        for first, last in zip(eligible[:-1], eligible[1:]):
            ratio = first["one_minus_overlap"] / last["one_minus_overlap"]
            expected = (last["steps"] / first["steps"]) ** 2
            pairs.append({"steps": [first["steps"], last["steps"]], "E_ratio": ratio,
                          "inverse_square_ratio": expected, "relative_ratio_deviation": ratio / expected - 1})
        slope = None
        if len(eligible) >= 2:
            first, last = eligible[-2:]
            slope = math.log(last["one_minus_overlap"] / first["one_minus_overlap"]) / math.log(last["steps"] / first["steps"])
        scaling.append({"case": case, "T": T, "ordering": ordering, "eligible_points":len(eligible),
                        "last_two_eligible_log_log_slope":slope, "successive_ratios":pairs})
    diagnostics = []
    for row in rows:
        E, I = row["one_minus_overlap"], row["infidelity"]
        old = row["naive_unnormalized_one_minus_overlap"]
        static_I = row["T"] ** 4 * checks[row["case"]]["initial_BCH"][row["ordering"]]["norm"] ** 2 / (4 * row["steps"] ** 2)
        diagnostics.append({"case":row["case"],"T":row["T"],"steps":row["steps"],"ordering":row["ordering"],
            "eligible_above_floor": E > args.error_floor,
            "naive_error_absolute_difference":old-E,
            "naive_error_relative_difference":(old-E)/E if E else None,
            "naive_floor_classification_changed":(old <= args.error_floor)!=(E <= args.error_floor),
            "static_HF_BCH_predicted_infidelity":static_I,
            "static_HF_BCH_relative_error":static_I/I-1 if I else None,
            "dynamic_prediction_relative_error":row.get("prediction_relative_error"),
            "previous_relative_difference":row.get("previous_relative_difference")})
    eligible = [r for r in diagnostics if r["eligible_above_floor"]]
    prediction_summary = []
    for steps in sorted({row["steps"] for row in rows}):
        d = [r for r in eligible if r["steps"] == steps]
        p = [r for r in prediction_pairs if r["steps"] == steps]
        common_pairs = [r for r in p if r["dynamic_correct"] is not None]
        common_rows = [r for r in d if r["dynamic_prediction_relative_error"] is not None]
        prediction_summary.append({"steps":steps,
            "static_abs_relative_error_on_dynamic_subset":distribution(abs(r["static_HF_BCH_relative_error"]) for r in common_rows if r["static_HF_BCH_relative_error"] is not None),
            "common_prediction_rank_pairs":len(common_pairs),"static_correct_on_dynamic_subset":sum(r["static_HF_BCH_correct"] for r in common_pairs),
            "dynamic_abs_relative_error":distribution(abs(r["dynamic_prediction_relative_error"]) for r in d if r["dynamic_prediction_relative_error"] is not None),
            "static_abs_relative_error":distribution(abs(r["static_HF_BCH_relative_error"]) for r in d if r["static_HF_BCH_relative_error"] is not None),
            "eligible_static_rank_pairs":len(p),"correct_static_rank_pairs":sum(r["static_HF_BCH_correct"] for r in p),
            "eligible_dynamic_rank_pairs":sum(r["dynamic_correct"] is not None for r in p),"correct_dynamic_rank_pairs":sum(r["dynamic_correct"] is True for r in p)})
    reversals = []
    anchor_comparisons = defaultdict(list)
    for item in comparisons:
        if item["steps"] == args.anchor_steps and item["outcome_for_fermionic"] in ("win", "loss"):
            anchor_comparisons[item["case"], item["baseline"]].append(item)
    for (case, baseline), group in sorted(anchor_comparisons.items()):
        if {r["outcome_for_fermionic"] for r in group} == {"win", "loss"}:
            reversals.append({"case":case,"baseline":baseline,"steps":args.anchor_steps,
                              "sampled_times":[{"T":r["T"],"fermionic_outcome":r["outcome_for_fermionic"],
                                                "baseline_E_over_fermionic_E":r["baseline_E_over_fermionic_E"]}
                                               for r in sorted(group,key=lambda r:r["T"])]})
    conservation = {name: {"max_particle_number_leakage":max(r["particle_number_leakage"] for r in rows if r["ordering"]==name),
                           "max_spin_sector_leakage":max(r["spin_sector_leakage"] for r in rows if r["ordering"]==name)} for name in ORDERS}
    return {"unique_physics_rows":len(rows),"unique_cases":len(checks),"unique_molecular_families":len({c.split("_s")[0].lower() for c in checks}),
        "matched_comparison_contract":"Equal retained Hamiltonian, initial HF state, time, step count and Pauli exponential count; no postselection. Compiled gate costs are not assessed.",
        "thresholds":{"error_floor":args.error_floor,"win_tie_relative_tolerance":args.tie_rtol,"win_tie_absolute_tolerance":args.tie_atol,
                      "rank_minimum_relative_gap":args.rank_gap,"anchor_steps":args.anchor_steps},
        "baseline_comparisons":comparisons,"anchor_rows":anchor,"time_rank_reversals":reversals,
        "anchor_comparison_counts":[{"T":T,"steps":args.anchor_steps,"baseline":baseline,"counts":dict(count)} for (T,baseline),count in sorted(time_outcomes.items())],
        "scaling":scaling,"slope_summary":distribution(r["last_two_eligible_log_log_slope"] for r in scaling),
        "prediction_summary_by_steps":prediction_summary,"prediction_pair_diagnostics":prediction_pairs,
        "row_diagnostics":diagnostics,"below_floor_rows":[{k:r[k] for k in ("case","T","steps","ordering","one_minus_overlap")} for r in rows if r["one_minus_overlap"]<=args.error_floor],
        "missing_dynamic_prediction_rows":sum("predicted_infidelity" not in r for r in rows),
        "previous_comparison":{"matched_rows":sum("previous_one_minus_overlap" in r for r in rows),
           "above_floor_absolute_relative_difference":distribution(abs(r["previous_relative_difference"]) for r in rows if r.get("previous_relative_difference") is not None and min(r["one_minus_overlap"],r["previous_one_minus_overlap"])>args.error_floor)},
        "normalization":{"max_exact_norm_error":max(abs(r["exact_norm"]-1) for r in rows),"max_approximate_norm_error":max(abs(r["approx_norm"]-1) for r in rows),
           "max_infidelity_E_identity_residual":max(abs(r["infidelity"]-r["one_minus_overlap"]*(2-r["one_minus_overlap"])) for r in rows)},
        "legacy_metric_effect":{"above_floor_abs_relative_change":distribution(abs(r["naive_error_relative_difference"]) for r in eligible),
            "floor_classification_changes":sum(r["naive_floor_classification_changed"] for r in diagnostics),
            "pair_classification_changes_including_floor":sum(r["old_naive_classification_disagrees"] for r in comparisons)},
        "conservation":conservation,
        "max_owned_block_number_commutator_l1":max(c["block_commutator_bounds"]["max_N_coefficient_l1"] for c in checks.values()),
        "max_full_space_slice_residual":max(v for c in checks.values() for v in c["full_space_slice_residuals"].values()),
        "max_exact_reference_phase_aligned_distance":max(c["old_sector_reference_vs_retained_H"]["phase_aligned_distance"] for c in checks.values()),
        "max_exact_reference_phase_sensitive_distance":max(c["old_sector_reference_vs_retained_H"]["normalized_phase_sensitive_distance"] for c in checks.values()),
        "interpretation":"Descriptive audit of frozen selected Hamiltonians. Times, steps, orderings and repeated runs are not independent molecular samples. Prediction validation does not establish a cheap ordering-selection algorithm."}


def case_label(case):
    if case in CASE_LABELS:
        return CASE_LABELS[case]
    molecule, separator, scale = case.partition("_s")
    molecule = {"CH2O": r"CH$_2$O", "H2S": r"H$_2$S", "SiH4": r"SiH$_4$"}.get(molecule, molecule)
    return molecule + (" · scale " + scale if separator else "")


def plot_style():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlesize":10,
                        "axes.labelsize":9,"legend.fontsize":8,"pdf.fonttype":42,"ps.fonttype":42,
                        "axes.spines.top":False,"axes.spines.right":False,"savefig.bbox":"tight"})


def plots(rows, checks, summary, destination, args):
    plot_style();files=[]
    cases=sorted(checks,key=lambda c:(c not in ("f2","li2","nh3"),c))
    for T in sorted({r["T"] for r in rows}):
        present=[c for c in cases if any(r["case"]==c and r["T"]==T for r in rows)]
        columns=min(3,len(present));nrows=math.ceil(len(present)/columns)
        fig,axes=plt.subplots(nrows,columns,figsize=(4.2*columns,3.05*nrows),squeeze=False)
        for ax,case in zip(axes.flat,present):
            first=None
            for name in ORDERS:
                group=sorted((r for r in rows if r["case"]==case and r["T"]==T and r["ordering"]==name),key=lambda r:r["steps"])
                x=np.array([r["steps"] for r in group]);y=np.array([r["one_minus_overlap"] for r in group]);positive=y>0
                ax.loglog(x[positive],y[positive],color=COLORS[name],marker=MARKERS[name],ms=4,lw=1.5,label=LABELS[name])
                if first is None and np.any(positive):first=(x[positive],y[positive])
            if first:
                x,y=first;ax.loglog(x,y[-1]*(x[-1]/x)**2,color="0.55",ls=":",lw=1,label=r"$r^{-2}$ guide")
            positive_values=[r["one_minus_overlap"] for r in rows if r["case"]==case and r["T"]==T and r["one_minus_overlap"]>0]
            if positive_values:
                lower,upper=min(positive_values)*.6,max(positive_values)*1.5
            else:
                lower,upper=max(args.error_floor,1e-16)/10,max(args.error_floor,1e-16)*10
                ax.text(.5,.5,"All errors at numerical zero",ha="center",transform=ax.transAxes)
            ax.set_ylim(lower,upper)
            if lower <= args.error_floor <= upper:
                ax.axhline(args.error_floor,color="0.7",ls="--",lw=.7)
            ax.set_title(case_label(case))
            ax.set_xlabel("Trotter steps r");ax.set_ylabel(r"$E=1-|\langle\psi_{\rm exact}|\psi_r\rangle|$")
            ax.grid(True,which="major",alpha=.18);ax.set_xticks(sorted({r['steps'] for r in rows}))
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_xaxis().set_minor_locator(matplotlib.ticker.NullLocator())
            ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        for ax in list(axes.flat)[len(present):]:ax.set_visible(False)
        handles,labels=axes.flat[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc="upper center",ncol=4,bbox_to_anchor=(.5,1.005),frameon=False)
        fig.suptitle(f"Frozen-input state-error convergence · total time T={T:g}",y=1.035,fontsize=13)
        fig.tight_layout();stem=f"convergence_T{T:g}".replace('.','p')
        for suffix in ("pdf","png"):
            name=f"{stem}.{suffix}";fig.savefig(destination/name,dpi=300);files.append(name)
        plt.close(fig)
    f2=[r for r in rows if r['case']=='f2' and r['steps']==args.anchor_steps]
    times=sorted({r['T'] for r in f2})
    if len(times)>=2 and all(any(r['T']==T and r['ordering']==o and r['infidelity']>0 for r in f2) for T in times for o in ORDERS):
        lookup={(r['T'],r['ordering']):r for r in f2}
        ratios=[lookup[T,'jw_signed']['infidelity']/lookup[T,'fermionic_signed']['infidelity'] for T in times]
        if min(ratios)<1<max(ratios):
            fig,axes=plt.subplots(1,2,figsize=(9,3.7))
            for name in ORDERS:
                axes[0].semilogy(times,[lookup[T,name]['one_minus_overlap'] for T in times],color=COLORS[name],marker=MARKERS[name],lw=1.6,label=LABELS[name])
            axes[0].set_ylabel(r"State error $E$");axes[0].legend(frameon=False)
            axes[1].semilogy(times,ratios,color='k',marker='o',lw=1.5,label="Measured infidelity ratio")
            prediction=[lookup[T,'jw_signed'].get('predicted_infidelity')/lookup[T,'fermionic_signed']['predicted_infidelity'] if 'predicted_infidelity' in lookup[T,'jw_signed'] and lookup[T,'fermionic_signed'].get('predicted_infidelity') else math.nan for T in times]
            if all(math.isfinite(v) for v in prediction):axes[1].semilogy(times,prediction,color=COLORS['fermionic_signed'],marker='x',ls='--',lw=1,label="Transported BCH ratio")
            norms=checks['f2']['initial_BCH'];static=(norms['jw_signed']['norm']/norms['fermionic_signed']['norm'])**2
            axes[1].axhline(static,color='0.45',ls=':',label="Initial HF-BCH ratio")
            axes[1].axhline(1,color='0.6',lw=.8);axes[1].set_ylabel(r"$I_{\rm JW\ signed}/I_{\rm fermionic\ signed}$")
            axes[1].legend(frameon=False,fontsize=7.5)
            for ax in axes:ax.set_xlabel("Total evolution time T");ax.set_xticks(times);ax.grid(True,alpha=.18)
            fig.suptitle(f"F₂ ordering reversal at r={args.anchor_steps}; lines join sampled times",fontsize=11)
            fig.tight_layout()
            for suffix in ('pdf','png'):
                name=f"f2_time_reversal.{suffix}";fig.savefig(destination/name,dpi=300);files.append(name)
            plt.close(fig)
    return files


def format_value(value,precision=4):
    return "not available" if value is None else f"{value:.{precision}g}"


def report(summary,files,provenance,repeats,args):
    s=summary;lines=["# Frozen molecular ordering audit", "",
        f"Recomputed {s['unique_physics_rows']} distinct settings for {s['unique_cases']} frozen Hamiltonians in {s['unique_molecular_families']} molecular families from {len(provenance)} completed run(s).",
        "The counts below describe these inputs; repeated times, step counts and orderings are not independent molecular samples.","",
        "All three methods use the same retained Hamiltonian, initial Hartree–Fock state, total time, step count and number of Pauli exponentials. Errors use normalized full states, retain leakage, and ignore a global phase. Compiled gate costs are not compared.","",
        f"An error at or below {args.error_floor:g} is marked below the analysis floor. A win/tie uses relative tolerance {args.tie_rtol:g} and absolute tolerance {args.tie_atol:g}. Predictor ranking checks require an actual infidelity gap greater than {100*args.rank_gap:g}% and both state errors above the floor.","",
        "## Anchor comparisons", "",f"At r={args.anchor_steps}; each baseline is assessed separately.","",
        "| T | Comparator | Fermionic wins | Ties | Losses | Below floor |", "|---:|---|---:|---:|---:|---:|"]
    for row in s['anchor_comparison_counts']:
        c=row['counts'];lines.append(f"| {row['T']:g} | {LABELS[row['baseline']]} | {c.get('win',0)} | {c.get('tie',0)} | {c.get('loss',0)} | {c.get('below_floor',0)} |")
    lines.extend(["","Individual anchor errors:","","| Case | T | Fermionic signed E | JW signed E | JW magnitude E | Best within tolerance |","|---|---:|---:|---:|---:|---|"])
    for row in s['anchor_rows']:
        errors=row['errors'];best=', '.join(LABELS[o] for o in row['best_orderings_with_tolerance'])+(' (below floor)' if row['winner_below_floor'] else '')
        lines.append(f"| {row['case']} | {row['T']:g} | {errors['fermionic_signed']:.6g} | {errors['jw_signed']:.6g} | {errors['jw_magnitude']:.6g} | {best} |")
    slope=s['slope_summary'];lines.extend(["","## Convergence and error prediction","",
        f"The last two above-floor points in each case/time/ordering series give log–log slopes from {format_value(slope['min'],6)} to {format_value(slope['max'],6)} (median {format_value(slope['median'],6)}); first-order state infidelity is expected to scale as r⁻² in the asymptotic regime. Exact series and successive step-error ratios are in summary.json.","",
        "Dynamic predictions use the time-propagated projected leading defect, computed before the corresponding finite-step states. The static diagnostic uses the initial-HF BCH norm. Relative prediction errors below refer to infidelity I, not E.","",
        "| r | Dynamic prediction rows | Max relative error | Static proxy max relative error | Dynamic correct rank pairs | Static correct on same pairs |","|---:|---:|---:|---:|---:|---:|"])
    for p in s['prediction_summary_by_steps']:
        d=p['dynamic_abs_relative_error'];a=p['static_abs_relative_error_on_dynamic_subset'];lines.append(f"| {p['steps']} | {d['count']} | {format_value(100*d['max'] if d['max'] is not None else None)}% | {format_value(100*a['max'] if a['max'] is not None else None)}% | {p['correct_dynamic_rank_pairs']}/{p['eligible_dynamic_rank_pairs']} | {p['static_correct_on_dynamic_subset']}/{p['common_prediction_rank_pairs']} |")
    lines.extend(["",f"Dynamic prediction is absent for {s['missing_dynamic_prediction_rows']} rows (the runner has a matrix-size limit). Omitted rows are not counted as successful predictions. The table compares both predictors on the same available rows and rank pairs; additional static-only results remain in summary.json. Pairwise ranking comparisons share Hamiltonians and are not independent statistical trials."])
    if s["time_rank_reversals"]:
        descriptions=[r["case"] + " versus " + LABELS[r["baseline"]] for r in s["time_rank_reversals"]]
        lines.extend(["", "Anchor rankings change across sampled times for: " + "; ".join(descriptions) + ". The exact ratios and times are recorded in summary.json."])
    if any(name.startswith('f2_time_reversal') for name in files):
        lines.extend(["", "F₂ reverses its JW-signed/fermionic ranking between the sampled total times while step convergence remains regular. The initial-HF ratio is time independent and therefore cannot describe that reversal; the transported prediction can be checked against it in the figure. F₂ is a previously identified diagnostic example, not a newly prospective discovery."])
    n=s['normalization'];old=s['legacy_metric_effect'];previous=s['previous_comparison'];lines.extend(["","## Numerical accuracy and old metric effects","",
        f"Maximum state-norm errors: exact {n['max_exact_norm_error']:.3g}, approximate {n['max_approximate_norm_error']:.3g}. Maximum full-space/closed-space slice residual: {s['max_full_space_slice_residual']:.3g}. Exact-reference phase-aligned discrepancy: {s['max_exact_reference_phase_aligned_distance']:.3g}; the phase-sensitive discrepancy reaches {s['max_exact_reference_phase_sensitive_distance']:.3g} because the historical scalar-energy convention differs.","",
        f"{len(s['below_floor_rows'])} rows fall at or below the analysis floor. Replacing the stable metric by the old clipped unnormalized-overlap subtraction changes {old['floor_classification_changes']} floor labels and {old['pair_classification_changes_including_floor']} pair classifications (including floor labels). The maximum above-floor relative metric change is {format_value(old['above_floor_abs_relative_change']['max'])}. These counts use the stated tie threshold; they do not estimate floating-point uncertainty independently.","",
        f"{previous['matched_rows']} rows have saved prior frozen-input results. The largest above-floor relative difference is {format_value(previous['above_floor_absolute_relative_difference']['max'])}. Agreement with these September frozen campaigns does not repair missing input provenance in the older August aggregate.","",
        "## Symmetry and fairness", "", "| Ordering | Maximum particle-number leakage | Maximum spin-sector leakage |", "|---|---:|---:|"])
    changed=[r for r in s["baseline_comparisons"] if r["old_naive_classification_disagrees"]]
    if changed:
        details=[f"{r['case']}, T={r['T']:g}, r={r['steps']} versus {LABELS[r['baseline']]}: stable E values {r['fermionic_E']:.6g} / {r['baseline_E']:.6g}, old values {r['legacy_naive_fermionic_E']:.6g} / {r['legacy_naive_baseline_E']:.6g} (fermionic / baseline)" for r in changed]
        # Put the examples before the symmetry heading, leaving its table intact.
        at=lines.index("## Symmetry and fairness")
        lines[at:at]=["Observed classification changes:", ""]+["- "+text for text in details]+[""]
    for name,c in s['conservation'].items():lines.append(f"| {LABELS[name]} | {c['max_particle_number_leakage']:.5g} | {c['max_spin_sector_leakage']:.5g} |")
    lines.extend(["",f"Maximum owned-block number-commutator coefficient-l1 residual: {s['max_owned_block_number_commutator_l1']:.3g}. Symmetry conservation and smaller total error are distinct properties: a leaking schedule can still have less total infidelity. Sector-conditional errors in the raw data are diagnostic and are not a free postselection benefit.","",
        "The BCH pair denominator is common to permutations of each fixed Pauli list. Its relative initial-state ratio is meaningful as a cancellation diagnostic; it is not the fixed-time fidelity coefficient. The present predictions use that time-dependent coefficient explicitly.","","## Repeated-run and artifact record","",
        f"Repeated-run rows: {len(repeats)}; disagreements beyond rtol=1e-8, atol=1e-20: {sum(not r['within_repeat_tolerance'] for r in repeats)}. Repeats are not included in unique-case performance counts.","",
        "summary.json contains the thresholds, per-row diagnostics, pair comparisons, input hashes, run manifests and generator hash. normalized_results.csv is a flat export of the accepted unique rows. Source run files were read only.","","Figures (convergence panels use separate vertical ranges; a horizontal dashed line marks the analysis floor when visible):",""])
    lines.extend(f"- [{name}]({name})" for name in files if name.endswith('.pdf'))
    return '\n'.join(lines)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path,help='New directory; existing directories are refused.')
    parser.add_argument('--error-floor',type=float,default=1e-12)
    parser.add_argument('--tie-rtol',type=float,default=1e-6)
    parser.add_argument('--tie-atol',type=float,default=1e-14)
    parser.add_argument('--rank-gap',type=float,default=.01)
    parser.add_argument('--anchor-steps',type=int,default=100)
    args=parser.parse_args()
    if any(not math.isfinite(v) or v<0 for v in (args.error_floor,args.tie_rtol,args.tie_atol,args.rank_gap)) or args.anchor_steps<=0:
        parser.error('Thresholds must be finite/nonnegative and anchor steps positive.')
    if args.output.exists():parser.error('Output directory already exists; select a new directory.')
    rows,checks,provenance,repeats=load_runs(args.runs)
    result=summarize(rows,checks,args)
    args.output.mkdir(parents=True,exist_ok=False)
    files=plots(rows,checks,result,args.output,args)
    result.update({'generator_sha256':sha(__file__),'runs':provenance,'physical_checks':checks,'repeat_comparisons':repeats,'figures':files})
    dump(args.output/'summary.json',result)
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (args.output/'normalized_results.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    with (args.output/'REPORT.md').open('x') as stream:stream.write(report(result,files,provenance,repeats,args))
    print(json.dumps({'output':str(args.output),'unique_cases':len(checks),'physics_rows':len(rows),'figures':files,'slope_summary':result['slope_summary'],'anchor_comparison_counts':result['anchor_comparison_counts']},indent=2))

if __name__=='__main__':main()
