"""Bounded synthetic analyzer contracts; no molecular data or scientific imports.

Run this file with Python3.11+ or select its unittest classes. The final class
contains two SUSPECTED EXPECTED-RED source-level receipt-identity regressions:
root must observe/fix them rather than decorate them as expectedFailure.
No external study directory is read; all artifacts are tiny temporary fixtures.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import statistics
import tempfile
import unittest
from unittest import mock

_SPEC = importlib.util.spec_from_file_location(
    "prospective_analysis_contract_target", Path(__file__).with_name("analyze_study.py"))
assert _SPEC is not None and _SPEC.loader is not None
A = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(A)

CASES = ("LiH_r1p487_631g_14m4e", "LiH_r2p113_631g_14m4e",
         "BeH2_r1p369_631g_14m4e", "BeH2_r1p907_631g_14m4e")
PATHS = (("canonical", 1), ("canonical", 2), ("magnitude", 2),
         ("current", 1), ("current", 2), ("first", 1), ("first", 2),
         ("pauli", 2), ("cubic", 2), ("window_sparse", 2),
         ("window_magnitude", 2))
TIMES = (.4, 1.6)
CAPS = (8, 24)


def identity(row):
    return (row["case"], row["arm"], row["repetition"], row["formula"],
            row["time"], row["budget_multiplier"])


def raw_panel():
    rows = []
    for case in CASES:
        for rep in (1, 2):
            for arm, formula in PATHS:
                for t in TIMES:
                    for cap in CAPS:
                        rows.append(dict(case=case, arm=arm, repetition=rep,
                            formula=formula, time=t, budget_multiplier=cap,
                            path=f"{arm}_f{formula}", status="physical_checked",
                            physical_checked=True, fully_validated=True,
                            objective=dict(status="verified", accepted=True),
                            metrics=dict(projective_error=.3, infidelity=.15,
                                norm_drift=0., number_leakage=0., spin_leakage=0.),
                            constructor_CPU_seconds=2., constructor_wall_seconds=8.))
    return rows


def select(rows, *, case=CASES[0], arm="window_sparse", rep=1,
           formula=2, t=.4, cap=8):
    key = (case, arm, rep, formula, t, cap)
    return next(row for row in rows if identity(row) == key)


def comparison(items, window="window_sparse", baseline="cubic"):
    return next(item for item in items
                if item["window"] == window and item["baseline"] == baseline)


def bound_json(root, relative, value):
    path = Path(root) / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    path.write_bytes(raw)
    return dict(path=relative, sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw))


def guard_receipt(unit, cpu, wall):
    # Independent literal profile: this fixture must not inherit a weakened
    # expected profile from the analyzer under test.
    profile = dict(CPUQuota="25%", CPUWeight="5", IOWeight="10", Nice="19",
        MemoryHigh="448M", MemoryMax="512M", MemorySwapMax="0",
        MemoryAccounting="yes", CPUAccounting="yes")
    return dict(unit=unit, settings=profile, final_accounting=dict(
        CPU_seconds=cpu, elapsed_seconds=wall, CPU_accounting_available=True,
        process_CPU_timing_available=True), guard_error=None, cleanup_errors=[],
        exit_code=0, termination=None, frozen_seconds=0.,
        budget=dict(cpu_seconds=4., wall_seconds=20.), command=["synthetic-constructor"])


def partial_graph(root):
    """Authenticated relocated protocol + four explicit input refusals, no models.

    This exercises analyze's actual graph reader and complete logical roster,
    without creating fake molecular coefficients or any physical state vectors.
    """
    registry_name = "numerical_audit/prospective_window_ordering/REGISTRY.json"
    registry = dict(cases=[dict(case=case, n=14, electrons=4) for case in CASES])
    snapshot = bound_json(root, "data-a2/source_snapshot/" + registry_name, registry)
    protocol = dict(source_commit="a" * 40, sources={registry_name: snapshot["sha256"]},
        registry=registry, internal_PySCF_max_memory_MB=256,
        configuration=dict(training_time=1, training_cap_times_original_words=30,
            fresh_repetitions=2,
            arms=["canonical", "magnitude", "current", "first", "pauli", "cubic",
                  "window_sparse", "window_magnitude"],
            physical=dict(times=[.4, 1.6], cap_times_original_words=[8, 24],
                formulas={"canonical": [1, 2], "magnitude": [2], "current": [1, 2],
                    "first": [1, 2], "pauli": [2], "cubic": [2],
                    "window_sparse": [2], "window_magnitude": [2]},
                win_tie_ambiguity=1e-12)))
    pd = bound_json(root, "data-a2/INPUT_PROTOCOL_FROZEN.json", protocol)
    freeze = dict(protocol_sha256=pd["sha256"], expected_cases=4, status="partial",
        successes={}, refusals={case: dict(status="unperformed_fixture",
            reason="Synthetic input refusal; no chemistry run") for case in CASES})
    fd = bound_json(root, "data-a2/INPUTS_FROZEN.json", freeze)
    index = dict(schema="prospective-window-analysis-index-v1", mode="main",
        expected_source_commit="a" * 40, input_protocol=pd, input_freeze=fd,
        plan=None, all_orders_freeze=None, cases=[], partial_attempts=[],
        phase_receipts=[], historical_artifacts=[])
    bound_json(root, "ANALYSIS_INDEX.main.json", index)
    return index, snapshot


def save_index(root, index):
    return bound_json(root, "ANALYSIS_INDEX.main.json", index)


class TestPairedContracts(unittest.TestCase):
    def test_all_32_observations_pair_before_global_median(self):
        rows = raw_panel()
        self.assertEqual(len(rows), 352)
        # 16(-.1),8(-.5),8(+.7): median paired delta=-.1.
        # Unpaired arm medians are .15 and .20, giving a different -.05.
        assignments = [( .1, .2)] * 16 + [(.2, .7)] * 8 + [(.8, .1)] * 8
        windows, baselines, expected = [], [], {}
        k = 0
        for case in CASES:
            for rep in (1, 2):
                for t in TIMES:
                    for cap in CAPS:
                        we, be = assignments[k]
                        window = select(rows, case=case, rep=rep, t=t, cap=cap)
                        baseline = select(rows, case=case, arm="cubic", rep=rep, t=t, cap=cap)
                        window["metrics"]["projective_error"] = we
                        baseline["metrics"]["projective_error"] = be
                        window["constructor_CPU_seconds"] = 10. + k
                        baseline["constructor_CPU_seconds"] = 2. + k / 2
                        expected[(case, rep, t, cap)] = (we - be, 8. + k / 2)
                        windows.append(we); baselines.append(be); k += 1
        pairs, medians, summaries = A.paired_derivations(rows, CASES, 2)
        self.assertEqual(len(pairs), 12 * 32)
        self.assertEqual(len(medians), 12 * 16)
        self.assertEqual(len(summaries), 12)
        chosen = [p for p in pairs if p["window"] == "window_sparse" and p["baseline"] == "cubic"]
        self.assertEqual(len(chosen), 32)
        for pair in chosen:
            delta, cpu_delta = expected[(pair["case"], pair["repetition"], pair["time"], pair["budget_multiplier"])]
            self.assertAlmostEqual(pair["error_delta"], delta)
            self.assertAlmostEqual(pair["CPU_delta"], cpu_delta)
        summary = comparison(summaries)
        self.assertEqual(summary["expected_pairs"], 32)
        self.assertEqual(summary["physically_checked_pairs"], 32)
        self.assertEqual(summary["fully_validated_pairs"], 32)
        self.assertEqual((summary["win"], summary["tie"], summary["loss"], summary["censored"]), (24, 0, 8, 0))
        self.assertAlmostEqual(summary["complete_panel_median_of_paired_deltas"], -.1)
        unrelated = statistics.median(windows) - statistics.median(baselines)
        self.assertAlmostEqual(unrelated, -.05)
        self.assertNotAlmostEqual(summary["complete_panel_median_of_paired_deltas"], unrelated)
        group = next(g for g in medians if g["window"] == "window_sparse" and g["baseline"] == "cubic"
                     and g["case"] == CASES[0] and g["time"] == .4 and g["budget_multiplier"] == 8)
        self.assertEqual(group["available_repetitions"], [1, 2])
        self.assertTrue(group["complete_repetitions"])
        self.assertAlmostEqual(group["median_paired_error_delta"], -.1)

    def test_first_and_pauli_formula2_pairs_keep_counts_costs_and_censors(self):
        rows = raw_panel()
        for case in CASES:
            for rep in (1, 2):
                for t in TIMES:
                    for cap in CAPS:
                        # Formula-1 first remains a distinct reported control;
                        # its deliberately different numbers must not be paired.
                        first1 = select(rows, case=case, arm="first", rep=rep,
                                        formula=1, t=t, cap=cap)
                        first1["metrics"]["projective_error"] = .9
                        first1.update(constructor_CPU_seconds=99., constructor_wall_seconds=199.)
                        first2 = select(rows, case=case, arm="first", rep=rep, t=t, cap=cap)
                        first2["metrics"]["projective_error"] = .4
                        first2.update(constructor_CPU_seconds=5., constructor_wall_seconds=20.)
                        pauli = select(rows, case=case, arm="pauli", rep=rep, t=t, cap=cap)
                        pauli["metrics"]["projective_error"] = .1
                        pauli.update(constructor_CPU_seconds=7., constructor_wall_seconds=28.)
        missing = select(rows, arm="first", rep=2)
        missing.update(physical_checked=False, fully_validated=False, status="censored",
                       constructor_CPU_seconds=None, constructor_wall_seconds=None)
        # A measured control with an unavailable cost is not a free constructor.
        select(rows, case=CASES[1], arm="pauli")["constructor_CPU_seconds"] = None
        pairs, medians, summaries = A.paired_derivations(rows, CASES, 2)
        self.assertEqual(len(pairs), 384)
        self.assertEqual(len(medians), 192)
        self.assertEqual(len(summaries), 12)
        self.assertEqual({s["baseline"] for s in summaries},
                         {"cubic", "current", "magnitude", "canonical", "first", "pauli"})
        for window in ("window_sparse", "window_magnitude"):
            for baseline, expected_delta, expected_cpu, expected_wall in (
                    ("first", -.1, -3., -12.), ("pauli", .2, -5., -20.)):
                chosen = [p for p in pairs if p["window"] == window and p["baseline"] == baseline]
                self.assertEqual(len(chosen), 32)
                self.assertTrue(all(p["formula"] == 2 for p in chosen))
                summary = comparison(summaries, window, baseline)
                self.assertEqual(summary["expected_pairs"], 32)
                if baseline == "first":
                    self.assertEqual((summary["win"], summary["loss"], summary["censored"]), (31, 0, 1))
                    self.assertIsNone(summary["complete_panel_median_of_paired_deltas"])
                    censored = next(p for p in chosen if p["classification"] == "censored")
                    self.assertIsNone(censored["error_delta"])
                    self.assertIsNone(censored["CPU_delta"])
                    self.assertIsNone(censored["wall_delta"])
                    group = next(g for g in medians if g["window"] == window and
                                 g["baseline"] == baseline and g["case"] == CASES[0] and
                                 g["time"] == .4 and g["budget_multiplier"] == 8)
                    self.assertFalse(group["complete_repetitions"])
                    self.assertIsNone(group["median_paired_error_delta"])
                else:
                    self.assertEqual((summary["win"], summary["loss"], summary["censored"]), (0, 32, 0))
                    self.assertAlmostEqual(summary["complete_panel_median_of_paired_deltas"], .2)
                    unknown_cost = next(p for p in chosen if p["case"] == CASES[1] and
                                        p["repetition"] == 1 and p["time"] == .4 and
                                        p["budget_multiplier"] == 8)
                    self.assertIsNone(unknown_cost["CPU_delta"])
                    self.assertEqual(unknown_cost["classification"], "loss")
                known = next(p for p in chosen if p["case"] == CASES[2] and p["repetition"] == 1)
                self.assertAlmostEqual(known["error_delta"], expected_delta)
                self.assertEqual(known["CPU_delta"], expected_cpu)
                self.assertEqual(known["wall_delta"], expected_wall)

    def test_missing_repetition_and_cost_remain_null_not_free(self):
        rows = raw_panel()
        for t in TIMES:
            for cap in CAPS:
                row = select(rows, rep=2, t=t, cap=cap)
                row.update(physical_checked=False, fully_validated=False, status="censored",
                    constructor_CPU_seconds=None, constructor_wall_seconds=None)
        # An otherwise measured pair with absent cost retains its physical result.
        select(rows, case=CASES[1])["constructor_CPU_seconds"] = None
        pairs, medians, summaries = A.paired_derivations(rows, CASES, 2)
        selected = [p for p in pairs if p["window"] == "window_sparse" and p["baseline"] == "cubic"]
        missing = next(p for p in selected if p["case"] == CASES[0] and p["repetition"] == 2)
        self.assertEqual(missing["classification"], "censored")
        self.assertIsNone(missing["error_delta"])
        self.assertIsNone(missing["CPU_delta"])
        self.assertIsNone(missing["baseline_over_window_error_ratio"])
        for group in medians:
            if group["window"] == "window_sparse" and group["baseline"] == "cubic" and group["case"] == CASES[0]:
                self.assertEqual(group["available_repetitions"], [1])
                self.assertFalse(group["complete_repetitions"])
                self.assertIsNone(group["median_paired_error_delta"])
                self.assertIsNone(group["fully_validated_median_paired_error_delta"])
                self.assertEqual(group["observed_only_median_paired_error_delta"], 0.)
                self.assertIsNone(group["median_paired_CPU_delta"])
                self.assertIsNone(group["median_paired_wall_delta"])
        known = next(p for p in selected if p["case"] == CASES[1] and p["repetition"] == 1
                     and p["time"] == .4 and p["budget_multiplier"] == 8)
        self.assertEqual(known["error_delta"], 0.)
        self.assertIsNone(known["CPU_delta"])
        summary = comparison(summaries)
        self.assertEqual(summary["censored"], 4)
        self.assertEqual(summary["physically_checked_pairs"], 28)
        self.assertIsNone(summary["complete_panel_median_of_paired_deltas"])

    def test_floors_ties_losses_objective_refusals_and_censors_survive(self):
        rows = raw_panel()
        for t, cap, we, be in ((.4, 8, 0., .2), (.4, 24, .1, .1 + 5e-13),
                              (1.6, 8, .4, .2), (1.6, 24, .1, .2)):
            select(rows, t=t, cap=cap)["metrics"]["projective_error"] = we
            select(rows, arm="cubic", t=t, cap=cap)["metrics"]["projective_error"] = be
        refusal = select(rows, t=1.6, cap=24)
        refusal.update(fully_validated=False, objective=dict(status="objective_floor", accepted=False))
        select(rows, case=CASES[1]).update(physical_checked=False, fully_validated=False, status="censored")
        pairs, _, summaries = A.paired_derivations(rows, CASES, 2)
        chosen = [p for p in pairs if p["window"] == "window_sparse" and p["baseline"] == "cubic"]
        def pair(t, cap):
            return next(p for p in chosen if p["case"] == CASES[0] and p["repetition"] == 1
                        and p["time"] == t and p["budget_multiplier"] == cap)
        floor = pair(.4, 8)
        self.assertTrue(floor["window_floor"])
        self.assertEqual(floor["classification"], "win")
        self.assertIsNone(floor["baseline_over_window_error_ratio"])
        self.assertIn("floor", floor["ratio_suppression_reason"])
        self.assertEqual(pair(.4, 24)["classification"], "tie")
        self.assertEqual(pair(1.6, 8)["classification"], "loss")
        self.assertEqual(pair(1.6, 8)["baseline_over_window_error_ratio"], .5)
        self.assertEqual(pair(1.6, 24)["classification"], "win")
        self.assertFalse(pair(1.6, 24)["fully_validated"])
        self.assertIsNone(pair(1.6, 24)["baseline_over_window_error_ratio"])
        self.assertIn("Objective", pair(1.6, 24)["ratio_suppression_reason"])
        self.assertEqual(comparison(summaries)["censored"], 1)
        self.assertEqual(A.classify(-1e-12), "tie")
        self.assertEqual(A.classify(1e-12), "tie")
        self.assertEqual(A.classify(-2e-12), "win")
        self.assertEqual(A.classify(2e-12), "loss")

    def test_pareto_keeps_objective_refusals_and_separates_cost_axes(self):
        rows = raw_panel()
        group = [r for r in rows if r["case"] == CASES[0] and r["repetition"] == 1
                 and r["time"] == .4 and r["budget_multiplier"] == 8]
        self.assertEqual(len(group), 11)
        for row in group:
            row["metrics"]["projective_error"] = .9
            row.update(constructor_CPU_seconds=9., constructor_wall_seconds=9.)
        def point(arm, formula, error, cpu, wall):
            row = select(rows, arm=arm, formula=formula)
            row["metrics"]["projective_error"] = error
            row.update(constructor_CPU_seconds=cpu, constructor_wall_seconds=wall)
            return row
        point("canonical", 1, .1, 1., 10.)
        point("canonical", 2, .2, 2., 9.)
        point("current", 1, .05, 3., 3.)
        point("current", 2, .1 + 5e-13, .5, 12.)
        point("magnitude", 2, .04, None, 1.)
        point("window_sparse", 2, .03, 4., 4.).update(
            fully_validated=False, objective=dict(status="replay_refused", accepted=False))
        point("window_magnitude", 2, 2e-13, 5., 5.)
        point("pauli", 2, 0., 0., 0.).update(physical_checked=False, fully_validated=False)
        output = A.pareto_sets(rows, CASES, 2)
        self.assertEqual(len(output), 32)
        selected = next(r for r in output if r["case"] == CASES[0] and r["repetition"] == 1
                        and r["time"] == .4 and r["budget_multiplier"] == 8)
        cpu = selected["sets"]["CPU_seconds"]
        wall = selected["sets"]["wall_seconds"]
        self.assertEqual(set(cpu["pareto_paths"]),
                         {"current_f1", "current_f2", "window_sparse_f2", "window_magnitude_f2"})
        self.assertEqual(set(wall["pareto_paths"]),
                         {"magnitude_f2", "window_sparse_f2", "window_magnitude_f2"})
        self.assertIn("current_f2", cpu["dominated_by"]["canonical_f1"])
        self.assertEqual({p["path"] for p in cpu["excluded"]}, {"magnitude_f2", "pauli_f2"})
        retained = next(p for p in cpu["points"] if p["path"] == "window_sparse_f2")
        self.assertFalse(retained["fully_validated"])
        self.assertEqual(retained["objective_status"], "replay_refused")
        self.assertTrue(next(p for p in cpu["points"] if p["path"] == "window_magnitude_f2")["physical_floor"])


class TestAuthenticatedPartialGraph(unittest.TestCase):
    def test_cli_versioned_index_binds_actual_bytes_and_preserves_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            original = (root / "ANALYSIS_INDEX.main.json").read_bytes()
            relative = "indices/ANALYSIS_INDEX.main.a2.json"
            versioned = root / relative
            versioned.parent.mkdir()
            raw = (json.dumps(index, indent=2) + "\n").encode()
            versioned.write_bytes(raw)
            digest = hashlib.sha256(raw).hexdigest()
            self.assertNotEqual(digest, hashlib.sha256(original).hexdigest())
            output = root / "summary.a2.json"
            self.assertEqual(A.main([str(root), str(output), "--index", relative]), 2)
            summary = json.loads(output.read_text())
            self.assertEqual(len(summary["raw_rows"]), 352)
            self.assertTrue(all(not row["physical_checked"] for row in summary["raw_rows"]))
            self.assertFalse(summary["completion"]["complete_352_main_claim_allowed"])
            self.assertEqual(summary["analysis_index_sha256"], digest)
            self.assertEqual(summary["source_artifact_SHA256_map"][relative]["sha256"], digest)
            self.assertNotIn("ANALYSIS_INDEX.main.json", summary["source_artifact_SHA256_map"])
            self.assertEqual((root / "ANALYSIS_INDEX.main.json").read_bytes(), original)
            default_output = root / "summary.default.json"
            self.assertEqual(A.main([str(root), str(default_output)]), 2)
            self.assertEqual(json.loads(default_output.read_text())["analysis_index_sha256"],
                             hashlib.sha256(original).hexdigest())
            with self.assertRaisesRegex(ValueError, "Append-only"):
                A.main([str(root), str(output), "--index", relative])

    def test_cli_unsafe_index_refused_before_artifact_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            root = parent / "root"
            root.mkdir()
            partial_graph(root)
            outside = parent / "outside.json"
            outside.write_text("not an analysis index\n")
            (root / "escape.json").symlink_to(outside)
            output = root / "new-summary.json"
            for spelling in (str(outside), "../outside.json", "escape.json"):
                with self.subTest(index=spelling), mock.patch.object(A.Artifacts, "read_path") as read:
                    with self.assertRaisesRegex(ValueError, "Artifact .*escape"):
                        A.main([str(root), str(output), "--index", spelling])
                    read.assert_not_called()
                self.assertFalse(output.exists())

    def test_analyze_retains_352_cells_and_partial_cost_lower_bounds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            r1 = bound_json(root, "receipts/first.json", guard_receipt("fixture-call-attempt1", 1.25, 5.))
            r2 = bound_json(root, "receipts/retry.json", guard_receipt("fixture-call-attempt2", .75, 3.))
            prep = bound_json(root, "receipts/prep.json", guard_receipt("fixture-preparation", .5, 2.))
            index["partial_attempts"] = [dict(case=CASES[0], arm="window_sparse", repetition=1,
                receipt=r, record=None) for r in (r1, r2)]
            # Main's documented owner is constructors; the old synthetic
            # singular construction spelling was not a frozen campaign phase.
            index["phase_receipts"] = [dict(phase="constructors", receipt=r1), dict(phase="preparation", receipt=prep)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertEqual(summary["schema"], "prospective-window-report-data-v1")
            self.assertEqual(summary["status"], "partial_or_censored")
            self.assertEqual(summary["failures"], [])
            self.assertEqual(len(summary["raw_rows"]), 352)
            self.assertEqual(len({identity(r) for r in summary["raw_rows"]}), 352)
            self.assertTrue(all(not r["physical_checked"] for r in summary["raw_rows"]))
            self.assertTrue(all(r["constructor_CPU_seconds"] is None for r in summary["raw_rows"]))
            self.assertEqual(summary["completion"]["expected_logical_rows"], 352)
            self.assertFalse(summary["completion"]["complete_352_main_claim_allowed"])
            self.assertEqual(len(summary["paired_observations"]), 384)
            self.assertTrue(all(p["classification"] == "censored" for p in summary["paired_observations"]))
            self.assertEqual(len(summary["comparisons"]), 12)
            self.assertTrue(all(c["expected_pairs"] == 32 and c["censored"] == 32 for c in summary["comparisons"]))
            totals = summary["construction_totals"]
            self.assertEqual(totals["expected_calls"], 64)
            self.assertEqual(totals["known_calls"], 0)
            self.assertIsNone(totals["CPU_seconds"])
            self.assertIsNone(totals["wall_seconds"])
            self.assertEqual(totals["known_CPU_seconds"], 2.)
            self.assertEqual(totals["known_wall_seconds"], 8.)
            self.assertEqual(summary["construction_by_arm"]["window_sparse"]["indexed_calls"], 1)
            self.assertIsNone(summary["construction_by_arm"]["window_sparse"]["CPU_seconds"])
            call = summary["constructor_costs"][0]
            self.assertFalse(call["accounting_complete"])
            self.assertTrue(call["provided_receipts_accounting_complete"])
            self.assertEqual(len(call["attempts"]), 2)
            self.assertTrue(summary["other_phase_costs"][0]["already_counted_in_construction"])
            self.assertEqual(summary["other_phase_totals"]["constructors"]["indexed_attempts"], 0)
            self.assertEqual(summary["other_phase_totals"]["constructors"]["CPU_seconds"], 0)
            self.assertEqual(summary["other_phase_totals"]["preparation"]["CPU_seconds"], .5)
            source = "data-a2/source_snapshot/numerical_audit/prospective_window_ordering/REGISTRY.json"
            self.assertIn(source, summary["source_artifact_SHA256_map"])
            self.assertEqual(summary["scientific_source"]["commit"], "a" * 40)
            self.assertLess(summary["bytes_hashed"], 100_000)
            output = root / "new-summary.json"
            A.write_new(output, summary)
            self.assertEqual(json.loads(output.read_text())["completion"]["expected_logical_rows"], 352)

    def test_tampered_source_or_descriptor_path_quarantines_full_roster(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, snapshot = partial_graph(root)
            (root / snapshot["path"]).write_text("tampered source snapshot\n")
            summary = A.analyze(root, "main")
            self.assertEqual(summary["status"], "partial_or_censored")
            self.assertEqual(len(summary["raw_rows"]), 352)
            self.assertTrue(summary["failures"])
            self.assertIn("hash differs", summary["failures"][0]["message"])
            self.assertTrue(all(not row["physical_checked"] for row in summary["raw_rows"]))
            self.assertFalse(summary["completion"]["complete_352_main_claim_allowed"])
            index["input_protocol"]["path"] = "../outside-protocol.json"
            save_index(root, index)
            escaped = A.analyze(root, "main")
            self.assertEqual(len(escaped["raw_rows"]), 352)
            self.assertIn("escape", escaped["failures"][0]["message"])

    def test_duplicate_same_hash_partial_attempt_rejected_without_invented_cost(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            rd = bound_json(root, "receipt.json", guard_receipt("one-unit", 1.25, 5.))
            item = dict(case=CASES[0], arm="cubic", repetition=1, receipt=rd, record=None)
            index["partial_attempts"] = [item, copy.deepcopy(item)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertTrue(any("receipt repeated" in f["message"] for f in summary["failures"]))
            self.assertEqual(summary["construction_totals"]["known_CPU_seconds"], 1.25)
            self.assertEqual(summary["construction_totals"]["known_wall_seconds"], 5.)
            self.assertEqual(sum(f.get("phase")=="partial_attempts" for f in summary["failures"]), 1)
            self.assertIsNone(summary["construction_totals"]["CPU_seconds"])


class TestPrimarySourceBinding(unittest.TestCase):
    """Tiny inventory contracts; no fake complete physics/constructor graph.

    Protocol-source authentication is exercised elsewhere. These tests pin the
    primary driver's exact runtime dependency inventory and frozen-value check.
    """

    def fixture(self):
        # Independent literal logical paths, not a production inventory constant.
        sources = {
            "numerical_audit/prospective_window_ordering/physical.jl": "a" * 64,
            "numerical_audit/prospective_window_ordering/src/WindowPhysical.jl": "b" * 64,
            "numerical_audit/prospective_window_ordering/src/WindowStudyIO.jl": "c" * 64,
            "numerical_audit/julia_cubic_feasibility/src/ResponseCore.jl": "d" * 64,
        }
        hashes = {"physical.jl": "a" * 64, "src/WindowPhysical.jl": "b" * 64,
                  "src/WindowStudyIO.jl": "c" * 64,
                  "../julia_cubic_feasibility/src/ResponseCore.jl": "d" * 64}
        return dict(source_hashes=hashes), dict(sources=sources)

    def test_exact_primary_inventory_matches_frozen_protocol_without_mutation(self):
        manifest, protocol = self.fixture()
        before = copy.deepcopy((manifest, protocol))
        A.primary_source_binding(manifest, protocol)
        self.assertEqual((manifest, protocol), before)

    def test_missing_or_extra_primary_dependency_rejected(self):
        manifest, protocol = self.fixture()
        for name in tuple(manifest["source_hashes"]):
            changed = copy.deepcopy(manifest)
            del changed["source_hashes"][name]
            with self.subTest(missing=name), self.assertRaisesRegex(ValueError, "Complete primary"):
                A.primary_source_binding(changed, protocol)
        changed = copy.deepcopy(manifest)
        changed["source_hashes"]["unused.jl"] = "e" * 64
        with self.assertRaisesRegex(ValueError, "Complete primary"):
            A.primary_source_binding(changed, protocol)
        for hashes in (None, [], {}):
            with self.subTest(inventory=hashes), self.assertRaisesRegex(ValueError, "Complete primary"):
                A.primary_source_binding(dict(source_hashes=hashes), protocol)
        with self.assertRaisesRegex(ValueError, "Complete primary"):
            A.primary_source_binding({}, protocol)

    def test_changed_driver_or_dependency_hash_cannot_inherit_input_order_identity(self):
        manifest, protocol = self.fixture()
        for name in manifest["source_hashes"]:
            changed = copy.deepcopy(manifest)
            changed["source_hashes"][name] = "f" * 64
            with self.subTest(changed=name), self.assertRaisesRegex(ValueError, "scientific source differs"):
                A.primary_source_binding(changed, protocol)
        changed = copy.deepcopy(manifest)
        changed["source_hashes"]["physical.jl"] = None
        del protocol["sources"]["numerical_audit/prospective_window_ordering/physical.jl"]
        # Missing expected and claimed values cannot authenticate as None==None.
        with self.assertRaisesRegex(ValueError, "scientific source differs"):
            A.primary_source_binding(changed, protocol)


class TestArtifactAndOutputContracts(unittest.TestCase):
    def test_hash_size_mutation_and_lexical_path_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            descriptor = bound_json(root, "artifact.json", dict(value=1))
            wrong = dict(descriptor, sha256="0" * 64)
            with self.assertRaisesRegex(ValueError, "hash differs"):
                A.Artifacts(root).read(wrong)
            with self.assertRaisesRegex(ValueError, "byte-size"):
                A.Artifacts(root).read(dict(descriptor, bytes=descriptor["bytes"] + 1))
            for spelling in ("../artifact.json", "/tmp/outside.json", "dir\\escape.json"):
                with self.assertRaisesRegex(ValueError, "escape"):
                    A.Artifacts(root).path(spelling)
            reader = A.Artifacts(root)
            reader.read(descriptor)
            (root / descriptor["path"]).write_text('{"value":2}\n')
            with self.assertRaisesRegex(ValueError, "changed after authentication"):
                reader.read(descriptor)

    def test_symlink_escape_and_append_only_publication(self):
        with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as outside:
            root = Path(tmp)
            foreign = Path(outside) / "foreign.json"
            foreign.write_text("{}\n")
            (root / "link.json").symlink_to(foreign)
            with self.assertRaisesRegex(ValueError, "symlink escape"):
                A.Artifacts(root).path("link.json")
            output = root / "summary.json"
            A.write_new(output, dict(status="synthetic"))
            original = output.read_bytes()
            self.assertFalse(output.with_name("summary.json.partial").exists())
            with self.assertRaisesRegex(ValueError, "Append-only"):
                A.write_new(output, dict(status="replacement"))
            self.assertEqual(output.read_bytes(), original)
            other = root / "other.json"
            partial = root / "other.json.partial"
            partial.write_bytes(b"retained unfinished output\n")
            with self.assertRaisesRegex(ValueError, "Append-only"):
                A.write_new(other, dict(status="new"))
            self.assertFalse(other.exists())
            self.assertEqual(partial.read_bytes(), b"retained unfinished output\n")

    def test_duplicate_json_keys_and_nonfinite_constants_rejected(self):
        for raw in ('{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                A.decode_json(raw)


class TestSuspectedReceiptIdentity(unittest.TestCase):
    """SUSPECTED EXPECTED-RED until root observes and fixes analyzer source.

    Distinct hashes do not make one guard service two independent attempts.
    Failure/quarantine is acceptable; summing the duplicated unit is not.
    """
    def test_suspected_red_same_unit_different_hash_partial_cost_not_double_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            first = bound_json(root, "first.json", guard_receipt("reused-unit", 1.25, 5.))
            second = bound_json(root, "second.json", guard_receipt("reused-unit", .75, 3.))
            self.assertNotEqual(first["sha256"], second["sha256"])
            index["partial_attempts"] = [dict(case=CASES[0], arm="cubic", repetition=1,
                receipt=rd, record=None) for rd in (first, second)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertLessEqual(summary["construction_totals"]["known_CPU_seconds"], 1.25,
                "One guard unit must not become two charged partial attempts solely through distinct receipt hashes")
            self.assertTrue(any("unit" in str(f).lower() for f in summary["failures"]),
                "Duplicate unit identity must be explicitly rejected/quarantined")

    def test_suspected_red_same_unit_different_hash_phase_cost_not_double_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            first = bound_json(root, "first.json", guard_receipt("reused-phase-unit", 1.25, 5.))
            second = bound_json(root, "second.json", guard_receipt("reused-phase-unit", .75, 3.))
            self.assertNotEqual(first["sha256"], second["sha256"])
            index["phase_receipts"] = [dict(phase="preparation", receipt=rd) for rd in (first, second)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertLessEqual(summary["other_phase_totals"]["preparation"]["known_CPU_seconds"], 1.25,
                "One guard unit must not become two phase attempts solely through distinct receipt hashes")
            self.assertTrue(any("unit" in str(f).lower() for f in summary["failures"]),
                "Duplicate phase unit identity must be explicitly rejected/quarantined")


class TestReceiptIdentityFix(unittest.TestCase):
    def test_partial_alias_across_calls_preserves_first_valid_known_cost(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            first = bound_json(root, "first.json", guard_receipt("same-constructor-unit", 1.25, 5.))
            second = bound_json(root, "alias.json", guard_receipt("same-constructor-unit", .75, 3.))
            index["partial_attempts"] = [
                dict(case=CASES[0], arm="cubic", repetition=1, receipt=first, record=None),
                dict(case=CASES[1], arm="current", repetition=2, receipt=second, record=None)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertEqual(summary["status"], "partial_or_censored")
            self.assertFalse(summary["completion"]["complete_352_main_claim_allowed"])
            totals = summary["construction_totals"]
            self.assertEqual(totals["known_CPU_seconds"], 1.25)
            self.assertEqual(totals["known_wall_seconds"], 5.)
            self.assertIsNone(totals["CPU_seconds"])
            self.assertIsNone(totals["wall_seconds"])
            self.assertEqual(len(summary["constructor_costs"]), 1)
            self.assertEqual(summary["constructor_costs"][0]["arm"], "cubic")
            rejected = next(f for f in summary["failures"] if f.get("phase") == "partial_attempts")
            self.assertTrue(rejected["quarantined"])
            self.assertEqual(rejected["measurement"]["CPU_seconds"], .75)
            self.assertEqual(rejected["receipt"], second)
            self.assertEqual(rejected["first_receipt_sha256"], first["sha256"])

    def test_construction_to_phase_changed_hash_alias_is_retained_but_not_charged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            first = bound_json(root, "construction.json", guard_receipt("construction-phase-unit", 1.25, 5.))
            alias = bound_json(root, "phase-alias.json", guard_receipt("construction-phase-unit", .75, 3.))
            index["partial_attempts"] = [dict(case=CASES[0], arm="cubic", repetition=1,
                receipt=first, record=None)]
            index["phase_receipts"] = [dict(phase="preparation", receipt=alias)]
            save_index(root, index)
            summary = A.analyze(root, "main")
            self.assertEqual(summary["construction_totals"]["known_CPU_seconds"], 1.25)
            self.assertEqual(summary["construction_totals"]["known_wall_seconds"], 5.)
            phase = summary["other_phase_totals"]["preparation"]
            self.assertFalse(phase["accounting_complete"])
            self.assertEqual(phase["quarantined_attempts"], 1)
            self.assertEqual(phase["indexed_attempts"], 0)
            self.assertIsNone(phase["CPU_seconds"])
            self.assertIsNone(phase["wall_seconds"])
            self.assertEqual(phase["known_CPU_seconds"], 0.)
            retained = summary["other_phase_costs"][0]
            self.assertTrue(retained["quarantined"])
            self.assertTrue(retained["already_counted_in_construction"])
            self.assertEqual(retained["measurement"]["CPU_seconds"], .75)
            self.assertEqual(retained["measurement"]["descriptor"], alias)
            self.assertTrue(any("aliases construction" in f["message"] for f in summary["failures"]))
            self.assertFalse(summary["completion"]["complete_352_main_claim_allowed"])

    def test_missing_guard_unit_is_not_a_fictitious_distinct_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for i, value in enumerate((None, "", 7)):
                invalid = guard_receipt("temporary", 1., 4.)
                invalid["unit"] = value
                rd = bound_json(root, f"missing-unit-{i}.json", invalid)
                with self.subTest(unit=value), self.assertRaisesRegex(ValueError, "unit identity"):
                    A.receipt(A.Artifacts(root), rd)


class PreflightRelocationContracts(unittest.TestCase):
    """Real hash-bound input graph; no chemistry or physical outcomes invented."""

    def graph(self, root):
        index, _ = partial_graph(root)
        index['mode'] = 'preflight'
        protocol = json.loads((root / index['input_protocol']['path']).read_text())
        altered = dict(expected_cases=1, preflight_only=True, status='complete',
            parent_inventory_sha256=index['input_freeze']['sha256'],
            protocol_sha256=index['input_protocol']['sha256'],
            source_protocol_path=str(root / index['input_protocol']['path']),
            refusals={}, successes={CASES[0] + '_preflight_7over5': {}})
        index['preflight_input_freeze'] = bound_json(root, 'preflight/INPUTS_FROZEN.json', altered)
        plan = dict(schema='prospective-window-constructor-plan-v1', mode='preflight',
            root=str(root), source_freeze=index['input_protocol'],
            source_snapshot=protocol['sources'], physical_error_read=False,
            input_freeze=index['preflight_input_freeze'])
        index['plan'] = bound_json(root, 'preflight-plan.json', plan)
        return index, altered, plan

    def read_graph(self, root, index):
        summary = {}
        result = A.source_and_inputs(A.Artifacts(root), index, 'preflight', summary)
        # The location check must succeed without upgrading four retained input
        # refusals into scientific completeness or reading absent model payloads.
        self.assertEqual(result[3], [CASES[0] + '_preflight_7over5'])
        self.assertFalse(result[4])
        self.assertFalse(summary['inputs']['original_complete'])

    def test_same_location_and_moved_archive_share_authenticated_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            original = Path(tmp) / 'original'
            index, _, _ = self.graph(original)
            self.read_graph(original, index)
            # Preserve the old same-location, no-plan partial-evidence behavior.
            self.read_graph(original, dict(index, plan=None))
            before = {p.relative_to(original): p.read_bytes() for p in original.rglob('*') if p.is_file()}
            relocated = Path(tmp) / 'published-copy'
            shutil.move(str(original), relocated)
            self.assertFalse(original.exists())
            self.read_graph(relocated, index)
            after = {p.relative_to(relocated): p.read_bytes() for p in relocated.rglob('*') if p.is_file()}
            self.assertEqual(before, after)
            with self.assertRaisesRegex(ValueError, 'needs an authenticated plan'):
                self.read_graph(relocated, dict(index, plan=None))

    def test_rehashed_forged_original_root_is_not_corroborated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'campaign'
            index, _, plan = self.graph(root)
            plan['root'] = '/unrelated/original-root'
            index['plan'] = bound_json(root, 'preflight-plan.json', plan)
            with self.assertRaisesRegex(ValueError, 'differs from plan root/source'):
                self.read_graph(root, index)

    def test_rehashed_forged_protocol_path_is_not_corroborated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'campaign'
            index, altered, plan = self.graph(root)
            altered['source_protocol_path'] = str(root / 'wrong/INPUT_PROTOCOL_FROZEN.json')
            index['preflight_input_freeze'] = bound_json(root, 'preflight/INPUTS_FROZEN.json', altered)
            plan['input_freeze'] = index['preflight_input_freeze']
            index['plan'] = bound_json(root, 'preflight-plan.json', plan)
            with self.assertRaisesRegex(ValueError, 'differs from plan root/source'):
                self.read_graph(root, index)

    def test_plan_hash_and_source_binding_are_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'campaign'
            index, _, plan = self.graph(root)
            plan['root'] = '/different-root'
            bound_json(root, 'preflight-plan.json', plan)  # Intentionally leave stale descriptor.
            with self.assertRaisesRegex(ValueError, 'Artifact hash differs'):
                self.read_graph(root, index)
            plan['root'] = str(root)
            plan['source_freeze'] = dict(index['input_protocol'], path='../outside.json')
            index['plan'] = bound_json(root, 'preflight-plan.json', plan)
            with self.assertRaisesRegex(ValueError, 'plan/source binding'):
                self.read_graph(root, index)

    def test_recorded_path_must_be_absolute_and_canonical(self):
        for spelling in ('relative/root', '/root/../elsewhere', '//host/root', '/root\\escape'):
            with self.subTest(path=spelling), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp) / 'campaign'
                index, _, plan = self.graph(root)
                plan['root'] = spelling
                index['plan'] = bound_json(root, 'preflight-plan.json', plan)
                with self.assertRaisesRegex(ValueError, 'Canonical recorded absolute path'):
                    self.read_graph(root, index)


def constructor_metadata(key, source="a" * 64, model="b" * 64):
    """Literal authenticated call metadata only, not a molecular result graph."""
    return dict(schema="prospective-window-constructor-v1", artifact_complete=True,
        case=key[0], arm=key[1], repetition=key[2], source_freeze_sha256=source,
        model_sha256=model, physical_error_read=False,
        score_kind="diagonal_reference_cubic_response",
        constructor=dict(order=[0, 1], score=.25, status="complete",
            stats=dict(window_status="completed")))


def panel_receipt_fixture(root, case, shared=None):
    """Eight cost-only censored ledgers; no candidate/physical completeness."""
    entries, retained = [], []
    for arm in ("canonical", "magnitude", "current", "first", "pauli", "cubic",
                "window_sparse", "window_magnitude"):
        rd = bound_json(root, f"{case}-{arm}.json",
                        guard_receipt(f"{case}-{arm}-unit", .25, 1.))
        attempts = [dict(receipt=rd, record=None)]
        cpu, wall = .25, 1.
        if shared is not None and arm == "cubic":
            attempts.insert(0, dict(receipt=shared, record=None))
            cpu += 1.25; wall += 5.
        entries.append(dict(arm=arm, repetition=1, attempts=attempts, selected_attempt=None))
        retained.append(dict(arm=arm, repetition=1, status="censored", constructor=None,
            selected_record=None, attempts=[dict(receipt=a["receipt"], constructor=None) for a in attempts],
            cost_complete=True, guard_CPU_seconds=cpu, guard_wall_seconds=wall,
            within_total_construction_ceiling=True))
    return dict(constructors=retained), dict(constructors=entries)


class TestCampaignReceiptOwnership(unittest.TestCase):
    def test_exact_hash_duplicate_within_and_across_partial_calls_keeps_prefix(self):
        for cross_call in (False, True):
            with self.subTest(cross_call=cross_call), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp); index, _ = partial_graph(root)
                rd = bound_json(root, "first.json", guard_receipt("prefix-unit", 1.25, 5.))
                alias = bound_json(root, "copy.json", guard_receipt("prefix-unit", 1.25, 5.))
                tail = bound_json(root, "tail.json", guard_receipt("tail-unit", .75, 3.))
                self.assertEqual(rd["sha256"], alias["sha256"])
                index["partial_attempts"] = [
                    dict(case=CASES[0], arm="cubic", repetition=1, receipt=rd, record=None),
                    dict(case=CASES[1] if cross_call else CASES[0], arm="cubic",
                         repetition=1, receipt=alias, record=None),
                    "malformed row",
                    dict(case=CASES[0], arm="current", repetition=2, receipt=tail, record=None)]
                save_index(root, index); result = A.analyze(root, "main")
                self.assertEqual(result["construction_totals"]["known_CPU_seconds"], 2.)
                self.assertEqual(result["construction_totals"]["known_wall_seconds"], 8.)
                self.assertIsNone(result["construction_totals"]["CPU_seconds"])
                errors = [f for f in result["failures"] if f.get("phase") == "partial_attempts"]
                self.assertEqual(len(errors), 2)
                self.assertTrue(all(f["quarantined"] for f in errors))
                self.assertIn("fully_validated_rows", result["completion"])
                self.assertNotIn("source_or_artifact_validation_failed", result["completion"])
                self.assertEqual(len(result["raw_rows"]), 352)
                self.assertEqual(len(result["constructor_costs"]), 2)

    def test_own_panel_and_global_reference_is_one_attempt_not_a_duplicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reader = A.Artifacts(root)
            rd = bound_json(root, "own.json", guard_receipt("own-unit", 1.25, 5.))
            measured = A.receipt(reader, rd); key = (CASES[0], "cubic", 1)
            ctor = dict(status="complete", score=.25, score_kind="diagonal_reference_cubic_response",
                        order=[0, 1], window_status="completed")
            existing = {key: dict(case=key[0], arm=key[1], repetition=1, status="retained_order",
                constructor=ctor, attempts=[measured], accounting_complete=True,
                CPU_seconds=1.25, wall_seconds=5., known_CPU_seconds=1.25, known_wall_seconds=5.,
                source="PANEL complete per-call attempt ledger")}
            global_map = {key: dict(receipt=measured, record=constructor_metadata(key))}
            row = dict(case=key[0], arm=key[1], repetition=1, receipt=rd, record=None)
            errors = []
            A.partial_costs(reader, dict(partial_attempts=[row]), CASES, 2,
                            global_map, existing, errors)
            self.assertEqual(errors, [])
            self.assertEqual(existing[key]["constructor"], ctor)
            self.assertEqual(existing[key]["known_CPU_seconds"], 1.25)
            self.assertEqual(len(existing[key]["attempts"]), 1)

    def test_partial_retry_union_preserves_distinct_global_final_and_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reader = A.Artifacts(root); key = (CASES[0], "cubic", 1)
            first = bound_json(root, "retry.json", guard_receipt("retry-unit", 1.25, 5.))
            last = bound_json(root, "final.json", guard_receipt("final-unit", .75, 3.))
            global_map = {key: dict(receipt=A.receipt(reader, last), record=constructor_metadata(key))}
            rows = [dict(case=key[0], arm=key[1], repetition=1, receipt=rd, record=None)
                    for rd in (first, last, last)]
            costs, errors = {}, []
            A.partial_costs(reader, dict(partial_attempts=rows), CASES, 2, global_map, costs, errors)
            cost = costs[key]
            self.assertEqual(cost["known_CPU_seconds"], 2.)
            self.assertEqual(cost["known_wall_seconds"], 8.)
            self.assertEqual({a["descriptor"]["sha256"] for a in cost["attempts"]},
                             {first["sha256"], last["sha256"]})
            self.assertEqual(cost["constructor"]["score"], .25)
            self.assertEqual(cost["constructor"]["window_status"], "completed")
            self.assertIsNone(cost["CPU_seconds"])
            self.assertFalse(cost["accounting_complete"])
            self.assertEqual(len(errors), 1)

    def test_partial_constructor_must_bind_source_model_call_and_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); key = (CASES[0], "cubic", 1)
            rd = bound_json(root, "guard.json", guard_receipt("identity-unit", 1.25, 5.))
            good = constructor_metadata(key)
            variants = [("valid", good)]
            for field, value in (("schema", "wrong"), ("artifact_complete", False),
                    ("case", CASES[1]), ("source_freeze_sha256", "c" * 64),
                    ("model_sha256", "d" * 64), ("physical_error_read", True)):
                bad = copy.deepcopy(good); bad[field] = value; variants.append((field, bad))
            for name, record in variants:
                with self.subTest(field=name):
                    desc = bound_json(root, name + ".json", record)
                    index = dict(input_protocol=dict(sha256="a" * 64),
                        partial_attempts=[dict(case=key[0], arm=key[1], repetition=1, receipt=rd, record=desc)])
                    costs, errors = {}, []
                    A.partial_costs(A.Artifacts(root), index, CASES, 2, {}, costs, errors,
                                    models={key[0]: dict(json_hash="b" * 64)})
                    if name == "valid":
                        self.assertEqual(errors, [])
                        self.assertEqual(costs[key]["known_CPU_seconds"], 1.25)
                    else:
                        self.assertEqual(costs, {})
                        self.assertEqual(len(errors), 1)
                        self.assertTrue(errors[0]["quarantined"])
                        self.assertEqual(errors[0]["measurement"]["CPU_seconds"], 1.25)

    def test_cross_case_nonfinal_panel_alias_is_rejected_and_first_costs_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reader = A.Artifacts(root); owners = A.ReceiptOwners()
            rd = bound_json(root, "shared.json", guard_receipt("shared-nonfinal-unit", 1.25, 5.))
            panel, inventory = panel_receipt_fixture(root, CASES[0], rd)
            first = {}
            A.panel_costs(reader, panel, inventory, CASES[0], 1, {}, owners=owners, retained=first)
            panel, inventory = panel_receipt_fixture(root, CASES[1], rd)
            second = {}
            with self.assertRaisesRegex(ValueError, "receipt/unit reused"):
                A.panel_costs(reader, panel, inventory, CASES[1], 1, {}, owners=owners, retained=second)
            costs = {(CASES[0], *key): value for key, value in first.items()}
            costs.update({(CASES[1], *key): value for key, value in second.items()})
            errors = []
            A.partial_costs(reader, {}, CASES, 1, {}, costs, errors, owners=owners)
            self.assertEqual(sum(c["known_CPU_seconds"] for c in costs.values()), 3.5)
            self.assertEqual(sum(a["descriptor"]["sha256"] == rd["sha256"]
                                 for c in costs.values() for a in c["attempts"]), 1)

    def test_global_final_owner_blocks_other_case_nonfinal_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reader = A.Artifacts(root); owners = A.ReceiptOwners()
            rd = bound_json(root, "global.json", guard_receipt("global-final-unit", 1.25, 5.))
            key = (CASES[1], "current", 1); measured = A.receipt(reader, rd)
            owners.claim(key, measured)
            panel, inventory = panel_receipt_fixture(root, CASES[0], rd)
            with self.assertRaisesRegex(ValueError, "receipt/unit reused"):
                A.panel_costs(reader, panel, inventory, CASES[0], 1, {}, owners=owners)
            self.assertEqual(owners.lookup(measured)["key"], key)
            self.assertEqual(len(owners.calls[key]), 1)

    def test_phase_duplicates_invalid_receipts_and_nullable_tombstones(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); index, _ = partial_graph(root)
            valid = bound_json(root, "prep.json", guard_receipt("prep-unit", 1.25, 5.))
            bad = guard_receipt("invalid-unit", .75, 3.); bad["settings"]["MemoryMax"] = "1G"
            invalid = bound_json(root, "invalid.json", bad)
            index["phase_receipts"] = [
                dict(phase="preparation", receipt=valid),
                dict(phase="preparation", receipt=valid),
                dict(phase="independent", receipt=invalid),
                dict(phase="physical", receipt=dict(path="missing.json", sha256="a" * 64))]
            save_index(root, index); result = A.analyze(root, "main")
            prep = result["other_phase_totals"]["preparation"]
            self.assertEqual(prep["known_CPU_seconds"], 1.25)
            self.assertEqual(prep["quarantined_attempts"], 1)
            self.assertFalse(prep["accounting_complete"])
            self.assertIsNone(prep["CPU_seconds"])
            for phase in ("physical", "independent"):
                total = result["other_phase_totals"][phase]
                self.assertEqual(total["quarantined_attempts"], 1)
                self.assertEqual(total["indexed_attempts"], 0)
                self.assertFalse(total["accounting_complete"])
                self.assertIsNone(total["CPU_seconds"])
                self.assertIsNone(next(p for p in result["other_phase_costs"] if p["phase"] == phase)["measurement"])

    def test_exact_constructor_receipt_only_deduplicates_in_owning_phase(self):
        for wrong in ("physical", "packaging", "construction"):
            with self.subTest(phase=wrong), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp); index, _ = partial_graph(root)
                rd = bound_json(root, "construction.json", guard_receipt("owner-unit", 1.25, 5.))
                index["partial_attempts"] = [dict(case=CASES[0], arm="cubic", repetition=1, receipt=rd, record=None)]
                index["phase_receipts"] = [dict(phase=wrong, receipt=rd), dict(phase="constructors", receipt=rd)]
                save_index(root, index); result = A.analyze(root, "main")
                self.assertEqual(result["construction_totals"]["known_CPU_seconds"], 1.25)
                self.assertFalse(result["other_phase_totals"][wrong]["accounting_complete"])
                self.assertIsNone(result["other_phase_totals"][wrong]["CPU_seconds"])
                own = result["other_phase_totals"]["constructors"]
                self.assertTrue(own["accounting_complete"])
                self.assertEqual(own["already_counted_in_construction_attempts"], 1)
                self.assertEqual(own["CPU_seconds"], 0.)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reader = A.Artifacts(root); owners = A.ReceiptOwners()
            rd = bound_json(root, "preflight.json", guard_receipt("preflight-owner", 1.25, 5.))
            owners.claim((CASES[0], "cubic", 1), A.receipt(reader, rd)); failures = []
            phase = A.phase_costs(reader, dict(phase_receipts=[dict(phase="preflight", receipt=rd)]),
                                  owners, "preflight", failures)
            self.assertEqual(failures, [])
            self.assertTrue(phase[0]["already_counted_in_construction"])
            self.assertFalse(phase[0]["quarantined"])


class TestObjectiveNumbersAndCaseFlags(unittest.TestCase):
    def fixture(self, actual=1.000000002, claimed=1.):
        # Existing producer compare_scores contract, stated independently in
        # literals. This fixture asserts no molecular/scorer result.
        difference = abs(actual - claimed); scale = max(actual, claimed)
        tolerance = 1e-24 + 1e-8 * scale
        within, floor = difference <= tolerance, scale <= 1e-24
        status = "mismatch" if not within else "floor_censored" if floor else "verified"
        candidate = dict(name="literal", partition="conserving", order=[1, 2],
            constructor_scored=True, constructor_score=claimed,
            score_kind="diagonal_reference_cubic_response")
        logical = dict(arm="cubic", repetition=1, formula=2, time=.4, budget_multiplier=8,
            name="literal", partition="conserving", candidate_index=1, check_index=1, status=status)
        record = dict(check_index=1, partition="conserving", order=[1, 2],
            score_kind=candidate["score_kind"], constructor_score=claimed, replayed_score=actual,
            absolute_difference=difference, relative_difference=difference / scale if scale else 0.,
            tolerance=tolerance, numerical_floor=floor, within_tolerance=within,
            verified=within and not floor, status=status)
        check = dict(relative_tolerance=1e-8, absolute_tolerance=1e-24, logical_rows=[logical], checks=[record])
        return check, {("cubic", 1, 2, .4, 8): candidate}

    def test_numeric_rule_verified_mismatch_and_boolean_floor(self):
        for actual, claimed, accepted, status in (
                (1.000000002, 1., True, "verified"), (1.1, 1., False, "mismatch"),
                (0., 0., False, "floor_censored"), (5e-25, 1e-25, False, "floor_censored")):
            with self.subTest(actual=actual, claimed=claimed):
                check, candidates = self.fixture(actual, claimed)
                row = next(iter(A.objective_rows(check, candidates, 2).values()))
                self.assertIs(row["accepted"], accepted)
                self.assertEqual(row["status"], status)
                self.assertIs(type(row["numerical_floor"]), bool)

    def test_missing_nonfinite_negative_and_forged_comparison_fields_fail(self):
        fields = ("replayed_score", "absolute_difference", "relative_difference", "tolerance",
                  "numerical_floor", "verified", "within_tolerance")
        for field in fields:
            check, candidates = self.fixture()
            del check["checks"][0][field]
            with self.subTest(missing=field), self.assertRaises((ValueError, KeyError)):
                A.objective_rows(check, candidates, 2)
        for field in ("replayed_score", "constructor_score", "absolute_difference", "relative_difference", "tolerance"):
            for value in (float("nan"), float("inf"), -1.):
                check, candidates = self.fixture()
                check["checks"][0][field] = value
                with self.subTest(field=field, value=value), self.assertRaises((ValueError, KeyError)):
                    A.objective_rows(check, candidates, 2)
        for field, value in (("absolute_difference", 0.), ("relative_difference", 0.),
                ("tolerance", 1.), ("numerical_floor", True), ("numerical_floor", 0),
                ("verified", False), ("within_tolerance", False)):
            check, candidates = self.fixture(); check["checks"][0][field] = value
            with self.subTest(forged=field, value=value), self.assertRaises(ValueError):
                A.objective_rows(check, candidates, 2)
        check, candidates = self.fixture(1.1, 1.)
        check["checks"][0].update(status="verified", verified=True, within_tolerance=True)
        check["logical_rows"][0]["status"] = "verified"
        with self.assertRaises(ValueError): A.objective_rows(check, candidates, 2)
        check, candidates = self.fixture(0., 0.)
        check["checks"][0].update(numerical_floor=False, verified=True, status="verified")
        check["logical_rows"][0]["status"] = "verified"
        with self.assertRaises(ValueError): A.objective_rows(check, candidates, 2)

    def test_refused_objective_stays_censored_without_invented_numbers(self):
        check, candidates = self.fixture()
        record = check["checks"][0]
        for field in ("replayed_score", "absolute_difference", "relative_difference", "tolerance", "numerical_floor"):
            del record[field]
        record.update(status="numerical_refusal", verified=False, within_tolerance=False)
        check["logical_rows"][0]["status"] = "numerical_refusal"
        output = next(iter(A.objective_rows(check, candidates, 2).values()))
        self.assertFalse(output["accepted"])
        self.assertIsNone(output["replayed_score"])

    def test_selected_constructor_identity_and_metadata_parity_without_whole_freeze(self):
        key = (CASES[0], "cubic", 1); record = constructor_metadata(key)
        presented = dict(order=[0, 1], score=.25, score_kind="diagonal_reference_cubic_response",
                         status="complete", window_status="completed")
        A.selected_constructor_binding(record, presented, key, "a" * 64, "b" * 64)
        for field, value in (("schema", "wrong"), ("artifact_complete", False),
                ("source_freeze_sha256", "c" * 64), ("model_sha256", "d" * 64),
                ("case", CASES[1]), ("physical_error_read", True)):
            changed = copy.deepcopy(record); changed[field] = value
            with self.subTest(record=field), self.assertRaises(ValueError):
                A.selected_constructor_binding(changed, presented, key, "a" * 64, "b" * 64)
        for field, value in (("score_kind", "retained_first_order_response"), ("status", "wrong"),
                              ("window_status", "cpu_refusal"), ("order", [1, 0]), ("score", .5)):
            changed = dict(presented, **{field: value})
            with self.subTest(presented=field), self.assertRaises(ValueError):
                A.selected_constructor_binding(record, changed, key, "a" * 64, "b" * 64)

    def test_case_completion_cannot_override_incomplete_row_or_graph_failure(self):
        check = dict(status="complete", failures=[], checked_logical_rows=1)
        candidates = {("literal",): {}}
        rows = {("literal",): dict(status="complete")}
        self.assertTrue(A.independent_complete(check, candidates, rows))
        rows[("literal",)]["status"] = "incomplete"
        self.assertFalse(A.independent_complete(check, candidates, rows))
        self.assertFalse(A.independent_complete(dict(check, checked_logical_rows=0), {}, {}))
        case = dict(primary_panel_complete=True, independent_panel_complete=True,
                    objective_claims_complete=True, graph_valid=True, failures=[], input_hash="retained")
        A.case_failure(case, ValueError("saved state hash differs"))
        self.assertFalse(any(case[key] for key in ("primary_panel_complete",
            "independent_panel_complete", "objective_claims_complete", "graph_valid")))
        self.assertEqual(case["input_hash"], "retained")
        self.assertEqual(len(case["failures"]), 1)


if __name__ == "__main__":
    unittest.main()
