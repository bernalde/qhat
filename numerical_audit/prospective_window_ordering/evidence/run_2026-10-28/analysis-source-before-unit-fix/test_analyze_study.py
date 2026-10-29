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
import statistics
import tempfile
import unittest

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
        self.assertEqual(len(pairs), 8 * 32)
        self.assertEqual(len(medians), 8 * 16)
        self.assertEqual(len(summaries), 8)
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
    def test_analyze_retains_352_cells_and_partial_cost_lower_bounds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = partial_graph(root)
            r1 = bound_json(root, "receipts/first.json", guard_receipt("fixture-call-attempt1", 1.25, 5.))
            r2 = bound_json(root, "receipts/retry.json", guard_receipt("fixture-call-attempt2", .75, 3.))
            prep = bound_json(root, "receipts/prep.json", guard_receipt("fixture-preparation", .5, 2.))
            index["partial_attempts"] = [dict(case=CASES[0], arm="window_sparse", repetition=1,
                receipt=r, record=None) for r in (r1, r2)]
            index["phase_receipts"] = [dict(phase="construction", receipt=r1), dict(phase="preparation", receipt=prep)]
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
            self.assertEqual(len(summary["paired_observations"]), 256)
            self.assertTrue(all(p["classification"] == "censored" for p in summary["paired_observations"]))
            self.assertEqual(len(summary["comparisons"]), 8)
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
            self.assertEqual(summary["other_phase_totals"]["construction"]["indexed_attempts"], 0)
            self.assertEqual(summary["other_phase_totals"]["construction"]["CPU_seconds"], 0)
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
            self.assertEqual(summary["construction_totals"]["known_CPU_seconds"], 0)
            self.assertIsNone(summary["construction_totals"]["CPU_seconds"])


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


if __name__ == "__main__":
    unittest.main()
