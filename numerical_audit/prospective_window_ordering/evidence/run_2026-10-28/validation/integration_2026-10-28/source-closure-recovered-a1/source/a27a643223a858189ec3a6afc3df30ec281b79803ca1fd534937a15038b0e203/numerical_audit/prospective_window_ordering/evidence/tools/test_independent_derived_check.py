"""UNRUN literal arithmetic/path/accounting contracts; no molecular successes.

These tests do not fabricate a complete scientific artifact graph. Root must
later exercise the checker on actual final artifacts and hash-updated tampering.
"""
import math
from pathlib import Path
import tempfile
import unittest

import independent_derived_check as D


class LiteralArithmetic(unittest.TestCase):
    def test_paired_subtraction_median_and_missing_repetition(self):
        differences = [D.subtract(.2, .5), D.subtract(.8, .6)]
        self.assertAlmostEqual(D.median(differences), -.05)
        self.assertIsNone(D.median([]))
        self.assertEqual(D.median([-.3]), -.3)
        self.assertEqual(D.outcome(-D.EPS), "tie")
        self.assertEqual(D.outcome(D.EPS), "tie")
        self.assertEqual(D.outcome(math.nextafter(D.EPS, math.inf)), "loss")

    def test_floor_and_objective_censor_suppress_ratio_without_hiding_loss(self):
        def operand(error, valid):
            return dict(values=dict(projective_error=error, infidelity=error), validated=valid)
        floor = D.paired_record(operand(0., True), operand(.1, True), None, None)
        self.assertEqual(floor["classification"], "win")
        self.assertIsNone(floor["baseline_over_window_error_ratio"])
        censored_objective = D.paired_record(operand(.5, False), operand(.1, True), 2., None)
        self.assertEqual(censored_objective["classification"], "loss")
        self.assertFalse(censored_objective["fully_validated"])
        self.assertIsNone(censored_objective["baseline_over_window_error_ratio"])
        missing = D.paired_record(None, operand(.1, True), None, 4.)
        self.assertEqual(missing["classification"], "censored")
        self.assertIsNone(missing["error_delta"])
        self.assertEqual(missing["wall_delta"], 4.)

    def test_all_missing_inventory_stays_censored(self):
        pairs, medians, summaries = D.derivations({}, {})
        self.assertEqual((len(pairs), len(medians), len(summaries)), (384, 192, 12))
        self.assertEqual({key[1] for key in summaries}, set(D.CONTROLS))
        self.assertTrue(all(x["classification"] == "censored" and x["CPU_delta"] is None for x in pairs.values()))
        self.assertTrue(all(x["median_paired_error_delta"] is None and not x["complete_repetitions"] for x in medians.values()))
        self.assertTrue(all(x["censored"] == 32 and x["complete_panel_median_of_paired_deltas"] is None for x in summaries.values()))

    def test_pareto_strictness_and_ambiguity(self):
        self.assertFalse(D.dominates(1., .2, 1., .2))
        self.assertTrue(D.dominates(.5, .2, 1., .2))
        self.assertFalse(D.dominates(2., .1, 1., .2))
        self.assertTrue(D.dominates(1., .1, 1., .2))
        self.assertFalse(D.dominates(1., .2 + 2*D.EPS, 2., .2))
        self.assertEqual(len(D.PATHS), 11)
        self.assertIn(("first", 1), D.PATHS)
        self.assertIn(("pauli", 2), D.PATHS)


class LiteralObjectiveComparison(unittest.TestCase):
    """Scalar compare_scores witnesses only, never molecular success graphs."""
    RTOL, ATOL = 1e-8, 1e-24

    def pass_record(self):
        return dict(status="verified", verified=True, within_tolerance=True,
                    numerical_floor=False, replayed_score=1., constructor_score=1.,
                    absolute_difference=0., relative_difference=0., tolerance=self.ATOL+self.RTOL)

    def accept(self, record):
        return D.objective_numeric_acceptance(record, self.RTOL, self.ATOL)

    def test_consistent_pass_floor_and_mismatch(self):
        self.assertTrue(self.accept(self.pass_record()))
        floor = dict(self.pass_record(), status="floor_censored", verified=False, numerical_floor=True,
                     replayed_score=0., constructor_score=0., tolerance=self.ATOL)
        self.assertFalse(self.accept(floor))
        # A floor is Boolean and does not erase the underlying numeric record.
        floor.update(replayed_score=self.ATOL, constructor_score=self.ATOL,
                     tolerance=self.ATOL+self.RTOL*self.ATOL)
        self.assertFalse(self.accept(floor))
        with self.assertRaises(ValueError):
            self.accept(dict(floor, tolerance=max(self.ATOL, self.RTOL*self.ATOL)))
        above = math.nextafter(self.ATOL, math.inf)
        floor.update(status="verified", verified=True, numerical_floor=False,
                     replayed_score=above, constructor_score=above,
                     tolerance=self.ATOL+self.RTOL*above)
        self.assertTrue(self.accept(floor))
        mismatch = dict(self.pass_record(), status="mismatch", verified=False, within_tolerance=False,
                        replayed_score=2., absolute_difference=1., relative_difference=.5,
                        tolerance=self.ATOL+2*self.RTOL)
        self.assertFalse(self.accept(mismatch))

    def test_missing_nonfinite_negative_and_boolean_numeric_fields(self):
        for field in ("replayed_score", "constructor_score", "absolute_difference", "relative_difference", "tolerance"):
            missing = self.pass_record()
            del missing[field]
            with self.subTest(field=field, value="missing"), self.assertRaises(ValueError):
                self.accept(missing)
            for bad in (math.nan, math.inf, -math.inf, -1., True):
                with self.subTest(field=field, value=bad), self.assertRaises(ValueError):
                    self.accept(dict(self.pass_record(), **{field: bad}))

    def test_flags_status_and_comparison_fields_are_reconstructed(self):
        changes = (dict(replayed_score=2.), dict(absolute_difference=1e-30),
                   dict(relative_difference=1e-30), dict(tolerance=1.),
                   dict(numerical_floor=True), dict(numerical_floor=0),
                   dict(verified=False), dict(verified=1), dict(within_tolerance=False),
                   dict(status="mismatch"))
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.accept(dict(self.pass_record(), **change))
        for field in ("numerical_floor", "verified", "within_tolerance"):
            missing = self.pass_record()
            del missing[field]
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.accept(missing)
        forged_floor = dict(self.pass_record(), replayed_score=0., constructor_score=0., tolerance=self.ATOL)
        with self.assertRaises(ValueError):
            self.accept(forged_floor)
        # Correct differences but false positive flags outside tolerance.
        forged_pass = dict(self.pass_record(), replayed_score=2., absolute_difference=1.,
                           relative_difference=.5, tolerance=self.ATOL+2*self.RTOL)
        with self.assertRaises(ValueError):
            self.accept(forged_pass)

    def test_refusal_and_missing_replay_remain_unaccepted(self):
        self.assertFalse(self.accept(dict(status="numerical_refusal", verified=False, constructor_score=1.)))
        self.assertFalse(self.accept(dict(status="budget_refusal", verified=False)))
        with self.assertRaises(ValueError):
            self.accept(dict(status="numerical_refusal", verified=True))
        for status in D.NUMERIC_OBJECTIVE_STATUSES:
            with self.subTest(status=status), self.assertRaises(ValueError):
                self.accept(dict(status=status, verified=False))

    def test_declared_tolerances_must_be_positive_and_bounded(self):
        for rtol, atol in ((0., self.ATOL), (self.RTOL, 0.), (2*self.RTOL, self.ATOL),
                           (self.RTOL, 2*self.ATOL), (math.nan, self.ATOL)):
            with self.subTest(rtol=rtol, atol=atol), self.assertRaises(ValueError):
                D.objective_numeric_acceptance(self.pass_record(), rtol, atol)


class LiteralIndependentReferences(unittest.TestCase):
    """Tiny scalar witness dictionaries; no states, cases, or artifact graph."""
    def operands(self, hard=False):
        primary = dict.fromkeys(D.METRICS, 0.)
        primary["projective_error"] = .25
        checked = dict(metrics=dict(primary), independent_reference_metrics=dict(primary),
            saved_reference_metric_discrepancies=dict.fromkeys(D.METRICS, 0.),
            independent_reference_metric_discrepancies=dict.fromkeys(D.METRICS, 0.),
            independent_unfused_replay=hard)
        if hard:
            checked.update(unfused_state_norm_residual=0., unfused_saved_reference_metrics=dict(primary),
                unfused_independent_reference_metrics=dict(primary),
                unfused_independent_reference_metric_discrepancies=dict.fromkeys(D.METRICS, 0.))
        return primary, checked

    def test_complete_scalar_witnesses_and_absolute_eps_boundary(self):
        for hard in (False, True):
            primary, checked = self.operands(hard)
            D.independent_metric_acceptance(primary, checked, hard)
        primary, checked = self.operands()
        checked["independent_reference_metrics"]["infidelity"] = D.EPS
        checked["independent_reference_metric_discrepancies"]["infidelity"] = D.EPS
        D.independent_metric_acceptance(primary, checked, False)
        checked["independent_reference_metrics"]["infidelity"] = math.nextafter(D.EPS, math.inf)
        checked["independent_reference_metric_discrepancies"]["infidelity"] = math.nextafter(D.EPS, math.inf)
        with self.assertRaises(ValueError):
            D.independent_metric_acceptance(primary, checked, False)

    def test_each_new_metric_surface_is_required_and_checked(self):
        for field in ("independent_reference_metrics", "unfused_saved_reference_metrics",
                      "unfused_independent_reference_metrics"):
            primary, checked = self.operands(True)
            del checked[field]
            with self.subTest(field=field, kind="missing"), self.assertRaises(KeyError):
                D.independent_metric_acceptance(primary, checked, True)
            for value in (2*D.EPS, math.nan, -1.):
                primary, checked = self.operands(True)
                checked[field]["infidelity"] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    D.independent_metric_acceptance(primary, checked, True)

    def test_stored_discrepancies_cannot_substitute_for_comparisons(self):
        for field in ("saved_reference_metric_discrepancies", "independent_reference_metric_discrepancies",
                      "unfused_independent_reference_metric_discrepancies"):
            primary, checked = self.operands(True)
            checked[field]["projective_error"] = 1e-30
            with self.subTest(field=field), self.assertRaises(ValueError):
                D.independent_metric_acceptance(primary, checked, True)
        primary, checked = self.operands(True)
        checked["independent_unfused_replay"] = False
        with self.assertRaises(ValueError):
            D.independent_metric_acceptance(primary, checked, True)


class LiteralAccounting(unittest.TestCase):
    def receipt(self, cpu):
        return dict(unit="literal-owned-job", elapsed=100., settings=dict(CPUQuota="25%", CPUWeight="5",
            IOWeight="10", Nice="19", MemoryHigh="448M", MemoryMax="512M", MemorySwapMax="0",
            CPUAccounting="yes", MemoryAccounting="yes"), final_accounting=dict(
            CPU_seconds=cpu, elapsed_seconds=100., CPU_accounting_available=cpu is not None,
            process_CPU_timing_available=cpu is not None), guard_error=None, cleanup_errors=[])

    def test_conservative_boundary_is_exposed_despite_raw_pass(self):
        plan = dict(attempt_CPU_seconds=64., attempt_wall_seconds=320.)
        final = dict(charged_CPU_seconds=64.002, absolute_guarded_wall_seconds=100.,
                     final_service_snapshot={"CPUUsageNSec": "64000500000"})
        result = D.final_charge(plan, self.receipt(63.999), final, "CPU time consumed: 1min 4.001s\n")
        self.assertFalse(result["raw_cap_exceeded"])
        self.assertTrue(result["conservative_cap_exceeded"])
        self.assertAlmostEqual(result["charge_minus_raw_CPU"], .003)

    def test_unmeasured_cpu_is_null_and_reserved_charge_is_not_measurement(self):
        result = D.final_charge(dict(attempt_CPU_seconds=64., attempt_wall_seconds=320.),
            self.receipt(None), dict(charged_CPU_seconds=64., absolute_guarded_wall_seconds=100.,
                                     final_service_snapshot={}), "")
        self.assertIsNone(result["raw_CPU_seconds"])
        self.assertFalse(result["raw_complete"])
        self.assertTrue(result["charge_is_reserved_fallback"])
        self.assertEqual(result["charged_CPU_seconds"], 64.)

    def test_forged_final_charge_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not reconcile"):
            D.final_charge(dict(attempt_CPU_seconds=64., attempt_wall_seconds=320.),
                self.receipt(2.), dict(charged_CPU_seconds=1., absolute_guarded_wall_seconds=100.,
                                       final_service_snapshot={}), "")

    def test_paths_cannot_leave_current_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = D.Store(Path(tmp))
            for path in ("../outside", "/absolute", "a\\b"):
                with self.subTest(path=path), self.assertRaises(ValueError):
                    store.path(path)


if __name__ == "__main__":
    unittest.main()
