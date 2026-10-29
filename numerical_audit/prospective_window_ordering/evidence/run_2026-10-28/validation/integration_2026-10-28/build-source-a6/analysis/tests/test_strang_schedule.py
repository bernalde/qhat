"""Budget and physical-evolution regressions for the maintained Strang planner."""
import contextlib
import io
import itertools
import json
import unittest

import numpy as np
from scipy.linalg import expm

from qhat.analysis.strang_schedule import (
    iter_strang_pulses, main, plan_strang_schedule,
)


def unmerged_units(lengths, steps):
    """Literal half-step palindrome, without using the production cost formula."""
    for _ in range(steps):
        yield from range(len(lengths))
        yield from reversed(range(len(lengths)))


def oracle_rotation_count(lengths, steps):
    # Merge identical consecutive group exponentials, then expand each group.
    cost, previous = 0, None
    for group in unmerged_units(lengths, steps):
        if group != previous:
            cost += lengths[group]
        previous = group
    return cost


def pulse_matrix(terms, plan, time):
    actual = np.eye(terms[0].shape[0], dtype=complex)
    for term, numerator in plan.iter_pulses():
        actual = expm(-1j * time * numerator / plan.duration_denominator * terms[term]) @ actual
    return actual


def paulis():
    identity = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.diag([1., -1.])
    return identity, x, y, z


class StrangScheduleTests(unittest.TestCase):
    def test_singleton_budget_boundary_and_surplus(self):
        for size in (1, 3, 19):
            for budget in (0, size - 1, size, size + 1000):
                with self.subTest(size=size, budget=budget):
                    plan = plan_strang_schedule([size], budget)
                    self.assertEqual(plan.feasible, budget >= size)
                    self.assertEqual(plan.steps, int(budget >= size))
                    self.assertEqual(plan.rotations, size if budget >= size else 0)
                    self.assertEqual(list(plan.iter_pulses()), [(i, 2) for i in range(size)] if budget >= size else [])
                    self.assertEqual(plan.duration_denominator, 2 if budget >= size else None)

    def test_singleton_direct_pulses_have_full_not_double_time(self):
        self.assertEqual(list(iter_strang_pulses([3], 7)), [(0, 14), (1, 14), (2, 14)])

    def test_small_partitions_match_independent_palindrome_counts(self):
        for groups in range(1, 5):
            for lengths in itertools.product(range(1, 4), repeat=groups):
                for budget in range(61):
                    plan = plan_strang_schedule(lengths, budget)
                    pulse_list = list(plan.iter_pulses())
                    self.assertLessEqual(plan.rotations, budget)
                    self.assertEqual(len(pulse_list), plan.rotations)
                    if not plan.feasible:
                        self.assertLess(budget, oracle_rotation_count(lengths, 1))
                        self.assertIsNone(plan.duration_denominator)
                        continue
                    self.assertEqual(plan.rotations, oracle_rotation_count(lengths, plan.steps))
                    if groups > 1:
                        self.assertGreater(oracle_rotation_count(lengths, plan.steps + 1), budget)
                    totals = [0] * sum(lengths)
                    for index, numerator in pulse_list:
                        totals[index] += numerator
                    self.assertEqual(totals, [plan.duration_denominator] * sum(lengths))

    def test_worked_orientation_example(self):
        forward = plan_strang_schedule([1, 2, 3], 60)
        reverse = plan_strang_schedule([3, 2, 1], 60)
        self.assertEqual((forward.steps, forward.rotations), (7, 57))
        self.assertEqual((reverse.steps, reverse.rotations), (7, 59))
        self.assertEqual(plan_strang_schedule([1, 2, 3], 8).steps, 0)
        self.assertEqual(plan_strang_schedule([1, 2, 3], 9).steps, 1)

    def test_singleton_equals_exact_dense_evolution(self):
        identity, x, y, z = paulis()
        # Non-diagonal commuting group: catches a factor-of-two duration error.
        terms = [0.3 * np.kron(x, x), -0.2 * np.kron(y, y), 0.7 * np.kron(z, z)]
        for selected in (terms[:1], terms):
            plan = plan_strang_schedule([len(selected)], 100)
            for time in (0., 0.43, -0.7):
                actual = pulse_matrix(selected, plan, time)
                expected = expm(-1j * time * sum(selected))
                np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0.)
        actual = pulse_matrix(terms, plan_strang_schedule([3], 3), 0.43)
        self.assertGreater(np.linalg.norm(actual - expm(-2j * 0.43 * sum(terms))), 0.1)

    def test_multiunit_fusion_equals_unmerged_group_product(self):
        identity, x, y, z = paulis()
        groups = [[0.3 * np.kron(x, x), -0.2 * np.kron(y, y)],
                  [0.4 * np.kron(z, identity)],
                  [0.6 * np.kron(identity, z), 0.7 * np.kron(z, z)]]
        outputs = []
        for chosen in (groups, list(reversed(groups))):
            lengths = [len(g) for g in chosen]
            terms = [t for group in chosen for t in group]
            plan = plan_strang_schedule(lengths, 16)
            self.assertEqual(plan.steps, 2)
            for time in (0., 0.43, -0.7):
                actual = pulse_matrix(terms, plan, time)
                expected = np.eye(4, dtype=complex)
                for group in unmerged_units(lengths, plan.steps):
                    expected = expm(-1j * time / (2 * plan.steps) * sum(chosen[group])) @ expected
                np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0.)
            outputs.append(actual)
        self.assertGreater(np.linalg.norm(outputs[0] - outputs[1]), 1e-4)

    def test_invalid_partition_and_budget_are_rejected(self):
        for lengths in ([], [0], [-1], [1, 0], [1.5], [True]):
            with self.subTest(lengths=lengths), self.assertRaises((ValueError, TypeError)):
                plan_strang_schedule(lengths, 10)
        for budget in (-1, 1.5, True):
            with self.subTest(budget=budget), self.assertRaises((ValueError, TypeError)):
                plan_strang_schedule([1], budget)
        for steps in (0, -1, 1.5, True):
            with self.subTest(steps=steps), self.assertRaises((ValueError, TypeError)):
                iter_strang_pulses([1], steps)

    def test_integer_like_values_and_one_pass_iterables(self):
        plan = plan_strang_schedule(iter([np.int64(1), np.int64(2)]), np.int64(10))
        self.assertEqual((plan.steps, plan.rotations), (3, 10))

    def test_cli_exposes_singleton_and_infeasible_plans(self):
        for budget in (2, 3):
            capture = io.StringIO()
            with contextlib.redirect_stdout(capture):
                self.assertEqual(main(['--unit-lengths', '3', '--budget', str(budget)]), 0)
            record = json.loads(capture.getvalue())
            self.assertEqual(record['feasible'], budget == 3)
            self.assertEqual(record['pulses'], [[0, 2], [1, 2], [2, 2]] if budget == 3 else [])


if __name__ == '__main__':
    unittest.main()
