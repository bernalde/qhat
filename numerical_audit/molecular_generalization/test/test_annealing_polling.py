"""Controlled-clock monitoring regressions; no real Julia child is launched."""
from pathlib import Path
import signal
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import call, patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_annealing_adapter as fixtures

a = fixtures.a


class ControlledClock:
    def __init__(self):
        self.cpu = 0.
        self.wall = 0.

    def process_time(self):
        return self.cpu

    def monotonic(self):
        return self.wall

    def sleep(self, seconds):
        self.wall += seconds


class PollingTests(unittest.TestCase):
    def test_wall_advances_during_poll_sleeps_then_only_owned_child_is_killed_and_reaped(self):
        clock = ControlledClock()
        pid = 876543
        process = SimpleNamespace(pid=pid, returncode=None)
        usage = SimpleNamespace(ru_utime=.125, ru_stime=.0625, ru_maxrss=73)
        stats, checked = {}, []

        def wait_owned(requested_pid, flags):
            self.assertEqual(requested_pid, pid)
            self.assertIn(flags, (a.os.WNOHANG, 0))
            return (0, 0, None) if flags == a.os.WNOHANG else (pid, signal.SIGKILL, usage)

        with patch.object(a.subprocess, 'Popen', return_value=process), \
             patch.object(a.os, 'wait4', side_effect=wait_owned) as wait, \
             patch.object(a.os, 'killpg') as kill, \
             patch.object(a, '_proc_cpu', return_value=.125), \
             patch.object(a._time, 'monotonic', side_effect=clock.monotonic), \
             patch.object(a._time, 'sleep', side_effect=clock.sleep) as sleep:
            outcome = a._run_worker('/fixture/julia', Path('/fixture/request'), Path('/fixture/reply'),
                1., .03, checked.append, stats)

        self.assertEqual(sleep.call_args_list, [call(a.OWNED_CHILD_POLL_SECONDS)]*2)
        self.assertEqual(wait.call_args_list, [call(pid, a.os.WNOHANG)]*3+[call(pid, 0)])
        kill.assert_called_once_with(pid, signal.SIGKILL)
        self.assertEqual(checked, [.125, .125])
        self.assertEqual(outcome, 'child_wall_limit')
        self.assertTrue(stats['child_reaped'])
        self.assertEqual(stats['child_returncode'], -signal.SIGKILL)
        self.assertEqual(stats['child_process_CPU_seconds'], usage.ru_utime+usage.ru_stime)
        self.assertEqual(stats['child_peak_RSS_KiB'], usage.ru_maxrss)
        self.assertEqual(stats['child_wall_seconds'], 2*a.OWNED_CHILD_POLL_SECONDS)

    def test_valid_checked_checkpoint_still_fails_final_complete_cpu_overrun_gate(self):
        # Reuse the existing tiny native/checkpoint fixture without inheriting
        # its test methods or launching a child process.
        fixture = fixtures.AdapterTests(methodName='runTest')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        clock = ControlledClock()
        budget_check = fixtures.native._EXTERNAL_BUDGET_CHECK
        previous = budget_check.get()

        class CompleteCallBoundary:
            resets = 0

            def get(self):
                return budget_check.get()

            def set(self, value):
                return budget_check.set(value)

            def reset(self, token):
                budget_check.reset(token)
                self.resets += 1
                # All checkpoint/native checks completed under the allowance.
                # Account an injected finalization CPU overrun at the boundary.
                clock.cpu = 10.25

        boundary = CompleteCallBoundary()
        with patch.object(fixtures.native, '_EXTERNAL_BUDGET_CHECK', boundary), \
             patch.object(a._time, 'process_time', side_effect=clock.process_time), \
             patch.object(a._time, 'monotonic', side_effect=clock.monotonic), \
             fixture.environment(), patch.object(a.shutil, 'which', return_value='/fixture/julia'):
            result = fixture.call(cpu_seconds=10., wall_seconds=50.)

        self.assertEqual(fixture.worker_calls, 1)
        self.assertEqual(boundary.resets, 1)
        self.assertIs(budget_check.get(), previous)
        self.assertTrue(result.stats['four_chain_coverage'])
        self.assertEqual(result.stats['native_verifications'], 1)
        self.assertIsNotNone(result.raw_julia_candidate['native_score'])
        self.assertFalse(result.eligible)
        self.assertEqual(result.status, 'whole_call_overrun')
        self.assertEqual(result.refusal['code'], 'whole_call_overrun')
        self.assertEqual(result.stats['parent_process_CPU_seconds'], 10.25)
        self.assertEqual(result.CPU_seconds, 10.25+result.stats['child_process_CPU_seconds'])
        self.assertGreater(result.CPU_seconds, 10.)


if __name__ == '__main__':
    unittest.main()
