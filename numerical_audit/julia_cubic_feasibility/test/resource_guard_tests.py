"""Pure filesystem/accounting tests and mocked supervisor tests; never run systemd."""
from contextlib import ExitStack, redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('ordering_resource_guard', Path(__file__).resolve().parents[1] / 'tools/resource_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class FakeChild:
    def __init__(self, done=False):
        self.done = done
        self.returncode = 0 if done else None

    def poll(self):
        return self.returncode if self.done else None

    def wait(self, timeout=None):
        self.done = True
        self.returncode = 0 if self.returncode is None else self.returncode
        return self.returncode

    def kill(self):
        self.done = True
        self.returncode = -9


class ResourceGuardTests(unittest.TestCase):
    def test_output_accounting_charges_files_and_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.assertEqual(guard.output_bytes(root / 'absent'), 0)
            output = root / 'one.toml'
            output.write_bytes(b'1234567')
            self.assertEqual(guard.output_bytes(output), 7)
            nested = root / 'subdir'
            nested.mkdir()
            (nested / 'two.toml').write_bytes(b'123')
            self.assertEqual(guard.output_bytes(root), 10)

    def test_missing_measurements_differ_from_observed_zero(self):
        self.assertEqual(guard.process_measurements('launch failed'), {})
        self.assertIsNone(guard.service_cpu_seconds('MemoryCurrent=0'))
        self.assertIsNone(guard.service_cpu_seconds('CPUUsageNSec=[not set]'))
        self.assertEqual(guard.service_cpu_seconds('CPUUsageNSec=0'), 0.0)

    def run_mocked(self, root, *, child=None, contents='', headroom=None, transition_race=False,
                   show_failure=False, launch_failure=False, cpu=10, wall=20):
        child = child or FakeChild()
        now = [0.0]
        calls = []

        def launch(invocation, *, stdout, stderr):
            calls.append(invocation)
            if launch_failure:
                raise OSError('mocked launch failure')
            stdout.write(contents)
            stdout.flush()
            return child

        def run(command, **kwargs):
            calls.append(command)
            self.assertIs(kwargs.get('check'), False)
            if 'show' in command:
                if show_failure:
                    raise OSError('mocked show failure')
                return subprocess.CompletedProcess(command, 0, stdout='CPUUsageNSec=0\nMemoryCurrent=0\n')
            if 'freeze' in command and transition_race:
                child.done = True
                child.returncode = 0
                return subprocess.CompletedProcess(command, 1)
            if 'stop' in command:
                child.done = True
                child.returncode = 0
            return subprocess.CompletedProcess(command, 0)

        def sleep(seconds):
            now[0] += seconds

        argv = ['--label', 'mock', '--logs', str(root / 'logs'), '--cpu-seconds', str(cpu),
                '--wall-seconds', str(wall), '--disk-mib', '1', '--output', str(root / 'output.toml'),
                '--', '/bin/true']
        with ExitStack() as stack:
            stack.enter_context(patch.object(guard, 'available', side_effect=headroom or [6 * 1024**3]))
            stack.enter_context(patch.object(guard.subprocess, 'Popen', side_effect=launch))
            stack.enter_context(patch.object(guard.subprocess, 'run', side_effect=run))
            stack.enter_context(patch.object(guard.time, 'monotonic', side_effect=lambda: now[0]))
            stack.enter_context(patch.object(guard.time, 'sleep', side_effect=sleep))
            stack.enter_context(redirect_stdout(io.StringIO()))
            with self.assertRaises(SystemExit) as result:
                guard.main(argv)
        receipt = json.loads((root / 'logs/mock.resources.json').read_text())
        self.assertEqual(result.exception.code, receipt['exit_code'])
        return receipt, calls

    def test_transition_child_exit_race_still_writes_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt, calls = self.run_mocked(Path(directory), transition_race=True,
                headroom=[6 * 1024**3, 4 * 1024**3])
            self.assertTrue(any('freeze' in command for command in calls))
            self.assertEqual(receipt['exit_code'], 0)
            self.assertIsNone(receipt['termination'])
            self.assertIsNone(receipt['guard_error'])

    def test_frozen_time_still_exhausts_absolute_wall_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt, calls = self.run_mocked(Path(directory), wall=2,
                headroom=[6 * 1024**3, 4 * 1024**3, 4 * 1024**3])
            self.assertEqual(receipt['termination'], 'Wall budget exhausted')
            self.assertEqual(receipt['settings']['RuntimeMaxSec'], '2')
            self.assertGreater(receipt['frozen_seconds'], 0)
            self.assertEqual(receipt['elapsed'], 3.0)
            self.assertNotEqual(receipt['exit_code'], 0)
            self.assertTrue(any('stop' in command for command in calls))
            self.assertEqual(receipt['settings']['MemoryMax'], '512M')
            self.assertEqual(receipt['settings']['MemorySwapMax'], '0')
            self.assertEqual(receipt['settings']['CPUQuota'], '25%')

    def test_polling_and_launch_errors_leave_best_effort_receipts(self):
        for failure in ('show_failure', 'launch_failure'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as directory:
                receipt, _ = self.run_mocked(Path(directory), headroom=[6 * 1024**3, 6 * 1024**3], **{failure: True})
                self.assertNotEqual(receipt['exit_code'], 0)
                self.assertIn('mocked', receipt['guard_error'])
                self.assertIsNone(receipt['final_accounting']['CPU_seconds'])
                self.assertFalse(receipt['final_accounting']['CPU_accounting_available'])
                self.assertFalse(receipt['final_accounting']['process_CPU_timing_available'])

    def test_short_process_final_measurements_are_charged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'output.toml').write_bytes(b'x' * 11)
            receipt, _ = self.run_mocked(root, child=FakeChild(done=True), cpu=4,
                contents='User time (seconds): 3.00\nSystem time (seconds): 2.00\nMaximum resident set size (kbytes): 17\n')
            self.assertEqual(receipt['samples'], [])
            self.assertEqual(receipt['final_accounting']['CPU_seconds'], 5.0)
            self.assertEqual(receipt['final_accounting']['output_bytes'], 11)
            self.assertTrue(receipt['final_accounting']['process_CPU_timing_available'])
            self.assertEqual(receipt['termination'], 'CPU budget exhausted')

    def test_existing_log_is_not_reused_as_fresh_measurement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            logs = root / 'logs'
            logs.mkdir()
            old = logs / 'mock.log'
            old.write_text('User time (seconds): 99\n')
            args = ['--label', 'mock', '--logs', str(logs), '--cpu-seconds', '10',
                    '--wall-seconds', '20', '--disk-mib', '1', '--output', str(root), '--', '/bin/true']
            with patch.object(guard, 'available', return_value=6 * 1024**3), patch.object(guard.subprocess, 'Popen') as launch:
                with self.assertRaises(FileExistsError):
                    guard.main(args)
                launch.assert_not_called()
            self.assertEqual(old.read_text(), 'User time (seconds): 99\n')
            self.assertFalse((logs / 'mock.resources.json').exists())


if __name__ == '__main__':
    unittest.main()
