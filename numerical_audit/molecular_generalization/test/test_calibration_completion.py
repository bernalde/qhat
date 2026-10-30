"""Recovery-accounting checks; synthetic receipts are not scientific evidence."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

path = Path(__file__).resolve().parents[1]/'calibration_completion.py'
spec = importlib.util.spec_from_file_location('molecular_completion_test', path)
c = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = c
spec.loader.exec_module(c)


class CompletionContracts(unittest.TestCase):
    def test_service_footer_rounding_and_failed_costs_are_retained(self):
        self.assertEqual(c.duration('3min 58.071s'), (238.071, .0005))
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); log = base/'attempt.log'
            log.write_text('CPU time consumed: 3min 58.071s\n')
            receipt = dict(command=['owned-fixture'], completed_utc='fixed-test-time',
                settings=c.PROFILE, environment=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
                exit_code=2, elapsed=952.4, guard_error=None, cleanup_errors=[], termination=None, log=str(log),
                process_measurements=dict(user_cpu_seconds=220.25, system_cpu_seconds=17.80),
                final_accounting=dict(process_CPU_timing_available=True, CPU_seconds=238.05, output_bytes=1024))
            target = base/'receipt.json'; target.write_text(json.dumps(receipt))
            _, cost = c.receipt_cost(target, ['owned-fixture'], allow_failure=True)
            self.assertEqual(cost['charged_CPU_seconds'], 238.072)
            self.assertTrue(cost['failed_attempt'])
            with self.assertRaisesRegex(c.Refusal, 'Successful completion'):
                c.receipt_cost(target, ['owned-fixture'], allow_failure=False)
            with self.assertRaisesRegex(c.Refusal, 'command receipt'):
                c.receipt_cost(target, ['other-worker'], allow_failure=True)
            receipt['process_measurements'] = {}
            receipt['final_accounting']['process_CPU_timing_available'] = False
            log.write_text('Only live samples, not terminal proof\n')
            target.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(c.Refusal, 'Unknown terminal CPU'):
                c.receipt_cost(target, ['owned-fixture'], allow_failure=True)

    def test_exact_cache_pins_refuse_changed_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'row.json'; path.write_text('{"score":1}\n')
            reference = c.ref(path)
            self.assertEqual(c.pin(reference), path)
            path.write_text('{"score":2}\n')
            with self.assertRaisesRegex(c.Refusal, 'Pinned artifact changed'):
                c.pin(reference)

    def test_endpoint_fusion_normalization_and_underflow(self):
        model = dict(blocks=[{1: 1, 2: 1}, {1: 1, 2: 1, 3: 1}, {1: 1}])
        order = [0, 1, 2]; steps = 19
        valid = dict(value=1/steps**4, cubic_response=1., steps=steps)
        c.score_valid(valid, order, model)
        with self.assertRaisesRegex(c.Refusal, 'endpoint fusion'):
            c.score_valid(dict(valid, steps=20), order, model)
        with self.assertRaisesRegex(c.Refusal, 'underflow'):
            c.score_valid(dict(valid, value=0.), order, model)


if __name__ == '__main__':
    unittest.main()
