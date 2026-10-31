"""Cost and identity regressions, with no molecular panel or timing thresholds."""
from pathlib import Path
import importlib.util
import math
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import common as k
import construct


class ContractTests(unittest.TestCase):
    def test_unmeasured_cpu_is_not_a_free_attempt(self):
        receipt = dict(elapsed=2., samples=[], process_measurements={})
        value = construct.cost_from_receipt(receipt)
        self.assertIsNone(value['raw_CPU_seconds'])
        self.assertIsNone(value['charged_CPU_seconds'])
        receipt['process_measurements'] = dict(user_cpu_seconds=0., system_cpu_seconds=0.)
        self.assertEqual(construct.cost_from_receipt(receipt)['charged_CPU_seconds'], .011)
        receipt['process_measurements'] = dict(user_cpu_seconds=1., system_cpu_seconds=.2)
        receipt['samples'] = [dict(service='CPUUsageNSec=1250000000')]
        self.assertEqual(construct.cost_from_receipt(receipt)['charged_CPU_seconds'], 1.25)

    def test_shared_parent_and_cheap_allowances(self):
        c = k.configuration()
        p = dict(costs=dict(total=dict(charged_CPU_seconds=7., wall_seconds=29.)))
        self.assertEqual(construct.allowances(c, 'full', 32, p), (25., 131))
        self.assertEqual(construct.allowances(c, 'annealing', 64, p), (57., 291))
        self.assertEqual(construct.allowances(c, 'magnitude', 0), (8, 40))
        self.assertEqual(construct.allowances(c, 'parent', 0), (8, 40))
        p['costs']['total']['charged_CPU_seconds'] = None
        with self.assertRaisesRegex(ValueError, 'Measured parent'):
            construct.allowances(c, 'window', 32, p)

    def test_same_four_starts_and_complete_nine_identity_roster(self):
        base = list(range(20))
        value = construct.restarts(base, 'old-case', 1)
        self.assertEqual(value, construct.restarts(base, 'old-case', 1))
        self.assertEqual(len({tuple(o) for o in value['orders']}), 4)
        for order in value['orders']:
            self.assertEqual(sorted(order), base)
            self.assertEqual((order[0], order[-1]), (base[0], base[-1]))
        self.assertEqual(len(set(k.identities(1))), 9)
        self.assertEqual(len(set(k.identities(2))), 18)

    def test_exact_target_partition_and_dynamic_reference_dimension(self):
        model = dict(n=12, electrons=8, denominator='3',
            coefficients=[[[[0, 'Z']], '1'], [[[1, 'Z']], '2']],
            blocks=[[[[[0, 'Z']], '1']], [[[[1, 'Z']], '2']]])
        blocks = k.blocks(model)
        self.assertEqual(sum(len(b) for b in blocks), 2)
        self.assertEqual(str(blocks[0][((0, 'Z'),)]), '1/3')
        model['blocks'].append(model['blocks'][0])
        with self.assertRaisesRegex(ValueError, 'partition'):
            k.blocks(model)
        model['n'], model['electrons'] = 14, 7
        with self.assertRaisesRegex(ValueError, 'fixed-N'):
            k.blocks(model)

    def test_main_is_refused_before_output_creation(self):
        disabled = dict(k.configuration(), main_execution_allowed=False,
                        protocol_state='old-case_preflight_development')
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'main'
            with patch.object(k, 'configuration', return_value=disabled):
                with self.assertRaisesRegex(ValueError, 'Main disabled'):
                    k.freeze(output, mode='main')
            self.assertFalse(output.exists())

    def test_contained_artifact_rejects_traversal_and_symlink_escape(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'data'; root.mkdir()
            (root/'escape').symlink_to(Path(temp), target_is_directory=True)
            for value in ('../other', '/absolute', 'a//b', 'a/./b', 'escape/file'):
                with self.assertRaises(ValueError):
                    k.contained(root, value)

    def test_observer_preserves_values_and_restores_on_failure(self):
        native = k.bootstrap()
        from qhat.analysis import window_ordering as window
        class Context:
            def __init__(self): self.calls = 0
            def score(self, order):
                self.calls += 1
                return native.Score(1., 16., 2)
        context = Context()
        original_score = context.score
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'reply'; path.write_text('exact raw checkpoint')
            reply = dict(score=1., cubic_response=16., steps=2)
            seen = []
            def fake_reply(*args): return (1, 0), reply
            with patch.object(window, '_reply', fake_reply):
                with self.assertRaisesRegex(RuntimeError, 'after check'):
                    with construct.observe_window(seen):
                        candidate, got = window._reply(path, 'digest', context)
                        self.assertIs(got, reply)
                        self.assertEqual(context.score(candidate).value, 1.)
                        raise RuntimeError('after check')
                self.assertIs(window._reply, fake_reply)
            self.assertEqual(context.calls, 1)
            self.assertEqual(context.score, original_score)
            self.assertEqual(seen[0]['native_score'], dict(value=1., cubic_response=16., steps=2))
            self.assertNotIn('native_refusal', seen[0])


if __name__ == '__main__':
    unittest.main()
