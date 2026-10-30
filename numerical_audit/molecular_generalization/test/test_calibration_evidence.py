"""Contained calibration copies, using metadata fixtures and no numerical work."""
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import common as k


class CalibrationEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base/'repo'
        self.here = self.root/'numerical_audit/molecular_generalization'
        self.here.mkdir(parents=True)
        self.prefix = 'evidence/calibration_fixture'
        self.report = self.raw(dict(complete=True, physical_errors_read=False,
                                    original_report='/historical/work/CALIBRATION.json'))
        self.expected = hashlib.sha256(self.report).hexdigest()
        self.members = {'resume/work/CALIBRATION.json': self.report,
                        'producer/source.py': b'# original producer\n', 'original/empty.log': b''}
        self.manifest = dict(schema='molecular-calibration-evidence-v1', prefix=self.prefix,
            accepted_report='resume/work/CALIBRATION.json', accepted_sha256=self.expected,
            scope='Synthetic metadata only', files={name: dict(bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest(), original_locator='/historical/'+name)
                for name, raw in self.members.items()})
        for name, raw in self.members.items():
            self.put(self.prefix+'/'+name, raw)
        self.put('CALIBRATION.json', self.report)
        self.save_manifest()
        self.put(self.prefix+'/undeclared/private.json', b'{"excluded":true}\n')
        self.config = dict(main_execution_allowed=False, protocol_state='development',
                           annealing=dict(calibration_sha256=self.expected))
        self.put('PROTOCOL.json', self.raw(self.config))
        self.put('ENVIRONMENT.json', b'{}\n')

    @staticmethod
    def raw(value):
        return (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()

    def put(self, name, raw):
        path = self.here/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
        return path

    def save_manifest(self):
        self.put('CALIBRATION_EVIDENCE.json', self.raw(self.manifest))

    def sources(self):
        return {str((self.here/name).relative_to(self.root)): k.sha(self.here/name)
                for name in ('CALIBRATION.json', 'CALIBRATION_EVIDENCE.json',
                             'PROTOCOL.json', 'ENVIRONMENT.json')}

    def context(self):
        stack = ExitStack()
        for name, value in dict(HERE=self.here, ROOT=self.root,
                configuration=lambda: self.config, registry=lambda: dict(cases=[]),
                source_map=self.sources, committed_sources=lambda: 'fixture-commit').items():
            stack.enter_context(patch.object(k, name, value))
        stack.enter_context(patch.object(k.prior, 'environment', return_value={}))
        return stack

    def test_freeze_copies_only_declared_bytes_and_replays_without_original_tree(self):
        data = self.base/'data'
        with self.context():
            k.freeze(data, mode='preflight')
            snapshot = data/'source_snapshot'/self.here.relative_to(self.root)
            for name, raw in self.members.items():
                self.assertEqual((snapshot/self.prefix/name).read_bytes(), raw)
            self.assertEqual((snapshot/'CALIBRATION_EVIDENCE.json').read_bytes(),
                             (self.here/'CALIBRATION_EVIDENCE.json').read_bytes())
            self.assertEqual((snapshot/'CALIBRATION.json').read_bytes(), self.report)
            self.assertFalse((snapshot/self.prefix/'undeclared').exists())
            shutil.rmtree(self.here)
            frozen = k.frozen(data)
            self.assertEqual(frozen['calibration_evidence']['declared_files'], 3)
            self.assertEqual(frozen['calibration_evidence']['accepted_sha256'], self.expected)
            copied = k.read(snapshot/'CALIBRATION_EVIDENCE.json')
            self.assertEqual(copied['files']['producer/source.py']['original_locator'],
                             '/historical/producer/source.py')

    def test_changed_member_size_or_hash_is_refused_before_output_creation(self):
        member = self.here/self.prefix/'producer/source.py'
        for raw in (b'changed size', b'# replaced producer\n'):
            with self.subTest(raw=raw), self.context():
                member.write_bytes(raw)
                output = self.base/'refused'
                with self.assertRaisesRegex(ValueError, 'Calibration evidence (size|hash/size) differs'):
                    k.freeze(output, mode='preflight')
                self.assertFalse(output.exists())

    def test_top_report_must_be_the_exact_accepted_report(self):
        self.put('CALIBRATION.json', self.report+b' ')
        with self.assertRaisesRegex(ValueError, 'Top calibration must equal'):
            k.calibration_evidence(self.here, self.expected)
        with self.assertRaisesRegex(ValueError, 'Accepted calibration report declaration'):
            k.calibration_evidence(self.here, '0'*64)

    def test_self_consistent_incomplete_or_physical_report_is_still_refused(self):
        for report in (dict(complete=False, physical_errors_read=False),
                       dict(complete=True, physical_errors_read=True)):
            with self.subTest(report=report):
                raw = self.raw(report)
                digest = hashlib.sha256(raw).hexdigest()
                name = self.manifest['accepted_report']
                self.put('CALIBRATION.json', raw)
                self.put(self.prefix+'/'+name, raw)
                self.manifest['accepted_sha256'] = digest
                self.manifest['files'][name].update(bytes=len(raw), sha256=digest)
                self.save_manifest()
                with self.assertRaisesRegex(ValueError, 'Complete old-only calibration'):
                    k.calibration_evidence(self.here, digest)

    def test_paths_reject_traversal_absolute_empty_and_backslash_components(self):
        entry = self.manifest['files']['producer/source.py']
        for name in ('../escape', '/absolute', 'a//b', 'a/./b', 'a\\b', 'a\x00b'):
            with self.subTest(name=name):
                self.manifest['files'][name] = entry
                self.save_manifest()
                with self.assertRaisesRegex(ValueError, 'Canonical calibration evidence path'):
                    k.calibration_evidence(self.here, self.expected)
                del self.manifest['files'][name]
        self.manifest['prefix'] = '../elsewhere'
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, 'Canonical calibration evidence path'):
            k.calibration_evidence(self.here, self.expected)

    def test_symlink_file_and_directory_are_refused_even_with_matching_bytes(self):
        member = self.here/self.prefix/'producer/source.py'
        target = self.put('same.py', self.members['producer/source.py'])
        member.unlink()
        member.symlink_to(target)
        with self.assertRaisesRegex(ValueError, 'Calibration evidence I/O refused'):
            k.calibration_evidence(self.here, self.expected)
        member.unlink()
        member.write_bytes(self.members['producer/source.py'])
        directory = member.parent
        directory.rename(directory.with_name('real_producer'))
        directory.symlink_to(directory.with_name('real_producer'), target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'Calibration evidence I/O refused'):
            k.calibration_evidence(self.here, self.expected)

    def test_nonregular_file_and_missing_member_are_refused(self):
        member = self.here/self.prefix/'original/empty.log'
        member.unlink()
        with self.assertRaisesRegex(ValueError, 'Calibration evidence I/O refused'):
            k.calibration_evidence(self.here, self.expected)
        os.mkfifo(member)
        with self.assertRaisesRegex(ValueError, 'Bounded regular calibration evidence'):
            k.calibration_evidence(self.here, self.expected)

    def test_duplicate_json_and_declared_resource_overflow_are_refused(self):
        self.put('CALIBRATION_EVIDENCE.json', b'{"schema":1,"schema":2}')
        with self.assertRaisesRegex(ValueError, 'Duplicate calibration JSON key'):
            k.calibration_evidence(self.here, self.expected)
        self.save_manifest()
        with patch.object(k, 'MAX_CALIBRATION_FILES', 1):
            with self.assertRaisesRegex(ValueError, 'Bounded calibration manifest schema'):
                k.calibration_evidence(self.here, self.expected)
        with patch.object(k, 'MAX_CALIBRATION_TOTAL_BYTES', 1):
            with self.assertRaisesRegex(ValueError, 'total byte cap'):
                k.calibration_evidence(self.here, self.expected)

    def test_copy_refuses_overwrite_and_changed_snapshot_is_detected(self):
        copied = self.base/'copied'; copied.mkdir()
        k.calibration_evidence(self.here, self.expected, copy_to=copied)
        with self.assertRaisesRegex(ValueError, 'Calibration evidence I/O refused'):
            k.calibration_evidence(self.here, self.expected, copy_to=copied)
        self.assertEqual((copied/'CALIBRATION_EVIDENCE.json').read_bytes(),
                         (self.here/'CALIBRATION_EVIDENCE.json').read_bytes())
        data = self.base/'data'
        with self.context():
            k.freeze(data, mode='preflight')
            snapshot = data/'source_snapshot'/self.here.relative_to(self.root)
            (snapshot/self.prefix/'producer/source.py').write_bytes(b'# replaced producer\n')
            with self.assertRaisesRegex(ValueError, 'Calibration evidence hash/size differs'):
                k.frozen(data)

    def test_copy_failure_never_publishes_freeze_marker(self):
        original = k._calibration_file
        def changed_during_copy(root, relative, **kwargs):
            if kwargs.get('copy_to') is not None and relative.endswith('/producer/source.py'):
                (self.here/relative).write_bytes(b'changed after initial validation')
            return original(root, relative, **kwargs)
        data = self.base/'interrupted'
        with self.context(), patch.object(k, '_calibration_file', changed_during_copy):
            with self.assertRaisesRegex(ValueError, 'Calibration evidence size differs'):
                k.freeze(data, mode='preflight')
        self.assertFalse((data/'INPUT_PROTOCOL_FROZEN.json').exists())

    def test_new_freeze_cannot_drop_its_evidence_binding(self):
        data = self.base/'data'
        with self.context():
            k.freeze(data, mode='preflight')
            path = data/'INPUT_PROTOCOL_FROZEN.json'
            value = k.read(path)
            del value['calibration_evidence']
            path.write_bytes(self.raw(value))
            with self.assertRaisesRegex(ValueError, 'Frozen calibration evidence binding changed'):
                k.frozen(data)

    def test_source_scan_excludes_evidence_and_test_subtrees(self):
        self.put('helper.py', b'# current source\n')
        self.put('test/fixture.py', b'# test\n')
        with patch.object(k, 'HERE', self.here), patch.object(k, 'ROOT', self.root), \
                patch.object(k.retained, 'source_paths', return_value=[]):
            paths = k.source_paths()
        self.assertIn(self.here/'helper.py', paths)
        self.assertIn(self.here/'CALIBRATION_EVIDENCE.json', paths)
        self.assertFalse(any('evidence' in p.relative_to(self.here).parts or
                             'test' in p.relative_to(self.here).parts
                             for p in paths if p.is_relative_to(self.here)))


if __name__ == '__main__':
    unittest.main()
