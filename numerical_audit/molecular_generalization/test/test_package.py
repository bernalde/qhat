"""Small evidence round-trip and tamper regression; no scientific numerics."""
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import package as evidence


class EvidenceArchiveTests(unittest.TestCase):
    def test_round_trip_and_existing_inventory_prevents_relabelling(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); data = base/'original'; data.mkdir()
            (data/'INPUT_PROTOCOL_FROZEN.json').write_text('{}\n')
            (data/'state.bin').write_bytes(bytes(range(256)))
            archive = base/'study.tar.gz'
            evidence.pack(data, archive)
            evidence.extract(archive, base/'relocated')
            self.assertEqual((base/'relocated/data/state.bin').read_bytes(), bytes(range(256)))
            self.assertTrue((base/'relocated/EXTRACTION_VERIFIED.json').is_file())
            (data/'state.bin').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'inventory differs'):
                evidence.pack(data, base/'replacement.tar.gz')

    def test_archive_traversal_and_tree_symlink_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); archive = base/'malformed.tar.gz'
            with tarfile.open(archive, 'w:gz') as tar:
                info = tarfile.TarInfo('data/../escaped'); info.size = 1
                tar.addfile(info, io.BytesIO(b'x'))
            with self.assertRaisesRegex(ValueError, 'relative'):
                evidence.extract(archive, base/'out')
            self.assertFalse((base/'out/escaped').exists())
            data = base/'data'; data.mkdir()
            (data/'INPUT_PROTOCOL_FROZEN.json').write_text('{}\n')
            (data/'link').symlink_to(archive)
            with self.assertRaisesRegex(ValueError, 'symlinks'):
                evidence.pack(data, base/'bad.tar.gz')


if __name__ == '__main__':
    unittest.main()
