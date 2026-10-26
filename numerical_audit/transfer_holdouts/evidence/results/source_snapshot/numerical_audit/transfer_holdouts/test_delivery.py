"""Archive safety, deterministic bytes and exact binary inventory checks."""
import io
import json
import os
from pathlib import Path
import sys
import tarfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pytest
import delivery


@pytest.mark.parametrize('name,kind', [('../escape', tarfile.REGTYPE), ('delivery/link', tarfile.SYMTYPE)])
def test_archive_traversal_and_links_refused_before_extraction(tmp_path, name, kind):
    archive = tmp_path / 'bad.tar.gz'
    with tarfile.open(archive, 'w:gz') as stream:
        member = tarfile.TarInfo(name)
        member.type = kind
        if kind == tarfile.REGTYPE:
            member.size = 1
            stream.addfile(member, io.BytesIO(b'x'))
        else:
            member.linkname = '../../escape'
            stream.addfile(member)
    with pytest.raises(ValueError):
        delivery.extract_check(archive, tmp_path / 'extracted')
    assert not (tmp_path / 'escape').exists()


def test_unlisted_binary_archive_is_not_silently_omitted_from_delivery(tmp_path):
    (tmp_path / 'MANIFEST.json').write_text(json.dumps(dict(files={})))
    (tmp_path / 'required-state.npz').write_bytes(b'synthetic-unlisted-binary')
    with pytest.raises(ValueError, match='Strict complete delivery inventory differs'):
        delivery.check(tmp_path)


def test_fixed_payload_archive_bytes_ignore_mtime_and_output_filename(tmp_path, monkeypatch):
    # Isolate deterministic container encoding on an explicitly synthetic tree;
    # full scientific verification remains mandatory in the unpatched API.
    root = tmp_path / 'synthetic-delivery'
    root.mkdir()
    (root / 'data.json').write_text('{"fixture":true}\n')
    (root / 'state.npz').write_bytes(b'synthetic-binary')
    def fixture_check(path, commit=None):
        assert Path(path) == root
        return dict(status='synthetic-container-only')
    monkeypatch.setattr(delivery, 'check', fixture_check)
    first, second = tmp_path / 'first.tar.gz', tmp_path / 'second.tar.gz'
    one = delivery.archive(root, first)
    os.utime(root / 'data.json', (100, 100))
    os.utime(root / 'state.npz', (200, 200))
    two = delivery.archive(root, second)
    assert one == two and first.read_bytes() == second.read_bytes()
