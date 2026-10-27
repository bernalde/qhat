"""Byte-identity delivery tests; archived code must remain inert evidence."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import pytest
from delivery import package


def test_delivery_keeps_exact_bytes_and_does_not_execute_archived_code(tmp_path):
    source = tmp_path/'producer'
    source.mkdir()
    payload = b'{"coefficient":"17/23"}\n'
    (source/'data.json').write_bytes(payload)
    (source/'never_execute.py').write_text('raise RuntimeError("archived source was executed")\n')
    files,archive,extracted = package(source,tmp_path/'delivery')
    assert set(files) == {'data.json','never_execute.py'}
    assert (extracted/'data.json').read_bytes() == payload
    assert (extracted/'never_execute.py').read_bytes() == (source/'never_execute.py').read_bytes()
    assert archive.is_file()


def test_delivery_refuses_nested_output_and_symlinks(tmp_path):
    source = tmp_path/'producer'
    source.mkdir()
    (source/'data').write_text('evidence')
    with pytest.raises(ValueError,match='separate'):
        package(source,source/'recursive')
    (source/'linked').symlink_to(source/'data')
    with pytest.raises(ValueError,match='symlink'):
        package(source,tmp_path/'delivery')
    assert not (tmp_path/'delivery').exists()


def test_delivery_refuses_unfinished_payload_before_packaging(tmp_path):
    source = tmp_path/'producer'
    source.mkdir()
    (source/'unfinished.json.partial').write_text('retained')
    with pytest.raises(ValueError,match='interrupted'):
        package(source,tmp_path/'delivery')
    assert (source/'unfinished.json.partial').read_text() == 'retained'
