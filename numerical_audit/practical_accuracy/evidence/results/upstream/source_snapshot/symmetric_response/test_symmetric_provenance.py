"""A current verifier may authenticate old source; propagation stays exact-source."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import pytest
import study


def fixture(tmp_path):
    archive=tmp_path/'source_snapshot';archive.mkdir()
    current=study.sources();historical={}
    for name in current:
        data=(study.HERE.parent/name).read_bytes()
        if name=='symmetric_response/PROTOCOL.md':data+=b'\nHistorical source fixture.\n'
        target=archive/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        historical[name]=hashlib.sha256(data).hexdigest()
    protocol=dict(source_commit='a'*40,environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')})
    frozen=dict(sources=historical,artifacts={},prior_models={})
    for name,value in [('PROTOCOL_FROZEN.json',protocol),('SCHEDULES_FROZEN.json',frozen)]:
        (tmp_path/name).write_text(json.dumps(value))
    (archive/'MANIFEST.json').write_text(json.dumps(dict(source_commit=protocol['source_commit'],sources=historical)))
    return archive,protocol


def test_propagation_rejects_changed_producing_sources_even_with_valid_archive(tmp_path):
    fixture(tmp_path)
    with pytest.raises(ValueError,match='source changed since freeze'):study.bind(tmp_path)


def test_independent_verifier_can_authenticate_complete_original_source(tmp_path):
    archive,protocol=fixture(tmp_path)
    assert study.bind(tmp_path,allow_archived_sources=True)==protocol


@pytest.mark.parametrize('corruption',['bytes','missing','extra','commit'])
def test_stale_or_incomplete_source_archive_is_rejected(tmp_path,corruption):
    archive,_=fixture(tmp_path);target=archive/'symmetric_response/PROTOCOL.md'
    if corruption=='bytes':target.write_text('changed')
    elif corruption=='missing':target.unlink()
    elif corruption=='extra':(archive/'unexpected.py').write_text('not imported')
    else:
        path=archive/'MANIFEST.json';record=json.loads(path.read_text());record['source_commit']='b'*40;path.write_text(json.dumps(record))
    with pytest.raises(ValueError):study.bind(tmp_path,allow_archived_sources=True)


def test_changed_source_domain_is_rejected_before_archive_paths_are_read(tmp_path):
    fixture(tmp_path);path=tmp_path/'SCHEDULES_FROZEN.json';record=json.loads(path.read_text())
    record['sources']['../outside']='0'*64;path.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='domain changed'):study.bind(tmp_path,allow_archived_sources=True)
