"""Execution gates, archive allocation guards and complete caller charges."""
from pathlib import Path
import io
import hashlib
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import zipfile
import numpy as np
import pytest
import transfer_common as k

study = k.load('transfer_full_panel_study_tests', k.HERE / 'study.py')


def reviewed_record():
    names = ['numerical_audit/practical_accuracy/' + name for name in
             ('common.py', 'construction.py', 'study.py', 'verify.py', 'PROTOCOL.md')]
    names.append('numerical_audit/compiled_accuracy/compiler.py')
    blob = b'authenticated interface fixture\n'
    return dict(schema=1, repository='bernalde/qhat', issue=16, state='reviewed_frozen_interface',
                automated_handoff='READY FOR HUMAN REVIEW',
                requested_reviewers=['gpt-6-astra/high', 'gemini-4-argon-high/high'],
                actionable_findings_remaining=0, controller_policy_version=2, head='a' * 40,
                cost_convention='caller_through_return_including_serialization_hash_cleanup',
                source_sha256={name: hashlib.sha256(blob).hexdigest() for name in names}), blob


def test_unreviewed_dependency_stops_before_reference_io_or_training(tmp_path, monkeypatch):
    path = tmp_path / 'interface.json'
    record, _ = reviewed_record()
    record['state'] = 'implementation_only'
    k.write(path, record)
    def forbidden(*args, **kwargs):
        pytest.fail('No work may start before the reviewed practical interface')
    monkeypatch.setattr(study.k, 'committed_sources', forbidden)
    monkeypatch.setattr(study, 'References', forbidden)
    with pytest.raises(ValueError, match='Reviewed practical interface required'):
        study.prepare(tmp_path / 'missing-references', tmp_path / 'output', path)
    assert not (tmp_path / 'output').exists()


def test_reviewed_attestation_also_requires_actual_matching_git_source(tmp_path, monkeypatch):
    path = tmp_path / 'interface.json'
    record, blob = reviewed_record()
    k.write(path, record)
    def git(args, **kwargs):
        assert args[:4] == ['git', '-C', str(k.ROOT), 'show']
        assert args[4].startswith(record['head'] + ':')
        return blob
    monkeypatch.setattr(study.subprocess, 'check_output', git)
    monkeypatch.setattr(study.k, 'sha', lambda path: hashlib.sha256(blob).hexdigest())
    assert study.practical_interface(path) == record
    record['source_sha256']['numerical_audit/practical_accuracy/construction.py'] = '0' * 64
    k.write(path, record)
    with pytest.raises(ValueError, match='Practical reviewed Git source differs'):
        study.practical_interface(path)


def minimal_library(tmp_path, raw):
    root = tmp_path / 'references'
    model_path = root / 'inputs_evidence/inputs/toy/model.json'
    k.write(model_path, dict(n=4))
    path = root / 'states/toy.npz'
    path.parent.mkdir()
    with zipfile.ZipFile(path, 'w') as stream:
        stream.writestr('vectors.npy', raw)
    reference = study.References.__new__(study.References)
    reference.root, reference.cases = root, ['toy']
    reference.records = {'toy': {'vectors_file': 'states/toy.npz', 'vectors_sha256': k.sha(path)}}
    reference.vectors = {}
    return reference


def test_declared_large_npy_shape_refused_before_numpy_load(tmp_path, monkeypatch):
    raw = io.BytesIO()
    np.lib.format.write_array_header_1_0(raw, dict(descr=np.dtype('complex128').str,
                                                fortran_order=False, shape=(1 << 40, 5)))
    reference = minimal_library(tmp_path, raw.getvalue())
    def forbidden(*args, **kwargs):
        pytest.fail('Array allocation must follow the bounded header check')
    monkeypatch.setattr(study.np, 'load', forbidden)
    with pytest.raises(ValueError, match='Reference array shape/dtype before allocation'):
        reference.vector('toy', 'HF')


def test_object_npy_dtype_refused_before_numpy_load(tmp_path, monkeypatch):
    raw = io.BytesIO()
    np.lib.format.write_array_header_1_0(raw, dict(descr=np.dtype('object').str,
                                                fortran_order=False, shape=(16, 5)))
    reference = minimal_library(tmp_path, raw.getvalue())
    def forbidden(*args, **kwargs):
        pytest.fail('Object arrays cannot reach NumPy load')
    monkeypatch.setattr(study.np, 'load', forbidden)
    with pytest.raises(ValueError, match='Reference array shape/dtype before allocation'):
        reference.vector('toy', 'HF')


def test_bounded_reference_archive_returns_independent_state_copy(tmp_path):
    values = np.zeros((16, 5), complex)
    values[3, :] = 1
    raw = io.BytesIO()
    np.save(raw, values, allow_pickle=False)
    reference = minimal_library(tmp_path, raw.getvalue())
    first = reference.vector('toy', 'HF')
    first[3] = 0
    assert reference.vector('toy', 'HF')[3] == 1
    assert np.array_equal(reference.vector('toy', 'coherent_imag'), values[:, 4])


def test_outer_training_charge_includes_final_constructor_return(monkeypatch):
    class Reference:
        inputs = {'successes': {'toy': 'expected-model-sha'}}
        def model_path(self, case): return Path('/not-read-by-fake-constructor')
        def vector(self, case, state): return np.eye(16, dtype=complex)[:, 3]
        def cost(self, case, state): return .5
    times = iter([0., 8.])
    monkeypatch.setattr(study.time, 'process_time', lambda: next(times))
    def first(path, digest, psi, partition, reference_CPU_charge):
        assert digest == 'expected-model-sha' and partition == 'conserving' and reference_CPU_charge == .5
        return dict(actual_local_CPU_seconds=1., record={'order': [0]}, phases={})
    monkeypatch.setattr(study.training, 'first', first)
    result = study.constructor(Reference(), dict(case='toy', state='HF', method='first_JD'))
    assert result['caller_CPU_seconds'] == 8.
    assert result['actual_local_CPU_seconds'] == 1.
    assert result['standalone_reference_and_ordering_CPU_charge'] == 8.5


def test_frozen_constructor_tampering_rejected_before_screening(tmp_path):
    plan, dependency = tmp_path / 'PLAN_FROZEN.json', tmp_path / 'PRACTICAL_INTERFACE.json'
    k.write(plan, {'logical': []})
    k.write(dependency, {'head': 'fixture'})
    constructor = tmp_path / 'constructors/order.json'
    k.write(constructor, {'order': [0, 1]})
    k.write(tmp_path / 'PROTOCOL_FROZEN.json', dict(plan_sha256=k.sha(plan),
            reviewed_practical_interface_sha256=k.sha(dependency)))
    k.write(tmp_path / 'SCHEDULES_FROZEN.json', dict(artifacts=study.inventory(tmp_path)))
    study.bind(tmp_path)
    k.write(constructor, {'order': [1, 0]})
    with pytest.raises(ValueError, match='Frozen transfer construction changed'):
        study.bind(tmp_path)


def test_symlink_evidence_refused_without_following_it(tmp_path):
    (tmp_path / 'link').symlink_to('/not-a-task-owned-artifact')
    with pytest.raises(ValueError, match='Evidence symlink refused'):
        study.inventory(tmp_path)
