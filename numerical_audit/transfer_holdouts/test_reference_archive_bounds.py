"""Regressions through the actual independent reference auditor entrypoints.

The original review target remains frozen in a separate worktree.
Malformed header fixtures never allocate the declared oversized arrays: a spy
fails if they reach NumPy's loader before controlled archive rejection.
"""
import hashlib
import io
import json
from pathlib import Path
import shutil
import struct
import zipfile
import numpy as np
import pytest
import verify_references as reference


def _guard_sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _guard_npy(descr,shape,payload=b''):
    text=repr(dict(descr=descr,fortran_order=False,shape=shape))
    header=(text+' '*((-10-len(text)-1)%64)+'\n').encode('ascii')
    return b'\x93NUMPY\x01\x00'+struct.pack('<H',len(header))+header+payload

def _guard_fixture(tmp_path):
    np.savez_compressed(tmp_path/'tensors.npz',constant=0.,one_body=np.zeros((6,6)),two_body=np.zeros((6,6,6,6)))
    np.savez_compressed(tmp_path/'orbitals.npz',mo_coeff=np.eye(6),ao_overlap=np.eye(6),mo_occ=np.array([2.,2.,0.,0.,0.,0.]),mo_energy=np.arange(6,dtype=float),ao_hcore=np.zeros((6,6)))
    model=dict(n=6,electrons=2,total_electrons=4,basis='sto-6g',family='LiH',case='LiH_r1p713',geometry=[['Li',[0.,0.,0.]],['H',[0.,0.,1.713]]],denominator='1',removed_scalar_exact='0',hf_energy=0.,fermionic_parents=[],policy=dict(canonical_hermitian_adjustment_fermionic_l1=0.,discarded_parent_fermionic_l1=0.,input_constant=0.,parent_cutoff=1e-12,removed_parent_pauli_scalar=0.))
    return model

def _guard_hashes(folder,model):
    model['tensors_sha256']=_guard_sha(folder/'tensors.npz');model['orbitals_sha256']=_guard_sha(folder/'orbitals.npz')

def _guard_corrupt(path,kind,member):
    with zipfile.ZipFile(path) as z:entries={i.filename:z.read(i) for i in z.infolist()}
    shape=(6,6,6,6) if member=='two_body.npy' else () if member=='constant.npy' else (6,) if member in ('mo_occ.npy','mo_energy.npy') else (6,6)
    if kind=='shape':bad=_guard_npy('<f8',(1_000_000_000,6))
    elif kind=='dtype':bad=_guard_npy('|O',shape)
    elif kind=='uncompressed':bad=_guard_npy('<f8',shape,b'\0'*40_000)
    elif kind=='payload':bad=_guard_npy('<f8',shape,b'\0')
    elif kind=='header':bad=b'\x93NUMPY\x02\x00'+struct.pack('<I',2**32-1)+b'{}'
    else:bad=_guard_npy('<f8',(),b'\0'*8)
    entries[member]=bad
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,data in entries.items():z.writestr(name,data)
    return bad

def _guard_loader_spy(monkeypatch,bad_path,bad_payload):
    real=np.load;hits=[];compressed=bad_path.read_bytes()
    def load(source,*args,**kwargs):
        raw=source.getvalue() if isinstance(source,io.BytesIO) else None
        path=Path(source) if isinstance(source,(str,Path)) else None
        if path==bad_path or raw==compressed or (bad_payload is not None and raw==bad_payload):
            hits.append(True);raise AssertionError('Unsafe NumPy allocation path reached before archive validation')
        return real(source,*args,**kwargs)
    monkeypatch.setattr(reference.np,'load',load)
    return hits

def test_actual_molecular_archive_auditor_accepts_bounded_normal_fixture(tmp_path):
    model=_guard_fixture(tmp_path);_guard_hashes(tmp_path,model)
    identity,_=reference.molecular_identity(tmp_path,model,np.zeros((64,64)))
    assert identity['raw_tensor_target_difference']==0.
    assert identity['orbital_orthogonality_error']==0.

@pytest.mark.parametrize('kind,member',[('shape','one_body.npy'),('shape','two_body.npy'),('shape','constant.npy'),
    ('dtype','two_body.npy'),('dtype','constant.npy'),('uncompressed','two_body.npy'),
    ('inventory','unexpected.npy'),('payload','one_body.npy'),('header','two_body.npy')])
def test_actual_tensor_auditor_refuses_before_numpy_load(tmp_path,monkeypatch,kind,member):
    model=_guard_fixture(tmp_path);path=tmp_path/'tensors.npz';bad=_guard_corrupt(path,kind,member);_guard_hashes(tmp_path,model);hits=_guard_loader_spy(monkeypatch,path,bad)
    with pytest.raises(ValueError,match='(?i)(archive|NPY|shape|dtype|byte|member)'):
        reference.molecular_identity(tmp_path,model,np.zeros((64,64)))
    assert not hits

@pytest.mark.parametrize('kind,member',[('shape','mo_coeff.npy'),('shape','ao_overlap.npy'),('dtype','mo_occ.npy'),
    ('shape','mo_energy.npy'),('shape','ao_hcore.npy')])
def test_actual_orbital_auditor_refuses_before_numpy_load(tmp_path,monkeypatch,kind,member):
    model=_guard_fixture(tmp_path);path=tmp_path/'orbitals.npz';bad=_guard_corrupt(path,kind,member);_guard_hashes(tmp_path,model);hits=_guard_loader_spy(monkeypatch,path,bad)
    with pytest.raises(ValueError,match='(?i)(archive|NPY|shape|dtype|byte|member)'):
        reference.molecular_identity(tmp_path,model,np.zeros((64,64)))
    assert not hits

@pytest.mark.parametrize('kind',['shape','dtype'])
def test_actual_reference_entry_rejects_before_vector_allocation(tmp_path,monkeypatch,kind):
    source=Path(reference.__file__).resolve().parent/'evidence/results/reference_evidence'
    assert source.is_dir(),'The complete committed numerical reference payload is required'
    root=tmp_path/'reference';shutil.copytree(source,root)
    registry=reference.k.read(reference.k.HERE/'REGISTRY.json');case=registry['cases'][0]['case']
    record_path=root/'references'/(case+'.json');record=reference.k.read(record_path);vectors=root/record['vectors_file']
    bad=_guard_npy('<c16',(64,1_000_000_000)) if kind=='shape' else _guard_npy('|O',(64,5))
    with zipfile.ZipFile(vectors,'w',compression=zipfile.ZIP_DEFLATED) as z:z.writestr('vectors.npy',bad)
    record['vectors_sha256']=_guard_sha(vectors);record_path.write_text(json.dumps(record,indent=2)+'\n')
    frozen_path=root/'REFERENCE_STATES_FROZEN.json';frozen=reference.k.read(frozen_path)
    for path in (vectors,record_path):frozen['artifacts'][str(path.relative_to(root))]=_guard_sha(path)
    frozen_path.write_text(json.dumps(frozen,indent=2)+'\n')
    monkeypatch.setattr(reference.geometry_identity,'audit',lambda registry:dict(fixture_only='Unrelatedbaselinegeometryauditomittedfromthisallocation-pathregression'))
    hits=_guard_loader_spy(monkeypatch,vectors,bad)
    with pytest.raises(ValueError,match='(?i)(archive|NPY|shape|dtype|byte|member)'):
        reference.verify(root,receipt_dir=tmp_path/'independent-audit')
    assert not hits


@pytest.mark.parametrize('kind', ['compressed_bytes', 'duplicate_member', 'missing_member'])
def test_actual_tensor_auditor_rejects_archive_container_before_load(tmp_path, monkeypatch, kind):
    model = _guard_fixture(tmp_path)
    path = tmp_path / 'tensors.npz'
    with zipfile.ZipFile(path) as archive:
        entries = {item.filename: archive.read(item) for item in archive.infolist()}
    if kind == 'compressed_bytes':
        path.write_bytes(path.read_bytes() + b'x' * 40_000)
    else:
        with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in entries.items():
                if kind != 'missing_member' or name != 'constant.npy':
                    archive.writestr(name, data)
            if kind == 'duplicate_member':
                with pytest.warns(UserWarning, match='Duplicate name'):
                    archive.writestr('one_body.npy', entries['one_body.npy'])
    _guard_hashes(tmp_path, model)
    hits = _guard_loader_spy(monkeypatch, path, None)
    with pytest.raises(ValueError, match='(?i)(archive|byte|member)'):
        reference.molecular_identity(tmp_path, model, np.zeros((64, 64)))
    assert not hits


def test_actual_auditor_checks_every_tensor_header_before_first_numpy_allocation(tmp_path, monkeypatch):
    model = _guard_fixture(tmp_path)
    _guard_corrupt(tmp_path / 'tensors.npz', 'shape', 'two_body.npy')
    _guard_hashes(tmp_path, model)
    calls = []
    def refuse(*args, **kwargs):
        calls.append(True)
        raise AssertionError('An array was loaded before all tensor headers passed')
    monkeypatch.setattr(reference.np, 'load', refuse)
    with pytest.raises(ValueError, match='(?i)(archive|shape)'):
        reference.molecular_identity(tmp_path, model, np.zeros((64, 64)))
    assert not calls


def test_actual_orbital_bounds_follow_sto6g_atom_inventory(tmp_path, monkeypatch):
    model = _guard_fixture(tmp_path)
    model.update(family='BeH2', total_electrons=6,
                 geometry=[['Be', [0., 0., 0.]], ['H', [0., 0., -1.2]], ['H', [0., 0., 1.2]]])
    # These otherwise well-formed arrays have six orbitals; BeH2 STO-6G has seven.
    path = tmp_path / 'orbitals.npz'
    np.savez_compressed(path, mo_coeff=np.eye(6), ao_overlap=np.eye(6),
                        mo_occ=np.array([2., 2., 2., 0., 0., 0.]),
                        mo_energy=np.arange(6, dtype=float), ao_hcore=np.zeros((6, 6)))
    _guard_hashes(tmp_path, model)
    with zipfile.ZipFile(path) as archive:
        bad = archive.read('mo_coeff.npy')
    hits = _guard_loader_spy(monkeypatch, path, bad)
    with pytest.raises(ValueError, match='(?i)(archive|shape)'):
        reference.molecular_identity(tmp_path, model, np.zeros((64, 64)))
    assert not hits


@pytest.mark.parametrize('change', ['basis', 'atoms'])
def test_actual_auditor_rejects_unsupported_orbital_domain_before_load(tmp_path, monkeypatch, change):
    model = _guard_fixture(tmp_path)
    if change == 'basis':
        model['basis'] = 'unbounded-user-supplied-basis'
    else:
        model['geometry'][0][0] = 'N'
    _guard_hashes(tmp_path, model)
    def refuse(*args, **kwargs):
        raise AssertionError('Unsupported chemistry reached NumPy allocation')
    monkeypatch.setattr(reference.np, 'load', refuse)
    with pytest.raises(ValueError, match='(?i)(archive|domain|inventory)'):
        reference.molecular_identity(tmp_path, model, np.zeros((64, 64)))
