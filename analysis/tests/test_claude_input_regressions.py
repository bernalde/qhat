"""Behavioral regressions for supported loaders, dispatch, and phase sizing."""
from types import SimpleNamespace
import builtins
import json
import math
import h5py
import networkx as nx
import numpy as np
import pytest
from qhat.analysis.hamiltonian import load_pauli, load_hamlib_hdf5, get_physical_hamiltonian, boson_to_qubit_operator
from qhat.analysis.phase_parameters import phase_qubits_for_time


@pytest.mark.parametrize('fmt', ['dense', 'sparse', 'json', 'hamlib'])
def test_repeated_pauli_inputs_accumulate_before_validation(tmp_path, fmt):
    if fmt == 'dense':
        p=tmp_path/'terms.txt';p.write_text('2+1j XZ\n3-1j XZ\n1 II\n-1 II\n')
    elif fmt == 'sparse':
        p=tmp_path/'terms.txt';p.write_text('(2+1j) [X0 Z1]\n(3-1j) [Z1 X0]\n(1) []\n(-1) []\n')
    elif fmt == 'json':
        p=tmp_path/'terms.json';p.write_text(json.dumps(dict(n_qubits=2, terms=[
            dict(ops=[[0,'X'],[1,'Z']], coeff=2),dict(ops=[[1,'Z'],[0,'X']], coeff=3),
            dict(ops=[],coeff=1),dict(ops=[],coeff=-1)])))
    else:
        p=tmp_path/'terms.h5'
        with h5py.File(p,'w') as f:
            f.create_dataset('H',data=b'(2+1j) [X0 Z1] +\n(3-1j) [Z1 X0] +\n(1+0j) [] +\n(-1+0j) []')
    cfg=SimpleNamespace(filename=str(p),hdf5_key='H')
    h=(load_hamlib_hdf5 if fmt=='hamlib' else load_pauli)(cfg)
    assert h.get_all_pauli_strings(return_as='strings') == {'XZ':5,'II':0}


def test_nonreal_accumulated_coefficient_is_rejected(tmp_path):
    p=tmp_path/'bad.txt';p.write_text('1+1j X\n1 X\n')
    with pytest.raises(ValueError,match='Hermitian'):
        load_pauli(SimpleNamespace(filename=str(p)))


def test_unsupported_lcps_and_boson_stub_fail_clearly():
    with pytest.raises(ValueError,match='supported sources'):
        get_physical_hamiltonian(SimpleNamespace(source='LCPS'))
    with pytest.raises(NotImplementedError,match='BosonicBinaryEncoding'):
        boson_to_qubit_operator(None, 2)


def test_optional_gurobi_import_has_clear_failure(monkeypatch):
    from qhat.analysis.ordering import chromatic_number_and_coloring_gurobi
    original=builtins.__import__
    def without_gurobi(name,*args,**kwargs):
        if name=='gurobipy':
            raise ImportError('deliberate missing optional dependency')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',without_gurobi)
    with pytest.raises(ImportError,match='optional gurobipy'):
        chromatic_number_and_coloring_gurobi(nx.complete_graph(3))


@pytest.mark.parametrize('time', [0, -1, float('inf'), float('nan'), math.pi])
def test_invalid_or_wrapping_time_is_rejected(time):
    with pytest.raises(ValueError):
        phase_qubits_for_time(time,-1,1,1.01,.1)


def test_precision_uses_actual_safe_time_and_margin():
    assert phase_qubits_for_time(1,-1,1,1.01,.1)==6
    assert phase_qubits_for_time(.5,-1,1,1.01,.1)==7
    assert phase_qubits_for_time(.5,-1,1,1.01,None) is None
    with pytest.raises(ValueError,match='margin'):
        phase_qubits_for_time(math.pi/1.005,-1,1,1.01,.1)


@pytest.mark.parametrize('method', ['', 'r', 'ramped', 'trotter'])
def test_substring_unitary_methods_are_rejected(method):
    from qhat.analysis.unitary import encode_as_unitary
    with pytest.raises(ValueError,match='Invalid unitary encoding method'):
        encode_as_unitary(SimpleNamespace(method=method),None,None)


def test_dispatch_normalizes_case_and_whitespace(monkeypatch):
    import qhat.analysis.unitary as unitary
    sentinel=object()
    monkeypatch.setattr(unitary,'encode_ramped_trotter',lambda *args:sentinel)
    assert unitary.encode_as_unitary(SimpleNamespace(method=' Ramped Trotter '),None,1) is sentinel


@pytest.mark.parametrize('ordering,expected', [('magnitude',['Z','X']),('lexicographical',['X','Z'])])
def test_estimator_and_builder_receive_identical_order(monkeypatch,ordering,expected):
    import qhat.analysis.unitary as unitary
    import qhat.analysis.trotter_coefficients_fast as fast
    import qhat.common.trotter_flattened as builder
    from qhat.analysis.config_types import UnitaryConfiguration
    from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings
    cfg=UnitaryConfiguration();cfg.encode_ramped_trotter(energy_error=1,
        ordering_method=ordering,trotter_order='second order')
    h=Hamiltonian(LinearCombinationOfPauliStrings(dense={'Z':2,'X':1},num_qubits=1))
    seen={}
    def estimate(terms,*args,**kwargs):
        seen['estimated']=[next(iter(t.terms)) for t in terms]
        return 2,5/3
    def build(terms,*args,**kwargs):
        seen['built']=list(terms)
        return seen
    monkeypatch.setattr(fast,'trotter_error_estimator_fast',estimate)
    monkeypatch.setattr(builder,'build_ramped_trotterized_unitary',build)
    result=unitary.encode_ramped_trotter(cfg,h,1)
    assert [word for word,_ in result['built']]==expected
    assert result['estimated']==[((0,word),) for word in expected]


def test_driver_carries_canonical_method_and_override_precision(monkeypatch):
    import qhat.analysis.driver as driver
    observed={}
    general=SimpleNamespace(get_output_path=lambda name:name,logfile='unused.log',
                            loglevel='warning',git_hash='fixture')
    unitary=SimpleNamespace(method=' Ramped Trotter ',phase_scale_factor=1.01,timestep=.5)
    algorithm=SimpleNamespace(energy_error=.1)
    state=SimpleNamespace(config_general=general,config_unitary=unitary,
        config_algorithm=algorithm,config_hamiltonian=object(),config_analysis=object(),
        store_results=lambda results:None,show_results=lambda:None,save_summary=lambda:None)
    h=SimpleNamespace(compute_initial_energy_bounds=lambda cfg:(-1,1),
                      energy_shift=lambda shift:None,get_energy_shift=lambda:0)
    monkeypatch.setattr(driver,'load_configuration',lambda:(state,'fixture',None))
    monkeypatch.setattr(driver,'configure_logging',lambda **kwargs:None)
    monkeypatch.setattr(driver,'get_physical_hamiltonian',lambda cfg:h)
    monkeypatch.setattr(driver,'compute_initial_phase_qubits',lambda *args:(3,-1,1))
    def encode(cfg,physical,time):
        observed.update(method=cfg.method,time=time)
        return object()
    def build(cfg,encoded,p0):
        observed['precision']=p0
        return object()
    def analyze(cfg,built,method,**kwargs):
        observed['analysis_method']=method
        return {}
    monkeypatch.setattr(driver,'encode_as_unitary',encode)
    monkeypatch.setattr(driver,'build_algorithm',build)
    monkeypatch.setattr(driver,'analyze_algorithm',analyze)
    driver.run()
    assert observed==dict(method='ramped trotter',analysis_method='ramped trotter',time=.5,precision=7)
