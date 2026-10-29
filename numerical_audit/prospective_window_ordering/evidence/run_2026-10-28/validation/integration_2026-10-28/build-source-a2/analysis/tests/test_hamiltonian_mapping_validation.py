"""Reject unsupported mappings before dispatch or filename construction."""

from types import SimpleNamespace

import pytest

from qhat.hamiltonian_generator import hamgen
from qhat.hamiltonian_generator.hamgen_types import HamiltonianConfiguration, State


@pytest.mark.parametrize('mapping', ['unsupported', 'JW_typo', None])
@pytest.mark.parametrize('accessor', ['fermion_to_qubit_name', 'f2q_tag'])
def test_unknown_mapping_is_rejected_by_configuration(mapping, accessor):
    config = HamiltonianConfiguration()
    config.f2q_mapping = mapping
    with pytest.raises(NotImplementedError, match='invalid fermion-to-qubit mapping'):
        getattr(config, accessor)()


@pytest.mark.parametrize('mapping,name,tag', [
    ('Jordan-Wigner', 'jordan-wigner', 'jw'),
    ('Jordan Wigner', 'jordan-wigner', 'jw'),
    ('JW', 'jordan-wigner', 'jw'),
    ('jordan-wigner', 'jordan-wigner', 'jw'),
    ('jordan_wigner', 'jordan-wigner', 'jw'),
    ('jordan wigner', 'jordan-wigner', 'jw'),
    ('jw', 'jordan-wigner', 'jw'),
    ('Bravyi-Kitaev', 'bravyi-kitaev', 'bk'),
    ('Bravyi Kitaev', 'bravyi-kitaev', 'bk'),
    ('BK', 'bravyi-kitaev', 'bk'),
    ('bravyi-kitaev', 'bravyi-kitaev', 'bk'),
    ('bravyi_kitaev', 'bravyi-kitaev', 'bk'),
    ('bravyi kitaev', 'bravyi-kitaev', 'bk'),
    ('bk', 'bravyi-kitaev', 'bk'),
])
def test_supported_mapping_aliases_keep_their_names_and_tags(mapping, name, tag):
    config = HamiltonianConfiguration()
    config.f2q_mapping = mapping
    assert config.fermion_to_qubit_name() == name
    assert config.f2q_tag() == tag


def test_mapping_dispatch_reports_the_configured_string():
    config = SimpleNamespace(
        fermion_to_qubit_name=lambda: 'unsupported', f2q_mapping='bad mapping')
    state = SimpleNamespace(config_hamiltonian=config)
    with pytest.raises(NotImplementedError, match='bad mapping'):
        hamgen.map_fermions_to_qubits(state, object())


def test_unknown_mapping_cannot_silently_create_a_none_filename_tag():
    config = HamiltonianConfiguration()
    config.f2q_mapping = 'unsupported'
    config.num_active_occupied = config.num_active_vacant = 1
    state = State.__new__(State)
    state.config_hamiltonian = config
    state.config_general = SimpleNamespace(file_stub='toy', ham3_ext=lambda: 'npz')
    with pytest.raises(NotImplementedError, match='invalid fermion-to-qubit mapping'):
        state.filename_ham3()
