"""Source/input contracts; no chemistry or molecular quality is measured."""
from fractions import Fraction
import json
from pathlib import Path
import sys
import tomllib

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import study_common as k


def model():
    # Integer larger than binary64's exact integer range: serialization must
    # preserve its last bit rather than round through a numerical coefficient.
    numerator = 2**100 + 1
    words = [[[[0, 'Z']], str(numerator)], [[[13, 'Z']], '3']]
    return dict(case='contract_only', n=14, electrons=4, denominator=str(2**107),
                coefficients=words, blocks=[[words[0]], [words[1]]])


def test_exact_json_to_toml_and_original_labels(tmp_path):
    data = model()
    target = tmp_path / 'model.toml'
    k.write_model_toml(target, data, source_json_sha256='a' * 64, order=[1, 0])
    loaded = tomllib.loads(target.read_text())
    assert loaded['baseline_order'] == [2, 1]
    assert loaded['denominator'] == str(2**107)
    assert loaded['units'][0]['numerators'] == [str(2**100 + 1)]
    assert loaded['units'][1]['words'] == ['13Z']
    assert k.blocks(data)[0][((0, 'Z'),)] == Fraction(2**100 + 1, 2**107)
    with pytest.raises(FileExistsError):
        k.write_model_toml(target, data, source_json_sha256='b' * 64)
    assert tomllib.loads(target.read_text())['source_json_sha256'] == 'a' * 64


def test_partition_and_register_identity_are_not_inferred_from_labels():
    data = model()
    data['blocks'][1] = data['blocks'][0]
    with pytest.raises(ValueError, match='partition'):
        k.blocks(data)
    data = model()
    data['n'] = 12
    with pytest.raises(ValueError, match='fourteen'):
        k.blocks(data)


def test_result_publication_cannot_replace_existing_evidence(tmp_path):
    target = tmp_path / 'result.json'
    k.write_new(target, dict(original=True))
    before = target.read_bytes()
    with pytest.raises(ValueError):
        k.write_new(target, dict(original=False))
    assert target.read_bytes() == before


def test_retained_control_has_unchanged_bounds_and_scoped_import_alias():
    sentinel = object()
    previous = sys.modules.get('common')
    sys.modules['common'] = sentinel
    try:
        control = k.retained_training()
        assert sys.modules['common'] is sentinel
        assert control.k.MAX_TERMS == 4096
        assert Path(control.current.__code__.co_filename).resolve() == (
            k.ROOT / 'numerical_audit/scalable_ordering/training.py')
        assert control.k.ROOT == k.ROOT
    finally:
        if previous is None:
            sys.modules.pop('common', None)
        else:
            sys.modules['common'] = previous


def test_predata_inventory_counts_include_pauli_and_first_order_controls():
    registry = json.loads((HERE / 'REGISTRY.json').read_text())
    assert len(registry['cases']) == 4 and {case['n'] for case in registry['cases']} == {14}
    assert len(k.ARMS) == 8 and 'pauli' in k.ARMS
    assert len(registry['cases']) * len(k.ARMS) * k.REPETITIONS == 64
    logical = len(registry['cases']) * k.REPETITIONS * len(k.TIMES) * len(k.CAPS) * sum(
        len(k.FORMULAS[arm]) for arm in k.ARMS)
    assert logical == 352
