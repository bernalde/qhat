"""Exact altered-input selection, preserving original source evidence."""
from fractions import Fraction
import json
from pathlib import Path
import sys
import tomllib

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import study_common as k
import preflight


def source_inventory(directory):
    registry = k.read(HERE / 'REGISTRY.json')
    directory.mkdir()
    k.write_new(directory / 'INPUT_PROTOCOL_FROZEN.json', dict(sources=k.source_map()))
    successes = {}
    for spec in registry['cases']:
        target = directory / 'inputs' / spec['case']
        target.mkdir(parents=True)
        words = [[[[0, 'Z']], str(2**100 + 1)], [[[13, 'Z']], '-3']]
        data = dict(spec, basis=registry['basis'], unit=registry['unit'], charge=0, spin=0,
                    denominator=str(2**107), coefficients=words, blocks=[[words[0]], [words[1]]],
                    removed_scalar_exact='-3/7')
        k.write_new(target / 'model.json', data)
        successes[spec['case']] = {'model.json': k.sha(target / 'model.json')}
    k.write_new(directory / 'INPUTS_FROZEN.json', dict(status='complete', expected_cases=4,
        successes=successes, refusals={}, protocol_sha256=k.sha(directory / 'INPUT_PROTOCOL_FROZEN.json')))
    return registry


def test_exact_scaling_and_lexical_tie_preserve_originals(tmp_path):
    source = tmp_path / 'source'
    registry = source_inventory(source)
    before = {path: k.sha(path) for path in source.rglob('*') if path.is_file()}
    output = tmp_path / 'altered'
    preflight.alter(source, output)
    inventory = k.read(output / 'INPUTS_FROZEN.json')
    selected = min(spec['case'] for spec in registry['cases'])
    assert inventory['selected_case'] == selected
    name = selected + '_preflight_7over5'
    model = k.read(output / 'inputs' / name / 'model.json')
    assert k.blocks(model)[0][((0, 'Z'),)] == Fraction(7 * (2**100 + 1), 5 * 2**107)
    assert k.blocks(model)[1][((13, 'Z'),)] == Fraction(-21, 5 * 2**107)
    assert model['removed_scalar_exact'] == '-3/5'
    assert inventory['parent_inventory_sha256'] == k.sha(source / 'INPUTS_FROZEN.json')
    encoded = tomllib.loads((output / 'inputs' / name / 'model.toml').read_text())
    assert encoded['units'][0]['numerators'] == [str(7 * (2**100 + 1))]
    assert encoded['source_json_sha256'] == k.sha(output / 'inputs' / name / 'model.json')
    assert all(k.sha(path) == value for path, value in before.items())
    with pytest.raises(ValueError, match='Fresh'):
        preflight.alter(source, output)


def test_missing_larger_case_does_not_create_easier_preflight(tmp_path):
    source = tmp_path / 'source'
    source_inventory(source)
    path = source / 'INPUTS_FROZEN.json'
    frozen = json.loads(path.read_text())
    frozen['status'] = 'incomplete'
    path.write_text(json.dumps(frozen))
    output = tmp_path / 'must_not_exist'
    with pytest.raises(ValueError, match='Complete four-case'):
        preflight.alter(source, output)
    assert not output.exists()
