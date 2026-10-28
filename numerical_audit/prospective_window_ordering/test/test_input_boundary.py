"""The live chemistry encoder must satisfy the strict reader before JSON IO."""
from fractions import Fraction
import json
from pathlib import Path
import sys
import tomllib

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import generate_inputs as generator
import study_common as k


def test_tuple_pauli_keys_survive_in_memory_json_and_toml(tmp_path):
    total = {((13, 'Z'),): -3, ((0, 'Z'),): 2**100 + 1}
    encoded = generator.encode_integer_terms(total)
    assert all(type(word) is list and all(type(pair) is list for pair in word)
               for word, coefficient in encoded)
    model = dict(case='encoder_boundary', n=14, electrons=4, denominator=str(2**107),
                 coefficients=encoded,
                 blocks=[generator.encode_integer_terms({word: total[word]}) for word in sorted(total)])
    expected = ({((0, 'Z'),): Fraction(2**100 + 1, 2**107)},
                {((13, 'Z'),): Fraction(-3, 2**107)})
    # The failing producer performs this check BEFORE its json.dump call.
    assert k.blocks(model) == expected
    assert k.blocks(json.loads(json.dumps(model))) == expected
    output = tmp_path / 'model.toml'
    k.write_model_toml(output, model, source_json_sha256='a' * 64)
    restored = tomllib.loads(output.read_text())
    assert restored['units'][0]['numerators'] == [str(2**100 + 1)]
    assert restored['units'][1]['words'] == ['13Z']
