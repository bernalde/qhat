"""Select the fixed largest-word input and encode its exact 7/5 perturbation.

This is data preparation only. No constructor, error, or quality is evaluated.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import study_common as k


def alter(inputs, output):
    frozen = k.read(inputs / 'INPUTS_FROZEN.json')
    registry = k.read(k.HERE / 'REGISTRY.json')
    cases = [row['case'] for row in registry['cases']]
    k.require(frozen['status'] == 'complete' and frozen['expected_cases'] == 4 and
              set(frozen['successes']) == set(cases) and not frozen['refusals'],
              'Complete four-case source inventory required; no substituted case')
    protocol_path = inputs / 'INPUT_PROTOCOL_FROZEN.json'
    k.require(frozen['protocol_sha256'] == k.sha(protocol_path), 'Input protocol hash differs')
    protocol = k.read(protocol_path)
    k.require(protocol['sources'] == k.source_map(), 'Preflight scientific sources changed')
    admitted = []
    for case in cases:
        path = inputs / 'inputs' / case / 'model.json'
        k.require(k.sha(path) == frozen['successes'][case]['model.json'], 'Frozen model differs')
        model = k.read(path)
        k.blocks(model)
        admitted.append((len(model['coefficients']), case, path, model))
    # Largest count first, lexical case name as the declared tie rule.
    count, case, path, parent = sorted(admitted, key=lambda row: (-row[0], row[1]))[0]
    k.require(not output.exists(), 'Fresh preflight directory required')
    output.mkdir(parents=True)
    name = case + '_preflight_7over5'
    target = output / 'inputs' / name
    target.mkdir(parents=True)
    model = {key: parent[key] for key in ('n', 'electrons', 'family', 'geometry',
             'basis', 'unit', 'charge', 'spin', 'active_spatial', 'frozen_spatial')}
    model.update(case=name, denominator=str(5 * int(parent['denominator'])),
        coefficients=[[word, str(7 * int(value))] for word, value in parent['coefficients']],
        blocks=[[[word, str(7 * int(value))] for word, value in block] for block in parent['blocks']],
        preflight_parent_case=case, preflight_parent_model_sha256=k.sha(path),
        preflight_scale='7/5', preflight_only=True,
        tensors_sha256=parent['tensors_sha256'], orbitals_sha256=parent['orbitals_sha256'],
        provenance='Exact rescaling of the retained nonscalar Hamiltonian; no new SCF/tensor result')
    if 'removed_scalar_exact' in parent:
        model['removed_scalar_exact'] = str(Fraction(parent['removed_scalar_exact']) * Fraction(7, 5))
    k.blocks(model)
    k.write_new(target / 'model.json', model)
    k.write_model_toml(target / 'model.toml', model, source_json_sha256=k.sha(target / 'model.json'))
    k.write_new(output / 'INPUTS_FROZEN.json', dict(status='complete', expected_cases=1,
        preflight_only=True, protocol_sha256=k.sha(protocol_path), parent_inventory_sha256=k.sha(inputs / 'INPUTS_FROZEN.json'),
        source_protocol_path=str(protocol_path.resolve()), refusals={},
        selected_case=case, retained_word_count=count, selection_rule='largest words, lexical case tie',
        successes={name: {file: k.sha(target / file) for file in ('model.json', 'model.toml')}}))
    print('ALTERED_INPUT_FROZEN', name, 'words', count, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=k.Path)
    parser.add_argument('output', type=k.Path)
    args = parser.parse_args()
    alter(args.inputs, args.output)
