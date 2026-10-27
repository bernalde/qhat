"""Construct a sparse response ordering from an authenticated-model JSON shape.

Example: python -m qhat.analysis.train_ordering model.json --time 1 \
    --rotation-budget 1000 --cpu-seconds 2 --error-tolerance 1e-3

The model must supply n, electrons, denominator and blocks in QHAT research-model
notation: each block is a list of [Pauli-word, integer numerator] pairs. A word is
a list of [qubit index, X/Y/Z] pairs, with qubit zero the least significant bit.
No chemistry calculation, propagation or hidden validation oracle is invoked.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

from qhat.analysis.sparse_response_ordering import (
    ResourceLimit, conservative_steps, train_ordering,
)


def load_model(path):
    path = Path(path)
    with path.open('rb') as handle:
        raw = handle.read(10 * 1024**2 + 1)
    if len(raw) > 10 * 1024**2:
        raise ResourceLimit('model JSON exceeds 10 MiB')
    model = json.loads(raw)
    if not isinstance(model, dict):
        raise ValueError('model must be a JSON object')
    n, electrons = model['n'], model['electrons']
    if type(n) is not int or not 1 <= n <= 4096:
        raise ValueError('n must be an integer in 1..4096')
    if type(electrons) is not int or not 0 <= electrons <= n:
        raise ValueError('invalid active-electron count')
    def integer(value):
        if type(value) is int:
            return value
        if isinstance(value, str) and re.fullmatch(r'[+-]?[0-9]+', value):
            return int(value)
        raise ValueError('coefficient integers must not be truncated or coerced from booleans')
    denominator = integer(model['denominator'])
    if denominator <= 0 or denominator.bit_length() > 4096:
        raise ValueError('invalid coefficient denominator')
    blocks = []
    if not isinstance(model['blocks'], list) or len(model['blocks']) > 4096:
        raise ResourceLimit('invalid or oversized block list')
    terms = 0
    for raw_block in model['blocks']:
        if not isinstance(raw_block, list):
            raise ValueError('block must be a list')
        terms += len(raw_block)
        if terms > 50000:
            raise ResourceLimit('too many Pauli terms')
        block = {}
        for raw_word, raw_coefficient in raw_block:
            word = tuple(tuple(letter) for letter in raw_word)
            if word in block:
                raise ValueError('duplicate Pauli word within block')
            coefficient = integer(raw_coefficient)
            if coefficient.bit_length() > 4096:
                raise ValueError('oversized coefficient numerator')
            block[word] = Fraction(coefficient, denominator)
        blocks.append(block)
    return model, tuple(blocks), hashlib.sha256(raw).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path)
    parser.add_argument('--time', default='1', help='evolution time with hbar=1')
    parser.add_argument('--rotation-budget', required=True, type=int,
                        help='ideal Pauli-rotation budget used by the response objective')
    parser.add_argument('--cpu-seconds', type=float, default=2.)
    parser.add_argument('--max-evaluations', type=int, default=3000)
    parser.add_argument('--max-entries', type=int, default=50000)
    parser.add_argument('--error-tolerance', help='optional conservative ideal-error target')
    parser.add_argument('--metric', choices=['operator', 'infidelity', 'squared_infidelity'],
                        default='operator')
    parser.add_argument('--max-certified-steps', type=int, default=10**9)
    parser.add_argument('--output', type=Path, help='write a new JSON file; existing files are refused')
    args = parser.parse_args(argv)
    model, blocks, digest = load_model(args.model)
    duration = Fraction(args.time)
    reference = (1 << model['electrons']) - 1
    trained = train_ordering(blocks, model['n'], reference, time=duration,
                             rotation_budget=args.rotation_budget, cpu_seconds=args.cpu_seconds,
                             max_evaluations=args.max_evaluations, max_entries=args.max_entries)
    output = dict(model_sha256=digest, n_qubits=model['n'], reference='Hartree-Fock determinant',
                  time=str(duration), training=asdict(trained),
                  score_scope='estimated diagonal-reference cubic response; not a finite-error bound',
                  rotation_budget=args.rotation_budget,
                  certificate_scope='ideal supplied coefficients and exact unit exponentials; excludes chemistry, synthesis and hardware error')
    if args.error_tolerance is not None:
        try:
            certificate = conservative_steps(blocks, model['n'], time=duration,
                                              error_tolerance=Fraction(args.error_tolerance),
                                              metric=args.metric, max_steps=args.max_certified_steps)
        except ResourceLimit as exc:
            output['certificate'] = dict(status='unresolved_within_step_limit', reason=str(exc))
        else:
            affordable = trained.score.steps if trained.score is not None else None
            output['certificate'] = dict(status='certified_ideal_bound', steps=certificate.steps,
                                         metric=certificate.metric,
                                         error_bound=str(certificate.error_bound),
                                         operator_error_bound=str(certificate.operator_error_bound),
                                         fits_rotation_budget=(certificate.steps <= affordable
                                                               if affordable is not None else None),
                                         order_dependent=False)
    text = json.dumps(output, indent=2) + '\n'
    if args.output is not None:
        with args.output.open('x') as handle:
            handle.write(text)
    else:
        print(text, end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
