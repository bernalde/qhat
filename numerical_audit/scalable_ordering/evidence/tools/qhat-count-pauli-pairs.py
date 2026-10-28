#!/usr/bin/env python3
"""Count singleton-Pauli commutator support without storing pairs or evolving states.

This stdlib-only resource diagnostic reads exactly its supplied frozen model.
It does not import/read project source, inspect validation errors, or change a
constructor cap. Run under the coordinator's bounded operating-system guard.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import resource
import time


EXPECTED_WORDS = 1722
EXPECTED_UNORDERED_PAIRS = 1481781
MAX_MODEL_BYTES = 16 * 1024**2


def require(condition, message):
    if not condition:
        raise ValueError(message)


def update(path, result, cpu, wall):
    result.update(CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
                  max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  RSS_units='KiB on Linux; lifetime process high-water mark')
    temporary = path.with_suffix(path.suffix + '.partial')
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    temporary.replace(path)


def main():
    cpu, wall = time.process_time(), time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Fresh JSON artifact')
    args = parser.parse_args()
    output, model_path = args.output.resolve(), args.model.resolve(strict=True)
    require(not output.exists() and not output.is_symlink() and
            not output.with_suffix(output.suffix + '.partial').exists(), 'Fresh output and partial paths required')
    require(model_path.is_file() and model_path.stat().st_size <= MAX_MODEL_BYTES, 'Bounded model JSON required')
    require(not output.is_relative_to(Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')),
            'Diagnostic output must be outside the project checkout')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        stream.write('{}\n')
    raw = model_path.read_bytes()
    result = dict(schema='qhat-singleton-pauli-pair-count-v1', status='initializing',
        started_UTC=datetime.now(timezone.utc).isoformat(),
        model=dict(path=str(model_path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()),
        counter=dict(path=str(Path(__file__).resolve()),
                     sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        expected_words=EXPECTED_WORDS, expected_unordered_pairs=EXPECTED_UNORDERED_PAIRS,
        processed_unordered_pairs=0, qualifying_pairs=0, same_flip_pairs=0,
        commuting_nonzero_flip_pairs=0, stored_pairs=0,
        predicate='x_i XOR x_j != 0 and parity(popcount(x_i AND z_j) + popcount(z_i AND x_j)) == 1',
        interpretation='For distinct nonzero real singleton Paulis, anticommutation makes the HF '
            'commutator column nonzero. A nonzero flip XOR excludes the HF destination. '
            'This counts support entries only; it computes no objective or physical error.',
        source_imports=False, project_source_read=False, physical_errors_read=False,
        constructor_cap_changed=False)
    update(output, result, cpu, wall)
    try:
        model = json.loads(raw, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
        require(type(model['n']) is int and model['n'] == 12, 'Expected 12 spin modes')
        require(isinstance(model['denominator'], str) and int(model['denominator']) > 0,
                'Positive exact denominator required')
        coefficients = model['coefficients']
        require(isinstance(coefficients, list) and len(coefficients) == EXPECTED_WORDS,
                'Frozen CH4 Pauli-word inventory differs')
        masks, seen = [], set()
        for word, coefficient in coefficients:
            require(isinstance(word, list) and 0 < len(word) <= 12 and
                    isinstance(coefficient, str) and int(coefficient) != 0 and
                    str(int(coefficient)) == coefficient, 'Nonzero exact real coefficient required')
            x = z = ys = 0
            sites = []
            for q, letter in word:
                require(type(q) is int and 0 <= q < 12 and letter in ('X', 'Y', 'Z'), 'Pauli support bound')
                sites.append(q)
                if letter in ('X', 'Y'):
                    x |= 1 << q
                if letter in ('Z', 'Y'):
                    z |= 1 << q
                ys += letter == 'Y'
            require(sites == sorted(set(sites)) and ys % 2 == 0, 'Unique ordered support and even Y required')
            identity = tuple((q, letter) for q, letter in word)
            require(identity not in seen, 'Duplicate Pauli word')
            seen.add(identity)
            masks.append((x, z))
        count = len(masks)
        require(count * (count - 1) // 2 == EXPECTED_UNORDERED_PAIRS, 'Unordered-pair inventory differs')
        # Discard model payloads; the loop retains only O(number of words) masks.
        del model, coefficients, raw, seen, identity, word
        result.update(status='counting', words=count, n=12, real_even_Y_words_validated=True)
        update(output, result, cpu, wall)
        qualifying = same_flip = commuting = processed = 0
        for i, (x, z) in enumerate(masks):
            for j in range(i + 1, count):
                other_x, other_z = masks[j]
                if x == other_x:
                    same_flip += 1
                elif ((x & other_z).bit_count() + (z & other_x).bit_count()) & 1:
                    qualifying += 1
                else:
                    commuting += 1
            processed += count - i - 1
            if i % 128 == 0 or i == count - 1:
                result.update(processed_unordered_pairs=processed, qualifying_pairs=qualifying,
                              same_flip_pairs=same_flip, commuting_nonzero_flip_pairs=commuting,
                              completed_outer_indices=i + 1)
                update(output, result, cpu, wall)
        require(processed == qualifying + same_flip + commuting == EXPECTED_UNORDERED_PAIRS,
                'Complete disjoint pair partition failed')
        result.update(status='passed', qualifying_pairs=qualifying,
                      processed_unordered_pairs=processed, exact_integer_support_count=True)
        update(output, result, cpu, wall)
        print(json.dumps(dict(status='passed', words=count, unordered_pairs=processed,
                              qualifying_pairs=qualifying, CPU_seconds=result['CPU_seconds'],
                              wall_seconds=result['wall_seconds'], max_RSS_KiB=result['max_RSS_KiB']),
                         sort_keys=True), flush=True)
    except BaseException as error:
        result.update(status='failed', exception=dict(type=type(error).__name__, message=str(error)))
        update(output, result, cpu, wall)
        raise


if __name__ == '__main__':
    main()
