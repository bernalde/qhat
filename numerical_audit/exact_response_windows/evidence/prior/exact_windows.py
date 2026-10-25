"""Bounded exhaustive eight-unit windows; exact for supplied floating weights."""
import argparse
import itertools
import json
import math
from pathlib import Path
import time

from campaign import frozen, helpers, require, write


def window_solve(order, pairs, weights, denominator, offset, width):
    common, optimizer, _, _ = helpers()
    started = time.monotonic()
    window = order[offset:offset+width]
    require(len(window) == width, 'Incomplete exact window')
    selected = set(window)
    full = optimizer.amplitudes(order, pairs)
    inside = {(i, j): value for (i, j), value in pairs.items() if i in selected and j in selected}
    original = optimizer.amplitudes(window, inside)
    base = dict(full)
    for dest, value in original.items():
        base[dest] -= value
    ratios = {dest: weight.as_integer_ratio() for dest, weight in weights.items()}
    scale = max((den for _, den in ratios.values()), default=1)
    integer_weights = {dest: num*(scale//den) for dest, (num, den) in ratios.items()}
    affected = sorted({dest for dest, _ in inside.values()})
    constant = sum(integer_weights[dest]*value*value for dest, value in base.items() if dest not in affected)
    initial_integer = sum(integer_weights[dest]*value*value for dest, value in full.items())
    best = initial_integer
    winner = list(window)
    evaluated = 0
    ties = 0
    for permutation in itertools.permutations(window):
        require(time.monotonic()-started <= 180, 'Exact window exceeds 180-second limit')
        positions = {unit: i for i, unit in enumerate(permutation)}
        amplitudes = {dest: base.get(dest, 0) for dest in affected}
        for (i, j), (dest, value) in inside.items():
            amplitudes[dest] += value if positions[i] < positions[j] else -value
        score = constant + sum(integer_weights[dest]*value*value for dest, value in amplitudes.items())
        evaluated += 1
        if score < best:
            best, winner, ties = score, list(permutation), 1
        elif score == best:
            ties += 1
            if tuple(permutation) < tuple(winner):
                winner = list(permutation)
    require(evaluated == math.factorial(width), 'Incomplete enumeration')
    optimum_order = order[:offset] + winner + order[offset+width:]
    exact_recomputed = sum(integer_weights[dest]*value*value
                           for dest, value in optimizer.amplitudes(optimum_order, pairs).items())
    require(exact_recomputed == best, 'Window objective does not match complete-order reconstruction')
    divisor = 4*denominator**4*scale
    return dict(offset=offset, width=width, evaluated=evaluated, minimizers=ties,
                initial_J=initial_integer/divisor, optimum_J=best/divisor,
                incumbent_integer=str(best), lower_bound_integer=str(best),
                objective_denominator=str(divisor), original_window=window, selected_window=winner,
                elapsed=time.monotonic()-started,
                scope='Complete enumeration of this window only; exact rational values of supplied binary64 weights; outside order fixed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    record = frozen(args)
    common, optimizer, _, _ = helpers()
    results = []
    started = time.monotonic()
    for case in ('CH2O_s1.00', 'nh3_coeff1p00'):
        design = json.loads((args.output/'designs'/(case+'.json')).read_text())
        blocks = [common.block(raw) for raw in design['blocks']['unsplit']]
        flips = [common.masks(next(iter(block)))[0] for block in blocks]
        denominator = int(design['denominator'])
        pairs, weights = optimizer.pair_setup(blocks, flips, (1<<design['electrons'])-1, denominator)
        order = design['methods']['unsplit_JD']['unit_order']
        for offset in (0, (len(order)-8)//2):
            result = window_solve(order, pairs, weights, denominator, offset, 8)
            result['case'] = case
            results.append(result)
            print('EXACT WINDOW', case, offset, result['initial_J'], result['optimum_J'], flush=True)
    write(args.output/'exact_windows.json', dict(status='passed', windows=results,
          seconds=time.monotonic()-started, global_order_certificate=False,
          real_sinc_weight_certificate=False, finite_step_error_certificate=False))


if __name__ == '__main__':
    main()
