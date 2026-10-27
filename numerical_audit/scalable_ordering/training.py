"""Six bounded constructors; no validation reference or physical error is read."""
import math
import time
import numpy as np
import common as k


class Deadline(RuntimeError):
    pass


def check(deadline):
    if time.process_time() >= deadline:
        raise Deadline('constructor CPU deadline')


def descriptor(block):
    rows = []
    flip = None
    for word, coefficient in block.items():
        x = z = ys = 0
        for q, letter in word:
            if letter in ('X', 'Y'):
                x |= 1 << q
            if letter in ('Y', 'Z'):
                z |= 1 << q
            ys += letter == 'Y'
        k.require(ys % 2 == 0 and (flip is None or flip == x), 'Real equal-flip complete unit')
        flip = x
        rows.append((z, ys // 2, coefficient))
    k.require(rows, 'Nonempty complete unit')
    return flip, rows


def value(desc, occupation):
    return sum(c * (-1 if ((z & occupation).bit_count() + y) % 2 else 1) for z, y, c in desc[1])


def current(model, deadline=float('inf')):
    """The retained exact supplied-weight diagonal insertion rule at T=1."""
    raw = [k.decode(block) for block in model['blocks']]
    desc = [descriptor(block) for block in raw]
    diagonal = sorted((i for i, d in enumerate(desc) if d[0] == 0),
                      key=lambda i: (-abs(next(iter(raw[i].values()))), next(iter(raw[i]))))
    off = sorted((i for i, d in enumerate(desc) if d[0]),
                 key=lambda i: (-sum(abs(c) for c in raw[i].values()), desc[i][0]))
    if not diagonal or not off:
        return diagonal + off
    hf, denominator = (1 << model['electrons']) - 1, int(model['denominator'])
    masks = [desc[i][0] for i in off]
    amplitude = {}
    for at, i in enumerate(off):
        check(deadline)
        for j in off[at + 1:]:
            x, y = desc[i][0], desc[j][0]
            b = value(desc[j], hf) * value(desc[i], hf ^ y) - value(desc[i], hf) * value(desc[j], hf ^ x)
            destination = x ^ y
            if destination and b:
                amplitude[destination] = amplitude.get(destination, 0) + b
    rows = [[value(desc[j], hf) * (value(desc[i], hf ^ desc[j][0]) - value(desc[i], hf))
             for j in off] for i in diagonal]
    D = descriptor({word: c for i in diagonal for word, c in raw[i].items()})
    ratios = {d: float(np.sinc((value(D, hf ^ d) - value(D, hf)) /
                             (2 * np.pi * denominator))**2).as_integer_ratio()
              for d in set(amplitude) | set(masks)}
    common = max((b for a, b in ratios.values()), default=1)
    weights = {d: a * (common // b) for d, (a, b) in ratios.items()}
    score = lambda a: sum(weights[d] * x * x for d, x in a.items() if x)
    seed_scores = []
    for slot in range(len(off) + 1):
        a = dict(amplitude)
        for row in rows:
            for j, (d, b) in enumerate(zip(masks, row)):
                a[d] = a.get(d, 0) + (b if j >= slot else -b)
        seed_scores.append(score(a))
    start = min(range(len(seed_scores)), key=lambda i: (seed_scores[i], i))
    slots = [start] * len(diagonal)
    a = dict(amplitude)
    for row in rows:
        for j, (d, b) in enumerate(zip(masks, row)):
            a[d] = a.get(d, 0) + (b if j >= start else -b)
    for index, row in enumerate(rows):
        check(deadline)
        trial = dict(a)
        for j in range(slots[index]):
            d = masks[j]
            trial[d] = trial.get(d, 0) + 2 * row[j]
        scores = [score(trial)]
        for d, b in zip(masks, row):
            trial[d] = trial.get(d, 0) - 2 * b
            scores.append(score(trial))
        chosen = min(range(len(scores)), key=lambda j: scores[j])
        for j, (d, b) in enumerate(zip(masks, row)):
            a[d] = a.get(d, 0) + ((1 if j >= chosen else -1) - (1 if j >= slots[index] else -1)) * b
        slots[index] = chosen
    order = []
    for slot in range(len(off) + 1):
        order.extend(i for i, position in zip(diagonal, slots) if position == slot)
        if slot < len(off):
            order.append(off[slot])
    return order


class FirstResponse:
    """HF response from sparse pair columns; no 2**n vectors or propagation."""
    def __init__(self, raw, denominator, hf, diagonal, deadline):
        self.m = len(raw)
        k.require(self.m * (self.m - 1) // 2 <= 500_000, 'First-response pair storage bound')
        desc = [descriptor(b) for b in raw]
        D = descriptor(diagonal or {(): 0})
        self.pairs, self.weights = [], {}
        for i in range(self.m):
            check(deadline)
            for j in range(i + 1, self.m):
                x, y = desc[i][0], desc[j][0]
                destination = hf ^ x ^ y
                if destination == hf:
                    continue
                b = value(desc[j], hf) * value(desc[i], hf ^ y) - value(desc[i], hf) * value(desc[j], hf ^ x)
                if not b:
                    continue
                self.pairs.append((i, j, destination, 0.5 * (b / denominator**2)))
                if destination not in self.weights:
                    gap = (value(D, destination) - value(D, hf)) / denominator
                    xphase = k.TRAIN_TIME * gap / 2
                    self.weights[destination] = k.TRAIN_TIME**2 * (math.sin(xphase) / xphase if xphase else 1.)**2
        self.lookup = {i * self.m + j: (d, c) for i, j, d, c in self.pairs}

    def amplitude(self, order, deadline):
        check(deadline)
        positions = {unit: at for at, unit in enumerate(order)}
        amplitude = {}
        for at, (i, j, destination, coefficient) in enumerate(self.pairs):
            if at % 1024 == 0:
                check(deadline)
            amplitude[destination] = amplitude.get(destination, 0.) + (coefficient if positions[i] < positions[j] else -coefficient)
        return amplitude

    def score(self, order, deadline):
        amplitude = self.amplitude(order, deadline)
        result = sum(self.weights[d] * a * a for d, a in amplitude.items())
        k.require(math.isfinite(result) and result >= 0, 'Finite nonnegative first response')
        return result

    def crossing(self, amplitude, score, left, right):
        pair = self.lookup.get(min(left, right) * self.m + max(left, right))
        if pair is None:
            return score
        destination, coefficient = pair
        old = amplitude.get(destination, 0.)
        new = old + (-2 if left < right else 2) * coefficient
        amplitude[destination] = new
        return score + self.weights[destination] * (new * new - old * old)


def first_search(model, arm, initial, deadline):
    raw = [k.decode(block) for block in model['blocks']]
    diagonal = {w: c for b in raw if descriptor(b)[0] == 0 for w, c in b.items()}
    if arm == 'pauli':
        index = {w: i for i, w in enumerate(sorted(k.decode(model['coefficients'])))}
        initial = [index[w] for i in initial for w in sorted(raw[i])]
        raw = [{w: c} for w, c in sorted(k.decode(model['coefficients']).items())]
    seeds = [tuple(initial), tuple(reversed(initial)), tuple(range(len(raw)))]
    best, best_score, evaluations, moves = seeds[0], None, 0, 0
    context, selected = None, None
    try:
        context = FirstResponse(raw, int(model['denominator']), (1 << model['electrons']) - 1, diagonal, deadline)
        scored = []
        for seed in dict.fromkeys(seeds):
            score = context.score(seed, deadline)
            evaluations += 1
            scored.append((seed, score))
            if best_score is None or (score, seed) < (best_score, best):
                best, best_score = seed, score
        for seed, initial_score in scored:
            current_order, score = seed, initial_score
            for _ in range(5):
                selected, selected_score = current_order, score
                base_amplitude = context.amplitude(current_order, deadline)
                for source in range(len(raw)):
                    for direction in (-1, 1):
                        candidate = list(current_order)
                        amplitude = dict(base_amplitude)
                        value_, position = score, source
                        while 0 <= position + direction < len(raw):
                            if evaluations >= 3000:
                                raise Deadline('candidate evaluation limit')
                            check(deadline - .02)
                            left = position - 1 if direction < 0 else position
                            value_ = context.crossing(amplitude, value_, candidate[left], candidate[left + 1])
                            candidate[left], candidate[left + 1] = candidate[left + 1], candidate[left]
                            position += direction
                            evaluations += 1
                            choice = tuple(candidate)
                            if (value_, choice) < (selected_score, selected):
                                selected, selected_score = choice, value_
                if selected == current_order:
                    break
                actual = context.score(selected, deadline)
                if actual >= score - max(1e-30, 1e-12 * abs(score)):
                    break
                current_order, score = selected, actual
                moves += 1
                if (score, current_order) < (best_score, best):
                    best, best_score = current_order, score
        status = 'complete'
    except Deadline:
        status = 'evaluation_limit' if evaluations >= 3000 else 'cpu_limit'
        if context is not None and selected is not None:
            try:
                actual = context.score(selected, deadline)
            except Deadline:
                pass
            else:
                if best_score is None or (actual, selected) < (best_score, best):
                    best, best_score = selected, actual
                    moves += 1
    return dict(order=list(best), score=best_score, evaluations=evaluations,
                accepted_moves=moves, status=status, partition='pauli' if arm == 'pauli' else 'conserving')


def construct(model_path, digest, arm):
    cpu, wall = time.process_time(), time.monotonic()
    k.require(arm in k.ARMS and k.sha(model_path) == digest, 'Unknown arm or input changed')
    model = k.read(model_path)
    units = k.blocks(model)
    canonical = list(range(len(units)))
    magnitude = sorted(canonical, key=lambda i: (-sum(abs(c) for c in units[i].values()), i))
    deadline = cpu + k.TRAIN_CPU_SECONDS
    result = dict(order=canonical, score=None, evaluations=0, accepted_moves=0,
                  status='complete', partition='conserving')
    if arm == 'pauli':
        result.update(partition='pauli', order=list(range(len(model['coefficients']))))
    try:
        if arm == 'canonical':
            result['order'] = canonical
        elif arm == 'magnitude':
            result['order'] = magnitude
        else:
            placed = current(model, deadline)
            if arm != 'pauli':
                result['order'] = placed
            if arm in ('first', 'pauli'):
                result = first_search(model, arm, placed, deadline)
            elif arm == 'cubic':
                remaining = deadline - time.process_time()
                if remaining > 0:
                    ordering = k.bootstrap()
                    trained = ordering.train_ordering(units, model['n'], (1 << model['electrons']) - 1,
                        time=k.TRAIN_TIME, rotation_budget=30 * len(model['coefficients']),
                        cpu_seconds=remaining, max_evaluations=3000, max_moves=5,
                        starts=[placed, list(reversed(placed)), canonical])
                    result.update(order=list(trained.order), score=None if trained.score is None else trained.score.value,
                        evaluations=trained.evaluations, accepted_moves=trained.accepted_moves,
                        status=trained.status, stats=trained.stats)
                else:
                    result['status'] = 'cpu_limit'
    except Deadline:
        result['status'] = 'cpu_limit'
    if result['partition'] == 'pauli':
        raw = [{word: c} for word, c in sorted(k.decode(model['coefficients']).items())]
    else:
        raw = [k.decode(block) for block in model['blocks']]
    order = result['order']
    k.require(sorted(order) == list(range(len(raw))), 'Constructor must return complete permutation')
    result.update(arm=arm, formula=2, time=k.TRAIN_TIME, model_sha256=digest,
        unit_lengths=[len(raw[i]) for i in order],
        schedule=[[[[q, p] for q, p in word], str(c)] for i in order for word, c in sorted(raw[i].items())],
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall,
        nominal_CPU_seconds=k.TRAIN_CPU_SECONDS if arm not in ('canonical', 'magnitude') else None,
        response_is_error_bound=False,
        timing='Caller through return; input authentication, setup, initialization, search and layout; imports/campaign startup separate')
    return result
