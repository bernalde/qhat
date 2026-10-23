"""Sparse occupation-response features and bounded, deterministic order searches.

Mode q is bit 1<<q. Coefficients are integer numerators with one positive
denominator. Search never calls an exact-Hamiltonian propagator.
"""
from __future__ import annotations

import math
import time
import numpy as np


def require(ok, message):
    if not ok:
        raise ValueError(message)


def masks(word):
    x = z = ny = 0
    for q, p in word:
        x |= (1 << q) if p in 'XY' else 0
        z |= (1 << q) if p in 'ZY' else 0
        ny += p == 'Y'
    return x, z, ny


def descriptor(block):
    rows = [(masks(w), int(c)) for w, c in block.items()]
    require(bool(rows), 'Empty unit')
    x = rows[0][0][0]
    require(all(m[0] == x and m[2] % 2 == 0 for m, c in rows),
            'A real equal-flip unit is required')
    return x, [(m[1], m[2] // 2, c) for m, c in rows]


def value(desc, v):
    return sum(c * (-1 if ((z & v).bit_count() + y) % 2 else 1)
               for z, y, c in desc[1])


def kernel(gap, T):
    return T * np.exp(.5j * T * gap) * np.sinc(T * gap / (2 * np.pi))


class Reference:
    """D or D plus one chosen determinant-to-determinant coupling.

    Only the supplied sparse C columns and their destination energies are
    evaluated. The corrected reference diagonalizes a 2x2 matrix, never H.
    """
    def __init__(self, diagonal, v, T, denominator, partner=None, coupling=0):
        require(T > 0 and denominator > 0, 'Invalid time/denominator')
        self.diagonal = descriptor(diagonal)
        self.v, self.T, self.den = v, T, denominator
        self.inputs = [v] if partner is None or coupling == 0 else [v, partner]
        self.energies = {}
        self.coupling = coupling / denominator
        if len(self.inputs) == 2:
            a, b = (self.energy(w) for w in self.inputs)
            self.eigenvalues, self.U = np.linalg.eigh([[a, self.coupling],
                                                       [self.coupling, b]])

    def energy(self, w):
        if w not in self.energies:
            self.energies[w] = value(self.diagonal, w) / self.den
        return self.energies[w]

    def response(self, columns):
        """Return (I-|v><v|) integral exp(itA) C exp(-itA)|v> dt."""
        v, T = self.v, self.T
        if len(self.inputs) == 1:
            return {w: kernel(self.energy(w) - self.energy(v), T) * c
                    for w, c in columns[0].items() if w != v and c}
        U, e = self.U, self.eigenvalues
        destinations = set(columns[0]) | set(columns[1])
        result = {}
        for w in destinations - set(self.inputs):
            row = np.array([col.get(w, 0j) for col in columns]) @ U
            result[w] = np.sum(kernel(self.energy(w) - e, T) * row * U[0])
        inside = np.array([[col.get(w, 0j) for col in columns]
                           for w in self.inputs])
        transformed = U.T @ inside @ U
        transported = kernel(e[:, None] - e[None, :], T) * transformed
        output = U @ (transported @ U[0])
        result[self.inputs[1]] = output[1]
        return {w: c for w, c in result.items() if c != 0}


def choose_partner(blocks, v):
    choices = []
    for block in blocks:
        desc = descriptor(block)
        if desc[0]:
            f = value(desc, v)
            if f:
                choices.append((-abs(f), desc[0], v ^ desc[0], f))
    if not choices:
        return None, 0
    _, _, partner, coupling = min(choices)
    return partner, coupling


class Features:
    """Linear response of every ordered unit pair; score is its squared norm."""
    def __init__(self, blocks, reference):
        began = time.monotonic()
        require(0<len(blocks)<=512,'Response unit limit is 512')
        desc = [descriptor(b) for b in blocks]
        maps, destinations = {}, set()
        den2 = reference.den ** 2
        for i in range(len(blocks)):
            for j in range(i + 1, len(blocks)):
                cols = []
                x, y = desc[i][0], desc[j][0]
                for v in reference.inputs:
                    b = (value(desc[j], v) * value(desc[i], v ^ y)
                         - value(desc[i], v) * value(desc[j], v ^ x))
                    cols.append({v ^ x ^ y: -.5j * (b / den2)} if b else {})
                response = reference.response(cols)
                if response:
                    maps[i, j] = response
                    destinations.update(response)
        self.destinations = sorted(destinations)
        columns = {w: k for k, w in enumerate(self.destinations)}
        self.pairs = {}
        for pair, response in maps.items():
            ids = np.array([columns[w] for w in response], dtype=np.int32)
            self.pairs[pair] = (ids, np.array(list(response.values()), dtype=complex))
        self.m = len(blocks)
        self.setup_seconds = time.monotonic() - began
        self.reference = reference

    def amplitude(self, order):
        require(sorted(order) == list(range(self.m)), 'Invalid unit permutation')
        pos = {i: k for k, i in enumerate(order)}
        a = np.zeros(len(self.destinations), dtype=complex)
        for (i, j), (ids, b) in self.pairs.items():
            a[ids] += b if pos[i] < pos[j] else -b
        return a

    @staticmethod
    def score(a):
        return float(np.vdot(a, a).real)

    def crossing(self, a, left, right):
        pair = min(left, right), max(left, right)
        if pair in self.pairs:
            ids, b = self.pairs[pair]
            a[ids] += -2 * b if left < right else 2 * b


def best_slot(order, moving, features):
    """Scan every insertion slot, recomputing each nonnegative score."""
    remaining = [i for i in order if i != moving]
    first = [moving] + remaining
    a = features.amplitude(first)
    scores = [features.score(a)]
    for crossed in remaining:
        features.crossing(a, moving, crossed)
        scores.append(features.score(a))
    slot = min(range(len(scores)), key=lambda k: (scores[k], k))
    return remaining[:slot] + [moving] + remaining[slot:], scores, slot


def placement(blocks, diagonal_ids, off_ids, features):
    began = time.monotonic()
    best = None
    evaluations = 0
    for slot in range(len(off_ids) + 1):
        order = off_ids[:slot] + diagonal_ids + off_ids[slot:]
        score = features.score(features.amplitude(order))
        evaluations += 1
        if best is None or score < best[0]:
            best = score, order
    initial_score, order = best
    initial_slot = sum(i in off_ids for i in order[:order.index(diagonal_ids[0])]) if diagonal_ids else 0
    slots = {i:initial_slot for i in diagonal_ids}
    def assemble():
        out=[]
        for slot in range(len(off_ids)+1):
            out.extend(i for i in diagonal_ids if slots[i]==slot)
            if slot<len(off_ids):out.append(off_ids[slot])
        return out
    history = []
    for moving in diagonal_ids:
        before = features.score(features.amplitude(order))
        a = features.amplitude(order)
        old = slots[moving]
        for j in reversed(range(old)):
            features.crossing(a,off_ids[j],moving)
        scores=[features.score(a)]
        for crossed in off_ids:
            features.crossing(a,moving,crossed)
            scores.append(features.score(a))
        slot=min(range(len(scores)),key=lambda k:(scores[k],k))
        slots[moving]=slot
        candidate=assemble()
        evaluations += len(scores)
        after = features.score(features.amplitude(candidate))
        # Cancellation can put incremental scores a few ulps apart. Always
        # retain the original position if the independently rescored candidate
        # would increase the objective.
        if after <= before:
            order = candidate
        else:
            slots[moving]=old
            after = before
        history.append(dict(unit=moving,slot=slot,before=before,after=after))
    return dict(order=order,score=features.score(features.amplitude(order)),
                initial_score=initial_score,evaluations=evaluations,history=history,
                seconds=time.monotonic()-began)


def insertion(starts, features, budget, max_moves=10):
    began = time.monotonic()
    require(bool(starts) and 0 <= budget <= 2_000_000 and 0 <= max_moves <= 10,
            'Invalid or oversized search allowance')
    candidates = []
    used = 0
    for start_id, start in enumerate(starts):
        cap = budget // len(starts) + (start_id < budget % len(starts))
        order = list(start)
        a = features.amplitude(order)
        score = features.score(a)
        count = 0
        accepted = []
        for sweep in range(max_moves):
            best, move = score, None
            exhausted = False
            for source in range(len(order)):
                for direction in (-1, 1):
                    trial = a.copy()
                    at = source
                    while 0 <= at + direction < len(order):
                        if count >= cap:
                            exhausted = True
                            break
                        crossed = order[at+direction]
                        left, right = ((crossed,order[source]) if direction < 0
                                       else (order[source],crossed))
                        features.crossing(trial,left,right)
                        at += direction
                        count += 1
                        candidate = features.score(trial)
                        if candidate < best-max(1e-30,1e-12*best):
                            best, move = candidate,(source,at)
                    if exhausted:
                        break
                if exhausted:
                    break
            if move is None:
                break
            candidate_order = list(order)
            item = candidate_order.pop(move[0]);candidate_order.insert(move[1],item)
            actual_a = features.amplitude(candidate_order)
            actual_score = features.score(actual_a)
            if actual_score >= score-max(1e-30,1e-12*score):
                break
            order,a,score = candidate_order,actual_a,actual_score
            accepted.append(dict(move=move,score=score,evaluations=count))
            if exhausted:
                break
        used += count
        candidates.append(dict(order=order,score=score,evaluations=count,accepted=accepted))
    selected = min(candidates,key=lambda x:(x['score'],tuple(x['order'])))
    return dict(order=selected['order'],score=selected['score'],evaluations=used,
                budget=budget,starts=candidates,seconds=time.monotonic()-began)


def pulses(unit_lengths, steps, symmetric=False):
    """Exact chronological word indices and fractions of total physical time."""
    require(steps >= 1 and bool(unit_lengths) and min(unit_lengths)>0,
            'Invalid formula dimensions')
    offsets = np.cumsum([0]+list(unit_lengths))
    if not symmetric:
        ids = np.tile(np.arange(offsets[-1],dtype=np.int32),steps)
        return ids,np.full(len(ids),1/steps)
    if len(unit_lengths)==1:
        return np.arange(offsets[-1],dtype=np.int32),np.ones(offsets[-1])
    sequence = []
    def append(unit, duration):
        if sequence and sequence[-1][0]==unit:
            sequence[-1] = unit,sequence[-1][1]+duration
        else:
            sequence.append((unit,duration))
    for _ in range(steps):
        for unit in range(len(unit_lengths)-1):append(unit,1)
        append(len(unit_lengths)-1,2)
        for unit in reversed(range(len(unit_lengths)-1)):append(unit,1)
    ids,fractions=[],[]
    for unit,duration in sequence:
        indices=range(offsets[unit],offsets[unit+1])
        ids.extend(indices);fractions.extend([duration/(2*steps)]*len(indices))
    return np.array(ids,dtype=np.int32),np.array(fractions)


def feasible_steps(lengths, budget, symmetric=False):
    require(bool(lengths) and min(lengths)>0 and budget>=0,'Invalid budget/units')
    U=sum(lengths)
    if not symmetric:return budget//U
    if len(lengths)==1:return 1 if budget>=U else 0
    return max(0,(budget-lengths[0])//(2*U-lengths[0]-lengths[-1]))
