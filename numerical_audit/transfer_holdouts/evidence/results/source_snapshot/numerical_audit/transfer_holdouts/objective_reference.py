"""Independent literal coherent objectives and exact legacy-placement audit.

No producer feature, insertion, placement, descriptor or cubic-polynomial helper
is imported. Literal sparse matrices supply both nested commutators and the
integer occupation-column actions used to audit the legacy fixed-T1 placement.
"""
from functools import lru_cache
import math
from pathlib import Path
import numpy as np
from scipy import sparse
import circuit_reference

if Path(circuit_reference.__file__).resolve() != Path(__file__).resolve().parent / 'circuit_reference.py':
    raise ValueError('Wrong literal reference source')
require = circuit_reference.require


class Reference(circuit_reference.Reference):
    def __init__(self, model, initial):
        super().__init__(model, initial)
        self.blocks, self.raw_columns, self.masks = [], [], []
        for raw in model['blocks']:
            matrix = sparse.csr_matrix(self.H.shape, dtype=complex)
            columns = [{} for _ in range(self.dimension)]
            for word, coefficient in raw:
                operator = self.pool[tuple((q, p) for q, p in word)]
                matrix += (int(coefficient) / self.den) * operator
                exact = operator.tocoo()
                for row, column, phase in zip(exact.row, exact.col, exact.data):
                    require(phase.imag == 0 and phase.real in (-1, 1), 'Real-unit exact column domain')
                    values = columns[int(column)]
                    values[int(row)] = values.get(int(row), 0) + int(coefficient) * int(phase.real)
            columns = [{row: value for row, value in values.items() if value} for values in columns]
            masks = {row ^ column for column, values in enumerate(columns) for row in values}
            require(len(masks) == 1, 'Declared equal-flip complete unit required')
            self.blocks.append(matrix)
            self.raw_columns.append(columns)
            self.masks.append(next(iter(masks)))
        self.lengths = [len(raw) for raw in model['blocks']]
        self.paulis = [(coefficient / self.den) * self.pool[word] for word, coefficient in sorted(self.words.items())]
        self.support = np.flatnonzero(self.initial != 0)
        require(0 < len(self.support) <= 9, 'Bounded support before column allocation')
        self.columns = np.eye(self.dimension)[:, self.support]
        self.diagonal = self.H.diagonal().real

    @lru_cache(maxsize=4)
    def kernel(self, T):
        require(math.isfinite(T) and T > 0, 'Positive finite reference time required')
        result = np.empty((self.dimension, len(self.support)), complex)
        for u in range(self.dimension):
            for column, v in enumerate(self.support):
                x = T * float(self.diagonal[u] - self.diagonal[v]) / 2
                result[u, column] = T * (math.sin(x) / x if x else 1.) * complex(math.cos(x), math.sin(x))
        return result

    def projected_value(self, columns, T):
        action = (self.kernel(T) * columns) @ self.initial[self.support]
        action -= self.initial * np.vdot(self.initial, action)
        value = float(np.vdot(action, action).real)
        require(math.isfinite(value) and value >= 0, 'Invalid independently computed coherent score')
        return value

    def cubic_columns(self, order):
        require(all(type(i) is int for i in order) and sorted(order) == list(range(len(self.blocks))),
                'Independent complete-unit permutation required')
        rest = sparse.csr_matrix(self.H.shape, dtype=complex)
        out = np.zeros_like(self.columns, dtype=complex)
        for index in reversed(order):
            outer = self.blocks[index]
            av, bv = outer @ self.columns, rest @ self.columns
            aab = outer @ (outer @ bv) - 2 * outer @ (rest @ av) + rest @ (outer @ av)
            bba = rest @ (rest @ av) - 2 * rest @ (outer @ bv) + outer @ (rest @ bv)
            out += -aab / 24 + bba / 12
            rest += outer
        require(float(np.max(abs(out.imag))) < 2e-12, 'Real cubic coefficient became complex')
        return out.real

    def cubic_score(self, order, T, rotation_cap):
        require(type(rotation_cap) is int and rotation_cap > 0, 'Positive rotation cap')
        columns = self.cubic_columns(order)
        value = self.projected_value(columns, T)
        lengths = [self.lengths[i] for i in order]
        U = sum(lengths)
        steps = (1 if rotation_cap >= U else 0) if len(lengths) == 1 else max(
            0, (rotation_cap - lengths[0]) // (2 * U - lengths[0] - lengths[-1]))
        require(steps >= 1, 'Independent symmetric affordability failure')
        return dict(JD2=value, score=value / steps**4, steps=steps,
                    column_sources=self.support.tolist(), cubic_columns=columns.tolist())

    def first_score(self, order, T, rotation_cap, partition):
        require(partition in ('conserving', 'pauli') and type(rotation_cap) is int and rotation_cap > 0,
                'Declared first-order partition/cap')
        matrices = self.blocks if partition == 'conserving' else self.paulis
        require(all(type(i) is int for i in order) and sorted(order) == list(range(len(matrices))),
                'Independent first-order permutation required')
        prefix = sparse.csr_matrix(self.H.shape, dtype=complex)
        columns = np.zeros_like(self.columns, dtype=complex)
        for index in order:
            operator = matrices[index]
            columns += -.5j * (prefix @ (operator @ self.columns) - operator @ (prefix @ self.columns))
            prefix += operator
        value = self.projected_value(columns, T)
        steps = rotation_cap // len(self.model['coefficients'])
        require(steps >= 1, 'Independent first-order affordability failure')
        return dict(JD1=value, score=value / steps**2, steps=steps)

    def _raw_action(self, unit, vector):
        result = {}
        for column, amplitude in vector.items():
            for row, value in self.raw_columns[unit][column].items():
                result[row] = result.get(row, 0) + value * amplitude
        return result

    def _raw_commutator_column(self, left, right):
        first = self._raw_action(left, self.raw_columns[right][self.v])
        second = self._raw_action(right, self.raw_columns[left][self.v])
        for row, value in second.items():
            first[row] = first.get(row, 0) - value
        return {row: value for row, value in first.items() if value and row != self.v}

    def _raw_bch_column(self, order):
        total = {}
        for at, left in enumerate(order):
            for right in order[at + 1:]:
                for row, value in self._raw_commutator_column(left, right).items():
                    total[row] = total.get(row, 0) + value
        return total

    def current_placement(self):
        """Re-enumerate every slot using literal integer BCH columns.

        Supplied binary64 sinc weights follow the pinned legacy policy; exact
        integers reproduce its comparisons and ties. Its incremental updater
        and descriptor formulas are not reused.
        """
        raw_blocks = [{tuple((q, p) for q, p in word): int(value) for word, value in block}
                      for block in self.model['blocks']]
        diagonal = [i for i, mask in enumerate(self.masks) if mask == 0]
        require(all(len(raw_blocks[i]) == 1 for i in diagonal), 'Legacy one-word diagonal unit domain')
        diagonal.sort(key=lambda i: (-abs(next(iter(raw_blocks[i].values()))), next(iter(raw_blocks[i]))))
        off = sorted((i for i, mask in enumerate(self.masks) if mask != 0),
                     key=lambda i: (-sum(abs(c) for c in raw_blocks[i].values()), self.masks[i]))
        destinations = {self.v ^ self.masks[i] for i in off}
        for at, left in enumerate(off):
            for right in off[at + 1:]:
                destinations.update(self._raw_commutator_column(left, right))
        energies = [sum(self.raw_columns[i][v].get(v, 0) for i in diagonal) for v in range(self.dimension)]
        rho = {row: float(np.sinc((energies[row] - energies[self.v]) / (2 * np.pi * self.den))**2)
               for row in destinations}
        ratios = {row: value.as_integer_ratio() for row, value in rho.items()}
        denominator = max((b for a, b in ratios.values()), default=1)
        require(all(denominator % b == 0 for a, b in ratios.values()), 'Binary64 weight denominator differs')
        weights = {row: a * (denominator // b) for row, (a, b) in ratios.items()}
        def order_at(slots):
            order = []
            for slot in range(len(off) + 1):
                order.extend(i for i, position in zip(diagonal, slots) if position == slot)
                if slot < len(off):
                    order.append(off[slot])
            return order
        def score(slots):
            return sum(weights[row] * value * value for row, value in self._raw_bch_column(order_at(slots)).items()
                       if value)
        initial = [score([slot] * len(diagonal)) for slot in range(len(off) + 1)]
        start = min(range(len(initial)), key=lambda slot: (initial[slot], slot))
        slots, history = [start] * len(diagonal), []
        for index in range(len(diagonal)):
            old, before = slots[index], score(slots)
            values = []
            for slot in range(len(off) + 1):
                trial = list(slots)
                trial[index] = slot
                values.append(score(trial))
            selected = min(range(len(values)), key=lambda slot: (values[slot], slot))
            history.append(dict(index=index, old_slot=old, new_slot=selected,
                                old_score=str(before), new_score=str(values[selected]),
                                slot_scores=[str(value) for value in values]))
            require(values[selected] <= before, 'Legacy placement increased its exact surrogate')
            slots[index] = selected
        order = order_at(slots)
        amplitude = self._raw_bch_column(order)
        physical_score = .25 * sum(rho[row] * (value / (self.den * self.den))**2
                                   for row, value in amplitude.items() if value)
        return dict(order=order, score=physical_score, start_slot=start, slots=slots, history=history,
                    initial_slot_scores=[str(value) for value in initial], exact_weight_denominator=str(denominator),
                    exact_score_numerator=str(score(slots)), evaluations=(len(diagonal) + 1) * (len(off) + 1),
                    diagonal_priority=diagonal, off_priority=off)
