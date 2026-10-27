"""Independent literal states, primitive Clifford-frame replay and DAG costs.

No sparse trainer, producer descriptor, search, layout or compiler helper is
imported. Full-space allocation belongs only to this twelve-mode validator.
"""
from array import array
from collections import OrderedDict
from fractions import Fraction
import gzip
import json
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import common as k


def multiply(left, right):
    """Multiply i**phase X**x Z**z Pauli frames exactly."""
    x, z, phase = left
    xx, zz, other = right
    return x ^ xx, z ^ zz, (phase + other + 2 * (z & xx).bit_count()) % 4


class Reference:
    def __init__(self, model):
        self.n = model['n']
        k.require(type(self.n) is int and 1 <= self.n <= 12, 'Reference register before allocation')
        k.require(0 < len(model['coefficients']) <= k.MAX_TERMS, 'Reference term count')
        self.dimension = 1 << self.n
        self.model = model
        self.den = int(model['denominator'])
        k.require(self.den > 0, 'Reference denominator')
        self.hf = (1 << model['electrons']) - 1
        k.require(0 <= model['electrons'] <= self.n, 'Reference occupation')
        self.basis = np.arange(self.dimension, dtype=np.int64)
        self.initial = np.zeros(self.dimension, dtype=np.clongdouble)
        self.initial[self.hf] = 1
        self.words = {}
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.max_cache_bytes = 64 * 1024**2
        flips = set()
        for raw, coefficient in model['coefficients']:
            word = tuple((q, p) for q, p in raw)
            k.require(word and word not in self.words and len({q for q, p in word}) == len(word) and
                      all(type(q) is int and 0 <= q < self.n and p in ('X', 'Y', 'Z') for q, p in word),
                      'Independent Pauli support')
            coefficient = int(coefficient)
            self.words[word] = coefficient
            flips.add(sum(1 << q for q, p in word if p in ('X', 'Y')))
        # At most one H entry per column and distinct flip mask. Reserve room
        # for CSR addition/Hermiticity temporaries before building any matrix.
        k.require(len(flips) <= 512 and self.dimension * len(flips) <= 2_097_152,
                  'Reference sparse matrix storage bound before allocation')
        self.H = sparse.csr_matrix((self.dimension, self.dimension), dtype=np.complex128)
        for word, coefficient in self.words.items():
            rows, values = self.literal_action(word)
            self.H += sparse.csr_matrix(((coefficient / self.den) * values, (rows, self.basis)),
                                        shape=self.H.shape)
        self.H.eliminate_zeros()
        k.require(self.H.nnz <= self.dimension * min(1024, len(self.words)), 'Reference sparse H storage')
        difference = self.H - self.H.getH()
        difference.eliminate_zeros()
        k.require(difference.nnz == 0, 'Independent Hermiticity')
        flat = [(tuple((q, p) for q, p in word), int(c)) for block in model['blocks'] for word, c in block]
        k.require(len(flat) == len(self.words) and dict(flat) == self.words, 'Reference exact partition')
        self.low = {q: self.basis[(self.basis & (1 << q)) == 0] for q in range(self.n)}
        self.high = {q: self.low[q] ^ (1 << q) for q in range(self.n)}
        self.diagonal = np.zeros(self.dimension, dtype=np.longdouble)
        for word, c in self.words.items():
            if all(p == 'Z' for q, p in word):
                _, values = self.literal_action(word)
                self.diagonal += (np.longdouble(c) / np.longdouble(self.den)) * values.real

    def literal_action(self, word):
        key = tuple(word)
        if key in self.cache:
            result = self.cache.pop(key)
            self.cache[key] = result
            return result
        destination = self.basis.copy()
        values = np.ones(self.dimension, dtype=np.complex128)
        for q, letter in key:
            bit = (self.basis >> q) & 1
            if letter in ('Z', 'Y'):
                values *= 1 - 2 * bit
            if letter == 'Y':
                values *= 1j
            if letter in ('X', 'Y'):
                destination ^= 1 << q
        size = destination.nbytes + values.nbytes
        while self.cache and self.cache_bytes + size > self.max_cache_bytes:
            _, old = self.cache.popitem(last=False)
            self.cache_bytes -= old[0].nbytes + old[1].nbytes
        if size <= self.max_cache_bytes:
            self.cache[key] = (destination, values)
            self.cache_bytes += size
        return destination, values

    def apply_word(self, word, state):
        destination, values = self.literal_action(word)
        out = np.empty_like(state)
        out[destination] = values * state
        return out

    def exact_action(self, T):
        k.require(type(T) in (int, float) and math.isfinite(T) and T > 0, 'Reference time')
        matrix = (-1j * T) * self.H
        result = expm_multiply(matrix, np.asarray(self.initial, np.complex128), traceA=matrix.diagonal().sum())
        k.require(np.all(np.isfinite(result)) and abs(np.linalg.norm(result) - 1) < 2e-10,
                  'Reference exponential finite/norm check')
        return result

    def series_action(self, T):
        """Separate scaled 24-term Taylor action for numerical reference parity.

        The truncation is tiny at step norm <=1/8; floating arithmetic is still
        not an interval certificate. Bound the scaling count before work.
        """
        k.require(type(T) in (int, float) and math.isfinite(T) and T > 0, 'Series reference time')
        norm_bound = float(np.max(np.asarray(abs(self.H).sum(axis=1)).ravel())) * T
        steps = max(1, math.ceil(8 * norm_bound))
        k.require(steps <= 4096, 'Series reference scaling-work bound')
        state = self.initial.copy()
        for _ in range(steps):
            term, out = state.copy(), state.copy()
            for degree in range(1, 25):
                term = (-1j * T / steps / degree) * (self.H @ term)
                out += term
            state = out
        return state

    def ideal_formula(self, record, steps, T):
        k.require(type(steps) is int and 1 <= steps <= 128 and record['formula'] == 2, 'Symmetric step domain')
        schedule = [(tuple((q, p) for q, p in word), int(c)) for word, c in record['schedule']]
        lengths = record['unit_lengths']
        k.require(sum(lengths) == len(schedule) and all(c == self.words[word] for word, c in schedule),
                  'Reference schedule target')
        offsets = np.cumsum([0] + lengths)
        state = self.initial.copy()
        def unit(index, duration):
            nonlocal state
            for at in range(offsets[index], offsets[index + 1]):
                word, coefficient = schedule[at]
                theta = np.longdouble(coefficient) / np.longdouble(self.den) * duration
                state = np.cos(theta) * state - 1j * np.sin(theta) * self.apply_word(word, state)
        for _ in range(steps):
            for i in range(len(lengths) - 1):
                unit(i, np.longdouble(T) / (2 * steps))
            unit(len(lengths) - 1, np.longdouble(T) / steps)
            for i in reversed(range(len(lengths) - 1)):
                unit(i, np.longdouble(T) / (2 * steps))
        return state

    def apply_frame(self, pauli, state):
        x, z, phase = pauli
        destination = self.basis ^ x
        values = np.full(self.dimension, 1j**phase, dtype=np.complex128)
        remaining = z
        while remaining:
            bit = remaining & -remaining
            values *= 1 - 2 * ((self.basis & bit) != 0)
            remaining ^= bit
        out = np.empty_like(state)
        out[destination] = values * state
        return out

    def objective(self, record, T=1., rotation_multiple=30, *, kind='cubic'):
        """Literal nested-commutator oracle for each final frozen order.

        This separate full-space calculation does not call SparseResponse,
        Taylor-column producer code, descriptors or insertion helpers. The
        first-response baseline stores its raw squared response; the cubic
        scorer additionally uses the actual symmetric affordability factor.
        """
        k.require(kind in ('first', 'cubic') and T > 0 and math.isfinite(T), 'Objective oracle domain')
        raw = ([{w: c} for w, c in sorted(self.words.items())] if record['partition'] == 'pauli'
               else [dict((tuple((q, p) for q, p in word), int(c)) for word, c in b)
                     for b in self.model['blocks']])
        order = record['order']
        k.require(all(type(i) is int for i in order) and sorted(order) == list(range(len(raw))),
                  'Independent objective complete permutation')
        def matrix_for(block):
            matrix = sparse.csr_matrix(self.H.shape, dtype=np.complex128)
            for word, coefficient in block.items():
                rows, values = self.literal_action(word)
                matrix += sparse.csr_matrix(((coefficient / self.den) * values, (rows, self.basis)),
                                            shape=self.H.shape)
            return matrix
        column = np.zeros(self.dimension, dtype=np.clongdouble)
        suffix = sparse.csr_matrix(self.H.shape, dtype=np.complex128)
        for index in reversed(order):
            # Keep only the current unit and accumulated suffix. The complete
            # Hamiltonian's guarded flip count bounds both sparse supports.
            outer = matrix_for(raw[index])
            av, bv = outer @ self.initial, suffix @ self.initial
            if kind == 'first':
                column += -.5j * (outer @ bv - suffix @ av)
            else:
                aab = outer @ (outer @ bv) - 2 * (outer @ (suffix @ av)) + suffix @ (outer @ av)
                bba = suffix @ (suffix @ av) - 2 * (suffix @ (outer @ bv)) + outer @ (suffix @ bv)
                column += -aab / 24 + bba / 12
            suffix += outer
        phase = np.longdouble(T) * (self.diagonal - self.diagonal[self.hf]) / 2
        sinc = np.ones_like(phase)
        np.divide(np.sin(phase), phase, out=sinc, where=phase != 0)
        action = np.longdouble(T) * sinc * np.exp(1j * phase) * column
        action[self.hf] = 0
        response = float(np.vdot(action, action).real)
        lengths = [len(raw[i]) for i in order]
        if len(lengths) == 1:
            steps = 1
        else:
            budget = rotation_multiple * len(self.words)
            steps = (budget - lengths[0]) // (2 * sum(lengths) - lengths[0] - lengths[-1])
        k.require(steps >= 1 and math.isfinite(response) and response >= 0, 'Independent response value')
        result = dict(value=response if kind == 'first' else response / steps**4,
            raw_response=response, steps=steps, kind=kind,
            method='Literal sparse matrices and nested commutators; independent of trainer')
        del suffix
        return result

    def audit_trace(self, path, identity, compiler_version, *, literal=False):
        """Inspect every primitive; numerically apply only its actual Rz rotations.

        Maintain images C†XqC/C†ZqC, so left-appended Rz(theta) equals C times
        exp(-i theta C†ZqC/2). The final Clifford is HF's X mask up to phase.
        """
        x_images = [(1 << q, 0, 0) for q in range(self.n)]
        z_images = [(0, 1 << q, 0) for q in range(self.n)]
        state = np.zeros(self.dimension, dtype=np.clongdouble)
        state[0] = 1
        direct = state.copy() if literal else None
        previous = [-1] * self.n
        left, right, twoq = array('i'), array('i'), bytearray()
        one = two = rotations = 0
        with gzip.open(path, 'rt') as stream:
            header_line = stream.readline(16385)
            k.require(0 < len(header_line) <= 16384, 'Independent header line bound')
            header = json.loads(header_line)
            k.require(header['n'] == self.n and header['compiler'] == compiler_version and
                      header['identity'] == identity, 'Trace header identity')
            index = -1
            while True:
                line = stream.readline(16385)
                if not line:
                    break
                index += 1
                k.require(index < 2_000_000 and len(line) <= 16384, 'Independent gate/line bound')
                gate = json.loads(line)
                k.require(isinstance(gate, list) and len(gate) >= 2 and type(gate[1]) is int and
                          0 <= gate[1] < self.n, 'Primitive gate support')
                kind, q = gate[:2]
                if kind == 'cx':
                    k.require(len(gate) == 3 and type(gate[2]) is int and 0 <= gate[2] < self.n and
                              abs(q - gate[2]) == 1, 'Primitive line CX')
                    t = gate[2]
                    left.append(previous[q]); right.append(previous[t]); twoq.append(1)
                    previous[q] = previous[t] = index
                    two += 1
                    x_images[q] = multiply(x_images[q], x_images[t])
                    z_images[t] = multiply(z_images[q], z_images[t])
                    if literal:
                        destination = self.basis ^ (((self.basis >> q) & 1) << t)
                        direct = direct[destination]
                    continue
                k.require(kind in ('h', 's', 'sdg', 'rz'), 'Primitive gate kind')
                left.append(previous[q]); right.append(-1); twoq.append(0)
                previous[q] = index
                one += 1
                if kind == 'h':
                    k.require(len(gate) == 2, 'H payload')
                    x_images[q], z_images[q] = z_images[q], x_images[q]
                    if literal:
                        a, b = direct[self.low[q]].copy(), direct[self.high[q]].copy()
                        direct[self.low[q]], direct[self.high[q]] = (a + b) / np.sqrt(np.longdouble(2)), (a - b) / np.sqrt(np.longdouble(2))
                elif kind in ('s', 'sdg'):
                    k.require(len(gate) == 2, 'S payload')
                    x, z, phase = multiply(x_images[q], z_images[q])
                    x_images[q] = x, z, (phase + (3 if kind == 's' else 1)) % 4
                    if literal:
                        direct[self.high[q]] *= 1j if kind == 's' else -1j
                else:
                    k.require(len(gate) == 4 and isinstance(gate[2], str) and isinstance(gate[3], str) and
                              int(gate[3]) > 0, 'Rz exact payload')
                    angle = Fraction(int(gate[2]), int(gate[3]))
                    k.require(angle and str(angle.numerator) == gate[2] and str(angle.denominator) == gate[3],
                              'Canonical nonzero Rz angle')
                    theta = np.longdouble(angle.numerator) / np.longdouble(angle.denominator) / 2
                    state = np.cos(theta) * state - 1j * np.sin(theta) * self.apply_frame(z_images[q], state)
                    rotations += 1
                    if literal:
                        direct[self.low[q]] *= np.cos(theta) - 1j * np.sin(theta)
                        direct[self.high[q]] *= np.cos(theta) + 1j * np.sin(theta)
        for q in range(self.n):
            k.require(x_images[q] == (1 << q, 0, 0) and
                      z_images[q] == (0, 1 << q, 2 * ((self.hf >> q) & 1)), 'Final frame must be HF bit flip')
        state = state[self.basis ^ self.hf]
        if literal:
            overlap = np.vdot(direct, state)
            k.require(abs(abs(overlap) - 1) < 2e-11 and
                      np.linalg.norm(state - direct * overlap / abs(overlap)) < 2e-10,
                      'Frame replay differs from literal primitive replay')
        depth, cx_depth = array('i'), array('i')
        for index, (a, b, cx) in enumerate(zip(left, right, twoq)):
            k.require(a < index and b < index, 'Independent DAG acyclicity')
            depth.append(max(depth[a] if a >= 0 else 0, depth[b] if b >= 0 else 0) + 1)
            cx_depth.append(max(cx_depth[a] if a >= 0 else 0, cx_depth[b] if b >= 0 else 0) + cx)
        return state, dict(gates=one + two, single_qubit_gates=one, two_qubit_gates=two,
            unit_duration_ASAP_depth=max(depth, default=0), CX_dependency_depth=max(cx_depth, default=0),
            actual_Rz_rotations=rotations)


def metrics(exact, actual, model):
    norm = float(np.linalg.norm(actual))
    k.require(math.isfinite(norm) and norm > 0, 'Finite nonzero state norm')
    a, z = actual / norm, exact / np.linalg.norm(exact)
    overlap = np.vdot(z, a)
    loss = float(np.sum(np.abs(a - z * overlap)**2))
    hf, n = (1 << model['electrons']) - 1, model['n']
    target = [sum((hf >> q) & 1 for q in range(p, n, 2)) for p in (0, 1)]
    number = joint = 0.
    for v, value in enumerate(a):
        if v.bit_count() != model['electrons']:
            number += abs(value)**2
        if [sum((v >> q) & 1 for q in range(p, n, 2)) for p in (0, 1)] != target:
            joint += abs(value)**2
    return dict(error=loss / (1 + float(abs(overlap))), squared_infidelity=loss,
                norm=norm, number_leakage=float(number), joint_sector_leakage=float(joint))
