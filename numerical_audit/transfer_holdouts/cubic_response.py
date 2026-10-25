"""Symmetric cubic diagonal response of a normalized complex initial state.

Columns of the signed real cubic defect are transported coherently before the
state projection. No full-H dynamics or validation error is used by this score.
"""
import importlib.util
import math
from pathlib import Path
import time
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'symmetric_response/cubic.py'
spec = importlib.util.spec_from_file_location('transfer_prior_cubic', OLD)
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
engine = prior.engine
require = prior.require
Interrupted = prior.Interrupted


def validate(model):
    n = model['n']
    require(type(n) is int and 1 <= n <= 6, 'Register bound before allocation')
    require(type(model['electrons']) is int and 0 <= model['electrons'] <= n, 'Occupation bound')
    require(0 < len(model['blocks']) <= 256 and 0 < len(model['coefficients']) <= 256,
            'Unit/term bound before allocation')
    require(type(model['denominator']) in (int, str) and int(model['denominator']) > 0, 'Positive denominator')
    blocks = [prior.decode(raw) for raw in model['blocks']]
    total = prior.decode(model['coefficients'])
    require(all(blocks) and sum(map(len, blocks)) == len(total) and
            {w: c for b in blocks for w, c in b.items()} == total, 'Complete target partition required')
    for block in blocks:
        for word in block:
            require(word and len(dict(word)) == len(word) and
                    all(type(q) is int and 0 <= q < n and p in 'XYZ' for q, p in word), 'Invalid support')
        engine.descriptor(block)
    return blocks


class Context:
    def __init__(self, model, psi):
        blocks = validate(model)
        dimension = 1 << model['n']
        state = np.asarray(psi, complex)
        require(state.shape == (dimension,) and np.all(np.isfinite(state)) and
                abs(np.vdot(state, state) - 1) < 1e-10, 'Normalized full-space state required')
        support = np.flatnonzero(state != 0)
        require(0 < len(support) <= 9, 'Declared state-support bound before column allocation')
        require(len(blocks) * dimension * 24 + 12 * dimension * len(support) * 16 <= 8 * 1024**2,
                'Action/column storage guard')
        self.model = model
        self.m = len(blocks)
        self.blocks = blocks
        self.dimension = dimension
        self.psi = state.copy()
        self.support = support
        self.lengths = [len(b) for b in blocks]
        denominator = int(model['denominator'])
        basis = np.arange(dimension)
        self.actions = []
        diagonal = {}
        for block in blocks:
            descriptor = engine.descriptor(block)
            indices = basis ^ descriptor[0]
            values = np.array([engine.value(descriptor, int(v)) / denominator for v in basis], float)
            require(np.all(np.isfinite(values)), 'Nonfinite action coefficients')
            self.actions.append((indices, values))
            if descriptor[0] == 0:
                diagonal.update(block)
        desc_D = engine.descriptor(diagonal or {(): 0})
        self.energies = np.array([engine.value(desc_D, int(v)) / denominator for v in basis])
        self.initial_columns = np.eye(dimension)[:, support]
        self.H1 = self.hamiltonian_action(self.initial_columns)
        self.H2 = self.hamiltonian_action(self.H1) / 2
        self.H3 = self.hamiltonian_action(self.hamiltonian_action(self.H1)) / 6
        self.cache_bytes = sum(i.nbytes + v.nbytes for i, v in self.actions) + sum(
            a.nbytes for a in (self.psi, self.support, self.energies, self.initial_columns, self.H1, self.H2, self.H3))
        require(self.cache_bytes <= 8 * 1024**2, 'Actual action/column cache bound')

    def apply(self, unit, columns):
        indices, values = self.actions[unit]
        return values[indices, None] * columns[indices]

    def hamiltonian_action(self, columns):
        result = np.zeros_like(columns)
        for unit in range(self.m):
            result += self.apply(unit, columns)
        return result

    def cubic_columns(self, order, deadline=float('inf'), clock=time.process_time):
        prior.permutation(order, self.m)
        coefficients = [self.initial_columns.copy()] + [np.zeros_like(self.initial_columns) for _ in range(3)]
        completed = 0
        for unit, duration in prior.symmetric_units(order):
            if clock() >= deadline:
                raise Interrupted(completed)
            updated = [v.copy() for v in coefficients]
            powers = list(coefficients)
            factor = 1.
            for degree in range(1, 4):
                factor *= duration / degree
                for k in range(4 - degree):
                    powers[k] = self.apply(unit, powers[k])
                    updated[k + degree] += factor * powers[k]
            coefficients = updated
            completed += 1
        require(np.max(abs(coefficients[1] - self.H1)) < 2e-11 and
                np.max(abs(coefficients[2] - self.H2)) < 2e-11, 'Symmetric lower coefficients differ')
        return coefficients[3] - self.H3

    def score(self, order, T, rotation_cap, deadline=float('inf'), clock=time.process_time):
        require(math.isfinite(T) and T > 0 and type(rotation_cap) is int and rotation_cap > 0,
                'Positive finite time and cap required')
        columns = self.cubic_columns(order, deadline, clock)
        gaps = self.energies[:, None] - self.energies[self.support][None, :]
        kernel = T * np.sinc(T * gaps / (2 * np.pi)) * np.exp(.5j * T * gaps)
        action = np.sum(kernel * columns * self.psi[self.support][None, :], axis=1)
        projected = action - self.psi * np.vdot(self.psi, action)
        value = float(np.vdot(projected, projected).real)
        require(math.isfinite(value) and value >= 0, 'Invalid generalized cubic variance')
        steps = engine.feasible_steps([self.lengths[i] for i in order], rotation_cap, True)
        require(steps >= 1, 'No affordable symmetric step')
        return dict(order=list(order), JD2=value, score=value / steps**4, steps=steps,
                    column_sources=self.support.tolist(), cubic_columns=columns.tolist())
