#!/usr/bin/env python3
"""Bounded exact structural-certificate and fixed-circuit compilation audit.

This does not select an ordering by its measured error. Its certified swaps
preserve the original product formula, and its resource counts use one stated
all-to-all parity-ladder synthesis without optimization or routing.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '2'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-followup-mpl')

import argparse
from dataclasses import dataclass
import hashlib
import importlib.metadata
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import resource
import shutil

_soft, _hard = resource.getrlimit(resource.RLIMIT_AS)
_cap = 6 * 1024**3
resource.setrlimit(resource.RLIMIT_AS,
                   (min(_cap, _hard) if _hard >= 0 else _cap, _hard))
import cirq
import numpy as np
from openfermion import FermionOperator, hermitian_conjugated, jordan_wigner
from scipy.linalg import expm


def require(value, message):
    if not value:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Parent:
    creators: tuple
    annihilators: tuple
    spectators: tuple

    @property
    def active(self):
        return frozenset(self.creators + self.annihilators)

    def name(self):
        return f'A{self.creators}_C{self.annihilators}_R{self.spectators}'

    def operator(self):
        require(len(self.creators) == len(self.annihilators), 'Unbalanced parent')
        require(len(set(self.creators + self.annihilators + self.spectators))
                == len(self.creators + self.annihilators + self.spectators),
                'Canonical site sets must be disjoint')
        spectator = FermionOperator((), 1)
        for p in self.spectators:
            spectator *= FermionOperator(((p, 1), (p, 0)), 1)
        if not self.active:
            return spectator
        word = tuple((p, 1) for p in self.creators)
        word += tuple((p, 0) for p in reversed(self.annihilators))
        excitation = FermionOperator(word, 1)
        return spectator * (excitation + hermitian_conjugated(excitation))


def structural_certificate(left, right):
    """Sufficient certificate; None means unknown, never noncommuting."""
    if left.active == right.active:
        return 'equal_active_mask'
    if (left.active.isdisjoint(right.active)
            and left.active.isdisjoint(right.spectators)
            and right.active.isdisjoint(left.spectators)):
        return 'disjoint_active_no_cross_spectators'
    return None


def commute_words(a, b):
    right = dict(b)
    return sum(p in right and right[p] != v for p, v in a) % 2 == 0


def ladder(p, creation, n):
    out = np.zeros((1 << n, 1 << n), dtype=np.int64)
    bit = 1 << (n - p - 1)
    for source in range(1 << n):
        if bool(source & bit) != bool(creation):
            out[source ^ bit, source] = (-1) ** ((source >> (n - p)).bit_count())
    return out


def car_matrix(op, n):
    out = np.zeros((1 << n, 1 << n), dtype=np.int64)
    for word, coeff in op.terms.items():
        require(complex(coeff).imag == 0 and float(complex(coeff).real).is_integer(),
                'Exact matrix oracle requires integral parent coefficients')
        term = np.eye(1 << n, dtype=np.int64)
        for p, create in word:
            term = term @ ladder(p, create, n)
        out += int(complex(coeff).real) * term
    return out


def pauli_matrix(word, n):
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    labels = dict(word)
    out = np.ones((1, 1), dtype=complex)
    for p in range(n):
        out = np.kron(out, matrices[labels.get(p, 'I')])
    return out


def powerset(values):
    return [group for k in range(len(values) + 1) for group in combinations(values, k)]


def canonical_parents(n):
    parents = [Parent((), (), r) for r in powerset(tuple(range(n)))]
    for size in range(1, n // 2 + 1):
        for active in combinations(range(n), 2 * size):
            for creators in combinations(active, size):
                annihilators = tuple(p for p in active if p not in creators)
                if creators > annihilators:
                    continue
                remaining = tuple(p for p in range(n) if p not in active)
                parents.extend(Parent(creators, annihilators, r) for r in powerset(remaining))
    return parents


def certificate_audit(n=4):
    parents = canonical_parents(n)
    matrices = [car_matrix(p.operator(), n) for p in parents]
    descendants = [jordan_wigner(p.operator()).terms for p in parents]
    for matrix, terms in zip(matrices, descendants):
        mapped = sum((c * pauli_matrix(p, n) for p, c in terms.items()),
                     np.zeros_like(matrix, dtype=complex))
        require(np.array_equal(mapped, matrix), 'Actual JW expansion mismatch')
        require(all(commute_words(a, b) for a, b in combinations(terms, 2)),
                'Canonical parent internal commutation failed')
    counts = {'parents': len(parents), 'pairs_with_replacement': 0,
              'structurally_certified': 0, 'unknown': 0,
              'zero_product_cross_commuting': 0, 'zero_product_cross_noncommuting': 0,
              'commuting_nonzero_product': 0}
    witnesses = {}
    for i, j in combinations_with_replacement(range(len(parents)), 2):
        counts['pairs_with_replacement'] += 1
        forward, reverse = matrices[i] @ matrices[j], matrices[j] @ matrices[i]
        cross = all(commute_words(a, b) for a in descendants[i] for b in descendants[j])
        certificate = structural_certificate(parents[i], parents[j])
        counts['structurally_certified' if certificate else 'unknown'] += 1
        if certificate:
            require(cross, 'False-positive structural certificate')
        if np.array_equal(forward, reverse) and np.any(forward):
            counts['commuting_nonzero_product'] += 1
            require(cross, 'Canonical nonzero-product implication failed')
        if not np.any(forward) and not np.any(reverse):
            name = 'zero_product_cross_commuting' if cross else 'zero_product_cross_noncommuting'
            counts[name] += 1
            witnesses.setdefault(name, [parents[i].name(), parents[j].name()])
    require(counts['zero_product_cross_commuting'] > 0
            and counts['zero_product_cross_noncommuting'] > 0, 'Missing zero-product controls')
    return {'counts': counts, 'witnesses': witnesses,
            'interpretation': 'Exhaustive finite four-mode check; general guarantee comes from Lean.'}


def rotation_operations(key, angle, qubits):
    require(bool(key), 'Identity rotation not needed in this fixed example')
    active = [qubits[p] for p, _ in key]
    before, after = [], []
    for p, label in key:
        if label == 'X':
            before.append(cirq.H(qubits[p])); after.append(cirq.H(qubits[p]))
        elif label == 'Y':
            before.extend([cirq.S(qubits[p])**-1, cirq.H(qubits[p])])
            after.extend([cirq.H(qubits[p]), cirq.S(qubits[p])])
    parity = [cirq.CNOT(a, b) for a, b in zip(active[:-1], active[1:])]
    return before + parity + [cirq.rz(2 * angle)(active[-1])] + parity[::-1] + after


def circuit(terms, time=0.17, n=4):
    qubits = cirq.LineQubit.range(n)
    out = cirq.Circuit()
    for key, coefficient in terms:
        for op in rotation_operations(key, time * float(complex(coefficient).real), qubits):
            out.append(op, strategy=cirq.InsertStrategy.NEW)
    return out, out.unitary(qubit_order=qubits)


def compilation_audit():
    n = 4
    parents = [Parent((0,), (1,), ()), Parent((0,), (1,), (2,)),
               Parent((0,), (1,), (3,))]
    weights = [1, 2, 3]
    require(all(structural_certificate(a, b) for a, b in combinations(parents, 2)),
            'Compilation region is not structurally certified')
    raw, coefficients = [], {}
    for parent, weight in zip(parents, weights):
        for key, value in sorted(jordan_wigner(weight * parent.operator()).terms.items()):
            raw.append((key, value))
            coefficients[key] = coefficients.get(key, 0) + value
    tail_parent = Parent((1,), (2,), ())
    tail = sorted(jordan_wigner(tail_parent.operator()).terms.items())
    merged = [(key, value) for key, value in sorted(coefficients.items()) if value != 0]
    raw_circuit, U = circuit(raw + tail)
    merged_circuit, V = circuit(merged + tail)
    expected = np.eye(1 << n, dtype=complex)
    for parent, weight in zip(parents + [tail_parent], weights + [1]):
        expected = expm(-.17j * weight * car_matrix(parent.operator(), n)) @ expected
    require(np.linalg.norm(U - V, 2) < 1e-12, 'Certified merged circuit changed the operator')
    require(np.linalg.norm(U - expected, 2) < 1e-12, 'Parity-ladder synthesis disagrees with CAR')
    def costs(c, terms):
        return {'pauli_rotations': len(terms),
                'CNOT': sum(op.gate == cirq.CNOT for op in c.all_operations()),
                'one_qubit_gates_including_Rz': sum(len(op.qubits) == 1 for op in c.all_operations()),
                'serialized_depth': len(c)}
    return {'time': .17, 'hardware': 'All-to-all connectivity, exact H/S/CNOT/Rz, no routing',
            'compiler': 'Explicit fixed parity ladder; every operation gets its own moment; no optimizer',
            'chronology': 'h01, 2 n2 h01, 3 n3 h01, h12; factors act in listed order',
            'before': costs(raw_circuit, raw + tail), 'after': costs(merged_circuit, merged + tail),
            'unitary_difference_operator_norm': float(np.linalg.norm(U - V, 2)),
            'independent_CAR_formula_difference': float(np.linalg.norm(U - expected, 2)),
            'raw_terms': [[key, complex(value).real] for key, value in raw + tail],
            'merged_terms': [[key, complex(value).real] for key, value in merged + tail],
            'before_cirq': str(raw_circuit), 'after_cirq': str(merged_circuit),
            'scope': 'One fixed illustrative circuit; no accuracy or general hardware-resource advantage claimed.'}


def ownership_graph_witness():
    parents = [Parent((), (), (0,)), Parent((0,), (1,), ()), Parent((), (), (0, 1))]
    terms = [jordan_wigner(p.operator()).terms for p in parents]
    full = {}
    for group in terms:
        for key, value in group.items():
            full[key] = full.get(key, 0) + value
    seen, owned = set(), []
    for group in terms:
        block = {k: full[k] for k in group if k and full.get(k, 0) != 0 and k not in seen}
        seen.update(block)
        owned.append(sum((v * pauli_matrix(k, 2) for k, v in block.items()), np.zeros((4, 4), complex)))
    source = [car_matrix(p.operator(), 2) for p in parents]
    require(np.array_equal(source[1] @ source[2], source[2] @ source[1]), 'Source witness changed')
    changed = owned[1] @ owned[2] - owned[2] @ owned[1]
    require(np.any(changed != 0), 'Ownership graph witness disappeared')
    return {'original_commutator_zero': True,
            'owned_commutator_operator_norm': float(np.linalg.norm(changed, 2)),
            'meaning': 'Original-parent graph certificates cannot be reused for redistributed owned blocks.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.output / 'executed_source.py')
    result = {'structural': certificate_audit(), 'compilation': compilation_audit(),
              'ownership_graph_control': ownership_graph_witness(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'packages': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'openfermion', 'cirq-core')}}
    (args.output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'structural': result['structural'],
                      'before': result['compilation']['before'], 'after': result['compilation']['after']}, indent=2))


if __name__ == '__main__':
    main()
