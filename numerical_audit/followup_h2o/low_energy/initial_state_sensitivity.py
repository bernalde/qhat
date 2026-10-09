#!/usr/bin/env python3
"""Bounded exploratory state-sensitivity control for eight frozen 10-qubit cases.

For each case, compare the same HF state and five predeclared Haar-random
complex states in its fixed-(N_alpha,N_beta) sector under three frozen Pauli
orders. Exact and approximate propagation retain the entire 2**10 space.
The independent reference sums unit-Pauli sparse matrices with the immutable
hexadecimal coefficients, without any additional coefficient pruning.
"""
from __future__ import annotations

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/qhat-audit-numba')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-audit-mpl')

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from openfermion import QubitOperator, get_fermion_operator, get_sparse_operator, jordan_wigner
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis.state_metrics import state_error_metrics

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).with_name('fixtures')
ORDERS = ('fermionic_signed', 'jw_signed', 'jw_magnitude')
BASE_SEEDS = (2026100901, 2026100902, 2026100903, 2026100904, 2026100905)
TOTAL_TIME = 1.0
STEPS = 100
TIE_ABSOLUTE = 1e-12
_UNIT_PAULIS = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def state_sha(state):
    return hashlib.sha256(np.asarray(state, dtype='<c16').tobytes(order='C')).hexdigest()


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def unit_pauli(key, n):
    """Construct a unit-coefficient Pauli independently of compiled bit actions."""
    cache_key = (n, key)
    if cache_key not in _UNIT_PAULIS:
        _UNIT_PAULIS[cache_key] = get_sparse_operator(QubitOperator(key, 1.), n_qubits=n).tocsr()
    return _UNIT_PAULIS[cache_key]


def independent_hamiltonian(coefficients, n):
    hamiltonian = sparse.csr_matrix((2**n, 2**n), dtype=complex)
    for key, coefficient in coefficients.items():
        hamiltonian = hamiltonian + coefficient * unit_pauli(key, n)
    hamiltonian.sum_duplicates()
    # No compress/threshold/eliminate-small-entry operation is used.
    return hamiltonian


def physical_sector(n, electrons):
    """Independent enumeration with interleaved spin orbitals, q0 leftmost."""
    alpha = (electrons + 1) // 2
    beta = electrons // 2
    indices = []
    for index in range(2**n):
        occupations = [(index >> (n - 1 - q)) & 1 for q in range(n)]
        if sum(occupations[::2]) == alpha and sum(occupations[1::2]) == beta:
            indices.append(index)
    return np.asarray(indices, dtype=int), alpha, beta


def prepare(case, case_index, manifest, output):
    fixture = FIXTURES / case
    for name, expected in manifest[case]['files'].items():
        if sha(fixture / name) != expected:
            raise RuntimeError(f'Changed input fixture: {case}/{name}')
    setup = json.loads((fixture / 'setup.json').read_text())
    snapshot = json.loads((fixture / 'operators.json').read_text())
    n, electrons = setup['n_qubits'], setup['n_electrons']
    if n != 10:
        raise ValueError('This bounded experiment accepts only 10-qubit cases.')
    coefficients = {tuple((int(q), p) for q, p in key): complex(float.fromhex(real), float.fromhex(imag))
                    for key, real, imag in snapshot['coefficients']}
    canonical = [[key, value.real.hex(), value.imag.hex()] for key, value in sorted(coefficients.items())]
    if object_sha(canonical) != setup['fingerprints']['identity_free_jw_sha256']:
        raise RuntimeError(f'Frozen retained Hamiltonian hash differs: {case}')
    if any(value.imag != 0 for value in coefficients.values()):
        raise RuntimeError('A Hermitian Pauli expansion must have real coefficients.')

    # Recover historical insertion order from the pinned preprocessing recipe,
    # then require exact equality to the immutable coefficient and order hashes.
    interaction, regenerated_n = baseline.load_interaction_operator(fixture / manifest[case]['tensor'])
    fermion = baseline.clean_fermion_operator(get_fermion_operator(interaction), 1e-12)
    jw = jordan_wigner(fermion)
    jw.compress(abs_tol=1e-12)
    regenerated = {key: value for key, value in jw.terms.items() if key and abs(value) > 1e-12}
    if regenerated_n != n or regenerated != coefficients:
        raise RuntimeError(f'Regenerated retained coefficients differ: {case}')
    raw_keys = list(regenerated)
    ordered_keys = {}
    order_hashes = {}
    for name in ORDERS:
        order = snapshot['orders'][name]
        if sorted(order) != list(range(len(raw_keys))):
            raise RuntimeError(f'Invalid frozen ordering: {case}/{name}')
        keys = [raw_keys[i] for i in order]
        digest = object_sha([[key, coefficients[key].real.hex(), coefficients[key].imag.hex()] for key in keys])
        if digest != setup['fingerprints']['physical_order_sha256'][name]:
            raise RuntimeError(f'Frozen schedule hash differs: {case}/{name}')
        ordered_keys[name], order_hashes[name] = keys, digest

    sector, alpha, beta = physical_sector(n, electrons)
    states = []
    hf = np.zeros(2**n, dtype=complex)
    hf[((1 << electrons) - 1) << (n - electrons)] = 1.
    states.append(('hf', None, hf))
    for seed in BASE_SEEDS:
        rng = np.random.default_rng(np.random.SeedSequence([seed, case_index]))
        amplitudes = rng.normal(size=len(sector)) + 1j * rng.normal(size=len(sector))
        amplitudes /= np.linalg.norm(amplitudes)
        state = np.zeros(2**n, dtype=complex)
        state[sector] = amplitudes
        states.append((f'haar_{seed}', [seed, case_index], state))
    predeclared = {
        'case': case, 'case_index': case_index, 'n_qubits': n, 'n_electrons': electrons,
        'N_alpha': alpha, 'N_beta': beta, 'physical_sector_dimension': len(sector),
        'full_hilbert_dimension': 2**n, 'time': TOTAL_TIME, 'steps': STEPS,
        'pauli_count': len(coefficients), 'nominal_exponential_count': len(coefficients) * STEPS,
        'source_fixture_hashes': manifest[case]['files'], 'hamiltonian_sha256': object_sha(canonical),
        'ordered_schedule_sha256': order_hashes,
        'initial_states': [{'name': name, 'seed_sequence_entropy': seed,
                            'complex128_little_endian_sha256': state_sha(state)} for name, seed, state in states],
        'state_sampling': 'Independent complex Gaussian amplitudes normalized in the full fixed-spin sector; Haar on that sector. SeedSequence([base_seed, sorted_case_index]).',
        'identity_convention': 'Exactly the same frozen identity-free Hamiltonian for reference and each approximation.',
    }
    # Record all state seeds, hashes and schedules BEFORE any tested evolution.
    write_json(output / f'{case}.design.json', predeclared)
    np.savez(output / f'{case}.initial_states.npz', sector_indices=sector,
             **{name: state for name, _, state in states})

    hamiltonian = independent_hamiltonian(coefficients, n)
    hermiticity_error = float(sparse.linalg.norm(hamiltonian - hamiltonian.getH()))
    mask = np.zeros(2**n, dtype=bool)
    mask[sector] = True
    sector_invariance_error = float(sparse.linalg.norm(hamiltonian[~mask][:, sector]))
    if max(hermiticity_error, sector_invariance_error) > 1e-10:
        raise RuntimeError(f'Invalid reference symmetry/Hermiticity: {case}')
    compiled = {name: evolution.compile_ordered_terms(keys, coefficients, n, 1e-12)
                for name, keys in ordered_keys.items()}
    if any(len(terms) != len(coefficients) for terms in compiled.values()):
        raise RuntimeError('Compilation changed the retained Pauli count.')
    return hamiltonian, compiled, states, sector, mask, predeclared, hermiticity_error, sector_invariance_error


def summarize(rows, cases, output, elapsed):
    summaries = []
    for case in cases:
        records = [row for row in rows if row['case'] == case]
        by_state = {}
        for row in records:
            by_state.setdefault(row['initial_state'], {})[row['ordering']] = row
        winners = {}
        for name, group in by_state.items():
            if set(group) != set(ORDERS):
                continue
            best = min(row['one_minus_overlap'] for row in group.values())
            winners[name] = sorted(order for order, row in group.items()
                                   if row['one_minus_overlap'] - best <= TIE_ABSOLUTE)
        hf = winners.get('hf', [])
        random_winners = {name: winner for name, winner in winners.items() if name != 'hf'}
        counts = Counter(order for winner in random_winners.values() for order in winner)
        summaries.append({
            'case': case, 'hf_winner_set': hf,
            'haar_winner_sets': random_winners,
            'haar_winner_counts_including_ties': {order: counts[order] for order in ORDERS},
            'haar_states_with_different_winner_set_from_hf': sum(winner != hf for winner in random_winners.values()),
            'haar_states_where_no_hf_winner_is_best': sum(not set(winner).intersection(hf) for winner in random_winners.values()),
            'max_spin_sector_leakage_by_order': {
                order: max((row['spin_sector_leakage'] for row in records if row['ordering'] == order), default=0.)
                for order in ORDERS},
        })
    summary = {'cases': summaries, 'completed_records': len(rows), 'expected_records': len(cases) * 6 * 3,
               'elapsed_seconds': elapsed, 'winner_tie_absolute_tolerance': TIE_ABSOLUTE,
               'scope': 'Exploratory initial-state control with five predeclared Haar states; not an estimate of molecular ground-state or experimental performance.'}
    write_json(output / 'summary.json', summary)
    lines = ['# Initial-state sensitivity control', '',
             'At T=1 and r=100, each ordering evolves the same initial state under the same retained Hamiltonian with the same nominal Pauli-exponential count. Exact reference evolution and Trotter propagation both use the full 1024-dimensional space; there is no sector postselection.', '',
             'Five seeds were declared before evolution for each case. Complex Haar states cover the full fixed-(N_alpha,N_beta) sector, including components outside the HF-connected affine space. These are exploratory controls, not chemically representative low-energy states.', '',
             '| Case | HF winner(s) | Haar wins: fermionic / JW signed / JW magnitude | Haar states with no HF winner best |',
             '|---|---|---|---|']
    for case in summaries:
        counts = case['haar_winner_counts_including_ties']
        lines.append(f"| {case['case']} | {', '.join(case['hf_winner_set'])} | {counts[ORDERS[0]]} / {counts[ORDERS[1]]} / {counts[ORDERS[2]]} | {case['haar_states_where_no_hf_winner_is_best']} / {len(case['haar_winner_sets'])} |")
    lines += ['', f"Completed {len(rows)} / {summary['expected_records']} ordering-state comparisons in {elapsed:.1f} seconds. Winners use a predeclared absolute overlap-error tie tolerance of {TIE_ABSOLUTE:g}; ties count for every tied ordering.", '',
              'Conservation and accuracy answer different questions: retaining a physical symmetry does not establish a uniformly better approximation. These data only test the specified states, times, budgets, and historical retained Hamiltonians. Input hashes, source hashes, state hashes, exact-reference crosschecks, norm defects, and leakage are stored in the adjacent JSON files. No ordering is selected separately for each metric after inspecting the outcome.', '']
    (output / 'report.md').write_text('\n'.join(lines))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=Path(__file__).with_name('results'))
    parser.add_argument('--max-seconds', type=float, default=900.)
    args = parser.parse_args()
    if not np.isfinite(args.max_seconds) or args.max_seconds <= 0:
        parser.error('--max-seconds must be finite and positive')
    manifest = json.loads((FIXTURES / 'manifest.json').read_text())
    cases = sorted(case for case, entry in manifest.items() if entry['qubits'] == 10)
    if len(cases) != 8:
        raise RuntimeError('Expected exactly eight frozen 10-qubit cases.')
    args.output_root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='initial-state-', dir=args.output_root))
    print('OUTPUT', output, flush=True)
    sources = [Path(__file__), ROOT / 'analysis/state_metrics.py'] + [ROOT / 'analysis' / name for name in (
        'benchmark_L_sweep_trotter.py', 'benchmark_b2_active_spaces_matrix_free.py',
        'benchmark_b2_signed_coefficient_baseline.py', 'benchmark_b2_coloring_robustness.py')]
    source_hashes = {str(path.relative_to(ROOT)): sha(path) for path in sources}
    write_json(output / 'manifest.json', {
        'git_head': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        'source_hashes': source_hashes, 'fixture_manifest_sha256': sha(FIXTURES / 'manifest.json'),
        'python': sys.version, 'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'openfermion', 'numba')},
        'cases': cases, 'orders': ORDERS, 'base_seeds': BASE_SEEDS, 'time': TOTAL_TIME, 'steps': STEPS,
        'initial_state_count_per_case': 6, 'threads': 1, 'max_seconds': args.max_seconds,
        'product_order': 'Chronological: subsequent Pauli exponentials act on the left.',
        'metric': 'normalized_full_state_v2; no postselection',
        'winner_tie_absolute_tolerance': TIE_ABSOLUTE,
        'reference': 'Independent sum of unit-Pauli sparse matrices and frozen hexadecimal coefficients; full-space expm_multiply crosschecked against independent sector eigh.',
    })
    started = time.monotonic()
    rows = []
    try:
        with (output / 'results.jsonl').open('x') as stream:
            for case_index, case in enumerate(cases):
                if time.monotonic() - started > args.max_seconds:
                    raise TimeoutError('Initial-state sensitivity time budget reached.')
                hamiltonian, compiled, states, sector, mask, design, herm_error, invariant_error = prepare(case, case_index, manifest, output)
                initial_matrix = np.column_stack([state for _, _, state in states])
                exact_matrix = expm_multiply(-1j * TOTAL_TIME * hamiltonian, initial_matrix,
                                            traceA=-1j * TOTAL_TIME * hamiltonian.diagonal().sum())
                energies, vectors = eigh(hamiltonian[sector][:, sector].toarray())
                independent_sector = vectors @ (np.exp(-1j * TOTAL_TIME * energies)[:, None] *
                                                (vectors.conj().T @ initial_matrix[sector]))
                embedded = np.zeros_like(initial_matrix)
                embedded[sector] = independent_sector
                reference_errors = np.linalg.norm(embedded - exact_matrix, axis=0)
                if np.max(reference_errors) > 1e-10:
                    raise RuntimeError(f'Independent exact references disagree: {case}')
                write_json(output / f'{case}.reference_checks.json', {
                    'hermiticity_frobenius': herm_error, 'sector_offblock_frobenius': invariant_error,
                    'full_expm_vs_sector_eigh_state_residuals': dict(zip([name for name, _, _ in states], reference_errors.tolist())),
                    'hamiltonian_sparse_nonzeros': hamiltonian.nnz,
                })
                for index, (state_name, seed, initial) in enumerate(states):
                    for order in ORDERS:
                        if time.monotonic() - started > args.max_seconds:
                            raise TimeoutError('Initial-state sensitivity time budget reached.')
                        before = time.monotonic()
                        approximate, count = evolution.evolve_trotter_state(initial, compiled[order], 1, STEPS, TOTAL_TIME, 2**30)
                        metrics = state_error_metrics(exact_matrix[:, index], approximate)
                        norm_defect = max(abs(metrics['exact_norm'] - 1), abs(metrics['approx_norm'] - 1))
                        if norm_defect > 1e-8:
                            raise RuntimeError(f'Unacceptable norm drift: {case}/{state_name}/{order}')
                        if count != design['nominal_exponential_count']:
                            raise RuntimeError('The actual nominal exponential budget differs.')
                        normalized = approximate / metrics['approx_norm']
                        spin_leakage = float(np.vdot(normalized[~mask], normalized[~mask]).real)
                        row = {'case': case, 'initial_state': state_name, 'seed_sequence_entropy': seed,
                               'initial_state_sha256': state_sha(initial), 'ordering': order,
                               'time': TOTAL_TIME, 'steps': STEPS, 'nominal_exponential_count': count,
                               'elapsed_seconds': time.monotonic() - before,
                               'spin_sector_leakage': spin_leakage, 'maximum_norm_defect': norm_defect,
                               **metrics}
                        rows.append(row)
                        stream.write(json.dumps(row, allow_nan=False) + '\n')
                        stream.flush()
                print('COMPLETE', case, 'elapsed', round(time.monotonic() - started, 1), flush=True)
        if any(sha(ROOT / path) != digest for path, digest in source_hashes.items()):
            raise RuntimeError('A recorded source file changed during this experiment.')
        write_json(output / 'results.json', rows)
        summary = summarize(rows, cases, output, time.monotonic() - started)
        print(json.dumps({'completed_records': summary['completed_records'], 'elapsed_seconds': summary['elapsed_seconds']}), flush=True)
    except Exception as error:
        write_json(output / 'failure.json', {'type': type(error).__name__, 'message': str(error),
                                           'completed_records': len(rows), 'elapsed_seconds': time.monotonic() - started})
        summarize(rows, cases, output, time.monotonic() - started)
        raise


if __name__ == '__main__':
    main()
