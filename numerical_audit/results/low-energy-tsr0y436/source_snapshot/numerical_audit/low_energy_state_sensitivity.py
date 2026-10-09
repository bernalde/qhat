#!/usr/bin/env python3
"""Predeclared low-energy-state control on eight corrected frozen Hamiltonians.

Read the completed controlled-block campaign's hexadecimal coefficients and
chronological schedules directly. Propagate all vectors in the full 1024-state
Hilbert space. State preparation fixes particle number and spin projection,
not total spin. No state or schedule is chosen using observed Trotter errors.
"""
from __future__ import annotations

import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/qhat-audit-numba')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-audit-mpl')

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis.state_metrics import state_error_metrics
from qhat.numerical_audit.initial_state_sensitivity import (
    independent_hamiltonian, physical_sector, sha, object_sha, state_sha, write_json,
)

ROOT = Path(__file__).resolve().parents[1]
AUDIT = Path(__file__).resolve().parent
SOURCE_RUN = AUDIT / 'results/controlled-blocks-b53wa_g3'
ORDERS = ('fermionic_signed', 'jw_signed', 'jw_magnitude')
SOURCE_ORDERS = dict(zip(ORDERS, ('block_signed', 'jw_signed', 'jw_magnitude')))
STATE_NAMES = ('sector_ground', 'sector_next', 'equal_coherent')
TOTAL_TIME, STEPS = 1.0, 100
TIE_ABSOLUTE = RATIO_FLOOR = 1e-12
DEGENERACY_ATOL = 1e-10
BASIS_ATOL = 1e-10
PHASE_ATOL = 1e-12
PROTOCOL = 'corrected-low-energy-states-v1'


def canonical_phase(vector):
    """Earliest component within 1e-12 of maximum magnitude is positive real."""
    vector = np.asarray(vector, dtype=complex).copy()
    magnitude = np.abs(vector)
    if np.max(magnitude) <= BASIS_ATOL:
        raise ValueError('Cannot phase a zero vector.')
    pivot = int(np.flatnonzero(magnitude >= np.max(magnitude) - PHASE_ATOL)[0])
    vector *= np.conjugate(vector[pivot]) / abs(vector[pivot])
    return vector


def canonical_subspace_basis(vectors):
    """Coordinate-projector basis, invariant under a unitary input basis change.

Project coordinate vectors in increasing sector basis order, twice modified
Gram--Schmidt, accept projected residual norms >1e-10, then fix phase. This
uses the eigenspace projector, not the eigensolver's arbitrary basis.
"""
    vectors = np.asarray(vectors, dtype=complex)
    projector = vectors @ vectors.conj().T
    chosen = []
    for column in range(projector.shape[1]):
        candidate = projector[:, column].copy()
        for _ in range(2):
            for previous in chosen:
                candidate -= np.vdot(previous, candidate) * previous
        norm = np.linalg.norm(candidate)
        if norm > BASIS_ATOL:
            chosen.append(canonical_phase(candidate / norm))
        if len(chosen) == vectors.shape[1]:
            break
    if len(chosen) != vectors.shape[1]:
        raise RuntimeError('Failed to resolve a numerical eigenspace.')
    result = np.column_stack(chosen)
    if np.linalg.norm(result.conj().T @ result - np.eye(len(chosen))) > 1e-10:
        raise RuntimeError('Canonical eigenspace basis is not orthonormal.')
    return result


def low_energy_states(matrix):
    """Return frozen canonical first two sector vectors and their superposition.

Sorted eigenvalues join a cluster only while its total span is <=1e-10 Ha.
Within such a numerical cluster, vectors are approximate eigenvectors with
residual bounded by the cluster span plus eigensolver roundoff. No replacement
is made if ground/next lie in the same cluster.
"""
    values, vectors = eigh(matrix)
    if len(values) < 2:
        raise ValueError('The sector must contain at least two states.')
    clusters = []
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[end] - values[start] <= DEGENERACY_ATOL:
            end += 1
        clusters.append({'start': start, 'stop': end, 'multiplicity': end-start,
                         'minimum_energy': float(values[start]), 'maximum_energy': float(values[end-1])})
        start = end
    canonical = []
    selected_clusters = []
    for cluster in clusters:
        basis = canonical_subspace_basis(vectors[:, cluster['start']:cluster['stop']])
        for j in range(basis.shape[1]):
            canonical.append(basis[:, j])
            selected_clusters.append({**cluster, 'canonical_index_within_cluster': j})
            if len(canonical) == 2:
                break
        if len(canonical) == 2:
            break
    first, second = canonical
    states = {'sector_ground': first, 'sector_next': second,
              'equal_coherent': (first + second) / np.sqrt(2.)}
    metadata = {'eigenvalues_identity_free_hartree': values.tolist(),
                'first_gap_hartree': float(values[1]-values[0]),
                'second_gap_hartree': float(values[2]-values[1]) if len(values) > 2 else None,
                'numerical_clusters': clusters, 'selected_eigenvector_clusters': selected_clusters,
                'selected_vectors_same_numerical_eigenspace': selected_clusters[0]['start'] == selected_clusters[1]['start'],
                'selected_vector_inner_product_abs': float(abs(np.vdot(first, second)))}
    return states, metadata, values, vectors


def decode_frozen_design(design):
    if design['coefficient_policy'] != 'molecular-unit-weight-v1' or design['coefficient_cutoff'] != 1e-12:
        raise RuntimeError('Expected the corrected explicit-cutoff coefficient policy.')
    coefficients = {tuple((int(q), p) for q, p in key): complex(float.fromhex(real), float.fromhex(imag))
                    for key, real, imag in design['identity_free_coefficients']}
    canonical = [[key, value.real.hex(), value.imag.hex()] for key, value in sorted(coefficients.items())]
    if object_sha(canonical) != design['hamiltonian_sha256']:
        raise RuntimeError('Frozen retained Hamiltonian hash differs.')
    if any(value.imag != 0 for value in coefficients.values()):
        raise RuntimeError('Nonreal Pauli coefficient.')
    raw_keys = [tuple((int(q), p) for q, p in key) for key in design['pauli_order_key_reference']]
    if len(raw_keys) != len(coefficients) or set(raw_keys) != set(coefficients):
        raise RuntimeError('Pauli key reference changes retained terms.')
    schedules, hashes = {}, {}
    for name, source_name in SOURCE_ORDERS.items():
        indices = design['pauli_orders'][source_name]
        if sorted(indices) != list(range(len(raw_keys))):
            raise RuntimeError('Invalid frozen schedule permutation.')
        keys = [raw_keys[i] for i in indices]
        digest = object_sha([[key, coefficients[key].real.hex(), coefficients[key].imag.hex()] for key in keys])
        if digest != design['schedule_sha256'][source_name]:
            raise RuntimeError('Frozen chronological schedule hash differs.')
        schedules[name], hashes[name] = keys, digest
    return coefficients, schedules, hashes


def prepare(case, output):
    source = SOURCE_RUN / f'{case}.design.json'
    design = json.loads(source.read_text())
    n, electrons = design['n_qubits'], design['n_electrons']
    if (n, electrons) != (10, 4):
        raise RuntimeError('This protocol requires ten qubits and four electrons.')
    coefficients, schedules, hashes = decode_frozen_design(design)
    for name, digest in design['input_fixture_hashes'].items():
        if sha(AUDIT / 'fixtures' / case / name) != digest:
            raise RuntimeError('Original input fixture changed.')
    shutil.copyfile(source, output / 'input_designs' / source.name)
    hamiltonian = independent_hamiltonian(coefficients, n)
    sector, alpha, beta = physical_sector(n, electrons)
    if len(sector) != 100 or (alpha, beta) != (2, 2):
        raise RuntimeError('Unexpected physical sector.')
    mask = np.zeros(2**n, dtype=bool)
    mask[sector] = True
    number_mask = np.array([index.bit_count() == electrons for index in range(2**n)])
    herm_error = float(sparse.linalg.norm(hamiltonian-hamiltonian.getH()))
    invariant_error = float(sparse.linalg.norm(hamiltonian[~mask][:, sector]))
    if max(herm_error, invariant_error) > 1e-10:
        raise RuntimeError('Invalid reference Hermiticity or sector invariance.')
    matrix = hamiltonian[sector][:, sector].toarray()
    small_states, spectrum, values, vectors = low_energy_states(matrix)
    states, state_info = {}, {}
    hf = np.zeros(2**n, complex)
    hf[((1 << electrons)-1) << (n-electrons)] = 1.
    for name, small in small_states.items():
        full = np.zeros(2**n, complex)
        full[sector] = small
        energy = float(np.vdot(full, hamiltonian @ full).real)
        residual = float(np.linalg.norm(hamiltonian @ full-energy*full))
        if name != 'equal_coherent' and residual > 2*DEGENERACY_ATOL:
            raise RuntimeError('Low-energy vector eigenvalue residual is too large.')
        states[name] = full
        state_info[name] = {'complex128_little_endian_sha256': state_sha(full), 'norm': float(np.linalg.norm(full)),
                            'energy_identity_free_hartree': energy, 'energy_standard_deviation_hartree': residual,
                            'hf_overlap_magnitude': float(abs(np.vdot(hf, full))),
                            'hf_fidelity': float(abs(np.vdot(hf, full))**2),
                            'eigenvalue_residual_interpretation': 'eigenvector residual' if name != 'equal_coherent' else 'energy standard deviation, generally nonzero'}
    prepared = {'case': case, 'protocol': PROTOCOL, 'source_design_sha256': sha(source),
                'coefficient_policy': design['coefficient_policy'], 'coefficient_cutoff': design['coefficient_cutoff'],
                'hamiltonian_sha256': design['hamiltonian_sha256'], 'ordered_schedule_sha256': hashes,
                'n_qubits': n, 'n_electrons': electrons, 'N_alpha': alpha, 'N_beta': beta,
                'sector_dimension': len(sector), 'full_hilbert_dimension': 2**n,
                'sector_convention': 'interleaved alpha/beta modes, q0 most significant; fixes particle number and spin projection, not total spin',
                'pauli_count': len(coefficients), 'nominal_exponential_count': len(coefficients)*STEPS,
                'time': TOTAL_TIME, 'steps': STEPS, 'state_preparation': state_info, 'spectrum': spectrum,
                'discarded_global_scalar': design['discarded_global_scalar'],
                'reference_hermiticity_frobenius': herm_error, 'reference_off_sector_frobenius': invariant_error}
    write_json(output / f'{case}.design.json', prepared)
    np.savez_compressed(output / f'{case}.initial_states.npz', sector_indices=sector, **states)
    compiled = {name: evolution.compile_ordered_terms(keys, coefficients, n, 1e-12) for name, keys in schedules.items()}
    if any(len(terms) != len(coefficients) for terms in compiled.values()):
        raise RuntimeError('Compilation changed retained Pauli count.')
    return hamiltonian, states, sector, mask, number_mask, compiled, prepared, values, vectors


def outcome(left, right):
    return 'tie' if abs(left-right) <= TIE_ABSOLUTE else ('win' if left < right else 'loss')


def summarize(rows, output, elapsed):
    comparisons = []
    for case in sorted({row['case'] for row in rows}):
        for state in STATE_NAMES:
            selected = {row['ordering']: row for row in rows if row['case'] == case and row['initial_state'] == state}
            if set(selected) != set(ORDERS):
                raise RuntimeError('Incomplete case/state comparison.')
            errors = {name: selected[name]['one_minus_overlap'] for name in ORDERS}
            best = min(errors.values())
            comparisons.append({'case': case, 'initial_state': state, 'E': errors,
                                'winner_set': [name for name in ORDERS if errors[name]-best <= TIE_ABSOLUTE],
                                'signed_vs_comparator': {name: outcome(errors[ORDERS[0]], errors[name]) for name in ORDERS[1:]},
                                'signed_over_comparator': {name: errors[ORDERS[0]]/errors[name] if errors[name] > RATIO_FLOOR else None for name in ORDERS[1:]}})
    counts = {}
    for state in STATE_NAMES + ('all_states',):
        selected = [c for c in comparisons if state == 'all_states' or c['initial_state'] == state]
        counts[state] = {'winner_counts_including_ties': dict(Counter(name for c in selected for name in c['winner_set'])),
                         'signed_win_tie_loss': {name: dict(Counter(c['signed_vs_comparator'][name] for c in selected)) for name in ORDERS[1:]}}
    summary = {'protocol': PROTOCOL, 'completed_records': len(rows), 'expected_records': 72,
               'elapsed_seconds': elapsed, 'counts_by_state': counts, 'comparisons': comparisons,
               'max_norm_defect': max(row['maximum_norm_defect'] for row in rows),
               'max_spin_projection_sector_leakage_by_order': {order: max(row['spin_sector_leakage'] for row in rows if row['ordering'] == order) for order in ORDERS},
               'max_particle_number_leakage_by_order': {order: max(row['particle_number_leakage'] for row in rows if row['ordering'] == order) for order in ORDERS},
               'winner_tie_absolute_tolerance': TIE_ABSOLUTE, 'ratio_floor': RATIO_FLOOR,
               'scope': 'Descriptive paired controls on eight frozen Hamiltonians from four molecular families. States and geometries are not independent molecular samples. Sector low-energy does not assert a total-spin quantum number or full-system chemical ground state.'}
    write_json(output / 'summary.json', summary)
    lines = ['# Corrected low-energy-state control', '',
             'Eight frozen ten-qubit Hamiltonians; T=1, r=100; three initial states and three fixed chronological Pauli schedules. Ground and next refer to the retained active-space Hamiltonian in the (N_alpha,N_beta)=(2,2) sector, not a specified total-spin sector. All vectors are propagated in the full 1024-dimensional space without postselection.', '',
             'The Hamiltonian and schedules come directly from controlled-blocks-b53wa_g3 with molecular-unit-weight-v1 coefficients. All initial vectors and designs were frozen before any tested Trotter evolution. Exact references use independent sparse full-space expm_multiply and sector diagonalization. The coherent state uses the predeclared canonical relative phase.', '',
             '| Initial state | Proposed vs JW signed W/T/L | Proposed vs JW magnitude W/T/L | Best: proposed / JW signed / JW magnitude |', '|---|---|---|---|']
    for state, data in counts.items():
        counts_text = []
        for name in ORDERS[1:]:
            c = data['signed_win_tie_loss'][name]
            counts_text.append(' / '.join(str(c.get(key, 0)) for key in ('win','tie','loss')))
        best = ' / '.join(str(data['winner_counts_including_ties'].get(name,0)) for name in ORDERS)
        lines.append(f'| {state} | {counts_text[0]} | {counts_text[1]} | {best} |')
    lines += ['', '| Case | State | Proposed E | JW signed E | JW magnitude E | Winner set |', '|---|---|---|---|---|---|']
    for c in comparisons:
        values = ' | '.join(f"{c['E'][name]:.9g}" for name in ORDERS)
        lines.append(f"| {c['case']} | {c['initial_state']} | {values} | {', '.join(c['winner_set'])} |")
    lines += ['', f'Completed {len(rows)}/72 Trotter evaluations in {elapsed:.2f} seconds; tie tolerance and denominator floor 1e-12. All values and final vectors are retained. Eigenvalue gaps, cluster decisions, state hashes, HF overlaps, norms, residuals, and reference checks are saved per case.', '', summary['scope'], '']
    (output / 'report.md').write_text('\n'.join(lines))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=AUDIT/'results')
    args = parser.parse_args()
    source_manifest = json.loads((SOURCE_RUN/'manifest.json').read_text())
    cases = sorted(case for case in source_manifest['cases'] if json.loads((SOURCE_RUN/f'{case}.design.json').read_text())['n_qubits'] == 10)
    if len(cases) != 8:
        raise RuntimeError('Expected exactly eight frozen ten-qubit cases.')
    args.output_root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='low-energy-', dir=args.output_root))
    for name in ('source_snapshot', 'input_designs'):
        (output/name).mkdir()
    print('OUTPUT', output, flush=True)
    sources = [Path(__file__), AUDIT/'test_low_energy_state_sensitivity.py', AUDIT/'initial_state_sensitivity.py'] + [ROOT/'analysis'/name for name in ('state_metrics.py', 'benchmark_b2_active_spaces_matrix_free.py', 'benchmark_b2_signed_coefficient_baseline.py', 'benchmark_L_sweep_trotter.py', 'benchmark_b2_coloring_robustness.py')]
    source_hashes = {str(path.relative_to(ROOT)): sha(path) for path in sources}
    for path in sources:
        target = output/'source_snapshot'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    write_json(output/'manifest.json', {
        'protocol': PROTOCOL, 'git_head': subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        'git_status_porcelain': subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True),
        'source_hashes': source_hashes, 'source_campaign': str(SOURCE_RUN.relative_to(ROOT)),
        'source_manifest_sha256': sha(SOURCE_RUN/'manifest.json'),
        'cases': cases, 'orders': ORDERS, 'source_order_mapping': SOURCE_ORDERS, 'states': STATE_NAMES,
        'time': TOTAL_TIME, 'steps': STEPS, 'expected_trotter_evaluations': 72,
        'degeneracy_cluster_max_span_hartree': DEGENERACY_ATOL, 'coordinate_basis_acceptance_norm': BASIS_ATOL,
        'phase_tie_absolute': PHASE_ATOL,
        'degeneracy_rule': 'Cluster sorted eigenvalues while total span <= 1e-10 Ha. In each cluster project ascending coordinate vectors with the eigenspace projector, twice modified Gram-Schmidt, accept norm > 1e-10. Keep all declared states, including ground/next in one cluster.',
        'phase_rule': 'Earliest coordinate within 1e-12 of largest magnitude made real-positive; coherent vector (ground+next)/sqrt(2).',
        'all_case_initial_vectors_frozen_before_any_trotter_evaluation': True,
        'metric': 'normalized_full_state_v2; E from phase-aligned residual; I=E(2-E); no postselection',
        'winner_tie_absolute_tolerance': TIE_ABSOLUTE, 'ratio_floor': RATIO_FLOOR,
        'product_order': 'chronological; subsequent exponentials act on the left', 'threads': 1,
        'python': sys.version, 'versions': {name:importlib.metadata.version(name) for name in ('numpy','scipy','openfermion','numba')},
    })
    started, rows = time.monotonic(), []
    try:
        prepared = {case: prepare(case, output) for case in cases}
        frozen = {path.name:sha(path) for path in sorted(output.glob('*.design.json'))}
        frozen.update({path.name:sha(path) for path in sorted(output.glob('*.initial_states.npz'))})
        write_json(output/'state_design_freeze.json', {'sha256': frozen, 'trotter_evaluations_completed': 0,
                                                     'elapsed_seconds': time.monotonic()-started})
        print('ALL_DESIGNS_FROZEN', output/'state_design_freeze.json', flush=True)
        with (output/'results.jsonl').open('x') as stream:
            for case in cases:
                H, states, sector, mask, number_mask, compiled, design, values, vectors = prepared[case]
                initial_matrix = np.column_stack([states[name] for name in STATE_NAMES])
                exact = expm_multiply(-1j*TOTAL_TIME*H, initial_matrix, traceA=-1j*TOTAL_TIME*H.diagonal().sum())
                small_exact = vectors @ (np.exp(-1j*TOTAL_TIME*values)[:,None]*(vectors.conj().T @ initial_matrix[sector]))
                embedded = np.zeros_like(exact)
                embedded[sector] = small_exact
                residuals = np.linalg.norm(embedded-exact, axis=0)
                if max(residuals) > 1e-10:
                    raise RuntimeError('Independent full-space and sector exact evolutions disagree.')
                write_json(output/f'{case}.reference_checks.json', {
                    'full_expm_vs_sector_eigh_residual': dict(zip(STATE_NAMES, residuals.tolist())),
                    'full_expm_final_state_sha256': {name:state_sha(exact[:,i]) for i,name in enumerate(STATE_NAMES)},
                    'exact_spin_projection_sector_leakage': {name:float(np.linalg.norm(exact[~mask,i])**2) for i,name in enumerate(STATE_NAMES)},
                    'hamiltonian_sparse_nonzeros': H.nnz})
                final_states = {f'exact__{name}':exact[:,i] for i,name in enumerate(STATE_NAMES)}
                for i, state in enumerate(STATE_NAMES):
                    for order in ORDERS:
                        before = time.monotonic()
                        approximate, count = evolution.evolve_trotter_state(states[state],compiled[order],1,STEPS,TOTAL_TIME,2**30)
                        metrics = state_error_metrics(exact[:,i],approximate)
                        norm_defect = max(abs(metrics['exact_norm']-1),abs(metrics['approx_norm']-1))
                        if norm_defect > 1e-8 or count != design['nominal_exponential_count']:
                            raise RuntimeError('Norm defect or exponential count failed.')
                        normalized = approximate/metrics['approx_norm']
                        row = {'case':case, 'initial_state':state, 'ordering':order, 'time':TOTAL_TIME, 'steps':STEPS,
                               'initial_state_sha256':state_sha(states[state]), 'approximate_final_state_sha256':state_sha(approximate),
                               'hamiltonian_sha256':design['hamiltonian_sha256'], 'schedule_sha256':design['ordered_schedule_sha256'][order],
                               'nominal_exponential_count':count, 'elapsed_seconds':time.monotonic()-before,
                               'particle_number_leakage':float(np.linalg.norm(normalized[~number_mask])**2),
                               'spin_sector_leakage':float(np.linalg.norm(normalized[~mask])**2), 'maximum_norm_defect':norm_defect, **metrics}
                        rows.append(row)
                        stream.write(json.dumps(row,allow_nan=False)+'\n')
                        stream.flush()
                        final_states[f'{order}__{state}'] = approximate
                np.savez_compressed(output/f'{case}.final_states.npz',**final_states)
                print('COMPLETE',case,'elapsed',round(time.monotonic()-started,2),flush=True)
        if any(sha(ROOT/path) != digest for path,digest in source_hashes.items()):
            raise RuntimeError('A recorded source changed during the experiment.')
        if any(sha(output/name) != digest for name,digest in frozen.items()):
            raise RuntimeError('A frozen state/design changed during the experiment.')
        if len(rows) != 72:
            raise RuntimeError('Incomplete experiment.')
        write_json(output/'results.json',rows)
        summary = summarize(rows,output,time.monotonic()-started)
        write_json(output/'completion.json',{'completed_records':len(rows),'elapsed_seconds':time.monotonic()-started,
                                            'all_frozen_state_and_source_hashes_unchanged':True,
                                            'result_file_sha256':sha(output/'results.json'), 'summary_file_sha256':sha(output/'summary.json')})
        print(json.dumps({'completed_records':summary['completed_records'],'elapsed_seconds':summary['elapsed_seconds'], 'counts':summary['counts_by_state']}),flush=True)
    except Exception as error:
        write_json(output/'failure.json',{'type':type(error).__name__,'message':str(error),
                                         'completed_records':len(rows),'elapsed_seconds':time.monotonic()-started})
        raise


if __name__ == '__main__':
    main()
