#!/usr/bin/env python3
"""Controlled permutations of one frozen first-owner block decomposition.

Ownership is computed exactly once per input using signed canonical-parent
priority. Every block treatment permutes those same completed owned blocks;
no treatment recomputes ownership. See controlled_block_protocol.md.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
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
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis.fermionic_coefficients import PRECISION_POLICY_VERSION, fermion_from_interaction
from qhat.analysis.fermionic_ownership import build_molecular_ownership
from qhat.numerical_audit import reference_tools as ref
from qhat.numerical_audit.run_frozen_audit import (
    canonical, evolve, number_commutator_bound, object_sha, physical_metrics, sha, write,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).with_name('fixtures')
PROTOCOL_VERSION = 'frozen-owned-block-permutations-v1'
RANDOM_SEEDS = tuple(range(2026100901, 2026100921))
TIMES = (0.25, 1.0)
PRIMARY_STEPS = 100
CONVERGENCE_CASES = ('HCN_s1.00', 'li2', 'f2')
CONVERGENCE_STEPS = (50, 200)
TIE_ABSOLUTE = 1e-12
ANALYSIS_FLOOR = 1e-12
ATOL = 1e-12
BLOCK_DETERMINISTIC = ('block_signed', 'block_parent_magnitude', 'block_l1', 'block_reverse_signed')
EXTERNAL = ('jw_signed', 'jw_magnitude')


def state_sha(vector):
    return hashlib.sha256(np.asarray(vector, dtype='<c16').tobytes(order='C')).hexdigest()


def family(case):
    return case.split('_s')[0].upper()


def case_seed_word(case):
    return int(hashlib.sha256(case.encode()).hexdigest()[:8], 16)


def plan_block_orders(owned, case, seeds=RANDOM_SEEDS):
    """Return metadata and permutations; never rebuild/reassign ownership."""
    blocks = owned.owned_blocks
    coefficients = owned.final_coefficients
    metadata = []
    for index, block in enumerate(blocks):
        parent = owned.parents[block.parent_index]
        tie_key = tuple(sorted(parent.component_keys))
        weight = complex(parent.operator.terms[tie_key[0]])
        if weight.imag != 0:
            raise ValueError('Signed priority requires exactly real canonical parent weights.')
        metadata.append({
            'block_index': index, 'parent_index': block.parent_index,
            'parent_component_keys': tie_key, 'signed_parent_weight': weight.real,
            'signed_parent_weight_hex': weight.real.hex(),
            'retained_block_l1': float(sum(abs(coefficients[key]) for key in block.pauli_keys)),
            'pauli_keys': block.pauli_keys,
        })
    indices = list(range(len(blocks)))
    signed = sorted(indices, key=lambda i: (metadata[i]['signed_parent_weight'], metadata[i]['parent_component_keys']))
    # The constructor's initial schedule is also the signed treatment.
    if signed != indices:
        raise ValueError('Ownership must be frozen from signed-ascending parent priority.')
    orders = {
        'block_signed': signed,
        'block_parent_magnitude': sorted(indices, key=lambda i: (-abs(metadata[i]['signed_parent_weight']), metadata[i]['parent_component_keys'])),
        'block_l1': sorted(indices, key=lambda i: (-metadata[i]['retained_block_l1'], metadata[i]['parent_component_keys'])),
        'block_reverse_signed': list(reversed(signed)),
    }
    seed_metadata = {}
    for seed in seeds:
        entropy = [int(seed), case_seed_word(case)]
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy)))
        name = f'block_random_{seed}'
        orders[name] = list(map(int, rng.permutation(len(blocks))))
        seed_metadata[name] = {'base_seed': int(seed), 'seed_sequence_entropy': entropy, 'bit_generator': 'PCG64'}
    expected = list(range(len(blocks)))
    if any(sorted(order) != expected for order in orders.values()):
        raise RuntimeError('Invalid whole-block permutation.')
    return metadata, orders, seed_metadata


def flatten_blocks(owned, order):
    return [key for i in order for key in owned.owned_blocks[i].pauli_keys]


def planned_settings(case):
    settings = [(total_time, PRIMARY_STEPS, 'primary') for total_time in TIMES]
    if case in CONVERGENCE_CASES:
        settings += [(1., steps, 'convergence_anchor') for steps in CONVERGENCE_STEPS]
    return settings


def source_files():
    return [Path(__file__), Path(ref.__file__), Path(__file__).with_name('run_frozen_audit.py'),
            Path(__file__).with_name('controlled_block_protocol.md'),
            Path(__file__).with_name('test_controlled_block_orders.py'),
            Path(__file__).with_name('plot_controlled_block_orders.py')] + [ROOT / 'analysis' / name for name in (
        'fermionic_coefficients.py', 'fermionic_ownership.py', 'state_metrics.py',
        'benchmark_b2_active_spaces_matrix_free.py', 'benchmark_b2_signed_coefficient_baseline.py',
        'benchmark_L_sweep_trotter.py', 'benchmark_b2_coloring_robustness.py',
    )]


def prepare(case, manifest, output):
    fixture = manifest[case]
    for name, digest in fixture['files'].items():
        if sha(FIXTURES / case / name) != digest:
            raise RuntimeError(f'Changed input: {case}/{name}')
    tensor = FIXTURES / case / fixture['tensor']
    interaction, n = baseline.load_interaction_operator(tensor)
    # This is the sole ownership call in the entire per-case preparation.
    owned = build_molecular_ownership(fermion_from_interaction(interaction), atol=ATOL, ordering='signed_ascending')
    coefficients = owned.final_coefficients
    if any(complex(c).imag != 0 for c in coefficients.values()):
        raise RuntimeError('Nonreal retained Pauli coefficient; no real projection is permitted.')
    if owned.fallback_pauli_keys or owned.diagnostics['noncommuting_owned_block_count']:
        raise RuntimeError('The fixed decomposition is not exhaustive commuting ownership.')
    block_info, block_orders, seeds = plan_block_orders(owned, case)
    schedules = {name: flatten_blocks(owned, order) for name, order in block_orders.items()}
    keys = list(coefficients)
    schedules.update({
        'jw_signed': baseline.signed_coefficient_lexicographic_order(keys, coefficients, n, ATOL),
        'jw_magnitude': baseline.magnitude_descending_lexicographic_order(keys, coefficients, n, ATOL),
    })
    if any(len(order) != len(keys) or set(order) != set(keys) for order in schedules.values()):
        raise RuntimeError('A schedule changes the common retained Hamiltonian.')
    if schedules['block_signed'] != list(owned.pauli_order):
        raise RuntimeError('Signed treatment differs from frozen ownership baseline.')
    raw = evolution.compile_ordered_terms(keys, coefficients, n, ATOL)
    if len(raw) != len(keys):
        raise RuntimeError('Pauli compilation discarded retained coefficients.')
    full_initial = evolution.build_hartree_fock_state(n, fixture['electrons'])
    hf_index = int(np.argmax(abs(full_initial)))
    states, reduced, _ = ref.affine_space(raw, hf_index)
    H, _ = ref.assemble_hamiltonian(reduced, len(states))
    initial = np.zeros(len(states), complex)
    initial[0] = 1.
    key_index = {key: i for i, key in enumerate(keys)}
    compiled = {name: [reduced[key_index[key]] for key in order] for name, order in schedules.items()}
    nm = np.array([int(x).bit_count() == fixture['electrons'] for x in states])
    physical = set(map(int, evolution.spin_sector_basis_indices(n, fixture['electrons'])))
    sm = np.array([int(x) in physical for x in states])
    membership = [[b.parent_index, b.pauli_keys] for b in owned.owned_blocks]
    design = {
        'case': case, 'family': family(case), 'protocol_version': PROTOCOL_VERSION,
        'n_qubits': n, 'n_electrons': fixture['electrons'],
        'pauli_terms': len(keys), 'nonempty_owned_blocks': len(owned.owned_blocks),
        'tensor_sha256': sha(tensor), 'input_fixture_hashes': fixture['files'],
        'coefficient_policy': PRECISION_POLICY_VERSION, 'coefficient_cutoff': ATOL,
        'hamiltonian_sha256': object_sha(canonical(coefficients)),
        'identity_free_coefficients': canonical(coefficients),
        'discarded_global_scalar': [complex(owned.scalar).real, complex(owned.scalar).imag],
        'ownership_policy': 'signed ascending canonical parent once; fixed thereafter',
        'membership_sha256': object_sha(membership), 'blocks': block_info,
        'block_orders': block_orders, 'random_seed_specifications': seeds,
        'pauli_orders': {name: [key_index[key] for key in order] for name, order in schedules.items()},
        'pauli_order_key_reference': keys,
        'schedule_sha256': {name: object_sha([[key, complex(coefficients[key]).real.hex(), complex(coefficients[key]).imag.hex()] for key in order]) for name, order in schedules.items()},
        'schedule_membership_sha256': {name: object_sha(membership) for name in block_orders},
        'affine_basis_indices': list(map(int, states)), 'closed_dimension': len(states),
        'full_hilbert_dimension': len(full_initial), 'physical_sector_in_closure_dimension': int(sum(sm)),
        'initial_state': 'Hartree-Fock; q0 is most significant bit',
        'initial_full_state_sha256': state_sha(full_initial), 'initial_closed_state_sha256': state_sha(initial),
        'planned_settings': [{'T': t, 'steps': r, 'role': role} for t, r, role in planned_settings(case)],
        'construction_diagnostics': owned.diagnostics,
    }
    # This complete schedule/state design is saved before any evolution.
    write(output / f'{case}.design.json', design)
    checks = {
        'hermiticity_frobenius': float(sparse.linalg.norm(H-H.getH())),
        'reference_off_spin_sector_frobenius': float(sparse.linalg.norm(H[~sm][:, sm])),
        'full_space_one_step_residuals': {},
        'max_block_number_commutator_l1': {},
    }
    block_indices = [[key_index[key] for key in b.pauli_keys] for b in owned.owned_blocks]
    for parity in (None, 0, 1):
        checks['max_block_number_commutator_l1'][str(parity)] = max(
            (number_commutator_bound(keys, coefficients, bucket, parity) for bucket in block_indices), default=0.)
    if checks['hermiticity_frobenius'] > 1e-10:
        raise RuntimeError('Reference Hamiltonian is not Hermitian.')
    for name in ('block_signed', f'block_random_{RANDOM_SEEDS[0]}', 'jw_magnitude'):
        full = evolve(full_initial, [raw[key_index[key]] for key in schedules[name]], .013, 1)
        small = evolve(initial, compiled[name], .013, 1)
        embedded = np.zeros_like(full)
        embedded[states] = small
        error = float(np.linalg.norm(embedded-full))
        checks['full_space_one_step_residuals'][name] = error
        if error > 1e-11:
            raise RuntimeError('Affine evolution does not match full Hilbert-space evolution.')
    write(output / f'{case}.checks.json', checks)
    return H, initial, compiled, nm, sm, design, checks


def random_comparison(error, random_errors):
    values = np.asarray(random_errors, dtype=float)
    better = int(np.sum(values < error-TIE_ABSOLUTE))
    tied = int(np.sum(abs(values-error) <= TIE_ABSOLUTE))
    worse = int(np.sum(values > error+TIE_ABSOLUTE))
    return {
        'random_count': len(values), 'random_better': better, 'random_tied': tied,
        'random_worse': worse, 'signed_midrank_among_signed_and_random': 1+better+tied/2,
        'signed_midpercentile_against_random': (better+tied/2)/len(values),
        'random_median_E': float(np.median(values)),
        'random_q25_E': float(np.quantile(values, .25)), 'random_q75_E': float(np.quantile(values, .75)),
        'random_min_E': float(min(values)), 'random_max_E': float(max(values)),
        'signed_over_random_median': error/float(np.median(values)) if np.median(values) > ANALYSIS_FLOOR else None,
        'signed_vs_random_median': 'tie' if abs(error-float(np.median(values))) <= TIE_ABSOLUTE else ('win' if error < np.median(values) else 'loss'),
    }


def paired_outcome(left, right):
    if abs(left-right) <= TIE_ABSOLUTE:
        return 'tie'
    return 'win' if left < right else 'loss'


def summarize(rows, output):
    groups = {}
    for row in rows:
        groups.setdefault((row['case'], row['T'], row['steps']), {})[row['ordering']] = row
    comparisons = []
    for (case, total_time, steps), group in groups.items():
        if len(group) != len(BLOCK_DETERMINISTIC)+len(RANDOM_SEEDS)+len(EXTERNAL):
            raise RuntimeError('Incomplete comparison setting.')
        signed = group['block_signed']['one_minus_overlap']
        random_errors = [group[f'block_random_{seed}']['one_minus_overlap'] for seed in RANDOM_SEEDS]
        comparisons.append({
            'case': case, 'family': family(case), 'T': total_time, 'steps': steps,
            'role': group['block_signed']['role'], 'signed_E': signed,
            'deterministic_E': {name: group[name]['one_minus_overlap'] for name in BLOCK_DETERMINISTIC+EXTERNAL},
            'signed_vs_comparator': {name: paired_outcome(signed, group[name]['one_minus_overlap']) for name in BLOCK_DETERMINISTIC[1:]+EXTERNAL},
            **random_comparison(signed, random_errors),
            'max_block_spin_leakage': max(row['spin_sector_leakage'] for name, row in group.items() if name.startswith('block_')),
            'max_block_number_leakage': max(row['particle_number_leakage'] for name, row in group.items() if name.startswith('block_')),
        })
    by_time = {}
    for total_time in TIMES:
        primary = [c for c in comparisons if c['role'] == 'primary' and c['T'] == total_time]
        counts = {name: dict(Counter(c['signed_vs_comparator'][name] for c in primary)) for name in BLOCK_DETERMINISTIC[1:]+EXTERNAL}
        counts['random_median'] = dict(Counter(c['signed_vs_random_median'] for c in primary))
        by_time[str(total_time)] = {'cases': len(primary), 'families': sorted({c['family'] for c in primary}), 'signed_win_tie_loss': counts}
    summary = {
        'protocol_version': PROTOCOL_VERSION, 'rows': len(rows), 'settings': len(comparisons),
        'cases': sorted({r['case'] for r in rows}), 'families': sorted({r['family'] for r in rows}),
        'tie_absolute': TIE_ABSOLUTE, 'analysis_floor': ANALYSIS_FLOOR,
        'random_permutations_per_case': len(RANDOM_SEEDS), 'primary_by_time': by_time,
        'comparisons': comparisons,
        'inference_scope': 'Descriptive paired case results. Geometries, times, steps and random schedules are not independent molecular samples. No population p-value and no best-random algorithm claim.',
    }
    write(output / 'summary.json', summary)
    lines = ['# Controlled whole-block ordering experiment', '',
        'Each Hamiltonian uses one signed-priority ownership decomposition. Its nonempty completed blocks, all retained Pauli coefficients, and each block\'s internal Pauli order are held fixed. Only whole-block chronology changes. The two JW comparators are external controls; they change contiguity as well as order.', '',
        'The primary error is stable normalized full-state E = 1 - |<exact, approximate>|, evaluated from the phase-aligned residual. No particle/spin postselection is used. Every comparison has the same Hamiltonian, HF state, total time, and Pauli-exponential count. Twenty predeclared uniform random block permutations supply a reference distribution, not a best-of-twenty algorithm.', '',
        f'Completed {len(rows)} evolutions across {len(summary["cases"])} Hamiltonians and {len(summary["families"])} molecular families. Absolute tie tolerance and ratio interpretation floor are both 1e-12. Values below this floor remain in the raw data; ratios with a denominator below the floor are omitted.', '',
        '| T | Comparator | Signed wins / ties / losses |', '|---|---|---|']
    for t, data in by_time.items():
        for comparator, counts in data['signed_win_tie_loss'].items():
            lines.append(f'| {t} | {comparator} | {counts.get("win",0)} / {counts.get("tie",0)} / {counts.get("loss",0)} |')
    lines += ['', 'Primary settings (r=100): lower E and lower signed rank are better. Random intervals are the empirical interquartile range across twenty fixed seeds. Rank is the signed method\'s midrank among itself and those twenty samples; ties use the declared absolute tolerance.', '',
        '| Case | T | Signed E | Parent magnitude E | Block l1 E | Reverse E | Random median [Q25,Q75] | Signed rank / 21 |', '|---|---|---|---|---|---|---|---|']
    for c in comparisons:
        if c['role'] != 'primary':
            continue
        d = c['deterministic_E']
        lines.append(f'| {c["case"]} | {c["T"]:g} | {c["signed_E"]:.5g} | {d["block_parent_magnitude"]:.5g} | {d["block_l1"]:.5g} | {d["block_reverse_signed"]:.5g} | {c["random_median_E"]:.5g} [{c["random_q25_E"]:.5g}, {c["random_q75_E"]:.5g}] | {c["signed_midrank_among_signed_and_random"]:g} |')
    lines += ['', 'Convergence anchors fixed before the experiment: HCN_s1.00, Li2, and F2 at T=1 and r=50,100,200. All twenty random permutations are reused at every time and step count; no schedule is selected after observing an error.', '',
        '| Case | Ordering | E(r=50)/E(r=100) | E(r=100)/E(r=200) |', '|---|---|---|---|']
    for case in CONVERGENCE_CASES:
        g = {r: groups.get((case, 1., r)) for r in (50, 100, 200)}
        if any(value is None for value in g.values()):
            continue
        for name in BLOCK_DETERMINISTIC+EXTERNAL:
            e = [g[r][name]['one_minus_overlap'] for r in (50, 100, 200)]
            ratios = [f'{e[i]/e[i+1]:.4g}' if e[i+1] > ANALYSIS_FLOOR else 'below floor' for i in (0, 1)]
            lines.append(f'| {case} | {name} | {ratios[0]} | {ratios[1]} |')
    lines += ['', f'Maximum particle-number leakage among all whole-block treatments: {max(c["max_block_number_leakage"] for c in comparisons):.6g}; maximum spin-sector leakage: {max(c["max_block_spin_leakage"] for c in comparisons):.6g}. These numerical residuals support conservation in this calculation; they are not machine-checked exact certificates.', '',
        'Interpretation is conditional on the frozen signed-priority ownership partition and the Hartree–Fock initial state. This experiment isolates sorting after that partition has been chosen; it does not compare different ownership assignments, establish an advantage for all states, measure compiled circuit costs, or treat repeated geometries as independent families. Source snapshots, manifest hashes, complete block membership and schedules, exact-reference state hashes, norm diagnostics and leakage are saved beside this report.', '']
    (output / 'REPORT.md').write_text('\n'.join(lines))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', nargs='+', help='Default: all 11 frozen cases; subset runs are separately recorded.')
    parser.add_argument('--output-root', type=Path, default=Path(__file__).with_name('results'))
    parser.add_argument('--max-seconds', type=float, default=3600.)
    parser.add_argument('--prepare-only', action='store_true', help='Save full designs/checks without the comparison campaign.')
    args = parser.parse_args()
    if not np.isfinite(args.max_seconds) or args.max_seconds <= 0:
        parser.error('--max-seconds must be finite and positive')
    manifest = json.loads((FIXTURES / 'manifest.json').read_text())
    cases = args.cases or sorted(manifest)
    if set(cases)-manifest.keys() or len(set(cases)) != len(cases):
        parser.error('Unknown or repeated case')
    args.output_root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='controlled-blocks-', dir=args.output_root))
    sources = source_files()
    hashes = {str(path.relative_to(ROOT)): sha(path) for path in sources}
    for source in sources:
        target = output / 'source_snapshot' / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    declared = {
        'protocol_version': PROTOCOL_VERSION, 'cases': cases, 'random_seeds': RANDOM_SEEDS,
        'times': TIMES, 'primary_steps': PRIMARY_STEPS, 'convergence_cases': CONVERGENCE_CASES,
        'convergence_times': [1.], 'convergence_steps': CONVERGENCE_STEPS,
        'expected_rows': sum(len(planned_settings(case)) for case in cases)*(len(BLOCK_DETERMINISTIC)+len(RANDOM_SEEDS)+len(EXTERNAL)),
        'deterministic_block_orders': BLOCK_DETERMINISTIC, 'external_pauli_comparators': EXTERNAL,
        'ownership_frozen_once_per_case': True, 'intrablock_pauli_order_frozen': True,
        'coefficient_cutoff': ATOL, 'coefficient_policy': PRECISION_POLICY_VERSION,
        'metric': 'normalized_full_state_v2; stable one_minus_overlap; no postselection',
        'tie_absolute': TIE_ABSOLUTE, 'analysis_floor': ANALYSIS_FLOOR,
        'primary_random_summary': 'median and empirical IQR; signed midrank using fixed absolute tie tolerance',
        'product_order': 'chronological; later entries multiply on the left', 'formula_order': 1,
        'fixture_manifest_sha256': sha(FIXTURES / 'manifest.json'), 'source_hashes': hashes,
        'python_version': sys.version,
        'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'openfermion', 'numba')},
        'git_head_before_run': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'prepare_only': args.prepare_only,
    }
    write(output / 'manifest.json', declared)
    print('OUTPUT', output, flush=True)
    start, rows = time.monotonic(), []
    try:
        for case in cases:
            case_start = time.monotonic()
            H, initial, schedules, nm, sm, design, checks = prepare(case, manifest, output)
            print('PREPARED', case, design['closed_dimension'], design['nonempty_owned_blocks'], flush=True)
            if args.prepare_only:
                continue
            exacts = {}
            for total_time in TIMES:
                exact = expm_multiply(-1j*total_time*H, initial, traceA=-1j*total_time*H.diagonal().sum())
                validation = {'T': total_time, 'closed_state_sha256': state_sha(exact), 'norm': float(np.linalg.norm(exact)), 'spin_leakage': ref.norm2(exact[~sm]), 'reference': 'scipy.sparse.linalg.expm_multiply of the same identity-free retained Hamiltonian'}
                if max(abs(validation['norm']-1), validation['spin_leakage']) > 1e-9:
                    raise RuntimeError('Exact reference failed normalization or physical-sector check.')
                if len(initial) <= 256:
                    values, vectors = eigh(H.toarray())
                    independent = vectors @ (np.exp(-1j*total_time*values)*(vectors.conj().T @ initial))
                    validation['dense_eigh_state_distance'] = float(np.linalg.norm(independent-exact))
                    if validation['dense_eigh_state_distance'] > 1e-10:
                        raise RuntimeError('Dense diagonalization disagrees with sparse exact reference.')
                write(output / f'{case}.T{total_time:g}.exact.json', validation)
                exacts[total_time] = exact
            for total_time, steps, role in planned_settings(case):
                for name, terms in schedules.items():
                    if time.monotonic()-start > args.max_seconds:
                        raise RuntimeError('Predeclared wall-time budget exhausted; incomplete run is not summarized as complete.')
                    began = time.monotonic()
                    approximate = evolve(initial, terms, total_time, steps)
                    row = {
                        'case': case, 'family': family(case), 'ordering': name,
                        'T': total_time, 'steps': steps, 'role': role,
                        'pauli_terms': design['pauli_terms'], 'nominal_exponential_count': steps*design['pauli_terms'],
                        'hamiltonian_sha256': design['hamiltonian_sha256'],
                        'schedule_sha256': design['schedule_sha256'][name],
                        'fixed_block_membership_sha256': design['membership_sha256'] if name.startswith('block_') else None,
                        'evolution_seconds': time.monotonic()-began,
                        **physical_metrics(exacts[total_time], approximate, nm, sm),
                    }
                    rows.append(row)
                    write(output / f'{case}.{name}.T{total_time:g}.r{steps}.json', row)
                    print(f'ROW {len(rows)}/{declared["expected_rows"]} {case} {name} T={total_time:g} r={steps} E={row["one_minus_overlap"]:.9g}', flush=True)
            write(output / f'{case}.completion.json', {'seconds': time.monotonic()-case_start})
        for rel, digest in hashes.items():
            if sha(ROOT / rel) != digest:
                raise RuntimeError(f'Source changed during campaign: {rel}')
        if args.prepare_only:
            write(output / 'prepared.json', {'status': 'prepared_only', 'cases': cases, 'seconds': time.monotonic()-start})
        else:
            if len(rows) != declared['expected_rows']:
                raise RuntimeError('Incomplete campaign')
            write(output / 'results.json', rows)
            summarize(rows, output)
            write(output / 'completion.json', {'status': 'complete', 'rows': len(rows), 'seconds': time.monotonic()-start})
    except Exception as exc:
        write(output / 'failure.json', {'type': type(exc).__name__, 'message': str(exc), 'completed_rows': len(rows)})
        raise
    print('COMPLETE', len(rows), output, flush=True)


if __name__ == '__main__':
    main()
