#!/usr/bin/env python3
"""Validate regenerated H2O tensors and run a bounded small-space integration pilot.

This is a reconstruction, not a replay of the unavailable historical tensors.
No fixed-particle-number projection is used in the evolution or error metric.
The independent reference implements occupation-basis CAR directly from arrays;
Pauli matrices use explicit Kronecker products, not OpenFermion's JW transform.
"""
from __future__ import annotations

import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '2'
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/qhat-h2o-numba')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-h2o-mpl')

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import resource
import shutil
import sys
import time

# Bound virtual memory before loading scientific libraries.
_cap = 12 * 1024**3
_soft, _hard = resource.getrlimit(resource.RLIMIT_AS)
resource.setrlimit(resource.RLIMIT_AS, (min(_cap, _hard) if _hard >= 0 else _cap, _hard))

import numpy as np
from scipy.linalg import expm
from openfermion import InteractionOperator

from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis.fermionic_coefficients import PRECISION_POLICY_VERSION, fermion_from_interaction
from qhat.analysis.fermionic_ownership import build_molecular_ownership
from qhat.analysis.state_metrics import state_error_metrics

ATOL = 1e-12
METHODS = ('fermionic_signed', 'jw_signed', 'jw_magnitude')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def array_sha(value):
    a = np.asarray(value)
    return object_sha({'shape': a.shape, 'dtype': a.dtype.str,
                       'raw_sha256': hashlib.sha256(a.tobytes(order='C')).hexdigest()})


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def hf_contraction(constant, one, two, electrons):
    """<HF|H|HF> for t[p,q,r,s] a_p^+ a_q^+ a_r a_s, no extra 1/2."""
    value = complex(constant)
    value += sum(one[p, p] for p in range(electrons))
    value += sum(two[p, q, q, p] - two[p, q, p, q]
                 for p in range(electrons) for q in range(electrons))
    return value


def occupation_action(source, word, n):
    """Apply rightmost operator first; q0 is the most significant bit."""
    state, sign = int(source), 1
    for mode, create in reversed(word):
        bit = 1 << (n - 1 - mode)
        if bool(state & bit) == bool(create):
            return None, 0
        if (state >> (n - mode)).bit_count() % 2:
            sign = -sign
        state ^= bit
    return state, sign


def dense_car(constant, one, two):
    n = len(one)
    result = complex(constant) * np.eye(2**n, dtype=complex)
    for indices, word_pattern, tensor in (
        (np.argwhere(one != 0), (1, 0), one),
        (np.argwhere(two != 0), (1, 1, 0, 0), two),
    ):
        for row in indices:
            key = tuple(map(int, row))
            word = tuple(zip(key, word_pattern))
            coefficient = tensor[key]
            for source in range(2**n):
                destination, phase = occupation_action(source, word, n)
                if phase:
                    result[destination, source] += coefficient * phase
    return result


def dense_pauli(key, n):
    matrices = {
        'I': np.eye(2, dtype=complex),
        'X': np.array([[0, 1], [1, 0]], dtype=complex),
        'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
        'Z': np.diag([1., -1.]).astype(complex),
    }
    word = dict(key)
    result = np.ones((1, 1), dtype=complex)
    for mode in range(n):
        result = np.kron(result, matrices[word.get(mode, 'I')])
    return result


def independent_metrics(exact, approx, n, electrons):
    norm_exact = float(np.linalg.norm(exact))
    norm_approx = float(np.linalg.norm(approx))
    require(abs(norm_exact - 1) < 1e-9 and abs(norm_approx - 1) < 1e-9,
            'Evolution normalization drift exceeds 1e-9.')
    e, a = exact / norm_exact, approx / norm_approx
    z = np.vdot(e, a)
    aligned = a * (np.conjugate(z) / abs(z)) if z != 0 else a
    E = float(np.vdot(e - aligned, e - aligned).real / 2)
    number_mask = np.array([s.bit_count() == electrons for s in range(2**n)])
    alpha_mask = sum(1 << (n - 1 - mode) for mode in range(0, n, 2))
    beta_mask = sum(1 << (n - 1 - mode) for mode in range(1, n, 2))
    spin_mask = np.array([(s & alpha_mask).bit_count() == electrons // 2
                          and (s & beta_mask).bit_count() == electrons // 2
                          for s in range(2**n)])
    result = {
        'one_minus_overlap': E, 'infidelity': E * (2 - E),
        'phase_aligned_distance': float(np.sqrt(2 * E)),
        'exact_norm': norm_exact, 'approx_norm': norm_approx,
        'particle_number_leakage': float(np.vdot(a[~number_mask], a[~number_mask]).real),
        'spin_sector_leakage': float(np.vdot(a[~spin_mask], a[~spin_mask]).real),
    }
    common = state_error_metrics(exact, approx)
    metric_difference = abs(common['phase_aligned_distance'] - result['phase_aligned_distance'])
    require(metric_difference <= 50 * np.finfo(float).eps * max(1., result['phase_aligned_distance']),
            'Independent and existing stable metrics disagree.')
    result['stable_metric_phase_distance_disagreement'] = metric_difference
    return result


def validate_tensor(case, tensor_path, electrons, hf_energy=None, hf_cutoff_bound=0.):
    with np.load(tensor_path, allow_pickle=False) as raw:
        require(set(raw.files) == {'constant', 'one_body', 'two_body'},
                f'{case}: unexpected tensor keys {raw.files}')
        constant = raw['constant'].item()
        one, two = raw['one_body'].copy(), raw['two_body'].copy()
    n = one.shape[0]
    require(one.shape == (n, n) and two.shape == (n, n, n, n),
            f'{case}: tensor shape mismatch')
    require(0 < electrons < n and electrons % 2 == 0 and n % 2 == 0,
            f'{case}: invalid restricted closed-shell active space')
    require(all(np.all(np.isfinite(a)) for a in (constant, one, two)),
            f'{case}: nonfinite tensor')
    imaginary = max(float(np.max(abs(np.imag(a)))) for a in (constant, one, two))
    one_res = float(np.max(abs(one - one.T.conj())))
    two_res = float(np.max(abs(two - two.transpose(3, 2, 1, 0).conj())))
    require(imaginary < 1e-12 and max(one_res, two_res) < 1e-10,
            f'{case}: real/Hermitian tensor check failed')
    hf = hf_contraction(constant, one, two, electrons)
    require(abs(hf.imag) < 1e-12, f'{case}: nonreal HF contraction')
    difference = None if hf_energy is None else float(hf.real - hf_energy)
    if difference is not None:
        require(abs(difference) < 1e-8 + hf_cutoff_bound,
                f'{case}: HF contraction differs from SCF beyond cutoff bound: {difference}')
    checks = {
        'case': case, 'tensor_path': str(tensor_path), 'tensor_sha256': sha(tensor_path),
        'n_qubits': n, 'n_electrons': electrons,
        'constant': float(np.real(constant)),
        'array_sha256': {name: array_sha(a) for name, a in
                         [('constant', np.asarray(constant)), ('one_body', one), ('two_body', two)]},
        'one_body_hermiticity_max': one_res, 'two_body_hermiticity_max': two_res,
        'tensor_imaginary_max': imaginary, 'hf_energy_from_occupation_contraction': float(hf.real),
        'scf_energy': hf_energy, 'hf_minus_scf_energy': difference,
        'hf_cutoff_bound': hf_cutoff_bound,
    }
    return constant, one, two, checks


def run_case(case, constant, one, two, electrons, output, steps=(50, 100, 200), tensor_policy='synthetic'):
    start = time.monotonic()
    n, dim = len(one), 2**len(one)
    require(n <= 6, 'Dense independent pilot is limited to <=6 qubits.')
    interaction = InteractionOperator(constant, one, two)
    owned = build_molecular_ownership(fermion_from_interaction(interaction), atol=ATOL,
                                      ordering='signed_ascending')
    require(not owned.fallback_pauli_keys and
            not owned.diagnostics['noncommuting_owned_block_count'],
            f'{case}: ownership is not an exhaustive commuting decomposition')
    coeff = owned.final_coefficients
    require(all(complex(c).imag == 0 for c in coeff.values()),
            f'{case}: retained Pauli coefficients are not exactly real')
    keys = list(coeff)
    schedules = {
        'fermionic_signed': list(owned.pauli_order),
        'jw_signed': baseline.signed_coefficient_lexicographic_order(keys, coeff, n, ATOL),
        'jw_magnitude': baseline.magnitude_descending_lexicographic_order(keys, coeff, n, ATOL),
    }
    require(all(len(order) == len(keys) and set(order) == set(keys)
                for order in schedules.values()), 'Schedules change the common Hamiltonian.')
    pauli = {key: dense_pauli(key, n) for key in keys}
    H = sum((coeff[key] * pauli[key] for key in keys), np.zeros((dim, dim), complex))
    car = dense_car(constant, one, two)
    residual = float(np.linalg.norm(car - owned.scalar * np.eye(dim) - H, 2))
    cutoff_bound = owned.diagnostics['discarded_pauli_coefficient_l1']
    require(residual <= cutoff_bound + 1e-10, 'CAR reference differs beyond final-cutoff bound.')
    car_herm = float(np.linalg.norm(car - car.conj().T, 2))
    require(car_herm < 1e-10, 'Independent CAR Hamiltonian is not Hermitian.')
    initial = np.zeros(dim, complex)
    initial[((1 << electrons) - 1) << (n - electrons)] = 1
    require(np.array_equal(initial, evolution.build_hartree_fock_state(n, electrons)),
            'HF index convention mismatch.')
    N = np.diag([index.bit_count() for index in range(dim)])
    block_comm = 0.
    for block in owned.owned_blocks:
        G = sum((coeff[k] * pauli[k] for k in block.pauli_keys), np.zeros_like(H))
        block_comm = max(block_comm, float(np.linalg.norm(N @ G - G @ N, 2)))
    require(block_comm < 1e-10, 'Owned block fails dense particle-number check.')
    exact = expm(-1j * H) @ initial
    exact_car = expm(-1j * car) @ initial
    states = {'initial': initial, 'exact_retained': exact, 'exact_input_tensor_car': exact_car}
    checks = {
        'case': case, 'n_qubits': n, 'n_electrons': electrons, 'tensor_policy': tensor_policy,
        'coefficient_cutoff': ATOL, 'coefficient_policy': PRECISION_POLICY_VERSION,
        'pauli_terms': len(keys), 'owned_blocks': len(owned.owned_blocks),
        'discarded_scalar': [complex(owned.scalar).real, complex(owned.scalar).imag],
        'input_tensor_car_hermiticity_operator_norm': car_herm,
        'car_minus_retained_and_scalar_operator_norm': residual,
        'final_cutoff_operator_norm_bound': cutoff_bound,
        'max_owned_block_number_commutator_operator_norm': block_comm,
        'retained_vs_input_tensor_exact_metrics': independent_metrics(exact_car, exact, n, electrons),
        'ownership_diagnostics': owned.diagnostics,
        'hamiltonian_sha256': object_sha([[key, complex(coeff[key]).real.hex(),
                                          complex(coeff[key]).imag.hex()] for key in sorted(keys)]),
        'pauli_orders': {name: [[key, complex(coeff[key]).real.hex()] for key in order]
                         for name, order in schedules.items()},
        'order_convention': 'chronological left-to-right list; first listed factor acts first',
        'state_convention': 'full Hilbert space, q0 most significant; HF first electrons occupied',
    }
    write(output / f'{case}.design.json', checks)
    rows = []
    identity = np.eye(dim, dtype=complex)
    for method, order in schedules.items():
        compiled = evolution.compile_ordered_terms(order, coeff, n, ATOL)
        for r in steps:
            row_start = time.monotonic()
            # Independent dense Pauli rotations; identity P^2=I is exact.
            one_step = identity.copy()
            for key in order:
                angle = float(complex(coeff[key]).real) / r
                one_step = (np.cos(angle) * identity - 1j * np.sin(angle) * pauli[key]) @ one_step
            dense = np.linalg.matrix_power(one_step, r) @ initial
            actual, count = evolution.evolve_trotter_state(initial, compiled, 1, r, 1., 2**30)
            absolute_difference = float(np.linalg.norm(dense - actual))
            require(absolute_difference < 1e-10, 'Dense rotations and compiled evolution differ.')
            require(count == r * len(keys), 'Unexpected Pauli exponential count.')
            require(float(np.linalg.norm(one_step.conj().T @ one_step - identity, 2)) < 1e-11,
                    'Dense product step is not unitary.')
            row = {
                'case': case, 'tensor_policy': tensor_policy, 'ordering': method, 'T': 1., 'steps': r,
                'n_qubits': n, 'n_electrons': electrons, 'nominal_exponentials': count,
                'compiled_vs_dense_state_distance': absolute_difference,
                'full_state_sha256': array_sha(actual),
                **independent_metrics(exact, actual, n, electrons),
                'vs_input_tensor_car': independent_metrics(exact_car, actual, n, electrons),
                'elapsed_seconds': time.monotonic() - row_start,
            }
            states[f'{method}_r{r}'] = actual
            states[f'{method}_r{r}_dense'] = dense
            rows.append(row)
    state_file = output / f'{case}.states.npz'
    np.savez_compressed(state_file, **states)
    result = {'case': case, 'rows': rows, 'elapsed_seconds': time.monotonic() - start,
              'state_file': state_file.name, 'state_file_sha256': sha(state_file),
              'state_array_sha256': {name: array_sha(value) for name, value in states.items()}}
    write(output / f'{case}.results.json', result)
    return result


def self_test(output):
    n = 4
    one = np.diag([-.9, -.9, .2, .2])
    one[0, 2] = one[2, 0] = .13
    one[1, 3] = one[3, 1] = .13
    two = np.zeros((n, n, n, n))
    two[0, 1, 1, 0] = .4
    two[1, 0, 0, 1] = .4
    constant = .7
    h = dense_car(constant, one, two)
    require(abs(h[12, 12] - hf_contraction(constant, one, two, 2)) < 1e-14,
            'Synthetic HF contraction failed.')
    for p in range(n):
        for q in range(n):
            pair = np.zeros_like(h)
            for word in (((p, 0), (q, 1)), ((q, 1), (p, 0))):
                for source in range(2**n):
                    destination, sign = occupation_action(source, word, n)
                    if sign:
                        pair[destination, source] += sign
            require(np.array_equal(pair, np.eye(2**n) if p == q else np.zeros_like(h)),
                    'Synthetic CAR anticommutator check failed.')
    return run_case('synthetic_four_mode', constant, one, two, 2, output, steps=(2, 10))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--qhat-root', type=Path, default=Path('/home/bernalde/repos/qhat'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    qhat_root = args.qhat_root.resolve()
    require(all(Path(module.__file__).resolve().is_relative_to(qhat_root)
                for module in (evolution, baseline)), 'Imported QHAT outside requested checkout.')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'source').mkdir()
    source_paths = [Path(__file__), Path(evolution.__file__), Path(baseline.__file__)]
    source_paths += [Path(baseline.__file__).parent / name for name in
                     ('fermionic_coefficients.py', 'fermionic_ownership.py',
                      'state_metrics.py', 'benchmark_L_sweep_trotter.py')]
    for path in source_paths:
        shutil.copyfile(path, args.output / 'source' / path.name)
    env = {
        'python': sys.version, 'executable': sys.executable,
        'script_sha256': sha(__file__),
        'memory_limit_bytes': resource.getrlimit(resource.RLIMIT_AS)[0],
        'thread_environment': {key: os.environ[key] for key in
                              ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS')},
        'packages': {name: importlib.metadata.version(name) for name in
                     ('numpy', 'scipy', 'openfermion', 'numba', 'pyscf', 'openfermionpyscf')},
        'qhat_import_sources': {name: {'path': str(module.__file__), 'sha256': sha(module.__file__)}
                               for name, module in [('evolution', evolution), ('baseline', baseline)]},
        'qhat_policy_source_sha256': {name: sha(Path(baseline.__file__).parent / name)
                                     for name in ('fermionic_coefficients.py', 'fermionic_ownership.py',
                                                  'state_metrics.py', 'benchmark_L_sweep_trotter.py')},
        'interpretation': 'Reconstructed inputs and bounded integration validation, not historical replication.',
    }
    write(args.output / 'environment.json', env)
    if args.self_test:
        self_test(args.output)
        print(json.dumps({'self_test': 'passed', 'output': str(args.output)}), flush=True)
        return
    require(args.manifest is not None, '--manifest required except for --self-test.')
    manifest = json.loads(args.manifest.read_text())
    require(set(manifest['bases']) == {'sto-6g', 'hgbs-5'},
            'Both completed reconstructed bases are required.')
    entries = []
    for basis, record in manifest['bases'].items():
        basis_dir = args.manifest.parent / basis
        for relative, digest in record['files_sha256'].items():
            require(sha(basis_dir / relative) == digest,
                    f'Changed reconstruction file: {basis}/{relative}')
        summary = json.loads((basis_dir / 'summary.json').read_text())
        require(summary['scf_converged'], f'{basis}: SCF did not converge.')
        for info in summary['cases']:
            for policy in ('unpruned', 'historical_policy'):
                relative = f'{policy}/{info["case_id"]}.tensors.npz'
                entries.append({
                    'case': info['case_id'] + '__' + policy,
                    'tensor_path': str(basis_dir / relative),
                    'n_electrons': info['active_electrons'],
                    'n_qubits': info['active_spin_orbitals'],
                    'hf_energy': summary['hf_energy_hartree'], 'policy': policy,
                    'tensor_sha256': record['files_sha256'][relative],
                    'hf_cutoff_bound': (info['cutoff_coefficient_l1_difference']
                                        if policy == 'historical_policy' else 0.),
                })
    require(len(entries) == 30, 'Expected 15 active spaces under each of two tensor policies.')
    validation, results = [], []
    for entry in entries:
        tensor = Path(entry['tensor_path'])
        if not tensor.is_absolute():
            tensor = args.manifest.parent / tensor
        case = entry['case']
        electrons = int(entry['n_electrons'])
        constant, one, two, checked = validate_tensor(case, tensor, electrons,
                                                     entry['hf_energy'], entry['hf_cutoff_bound'])
        checked['policy'] = entry['policy']
        require(len(one) == entry['n_qubits'], f'{case}: manifest qubit mismatch')
        if 'tensor_sha256' in entry:
            require(checked['tensor_sha256'] == entry['tensor_sha256'], f'{case}: tensor hash mismatch')
        validation.append(checked)
        if not args.validate_only and len(one) <= 6:
            results.append(run_case(case, constant, one, two, electrons, args.output,
                                    tensor_policy=entry['policy']))
        print(json.dumps({'validated': case, 'qubits': len(one),
                          'pilot': len(one) <= 6 and not args.validate_only}), flush=True)
    write(args.output / 'tensor_validation.json', {
        'manifest_sha256': sha(args.manifest), 'cases': validation,
        'independent_HF_contraction_matches_SCF_with_declared_cutoff_bound': all(
            v['hf_minus_scf_energy'] is not None
            and abs(v['hf_minus_scf_energy']) < 1e-8 + v['hf_cutoff_bound'] for v in validation),
    })
    write(args.output / 'pilot_summary.json', {
        'interpretation': env['interpretation'],
        'all_validated_cases': len(validation), 'piloted_cases': len(results),
        'evolutions': sum(len(r['rows']) for r in results),
        'rows': [row for result in results for row in result['rows']],
        'limit': 'Small-space integration checks do not establish a benefit across the full H2O sweep.',
    })


if __name__ == '__main__':
    main()
