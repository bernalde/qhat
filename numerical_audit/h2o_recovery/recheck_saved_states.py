#!/usr/bin/env python3
"""Independently audit saved H2O pilot states; never run Hamiltonian evolution.

E is recomputed from the orthogonal residual, not the pilot's phase alignment:
E = ||a - e<e,a>||^2 / (1 + |<e,a>|), for normalized e and a.
This audit certifies completeness/integration evidence, not an accuracy benefit.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '2'
import numpy as np


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(a):
    record = {'shape': a.shape, 'dtype': a.dtype.str,
              'raw_sha256': hashlib.sha256(a.tobytes(order='C')).hexdigest()}
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    folder = args.validation.resolve()
    env = json.loads((folder / 'environment.json').read_text())
    summary = json.loads((folder / 'pilot_summary.json').read_text())
    validation = json.loads((folder / 'tensor_validation.json').read_text())
    require(sha(args.manifest) == validation['manifest_sha256'], 'Input manifest changed.')
    manifest = json.loads(args.manifest.read_text())
    require(set(manifest['bases']) == {'sto-6g', 'hgbs-5'}, 'Expected exactly the two H2O bases.')
    tensor_root = args.manifest.resolve().parent
    traces = {str(args.manifest.resolve()): sha(args.manifest)}
    for name in ('environment.json', 'pilot_summary.json', 'tensor_validation.json'):
        traces[str(folder / name)] = sha(folder / name)

    expected_tensors = {
        f'H2O_s-1.00_{basis}_as-{occupied:03d}-{vacant:03d}__{policy}'
        for basis, pairs in {
            'sto-6g': [(2, 2), (2, 4), (4, 4), (6, 4), (8, 4), (10, 4)],
            'hgbs-5': [(2, 2), (2, 4), (4, 4), (4, 6), (4, 8), (6, 8),
                       (6, 10), (8, 10), (8, 12)],
        }.items() for occupied, vacant in pairs
        for policy in ('unpruned', 'historical_policy')
    }
    expected_cases = {
        f'H2O_s-1.00_{basis}_as-002-{vacant:03d}__{policy}'
        for basis in ('sto-6g', 'hgbs-5') for vacant in (2, 4)
        for policy in ('unpruned', 'historical_policy')
    }
    require(len(validation['cases']) == 30
            and {c['case'] for c in validation['cases']} == expected_tensors,
            'Expected exactly 30 distinct validated tensor inputs.')
    for checked in validation['cases']:
        # The saved absolute tensor_path is provenance metadata only. Resolve
        # content from the supplied manifest so the artifact can be relocated.
        stem, policy = checked['case'].rsplit('__', 1)
        basis = next(b for b in ('sto-6g', 'hgbs-5') if stem.startswith(f'H2O_s-1.00_{b}_as-'))
        relative = f'{policy}/{stem}.tensors.npz'
        require(manifest['bases'][basis]['files_sha256'].get(relative) == checked['tensor_sha256'],
                'Validated tensor is not certified by the supplied generation manifest.')
        path = tensor_root / basis / relative
        require(sha(path) == checked['tensor_sha256'], 'Validated tensor changed.')
        traces[str(path)] = checked['tensor_sha256']
        require(abs(checked['hf_minus_scf_energy']) < 1e-8 + checked['hf_cutoff_bound'],
                'Recorded independent HF contraction failed.')
    for basis, entry in manifest['bases'].items():
        for relative, digest in entry['files_sha256'].items():
            path = (tensor_root / basis / relative).resolve()
            require(path.is_relative_to(tensor_root / basis), 'Generation manifest path escapes its basis directory.')
            require(sha(path) == digest, f'Generation evidence changed: {path}')
            traces[str(path.resolve())] = digest

    expected_sources = {'validate_h2o_reconstruction.py': env['script_sha256']}
    expected_sources.update({Path(record['path']).name: record['sha256']
                            for record in env['qhat_import_sources'].values()})
    expected_sources.update(env['qhat_policy_source_sha256'])
    require({p.name for p in (folder / 'source').iterdir()} == set(expected_sources),
            'Frozen pilot source set differs from recorded evidence.')
    for name, digest in expected_sources.items():
        path = folder / 'source' / name
        require(sha(path) == digest, f'Frozen source hash mismatch: {name}')
        traces[str(path)] = digest

    require(summary['all_validated_cases'] == 30 and summary['piloted_cases'] == 8
            and summary['evolutions'] == 72 and len(summary['rows']) == 72,
            'Expected 30 tensors, 8 pilot inputs, and 72 evolutions.')
    require({row['case'] for row in summary['rows']} == expected_cases,
            'Pilot case coverage differs from the declared design.')
    expected_settings = {(method, steps) for method in
                         ('fermionic_signed', 'jw_signed', 'jw_magnitude')
                         for steps in (50, 100, 200)}
    summary_rows = {(r['case'], r['ordering'], r['steps']): r for r in summary['rows']}
    require(len(summary_rows) == 72, 'Duplicate rows in pilot summary.')
    max_metric = max_norm = max_leak_difference = max_dense = 0.
    counts, states_checked = [], 0
    for case in sorted(expected_cases):
        result_path = folder / f'{case}.results.json'
        design_path = folder / f'{case}.design.json'
        result = json.loads(result_path.read_text())
        design = json.loads(design_path.read_text())
        require(result['case'] == case and design['case'] == case, 'Case identity mismatch.')
        require(len(result['rows']) == 9
                and {(r['ordering'], r['steps']) for r in result['rows']} == expected_settings,
                f'{case}: incomplete/duplicate settings.')
        require(result['state_file'] == f'{case}.states.npz',
                'Saved state filename differs from the expected local case filename.')
        state_path = (folder / result['state_file']).resolve()
        require(state_path.is_relative_to(folder), 'Saved state path escapes validation directory.')
        require(sha(state_path) == result['state_file_sha256'], 'Saved NPZ hash mismatch.')
        for path in (result_path, design_path, state_path):
            traces[str(path)] = sha(path)
        with np.load(state_path, allow_pickle=False) as saved:
            require(set(saved.files) == set(result['state_array_sha256']) and len(saved.files) == 21,
                    'Expected exactly 21 saved arrays per pilot input.')
            for name in saved.files:
                require(array_digest(saved[name]) == result['state_array_sha256'][name],
                        f'Array hash mismatch: {case}/{name}')
                require(np.all(np.isfinite(saved[name])), 'Nonfinite saved state.')
                states_checked += 1
            exact = saved['exact_retained']
            exact = exact / np.linalg.norm(exact)
            for row in result['rows']:
                require(row == summary_rows[(case, row['ordering'], row['steps'])],
                        'Case-level and summary rows disagree.')
                require(row['T'] == 1. and row['tensor_policy'] == design['tensor_policy'],
                        'Time or input policy mismatch.')
                method, steps = row['ordering'], row['steps']
                raw = saved[f'{method}_r{steps}']
                norm = np.linalg.norm(raw)
                state = raw / norm
                overlap = np.vdot(exact, state)
                perpendicular = state - exact * overlap
                E = float(np.vdot(perpendicular, perpendicular).real / (1 + abs(overlap)))
                # Floating error enters through ~eps perturbations of a state;
                # error in its squared residual scales as eps*sqrt(E)+eps^2.
                metric_bound = 100 * np.finfo(float).eps * np.sqrt(max(E, row['one_minus_overlap']))
                metric_bound += 100 * np.finfo(float).eps**2
                difference = abs(E - row['one_minus_overlap'])
                require(difference < metric_bound, 'Alternative stable overlap metric disagrees.')
                bad_number = np.array([s.bit_count() != row['n_electrons'] for s in range(state.size)])
                leakage = float(np.sum(abs(state[bad_number])**2))
                leakage_diff = abs(leakage - row['particle_number_leakage'])
                require(leakage_diff <= 1e-30 + 1e-12 * leakage, 'Leakage recomputation disagrees.')
                dense_diff = float(np.linalg.norm(raw - saved[f'{method}_r{steps}_dense']))
                require(abs(dense_diff - row['compiled_vs_dense_state_distance']) < 1e-25,
                        'Dense-vs-compiled diagnostic differs from saved vectors.')
                require(dense_diff < 1e-10 and abs(norm - 1) < 1e-9,
                        'Saved state integration or normalization failure.')
                max_metric = max(max_metric, difference)
                max_norm = max(max_norm, abs(float(norm) - 1))
                max_leak_difference = max(max_leak_difference, leakage_diff)
                max_dense = max(max_dense, dense_diff)
        counts.append({'case': case, 'settings': 9, 'saved_arrays': 21})
    report = {
        'status': 'passed', 'python': sys.version, 'numpy': np.__version__,
        'checker_path': str(Path(__file__).resolve()), 'checker_sha256': sha(__file__),
        'input_tensor_files': 30, 'pilot_inputs': 8, 'evolutions': 72,
        'saved_state_arrays_checked': states_checked,
        'max_projector_vs_phase_alignment_E_difference': max_metric,
        'max_norm_drift': max_norm,
        'max_number_leakage_recompute_difference': max_leak_difference,
        'max_compiled_vs_dense_state_distance': max_dense,
        'max_HF_contraction_residual': max(abs(c['hf_minus_scf_energy']) for c in validation['cases']),
        'scope': 'Reconstruction and small-space integration evidence only; no accuracy advantage established.',
        'generation_manifest_sha256': sha(args.manifest),
        'cases': counts, 'evidence_file_sha256': traces,
    }
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: report[k] for k in ('status', 'input_tensor_files', 'pilot_inputs',
                                           'evolutions', 'saved_state_arrays_checked',
                                           'max_projector_vs_phase_alignment_E_difference')}, indent=2))


if __name__ == '__main__':
    main()
