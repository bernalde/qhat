#!/usr/bin/env python3
"""Re-evaluate immutable molecular schedules with full-state, stable metrics.

Run: python -m qhat.numerical_audit.run_frozen_audit --cases HCN_s1.00
No molecular integrals are regenerated. Chronological order is recorded
explicitly. Every output directory is new; comparisons require exact hashes.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qhat-audit-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/qhat-audit-numba')
import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from openfermion import get_fermion_operator, jordan_wigner
from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
from qhat.analysis import benchmark_b2_active_spaces_matrix_free as evolution
from qhat.analysis import benchmark_fermionic_structure_ablation as ablation
from qhat.analysis.state_metrics import state_error_metrics
from qhat.numerical_audit import reference_tools as ref

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).with_name('fixtures')
BASE = {'fermionic_signed': 'fermionic_signed_coefficient_lexicographic',
        'jw_signed': 'signed_coefficient_lexicographic',
        'jw_magnitude': 'jw_magnitude_descending_lexicographic'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def object_sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()

def write(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')

def canonical(coefficients):
    return [[k, complex(v).real.hex(), complex(v).imag.hex()]
            for k, v in sorted(coefficients.items())]

def evolve(initial, terms, total_time, steps):
    return evolution.evolve_trotter_state(initial, terms, 1, steps, total_time, 2**30)[0]

def number_commutator_bound(keys, coefficients, bucket, parity=None):
    """Unpruned Pauli coefficient l1 bound on [G,N] on the entire Hilbert space."""
    terms = {}
    for index in bucket:
        key = keys[index]
        for at, (q, p) in enumerate(key):
            if p not in ('X', 'Y') or (parity is not None and q % 2 != parity):
                continue
            changed = key[:at] + ((q, 'Y' if p == 'X' else 'X'),) + key[at+1:]
            value = (1j if p == 'X' else -1j) * coefficients[key]
            terms[changed] = terms.get(changed, 0) + value
    return float(sum(abs(v) for v in terms.values()))

def pair_denominator(compiled):
    """For Pauli summands, ||[cP,dQ] psi||=2|cd| for each anticommuting pair."""
    xs = np.array([int(t[1]) for t in compiled], dtype=np.uint64)
    zs = np.array([int(t[2]) for t in compiled], dtype=np.uint64)
    cs = np.array([abs(t[0]) for t in compiled])
    total = 0.0
    for i, (_, x, z, _) in enumerate(compiled):
        masks = (np.uint64(x) & zs[i+1:]) ^ (np.uint64(z) & xs[i+1:])
        # No dependency on a NumPy-version-specific population-count ufunc.
        parity = np.zeros(len(masks), dtype=np.uint8)
        while np.any(masks):
            parity ^= np.asarray(masks & np.uint64(1), dtype=np.uint8)
            masks >>= np.uint64(1)
        total += 2 * cs[i] * float(np.sum(cs[i+1:][parity.astype(bool)]))
    return total

def physical_metrics(exact, approx, number_mask, spin_mask):
    metrics = state_error_metrics(exact, approx)
    if max(abs(metrics['exact_norm']-1), abs(metrics['approx_norm']-1)) > 1e-8:
        raise RuntimeError('Unacceptable state normalization drift')
    phi = exact / metrics['exact_norm']; psi = approx / metrics['approx_norm']
    residual = psi - phi * np.vdot(phi, psi)
    metrics.update({
        'particle_number_leakage': ref.norm2(psi[~number_mask]),
        'spin_sector_leakage': ref.norm2(psi[~spin_mask]),
        'inside_sector_infidelity': ref.norm2(residual[spin_mask]),
        'outside_sector_infidelity': ref.norm2(residual[~spin_mask]),
        'naive_unnormalized_one_minus_overlap': max(0., 1.-min(1., float(abs(np.vdot(exact, approx))))),
        'naive_normalized_one_minus_overlap': max(0., 1.-float(abs(np.vdot(phi, psi)))),
        'below_analysis_floor': metrics['one_minus_overlap'] <= 1e-12,
    })
    weight = ref.norm2(psi[spin_mask])
    metrics['conditional_sector_infidelity'] = (
        state_error_metrics(phi[spin_mask], psi[spin_mask])['infidelity'] if weight > 0 else None)
    return metrics

def prepare(case, manifest, out):
    fixture = FIXTURES / case
    for name, digest in manifest[case]['files'].items():
        if sha(fixture/name) != digest:
            raise RuntimeError(f'Changed fixture: {case}/{name}')
    info = json.loads((fixture/'setup.json').read_text())
    snapshot = json.loads((fixture/'operators.json').read_text())
    tensor = fixture/manifest[case]['tensor']
    interaction, n = baseline.load_interaction_operator(tensor)
    fermion = baseline.clean_fermion_operator(get_fermion_operator(interaction), 1e-12)
    jw = jordan_wigner(fermion); jw.compress(abs_tol=1e-12)
    coefficients = {k:v for k,v in jw.terms.items() if k and abs(v)>1e-12}
    digest = object_sha(canonical(coefficients))
    if digest != info['fingerprints']['identity_free_jw_sha256']:
        raise RuntimeError(f'Frozen Hamiltonian hash mismatch: {case}/{digest}')
    keys = list(coefficients)
    indices = snapshot['orders']
    for name, order in indices.items():
        if sorted(order) != list(range(len(keys))):
            raise RuntimeError(f'Not a permutation: {case}/{name}')
        digest = object_sha([[keys[i], complex(coefficients[keys[i]]).real.hex(),
                              complex(coefficients[keys[i]]).imag.hex()] for i in order])
        if digest != info['fingerprints']['physical_order_sha256'][name]:
            raise RuntimeError(f'Frozen schedule mismatch: {case}/{name}')
    regenerated = baseline.build_deterministic_orderings(fermion, keys, coefficients, n, 1e-12)
    for name, full in BASE.items():
        if regenerated[full] != [keys[i] for i in indices[name]]:
            raise RuntimeError(f'Regenerated schedule differs: {case}/{name}')
    raw = evolution.compile_ordered_terms(keys, coefficients, n, 1e-12)
    full_initial = evolution.build_hartree_fock_state(n, info['n_electrons'])
    hf = int(np.argmax(abs(full_initial)))
    states, reduced, _ = ref.affine_space(raw, hf)
    np.testing.assert_array_equal(states, snapshot['closed_basis_indices'])
    H, parity = ref.assemble_hamiltonian(reduced, len(states))
    initial = np.zeros(len(states), complex); initial[0] = 1
    nm = np.array([int(x).bit_count() == info['n_electrons'] for x in states])
    spin = set(map(int, evolution.spin_sector_basis_indices(n, info['n_electrons'])))
    sm = np.array([int(x) in spin for x in states])
    sector = np.flatnonzero(sm)
    columns = {name: [reduced[i] for i in indices[name]] for name in BASE}
    checks = {'n_qubits':n, 'n_electrons':info['n_electrons'],
              'pauli_terms':len(keys), 'closed_dimension':len(states),
              'physical_sector_dimension':len(sector),
              'input_sha256':sha(tensor), 'fingerprints':info['fingerprints'],
              'hamiltonian_and_all_schedule_hashes_match':True,
              'regenerated_baseline_orders_match':True,
              'full_space_slice_residuals':{}, 'block_commutator_bounds':{}}
    for name in BASE:
        full = evolve(full_initial, [raw[i] for i in indices[name]], .013, 1)
        small = evolve(initial, columns[name], .013, 1)
        embedded = np.zeros_like(full); embedded[states] = small
        checks['full_space_slice_residuals'][name] = float(np.linalg.norm(full-embedded))
        if checks['full_space_slice_residuals'][name] > 1e-11:
            raise RuntimeError('Closed-space evolution failed full-space validation')
    bounds = [number_commutator_bound(keys, coefficients, b) for _, b in snapshot['parent_buckets']]
    checks['block_commutator_bounds']['max_N_coefficient_l1'] = max(bounds, default=0.)
    for p in (0,1):
        checks['block_commutator_bounds'][f'max_N{p}_coefficient_l1'] = max(
            (number_commutator_bound(keys, coefficients, b, p) for _, b in snapshot['parent_buckets']), default=0.)
    checks['fallback_terms'] = len(snapshot['fallback'])
    checks['off_sector_H_frobenius'] = float(sparse.linalg.norm(H[~sm][:,sm]))
    checks['BCH_pair_norm_denominator'] = pair_denominator(raw)
    checks['denominator_same_for_all_Pauli_permutations'] = True
    # Old sector reference has a different scalar convention. Compare rays.
    exact_old, basis, _, _ = evolution.exact_reference_state(fermion,n,info['n_electrons'],.37,1e-12,True)
    exact_new = expm_multiply(-.37j*H,initial,traceA=-.37j*H.diagonal().sum())
    full_exact = np.zeros(2**n,complex);full_exact[basis] = exact_old
    checks['old_sector_reference_vs_retained_H'] = state_error_metrics(full_exact[states],exact_new)
    if checks['old_sector_reference_vs_retained_H']['phase_aligned_distance'] > 1e-8:
        raise RuntimeError('Exact reference is not phase-equivalent to retained Hamiltonian')
    if 'within_blocks_20260917' in indices:
        shuffled=evolve(initial,[reduced[i] for i in indices['within_blocks_20260917']],.1,2)
        reference=evolve(initial,columns['fermionic_signed'],.1,2)
        checks['within_block_shuffle_distance']=float(np.linalg.norm(shuffled-reference))
    for name, terms in columns.items():
        v=ref.leading_defect(terms,initial,H,parity)
        checks.setdefault('initial_BCH',{})[name]={'norm':2*float(np.linalg.norm(v)),
            'cancellation_ratio':2*float(np.linalg.norm(v))/checks['BCH_pair_norm_denominator'] if checks['BCH_pair_norm_denominator'] else None}
    write(out/f'{case}.checks.json',checks)
    return H,initial,columns,nm,sm,sector,parity,checks

def dynamic_prediction(H,initial,terms,sector,parity,total_time):
    D=ref.defect_columns(terms,sector,H,parity)
    s0=initial[sector]
    error=float(np.linalg.norm(D@s0-ref.leading_defect(terms,initial,H,parity)))
    if error>1e-9:raise RuntimeError('Leading-defect columns disagree with direct vector')
    Hs=H[sector][:,sector]
    augmented=sparse.bmat([[-1j*H,sparse.csr_matrix(D)],[None,-1j*Hs]],format='csr')
    y=expm_multiply(total_time*augmented,np.r_[np.zeros(len(initial),complex),s0],
                    traceA=total_time*augmented.diagonal().sum())
    phi=np.zeros(len(initial),complex);phi[sector]=y[len(initial):]
    phi/=np.linalg.norm(phi)
    eta=y[:len(initial)];residual=eta-phi*np.vdot(phi,eta)
    return {'integrated_projected_defect_norm2':ref.norm2(residual),
            'direct_column_residual':error}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+')
    parser.add_argument('--times',type=float,nargs='+',default=[.25,1.])
    parser.add_argument('--steps',type=int,nargs='+',default=[25,50,100,200])
    parser.add_argument('--output-root',type=Path,default=Path(__file__).with_name('results'))
    parser.add_argument('--max-seconds',type=float,default=10800)
    args=parser.parse_args()
    if any(t<=0 or not math.isfinite(t) for t in args.times) or any(r<=0 for r in args.steps):
        parser.error('Times and step counts must be positive and finite')
    fixture_manifest=json.loads((FIXTURES/'manifest.json').read_text())
    cases=args.cases or list(fixture_manifest)
    if set(cases)-fixture_manifest.keys():parser.error('Unknown frozen case')
    args.output_root.mkdir(parents=True,exist_ok=True)
    out=Path(tempfile.mkdtemp(prefix='frozen-',dir=args.output_root))
    print('OUTPUT',out,flush=True)
    paths=[ROOT/'analysis'/name for name in
           ('benchmark_L_sweep_trotter.py','benchmark_b2_active_spaces_matrix_free.py',
            'benchmark_b2_signed_coefficient_baseline.py','benchmark_b2_coloring_robustness.py',
            'benchmark_fermionic_structure_ablation.py','state_metrics.py')]
    paths += [Path(__file__),Path(ref.__file__)]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
    manifest={'git_head':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
              'source_hashes':hashes,'fixture_manifest_sha256':sha(FIXTURES/'manifest.json'),
              'python':sys.version,'versions':{p:importlib.metadata.version(p) for p in ['numpy','scipy','openfermion','numba']},
              'cases':cases,'times':args.times,'steps':args.steps,
              'product_order':'chronological; later entries multiply on the left',
              'metric':'normalized_full_state_v2; no postselection',
              'budget':'equal final Pauli set, time, initial state and r; not compiled gate-cost equality',
              'input_scope':'historical retained Hamiltonians; original OpenFermion pruning semantics preserved'}
    write(out/'manifest.json',manifest)
    start=time.monotonic();rows=[]
    try:
        for case in cases:
            if time.monotonic()-start>args.max_seconds:raise TimeoutError('Audit time budget reached')
            began=time.monotonic()
            H,initial,columns,nm,sm,sector,parity,checks=prepare(case,fixture_manifest,out)
            print('PREPARED',case,'dim',len(initial),'terms',checks['pauli_terms'],flush=True)
            old=json.loads((FIXTURES/case/'previous_results.json').read_text())
            for total_time in args.times:
                exact=expm_multiply(-1j*total_time*H,initial,traceA=-1j*total_time*H.diagonal().sum())
                for name,terms in columns.items():
                    pred=(dynamic_prediction(H,initial,terms,sector,parity,total_time)
                          if len(initial)*len(sector)<=2_000_000 else None)
                    # Save coefficient prediction before observing finite-step errors.
                    write(out/f'{case}.{name}.T{total_time:g}.prediction.json',pred)
                    for steps in args.steps:
                        if time.monotonic()-start>args.max_seconds:raise TimeoutError('Audit time budget reached')
                        approximate=evolve(initial,terms,total_time,steps)
                        row={'case':case,'ordering':name,'T':total_time,'steps':steps,
                             'nominal_exponential_count':steps*checks['pauli_terms'],
                             **physical_metrics(exact,approximate,nm,sm)}
                        if pred:
                            row['predicted_infidelity']=(total_time/steps)**2*pred['integrated_projected_defect_norm2']
                            row['prediction_relative_error']=row['predicted_infidelity']/row['infidelity']-1 if row['infidelity'] else None
                        previous=next((r for r in old if r['ordering']==name and r['T']==total_time and r['steps']==steps),None)
                        if previous:
                            row['previous_one_minus_overlap']=previous['one_minus_overlap']
                            row['previous_relative_difference']=row['one_minus_overlap']/previous['one_minus_overlap']-1 if previous['one_minus_overlap'] else None
                        rows.append(row)
                        write(out/f'{case}.{name}.T{total_time:g}.r{steps}.json',row)
                        print(f"ROW {case} {name} T={total_time:g} r={steps} E={row['one_minus_overlap']:.9g} leakage={row['particle_number_leakage']:.3g}",flush=True)
            write(out/f'{case}.completion.json',{'seconds':time.monotonic()-began,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
        for rel,digest in hashes.items():
            if sha(ROOT/rel)!=digest:raise RuntimeError('Source changed during run: '+rel)
        write(out/'results.json',rows)
        write(out/'completion.json',{'status':'complete','rows':len(rows),'seconds':time.monotonic()-start,
              'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    except Exception as exc:
        write(out/'failure.json',{'type':type(exc).__name__,'message':str(exc),'completed_rows':len(rows)})
        raise
    print('COMPLETE',len(rows),out,flush=True)

if __name__=='__main__':main()
