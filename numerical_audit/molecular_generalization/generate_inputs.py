"""Frozen CAS(8e,6o) inputs; chemistry kernels follow the unchanged prior producer."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
from fractions import Fraction
import shutil
import sys
import time
from pathlib import Path
import common as k

def encode_integer_terms(block):
    """Canonical JSON containers also work before the JSON file is written."""
    return [[[[site, letter] for site, letter in word], str(coefficient)]
            for word, coefficient in sorted(block.items())]


def generate_one(spec, registry, destination):
    import numpy as np
    from pyscf import gto, scf, lib
    from openfermionpyscf._run_pyscf import compute_integrals
    from openfermion.ops.representations import get_active_space_integrals
    from openfermion import InteractionOperator
    k.bootstrap()
    from qhat.analysis import fermionic_coefficients as fc
    from qhat.numerical_audit.followup_classes import exact_blocks as eb
    k.require(k.Path(fc.__file__).resolve() == k.ROOT / 'analysis/fermionic_coefficients.py' and
              k.Path(eb.__file__).resolve() == k.ROOT / 'numerical_audit/followup_classes/exact_blocks.py',
              'Wrong chemistry conversion source')
    k.require(spec['n'] == 12 and spec['electrons'] == 8 and len(spec['active_spatial']) == 6,
              'Frozen CAS(8e,6o) before allocation')
    frozen = spec['frozen_spatial']
    k.require(frozen == list(range(len(frozen))) and
              spec['active_spatial'] == list(range(len(frozen), len(frozen) + 6)), 'Contiguous active window')
    destination.mkdir(parents=True)
    lib.num_threads(1)
    mol = gto.M(atom=spec['geometry'], basis=registry['basis'], unit=registry['unit'],
                charge=registry['charge'], spin=registry['spin'], symmetry=registry['symmetry'],
                verbose=0, max_memory=256)
    k.require(mol.nelectron == spec['total_electrons'] == 2 * len(frozen) + spec['electrons'], 'Electron/core count')
    k.require(mol.nao_nr() == spec['expected_AOs'], 'Declared basis AO count')
    mean = scf.RHF(mol)
    mean.conv_tol = 1e-11
    mean.max_cycle = 100
    mean.max_memory = 256
    energy = mean.kernel()
    k.require(mean.converged, 'RHF did not converge; preserve refusal')
    k.require(np.array_equal(mean.mo_occ, np.r_[np.full(mol.nelectron // 2, 2.),
              np.zeros(len(mean.mo_occ) - mol.nelectron // 2)]), 'Unexpected RHF occupation')
    one, two = compute_integrals(mol, mean)
    k.require(max(spec['active_spatial']) < len(one), 'Insufficient spatial orbitals')
    core, h, g = get_active_space_integrals(one, two, frozen, spec['active_spatial'])
    n = spec['n']
    h1, h2 = np.zeros((n, n)), np.zeros((n, n, n, n))
    for sigma in (0, 1):
        h1[sigma::2, sigma::2] = h
        for tau in (0, 1):
            h2[sigma::2, tau::2, tau::2, sigma::2] = g / 2
    constant = float(mol.energy_nuc() + core)
    occupied = spec['electrons']
    contraction = constant + sum(h1[p, p] for p in range(occupied)) + sum(
        h2[p, q, q, p] - h2[p, q, p, q] for p in range(occupied) for q in range(occupied))
    k.require(abs(contraction - energy) < 1e-8, 'Active-space HF contraction differs')
    ortho = float(np.max(abs(mean.mo_coeff.T @ mean.get_ovlp() @ mean.mo_coeff - np.eye(len(one)))))
    k.require(ortho < 1e-8, 'MO orthogonality failed')
    np.savez_compressed(destination / 'tensors.npz', constant=constant, one_body=h1, two_body=h2)
    np.savez_compressed(destination / 'orbitals.npz', mo_coeff=mean.mo_coeff, mo_energy=mean.mo_energy,
                        mo_occ=mean.mo_occ, ao_overlap=mean.get_ovlp(), ao_hcore=mean.get_hcore())
    raw = fc.fermion_from_interaction(InteractionOperator(constant, h1, h2))
    parents, total, denominator, policy = eb.exact_parents(raw, 1e-12)
    _, _, classes, masks = eb.decompose(parents, total)
    blocks = []
    for block, mask in zip(classes, masks):
        blocks.extend([{w: c} for w, c in sorted(block.items())] if mask == 0 else [block])
    k.require(0 < len(total) <= k.MAX_TERMS and 0 < len(blocks) <= k.MAX_UNITS, 'Molecular term/unit count')
    k.require(all(eb.internally_commuting(b) and not eb.number_commutator(b) for b in blocks),
              'Complete units must commute internally and conserve number')
    alpha_z = {((q, 'Z'),): 1 for q in range(0, n, 2)}
    k.require(all(not eb.commutator(alpha_z, block) for block in blocks),
              'Complete units must exactly conserve interleaved alpha number')
    retained = []
    for parent in parents:
        first = parent['components'][0]
        components = [[first, str(parent['weight'])]]
        if len(parent['components']) > 1:
            partner, phase = eb.adjoint(first)
            k.require(partner == parent['components'][1], 'Parent reconstruction differs')
            components.append([partner, str(phase * parent['weight'])])
        retained.append(dict(components=components, removed_scalar=str(parent['scalar'])))
    scalar = Fraction.from_float(constant) + Fraction(sum(p['scalar'] for p in parents), denominator)
    result = dict(spec)
    result.update(denominator=str(denominator), coefficients=encode_integer_terms(total),
        blocks=[encode_integer_terms(b) for b in blocks], fermionic_parents=retained, policy=policy,
        removed_scalar_exact=str(scalar), basis=registry['basis'], unit=registry['unit'],
        charge=registry['charge'], spin=registry['spin'], total_electrons=int(mol.nelectron),
        hf_energy=float(energy), hf_contraction=float(contraction), orthogonality_error=ortho,
        scf_converged=True, scf_cycles=int(mean.cycles), exact_number_spin_checked=True,
        tensors_sha256=k.sha(destination / 'tensors.npz'),
        orbitals_sha256=k.sha(destination / 'orbitals.npz'))
    return result


def generate_case(data, name):
    h, config = k.header(data, producer=True)
    k.require(h['mode'] == 'main' and name in h['cases'], 'Frozen main input identity')
    freeze = k.frozen(data, producer=True)
    spec = next(v for v in freeze['registry']['cases'] if v['case'] == name)
    target = data/'inputs'/name
    k.require(not target.exists(), 'One fixed input attempt; retain partial outputs')
    cpu, wall = time.process_time(), time.monotonic()
    try:
        model = generate_one(spec, freeze['registry'], target)
        model.update(generation_CPU_seconds=time.process_time()-cpu,
                     generation_wall_seconds=time.monotonic()-wall)
        k.blocks(model)
        k.require(len(json.dumps(model, allow_nan=False).encode()) <= k.MAX_MODEL_BYTES, 'Model byte bound')
        k.write_new(target/'model.json', model)
        k.write_model_toml(target/'model.toml', model, source_json_sha256=k.sha(target/'model.json'))
    except Exception as error:
        target.mkdir(parents=True, exist_ok=True)
        k.write_new(target/'REFUSAL.json', dict(type=type(error).__name__, message=str(error),
            CPU_seconds=time.process_time()-cpu, wall_seconds=time.monotonic()-wall))
        return False
    print('GENERATED', name, 'words', len(model['coefficients']), 'units', len(model['blocks']), flush=True)
    return True


def old_input(data):
    h, _ = k.header(data, producer=True)
    k.require(h['mode'] == 'preflight', 'Old input only in preflight')
    name = h['cases'][0]
    source = k.ROOT/'numerical_audit/prospective_window_ordering/evidence/run_2026-10-28/data-a6'
    index = k.read(source/'INPUTS_FROZEN.json')
    k.require(index['status'] == 'complete' and not index['refusals'], 'Authenticated retained old input')
    dest = data/'inputs'/name
    k.require(not dest.exists(), 'Fresh old input copy')
    dest.mkdir(parents=True)
    for item, digest in index['successes'][name].items():
        raw = source/'inputs'/name/item
        k.require(k.sha(raw) == digest, 'Retained old input digest')
        shutil.copyfile(raw, dest/item)
    k.blocks(k.read(dest/'model.json'))
    shutil.copyfile(source/'INPUTS_FROZEN.json', dest/'RETAINED_INPUTS_FROZEN.json')
    k.write_new(dest/'RETAINED_INPUT.json', dict(source=str(source.relative_to(k.ROOT)),
        source_index_sha256=k.sha(source/'INPUTS_FROZEN.json'), files=index['successes'][name],
        source_model_recomputed=False, scope='Old development input; no new chemistry'))


def finalize(data):
    h, _ = k.header(data, producer=True)
    successes, refusals = {}, {}
    for name in h['cases']:
        dest = data/'inputs'/name
        if (dest/'REFUSAL.json').exists():
            refusals[name] = dict(path=str((dest/'REFUSAL.json').relative_to(data)), sha256=k.sha(dest/'REFUSAL.json'))
        elif all((dest/p).is_file() for p in ('model.json', 'model.toml', 'orbitals.npz', 'tensors.npz')):
            model = k.read(dest/'model.json')
            k.require(model['case'] == name, 'Generated case identity')
            if h['mode'] == 'main':
                spec = next(v for v in k.frozen(data)['registry']['cases'] if v['case'] == name)
                k.require(all(model[key] == value for key, value in spec.items()), 'Frozen case specification')
            k.blocks(model)
            successes[name] = {p: k.sha(dest/p) for p in ('model.json', 'model.toml', 'orbitals.npz', 'tensors.npz')}
        else:
            refusals[name] = dict(reason='Missing output; retain owning receipt and partial files')
    k.write_new(data/'INPUTS_FROZEN.json', dict(h, successes=successes, refusals=refusals,
        status='complete' if not refusals else 'incomplete'))
    return not refusals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--freeze', choices=('preflight', 'main'))
    group.add_argument('--old-input', action='store_true')
    group.add_argument('--case')
    group.add_argument('--finalize', action='store_true')
    args = parser.parse_args()
    success = True
    if args.freeze:
        k.freeze(args.data, mode=args.freeze)
    elif args.old_input:
        old_input(args.data)
    elif args.case:
        success = generate_case(args.data, args.case)
    else:
        success = finalize(args.data)
    sys.exit(0 if success else 2)
