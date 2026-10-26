"""Authenticated new chemistry inputs; no ordering or validation error is computed."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from fractions import Fraction
import hashlib
import importlib.machinery
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def source_paths():
    paths = [HERE / n for n in ('generate_inputs.py', 'REGISTRY.json', 'PROTOCOL.md')]
    paths += [ROOT / n for n in ('analysis/fermionic_coefficients.py', 'analysis/__init__.py',
                                'numerical_audit/__init__.py', 'numerical_audit/followup_classes/exact_blocks.py')]
    return paths


def bootstrap():
    # This repository is a namespace package, and an isolated worktree need not
    # be named qhat. Bind its source root explicitly; never import another tree.
    if 'qhat' not in sys.modules:
        spec = importlib.machinery.ModuleSpec('qhat', loader=None, is_package=True)
        spec.submodule_search_locations = [str(ROOT)]
        sys.modules['qhat'] = importlib.util.module_from_spec(spec)
    from qhat.analysis import fermionic_coefficients as fc
    from qhat.numerical_audit.followup_classes import exact_blocks as eb
    require(Path(fc.__file__).resolve() == ROOT / 'analysis/fermionic_coefficients.py', 'Wrong fermion source')
    require(Path(eb.__file__).resolve() == ROOT / 'numerical_audit/followup_classes/exact_blocks.py', 'Wrong class source')
    return fc, eb


def generate_one(spec, registry, destination):
    import numpy as np
    from pyscf import gto, scf, lib
    from openfermionpyscf._run_pyscf import compute_integrals
    from openfermion.ops.representations import get_active_space_integrals
    from openfermion import InteractionOperator
    fc, eb = bootstrap()
    require(spec['n'] == 6 and spec['electrons'] in (2, 4) and len(spec['active_spatial']) == 3,
            'Declared active-space bounds before allocation')
    require(spec['frozen_spatial'] == list(range(len(spec['frozen_spatial']))) and
            spec['active_spatial'] == list(range(len(spec['frozen_spatial']), len(spec['frozen_spatial']) + 3)),
            'Fixed contiguous core/active windows required')
    destination.mkdir(parents=True)
    lib.num_threads(1)
    mol = gto.M(atom=spec['geometry'], basis=registry['basis'], unit=registry['unit'],
                charge=registry['charge'], spin=registry['spin'], symmetry=registry['symmetry'],
                verbose=0, max_memory=256)
    require(mol.nelectron == 2 * len(spec['frozen_spatial']) + spec['electrons'], 'Electron/core count differs')
    mean = scf.RHF(mol)
    mean.conv_tol = 1e-11
    mean.max_memory = 256
    mean.max_cycle = 100
    energy = mean.kernel()
    require(mean.converged, 'RHF did not converge; retain refusal')
    require(np.array_equal(mean.mo_occ, np.r_[np.full(mol.nelectron // 2, 2.),
                    np.zeros(len(mean.mo_occ) - mol.nelectron // 2)]), 'Unexpected closed-shell occupation')
    one, two = compute_integrals(mol, mean)
    require(max(spec['active_spatial']) < len(one), 'Insufficient orbitals for declared active window')
    core, h, g = get_active_space_integrals(one, two, spec['frozen_spatial'], spec['active_spatial'])
    n, occupied = spec['n'], spec['electrons']
    h1 = np.zeros((n, n))
    h2 = np.zeros((n, n, n, n))
    for sigma in (0, 1):
        h1[sigma::2, sigma::2] = h
        for tau in (0, 1):
            h2[sigma::2, tau::2, tau::2, sigma::2] = g / 2
    constant = float(mol.energy_nuc() + core)
    contraction = constant + sum(h1[p, p] for p in range(occupied)) + sum(
        h2[p, q, q, p] - h2[p, q, p, q] for p in range(occupied) for q in range(occupied))
    require(abs(contraction - energy) < 1e-8, 'Active-space HF contraction differs')
    orthogonality = float(np.max(abs(mean.mo_coeff.T @ mean.get_ovlp() @ mean.mo_coeff - np.eye(len(one)))))
    require(orthogonality < 1e-8, 'MO S-orthonormality failed')
    np.savez_compressed(destination / 'tensors.npz', constant=constant, one_body=h1, two_body=h2)
    np.savez_compressed(destination / 'orbitals.npz', mo_coeff=mean.mo_coeff, mo_energy=mean.mo_energy,
                        mo_occ=mean.mo_occ, ao_overlap=mean.get_ovlp(), ao_hcore=mean.get_hcore())
    raw = fc.fermion_from_interaction(InteractionOperator(constant, h1, h2))
    parents, total, denominator, policy = eb.exact_parents(raw, 1e-12)
    _, _, classes, masks = eb.decompose(parents, total)
    blocks = []
    for block, mask in zip(classes, masks):
        blocks.extend([{w: c} for w, c in sorted(block.items())] if mask == 0 else [block])
    require(0 < len(total) <= 256 and 0 < len(blocks) <= 256, 'Declared Pauli/unit storage bound exceeded')
    require(all(eb.internally_commuting(b) and not eb.number_commutator(b) for b in blocks),
            'Complete units do not commute internally/conserve number')
    retained = []
    for parent in parents:
        first = parent['components'][0]
        terms = [[first, str(parent['weight'])]]
        if len(parent['components']) > 1:
            partner, phase = eb.adjoint(first)
            require(partner == parent['components'][1], 'Parent reconstruction differs')
            terms.append([partner, str(phase * parent['weight'])])
        retained.append(dict(components=terms, removed_scalar=str(parent['scalar'])))
    encode = lambda block: [[word, str(c)] for word, c in sorted(block.items())]
    scalar = Fraction.from_float(constant) + Fraction(sum(p['scalar'] for p in parents), denominator)
    return dict(**spec, denominator=str(denominator), coefficients=encode(total), blocks=[encode(b) for b in blocks],
                fermionic_parents=retained, policy=policy, removed_scalar_exact=str(scalar),
                basis=registry['basis'], unit=registry['unit'], charge=registry['charge'], spin=registry['spin'],
                total_electrons=int(mol.nelectron), hf_energy=float(energy), hf_contraction=float(contraction),
                orthogonality_error=orthogonality, scf_converged=True, scf_cycles=int(mean.cycles),
                tensors_sha256=sha(destination / 'tensors.npz'), orbitals_sha256=sha(destination / 'orbitals.npz'))


def generate(output, preflight=False):
    require(not output.exists(), 'Fresh output required; retain every partial input attempt')
    commit = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    source_map = {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}
    for relative in source_map:
        require(subprocess.check_output(['git', '-C', str(ROOT), 'show', commit + ':' + relative]) ==
                (ROOT / relative).read_bytes(), 'Commit all generation sources before execution')
    registry = json.loads((HERE / 'REGISTRY.json').read_text())
    cases = [registry['chemistry_preflight']] if preflight else registry['cases']
    require(len(cases) == (1 if preflight else 6) and len({c['case'] for c in cases}) == len(cases), 'Input panel differs')
    output.mkdir(parents=True)
    env = {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'openfermion', 'openfermionpyscf', 'pyscf', 'h5py')}
    for relative in source_map:
        target = output / 'source_snapshot' / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    write(output / 'INPUT_PROTOCOL_FROZEN.json', dict(source_commit=commit, sources=source_map,
          environment=env, registry=registry, preflight=preflight, cases=cases,
          scope='Chemistrygenerationonly;nostate/order/errororcompiledcostclaim', internal_PySCF_max_memory_MiB=256))
    successes, refusals = {}, {}
    for spec in cases:
        cpu, wall = time.process_time(), time.monotonic()
        target = output / 'inputs' / spec['case']
        try:
            model = generate_one(spec, registry, target)
        except Exception as error:
            refusals[spec['case']] = dict(error_type=type(error).__name__, message=str(error),
                  CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall)
            write(target / 'REFUSAL.json', refusals[spec['case']])
            print('INPUT_REFUSAL', spec['case'], type(error).__name__, flush=True)
            continue
        model.update(generation_CPU_seconds=time.process_time() - cpu, generation_wall_seconds=time.monotonic() - wall)
        write(target / 'model.json', model)
        successes[spec['case']] = sha(target / 'model.json')
        print('GENERATED_NEW_HOLDOUT', spec['case'], 'words', len(model['coefficients']), 'units', len(model['blocks']), flush=True)
    write(output / 'INPUTS_FROZEN.json', dict(successes=successes, refusals=refusals,
          protocol_sha256=sha(output / 'INPUT_PROTOCOL_FROZEN.json'), expected_cases=len(cases)))
    require(not refusals, 'Input panel has retained refusals; do not silently substitute cases')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    generate(args.output.resolve(), args.preflight)
