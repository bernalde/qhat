"""Generate the frozen twelve-mode molecular inputs before any training/error."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from fractions import Fraction
import shutil
import sys
import time
import common as k
import identity


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
    k.require(spec['n'] == 12 and spec['electrons'] == 6 and len(spec['active_spatial']) == 6,
              'Declared active space before allocation')
    frozen = spec['frozen_spatial']
    k.require(frozen == list(range(len(frozen))) and
              spec['active_spatial'] == list(range(len(frozen), len(frozen) + 6)), 'Contiguous active window')
    destination.mkdir(parents=True)
    lib.num_threads(1)
    mol = gto.M(atom=spec['geometry'], basis=registry['basis'], unit=registry['unit'],
                charge=registry['charge'], spin=registry['spin'], symmetry=registry['symmetry'],
                verbose=0, max_memory=256)
    k.require(mol.nelectron == 2 * len(frozen) + 6, 'Electron/core count')
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
    retained = []
    for parent in parents:
        first = parent['components'][0]
        components = [[first, str(parent['weight'])]]
        if len(parent['components']) > 1:
            partner, phase = eb.adjoint(first)
            k.require(partner == parent['components'][1], 'Parent reconstruction differs')
            components.append([partner, str(phase * parent['weight'])])
        retained.append(dict(components=components, removed_scalar=str(parent['scalar'])))
    encode = lambda block: [[word, str(c)] for word, c in sorted(block.items())]
    scalar = Fraction.from_float(constant) + Fraction(sum(p['scalar'] for p in parents), denominator)
    return dict(**spec, denominator=str(denominator), coefficients=encode(total),
        blocks=[encode(b) for b in blocks], fermionic_parents=retained, policy=policy,
        removed_scalar_exact=str(scalar), basis=registry['basis'], unit=registry['unit'],
        charge=registry['charge'], spin=registry['spin'], total_electrons=int(mol.nelectron),
        hf_energy=float(energy), hf_contraction=float(contraction), orthogonality_error=ortho,
        scf_converged=True, scf_cycles=int(mean.cycles), tensors_sha256=k.sha(destination / 'tensors.npz'),
        orbitals_sha256=k.sha(destination / 'orbitals.npz'))


def generate(output):
    k.require(not output.exists(), 'Fresh output required; retain partial attempts')
    head = k.committed_sources()
    registry = k.read(k.HERE / 'REGISTRY.json')
    k.require(len(registry['cases']) == 4 and len({c['case'] for c in registry['cases']}) == 4,
              'Frozen four-case registry')
    prior = identity.audit(registry)
    output.mkdir(parents=True)
    sources = k.source_map()
    for name in sources:
        target = output / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT / name, target)
    k.write(output / 'INPUT_PROTOCOL_FROZEN.json', dict(source_commit=head, sources=sources,
        environment=k.environment(), registry=registry, geometry_identity=prior,
        scope='Input generation only; no ordering or validation errors', internal_PySCF_memory_MiB=256))
    successes, refusals = {}, {}
    for spec in registry['cases']:
        cpu, wall = time.process_time(), time.monotonic()
        target = output / 'inputs' / spec['case']
        try:
            model = generate_one(spec, registry, target)
        except Exception as error:
            refusals[spec['case']] = dict(type=type(error).__name__, message=str(error),
                CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall)
            k.write(target / 'REFUSAL.json', refusals[spec['case']])
            print('INPUT_REFUSAL', spec['case'], type(error).__name__, flush=True)
            continue
        model.update(generation_CPU_seconds=time.process_time() - cpu,
                     generation_wall_seconds=time.monotonic() - wall)
        k.write(target / 'model.json', model)
        successes[spec['case']] = k.sha(target / 'model.json')
        print('GENERATED', spec['case'], 'words', len(model['coefficients']), 'units', len(model['blocks']), flush=True)
    k.write(output / 'INPUTS_FROZEN.json', dict(successes=successes, refusals=refusals,
            protocol_sha256=k.sha(output / 'INPUT_PROTOCOL_FROZEN.json'), expected_cases=4))
    k.require(not refusals, 'Complete panel contains retained chemistry refusals')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=k.Path, required=True)
    generate(parser.parse_args().output)
