"""Generate the declared fourteen-mode molecular inputs before any training/error."""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
from fractions import Fraction
import shutil
import sys
import time
import study_common as k
import identity


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
    k.require(spec['n'] == 14 and spec['electrons'] == 4 and len(spec['active_spatial']) == 7,
              'Declared active space before allocation')
    frozen = spec['frozen_spatial']
    k.require(frozen == list(range(len(frozen))) and
              spec['active_spatial'] == list(range(len(frozen), len(frozen) + 7)), 'Contiguous active window')
    destination.mkdir(parents=True)
    lib.num_threads(1)
    mol = gto.M(atom=spec['geometry'], basis=registry['basis'], unit=registry['unit'],
                charge=registry['charge'], spin=registry['spin'], symmetry=registry['symmetry'],
                verbose=0, max_memory=256)
    k.require(mol.nelectron == 2 * len(frozen) + 4, 'Electron/core count')
    k.require(mol.nao_nr() == {'LiH': 11, 'BeH2': 13}[spec['family']], 'Declared basis AO count')
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
    return dict(**spec, denominator=str(denominator), coefficients=encode_integer_terms(total),
        blocks=[encode_integer_terms(b) for b in blocks], fermionic_parents=retained, policy=policy,
        removed_scalar_exact=str(scalar), basis=registry['basis'], unit=registry['unit'],
        charge=registry['charge'], spin=registry['spin'], total_electrons=int(mol.nelectron),
        hf_energy=float(energy), hf_contraction=float(contraction), orthogonality_error=ortho,
        scf_converged=True, scf_cycles=int(mean.cycles), exact_number_spin_checked=True,
        tensors_sha256=k.sha(destination / 'tensors.npz'),
        orbitals_sha256=k.sha(destination / 'orbitals.npz'))


def prepare(output):
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
    k.write_new(output / 'INPUT_PROTOCOL_FROZEN.json', dict(source_commit=head, sources=sources,
        environment=k.environment(), registry=registry, configuration=k.configuration(), geometry_identity=prior,
        scope='Input generation only; no ordering or validation errors; no new-family claim', internal_PySCF_max_memory_MB=256))


def frozen_protocol(output):
    protocol = k.read(output / 'INPUT_PROTOCOL_FROZEN.json')
    k.require(protocol['sources'] == k.source_map() and
              protocol['registry'] == k.read(k.HERE / 'REGISTRY.json'),
              'Generation source/registry changed after freeze')
    k.require(protocol['environment'] == k.environment(), 'Generation environment/basis changed after freeze')
    k.require(protocol.get('configuration') == k.configuration(), 'Frozen dispatch configuration changed')
    return protocol


def generate_case(output, name):
    protocol = frozen_protocol(output)
    specs = [case for case in protocol['registry']['cases'] if case['case'] == name]
    k.require(len(specs) == 1, 'Unknown or duplicate frozen case')
    target = output / 'inputs' / name
    k.require(not target.exists(), 'Fresh case attempt required; preserve old outputs')
    cpu, wall = time.process_time(), time.monotonic()
    try:
        model = generate_one(specs[0], protocol['registry'], target)
        model.update(generation_CPU_seconds=time.process_time() - cpu,
                     generation_wall_seconds=time.monotonic() - wall)
        k.blocks(model)  # Complete original-target/partition authentication.
        k.require(len(json.dumps(model, indent=2, allow_nan=False).encode()) <= k.MAX_MODEL_BYTES,
                  'Generated JSON model byte bound')
        k.write_new(target / 'model.json', model)
        k.write_model_toml(target / 'model.toml', model,
                           source_json_sha256=k.sha(target / 'model.json'))
    except Exception as error:
        target.mkdir(parents=True, exist_ok=True)
        k.write_new(target / 'REFUSAL.json', dict(type=type(error).__name__, message=str(error),
                    CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall))
        print('INPUT_REFUSAL', name, type(error).__name__, flush=True)
        return False
    print('GENERATED', name, 'words', len(model['coefficients']), 'units', len(model['blocks']), flush=True)
    return True


def finalize(output):
    protocol = frozen_protocol(output)
    successes, refusals = {}, {}
    for spec in protocol['registry']['cases']:
        target = output / 'inputs' / spec['case']
        if (target / 'REFUSAL.json').exists():
            refusals[spec['case']] = dict(path=str((target / 'REFUSAL.json').relative_to(output)),
                                        sha256=k.sha(target / 'REFUSAL.json'))
        elif (target / 'model.json').exists() and (target / 'model.toml').exists():
            model = k.read(target / 'model.json')
            k.require(all(model[key] == value for key, value in spec.items()), 'Frozen case identity')
            k.blocks(model)
            successes[spec['case']] = {name: k.sha(target / name) for name in ('model.json', 'model.toml', 'tensors.npz', 'orbitals.npz')}
        else:
            refusals[spec['case']] = dict(reason='Missing case output; inspect retained process receipt')
    k.write_new(output / 'INPUTS_FROZEN.json', dict(successes=successes, refusals=refusals,
        protocol_sha256=k.sha(output / 'INPUT_PROTOCOL_FROZEN.json'), expected_cases=4,
        status='complete' if not refusals else 'incomplete'))
    return not refusals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=k.Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true')
    group.add_argument('--case')
    group.add_argument('--finalize', action='store_true')
    args = parser.parse_args()
    success = True
    if args.prepare:
        prepare(args.output)
    elif args.case:
        success = generate_case(args.output, args.case)
    else:
        success = finalize(args.output)
    sys.exit(0 if success else 2)
