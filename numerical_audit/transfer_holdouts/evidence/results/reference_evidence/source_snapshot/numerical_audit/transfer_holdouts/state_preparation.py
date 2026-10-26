"""Five prospectively declared, fully charged reference-state constructions.

Each standalone state constructor rereads/validates its input. Non-HF states
include a fresh sector eigensystem; repeated standalone preparation is actual
campaign work here, never a free shared oracle or a quantum-preparation claim.
No Trotter circuit or ordering-validation error is computed.
"""
import argparse
import re
from pathlib import Path
import shutil
import time
import transfer_common as k
import numpy as np
from cubic_response import validate


def exact_conservation(model, blocks):
    n, occupied = model['n'], model['electrons']
    k.require(type(n) is int and 1 <= n <= 6 and 0 < len(blocks) <= 256, 'Conservation domain bound')
    target = (occupied // 2, occupied // 2)
    charges = lambda v: tuple(sum((v >> q) & 1 for q in range(p, n, 2)) for p in (0, 1))
    for block in blocks:
        descriptor = k.engine.descriptor(block)
        for source in range(1 << n):
            if charges(source) != charges(source ^ descriptor[0]):
                k.require(k.engine.value(descriptor, source) == 0, 'Exact complete-unit spin-sector conservation failed')
    return [v for v in range(1 << n) if charges(v) == target]


def hamiltonian(model):
    validate(model)
    dimension = 1 << model['n']
    basis = np.arange(dimension)
    result = np.zeros((dimension, dimension), complex)
    denominator = int(model['denominator'])
    for word, c in model['coefficients']:
        descriptor = k.engine.descriptor({tuple((q, p) for q, p in word): int(c)})
        result[basis ^ descriptor[0], basis] += [k.engine.value(descriptor, int(v)) / denominator for v in basis]
    k.require(np.max(abs(result - result.conj().T)) < 1e-12, 'Full input Hamiltonian is not Hermitian')
    return result


def occupations(psi, n):
    """Explicit CAR contraction; no determinant-participation correlation proxy."""
    k.require(type(n) is int and 1 <= n <= 6 and np.shape(psi) == (1 << n,), 'RDM size guard')
    gamma = np.zeros((n, n), complex)
    for source, amplitude in enumerate(psi):
        for q in range(n):
            if not (source >> q) & 1:
                continue
            intermediate = source ^ (1 << q)
            sign = (-1) ** ((source & ((1 << q) - 1)).bit_count())
            for p in range(n):
                if not (intermediate >> p) & 1:
                    factor = sign * (-1) ** ((intermediate & ((1 << p) - 1)).bit_count())
                    gamma[p, q] += np.conj(psi[intermediate ^ (1 << p)]) * amplitude * factor
    natural = np.linalg.eigvalsh(gamma)
    defect = float(np.sum(natural * (1 - natural)))
    label = ('correlated' if defect > k.CORRELATED_THRESHOLD else
             'Slater-like' if defect <= k.SLATER_THRESHOLD else 'indeterminate')
    return gamma, natural, defect, label


def _construct(model_path, digest, label, clock):
    phases = {}
    began = clock()
    k.require(label in k.STATES, 'Undeclared state label')
    k.require(model_path.stat().st_size <= 4 * 1024**2 and k.sha(model_path) == digest, 'Input size/digest differs')
    model = k.read(model_path)
    blocks = validate(model)
    k.require(model['n'] % 2 == 0 and model['electrons'] % 2 == 0, 'Balanced closed-shell reference domain')
    sector = exact_conservation(model, blocks)
    k.require(2 <= len(sector) <= 9, 'Declared nontrivial sector bound before eigensystem allocation')
    phases['model_load_validation_CPU_seconds'] = clock() - began
    began = clock()
    dimension = 1 << model['n']
    psi = np.zeros(dimension, complex)
    metadata = dict(label=label, sector=sector, construction='HF determinant' if label == 'HF' else
                    'Fresh exact sector eigensystem plus stated state construction',
                    quantum_preparation_circuit_cost='unknown; no efficient quantum preparation protocol supplied')
    H = values = vectors = None
    if label == 'HF':
        psi[(1 << model['electrons']) - 1] = 1
        phases['Hamiltonian_CPU_seconds'] = 0.
        phases['eigensystem_CPU_seconds'] = 0.
    else:
        H = hamiltonian(model)
        phases['Hamiltonian_CPU_seconds'] = clock() - began
        began = clock()
        values, vectors = k.response.canonical_eigensystem(H[np.ix_(sector, sector)])
        phases['eigensystem_CPU_seconds'] = clock() - began
        began = clock()
        if label == 'ground':
            psi[sector] = vectors[:, 0]
        elif label == 'excited':
            psi[sector] = vectors[:, 1]
        else:
            phase = 1. if label == 'coherent_real' else 1j
            psi[sector] = (vectors[:, 0] + phase * vectors[:, 1]) / np.sqrt(2)
        psi /= np.linalg.norm(psi)
        metadata.update(energies=values.tolist(), ground_excited_gap=float(values[1] - values[0]),
                        near_degenerate=bool(values[1] - values[0] <= k.DEGENERACY_LABEL),
                        canonical_cluster_tolerance=1e-11, degeneracy_label_threshold=k.DEGENERACY_LABEL)
        if label in ('ground', 'excited'):
            index = 0 if label == 'ground' else 1
            metadata['eigen_residual'] = float(np.linalg.norm(H @ psi - values[index] * psi))
    phases['state_embedding_CPU_seconds'] = clock() - began
    began = clock()
    gamma, natural, defect, correlation = occupations(psi, model['n'])
    weights = abs(psi)**2
    k.require(abs(np.linalg.norm(psi) - 1) < 1e-12 and abs(sum(natural) - model['electrons']) < 1e-10,
              'State norm/particle count differs')
    k.require(min(natural) >= -1e-10 and max(natural) <= 1 + 1e-10, 'Invalid natural occupations')
    metadata.update(norm=float(np.linalg.norm(psi)), natural_occupations=natural.tolist(),
                    one_RDM_idempotency_defect=defect, correlation_label=correlation,
                    participation_ratio=float(1 / np.sum(weights**2)),
                    maximum_determinant_weight=float(max(weights)),
                    support=np.flatnonzero(psi != 0).tolist(),
                    correlation_rule='defect>1e-10 correlated;defect<=1e-12 Slater-like;otherwise indeterminate')
    phases['labels_RDM_CPU_seconds'] = clock() - began
    began = clock()
    del H, values, vectors, gamma, natural, weights, blocks, model
    phases['cleanup_CPU_seconds'] = clock() - began
    return psi, dict(metadata=metadata, phases=phases)


def construct(model_path, digest, label, *, clock=time.process_time):
    began, wall = clock(), time.monotonic()
    psi, record = _construct(model_path, digest, label, clock)
    elapsed = clock() - began
    record.update(CPU_seconds=elapsed, wall_seconds=time.monotonic() - wall,
                  accounting_overhead_CPU_seconds=elapsed - sum(record['phases'].values()),
                  timing='caller through full state-constructor return,diagnostics and context cleanup')
    return psi, record


def authenticate_inputs(root, preflight):
    frozen = k.read(root / 'INPUTS_FROZEN.json')
    protocol = k.read(root / 'INPUT_PROTOCOL_FROZEN.json')
    k.require(not frozen['refusals'] and frozen['protocol_sha256'] == k.sha(root / 'INPUT_PROTOCOL_FROZEN.json'),
              'Input panel refused or protocol changed')
    registry = k.read(k.HERE / 'REGISTRY.json')
    expected = [registry['chemistry_preflight']] if preflight else registry['cases']
    k.require(protocol['registry'] == registry and protocol['preflight'] == preflight and protocol['cases'] == expected,
              'Frozen prospective input registry differs')
    k.require(set(frozen['successes']) == {c['case'] for c in expected}, 'Input cohort differs')
    for relative, digest in protocol['sources'].items():
        k.require(k.sha(root / 'source_snapshot' / k.safe(relative)) == digest, 'Input producing source snapshot changed')
    for case, digest in frozen['successes'].items():
        k.require(re.fullmatch('[A-Za-z0-9_]+', case) is not None, 'Unsafe input case')
        folder = root / 'inputs' / case
        model = k.read(folder / 'model.json')
        k.require(k.sha(folder / 'model.json') == digest and model['scf_converged'], 'Hamiltonian input changed')
        for name, field in (('tensors.npz', 'tensors_sha256'), ('orbitals.npz', 'orbitals_sha256')):
            k.require(k.sha(folder / name) == model[field], 'Molecular tensor/orbital source changed')
    return frozen, protocol


def prepare(input_root, output, preflight=False):
    began, wall = time.process_time(), time.monotonic()
    k.require(not output.exists(), 'Fresh reference-preparation output required; retain partial attempts')
    commit = k.committed_sources()
    frozen_inputs, input_protocol = authenticate_inputs(input_root, preflight)
    output.mkdir(parents=True)
    shutil.copytree(input_root, output / 'inputs_evidence')
    for relative in k.sources():
        destination = output / 'source_snapshot' / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT / relative, destination)
    k.write(output / 'PROTOCOL_FROZEN.json', dict(source_commit=commit, sources=k.sources(),
            environment=k.environment(), preflight=preflight, states=k.STATES, training_time=k.TRAIN_TIME,
            validation_times=k.TIMES, input_manifest_sha256=k.sha(input_root / 'INPUTS_FROZEN.json'),
            expected_models=len(frozen_inputs['successes']), expected_states=5 * len(frozen_inputs['successes']),
            scope='Referencepreparationonly;nocircuitvalidationerrorsororderingselection',
            repeated_standalone_constructors_are_actual_work=True))
    total_constructors = 0.
    for case in frozen_inputs['successes']:
        source = output / 'inputs_evidence/inputs' / case / 'model.json'
        vectors, records = [], {}
        for label in k.STATES:
            psi, record = construct(source, frozen_inputs['successes'][case], label)
            vectors.append(psi)
            records[label] = record
            total_constructors += record['CPU_seconds']
        path = output / 'states' / (case + '.npz')
        path.parent.mkdir(exist_ok=True)
        cpu = time.process_time()
        np.savez_compressed(path, vectors=np.column_stack(vectors))
        serialization = time.process_time() - cpu
        k.write(output / 'references' / (case + '.json'), dict(case=case, states=records,
                model_sha256=frozen_inputs['successes'][case], vectors_sha256=k.sha(path),
                vectors_file=str(path.relative_to(output)), archive_serialization_CPU_seconds=serialization))
        print('REFERENCE_STATES_PREPARED', case, flush=True)
    k.write(output / 'PREPARATION_COSTS.json', dict(CPU_seconds=time.process_time() - began,
            wall_seconds=time.monotonic() - wall, actual_constructor_CPU_seconds=total_constructors,
            scope='Actualfiveconstructorspermodel,source/inputauthentication,copies,allstatearchivesandmetadata;startup/finalfreezewritesinprocessreceipt',
            input_SCF_generation_is_separate=True, quantum_preparation_circuit_cost='unknown'))
    hashes = {str(p.relative_to(output)): k.sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    k.write(output / 'REFERENCE_STATES_FROZEN.json', dict(artifacts=hashes, models=len(frozen_inputs['successes']),
            states=5 * len(frozen_inputs['successes']), scope='Allreferencevectorsandlabelsfrozenbeforeneworderingvalidationerrors'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    prepare(args.inputs.resolve(), args.output.resolve(), args.preflight)
