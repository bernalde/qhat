"""Independent literal Pauli/CAR/tensor and full reference-state verification."""
import argparse
from fractions import Fraction
from pathlib import Path
import time
import transfer_common as k
import numpy as np
from scipy import sparse
import geometry_identity


def literal(word, n):
    k.require(type(n) is int and 1 <= n <= 6, 'Independent register bound')
    k.require(len(dict(word)) == len(word) and all(type(q) is int and 0 <= q < n and p in 'XYZ' for q, p in word), 'Literal support')
    letters = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
               'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    result = sparse.csr_matrix([[1. + 0j]])
    labels = dict(word)
    for q in reversed(range(n)):
        result = sparse.kron(result, letters[labels.get(q, 'I')], format='csr')
    return result


def annihilation(q, n):
    k.require(type(n) is int and 1 <= n <= 6 and type(q) is int and 0 <= q < n, 'Independent CAR bound')
    result = sparse.csr_matrix([[1. + 0j]])
    for p in reversed(range(n)):
        factor = (np.array([[0, 1], [0, 0]]) if p == q else
                  np.diag([1, -1]) if p < q else np.eye(2))
        result = sparse.kron(result, factor, format='csr')
    return result


def car(word, pool, dimension):
    result = sparse.eye(dimension, dtype=complex, format='csr')
    for q, create in word:
        result = result @ (pool[q].conj().T if create else pool[q])
    return result


def tensor_hamiltonian(one, two, constant):
    n = len(one)
    k.require(1 <= n <= 6 and one.shape == (n, n) and two.shape == (n, n, n, n), 'Tensor bound before CAR allocation')
    k.require(np.all(np.isfinite(one)) and np.all(np.isfinite(two)) and np.isfinite(constant), 'Nonfinite tensors')
    pool = [annihilation(q, n) for q in range(n)]
    dimension = 1 << n
    result = constant * sparse.eye(dimension, dtype=complex, format='csr')
    for p, q in np.argwhere(one != 0):
        result += one[p, q] * (pool[p].conj().T @ pool[q])
    for p, q, r, s in np.argwhere(two != 0):
        result += two[p, q, r, s] * car(((p, 1), (q, 1), (r, 0), (s, 0)), pool, dimension)
    return result.toarray(), pool


def molecular_identity(folder, model, H):
    n, occupied = model['n'], model['electrons']
    den = int(model['denominator'])
    for filename, key in (('tensors.npz', 'tensors_sha256'), ('orbitals.npz', 'orbitals_sha256')):
        k.require(k.sha(folder / filename) == model[key], 'Molecular source archive differs')
    with np.load(folder / 'tensors.npz', allow_pickle=False) as z:
        h, g, constant = z['one_body'], z['two_body'], float(z['constant'])
        raw, pool = tensor_hamiltonian(h, g, constant)
    policy = model['policy']
    allowance = policy['canonical_hermitian_adjustment_fermionic_l1'] + policy['discarded_parent_fermionic_l1']
    k.require(constant == policy['input_constant'] and policy['parent_cutoff'] == 1e-12 and 0 <= allowance <= 2e-8, 'Target-policy convention differs')
    scalar = constant + policy['removed_parent_pauli_scalar']
    difference = float(np.max(abs(raw - H - scalar * np.eye(1 << n))))
    k.require(difference <= allowance + 2e-10, 'Saved tensors do not produce the supplied target within declared policy')
    k.require(abs(float(Fraction(model['removed_scalar_exact'])) - scalar) < 2e-12, 'Removed scalar convention differs')
    hf_difference = float(abs(raw[(1 << occupied) - 1, (1 << occupied) - 1] - model['hf_energy']))
    k.require(hf_difference < 1e-8, 'Independent tensor/CAR HF contraction differs')
    parents = sparse.csr_matrix(H.shape, dtype=complex)
    for parent in model['fermionic_parents']:
        for word, c in parent['components']:
            parents += int(c) / den * car(word, pool, len(H))
        parents -= int(parent['removed_scalar']) / den * sparse.eye(len(H), dtype=complex, format='csr')
    parent_difference = float(np.max(abs(parents.toarray() - H)))
    k.require(parent_difference < 2e-11, 'Retained exact parent/Pauli identity differs')
    with np.load(folder / 'orbitals.npz', allow_pickle=False) as z:
        C, S, occ = z['mo_coeff'], z['ao_overlap'], z['mo_occ']
        orthogonality = float(np.max(abs(C.T @ S @ C - np.eye(C.shape[1]))))
        k.require(orthogonality < 1e-8 and sum(occ) == model['total_electrons'], 'Saved orbital identity differs')
        occupied_spatial = model['total_electrons'] // 2
        k.require(np.array_equal(occ, np.r_[np.full(occupied_spatial, 2.), np.zeros(len(occ) - occupied_spatial)]), 'Saved RHF occupations differ')
    return dict(raw_tensor_target_difference=difference, policy_L1_allowance=allowance,
                retained_parent_difference=parent_difference, HF_CAR_difference=hf_difference,
                orbital_orthogonality_error=orthogonality), pool


def verify(output, receipt_dir=None):
    receipt = output if receipt_dir is None else receipt_dir
    k.require(not (receipt / 'REFERENCE_VERIFIED.json').exists(), 'Preserve verification receipts')
    if receipt_dir is not None:
        k.require(not receipt.exists(), 'Fresh receipt directory required')
        receipt.mkdir(parents=True)
    began, wall = time.process_time(), time.monotonic()
    protocol = k.read(output / 'PROTOCOL_FROZEN.json')
    frozen = k.read(output / 'REFERENCE_STATES_FROZEN.json')
    k.require(protocol['states'] == list(k.STATES) and protocol['training_time'] == k.TRAIN_TIME and
              protocol['validation_times'] == list(k.TIMES), 'State/time registry changed')
    for name, digest in frozen['artifacts'].items():
        k.require(k.sha(output / k.safe(name)) == digest, 'Frozen reference artifact changed')
    source_files = {str(p.relative_to(output / 'source_snapshot')) for p in (output / 'source_snapshot').rglob('*') if p.is_file()}
    k.require(source_files == set(protocol['sources']), 'Reference producer source inventory differs')
    for name, digest in protocol['sources'].items():
        k.require(k.sha(output / 'source_snapshot' / k.safe(name)) == digest, 'Reference producing source changed')
    inputs = output / 'inputs_evidence'
    original = k.read(inputs / 'INPUTS_FROZEN.json')
    input_protocol = k.read(inputs / 'INPUT_PROTOCOL_FROZEN.json')
    k.require(protocol['input_manifest_sha256'] == k.sha(inputs / 'INPUTS_FROZEN.json') and
              original['protocol_sha256'] == k.sha(inputs / 'INPUT_PROTOCOL_FROZEN.json'), 'Input protocol/manifest identity differs')
    for name, digest in input_protocol['sources'].items():
        k.require(k.sha(inputs / 'source_snapshot' / k.safe(name)) == digest, 'Chemistry producer archive differs')
    registry = k.read(k.HERE / 'REGISTRY.json')
    specs = [registry['chemistry_preflight']] if protocol['preflight'] else registry['cases']
    k.require(not original['refusals'] and input_protocol['registry'] == registry and
              set(original['successes']) == {s['case'] for s in specs} and
              input_protocol['preflight'] == protocol['preflight'] and input_protocol['cases'] == specs, 'Input domain differs')
    geometry = geometry_identity.audit(registry)
    independent, total_constructors = [], 0.
    models = []
    for spec in specs:
        case = spec['case']
        folder = inputs / 'inputs' / case
        path = folder / 'model.json'
        model = k.read(path)
        k.require(k.sha(path) == original['successes'][case] and all(model[key] == value for key, value in spec.items()), 'Actual target geometry/active space differs')
        for key in ('basis', 'unit', 'charge', 'spin'):
            k.require(model[key] == registry[key], 'Molecular convention differs')
        n, electrons = model['n'], model['electrons']
        k.require(n == 6 and electrons in (2, 4) and 0 < len(model['coefficients']) <= 256 and
                  0 < len(model['blocks']) <= 256, 'Independent model bounds')
        target_terms = {tuple((q, p) for q, p in w): int(c) for w, c in model['coefficients']}
        supplied = [(tuple((q, p) for q, p in w), int(c)) for raw in model['blocks'] for w, c in raw]
        k.require(len(target_terms) == len(model['coefficients']) == len(supplied) and
                  dict(supplied) == target_terms and all(target_terms), 'Independent complete target partition differs')
        dimension = 1 << n
        pool = {tuple((q, p) for q, p in w): literal(w, n) for w, c in model['coefficients']}
        H = sum(int(c) / int(model['denominator']) * pool[tuple((q, p) for q, p in w)]
                for w, c in model['coefficients']).toarray()
        identity, destruction = molecular_identity(folder, model, H)
        numbers = [np.array([sum((v >> q) & 1 for q in range(p, n, 2)) for v in range(dimension)]) for p in (0, 1)]
        sector = [v for v in range(dimension) if all(count[v] == electrons // 2 for count in numbers)]
        for raw in model['blocks']:
            words = [tuple((q, p) for q, p in w) for w, c in raw]
            for at, word in enumerate(words):
                for other in words[at + 1:]:
                    commutator = pool[word] @ pool[other] - pool[other] @ pool[word]
                    commutator.eliminate_zeros()
                    k.require(commutator.nnz == 0, 'Complete unit does not factor into commuting Pauli words')
            G = sum(int(c) / int(model['denominator']) * pool[tuple((q, p) for q, p in w)] for w, c in raw).toarray()
            k.require(np.max(abs(G - G.conj().T)) < 1e-12, 'Non-Hermitian complete unit')
            for charge in numbers:
                k.require(np.max(abs(G * (charge[:, None] - charge[None, :]))) < 2e-12, 'Complete-unit sector conservation differs')
        ref = k.read(output / 'references' / (case + '.json'))
        vectors_path = output / k.safe(ref['vectors_file'])
        k.require(ref['model_sha256'] == k.sha(path) and ref['vectors_sha256'] == k.sha(vectors_path), 'State archive identity differs')
        with np.load(vectors_path, allow_pickle=False) as z:
            vectors = z['vectors']
        k.require(vectors.shape == (dimension, 5), 'Full reference-state cohort missing')
        values = np.linalg.eigvalsh(H[np.ix_(sector, sector)])
        g, e = vectors[:, 1], vectors[:, 2]
        k.require(abs(np.vdot(g, e)) < 1e-12 and
                  np.max(abs(vectors[:, 3] - (g + e) / np.sqrt(2))) < 1e-12 and
                  np.max(abs(vectors[:, 4] - (g + 1j * e) / np.sqrt(2))) < 1e-12, 'Coherent references lost their stated phase')
        for index, label in enumerate(k.STATES):
            psi = vectors[:, index]
            record = ref['states'][label]
            meta = record['metadata']
            k.require(np.max(abs(psi[[v for v in range(dimension) if v not in sector]])) < 1e-14 and
                      abs(np.linalg.norm(psi) - 1) < 1e-12 and meta['sector'] == sector, 'Preparation sector/norm differs')
            k.require(meta['support'] == np.flatnonzero(psi != 0).tolist() and
                      abs(meta['norm'] - np.linalg.norm(psi)) < 1e-12 and
                      abs(meta['participation_ratio'] - 1 / sum(abs(psi)**4)) < 1e-10, 'State support/participation metadata differs')
            if label == 'HF':
                expected = np.zeros(dimension, complex)
                expected[(1 << electrons) - 1] = 1
                k.require(np.array_equal(psi, expected) and record['phases']['eigensystem_CPU_seconds'] == 0, 'HF preparation used an uncharged eigensystem')
            else:
                k.require(record['phases']['eigensystem_CPU_seconds'] > 0, 'Non-HF eigensystem preparation is uncharged')
                k.require(np.max(abs(np.array(meta['energies']) - values)) < 2e-10, 'Sector energy labels differ')
                k.require(meta['near_degenerate'] == (meta['ground_excited_gap'] <= k.DEGENERACY_LABEL), 'Degeneracy label differs')
            if label in ('ground', 'excited'):
                at = index - 1
                residual = float(np.linalg.norm(H @ psi - values[at] * psi))
                k.require(residual < 2e-10 and abs(residual - meta['eigen_residual']) < 2e-10, 'Claimed eigenstate is not a sector eigenstate')
                pivot = psi[np.argmax(abs(psi))]
                k.require(abs(pivot.imag) < 1e-12 and pivot.real > 0, 'Eigenstate phase anchor changed')
            annihilated = [A @ psi for A in destruction]
            gamma = np.array([[np.vdot(a, b) for b in annihilated] for a in annihilated])
            natural = np.linalg.eigvalsh(gamma)
            defect = float(np.trace(gamma - gamma @ gamma).real)
            correlation = 'correlated' if defect > 1e-10 else 'Slater-like' if defect <= 1e-12 else 'indeterminate'
            k.require(np.max(abs(natural - meta['natural_occupations'])) < 2e-12 and
                      abs(defect - meta['one_RDM_idempotency_defect']) < 2e-12 and
                      correlation == meta['correlation_label'], 'RDM correlation metadata differs')
            k.require(all(np.isfinite(x) and x >= 0 for x in record['phases'].values()) and
                      abs(sum(record['phases'].values()) + record['accounting_overhead_CPU_seconds'] - record['CPU_seconds']) < 1e-12,
                      'Preparation timing does not close')
            total_constructors += record['CPU_seconds']
            independent.append(dict(case=case, state=label, correlation_label=correlation,
                       one_RDM_idempotency_defect=defect, natural_occupations=natural.tolist(),
                       preparation_CPU_seconds=record['CPU_seconds'], full_dimension=dimension, sector_dimension=len(sector)))
        models.append(dict(case=case, identity=identity))
        print('INDEPENDENT_REFERENCES_CHECKED', case, flush=True)
    k.require(len(independent) == protocol['expected_states'] == frozen['states'] and
              len(models) == protocol['expected_models'] == frozen['models'], 'Incomplete reference panel')
    costs = k.read(output / 'PREPARATION_COSTS.json')
    k.require(abs(costs['actual_constructor_CPU_seconds'] - total_constructors) < 1e-10, 'Actual preparation campaign charge differs')
    k.write(receipt / 'INDEPENDENT_REFERENCE_LABELS.json', independent)
    report = dict(status='passed', models=models, states=len(independent), geometry_audit=geometry,
                  producer_sources=protocol['sources'], verifier_sources=k.sources(),
                  reference_freeze_sha256=k.sha(output / 'REFERENCE_STATES_FROZEN.json'),
                  independent_labels_sha256=k.sha(receipt / 'INDEPENDENT_REFERENCE_LABELS.json'),
                  CPU_seconds=time.process_time() - began, wall_seconds=time.monotonic() - wall,
                  scope='LiteralPauli,CAR/tensorandstate/sector/RDM checks;notTrotter-errorvalidationorLeanproof')
    k.write(receipt / 'REFERENCE_VERIFIED.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt-dir', type=Path)
    args = parser.parse_args()
    verify(args.output.resolve(), args.receipt_dir.resolve() if args.receipt_dir else None)
