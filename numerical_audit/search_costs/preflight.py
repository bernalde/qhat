"""Representative largest-shape end-to-end preflight on altered coefficients."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tarfile
import time
import numpy as np
from scipy.sparse.linalg import expm_multiply
from bounded import run_arm, engine as e
from models import synthetic
from study import common, initializer, decorate, load_previous_placement, write, sha, CircuitTarget, physical_metrics, ARMS
from verify_costs import Reference, metrics, audit_ledger


def run(output):
    if output.exists():
        raise ValueError('Fresh preflight directory required; retain failed attempts')
    output.mkdir(parents=True)
    began_cpu, began_wall = time.process_time(), time.monotonic()
    model = synthetic(12)
    model['case'] = 'preflight_altered_extended_hubbard_12'
    for block in model['blocks']:
        for term in block:
            term[1] = str(7*int(term[1]))
    for term in model['coefficients']:
        term[1] = str(7*int(term[1]))
    model['denominator'] = str(5*int(model['denominator']))
    model['removed_scalar_exact'] = str(Fraction(model['removed_scalar_exact'])*Fraction(7,5))
    mp = output/'input.json'
    write(mp, model)
    context = common(mp)
    warm = initializer(context, load_previous_placement())
    records, warmups = [], []
    for partition, objective, initial in ARMS:
        bs = context['paulis'] if partition == 'pauli' else context['blocks']
        cap = 20*len(bs)*(len(bs)-1)
        warmups.append(run_arm(bs, context['D'], context['v'], context['den'], objective, initial,
                              warm[partition], .5, context['cpu_seconds'], warm['charges'][partition], cap))
        for budget in (.5, 2.):
            record = run_arm(bs, context['D'], context['v'], context['den'], objective, initial,
                             warm[partition], budget, context['cpu_seconds'], warm['charges'][partition], cap)
            record.update(arm=partition+'_'+objective+'_'+initial, repetition=0)
            records.append(decorate(record, context, partition))
    bs = context['paulis']
    larger = run_arm(bs, context['D'], context['v'], context['den'], 'JD', 'cold', warm['pauli'],
                     8., context['cpu_seconds'], warm['charges']['pauli'], 20*len(bs)*(len(bs)-1))
    larger.update(arm='pauli_JD_larger_CPU', repetition=0)
    records.append(decorate(larger, context, 'pauli'))
    baseline = dict(arm='current_placement', order=warm['conserving'], score=warm['result']['score'], objective='JD')
    records.append(decorate(baseline, context, 'conserving'))
    write(output/'construction.json', dict(records=records, warmups=warmups, common=context['cpu_seconds'], initializer=warm))
    construction_cpu = time.process_time()-began_cpu
    frozen_design = sha(output/'construction.json')
    target = CircuitTarget(model)
    initial = np.zeros(target.dim, complex)
    initial[target.v] = 1
    exact = expm_multiply(-1j*target.hamiltonian(), initial)
    rows = []
    for index, record in enumerate(records):
        for formula in (1, 2):
            steps = e.feasible_steps(record['unit_lengths'], 10*len(model['coefficients']), formula == 2)
            state, rotations = target.propagate(record, steps, formula)
            label = str(index)+'.f'+str(formula)
            path = output/'states'/(label+'.npz')
            path.parent.mkdir(exist_ok=True)
            np.savez_compressed(path, actual=state)
            rows.append(dict(index=index, formula=formula, steps=steps, rotations=rotations,
                             state_file=str(path.relative_to(output)), **physical_metrics(exact, state, model)))
    write(output/'rows.json', rows)
    propagation_cpu = time.process_time()-began_cpu-construction_cpu
    before_verify = time.process_time()
    reference = Reference(model)
    independent_exact = reference.exact_action()
    reference_delta = float(np.max(abs(exact-independent_exact)))
    assert reference_delta < 2e-10
    max_state_delta = max_score_delta = max_metric_delta = 0.
    for record in records:
        if record['arm'] != 'current_placement':
            audit_ledger(record, len(record['order']))
        score = reference.score(record)
        if record['score'] is not None:
            max_score_delta = max(max_score_delta, abs(score-record['score']))
            assert np.isclose(score, record['score'], rtol=2e-9, atol=5e-18)
    for row in rows:
        record = records[row['index']]
        independent = reference.replay(record, row['steps'], row['formula'])
        with np.load(output/row['state_file']) as data:
            state = data['actual']
        delta = float(np.max(abs(state-independent)))
        max_state_delta = max(max_state_delta, delta)
        assert delta < 2e-10
        for key, value in metrics(independent_exact, independent, model).items():
            delta = abs(value-row[key])
            max_metric_delta = max(max_metric_delta, delta)
            assert delta < 2e-12
    verification_cpu = time.process_time()-before_verify
    before_package = time.process_time()
    reloaded = json.loads((output/'rows.json').read_text())
    assert len(reloaded) == 36 and frozen_design == sha(output/'construction.json')
    archive = output/'preflight-evidence.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for path in sorted(output.rglob('*')):
            if path.is_file() and path != archive:
                tar.add(path, arcname=str(path.relative_to(output)), recursive=False)
    extraction = output/'extracted'
    with tarfile.open(archive) as tar:
        tar.extractall(extraction, filter='data')
    for path in extraction.rglob('*'):
        if path.is_file():
            assert sha(path) == sha(output/path.relative_to(extraction))
    package_cpu = time.process_time()-before_package
    stats = dict(status='passed', scope='Altered largest synthetic target, not prospective panel errors',
                 modes=12, Pauli_words=len(model['coefficients']), records=len(records), circuits=len(rows),
                 construction_CPU_seconds=construction_cpu, propagation_CPU_seconds=propagation_cpu,
                 verification_CPU_seconds=verification_cpu, extraction_packaging_replay_CPU_seconds=package_cpu,
                 CPU_seconds=time.process_time()-began_cpu, wall_seconds=time.monotonic()-began_wall,
                 max_state_delta=max_state_delta, max_score_delta=max_score_delta, max_metric_delta=max_metric_delta,
                 independent_reference_delta=reference_delta, cache_bytes=target.payload_bytes,
                 peak_pair_payload_bytes=max(r.get('array_payload_bytes',0) for r in records),
                 conservative_pair_bound_bytes=max(r.get('storage_upper_bound_bytes',0) for r in records),
                 retained_bytes=sum(p.stat().st_size for p in output.rglob('*') if p.is_file()),
                 archive_sha256=sha(archive))
    write(output/'PREFLIGHT.json', stats)
    print(json.dumps(stats), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
