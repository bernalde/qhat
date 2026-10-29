"""Freeze and generate the two declared fresh inputs without observing errors."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import shutil
import time
import common as k


def prepare(output, preflight):
    k.require(not output.exists(), 'Fresh input output required')
    config = k.validate_protocol(k.configuration())
    registry = k.registry()
    k.require(len(registry['cases']) == 2 and len({v['case'] for v in registry['cases']}) == 2,
              'Exactly two declared new cases')
    k.require(registry['base_commit'] == k.BASE, 'Current merged geometry-audit baseline')
    head = k.committed_sources()
    preflight = preflight.resolve(strict=True)
    proof_frozen = k.read(preflight/'INPUT_PROTOCOL_FROZEN.json')
    k.require(proof_frozen.get('mode') == 'preflight', 'A representative development preflight is required')
    sources = k.source_map()
    original_sources = proof_frozen['sources']
    consumer = 'numerical_audit/cost_matched_continuation/analyze.py'
    k.require(original_sources.keys() == sources.keys(), 'Preflight producer source inventory differs')
    delta = {name:dict(original=original_sources[name], replayed=digest)
             for name, digest in sources.items() if original_sources[name] != digest}
    k.require(set(delta) <= {consumer}, 'Changed numerical producers require a fresh representative preflight')
    k.require(proof_frozen['configuration'] == config and proof_frozen['registry'] == registry and
              proof_frozen['environment'] == k.prior.environment(), 'Preflight protocol/environment differs')
    checker = k._load('_continuation_preflight_gate', k.HERE/'analyze.py')
    replayed = checker.analyze(preflight)
    k.require(replayed['study_complete'] and replayed['mode'] == 'preflight' and
              replayed['primary']['decision'] == 'preflight_only', 'Preflight verification/admission failed')
    _, identity = k.chemistry_helpers()
    novelty = identity.audit(registry)
    output.mkdir(parents=True)
    for name in sources:
        target = output/'source_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT/name, target)
    shutil.copytree(preflight, output/'preflight')
    consumer_sha = k.sha(k.HERE/'analyze.py')
    source_path = 'preflight/validation/analyze.'+consumer_sha+'.py'
    analysis_path = 'preflight/validation/ANALYSIS.'+consumer_sha+'.json'
    (output/source_path).parent.mkdir(parents=True, exist_ok=True)
    if (output/source_path).exists():
        k.require(k.sha(output/source_path) == consumer_sha, 'Existing preflight consumer differs')
    else:
        with (output/source_path).open('xb') as stream:
            stream.write((k.HERE/'analyze.py').read_bytes())
    if (output/analysis_path).exists():
        k.require(k.read(output/analysis_path) == replayed, 'Existing preflight replay differs')
    else:
        k.write_new(output/analysis_path, replayed)
    proof = dict(path='preflight', analysis_path=analysis_path, analysis_sha256=k.sha(output/analysis_path),
        source_freeze_sha256=k.sha(preflight/'INPUT_PROTOCOL_FROZEN.json'),
        consumer_source_path=source_path, consumer_sha256=consumer_sha, consumer_only_source_delta=delta)
    k.write_new(output/'INPUT_PROTOCOL_FROZEN.json', dict(mode='main', source_commit=head, sources=sources,
        configuration_sha256=k.sha(k.HERE/'PROTOCOL.json'),
        registry=registry, configuration=config, environment=k.prior.environment(),
        geometry_identity=novelty, preflight=proof,
        scope='New geometry inputs; no ordering or physical quality observed'))


def generate_case(output, name):
    frozen = k.frozen_protocol(output)
    specs = [v for v in frozen['registry']['cases'] if v['case'] == name]
    k.require(len(specs) == 1, 'Exactly one declared case')
    target = output/'inputs'/name
    k.require(not target.exists(), 'Retain old input attempts')
    producer, _ = k.chemistry_helpers()
    cpu, wall = time.process_time(), time.monotonic()
    try:
        model = producer.generate_one(specs[0], frozen['registry'], target)
        model.update(generation_CPU_seconds=time.process_time()-cpu,
                     generation_wall_seconds=time.monotonic()-wall)
        k.prior.blocks(model)
        k.require(len((json.dumps(model, indent=2, allow_nan=False)+'\n').encode()) <= k.prior.MAX_MODEL_BYTES,
                  'Generated model byte bound')
        k.write_new(target/'model.json', model)
        k.prior.write_model_toml(target/'model.toml', model, source_json_sha256=k.sha(target/'model.json'))
    except Exception as exc:
        target.mkdir(parents=True, exist_ok=True)
        k.write_new(target/'REFUSAL.json', dict(type=type(exc).__name__, message=str(exc),
            CPU_seconds=time.process_time()-cpu, wall_seconds=time.monotonic()-wall))
        print('INPUT_REFUSAL', name, type(exc).__name__, flush=True)
        return False
    print('GENERATED', name, 'words', len(model['coefficients']), 'units', len(model['blocks']), flush=True)
    return True


def finalize(output):
    frozen = k.frozen_protocol(output)
    successes, refusals = {}, {}
    for spec in frozen['registry']['cases']:
        folder = output/'inputs'/spec['case']
        files = ('model.json', 'model.toml', 'orbitals.npz', 'tensors.npz')
        if (folder/'REFUSAL.json').exists():
            refusals[spec['case']] = dict(sha256=k.sha(folder/'REFUSAL.json'))
        elif all((folder/name).is_file() for name in files):
            model = k.read(folder/'model.json')
            k.require((folder/'model.json').stat().st_size <= k.prior.MAX_MODEL_BYTES,
                      'Serialized model exceeds byte cap')
            k.require(model['tensors_sha256'] == k.sha(folder/'tensors.npz') and
                      model['orbitals_sha256'] == k.sha(folder/'orbitals.npz'), 'Model tensor/orbital hashes differ')
            k.require(all(model[key] == value for key, value in spec.items()), 'Declared geometry identity')
            k.prior.blocks(model)
            successes[spec['case']] = {name:k.sha(folder/name) for name in files}
        else:
            refusals[spec['case']] = dict(reason='Missing output; inspect retained job receipt')
    k.write_new(output/'INPUTS_FROZEN.json', dict(status='complete' if not refusals else 'incomplete',
        expected_cases=2, successes=successes, refusals=refusals,
        protocol_sha256=k.sha(k.HERE/'PROTOCOL.json'),
        source_freeze_sha256=k.sha(output/'INPUT_PROTOCOL_FROZEN.json')))
    return not refusals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=k.Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true')
    group.add_argument('--case')
    group.add_argument('--finalize', action='store_true')
    parser.add_argument('--preflight', type=k.Path, help='Verified development DATA_ROOT; required with --prepare')
    args = parser.parse_args()
    if args.prepare and args.preflight is None:
        parser.error('--prepare requires --preflight')
    result = prepare(args.output, args.preflight) if args.prepare else generate_case(args.output, args.case) if args.case else finalize(args.output)
    raise SystemExit(2 if result is False else 0)
