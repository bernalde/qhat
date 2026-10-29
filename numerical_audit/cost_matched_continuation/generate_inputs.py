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


def prepare(output):
    k.require(not output.exists(), 'Fresh input output required')
    config = k.validate_protocol(k.configuration())
    registry = k.registry()
    k.require(len(registry['cases']) == 2 and len({v['case'] for v in registry['cases']}) == 2,
              'Exactly two declared new cases')
    k.require(registry['base_commit'] == k.BASE, 'Current merged geometry-audit baseline')
    head = k.committed_sources()
    _, identity = k.chemistry_helpers()
    novelty = identity.audit(registry)
    sources = k.source_map()
    output.mkdir(parents=True)
    for name in sources:
        target = output/'source_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT/name, target)
    k.write_new(output/'INPUT_PROTOCOL_FROZEN.json', dict(source_commit=head, sources=sources,
        configuration_sha256=k.sha(k.HERE/'PROTOCOL.json'),
        registry=registry, configuration=config, environment=k.prior.environment(),
        geometry_identity=novelty, scope='New geometry inputs; no ordering or physical quality observed'))


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
        k.require(len(json.dumps(model, allow_nan=False).encode()) <= k.prior.MAX_MODEL_BYTES,
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
    args = parser.parse_args()
    result = prepare(args.output) if args.prepare else generate_case(args.output, args.case) if args.case else finalize(args.output)
    raise SystemExit(2 if result is False else 0)
