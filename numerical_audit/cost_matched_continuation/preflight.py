"""Prepare the declared old development input, without generating holdouts."""
import argparse
import shutil
import common as k


def prepare(output):
    k.require(not output.exists(), 'Fresh preflight required')
    config = k.validate_protocol(k.configuration())
    case = config['preflight']['case']
    retained = k.OLD/'evidence/run_2026-10-28/data-a6'
    old_index = k.read(retained/'INPUTS_FROZEN.json')
    files = old_index['successes'][case]
    head = k.committed_sources()
    sources = k.source_map()
    output.mkdir(parents=True)
    for name in sources:
        target = output/'source_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT/name, target)
    k.write_new(output/'INPUT_PROTOCOL_FROZEN.json', dict(mode='preflight', source_commit=head,
        sources=sources, configuration_sha256=k.sha(k.HERE/'PROTOCOL.json'), registry=k.registry(),
        configuration=config, environment=k.prior.environment(),
        reuse=dict(path=str(retained.relative_to(k.ROOT)), input_index_sha256=k.sha(retained/'INPUTS_FROZEN.json')),
        scope='Reused merged development model; neither fresh geometry nor new SCF'))
    target = output/'inputs'/case
    target.mkdir(parents=True)
    for name in ('model.json', 'model.toml', 'orbitals.npz', 'tensors.npz'):
        k.require(k.sha(retained/'inputs'/case/name) == files[name], 'Retained input digest differs')
        shutil.copyfile(retained/'inputs'/case/name, target/name)
    model = k.read(target/'model.json')
    k.require((target/'model.json').stat().st_size <= k.prior.MAX_MODEL_BYTES and
              model['tensors_sha256'] == files['tensors.npz'] and
              model['orbitals_sha256'] == files['orbitals.npz'], 'Retained model/tensor/orbital binding')
    k.write_new(output/'INPUTS_FROZEN.json', dict(status='complete', expected_cases=1,
        successes={case: files}, refusals={}, protocol_sha256=k.sha(k.HERE/'PROTOCOL.json'),
        source_freeze_sha256=k.sha(output/'INPUT_PROTOCOL_FROZEN.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=k.Path, required=True)
    prepare(parser.parse_args().output)
