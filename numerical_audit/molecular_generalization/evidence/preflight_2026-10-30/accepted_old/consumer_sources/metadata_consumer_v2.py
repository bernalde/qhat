"""DRAFT later-consumer gate; invoke only inside the declared metadata service.

Reads sealed producer/extracted DATA and fresh external numerical reports.
Uses the existing inspect_case report-argument API; never substitutes artifacts
inside DATA or changes an original producer report or identity.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tomllib

if not __debug__:
    raise RuntimeError('Metadata gate requires normal Python assertion checking')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path, *, toml=False):
    path = Path(path)
    size = path.stat().st_size
    assert 0 <= size <= 8*1024**2, 'Bounded metadata artifact'
    with path.open('rb') as stream:
        raw = stream.read(size+1)
    assert len(raw) == size, 'Metadata changed while being read'
    return tomllib.loads(raw.decode()) if toml else json.loads(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-data', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--consumer-root', type=Path, required=True)
    parser.add_argument('--numerical-root', type=Path, required=True)
    parser.add_argument('--producer-analysis', type=Path, required=True)
    parser.add_argument('--metadata-plan', type=Path, required=True)
    args = parser.parse_args()
    data, original, output = (p.resolve(strict=True) for p in
                             (args.data, args.producer_data, args.consumer_root))
    assert not output.is_relative_to(data) and not output.is_relative_to(original)
    study = data/'source_snapshot/numerical_audit/molecular_generalization'
    sys.path.insert(0, str(study))
    import common as k
    import panels
    import analyze
    import package

    assert Path(k.__file__).resolve() == study/'common.py'
    assert Path(analyze.__file__).resolve() == study/'analyze.py'
    declared = read(data/'EVIDENCE_INVENTORY.json')
    assert package.inventory(data) == declared
    assert read(original/'EVIDENCE_INVENTORY.json') == declared
    assert package.inventory(original) == declared, 'Sealed producer changed'
    validation = panels.validate(data, preflight=True)
    assert validation['complete'] is True
    k.write_new(output/'PANELS_REPLAY.json', panels.public_validation(validation))

    # This replay consumes ORIGINAL reports and their ORIGINAL owning receipts.
    # It does not assign original provenance to the fresh numerical reports.
    replay = analyze.analyze(data, preflight=True)
    assert replay['complete'] is True and replay['status'] == 'complete_preflight'
    producer_analysis = read(args.producer_analysis)
    # The figure PLAN was created after ANALYSIS.json. Authenticate that
    # single additive ancillary digest; every original digest and all
    # scientific fields must still match exactly.
    original_hashes=producer_analysis['artifact_hashes']
    replay_hashes=replay['artifact_hashes']
    assert all(replay_hashes.get(key)==value for key,value in original_hashes.items())
    added=set(replay_hashes)-set(original_hashes)
    assert added=={'jobs/mg-old-preflight-figures-a4.PLAN.json'}
    for key in added:
        ancillary=read(data/key)
        assert ancillary['phase']=='analysis' and ancillary['source_freeze_sha256']==validation['source_freeze_sha256']
        assert ancillary['command'][1].endswith('/numerical_audit/molecular_generalization/plot.py')
        assert replay_hashes[key]==declared['files'][key]['sha256']
    comparison=json.loads(json.dumps(replay,allow_nan=False))
    comparison['artifact_hashes']=original_hashes
    assert comparison == producer_analysis, 'Original scientific report semantics changed'
    k.write_new(output/'ORIGINAL_ANALYSIS_REPLAY.json', replay)

    a, freeze = validation['_artifacts'], validation['_freeze']
    info, _ = k.header(data)
    assert info['cases'] == ['LiH_r1p487_631g_14m4e'] and info['repetitions'] == 1
    case = info['cases'][0]
    inventory = validation['_inventories'][case]
    panel = a.toml(f'panels/{case}.toml')
    witnesses = a.toml(f'witnesses/{case}.toml')
    numerical=args.numerical_root.resolve(strict=True)
    fresh_paths = {'independent': numerical/'independent'/case/'VERIFY.toml',
                   'objective': numerical/'objective'/case/'OBJECTIVE.toml'}
    reports = {'physical': a.toml(f'physical/{case}/PHYSICAL.toml')}
    reports.update({name: read(path, toml=True) for name, path in fresh_paths.items()})
    # Existing metadata join validates exact source/model/panel/witness hashes,
    # all selected rows and persisted raw claims, native/floor status, original
    # binary/circuit joins, reference residuals and every fixed tolerance.
    joined = analyze.inspect_case(a, freeze, info, inventory, panel, witnesses, reports)
    assert joined['complete'] and joined['raw_claims_verified'] and not joined['issues']
    assert len(joined['good']) == 36 and len(panel['candidates']) == 36
    assert len(reports['physical']['rows']) == reports['physical']['completed_logical_rows'] == 36
    assert len(reports['independent']['rows']) == reports['independent']['completed_logical_rows'] == 36
    assert reports['independent']['reference_check_count'] == len(reports['independent']['reference_checks']) == 2
    assert reports['objective']['processed_logical_rows'] == len(reports['objective']['logical_rows']) == 36
    assert reports['objective']['complete'] is True and reports['objective']['witness_inventory_complete'] is True
    assert reports['objective']['all_scored_objectives_verified'] is True
    assert len(reports['objective']['raw_witnesses']) == len(witnesses['witnesses']) == reports['objective']['processed_raw_witnesses']
    assert all(reports[name]['status'] == 'complete' and not reports[name]['failures'] for name in reports)
    assert not reports['independent']['censored_rows']
    assert reports['independent']['julia_version'] == '1.12.6'
    assert reports['independent']['julia_threads'] == reports['independent']['blas_threads'] == 1
    assert package.inventory(data) == declared and package.inventory(original) == declared
    # Own receipt/FINAL will be bound by the external dispatcher after teardown.
    k.write_new(output/'REPLAY_CHECKS.json', dict(schema='molecular-relocated-consumer-checks-v1',
        complete=True, role='later_relocated_numerical_consumer', preflight_only=True,
        source_freeze_sha256=k.sha(data/'INPUT_PROTOCOL_FROZEN.json'),
        original_producer_commit=freeze['source_commit'], original_producer_root=freeze['producer_root'],
        original_producer_data_root=freeze['producer_data_root'],
        original_frozen_sources=freeze['sources'], current_consumer_sources=k.source_map(),
        current_consumer_environment=k.prior.environment(),
        executing_metadata_consumer=dict(path=str(Path(__file__).resolve()), sha256=sha(__file__)),
        metadata_plan=dict(path=str(args.metadata_plan.resolve()), sha256=sha(args.metadata_plan)),
        producer_analysis=dict(path=str(args.producer_analysis.resolve()), sha256=sha(args.producer_analysis)),
        original_analysis_semantically_equal=True, sealed_trees_unchanged=True,
        additive_post_analysis_plot_plan_hashes={key:replay_hashes[key] for key in sorted(added)},
        comparison_scope="All original digests and scientific fields exact; one later plot PLAN digest separately authenticated",
        fresh_reports={name:dict(path=str(path),sha256=sha(path)) for name,path in fresh_paths.items()},
        model_sha256=sha(data/f'inputs/{case}/model.toml'), panel_sha256=sha(data/f'panels/{case}.toml'),
        witnesses_sha256=sha(data/f'witnesses/{case}.toml'),
        physical_sha256=sha(data/f'physical/{case}/PHYSICAL.toml'),
        logical_identities=[list(value) for value in sorted(joined['good'])],
        logical_rows=36, reference_checks=2, raw_claims=len(witnesses['witnesses']),
        raw_claims_verified=True, receipt_gate='pending_external_post_teardown_FINAL'))


if __name__ == '__main__':
    main()
