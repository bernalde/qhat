"""DRAFT later MAIN consumer; execute only inside its declared metadata guard.

Adapted from the adopted old metadata_consumer_v2.py. Original reports and fresh
external numerical reports have separate identities. No DATA/report overrides.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import censor_contract as cc

if not __debug__:
    raise RuntimeError('Metadata gate requires normal Python assertion checking')

CASES = ['HF_r1p00_631g_12m8e', 'H2S_r1p40_a96_631g_12m8e']
FIGURE_PLAN = 'jobs/mg-main-figures-a1.PLAN.json'


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


def ref(path):
    return dict(path=str(Path(path).resolve()), sha256=sha(path))


def observed_counters(record):
    """Preserve reported counters; absence is not a zero or inferred work count."""
    result = record.get('result') or {}
    stats = result.get('stats')
    stats = stats if isinstance(stats, dict) else {}
    scopes = [('result', result), ('result.stats', stats)]
    for prefix, parent in (('result', result), ('result.stats', stats)):
        raw = parent.get('raw_julia_candidate')
        if isinstance(raw, dict):
            scopes.append((prefix+'.raw_julia_candidate', raw))
            witness = raw.get('witness')
            if isinstance(witness, dict):
                scopes.append((prefix+'.raw_julia_candidate.witness', witness))
    values = {}
    for scope, source in scopes:
        if not isinstance(source, dict):
            continue
        for key, value in source.items():
            if ('eval' in key or key in ('full_scores', 'start_scores_completed', 'proposals',
                    'chain_proposals', 'accepted_moves', 'chain_accepted')) and (
                    type(value) in (int, float) or isinstance(value, list)):
                values[scope+'.'+key] = value
    return dict(constructor=dict(path=record['artifact_path'], sha256=record['artifact_sha256']),
                reported_counters=values, scope='Only actual persisted counter fields; no absent-to-zero conversion')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('producer-data', 'data', 'consumer-root', 'numerical-root',
                 'producer-analysis', 'producer-visual-gate', 'metadata-plan'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    data, original, output, numerical = (p.resolve(strict=True) for p in
        (args.data, args.producer_data, args.consumer_root, args.numerical_root))
    assert not output.is_relative_to(data) and not output.is_relative_to(original)
    assert not numerical.is_relative_to(data) and not numerical.is_relative_to(original)
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
    assert package.inventory(original) == declared, 'Sealed original producer changed'
    validation = panels.validate(data, preflight=False)
    assert validation['complete'] is True
    k.write_new(output/'PANELS_REPLAY.json', panels.public_validation(validation))

    # Replay ORIGINAL report semantics against ORIGINAL owning receipts.
    replay = analyze.analyze(data, preflight=False)
    cc.audit_analysis(replay)
    producer_analysis = read(args.producer_analysis)
    cc.audit_analysis(producer_analysis)
    original_hashes, replay_hashes = producer_analysis['artifact_hashes'], replay['artifact_hashes']
    assert all(replay_hashes.get(key) == value for key, value in original_hashes.items())
    gate = read(args.producer_visual_gate)
    post = gate['post_analysis_figure_plan']
    admitted = {}
    if post is not None:
        assert set(post) == {'path', 'sha256'} and post['path'] == FIGURE_PLAN
        assert sha(data/FIGURE_PLAN) == post['sha256'] == declared['files'][FIGURE_PLAN]['sha256']
        ancillary = read(data/FIGURE_PLAN)
        assert ancillary['label'] == 'mg-main-figures-a1' and ancillary['phase'] == 'analysis'
        assert ancillary['source_freeze_sha256'] == validation['source_freeze_sha256']
        frozen = validation['_freeze']
        command = ancillary['command']
        assert len(command) == 4 and command[:3] == [frozen['runtime']['executable'],
            frozen['producer_root']+'/numerical_audit/molecular_generalization/plot.py',
            str(args.producer_analysis.resolve())]
        assert replay_hashes[FIGURE_PLAN] == post['sha256']
        admitted[FIGURE_PLAN] = post['sha256']
    added = set(replay_hashes)-set(original_hashes)
    assert added == set(admitted)-set(original_hashes), 'No other additive digest exception'
    comparison = json.loads(json.dumps(replay, allow_nan=False))
    comparison['artifact_hashes'] = original_hashes
    assert comparison == producer_analysis, 'Every original scientific field/digest must remain identical'
    k.write_new(output/'ORIGINAL_ANALYSIS_REPLAY.json', replay)

    a, freeze = validation['_artifacts'], validation['_freeze']
    info, _ = k.header(data)
    assert info['mode'] == 'main' and info['cases'] == CASES and info['repetitions'] == 2
    case_checks, identities, constructors = [], [], []
    for case in CASES:
        inventory = validation['_inventories'][case]
        assert len(inventory['records']) == 18
        panel = a.toml(f'panels/{case}.toml')
        witnesses = a.toml(f'witnesses/{case}.toml')
        fresh_paths = {'independent': numerical/'independent'/case/'VERIFY.toml',
                       'objective': numerical/'objective'/case/'OBJECTIVE.toml'}
        reports = {'physical': a.toml(f'physical/{case}/PHYSICAL.toml')}
        reports.update({name: read(path, toml=True) for name, path in fresh_paths.items()})
        # Existing report-argument API; all binary/model/source/order/raw-claim
        # and tolerance joins remain unchanged. Fresh files stay external.
        joined = analyze.inspect_case(a, freeze, info, inventory, panel, witnesses, reports)
        cc.audit_inspection(case,joined)
        objective_summary=cc.audit_objective(case,reports['objective'])
        assert len(joined['good']) == len(panel['candidates']) == 72
        assert len(reports['physical']['rows']) == reports['physical']['completed_logical_rows'] == 72
        assert len(reports['independent']['rows']) == reports['independent']['completed_logical_rows'] == 72
        assert reports['independent']['reference_check_count'] == len(reports['independent']['reference_checks']) == 2
        assert reports['objective']['processed_logical_rows'] == len(reports['objective']['logical_rows']) == 72
        assert reports['objective']['witness_inventory_complete'] is True
        assert len(reports['objective']['raw_witnesses']) == len(witnesses['witnesses']) == reports['objective']['processed_raw_witnesses']
        assert reports['objective']['unique_objectives'] == len(reports['objective']['checks'])
        assert all(reports[name]['status']=='complete' and not reports[name]['failures'] for name in ('physical','independent'))
        assert not reports['objective']['failures']
        assert not reports['independent']['censored_rows']
        assert reports['independent']['julia_version'] == '1.12.6'
        assert reports['independent']['julia_threads'] == reports['independent']['blas_threads'] == 1
        identities.extend([[case, *value] for value in sorted(joined['good'])])
        constructors.extend(observed_counters(record) for record in inventory['records'].values())
        case_checks.append(dict(case=case,report_processing_complete=True,
            scientific_gate_complete=joined['complete'],objective_status=reports['objective']['status'],
            selected_objectives_verified=objective_summary['selected_verified'],logical_rows=len(joined['good']),
            reference_checks=len(reports['independent']['reference_checks']),
            raw_claims=len(witnesses['witnesses']),raw_claims_verified=joined['raw_claims_verified'],
            raw_proposed_comparisons_verified=True,native_refused_raw_names=objective_summary['native_refused_raw_names'],
            unique_objective_evaluations=len(reports['objective']['checks']),
            reported_objective_work=reports['objective'].get('charged_work_compile'),
            unique_unfused_circuit_replays=reports['independent']['independent_unfused_replay_count'],
            fresh_reports={name:ref(path) for name,path in fresh_paths.items()},
            model_json=ref(data/f'inputs/{case}/model.json'), model_toml=ref(data/f'inputs/{case}/model.toml'),
            panel=ref(data/f'panels/{case}.toml'), witnesses=ref(data/f'witnesses/{case}.toml'),
            original_physical=ref(data/f'physical/{case}/PHYSICAL.toml')))
    assert len(identities) == len({tuple(value) for value in identities}) == 144
    assert len(constructors) == 36 and sum(value['reference_checks'] for value in case_checks) == 4
    assert package.inventory(data) == declared and package.inventory(original) == declared
    k.write_new(output/'REPLAY_CHECKS.json', dict(schema='molecular-censored-relocated-consumer-checks-v1',
        delivery_complete=True,replay_complete=True,scientific_gate_complete=False,source_analysis_complete=False,
        selected_objectives_verified=144,raw_proposed_comparisons_verified=True, role='later_relocated_numerical_consumer', preflight_only=False, mode='main',
        source_freeze_sha256=k.sha(data/'INPUT_PROTOCOL_FROZEN.json'),
        original_producer_commit=freeze['source_commit'], original_producer_root=freeze['producer_root'],
        original_producer_data_root=freeze['producer_data_root'], original_frozen_sources=freeze['sources'],
        current_consumer_sources=k.source_map(), current_consumer_environment=k.prior.environment(),
        executing_metadata_consumer=ref(__file__), metadata_plan=ref(args.metadata_plan),
        producer_visual_gate=ref(args.producer_visual_gate), producer_analysis=ref(args.producer_analysis),
        original_analysis_semantically_equal=True, sealed_trees_unchanged=True,
        additive_post_analysis_plot_plan_hashes={key:replay_hashes[key] for key in sorted(added)},
        comparison_scope='All original digests and scientific fields exact; only exact root-bound main figure PLAN may be additive',
        cases=case_checks, logical_identities=identities, logical_rows=len(identities),
        reference_checks=sum(value['reference_checks'] for value in case_checks),
        raw_claims=sum(value['raw_claims'] for value in case_checks),raw_claims_verified=False,
        known_native_refusal=dict(case=cc.HF,name='rep2.window.B64.raw1',selected=False,native_score_absent=True),
        censor_contract=ref(cc.__file__),expected_outcomes=ref(Path(cc.__file__).with_name('EXPECTED_OUTCOMES.json')),
        source_primary_decisions=dict(window_minus_full='incomplete',annealing_minus_full='incomplete'),
        constructor_reported_counters=constructors,
        receipt_gate='pending_external_post_teardown_FINAL'))


if __name__ == '__main__':
    main()
