"""Retain every measured repetition, paired factorial comparison and CPU failure."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import statistics


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def build(results, output):
    if output.exists():
        raise ValueError('Fresh report directory required')
    verified = read(results/'VERIFIED.json')
    assert verified['status'] == 'passed' and verified['decision_agreement'] and verified['rows'] == 864
    assert verified['schedule_freeze_sha256'] == sha(results/'SCHEDULES_FROZEN.json')
    for name, digest in verified['artifacts'].items():
        assert sha(results/name) == digest, name
    freeze = read(results/'SCHEDULES_FROZEN.json')
    for name, digest in freeze['artifacts'].items():
        assert sha(results/name) == digest, name
    rows = [read(p) for p in sorted((results/'rows').glob('*.json'))]
    assert len(rows) == 864
    designs = {x['case']: x for x in [read(p) for p in sorted((results/'designs').glob('*.json'))]}
    assert len(designs) == 8
    physical = {(r['case'], r['record_index'], r['formula']): r for r in rows}
    lookup, searches, setup_rows, warmup_rows, campaign_rows = {}, [], [], [], []
    for case, design in designs.items():
        model = read(results/'inputs'/(case+'.json'))
        for warmup in design['warmups']:
            record=warmup['record']
            warmup_rows.append(dict(case=case,partition=warmup['partition'],objective=warmup['objective'],
                initialization=warmup['initialization'],CPU_budget=record['CPU_budget'],status=record['status'],
                actual_local_CPU_seconds=record['actual_local_CPU_seconds'],setup_CPU_seconds=record['setup_CPU_seconds'],
                search_CPU_seconds=record['search_CPU_seconds'],evaluations=record['evaluations'],
                within_CPU_budget=record['within_CPU_budget']))
        assert len(design['warmups'])==8
        shared_warmup=design['warmups'][0]['record']['common_setup_CPU_charge']
        initializer_warmup=next(x['record']['inherited_placement_CPU_charge'] for x in design['warmups']
            if x['partition']=='pauli' and x['initialization']=='placement')
        common_CPU=shared_warmup+sum(x['common_setup_cpu_seconds'] for x in design['repetitions'])
        initializer_CPU=initializer_warmup+sum(x['placement']['cpu_seconds'] for x in design['repetitions'])
        warmup_CPU=sum(x['record']['actual_local_CPU_seconds'] for x in design['warmups'])
        measured_optimizer_CPU=sum(x['actual_local_CPU_seconds'] for x in design['records'] if x['arm']!='current_placement')
        campaign_rows.append(dict(case=case,actual_common_setup_CPU_seconds=common_CPU,
            actual_initializer_CPU_seconds=initializer_CPU,actual_warmup_optimizer_CPU_seconds=warmup_CPU,
            actual_measured_optimizer_CPU_seconds=measured_optimizer_CPU,
            recorded_algorithm_CPU_seconds=common_CPU+initializer_CPU+warmup_CPU+measured_optimizer_CPU))
        for setup in design['repetitions']:
            warm = setup['placement']
            setup_rows.append(dict(case=case, repetition=setup['repetition'],
                                   common_setup_CPU_seconds=setup['common_setup_cpu_seconds'],
                                   placement_CPU_seconds=warm['placement_CPU_seconds'],
                                   pauli_conversion_CPU_seconds=warm['pauli_conversion_CPU_seconds'],
                                   actual_initializer_CPU_seconds=warm['cpu_seconds']))
        for index, record in enumerate(design['records']):
            lookup[case, record['arm'], record['CPU_budget'], record['repetition']] = (index, record)
            row = dict(case=case, modes=model['n'], Pauli_words=len(model['coefficients']),
                       conserving_units=len(model['blocks']), record_index=index)
            for field in ('arm','repetition','CPU_budget','partition','objective','initialization','status',
                          'total_charged_CPU_seconds','common_setup_CPU_charge','inherited_placement_CPU_charge',
                          'initializer_CPU_seconds','setup_CPU_seconds','search_CPU_seconds','actual_local_CPU_seconds',
                          'accounting_overhead_CPU_seconds','CPU_overrun_seconds','within_CPU_budget',
                          'evaluation_cap','evaluations','setup_pairs_visited','nonzero_pair_count','destinations',
                          'array_payload_bytes','storage_upper_bound_bytes','rescore_reservation_violated'):
                row[field] = record.get(field)
            row.update(starts=len(record.get('starts',[])),
                       sweeps=sum(len(x['sweeps']) for x in record.get('starts',[])),
                       incomplete_sweeps=sum(not s['complete'] for x in record.get('starts',[]) for s in x['sweeps']),
                       accepted_moves=sum(len(x['accepted']) for x in record.get('starts',[])),
                       max_full_rescore_CPU_seconds=max(record.get('full_rescore_CPU_seconds',[]),default=None),
                       peak_RSS_KiB_at_input_end=design.get('process_peak_RSS_KiB_at_input_end'))
            searches.append(row)
    assert len(searches) == len(lookup) == 432

    def pair(case, rep, formula, left_key, right_key, comparison, **factors):
        li, left = lookup[(case,)+left_key+(rep,)]
        ri, right = lookup[(case,)+right_key+(rep,)]
        a, b = physical[case,li,formula]['error'], physical[case,ri,formula]['error']
        return dict(case=case, repetition=rep, formula=formula, comparison=comparison, **factors,
                    left_arm=left['arm'], right_arm=right['arm'], left_budget=left['CPU_budget'], right_budget=right['CPU_budget'],
                    left_record=li,right_record=ri,left_error=a,right_error=b,error_difference=a-b,
                    classification='win' if a-b < -1e-12 else 'loss' if a-b > 1e-12 else 'tie',
                    denominator_above_floor=b>1e-12,error_ratio=a/b if b>1e-12 else None,
                    left_charged_CPU=left['total_charged_CPU_seconds'],right_charged_CPU=right['total_charged_CPU_seconds'],
                    left_within_CPU_budget=left.get('within_CPU_budget'),right_within_CPU_budget=right.get('within_CPU_budget'),
                    both_within_declared_cap=bool(left.get('within_CPU_budget') and right.get('within_CPU_budget')),
                    left_status=left['status'],right_status=right['status'])

    comparisons = []
    for case, rep, formula, cap in itertools.product(designs,range(3),(1,2),(.5,2.)):
        for objective, initial in itertools.product(('JD','J0'),('cold','placement')):
            comparisons.append(pair(case,rep,formula,('conserving_'+objective+'_'+initial,cap),('pauli_'+objective+'_'+initial,cap),
                                    'conserving_vs_pauli',objective=objective,initialization=initial,CPU_cap=cap))
        for partition, objective in itertools.product(('conserving','pauli'),('JD','J0')):
            comparisons.append(pair(case,rep,formula,(partition+'_'+objective+'_placement',cap),(partition+'_'+objective+'_cold',cap),
                                    'placement_vs_cold',partition=partition,objective=objective,CPU_cap=cap))
        for partition, initial in itertools.product(('conserving','pauli'),('cold','placement')):
            comparisons.append(pair(case,rep,formula,(partition+'_JD_'+initial,cap),(partition+'_J0_'+initial,cap),
                                    'JD_vs_J0',partition=partition,initialization=initial,CPU_cap=cap))
    assert len(comparisons) == 1152
    # Strong controls include a different formula. Keep that distinction explicit.
    controls = []
    for case, design in designs.items():
        for index, record in enumerate(design['records']):
            for formula, control_formula, (control, cap) in itertools.product((1,2),(1,2),
                    (('current_placement',None),('pauli_JD_larger_CPU',8.))):
                control_index, control_record = lookup[case,control,cap,record['repetition']]
                a = physical[case,index,formula]['error']
                b = physical[case,control_index,control_formula]['error']
                controls.append(dict(case=case,record_index=index,repetition=record['repetition'],arm=record['arm'],
                                     CPU_budget=record['CPU_budget'],formula=formula,control=control,control_formula=control_formula,
                                     error=a,control_error=b,error_difference=a-b,
                                     classification='win' if a-b < -1e-12 else 'loss' if a-b > 1e-12 else 'tie',
                                     error_ratio=a/b if b>1e-12 else None,denominator_above_floor=b>1e-12,
                                     same_formula=formula==control_formula,control_CPU_budget=cap))
    assert len(controls) == 3456
    aggregates = []
    for comparison, formula, cap in itertools.product(('conserving_vs_pauli','placement_vs_cold','JD_vs_J0'),(1,2),(.5,2.)):
        group = [p for p in comparisons if p['comparison']==comparison and p['formula']==formula and p['CPU_cap']==cap]
        ratios = [p['error_ratio'] for p in group if p['error_ratio'] is not None]
        aggregates.append(dict(comparison=comparison,formula=formula,CPU_cap=cap,pairs=len(group),
                               wins=sum(p['classification']=='win' for p in group),ties=sum(p['classification']=='tie' for p in group),
                               losses=sum(p['classification']=='loss' for p in group),floor_exclusions=len(group)-len(ratios),
                               median_of_paired_ratios=statistics.median(ratios) if ratios else None,
                               both_within_declared_cap=sum(p['both_within_declared_cap'] for p in group)))
    output.mkdir(parents=True)
    flat_rows = []
    for r in rows:
        flat = {k:v for k,v in r.items() if not isinstance(v,dict)}
        for name, values in r.items():
            if isinstance(values,dict):
                flat.update({name+'_'+k:v for k,v in values.items()})
        flat_rows.append(flat)
    for name, data in [('all_rows.csv',flat_rows),('searches.csv',searches),('shared_setup.csv',setup_rows),
                       ('paired_factorial.csv',comparisons),('strong_controls.csv',controls),('aggregate.csv',aggregates),('warmups.csv',warmup_rows),
                       ('campaign_construction.csv',campaign_rows)]:
        csv_write(output/name,data)
    over = sum(r['within_CPU_budget'] is False for r in searches)
    aggregate_table=['| Paired factor (left/right) | Formula | CPU cap | Wins/ties/losses | Median paired error ratio | Both within cap |',
        '| --- | --- | --- | --- | --- | --- |']
    for r in aggregates:
        aggregate_table.append('| '+' | '.join([r['comparison'],str(r['formula']),str(r['CPU_cap']),
            str(r['wins'])+'/'+str(r['ties'])+'/'+str(r['losses']),
            'floor-limited' if r['median_of_paired_ratios'] is None else format(r['median_of_paired_ratios'],'.6g'),
            str(r['both_within_declared_cap'])+'/'+str(r['pairs'])])+' |')
    aggregate_text='\n'.join(aggregate_table)
    first_order='; '.join(str(r['CPU_cap'])+' CPU seconds: '+str(r['wins'])+'/'+str(r['ties'])+'/'+str(r['losses'])
        for r in aggregates if r['comparison']=='conserving_vs_pauli' and r['formula']==1)
    (output/'README.md').write_text(f'''# Controlled search costs and declared synthetic scaling

This fixed panel retains 432 constructed orders and 864 full-space circuit comparisons: three retained molecular targets and five synthetic extended-Hubbard targets with 4–12 modes. The synthetic family extends the measured mode/term range; it supplies no larger-molecule validation. All repetitions and incomplete/out-of-budget records remain visible. {over} bounded-arm records exceeded their declared stand-alone CPU cap.

Conserving versus individual-Pauli first-order search has the following pooled wins/ties/losses: {first_order}. The complete table below includes both formulas/caps and every partition, initializer and objective comparison. These summaries pool the other declared factors and related repetitions; they do not establish a universal partition advantage or isolate conservation as a mechanism.

{aggregate_text}

Eight factorial arms separate conserving/individual-Pauli partition, J_D/J_0 objective and cold/placement initialization. Caps are 0.5 and 2 process-CPU seconds, with three repetitions, cyclic arm order and symmetric warmups. [Search records](searches.csv) count common setup, charged warm initializer/conversion, method setup, search, positions, partial sweeps, rescoring, storage and actual overruns. [Shared costs](shared_setup.csv) record actual common/initializer execution once per input/repetition; charging an inherited initializer to several standalone methods does not mean the campaign recomputed it each time. [All warmups](warmups.csv) retain actual optimizer work and CPU failures. [Actual campaign construction](campaign_construction.csv) counts shared setup and initializer once for the warmup phase and each measured repetition, then adds actual local optimizer CPU; it excludes file serialization, import overhead and original input generation, which remain in the full process resource ledger. Original input generation remains separate from stored-input search setup.

[Paired factorial comparisons](paired_factorial.csv) distinguish partition, objective and initializer effects at each input/repetition/budget/formula. Aggregate counts are descriptive summaries of related pairs, not independent statistical replicates. Ratios are formed per pair before aggregation, excluding and counting only denominator errors at or below 1e-12. Losses remain in every count. The both-within-cap field identifies failures to satisfy a nominal matched budget; such rows are retained and never credited as successful in-budget comparisons. Incomplete sweeps remain eligible executed policies only under their explicitly retained stopping records.

[Strong controls](strong_controls.csv) include separately timed current placement and the 8-CPU-second individual-Pauli search, with both first-order and symmetric reuse. The larger Pauli control has an unequal classical budget. Comparisons changing formula do not isolate ordering. Symmetric controls reuse first-order objectives, with exact central/boundary fusion and actual affordable steps at 10U ideal rotations; they do not complete second-order-specific optimization.

Full leakage-preserving propagation, independent literal-Pauli surrogate scoring, an independently scaled Taylor reference and unfused replay validate every retained circuit. Preparation, search storage and circuit cache limits are declared before execution. The process peak RSS at each input end is cumulative for the worker, not an independently reset per-input peak. Sparse array payload estimates omit Python/container overhead; the OS memory ceiling and measured RSS remain separate controls.

Conclusions are limited to this implementation, target family, budget/stopping policy and declared scaling range. Quarter-core/shared-machine CPU observations do not establish dedicated-machine wall speedups, an asymptotic advantage, conservation as the accuracy mechanism, compiled resources or hardware results. The rejected two-level correction is unchanged; no failed holdout gate is retuned. Full absolute errors and sector leakage are in [all rows](all_rows.csv).
''')
    write_manifest = dict(files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()},
                          schedule_freeze_sha256=sha(results/'SCHEDULES_FROZEN.json'),verification_sha256=sha(results/'VERIFIED.json'),
                          rows=864,orders=432,factorial_pairs=1152,strong_control_pairs=3456)
    (output/'MANIFEST.json').write_text(json.dumps(write_manifest,indent=2)+'\n')
    print(json.dumps(dict(status='passed',rows=864,orders=432,factorial_pairs=1152,strong_control_pairs=3456)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.results, args.output)
