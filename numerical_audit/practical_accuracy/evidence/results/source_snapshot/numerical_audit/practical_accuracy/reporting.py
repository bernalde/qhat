"""Deterministic complete tables from the retained full panel, with all losses."""
import argparse
from collections import Counter
import csv
import io
from pathlib import Path
import statistics
import common as k
import study


def csv_file(path, fields, rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(stream.getvalue())


def report(output, destination):
    k.require(not destination.exists(), 'Preserve reports; use a fresh destination')
    verified = k.read(output / 'VERIFIED.json')
    k.require(verified['status'] == 'passed' and verified['decision_agreement'], 'Independent audit required')
    destination.mkdir(parents=True)
    choices = k.read(output / 'SELECTIONS.json')
    pairs = k.read(output / 'COMPARISONS.json')
    rows = {p.stem: k.read(p) for p in sorted((output / 'rows').glob('*.json'))}
    candidate_rows = []
    for path in sorted((output / 'candidates').glob('*.json')):
        candidate = k.read(path)
        row = {key: candidate[key] for key in ('candidate', 'curve', 'case', 'method', 'cohort', 'repetition',
               'steps', 'status', 'reason', 'ordering_CPU_charge', 'compilation_CPU_seconds',
               'total_construction_CPU_seconds', 'total_CPU_budget', 'CPU_overrun_seconds', 'eligible')}
        cost = candidate['compiled_resources'] or candidate['partial_compiled_resources']
        row.update(cost)
        if candidate['candidate'] in rows:
            physical = rows[candidate['candidate']]
            row.update({key: physical[key] for key in ('error', 'infidelity', 'norm', 'number_leakage', 'joint_sector_leakage')})
        candidate_rows.append(row)
    csv_file(destination / 'candidates.csv', ['candidate', 'curve', 'case', 'method', 'cohort', 'repetition',
             'steps', 'status', 'reason', 'eligible', 'total_CPU_budget', 'ordering_CPU_charge',
             'compilation_CPU_seconds', 'total_construction_CPU_seconds', 'CPU_overrun_seconds',
             'two_qubit_gates', 'unit_duration_ASAP_depth', 'gates', 'single_qubit_gates', 'CX_dependency_depth',
             'error', 'infidelity', 'norm', 'number_leakage', 'joint_sector_leakage'], candidate_rows)
    target_rows = []
    for choice in choices:
        row = dict(choice)
        row['ambiguous_candidates'] = ';'.join(choice['ambiguous_candidates'])
        if choice['selected_candidate'] is not None:
            selected = rows[choice['selected_candidate']]
            row.update({key: selected[key] for key in ('error', 'steps', 'ordering_CPU_charge',
                       'compilation_CPU_seconds', 'total_construction_CPU_seconds')})
            row.update(selected['compiled_resources'])
        target_rows.append(row)
    csv_file(destination / 'targets.csv', ['curve', 'case', 'method', 'cohort', 'budget', 'repetition', 'target',
             'status', 'selected_candidate', 'qualifying_candidates', 'ambiguous_candidates', 'complete_candidates',
             'censored_candidates', 'construction_ineligible_candidates', 'steps', 'error', 'two_qubit_gates',
             'unit_duration_ASAP_depth', 'gates', 'ordering_CPU_charge', 'compilation_CPU_seconds',
             'total_construction_CPU_seconds'], target_rows)
    csv_file(destination / 'pairs.csv', ['case', 'budget', 'repetition', 'target', 'method', 'control',
             'control_cohort', 'method_status', 'control_status', 'two_qubit_gates', 'unit_duration_ASAP_depth',
             'method_selected', 'control_selected'], pairs)
    groups = []
    for target in k.TARGETS:
        for budget in k.BUDGETS:
            for method in ('cubic_JD', 'cubic_warm_JD'):
                for control in ('first_JD', 'pauli_JD', 'cubic_JD'):
                    panel = [p for p in pairs if (p['target'], p['budget'], p['method'], p['control']) ==
                             (target, budget, method, control)]
                    if panel:
                        groups.append(dict(target=target, budget=budget, method=method, control=control,
                                           paired_observations=len(panel),
                                           CX=dict(Counter(p['two_qubit_gates'] for p in panel)),
                                           depth=dict(Counter(p['unit_duration_ASAP_depth'] for p in panel))))
    summary = dict(curves=verified['curves'], candidates=verified['candidates'], complete=verified['complete'],
                   censored=verified['censored'], target_rows=verified['target_rows'], comparisons=verified['comparisons'],
                   target_status=dict(Counter(s['status'] for s in choices)),
                   matched_target_status=dict(Counter(s['status'] for s in choices if s['cohort'] == 'matched')),
                   complete_construction_overruns=sum(c['status'] == 'complete' and not c['eligible'] for c in candidate_rows),
                   groups=groups, actual_preparation=k.read(output / 'PREPARE_COSTS.json'),
                   actual_screening=k.read(output / 'RUN_COSTS.json'),
                   independent_CPU_seconds=verified['CPU_seconds'],
                   source_schedule_sha256=k.sha(output / 'SCHEDULES_FROZEN.json'))
    k.write(destination / 'SUMMARY.json', summary)
    text = ('# Practical accuracy under total classical construction caps\n\n'
            + ('Altered-model full-path preflight; not main-panel scientific evidence.\n\n'
               if k.read(output / 'PROTOCOL_FROZEN.json')['preflight'] else '')
            + f"All {summary['curves']} curves, {summary['candidates']} candidates, {summary['target_rows']} target rows "
              f"and {summary['comparisons']} paired decisions are retained. Complete/censored circuits: "
              f"{summary['complete']}/{summary['censored']}. Target status: {summary['target_status']}.\n\n"
            + 'The matched cohort charges full measured ordering plus one-circuit compilation through '
              'serialization/hash/cleanup. Frozen references retain inherited costs and are outside the '
              'matched-cap cohort. Repeated standalone charges are not campaign CPU. Full-grid screening '
              'and exact reference calculations are offline costs, reported separately.\n\n'
            + 'Candidates.csv includes every complete circuit and censored prefix; targets.csv includes '
              'every repetition and unmet/ambiguous target; pairs.csv retains every planned comparison '
              'including unavailable pairs. SUMMARY.json groups paired outcomes before aggregation.\n\n'
            + 'This is a three-model fixed-basis development panel at T=2.11 with HF input and one '
              'prototype exact-angle line compiler. It establishes neither transfer to new molecular '
              'families nor hardware execution performance, a causal conservation mechanism, global '
              'ordering/compiler optimality or efficient correlated-state preparation.\n')
    (destination / 'README.md').write_text(text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    report(args.output.resolve(), args.report.resolve())
