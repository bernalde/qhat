"""Complete resource/error grid, target decisions and classical-cost accounting."""
import argparse,csv,itertools,json,statistics
from pathlib import Path
import study
from compiler import require

ARMS=study.ARMS

def csvwrite(path,rows):
    require(rows,'Empty report table')
    fields=sorted({key for row in rows for key in row})
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)

def verdict(left,right):return 'win' if left<right else 'loss' if left>right else 'tie'

def report(output,destination):
    protocol,frozen=study.bind(output,independent=True);verified=study.read(output/'VERIFIED.json')
    require(verified['status']=='passed' and verified['decision_agreement'] and verified['schedule_freeze_sha256']==study.sha(output/'SCHEDULES_FROZEN.json'),'Unverified or stale compiled grid')
    require(verified['verification_sources']==study.verifier_sources(),'Current independent verifier differs')
    for name,digest in verified['artifacts'].items():require(study.sha(output/study.relative(name))==digest,'Verified physical artifact changed')
    require(not destination.exists(),'Fresh report directory required');destination.mkdir(parents=True)
    candidates={p.stem:study.read(p) for p in sorted((output/'candidates').glob('*.json'))}
    physical={p.stem:study.read(p) for p in sorted((output/'rows').glob('*.json'))}
    choices=study.read(output/'SELECTIONS.json');independent=study.read(output/'INDEPENDENT_SELECTIONS.json');require(choices==independent,'Target decision agreement differs')
    points=[];constructor=[];selected=[];curves=[]
    for label,c in candidates.items():
        common={k:c[k] for k in ('candidate','curve','case','T','arm','formula','steps','status','reason','generated','stored','inherited_ordering_CPU_charge')}
        row=dict(common,ideal_Pauli_rotations=c.get('ideal_Pauli_rotations'),completed_ideal_rotations=c.get('completed_ideal_rotations'),
            actual_local_constructor_CPU_seconds=c['actual_local_constructor_CPU_seconds'],constructor_wall_seconds=c['constructor_wall_seconds'])
        for key in ('gates','single_qubit_gates','two_qubit_gates','unit_duration_ASAP_depth','CX_dependency_depth'):
            row[key]=c['compiled_resources'][key] if c['status']=='complete' else None
            row['partial_'+key]=c['partial_compiled_resources'][key] if c['status']=='censored' else None
        if label in physical:
            row.update({k:physical[label][k] for k in ('error','infidelity','norm','number_leakage','joint_sector_leakage')})
            for target in study.TARGETS:row['target_'+str(target)+'_classification']=study.threshold(physical[label]['error'],target)
        else:
            for target in study.TARGETS:row['target_'+str(target)+'_classification']='compiler_censored'
        points.append(row)
        constructor.append(dict(common,lowering_CPU_seconds=c['lowering_CPU_seconds'],lowering_wall_seconds=c['lowering_wall_seconds'],
            resource_counting_CPU_seconds=c['resource_counting']['CPU_seconds'],resource_counting_wall_seconds=c['resource_counting']['wall_seconds'],
            serialization_CPU_seconds=c['serialization']['CPU_seconds'],serialization_wall_seconds=c['serialization']['wall_seconds'],
            hashing_CPU_seconds=c['hashing']['CPU_seconds'],hashing_wall_seconds=c['hashing']['wall_seconds'],
            bookkeeping_CPU_seconds=c['bookkeeping_CPU_seconds'],actual_local_constructor_CPU_seconds=c['actual_local_constructor_CPU_seconds'],constructor_wall_seconds=c['constructor_wall_seconds'],
            retained_trace_bytes=c.get('trace_bytes',c.get('partial_trace_bytes'))))
    for choice in choices:
        row={k:v for k,v in choice.items() if k!='ambiguous_candidates'}
        row['ambiguous_candidates']=json.dumps(choice['ambiguous_candidates'],separators=(',',':'))
        if choice['status']=='met':
            label=choice['selected_candidate'];p=physical[label];c=candidates[label]
            row.update(steps=p['steps'],ideal_Pauli_rotations=p['ideal_Pauli_rotations'],**p['compiled_resources'],error=p['error'],infidelity=p['infidelity'],
                number_leakage=p['number_leakage'],joint_sector_leakage=p['joint_sector_leakage'],
                actual_local_constructor_CPU_seconds=c['actual_local_constructor_CPU_seconds'],inherited_ordering_CPU_charge=c['inherited_ordering_CPU_charge'],
                standalone_ordering_plus_one_circuit_CPU_charge=c['inherited_ordering_CPU_charge']+c['actual_local_constructor_CPU_seconds'])
        selected.append(row)
    grouped={}
    for row in selected:grouped.setdefault((row['case'],row['T'],row['target']),{})[(row['arm'],row['formula'])]=row
    comparisons=[]
    for (case,T,target),cell in sorted(grouped.items()):
        require(set(cell)==set(ARMS),'Incomplete target controls')
        for i,j in itertools.combinations(range(8),2):
            left,right=cell[ARMS[i]],cell[ARMS[j]];both=left['status']==right['status']=='met'
            status='both_met' if both else 'left_only_met' if left['status']=='met' else 'right_only_met' if right['status']=='met' else 'both_unmet'
            row=dict(case=case,T=T,target=target,left_index=i,right_index=j,left_arm=ARMS[i][0],left_formula=ARMS[i][1],right_arm=ARMS[j][0],right_formula=ARMS[j][1],
                status=status,left_selected=left['selected_candidate'],right_selected=right['selected_candidate'],
                left_censored=left['censored_candidates'],right_censored=right['censored_candidates'],left_ambiguous=left['ambiguous_candidates'],right_ambiguous=right['ambiguous_candidates'])
            if both:
                for key in ('two_qubit_gates','unit_duration_ASAP_depth','gates','CX_dependency_depth','ideal_Pauli_rotations'):
                    row['left_'+key]=left[key];row['right_'+key]=right[key];row[key+'_difference']=left[key]-right[key];row[key+'_result']=verdict(left[key],right[key])
                row.update(left_error=left['error'],right_error=right['error'],left_constructor_CPU=left['actual_local_constructor_CPU_seconds'],right_constructor_CPU=right['actual_local_constructor_CPU_seconds'])
            comparisons.append(row)
    summaries=[]
    for target,(arm,formula) in itertools.product(study.TARGETS,ARMS[:-1]):
        pairs=[(cell[ARMS[-1]],cell[arm,formula]) for (case,T,t),cell in sorted(grouped.items()) if t==target]
        both=[(a,b) for a,b in pairs if a['status']==b['status']=='met']
        row=dict(target=target,control=arm,control_formula=formula,pairs=len(pairs),both_met=len(both),
            cubic_only_met=sum(a['status']=='met' and b['status']!='met' for a,b in pairs),control_only_met=sum(a['status']!='met' and b['status']=='met' for a,b in pairs),both_unmet=sum(a['status']!= 'met' and b['status']!='met' for a,b in pairs))
        for key in ('two_qubit_gates','unit_duration_ASAP_depth'):
            for label,plural in (('win','wins'),('tie','ties'),('loss','losses')):row[key+'_'+plural]=sum(verdict(a[key],b[key])==label for a,b in both)
            ratios=[a[key]/b[key] for a,b in both if b[key]>0]
            row[key+'_zero_denominator_exclusions']=len(both)-len(ratios);row[key+'_median_paired_ratio']=statistics.median(ratios) if ratios else None
        summaries.append(row)
    for path in sorted((output/'curves').glob('*.json')):
        curve=study.read(path);panel=[c for c in candidates.values() if c['curve']==curve['curve']]
        curves.append(dict(curve=curve['curve'],case=curve['case'],T=curve['T'],arm=curve['arm'],formula=curve['formula'],
            candidates=len(panel),feasible=sum(c['status']=='complete' for c in panel),censored=sum(c['status']=='censored' for c in panel),
            full_grid_constructor_CPU_seconds=sum(c['actual_local_constructor_CPU_seconds'] for c in panel),
            inherited_standalone_ordering_CPU_charge=curve['record']['standalone_construction_CPU_charge']))
    expected_comparisons=56 if protocol['preflight'] else 336
    require(len(points)==protocol['expected_candidates'] and len(choices)==2*protocol['expected_curves'] and len(comparisons)==expected_comparisons,'Incomplete compiled comparison report')
    for name,rows in [('all_candidates.csv',points),('constructors.csv',constructor),('targets.csv',selected),('comparisons.csv',comparisons),('summary.csv',summaries),('curves.csv',curves)]:csvwrite(destination/name,rows)
    costs=dict(prepare=study.read(output/'PREPARE_COSTS.json'),physical_grid=study.read(output/'RUN_COSTS.json'),
        measured_all_candidate_constructor_CPU_seconds=sum(c['actual_local_constructor_CPU_seconds'] for c in candidates.values()),
        complete_candidates=len(physical),compiler_censored_candidates=len(candidates)-len(physical),
        note='Inherited standalone ordering charges and per-row shared reference timings repeat across alternatives. They are not summed as actual campaign CPU. The full grid construction and physical screening costs are separate from construction of one selected circuit. No hardware time measurement.')
    require(abs(costs['measured_all_candidate_constructor_CPU_seconds']-costs['prepare']['constructor_CPU_seconds'])<1e-8,'Campaign constructor accounting differs')
    study.write(destination/'COSTS.json',costs)
    lines=['# Compiled circuit resources at fixed accuracy','',
        ('This is a representative altered-model preflight, not the main experimental result.' if protocol['preflight'] else 'This is the complete prospectively frozen compiled-resource study.'),'',
        f"Compiler `{protocol['compiler']}` lowers to H/S/Sdg/Rz/CX on a nearest-neighbor line, with exact rational Rz angles, no ancillas, restored routing and only globally adjacent cancellation. This prototype supplies ideal compiled counts, not hardware execution or fault-tolerant synthesis.",'',
        f"All {len(curves)} order/formula curves and {len(candidates)} candidates at steps {list(study.GRID)} were frozen before compiled-state errors. Orders are the eight fixed 30U-trained arms from numerical issue8; no order is retuned on compiled resources or these errors. There are {len(physical)} feasible circuits and {len(candidates)-len(physical)} compiler-censored candidates.",'',
        'The projective state-error targets are0.01 and0.001. Values within2e-12 of a target remain ambiguous and cannot establish attainment. Each curve selects the qualifying grid point with minimum CX, then unit-duration ASAP depth, total gates and step count; error need not be monotone. Unmet and censored cases remain explicit. Exact primitive traces include HF preparation and all allowed boundary cancellations.','',
        '| Target | Cubic insertion versus | Formula | Both met | CX W/T/L | Depth W/T/L | Cubic-only met | Control-only met | Both unmet |',
        '| ---: | --- | ---: | ---: | --- | --- | ---: | ---: | ---: |']
    for row in summaries:
        wt=lambda key:'/'.join(str(row[key+'_'+x]) for x in ('wins','ties','losses'))
        lines.append(f"| {row['target']} | {row['control']} | {row['control_formula']} | {row['both_met']} | {wt('two_qubit_gates')} | {wt('unit_duration_ASAP_depth')} | {row['cubic_only_met']} | {row['control_only_met']} | {row['both_unmet']} |")
    lines+=['','`all_candidates.csv` retains the entire error-versus-resource grid; `targets.csv` retains every attainment decision and selected one-circuit cost; `comparisons.csv` retains every control pair. A CX tie is distinct from an ASAP-depth tie. Both-met count/depth outcomes do not silently turn an unmet curve into a resource win.','',
        '`constructors.csv`, `curves.csv` and `COSTS.json` separate lowering, routing/cancellation, resource counting, exact-trace serialization/hashing, bookkeeping, source/upstream authentication, shared model setup and full-grid physical screening. Inherited ordering charges are standalone scenarios from issue8, not newly executed searches or a sum of campaign work. Per-row reference timings repeat the same shared solve. One chosen circuit cost is not the cost of screening the entire grid for target attainment.','',
        f"Independent verification checks literal sparse primitive replay against unfused Pauli formulas and a separate scaled-Taylor target action, with state tolerance2e-10 and metric tolerance2e-12. Actual maximum differences: state{verified['max_state_delta']:.9g}, metric{verified['max_metric_delta']:.9g}. An independent two-pass dependency DAG checks counts/depth. Target decisions, ambiguity and censoring agree. Exact canonical trace reproduction reuses the pinned producer compiler and is labelled separately from the independent replay/DAG oracle.",'',
        'Results apply to these related finite molecular inputs, frozen times, exact compiler and supplied grid. They do not establish global compilation/search optimality, matched classical CPU, conservation as the mechanism, larger-molecule validity, hardware depth/timing, noise performance or a Lean proof.']
    (destination/'README.md').write_text('\n'.join(lines)+'\n')
    study.write(destination/'MANIFEST.json',dict(files={p.name:study.sha(p) for p in sorted(destination.iterdir()) if p.is_file()},
        verified_sha256=study.sha(output/'VERIFIED.json'),schedule_freeze_sha256=study.sha(output/'SCHEDULES_FROZEN.json'),scope='All candidates, all target decisions and all pairwise resource comparisons'))
    print(json.dumps(dict(status='passed',candidates=len(candidates),feasible=len(physical),target_decisions=len(choices),comparisons=len(comparisons))))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--report',type=Path,required=True)
    a=parser.parse_args();report(a.output.resolve(),a.report.resolve())
