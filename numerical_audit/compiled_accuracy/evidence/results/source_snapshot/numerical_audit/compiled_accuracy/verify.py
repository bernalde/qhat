"""Offline independent replay, circuit DAG audit and accuracy-target selection.

The deterministic trace reproduction explicitly reuses the producing compiler;
it is distinct from literal sparse primitive/Pauli replay and independent depth.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import argparse,gzip,hashlib,itertools,json,math,shutil,tempfile,time
from pathlib import Path
import numpy as np
import study
from compiler import VERSION,write_trace
from reference import Reference,metrics,require

STATE_TOL=2e-10
METRIC_TOL=2e-12
ARMS=(('current',1),('first_JD_sameT',1),('pauli_JD_large_sameT',1),('current',2),
      ('first_JD_sameT',2),('pauli_JD_large_sameT',2),('cubic_orientation',2),('cubic_insertion',2))
GRID=(1,2,4,8,16,32,64,128)
TARGETS=(.01,.001)
AMBIGUITY=2e-12

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def canonical_trace_sha(path):
    h=hashlib.sha256()
    with gzip.open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def read(path):return json.loads(Path(path).read_text())
def inventory(root):
    paths=list(root.rglob('*'));require(not any(p.is_symlink() for p in paths),'Symlink refused')
    return {str(p.relative_to(root)) for p in paths if p.is_file()}
def nonnegative_costs(obj):
    for key,value in obj.items():
        if isinstance(value,dict):nonnegative_costs(value)
        elif key.endswith(('CPU_seconds','wall_seconds')):
            require(type(value) in (int,float) and math.isfinite(value) and value>=0,'Invalid retained phase cost: '+key)

def authenticate_upstream(output,protocol):
    root=output/'upstream';evidence=root/'evidence';manifest=read(evidence/'MANIFEST.json')
    require(inventory(evidence)==set(manifest['files'])|{'MANIFEST.json'},'Offline upstream inventory differs')
    require(sha(evidence/'MANIFEST.json')==protocol['upstream']['manifest_sha256'],'Upstream manifest identity differs')
    for name,digest in manifest['files'].items():require(sha(evidence/study.relative(name))==digest,'Offline upstream artifact differs')
    require(inventory(root/'source_snapshot')==set(manifest['sources']),'Offline upstream source inventory differs')
    for name,digest in manifest['sources'].items():require(sha(root/'source_snapshot'/study.relative(name))==digest,'Offline upstream source differs')
    require(read(root/'PROVENANCE.json')==protocol['upstream'],'Upstream provenance differs')
    prior=read(evidence/'results/PROTOCOL_FROZEN.json');verified=read(evidence/'results/VERIFIED.json')
    require(prior['preflight']==protocol['preflight'] and verified['status']=='passed' and verified['decision_agreement'],'Wrong or unverified upstream study')
    require(verified['rows']==(16 if protocol['preflight'] else 96) and verified['comparisons']==(56 if protocol['preflight'] else 336),'Incomplete upstream panel')
    require(verified['schedule_freeze_sha256']==sha(evidence/'results/SCHEDULES_FROZEN.json'),'Upstream schedule identity differs')
    require(protocol['upstream']['dependency_kind']==('development_preflight' if protocol['preflight'] else 'merged'),'Wrong upstream dependency class')
    require(protocol['upstream']['payload_Git_verified']==(not protocol['preflight']),'Wrong upstream Git authentication class')
    return evidence

def independent_choices(curves,candidates,rows,errors):
    choices=[]
    for curve in sorted(curves):
        c=curves[curve];panel=[candidates[curve+'.r'+format(r,'03d')] for r in GRID]
        for target in TARGETS:
            good=[];ambiguous=[]
            for candidate in panel:
                if candidate['status']!='complete':continue
                label=candidate['candidate'];error=errors[label]['error']
                if abs(error-target)<=AMBIGUITY:ambiguous.append(label)
                elif error<target:good.append(rows[label])
            winner=min(good,key=lambda row:(row['compiled_resources']['two_qubit_gates'],row['compiled_resources']['unit_duration_ASAP_depth'],row['compiled_resources']['gates'],row['steps'])) if good else None
            choices.append(dict(curve=curve,case=c['case'],T=c['T'],arm=c['arm'],formula=c['formula'],target=target,
                status='met' if winner else 'unmet',selected_candidate=winner['candidate'] if winner else None,
                ambiguous_candidates=ambiguous,censored_candidates=sum(x['status']=='censored' for x in panel),
                feasible_candidates=sum(x['status']=='complete' for x in panel),qualifying_candidates=len(good)))
    return choices

def verify(output):
    cpu=time.process_time();wall=time.monotonic();protocol,frozen=study.bind(output,independent=True)
    require(protocol['arms']==[list(x) for x in ARMS] and protocol['step_grid']==list(GRID) and protocol['accuracy_targets']==list(TARGETS),'Compiled comparison domain differs')
    require(protocol['threshold_ambiguity']==AMBIGUITY and protocol['training_rotation_cap_units']==30 and protocol['compiler']==VERSION,'Compiler or accuracy policy differs')
    expected_cases=['H2O_scale1p10'] if protocol['preflight'] else ['Li2_r2p95','F2_r1p55','H2O_scale1p10']
    expected_times=[2.11] if protocol['preflight'] else [.83,2.11]
    require(protocol['cases']==expected_cases and protocol['times']==expected_times,'Prospective model/time domain differs')
    require(protocol['compiler_sha256']==sha(study.HERE/'compiler.py'),'Trace reproduction requires the unchanged producing compiler')
    require(protocol['producer_sources']['numerical_audit/compiled_accuracy/study.py']==sha(study.HERE/'study.py'),'Trace reproduction requires unchanged producing controller')
    evidence=authenticate_upstream(output,protocol)
    cases=[('preflight_altered_' if protocol['preflight'] else '')+c for c in expected_cases]
    require(inventory(output/'models')=={c+'.json' for c in cases},'Compiled model inventory differs')
    expected_curve_names={c+'.T'+format(T,'.2f')+'.'+arm+'.f'+str(formula) for c,T,(arm,formula) in itertools.product(cases,expected_times,ARMS)}
    require(inventory(output/'curves')=={c+'.json' for c in expected_curve_names},'Curve inventory differs')
    curves={};models={};references={}
    for case in cases:
        model=read(output/'models'/(case+'.json'));require(model==read(evidence/'results/inputs'/(case+'.json')),'Compiled target differs from ordering target')
        models[case]=model;references[case]=Reference(model)
        design_path=evidence/'results/designs'/(case+'.json');design=read(design_path)
        for T,(arm,formula) in itertools.product(expected_times,ARMS):
            name=case+'.T'+format(T,'.2f')+'.'+arm+'.f'+str(formula);curve=read(output/'curves'/(name+'.json'))
            expected=[r for r in design['records'] if (r['T'],r['budget_units'],r['arm'],r['formula'])==(T,30,arm,formula)]
            require(len(expected)==1 and curve['record']==expected[0],'Order is not the fixed upstream30U record')
            require((curve['curve'],curve['case'],curve['T'],curve['arm'],curve['formula'])==(name,case,T,arm,formula),'Curve labels differ')
            require(curve['upstream_design_sha256']==sha(design_path) and curve['model_sha256']==sha(output/'models'/(case+'.json')),'Curve provenance differs')
            study.validate_record(model,curve['record']);curves[name]=curve
    expected_candidates={curve+'.r'+format(r,'03d') for curve,r in itertools.product(curves,GRID)}
    require(inventory(output/'candidates')=={c+'.json' for c in expected_candidates},'Candidate grid differs')
    candidates={name:read(output/'candidates'/(name+'.json')) for name in sorted(expected_candidates)}
    feasible={name for name,c in candidates.items() if c['status']=='complete'}
    require(all(c['status'] in ('complete','censored') for c in candidates.values()),'Unknown candidate termination')
    require(inventory(output/'rows')=={c+'.json' for c in feasible} and inventory(output/'states')=={c+'.npz' for c in feasible},'Physical grid missing or contains censored measurements')
    expected_traces={c+('.ndjson.gz' if c in feasible else '.partial.ndjson.gz') for c in candidates}
    require(inventory(output/'traces')==expected_traces,'Exact circuit inventory differs')
    require(frozen['candidates']==len(candidates)==protocol['expected_candidates'] and frozen['feasible']==len(feasible) and frozen['censored']==len(candidates)-len(feasible),'Compilation totals differ')
    require(len(curves)==protocol['expected_curves'],'Curve total differs')
    rows={};errors={};checks=[];maximum_state=maximum_metric=0.;exact={}
    scratch=Path(tempfile.mkdtemp(prefix='compiled-verification-',dir=output))
    try:
        for label,candidate in candidates.items():
            curve=curves[candidate['curve']];model=models[curve['case']];reference=references[curve['case']];record=curve['record']
            require(candidate['candidate']==label and candidate['steps'] in GRID and label==curve['curve']+'.r'+format(candidate['steps'],'03d'),'Candidate identity differs')
            for key in ('case','T','arm','formula'):require(candidate[key]==curve[key],'Candidate cell differs')
            identity=dict(candidate=label,curve=curve['curve'],case=curve['case'],T=curve['T'],T_ratio=list(float(curve['T']).as_integer_ratio()),steps=candidate['steps'],
                formula=curve['formula'],arm=curve['arm'],compilation_status=candidate['status'],model_sha256=curve['model_sha256'],curve_sha256=sha(output/'curves'/(curve['curve']+'.json')),protocol_sha256=sha(output/'PROTOCOL_FROZEN.json'))
            require(candidate['identity']==identity and candidate['inherited_ordering_CPU_charge']==record['standalone_construction_CPU_charge'],'Candidate provenance/cost attribution differs')
            nonnegative_costs(candidate)
            pieces=candidate['lowering_CPU_seconds']+candidate['bookkeeping_CPU_seconds']+sum(candidate[k]['CPU_seconds'] for k in ('resource_counting','serialization','hashing'))
            require(abs(pieces-candidate['actual_local_constructor_CPU_seconds'])<1e-9,'Constructor cost closure differs')
            complete=label in feasible;prefix='' if complete else 'partial_';trace=output/candidate[prefix+'trace_file']
            require(trace==output/'traces'/(label+('.ndjson.gz' if complete else '.partial.ndjson.gz')) and sha(trace)==candidate[prefix+'trace_sha256'] and trace.stat().st_size==candidate[prefix+'trace_bytes'],'Trace binding differs')
            state,costs=reference.audit_compiled_trace(trace,identity,VERSION,simulate=complete)
            require(costs==candidate['compiled_resources' if complete else 'partial_compiled_resources'],'Independent DAG/count audit differs')
            reproduced,termination=study.compile_candidate(model,record,candidate['steps'])
            require(termination['status']==candidate['status'] and termination['reason']==candidate['reason'] and termination['generated']==candidate['generated'] and termination['stored']==candidate['stored'],'Compiler deterministic termination differs')
            rotation_key='ideal_Pauli_rotations' if complete else 'completed_ideal_rotations'
            require(termination[rotation_key]==candidate[rotation_key],'Ideal-rotation count differs')
            replay_path=scratch/'reproduced.ndjson.gz';write_trace(replay_path,model['n'],reproduced.gates,identity)
            require(canonical_trace_sha(replay_path)==canonical_trace_sha(trace),'Deterministic exact canonical trace reproduction differs')
            replay_path.unlink();del reproduced
            if complete:
                row=read(output/'rows'/(label+'.json'));rows[label]=row
                for key in ('candidate','curve','case','T','arm','formula','steps','ideal_Pauli_rotations','compiled_resources','inherited_ordering_CPU_charge'):
                    require(row[key]==candidate[key],'Physical row attribution differs')
                require(row['candidate_sha256']==sha(output/'candidates'/(label+'.json')) and row['trace_sha256']==sha(trace),'Physical row source differs')
                state_path=output/study.relative(row['state_file']);require(state_path==output/'states'/(label+'.npz') and sha(state_path)==row['state_sha256'],'Physical state identity differs')
                with np.load(state_path,allow_pickle=False) as stored:
                    require(set(stored.files)=={'actual'},'Unexpected physical state payload');producer=stored['actual'].copy()
                ideal=reference.ideal_formula(record,candidate['steps'],curve['T'])
                state_delta=max(float(np.linalg.norm(state-producer)),float(np.linalg.norm(state-ideal)))
                require(state_delta<=STATE_TOL,'Compiled/literal/unfused state disagreement')
                key=(curve['case'],curve['T'])
                if key not in exact:exact[key]=reference.exact_action(curve['T'])
                independent=metrics(exact[key],state,model);ideal_metrics=metrics(exact[key],ideal,model)
                metric_delta=max(max(abs(independent[k]-row[k]),abs(independent[k]-ideal_metrics[k])) for k in independent)
                require(metric_delta<=METRIC_TOL,'Independent physical metric disagreement')
                errors[label]=independent;maximum_state=max(maximum_state,state_delta);maximum_metric=max(maximum_metric,metric_delta)
                checks.append(dict(candidate=label,state_delta=state_delta,metric_delta=metric_delta,**independent))
            print('VERIFIED_COMPILED',label,flush=True)
    finally:shutil.rmtree(scratch)
    choices=independent_choices(curves,candidates,rows,errors)
    require(choices==read(output/'SELECTIONS.json'),'Independent target selection/ambiguity/censor decisions differ')
    require(read(output/'RUN_COMPLETE.json')==dict(candidates=len(candidates),feasible=len(feasible),censored=len(candidates)-len(feasible),target_rows=2*len(curves)),'Run completion record differs')
    for file in ('PREPARE_COSTS.json','RUN_COSTS.json'):nonnegative_costs(read(output/file))
    study.write(output/'INDEPENDENT_ERRORS.json',checks);study.write(output/'INDEPENDENT_SELECTIONS.json',choices)
    result=dict(status='passed',candidates=len(candidates),feasible=len(feasible),censored=len(candidates)-len(feasible),curves=len(curves),target_decisions=len(choices),
        decision_agreement=True,max_state_delta=maximum_state,max_metric_delta=maximum_metric,schedule_freeze_sha256=sha(output/'SCHEDULES_FROZEN.json'),
        artifacts={str(p.relative_to(output)):sha(p) for folder in ('rows','states') for p in sorted((output/folder).iterdir())} |
            {n:sha(output/n) for n in ('SELECTIONS.json','RUN_COMPLETE.json','RUN_COSTS.json')},
        producing_source_commit=protocol['source_commit'],producing_sources=protocol['producer_sources'],verification_sources=study.verifier_sources(),
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        scope='Independent literal primitive and unfused Pauli replay, scaled Taylor target action, two-pass dependency DAG costs and independent all-grid target selection; deterministic exact canonical NDJSON reproduction reuses the producer compiler, with retained gzip bytes separately integrity-authenticated; no hardware execution')
    study.write(output/'VERIFIED.json',result);print(json.dumps(result))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    verify(parser.parse_args().output.resolve())
