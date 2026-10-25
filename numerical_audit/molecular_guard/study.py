"""Freeze the fixed LiH comparison, then retain every finite-step enclosure.

This producer does not supply the independent audit or a Lean molecular proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from intervals import Complex, require
import exact as ex
import propagation as prop

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODEL_SHA = 'ea38b9730eed2faf55fbd1fe22722d93b858b95fe4d804ecb72855d10176c134'
TIME = F(3,5)
EPSILONS = [F(1,1 << k) for k in range(0,21,2)]
STEPS = [1,2,4,8,16,32]
BITS,TERMS = 192,96


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, obj):
    path = Path(path)
    require(not path.exists(), 'Retain existing attempt artifacts: '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.partial')
    require(not temporary.exists(), 'Retain interrupted write: '+str(temporary))
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def source_paths():
    return sorted([p for p in HERE.iterdir() if p.suffix in ('.py','.md')] +
        [ROOT/'numerical_audit/prospective_reference'/name for name in ('engine.py','campaign.py','one_pass.py')])


def sources():
    return {str(path.relative_to(ROOT)):sha(path) for path in source_paths()}


def committed_sources():
    head = subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    for name in sources():
        committed = subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+name])
        require(committed == (ROOT/name).read_bytes(), 'Commit scientific sources before molecular freeze')
    return head


def legacy_module(name, path):
    if name in sys.modules:
        require(Path(sys.modules[name].__file__).resolve() == path.resolve(), 'Wrong legacy module '+name)
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def legacy_current(data, scale):
    folder = ROOT/'numerical_audit/prospective_reference'
    legacy_module('engine',folder/'engine.py')
    legacy_module('one_pass',folder/'one_pass.py')
    campaign = legacy_module('molecular_legacy_campaign',folder/'campaign.py')
    blocks = [{w:c*scale.numerator for w,c in block.items()} for block in data['blocks']]
    return campaign.current_placement(blocks,data['diag'],data['off'],data['hf'],data['den']*scale.denominator)


def freeze(model_path, output, *, preflight=False):
    head = committed_sources()
    require(not output.exists(), 'Use a fresh molecular output directory')
    require(sha(model_path) == MODEL_SHA, 'Wrong predeclared LiH model')
    model = read(model_path)
    require((model['case'],model['n'],model['electrons']) == ('LiH_r1p713',6,2), 'Wrong molecular input')
    scale = F(7,5) if preflight else F(1)
    data = ex.from_model(model,scale)
    cpu,wall = time.process_time(),time.monotonic()
    current = legacy_current(data,scale)
    construction = dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
                        boundary='Legacy current-placement imports/setup/search through return;input reconstruction and packaging separate')
    orders = dict(current=current['order'],diagonal_prefix=data['diag']+data['off'])
    for order in orders.values():
        require(sorted(order) == list(range(len(data['matrices']))), 'Incomplete frozen order')
    output.mkdir(parents=True)
    shutil.copyfile(model_path,output/'model.json')
    for path in source_paths():
        destination = output/'source_snapshot'/path.relative_to(ROOT)
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,destination)
    write(output/'FROZEN.json',dict(source_commit=head,source_sha256=sources(),model_sha256=MODEL_SHA,
        mode='preflight' if preflight else 'main',coefficient_scale=str(scale),time=str(TIME),
        epsilons=[str(x) for x in ([F(1),F(1,256)] if preflight else EPSILONS)],
        steps=[2] if preflight else STEPS,bits=BITS,terms=TERMS,orders=orders,
        legacy_current=current,construction_cost=construction,sector=data['sector'],hf=data['hf'],
        diagonal_ids=data['diag'],off_diagonal_ids=data['off'],
        exact_unit_matrices=[ex.encoded(a) for a in data['matrices']],
        full_space_invariance_columns_checked=(1 << data['n'])*len(data['matrices']),
        environment=dict(python=sys.version,executable=str(Path(sys.executable).resolve()),platform=platform.platform(),
                         packages={'numpy':importlib.metadata.version('numpy')}),
        dependency_status='Pinned development input;actual issue17 merge required before publication',
        product_convention='Written left-to-right exponentials;reverse chronological unit application'))
    print('FROZEN',output,flush=True)


def load_frozen(output):
    frozen = read(output/'FROZEN.json')
    require(frozen['source_sha256'] == sources(), 'Producer source changed after freeze')
    for name,digest in frozen['source_sha256'].items():
        require(sha(output/'source_snapshot'/name) == digest, 'Frozen source archive changed')
    require(sha(output/'model.json') == frozen['model_sha256'] == MODEL_SHA, 'Frozen input changed')
    require(frozen['time'] == str(TIME) and (frozen['bits'],frozen['terms']) == (BITS,TERMS), 'Changed time/precision')
    preflight = frozen['mode'] == 'preflight'
    require(frozen['mode'] in ('preflight','main'), 'Unknown run mode')
    require(frozen['coefficient_scale'] == str(F(7,5) if preflight else F(1)) and
            frozen['epsilons'] == [str(x) for x in ([F(1),F(1,256)] if preflight else EPSILONS)] and
            frozen['steps'] == ([2] if preflight else STEPS), 'Panel changed after freeze')
    data = ex.from_model(read(output/'model.json'),F(frozen['coefficient_scale']))
    require(frozen['exact_unit_matrices'] == [ex.encoded(a) for a in data['matrices']] and
            frozen['sector'] == data['sector'] and frozen['diagonal_ids'] == data['diag'] and
            frozen['off_diagonal_ids'] == data['off'], 'Frozen action identities changed')
    require(frozen['orders']['diagonal_prefix'] == data['diag']+data['off'] and
            frozen['orders']['current'] == frozen['legacy_current']['order'], 'Frozen order changed')
    for order in frozen['orders'].values():
        require(sorted(order) == list(range(len(data['matrices']))), 'Invalid frozen order')
    return frozen,data


def strict_winner(intervals):
    labels = list(intervals)
    require(len(labels) == 2, 'Two fixed schedules required')
    a,b = labels
    if intervals[a].hi < intervals[b].lo:
        return a
    if intervals[b].hi < intervals[a].lo:
        return b
    return 'unresolved'


def evaluate(output):
    cpu,wall = time.process_time(),time.monotonic()
    frozen,data = load_frozen(output)
    write(output/'EVALUATION_STARTED.json',dict(freeze_sha256=sha(output/'FROZEN.json'),
                                               source_sha256=sources()))
    leading = {}
    slopes = {}
    for label,order in frozen['orders'].items():
        k1,k2 = ex.bch(data,order)
        a,c,details = prop.leading_response(data,k1,k2,TIME,bits=BITS)
        leading[label] = dict(**details,imaginary_K1=ex.encoded(k1),imaginary_K2=ex.encoded(k2))
        slopes[label] = a,c
    write(output/'leading_response.json',leading)
    initial = [Complex.exact(int(j==data['hf_index']),bits=BITS) for j in range(len(data['sector']))]
    guards,comparisons,rows = [],[],[]
    for at,epsilon_string in enumerate(frozen['epsilons']):
        epsilon = F(epsilon_string)
        guard = {label:prop.response_guard(a,c,epsilon) for label,(a,c) in slopes.items()}
        guards.append(dict(epsilon=epsilon_string,intervals={k:v.witness() for k,v in guard.items()},
                           strict_smaller_response=strict_winner(guard)))
        matrices = [ex.times(h,1 if i in data['diag'] else epsilon) for i,h in enumerate(data['matrices'])]
        h = ex.zero(len(data['sector']))
        for matrix in matrices:
            h = ex.add(h,matrix)
        rcpu,rwall = time.process_time(),time.monotonic()
        reference_operator,reference_radius,reference_details = prop.exp_matrix(h,TIME,terms=TERMS,bits=BITS)
        reference = prop.matvec(reference_operator,initial)
        refname = f'references/epsilon-{at:02d}.json'
        write(output/refname,dict(epsilon=epsilon_string,state=[x.witness() for x in reference],
            analytic_radius=str(reference_radius),exponential=reference_details,
            CPU_seconds=time.process_time()-rcpu,wall_seconds=time.monotonic()-rwall))
        for steps in frozen['steps']:
            operators = [prop.exp_unit(matrix,TIME/steps,terms=TERMS,bits=BITS) for matrix in matrices]
            intervals = {}
            for label,order in frozen['orders'].items():
                scpu,swall = time.process_time(),time.monotonic()
                state,radius = initial,F(0)
                for _ in range(steps):
                    for index in reversed(order):
                        operator,tail = operators[index]
                        state = prop.matvec(operator,state)
                        radius = prop.combine(radius,tail,BITS)
                error,details = prop.infidelity(reference,reference_radius,state,radius)
                name = f'rows/epsilon-{at:02d}-r{steps:02d}-{label}.json'
                write(output/name,dict(epsilon=epsilon_string,steps=steps,order=label,
                    squared_infidelity=error.witness(),state=[x.witness() for x in state],
                    state_analytic_radius=str(radius),**details,
                    CPU_seconds=time.process_time()-scpu,wall_seconds=time.monotonic()-swall))
                rows.append(name)
                intervals[label] = error
            comparisons.append(dict(epsilon=epsilon_string,steps=steps,strict_smaller_error=strict_winner(intervals)))
            print('ENCLOSED',epsilon_string,'r',steps,{k:(str(v.bounds()[0]),str(v.bounds()[1])) for k,v in intervals.items()},flush=True)
    expected = len(frozen['epsilons'])*len(frozen['steps'])*2
    require(len(rows) == expected, 'Incomplete predeclared panel')
    write(output/'PRODUCER_COMPLETE.json',dict(status='producer_complete_independent_audit_pending',
        freeze_sha256=sha(output/'FROZEN.json'),circuits=len(rows),reference_evolutions=len(frozen['epsilons']),
        rows=rows,response_guards=guards,finite_rankings=comparisons,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        cost_boundary='Complete evaluation incl reconstruction,BCH,reference/unit polynomial construction,propagation and writes',
        certificate_scope='Outward software intervals;not independently audited or imported into Lean'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=('freeze','evaluate'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--model',type=Path)
    parser.add_argument('--preflight',action='store_true')
    args = parser.parse_args()
    if args.phase == 'freeze':
        require(args.model is not None,'Input model required')
        freeze(args.model,args.output,preflight=args.preflight)
    else:
        require(args.model is None and not args.preflight,'Evaluate the frozen panel without replacements')
        evaluate(args.output)


if __name__ == '__main__':
    main()
