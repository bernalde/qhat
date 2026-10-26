"""Report the complete fixed panel only after independent extracted verification.

All decisions use exact rational interval endpoints. Decimal approximations are
presentation columns and never select a schedule,epsilon,time or finite step.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LABELS = ('current','diagonal_prefix')
EPSILONS = [str(F(1,1 << k)) for k in range(0,21,2)]
STEPS = (1,2,4,8,16,32)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interval(raw):
    bits = raw['bits']
    require(type(bits) is int and bits == 192,'Expected producer192-bit enclosure')
    low,high = F(int(raw['lo']),1 << bits),F(int(raw['hi']),1 << bits)
    require(low <= high,'Reversed interval')
    return low,high


def winner(a, b):
    require(len(a) == len(b) == 2 and all(type(x) in (int,F) for x in (*a,*b)),
            'Exact rational interval endpoints required')
    return LABELS[0] if a[1] < b[0] else (LABELS[1] if b[1] < a[0] else 'unresolved')


def sufficient_epsilon_domain(a, b, c, d):
    """Positive-epsilon strict domain from fixed slope enclosures and constants.

    This is an algebraic consequence of the same fixed response bounds,not an
    additional validation panel or a replacement of any failed physical point.
    """
    require(len(a) == len(b) == 2 and all(type(x) in (int,F) for x in (*a,*b,c,d)),
            'Exact rational slope endpoints/constants required')
    a,b = tuple(map(F,a)),tuple(map(F,b))
    c,d = F(c),F(d)
    require(c >= 0 and d >= 0,'Nonnegative response remainder constants required')
    gap = b[0]-a[1]
    if gap <= 0:
        return dict(status='no_certified_slope_gap',gap_lower=str(gap))
    if c+d == 0:
        return dict(status='all_positive_epsilon_at_most_one',gap_lower=str(gap),
                    epsilon_upper='1',upper_inclusive=True)
    threshold = gap/(c+d)
    return dict(status='positive_epsilon_strict_upper',gap_lower=str(gap),
                epsilon_upper=str(min(F(1),threshold)),upper_inclusive=threshold > 1,
                raw_gap_over_constant_sum=str(threshold))


def validate_resources(record):
    ceilings = dict(tests=(120,900,64),preflight=(240,1800,256),main=(1200,7200,256),
                    independent=(1200,7200,256),reporting=(300,1800,256))
    require(set(record['phases']) == set(ceilings),'Incomplete numerical phase ledger')
    used_receipts = set()
    for phase,(cpu_cap,wall_cap,disk_cap) in ceilings.items():
        pool = record['phases'][phase]
        require((pool['CPU_seconds'],pool['wall_seconds'],pool['disk_MiB']) ==
                (cpu_cap,wall_cap,disk_cap),'Frozen phase ceiling changed')
        cpu,wall,disk = 0.,0.,0
        for attempt in pool['attempts']:
            name = Path(attempt['receipt']).name
            require(name not in used_receipts,'Duplicate charged resource receipt')
            used_receipts.add(name)
            entry = record['receipts'][name]
            require(hashlib.sha256(entry['text'].encode('utf8')).hexdigest() == entry['sha256'],
                    'Raw resource receipt hash differs')
            raw = json.loads(entry['text'])
            actual = raw['final_accounting']
            values = (actual['CPU_seconds'],actual['elapsed_seconds'])
            require(all(type(x) in (int,float) and math.isfinite(x) and x >= 0 for x in values),
                    'Invalid raw resource accounting')
            require(abs(attempt['CPU_seconds']-values[0]) < 1e-8 and
                    abs(attempt['wall_seconds']-values[1]) < 1e-8 and
                    attempt['exit_code'] == raw['exit_code'],'Attempt does not match raw receipt')
            budget = raw['budget']
            require(budget['cpu_seconds'] <= cpu_cap-cpu+1e-8 and
                    budget['wall_seconds'] <= wall_cap-wall+1e-8 and budget['disk_mib'] == disk_cap,
                    'Attempt reset/enlarged a cumulative phase allowance')
            expected = dict(CPUQuota='25%',Nice='19',CPUWeight='5',IOWeight='10',
                            MemoryHigh='448M',MemoryMax='512M',MemorySwapMax='0')
            require(all(raw['settings'].get(k) == v for k,v in expected.items()),'Wrong numerical resource profile')
            require(all(raw['environment'].get(k) == '1' for k in
                        ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS')),
                    'Numerical thread limits differ')
            require(type(actual['output_bytes']) is int and 0 <= actual['output_bytes'] <= disk_cap*1024**2,
                    'Output disk ceiling exceeded')
            disk += actual['output_bytes']
            cpu += values[0]
            wall += values[1]
        require(abs(cpu-pool['consumed_CPU_seconds']) < 1e-8 and
                abs(wall-pool['consumed_wall_seconds']) < 1e-8,'Cumulative resource totals differ')
        require(cpu <= cpu_cap+1e-8 and wall <= wall_cap+1e-8 and disk <= disk_cap*1024**2,
                'Cumulative resource ceiling exceeded')
    require(used_receipts == set(record['receipts']),'Uncharged or missing raw resource receipt')


def write_json(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    require(not path.exists(),'Preserve previous report outputs')
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def write_csv(path, rows):
    require(rows and not path.exists(),'Fresh nonempty CSV required')
    with path.open('x',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def authenticated_inputs(producer, delivery):
    complete = read(producer/'PRODUCER_COMPLETE.json')
    frozen = read(producer/'FROZEN.json')
    audit = read(delivery/'INDEPENDENT_VERIFIED.json')
    delivered = read(delivery/'DELIVERY_VERIFIED.json')
    require(audit['status'] == delivered['status'] == 'passed','Independent delivery/audit is unfinished')
    digest = sha(producer/'FROZEN.json')
    require(digest == complete['freeze_sha256'] == audit['producer_freeze_sha256'] ==
            delivered['producer_freeze_sha256'],'Producer/audit freeze identity differs')
    require(frozen['mode'] == 'main' and frozen['coefficient_scale'] == '1' and frozen['time'] == '3/5' and
            frozen['epsilons'] == EPSILONS and frozen['steps'] == list(STEPS),'Not the complete fixed main panel')
    require(complete['circuits'] == audit['circuits'] == delivered['circuits'] == 132 and
            complete['reference_evolutions'] == audit['reference_evolutions'] == 11,'Incomplete main audit')
    files = read(delivery/'DELIVERY_MANIFEST.json')
    actual = {p.relative_to(producer).as_posix() for p in producer.rglob('*') if p.is_file()}
    require(actual == set(files),'Original delivered file inventory changed')
    for name,row in files.items():
        path = producer/name
        require(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
                'Original delivered bytes changed: '+name)
    require(sha(delivery/'certificate.zip') == delivered['archive_sha256'],'Delivered archive changed')
    return frozen,complete,audit,delivered


def tables(producer, frozen, complete, audit):
    expected = [f'rows/epsilon-{at:02d}-r{r:02d}-{label}.json'
                for at in range(len(EPSILONS)) for r in STEPS for label in LABELS]
    require(complete['rows'] == expected and [x['path'] for x in audit['rows']] == expected,
            'Incomplete/duplicate/changed expected panel order')
    independent = {x['path']:x['independent_squared_infidelity'] for x in audit['rows']}
    errors,pairs,guards = [],[],[]
    for at,epsilon in enumerate(EPSILONS):
        for r in STEPS:
            intervals,independent_intervals = {},{}
            for label in LABELS:
                name = f'rows/epsilon-{at:02d}-r{r:02d}-{label}.json'
                row = read(producer/name)
                require((row['epsilon'],row['steps'],row['order']) == (epsilon,r,label),'Changed circuit key')
                low,high = interval(row['squared_infidelity'])
                require(0 <= low <= high <= 1,'True error interval outside[0,1]')
                check = independent[name]
                ilo,ihi = F(check['lo']),F(check['hi'])
                require(check['precision_bits'] == 160 and low <= ilo <= ihi <= high,
                        'Independent interval is not contained in producer enclosure')
                intervals[label] = low,high
                independent_intervals[label] = ilo,ihi
                errors.append(dict(epsilon=epsilon,physical_epsilon_one=epsilon=='1',time='3/5',steps=r,
                    order=label,squared_infidelity_lower=str(low),squared_infidelity_upper=str(high),
                    lower_approx=float(low),upper_approx=float(high),
                    independent_lower=str(ilo),independent_upper=str(ihi),
                    state_analytic_radius=row['state_analytic_radius'],
                    overlap_analytic_radius=row['analytic_overlap_radius'],
                    infidelity_analytic_radius=row['analytic_infidelity_radius']))
            decision = winner(*[intervals[x] for x in LABELS])
            idecision = winner(*[independent_intervals[x] for x in LABELS])
            original = next(x for x in complete['finite_rankings'] if x['epsilon']==epsilon and x['steps']==r)
            checked = next(x for x in audit['finite_rankings'] if x['epsilon']==epsilon and x['steps']==r)
            require(decision == original['strict_smaller_error'] and
                    idecision == checked['independent_strict_smaller_error'] and
                    decision in ('unresolved',idecision),'Changed finite ranking')
            pairs.append(dict(epsilon=epsilon,steps=r,strict_smaller_error=decision,
                              independent_strict_smaller_error=idecision,
                              prefix_minus_current_certified_gap_lower=str(intervals['diagonal_prefix'][0]-intervals['current'][1])))
        row = complete['response_guards'][at]
        require(row['epsilon'] == epsilon,'Guard panel identity changed')
        a,b = [interval(row['intervals'][x]) for x in LABELS]
        require(row['strict_smaller_response'] == winner(a,b),'Guard decision differs from its intervals')
        guards.append(dict(epsilon=epsilon,physical_epsilon_one=epsilon=='1',
            current_stddev_lower=str(a[0]),current_stddev_upper=str(a[1]),
            prefix_stddev_lower=str(b[0]),prefix_stddev_upper=str(b[1]),
            strict_smaller_response=row['strict_smaller_response']))
    return errors,pairs,guards


def dependency_record(path, model_sha, *, allow_development=False):
    if path is None:
        require(allow_development,'Actual merged input dependency required;development preparation needs an explicit flag')
        return dict(status='pinned_development_input_actual_issue17_merge_pending',
                    chemistry_producer='7c8c64acf0f13703df10b0335e0017da84e94f9b',model_sha256=model_sha)
    record = read(path)
    require(record.get('status') == 'verified_actual_merge' and record.get('repository') == 'bernalde/qhat',
            'Actual merged input provenance required')
    commit = record.get('merge_sha','')
    require(type(commit) is str and len(commit) == 40 and all(c in '0123456789abcdef' for c in commit),
            'Exact actual-merge SHA required')
    name = record.get('model_path','')
    relative = PurePosixPath(name)
    require(name.startswith('numerical_audit/transfer_holdouts/') and not relative.is_absolute() and
            '..' not in relative.parts and str(relative) == name,'Unsafe/wrong merged model path')
    blob = subprocess.check_output(['git','-C',str(ROOT),'show',commit+':'+name])
    require(hashlib.sha256(blob).hexdigest() == record.get('model_sha256') == model_sha,
            'Merged dependency model bytes differ from certified input')
    require(subprocess.run(['git','-C',str(ROOT),'merge-base','--is-ancestor',commit,'HEAD']).returncode == 0,
            'Verified numerical input merge is not in the current branch ancestry')
    return dict(record,offline_check='Exact Git model bytes and branch ancestry;live remote merge status established separately')


def trusted_verifier_sources(delivery):
    """Bind an executed audit revision to exact live and committed verifier bytes."""
    record = read(delivery/'TRUSTED_VERIFIER_SOURCES.json')
    delivered = read(delivery/'DELIVERY_VERIFIED.json')
    names = {'numerical_audit/molecular_guard/delivery.py',
             'numerical_audit/molecular_guard/independent.py'}
    require(record.get('head') == delivered.get('trusted_checker_head') and
            set(record.get('sha256',{})) == names, 'Verifier source binding is incomplete')
    commit = record['head']
    require(isinstance(commit,str) and len(commit) == 40 and
            all(c in '0123456789abcdef' for c in commit), 'Verifier source binding revision is invalid')
    for name,digest in record['sha256'].items():
        require(sha(ROOT/name) == digest, 'Verifier source binding differs from maintained source: '+name)
        blob = subprocess.check_output(['git','-C',str(ROOT),'show',commit+':'+name])
        require(hashlib.sha256(blob).hexdigest() == digest,
                'Verifier source binding differs from executed revision: '+name)
    return record['sha256']


def compare_reports(kept, replayed):
    for folder in (kept,replayed):
        manifest = read(folder/'REPORT_MANIFEST.json')
        require(set(manifest) == {p.name for p in folder.iterdir() if p.is_file() and p.name != 'REPORT_MANIFEST.json'},
                'Report manifest inventory differs')
        for name,digest in manifest.items():
            require(sha(folder/name) == digest,'Report bytes changed: '+name)
    for name in ('all_finite_intervals.csv','finite_comparisons.csv','response_guards.csv','README.md',
                 'EXECUTION_RESOURCES_BEFORE_REPORT.json'):
        require((kept/name).read_bytes() == (replayed/name).read_bytes(),'Report replay differs: '+name)
    a,b = read(kept/'SUMMARY.json'),read(replayed/'SUMMARY.json')
    # Execution revisions advance during publication/replay. Their independently
    # authenticated source hashes, all scientific data, and the input merge stay exact.
    names = {'numerical_audit/molecular_guard/delivery.py',
             'numerical_audit/molecular_guard/independent.py'}
    for row in (a,b):
        provenance = row['source_provenance']
        binding = provenance.get('verifier_sha256',{})
        require(set(binding) == names and all(isinstance(v,str) and len(v) == 64 and
                all(c in '0123456789abcdef' for c in v) for v in binding.values()),
                'Verifier source binding is incomplete')
        provenance.pop('reporter')
        provenance.pop('verifier')
    require(a == b,'Scientific/provenance report replay differs')


def build(producer, delivery, output, resources, dependency_path=None, *, allow_development=False):
    require(not output.exists(),'Use a fresh report directory')
    head = subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    require(subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+str(Path(__file__).relative_to(ROOT))])
            == Path(__file__).read_bytes(),'Commit reporting source before report generation')
    frozen,complete,audit,delivered = authenticated_inputs(producer,delivery)
    verifier_sources = trusted_verifier_sources(delivery)
    resource_record = read(resources)
    validate_resources(resource_record)
    require(resource_record['phases']['main']['attempts'] and resource_record['phases']['independent']['attempts'],
            'Completed main/independent resource receipts required')
    dependency = dependency_record(dependency_path,frozen['model_sha256'],allow_development=allow_development)
    errors,pairs,guards = tables(producer,frozen,complete,audit)
    leading = read(producer/'leading_response.json')
    domain = sufficient_epsilon_domain(interval(leading['current']['a']),interval(leading['diagonal_prefix']['a']),
                                       F(leading['current']['C']),F(leading['diagonal_prefix']['C']))
    finite_counts = dict(Counter(x['strict_smaller_error'] for x in pairs))
    physical_counts = dict(Counter(x['strict_smaller_error'] for x in pairs if x['epsilon']=='1'))
    guard_counts = dict(Counter(x['strict_smaller_response'] for x in guards))
    summary = dict(status='independent_software_intervals_verified',
        case='LiH_r1p713',time='3/5',state='HF',spin_modes=6,active_electrons=2,
        full_space_dimension=64,verified_invariant_sector_dimension=9,
        circuits=132,reference_evolutions=11,finite_comparisons=66,
        finite_counts=finite_counts,physical_epsilon_one_finite_counts=physical_counts,
        leading_guard_counts=guard_counts,physical_epsilon_one_leading_guard=guards[0]['strict_smaller_response'],
        sufficient_positive_epsilon_domain=domain,
        orders=frozen['orders'],response_constants={label:dict(a=leading[label]['a'],C=leading[label]['C']) for label in LABELS},
        source_provenance=dict(producer=frozen['source_commit'],verifier=delivered['trusted_checker_head'],
                               reporter=head,reporter_sha256=sha(Path(__file__)),verifier_sha256=verifier_sources,model_sha256=frozen['model_sha256']),
        input_dependency=dependency,
        formal_scope='Exact/rational software certificate;no concrete Lean instantiation or full numerical-search proof',
        limits=['One fixed active-space input,HF state,time and two predeclared first-order schedules.',
                'Physical leading-response guard and actual finite-step error rankings are separate conclusions.',
                'Squared infidelity values differ from1-abs(overlap),although their strict rankings agree.',
                'No SCF/integral/basis-set error bound,global optimality,conservation-only mechanism orhardware advantage.',
                'Source reviews,live merged-input checks and publication/merge status are recorded separately;software intervals do not establish these gates.'])
    output.mkdir(parents=True)
    write_csv(output/'all_finite_intervals.csv',errors)
    write_csv(output/'finite_comparisons.csv',pairs)
    write_csv(output/'response_guards.csv',guards)
    write_json(output/'SUMMARY.json',summary)
    write_json(output/'EXECUTION_RESOURCES_BEFORE_REPORT.json',resource_record)
    text=f'''# Fixed LiH finite-step interval certificate

All132circuits,11reference evolutions and66fixed schedule comparisons are retained.
The initial state isHF,time3/5,and the supplied rational Hamiltonian has six spin
modes/two active electrons. Complete-unit conservation was checked on the full
64-dimensional occupation space before using its invariant nine-dimensional sector.

Finite error comparisons(current/diagonal_prefix/unresolved):
{finite_counts.get('current',0)}/{finite_counts.get('diagonal_prefix',0)}/{finite_counts.get('unresolved',0)}.
At physicalepsilon1,the six finite-step comparisons give
{physical_counts.get('current',0)}/{physical_counts.get('diagonal_prefix',0)}/{physical_counts.get('unresolved',0)}.
These concern rigorous software intervals for squared infidelity1-|overlap|².
Their strict ordering agrees with1-|overlap|;the numerical error values differ.

The leading-response guard at physicalepsilon1 is **{guards[0]['strict_smaller_response']}**.
All eleven guard attempts remain in response_guards.csv,including vacuous cases.
The sufficient coupling-domain result in SUMMARY.json is an algebraic consequence
of the fixed slope bounds and constants;it does not replace any tested failure
or assert the physical Hamiltonian satisfies its leading-response guard.

A separate160-bit Fraction implementation reconstructed literal Pauli actions,
BCH corrections and generic matrix exponentials. Its rigorous state/error
intervals fit inside the producer's192-bit enclosures. One complete extracted
audit covers the byte-identical original payload. Archived Python was never
executed. Analytic tails and arithmetic rounding are separate;approximate states
were not normalized. Every interval and strict/unresolved pair is in the CSVs.

Resource boundaries are retained in EXECUTION_RESOURCES_BEFORE_REPORT.json;
completed reporting usage is recorded separately after this process exits.
Standalone propagation timers exclude reference/unit-operator preparation;
the complete phase receipts include these operations,imports,serialization and
retained failed/repeated preflight work. No timing is a hardware speedup claim.

Input dependency status:{dependency['status']}. Exact details are in SUMMARY.json.
Source review,publication and merge checks are separate gates. Generic Lean work is separate;no concrete molecular
witness or Python interval-arithmetic correctness proof has been checked in Lean.
The certificate concerns supplied exact rational coefficients,with no rigorous
SCF,integral or basis-set error estimate. It establishes no global order optimum,
conservation-only accuracy mechanism,larger-molecule validity or hardware result.
'''
    (output/'README.md').write_text(text)
    write_json(output/'REPORT_MANIFEST.json',{p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()})
    print(json.dumps({k:summary[k] for k in ('status','circuits','finite_counts','physical_epsilon_one_finite_counts',
                                           'physical_epsilon_one_leading_guard','sufficient_positive_epsilon_domain')}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer',type=Path,required=True)
    parser.add_argument('--delivery',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resources',type=Path,required=True)
    parser.add_argument('--input-dependency',type=Path)
    parser.add_argument('--allow-development-input',action='store_true')
    parser.add_argument('--compare-report',type=Path)
    args = parser.parse_args()
    build(args.producer,args.delivery,args.output,args.resources,args.input_dependency,
          allow_development=args.allow_development_input)
    if args.compare_report is not None:
        compare_reports(args.compare_report,args.output)
