"""Derive scoped bound ratios and cumulative costs from retained raw evidence."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(value):
    pieces = str(value).split('//')
    if len(pieces) > 2 or any(len(x) > 2500 for x in pieces):
        raise ValueError('Rational spelling bound')
    return Fraction(int(pieces[0]), int(pieces[1]) if len(pieces) == 2 else 1)


def metric(value):
    return dict(exact=str(value), decimal=float(value))


def main(root, output):
    if output.exists():
        raise ValueError('Fresh summary required')
    decision = json.loads((root / 'PREFLIGHT_DECISION.a1.json').read_text())
    if decision['status'] != 'failed_feasibility_gate' or decision['quality_used_in_decision']:
        raise ValueError('Original resource-only preflight disposition differs')
    rows = []
    for short in ('CH4', 'CO'):
        bound_path = root / 'windows-a1' / (short + '.BOUND.a1.toml')
        check_path = root / 'windows-a1' / (short + '.ORACLE_CHECK.a1.toml')
        bound = tomllib.loads(bound_path.read_text())
        checked = tomllib.loads(check_path.read_text())
        if not bound['complete'] or not bound['independent_exact_arithmetic_verified'] or checked['bound_sha256'] != sha(bound_path):
            raise ValueError('Incomplete or stale fixed-window check')
        if bound['flat_sha256'] != checked['flat_sha256'] or checked['full_molecule_ratio_claimed']:
            raise ValueError('Window target or scope differs')
        lower = rational(bound['L_Q'])
        optimum = rational(checked['optimum_Q'])
        incumbent = rational(bound['Q_inc'])
        if not 0 <= lower <= optimum <= incumbent:
            raise ValueError('Exact lower/optimum/upper ordering differs')
        r = int(bound['r_inc'])
        if rational(bound['L_J']) != lower / r**4 or rational(bound['UB_J']) != incumbent / r**4:
            raise ValueError('Fixed-step J normalization differs')
        rows.append(dict(case=short, domain='fixed eight-label interior window; all outside labels and endpoints frozen',
            steps=r, positive_lower_bound=lower > 0,
            retained_parent_certified_ratio=metric(incumbent / lower) if lower else None,
            archived_optimum_certified_ratio=metric(optimum / lower) if lower else None,
            lower_fraction_of_archived_optimum=metric(lower / optimum) if optimum else None,
            retained_parent_over_archived_optimum=metric(incumbent / optimum) if optimum else None,
            bound_sha256=sha(bound_path), oracle_check_sha256=sha(check_path),
            oracle_sha256=checked['oracle_sha256'], flat_sha256=bound['flat_sha256'],
            new_enumeration=False, full_domain_ratio_claimed=False))
    attempts = []
    factors = {'h': 3600, 'min': 60, 's': 1, 'ms': .001, 'us': .000001, 'µs': .000001, 'ns': .000000001}
    for path in sorted((root / 'logs').glob('*.resources.json')):
        receipt = json.loads(path.read_text())
        label = path.name.removesuffix('.resources.json')
        cpu = receipt['final_accounting']['CPU_seconds']
        if cpu is None:
            cpu = receipt['budget']['cpu_seconds']
        final = root / ('FINAL_ACCOUNTING.' + label.removeprefix('molecular-bounds-') + '.json')
        if final.is_file():
            cpu = max(cpu, json.loads(final.read_text())['charged_CPU_seconds'])
        log = root / 'logs' / (label + '.log')
        for text in re.findall(r'^CPU time consumed: (.+)$', log.read_text(), re.MULTILINE):
            tokens = re.findall(r'([0-9.]+)(h|min|ms|us|µs|ns|s)', text)
            if tokens:
                cpu = max(cpu, sum(float(value) * factors[unit] for value, unit in tokens) + .001)
        attempts.append(dict(label=label, exit_code=receipt['exit_code'],
            charged_CPU_seconds=cpu, guarded_wall_seconds=receipt['elapsed'], receipt_sha256=sha(path),
            log_sha256=sha(log), scope='All failed and completed attempts retained; rounded service summary charged conservatively'))
    summary = dict(schema='molecular-bounds-analysis-v1', analysis_sha256=sha(Path(__file__)),
        full_domain_status='No full-domain lower bound: complete projection refused at50k scopes in altered12-mode preflight',
        resource_decision_sha256=sha(root / 'PREFLIGHT_DECISION.a1.json'),
        main_full_cases_performed=0, main_full_cases_unperformed=6,
        refusal_check_scope='Metadata/caps authenticated only; scientific_verification_performed=false',
        restricted_windows=rows, attempts=attempts,
        cumulative_CPU_seconds=sum(row['charged_CPU_seconds'] for row in attempts),
        cumulative_guarded_wall_seconds=sum(row['guarded_wall_seconds'] for row in attempts),
        analysis_attempt_cost_included=False,
        limitations=['No global polynomial-time approximation ratio', 'No physical-error guarantee',
                     'No exact-transcendental kernel claim', 'No cold/warm speedup comparison'])
    with output.open('x') as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(windows=[dict(case=row['case'], retained_parent_ratio=row['retained_parent_certified_ratio']['decimal'],
        archived_optimum_ratio=row['archived_optimum_certified_ratio']['decimal']) for row in rows],
        prior_CPU_seconds=summary['cumulative_CPU_seconds'], prior_wall_seconds=summary['cumulative_guarded_wall_seconds'])), flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: analyze.py RUN_ROOT NEW_SUMMARY_JSON')
    main(Path(sys.argv[1]), Path(sys.argv[2]))
