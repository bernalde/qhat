"""Authenticate retained models and prepare exact altered/development registries."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import tomllib

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / 'julia_cubic_feasibility/evidence/inputs/cubic'
EXPECTED = {
    'CH4_r1p09': '00a283f3736e331b2c3bfc140213e1d5d1d3a62f119b34ab1450a1b091afbc97',
    'CO_r1p13': 'b78921326c04f331a0d95b8fadf1be600bb6bfd0caa65e464d24c7aa6de5fe52',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, text):
    with path.open('x') as stream:
        stream.write(text)


def registry(path, case, input_path, *, preflight=False):
    cpu, wall = (170.0, 850.0) if preflight else (180.0, 900.0)
    write(path, '\n'.join([
        'schema = "molecular-bounds-registry-v1"', '[caps]', 'profile = "default"',
        '[options]', 'iterations = 64', 'max_witnesses = 8', 'target_count = 8',
        'comparison_count = 1', f'attempt_CPU_seconds = {cpu}', f'attempt_wall_seconds = {wall}',
        '[[cases]]', 'id = ' + json.dumps(case), 'input = ' + json.dumps(str(input_path.resolve())),
        'freeze_kernel = true', 'time = 1.0', 'order_field = "cubic_order"',
        'budget_multiplier = 30', '']))


def prepare(output):
    if output.exists():
        raise ValueError('Fresh output directory required; preserve earlier attempts')
    models = {}
    for case, expected in EXPECTED.items():
        path = ARCHIVE / (case + '.toml')
        if not path.is_file() or path.stat().st_size > 8 * 1024**2 or sha(path) != expected:
            raise ValueError('Archived source identity differs: ' + case)
        model = tomllib.loads(path.read_text())
        m = len(model['units'])
        if model['n'] != 12 or model['electrons'] != 6 or sorted(model['cubic_order']) != list(range(1, m + 1)):
            raise ValueError('Original model register/order differs')
        models[case] = model
    selected = min(models, key=lambda case: (-sum(len(u['words']) for u in models[case]['units']), case))
    parent = models[selected]
    output.mkdir(parents=True)
    altered = output / (selected + '.coeff7over5.toml')
    denominator = 5 * int(parent['denominator'])
    lines = ['case = ' + json.dumps(selected + '_preflight_7over5'), 'n = 12', 'electrons = 6',
             'denominator = ' + json.dumps(str(denominator)),
             'cubic_order = ' + json.dumps(parent['cubic_order']),
             'parent_input_sha256 = ' + json.dumps(EXPECTED[selected]), 'coefficient_scale = "7//5"']
    for unit in parent['units']:
        lines += ['', '[[units]]', 'words = ' + json.dumps(unit['words']),
                  'numerators = ' + json.dumps([str(7 * int(x)) for x in unit['numerators']])]
    write(altered, '\n'.join(lines) + '\n')
    restored = tomllib.loads(altered.read_text())
    for old, new in zip(parent['units'], restored['units'], strict=True):
        if new['words'] != old['words'] or any(
            Fraction(int(b), denominator) != Fraction(7, 5) * Fraction(int(a), int(parent['denominator']))
            for a, b in zip(old['numerators'], new['numerators'], strict=True)):
            raise ValueError('Exact coefficient scaling or original word inventory differs')
    registry(output / 'PREFLIGHT.registry.toml', 'CH4-altered-preflight', altered, preflight=True)
    for case in models:
        registry(output / (case + '.registry.toml'), case + '-development', ARCHIVE / (case + '.toml'))
    manifest = dict(selected=selected, scale='7//5', source_inputs=EXPECTED,
        altered_input_sha256=sha(altered), altered_input=altered.name,
        registry_sha256={p.name: sha(p) for p in output.glob('*.registry.toml')},
        coefficient_exactness_verified=True, no_ordering_or_physical_error_computed=True,
        full_incumbent_policy='Original archived cubic_order; historical construction costs remain historical',
        prospective_incumbent_policy='Conserving full-cubic arm repetition1, frozen before physical error; missing artifacts censor rather than substitute')
    write(output / 'PREPARATION.json', json.dumps(manifest, indent=2) + '\n')
    print('DEVELOPMENT_REGISTRIES_FROZEN', selected, flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: prepare_development.py NEW_OUTPUT_DIRECTORY')
    prepare(Path(sys.argv[1]))
