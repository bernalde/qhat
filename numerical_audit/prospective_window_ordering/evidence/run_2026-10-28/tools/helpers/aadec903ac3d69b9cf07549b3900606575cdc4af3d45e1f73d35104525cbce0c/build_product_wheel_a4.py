"""Build/install the complete declared package payload in a bounded scratch tree."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile

SOURCE = Path('/tmp/qhat-normal-window-2026-10-28/worktree')
STAGE = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/integration_2026-10-28')
BUILD = STAGE / 'build-source-a4'
WHEELS = STAGE / 'wheels-a4'
INSTALLED = STAGE / 'installed-a4'
FREEZE = json.loads((STAGE / 'PRODUCT_SOURCE_FREEZE.a7.json').read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert not any(path.exists() for path in (BUILD, WHEELS, INSTALLED))
    tracked = subprocess.check_output(['git', '-C', str(SOURCE), 'ls-files', '-z']).decode().split('\0')
    names = {name for name in tracked if name and
             (len(Path(name).parts) == 1 or Path(name).parts[0] in ('analysis', 'common', 'hamiltonian_generator'))}
    names.update(name for name in FREEZE['source_hashes'] if name.startswith('analysis/'))
    names.add('pyproject.toml')
    total = sum((SOURCE / name).stat().st_size for name in names)
    assert total <= 16 * 1024**2, 'Declared complete build payload exceeds scratch reservation'
    for name, expected in FREEZE['source_hashes'].items():
        assert digest(SOURCE / name) == expected, name
    payload = {}
    for name in sorted(names):
        original = SOURCE / name
        assert original.is_file() and not original.is_symlink(), name
        destination = BUILD / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, destination)
        assert destination.read_bytes() == original.read_bytes(), name
        payload[name] = digest(original)
    (STAGE / 'BUILD_INPUT_MANIFEST.a4.json').write_text(json.dumps(dict(
        scope='Complete current declared qhat/analysis/common/hamiltonian_generator build inputs; numerical audit archives are not installed packages',
        bytes=total, files=payload), indent=2) + '\n')
    env = os.environ.copy()
    env['PIP_CACHE_DIR'] = str(STAGE / 'pip-cache')
    subprocess.run([sys.executable, '-m', 'pip', 'wheel', '.', '--no-deps', '--no-build-isolation',
                    '--wheel-dir', str(WHEELS)], cwd=BUILD, env=env, check=True)
    wheels = list(WHEELS.glob('*.whl'))
    assert len(wheels) == 1
    wheel = wheels[0]
    runtime = [name for name in FREEZE['source_hashes'] if
               name.startswith('analysis/julia_window/') or
               name in ('analysis/window_ordering.py', 'analysis/sparse_response_ordering.py',
                        'analysis/strang_schedule.py', 'analysis/train_ordering.py')]
    with zipfile.ZipFile(wheel) as archive:
        for name in runtime:
            assert archive.read('qhat/' + name) == (SOURCE / name).read_bytes(), name
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--no-deps', '--no-index',
                    '--target', str(INSTALLED), str(wheel)], env=env, check=True)
    for name in runtime:
        assert (INSTALLED / 'qhat' / name).read_bytes() == (SOURCE / name).read_bytes(), name
    record = dict(status='passed', wheel=str(wheel), wheel_sha256=digest(wheel),
                  installed=str(INSTALLED), complete_build_payload_bytes=total,
                  package_payload_files=len(payload), runtime_files_verified=runtime,
                  runtime_source_hashes={name: digest(SOURCE / name) for name in runtime},
                  source_freeze='PRODUCT_SOURCE_FREEZE.a7.json',
                  environment_modified=False)
    (STAGE / 'BUILD_WHEEL.a4.json').write_text(json.dumps(record, indent=2) + '\n')
    print('WHEEL_AND_ISOLATED_INSTALL_PASSED', wheel.name, flush=True)


if __name__ == '__main__':
    main()
