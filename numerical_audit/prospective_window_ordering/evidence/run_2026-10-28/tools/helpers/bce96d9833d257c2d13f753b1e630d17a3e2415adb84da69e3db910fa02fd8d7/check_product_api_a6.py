"""Run focused contracts against the verified isolated installed product."""
from pathlib import Path
import hashlib
import importlib
import json
import os
import sys

ROOT = Path('/tmp/qhat-normal-window-2026-10-28/worktree')
STAGE = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/integration_2026-10-28')
build = json.loads((STAGE / 'BUILD_WHEEL.a6.json').read_text())
freeze_path = STAGE / 'PRODUCT_SOURCE_FREEZE.a9.json'
freeze = json.loads(freeze_path.read_text())
for name, expected in freeze['source_hashes'].items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
installed = Path(build['installed'])
assert hashlib.sha256(Path(build['wheel']).read_bytes()).hexdigest() == build['wheel_sha256']
sys.path.insert(0, str(installed))
loaded = {}
for name in ('sparse_response_ordering', 'strang_schedule', 'window_ordering', 'train_ordering'):
    module = importlib.import_module('qhat.analysis.' + name)
    path = Path(module.__file__).resolve()
    expected = 'analysis/' + name + '.py'
    assert path.is_relative_to(installed.resolve()), path
    assert hashlib.sha256(path.read_bytes()).hexdigest() == build['runtime_source_hashes'][expected]
    loaded[name] = str(path)
os.environ['QHAT_WINDOW_WHEEL'] = build['wheel']

import pytest


class Reports:
    def __init__(self):
        self.rows = []

    def pytest_runtest_logreport(self, report):
        self.rows.append(dict(nodeid=report.nodeid, when=report.when, outcome=report.outcome,
                              duration_seconds=report.duration))


reports = Reports()
files = ['analysis/tests/test_window_ordering.py']
selectors = 'installed_asset_manifest_binds_the_product_target or opt_in_roundtrip_keeps_labels_lengths_cost_and_native_authority or real_julia_normal_api_and_native_parity or installed_wheel_normal_cli_and_asset_discovery'
status = int(pytest.main([str(ROOT / path) for path in files] + ['-q','-ra','-k',selectors], plugins=[reports]))
record = dict(status='passed' if status == 0 else 'failed', exit_code=status,
              driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              product_source_freeze_sha256=hashlib.sha256(freeze_path.read_bytes()).hexdigest(),
              loaded_modules=loaded, wheel_sha256=build['wheel_sha256'],
              test_files=files, phase_reports=reports.rows,
              test_source_hashes={name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                  for name in files},
              calls_passed=sum(row['when'] == 'call' and row['outcome'] == 'passed' for row in reports.rows),
              skipped=sum(row['outcome'] == 'skipped' for row in reports.rows),
              scope='Parent test-process product imports pinned to isolated installed wheel; dedicated -S wheel CLI test independently exercises real installed Julia assets. No molecular quality measured.')
(STAGE / 'PRODUCT_API_CHECK.a6.json').write_text(json.dumps(record, indent=2) + '\n')
raise SystemExit(status)
