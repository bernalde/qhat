"""Guarded metadata contracts; no molecular quality observations."""
from pathlib import Path
import hashlib
import json
import pytest

SOURCE = Path('/tmp/qhat-normal-window-2026-10-28/worktree')
STAGE = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/ordering_completion_2026-10-28/runs/integration_2026-10-28')
BASE = 'numerical_audit/prospective_window_ordering/'
FILES = [BASE+'test/test_panel.py', BASE+'test/test_preflight.py', BASE+'evidence/tools/test_analyze_study.py',
         BASE+'evidence/tools/test_independent_derived_check.py']
RUNTIME = [BASE+'panel.py', BASE+'preflight.py', BASE+'study_common.py',
           BASE+'evidence/tools/analyze_study.py', BASE+'evidence/tools/independent_derived_check.py']
hashes = {name: hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in FILES+RUNTIME}


class Reports:
    def __init__(self):
        self.rows = []

    def pytest_runtest_logreport(self, report):
        self.rows.append(dict(nodeid=report.nodeid, when=report.when,
                             outcome=report.outcome, duration_seconds=report.duration, report_type=type(report).__name__))


reports = Reports()
status = int(pytest.main([str(SOURCE/name) for name in FILES]+['-q', '-ra'], plugins=[reports]))
assert all(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() == expected
           for name, expected in hashes.items()), 'Metadata sources changed during tests'
record = dict(status='passed' if status == 0 else 'failed', exit_code=status,
              source_commit=None, base_commit='722013aa6d9a0a48776e26de5f356b5a441c3f61', source_state='uncommitted_tests_only_regression_delta_exact_hashes', source_hashes=hashes, driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              calls_passed=sum(r['when'] == 'call' and r['outcome'] == 'passed' for r in reports.rows),
              skipped=sum(r['outcome'] == 'skipped' for r in reports.rows), distinct_test_nodeids=len({r['nodeid'] for r in reports.rows if r['when']=='call'}), call_report_types={t:sum(r['when']=='call' and r['report_type']==t for r in reports.rows) for t in {r['report_type'] for r in reports.rows}}, phase_reports=reports.rows,
              scope='Literal source/receipt/metric contracts only; no molecular success graph or physical quality measured')
with (STAGE/'SOURCE_REVIEW_METADATA_CHECK.a4.red.json').open('x') as stream:
    json.dump(record, stream, indent=2, allow_nan=False)
    stream.write('\n')
raise SystemExit(status)
