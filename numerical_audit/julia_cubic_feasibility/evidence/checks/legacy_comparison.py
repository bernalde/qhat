"""Run the existing archived Python solver unchanged as an independent oracle.

The new implementation lives in Julia. This adapter only invokes the retained
baseline and preserves its measurements and exact certificate comparisons.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import statistics
import sys
import time

worktree=Path('/tmp/qhat-julia-cubic-2026-10-27/worktree')
source=worktree/'numerical_audit/exact_response_windows/exact.py'
spec=importlib.util.spec_from_file_location('retained_exact_window_baseline',source)
baseline=importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)
assert Path(baseline.__file__).resolve()==source.resolve()
root=source.parent/'evidence/results'
out=Path(sys.argv[1]);assert not out.exists()
records=[]
for case in ('F2_r1p55','H2O_scale1p10','Li2_r2p95'):
    domain=json.loads((root/'domains'/f'{case}.json').read_text())
    expected=json.loads((root/'certificates'/f'{case}.json').read_text())
    samples=[]
    for repetition in range(3):
        start=time.process_time()
        result=baseline.solve(domain,cpu_limit=10)
        total=time.process_time()-start
        for key in ('incumbent_integer','selected_window','minimizers','evaluated'):
            assert result[key]==expected[key],(case,key)
        assert result['complete']
        item=dict(repetition=repetition+1,total_CPU_seconds=total,solver=result)
        samples.append(item)
        checkpoint=out.with_name(out.name+f'.{case}.r{repetition+1}.json')
        assert not checkpoint.exists();checkpoint.write_text(json.dumps(item,indent=2)+'\n')
    records.append(dict(case=case,samples=samples,median_total_CPU_seconds=statistics.median(v['total_CPU_seconds'] for v in samples)))
out.write_text(json.dumps(dict(status='passed',source_path=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cases=records,scope='Fresh executions of the unchanged retained Python enumerator; original archived first-order JD exact integer objective.'),indent=2)+'\n')
print(json.dumps(dict(status='passed',cases=[dict(case=r['case'],median_CPU_seconds=r['median_total_CPU_seconds']) for r in records])))
