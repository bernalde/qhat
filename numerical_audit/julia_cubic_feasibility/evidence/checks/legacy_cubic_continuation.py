"""Invoke the unchanged retained optimizer from the same frozen parent.

No new Python numerical implementation. This is the same-budget legacy control
and independent scorer for the Julia window pilot.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import time

root=Path('/tmp/qhat-julia-cubic-2026-10-27/worktree')
common_path=root/'numerical_audit/scalable_ordering/common.py'
spec=importlib.util.spec_from_file_location('common',common_path)
common=importlib.util.module_from_spec(spec);sys.modules['common']=common;spec.loader.exec_module(common)
ordering=common.bootstrap()
assert Path(ordering.__file__).resolve()==root/'analysis/sparse_response_ordering.py'
case=sys.argv[1];out=Path(sys.argv[2]);assert not out.exists()
base=Path('/tmp/qhat-scalable-ordering-2026-10-27')
model_path=base/'inputs-a2/inputs'/case/'model.json'
parent_path=base/'study-a2/constructors'/f'{case}.cubic.json'
began=time.process_time();wall=time.monotonic()
model=common.read(model_path);parent=common.read(parent_path)
assert hashlib.sha256(model_path.read_bytes()).hexdigest()==parent['model_sha256']
units=common.blocks(model);budget=30*len(model['coefficients']);hf=(1<<model['electrons'])-1
allowance=32.0-(time.process_time()-began)
assert allowance>0
trained=ordering.train_ordering(units,model['n'],hf,time=1.,rotation_budget=budget,
    cpu_seconds=allowance,max_evaluations=1_000_000,max_moves=5,starts=[parent['order']])
construction_cpu=time.process_time()-began
check_start=time.process_time()
context=ordering.SparseResponse(units,model['n'],hf,time=1.)
verified=context.score(trained.order,budget)
parent_score=context.score(parent['order'],budget)
assert sorted(trained.order)==list(range(len(units)))
record=dict(case=case,status=trained.status,order=[i+1 for i in trained.order],
    score='unavailable' if trained.score is None else trained.score.value,
    independently_checked_score=verified.value,parent_score=parent_score.value,
    constructor_CPU_seconds=construction_cpu,validation_CPU_seconds=time.process_time()-check_start,
    wall_seconds=time.monotonic()-wall,evaluations=trained.evaluations,accepted_moves=trained.accepted_moves,
    additional_CPU_allowance_seconds=32.0,starts=1,stats_json=json.dumps(trained.stats,sort_keys=True),
    model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
    parent_sha256=hashlib.sha256(parent_path.read_bytes()).hexdigest(),
    implementation_sha256=hashlib.sha256(Path(ordering.__file__).read_bytes()).hexdigest(),
    scope='Unchanged legacy whole-order insertion continuation from frozen cubic parent; Julia indices in returned order.')
out.write_text('\n'.join(k+' = '+json.dumps(v,allow_nan=False) for k,v in record.items())+'\n')
print(json.dumps({k:record[k] for k in ('case','status','constructor_CPU_seconds','evaluations','independently_checked_score')}))
