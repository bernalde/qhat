"""Old-only two-version Julia bridge followed by unchanged native checks."""
import hashlib,json,math,subprocess,sys,tomllib
from pathlib import Path
p=Path(sys.argv[1]);plan=json.loads(p.read_text())
for key in ('old_kernel','new_kernel','request','model','julia_operator','python_operator'):
 assert hashlib.sha256(Path(plan[key]).read_bytes()).hexdigest()==plan[key+'_sha256'],key
command=[plan['julia'],'--startup-file=no','--history-file=no','--threads=1','--gcthreads=1','--heap-size-hint=32M','--optimize=2','--project='+plan['project'],plan['julia_operator'],plan['old_kernel'],plan['new_kernel'],plan['request'],plan['fresh_output']]
status=subprocess.run(command,check=False).returncode
work=Path(plan['fresh_output']);report=dict(schema='molecular-old-native-cache-bridge-v1',complete=False,Julia_return_code=status,plan_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),native_rows=[])
if (work/'BRIDGE.toml').is_file():
 report['Julia_report_sha256']=hashlib.sha256((work/'BRIDGE.toml').read_bytes()).hexdigest()
 julia=tomllib.loads((work/'BRIDGE.toml').read_text())
 sys.path.insert(0,str(Path(plan['old_snapshot_root'])/'numerical_audit/molecular_generalization'))
 import common as k
 model=k.read(Path(plan['model']));units=k.blocks(model);native=k.bootstrap()
 assert Path(native.__file__).resolve()==Path(plan['old_snapshot_root'])/'analysis/sparse_response_ordering.py'
 context=native.SparseResponse(units,model['n'],(1<<model['electrons'])-1,time=1.,max_entries=50000,max_cached_actions=100000)
 close=lambda a,b:math.isfinite(a)and math.isfinite(b)and abs(a-b)<=max(1e-28,abs(a)*1e-8)
 for row in julia['rows']:
  order=k.validate_order(row['order_zero_based'],len(units));score=context.score(order,30*len(model['coefficients']))
  actual=dict(value=score.value,cubic_response=score.cubic_response,steps=score.steps)
  good=all(actual['steps']==row[which]['steps']and close(actual['value'],row[which]['value'])and close(actual['cubic_response'],row[which]['cubic_response'])for which in ('old_score','new_score'))
  report['native_rows'].append(dict(label=row['label'],native_score=actual,verified=good))
 report['complete']=status==0 and julia['complete']is True and len(report['native_rows'])==2 and all(row['verified']for row in report['native_rows'])
for key in ('old_kernel','new_kernel','request','model'):
 assert hashlib.sha256(Path(plan[key]).read_bytes()).hexdigest()==plan[key+'_sha256'],'changed '+key
with (work/'NATIVE_BRIDGE.json').open('x')as out:json.dump(report,out,indent=2,allow_nan=False);out.write('\n')
print('OLD_NATIVE_CACHE_BRIDGE_COMPLETE',report['complete'],flush=True)
sys.exit(0 if report['complete']else 2)
