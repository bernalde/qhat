"""Read-only original-analysis representation diagnosis; no numerical replay."""
import json,sys
from pathlib import Path
r=Path('/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30/runs/preflight-a4')
data=r/'relocated-replay-a1/extracted/data'
sys.path.insert(0,str(data/'source_snapshot/numerical_audit/molecular_generalization'))
import analyze
actual=analyze.analyze(data,preflight=True)
expected=json.loads((r/'data/ANALYSIS.json').read_text())
diffs=[]
def compare(a,b,path=''):
    if type(a)is not type(b):
        diffs.append(dict(path=path,actual_type=type(a).__name__,expected_type=type(b).__name__,serialized_equal=json.loads(json.dumps(a,allow_nan=False))==b));return
    if isinstance(a,dict):
        if set(a)!=set(b):diffs.append(dict(path=path,keys_actual=sorted(a),keys_expected=sorted(b)))
        for key in sorted(set(a)&set(b)):compare(a[key],b[key],path+'/'+str(key))
    elif isinstance(a,(list,tuple)):
        if len(a)!=len(b):diffs.append(dict(path=path,length_actual=len(a),length_expected=len(b)))
        for i,(x,y)in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif a!=b:diffs.append(dict(path=path,actual=a,expected=b))
compare(actual,expected)
report=dict(schema='molecular-analysis-representation-diagnosis-v1',analysis_complete=actual['complete'],in_memory_equal=actual==expected,serialized_JSON_equal=json.loads(json.dumps(actual,allow_nan=False))==expected,differences=diffs,no_new_numerical_calculation=True)
output=r/'relocated-replay-a1/metadata-diagnosis-a1/DIAGNOSIS.json'
with output.open('x')as stream:json.dump(report,stream,indent=2,allow_nan=False)
print(json.dumps(report),flush=True)
