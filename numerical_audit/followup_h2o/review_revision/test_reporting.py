import resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
import copy,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
q=Path(__file__).resolve().parents[3];root=q/'numerical_audit/followup_h2o';rev=root/'review_revision'
spec=importlib.util.spec_from_file_location('revised',rev/'summarize_h2o.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s=json.loads((rev/'presentation/summary.json').read_text())
old=json.loads((root/'presentation/summary.json').read_text())
compared=['policy_summary','comparisons','convergence_anchors','low_energy_outcomes','validation','max_leakage_by_order']
for key in compared:
 if s[key]!=old[key]:raise RuntimeError('Scientific summary changed: '+key)
base=m.interpretation(s['policy_summary'],s['coefficient_policy_sensitivity'])
ps=copy.deepcopy(s['policy_summary']);ps['historical_policy']['best_method_sets']['fermionic_signed']=3;ps['historical_policy']['median_ratio']['jw_signed']=9.25
changed=m.interpretation(ps,s['coefficient_policy_sensitivity'])
if base==changed or '3 of 15' not in changed or '9.25' not in changed:raise RuntimeError('Outcome facts hard-coded')
pd=copy.deepcopy(s['coefficient_policy_sensitivity']);pd[0]['same_retained_H']=False;pd[0]['same_outcomes']['jw_signed']=False
changed=m.interpretation(s['policy_summary'],pd)
if '13 of 15' not in changed or '14 of 15' not in changed:raise RuntimeError('Policy facts hard-coded')
designs=[json.loads(p.read_text()) for p in (root/'low_energy').glob('H2O*/design.json')]
designs[0]['states']['sector_next']['S2_expectation']=0.5
if m.spin_summary(designs)['sector_next']['spin_label']!='not assigned':raise RuntimeError('Spin label hard-coded')
for relative in ['summarize_h2o.py','check_guards.py','README.md','artifact_manifest.json']:
 f=root/relative;blob=subprocess.check_output(['git','-C',str(q),'show','HEAD:'+str(f.relative_to(q))])
 if f.read_bytes()!=blob:raise RuntimeError('Historical source modified: '+relative)
# CLI rejects missing explicit output before any imports/runs beyond stdlib in guard.
a=subprocess.run([sys.executable,str(rev/'check_guards.py'),'--qhat-root',str(q)],capture_output=True,text=True)
if a.returncode==0 or '--output' not in a.stderr:raise RuntimeError('Output not required')
print(json.dumps({'status':'passed','mutations_recomputed':['unique_best_count_and_median','policy_identity_and_outcome_counts','spin_label'],'historical_files_byte_identical_to_HEAD':4,'guard_missing_output_rejected':True,'unchanged_scientific_summary_fields':compared},indent=2))
