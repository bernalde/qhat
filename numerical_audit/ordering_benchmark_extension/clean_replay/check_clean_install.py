from pathlib import Path
import hashlib,json,os,subprocess,sys,time,zipfile
import numpy as np

root=Path('/home/bernalde/repos/fermionic-ordering-research')
archive=root/'artifact-staging/extension_validation.zip'
out=root/'extension-clean-replay'
out.mkdir()
with zipfile.ZipFile(archive) as z:z.extractall(out)
fresh=out/'fresh';fresh.mkdir()
ext=out/'extension';qhat=out/'qhat'
env=os.environ.copy()
env.update({'QHAT_ROOT':str(qhat),'PYTHONPATH':str(out),'NUMBA_CACHE_DIR':str(fresh/'numba-cache'),'MPLCONFIGDIR':str(fresh/'mpl-cache'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
commands=[
 ['-O','extension/independent/check_references.py','--qhat-root','qhat','--output','fresh/references'],
 ['-O','extension/independent/check_partitions.py','--qhat-root','qhat','--output','fresh/partitions.json'],
 ['-O','extension/independent/objective_check.py','--campaign','extension/optimized','--output','fresh/objectives'],
]
replays=[('fixed','CH2O_s1.00','jw_magnitude',2),('optimized','CH2O_s1.00','strang_class_refined_response',2),('optimized','nh3_coeff1p00','class_refined_response',10),('published/v2','f2_coeff1p00','tranter_v2_l1_first_order',10),('published/v2','nh3_coeff1p00','tranter_v2_l1_strang',2)]
for i,(campaign,case,order,budget) in enumerate(replays):
 commands.append(['-O','extension/packaging/replay_one.py','--extension','extension','--qhat-root','qhat','--campaign',campaign,'--case',case,'--ordering',order,'--budget',str(budget),'--output',f'fresh/production_{i}.json'])
commands.extend([['extension/analysis/summarize.py','--root','extension','--output','fresh/analysis'],['extension/analysis/timings.py','--root','extension','--output','fresh/analysis'],['extension/analysis/render.py','--analysis','fresh/analysis','--output','fresh/figures']])
steps=[]
for i,cmd in enumerate(commands):
 started=time.monotonic();r=subprocess.run([sys.executable,*cmd],cwd=out,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (fresh/f'command_{i}.log').write_text(r.stdout)
 steps.append({'arguments':cmd,'returncode':r.returncode,'seconds':time.monotonic()-started})
 print(i,r.returncode,r.stdout[-800:],flush=True)
 if r.returncode:raise RuntimeError(f'Clean replay command {i} failed')
reference_differences=[]
for p in sorted((fresh/'references').glob('*.npz')):
 with np.load(p) as a,np.load(ext/'independent/references'/p.name) as b:
  differences={k:float(np.linalg.norm(a[k]-b[k])) for k in ('initial','exact')}
 if max(differences.values())>2e-11:raise RuntimeError('Archived independent reference cache differs: '+p.name)
 reference_differences.append({'case':p.stem,'differences':differences})
for name in ('summary.json','timings.json','optimizer_runtime.csv','all_rows.csv'):
 if (fresh/'analysis'/name).read_bytes()!=(ext/'analysis'/name).read_bytes():raise RuntimeError('Analysis reconstruction differs: '+name)
for name in ('headline_table.tex','all_comparisons.tex'):
 if (fresh/'figures'/name).read_bytes()!=(ext/'analysis/figures'/name).read_bytes():raise RuntimeError('Table reconstruction differs: '+name)
records=[json.loads((fresh/f'production_{i}.json').read_text()) for i in range(len(replays))]
result={'status':'passed','archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'python':sys.version,'python_executable':sys.executable,'clean_environment':json.loads((ext/'packaging/clean_environment.json').read_text()),'commands':steps,'fresh_reference_cases':len(reference_differences),'reference_cache_comparisons':reference_differences,'fresh_production_replays':records,'maximum_production_state_difference':max(x['state_distance'] for x in records),'analysis_and_tables_byte_identical':True,'scope':'Fresh exact pinned installation, archived-code import containment, 14 independent target references and partitions, 56 independent objectives, five fresh production circuits, and full non-propagating analysis/table regeneration. This is not full regeneration of every molecular evolution.'}
(fresh/'clean_replay.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','fresh_reference_cases','maximum_production_state_difference','analysis_and_tables_byte_identical')}),flush=True)
