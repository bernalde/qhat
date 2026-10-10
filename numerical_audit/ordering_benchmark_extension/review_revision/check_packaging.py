import hashlib,json,os,subprocess,tempfile,zipfile
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
Q=Path('/home/bernalde/repos/fermionic-ordering-research/qhat');P=Q.parent/'paper';D=Q/'numerical_audit/ordering_benchmark_extension/review_revision';py=Q.parent/'benchmark-venv-py312/bin/python';out=D/'packaging_check.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(v,m):
 if not v:raise RuntimeError(m)
record=json.loads((D/'archive_record.json').read_text());require(sha(D/'review_supplement.zip')==record['archive_sha256'],'Archive identity')
with tempfile.TemporaryDirectory(prefix='review-supplement-extract-') as t:
 x=Path(t)
 with zipfile.ZipFile(D/'review_supplement.zip') as z:z.extractall(x)
 cmds=[([str(py),'-O',str(x/'archive_support/check_historical_states.py'),'--archive-root',str(Q.parent/'extension-clean-replay'),'--output',str(x/'fresh_historical.json')],None),([str(py),str(x/'strang_controls/render_table.py')],x),([str(py),str(x/'analysis/render_formula_table.py'),'--checks',str(x/'analysis/framing_checks.json'),'--output',str(x/'formula_order_table.tex')],x)]
 logs=[]
 for cmd,cwd in cmds:
  r=subprocess.run(cmd,cwd=cwd or '/tmp',env=os.environ.copy(),text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90);require(r.returncode==0,r.stdout);logs.append({'command':cmd,'returncode':r.returncode,'output':r.stdout[-1200:]})
 require((x/'strang_controls/controls_table.tex').read_bytes()==(P/'theory/review_revision_evidence/strang_controls/controls_table.tex').read_bytes(),'Control table packaging')
 require((x/'formula_order_table.tex').read_bytes()==(P/'theory/review_revision_evidence/formula_order_table.tex').read_bytes(),'Formula table packaging')
 result={'status':'passed','archive_sha256':record['archive_sha256'],'historical_check':json.loads((x/'fresh_historical.json').read_text()),'both_manuscript_tables_regenerated_byte_identical':True,'commands':logs,'scope':'Fresh supplement extraction, all30savedhistorical arrays and bothnewmanuscript tables; independent molecular propagation is the separately recorded504-state check.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS: fresh supplement extraction, all30 historical arrays, two byte-identical manuscript tables')
