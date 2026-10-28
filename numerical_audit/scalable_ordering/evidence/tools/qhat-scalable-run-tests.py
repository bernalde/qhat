import hashlib,json,os,subprocess,sys
from pathlib import Path
n=Path('/tmp/qhat-scalable-ordering-2026-10-27');q=n/'worktree'
label=sys.argv[1] if len(sys.argv)>1 else 'unlabelled'
files=subprocess.check_output(['git','status','--porcelain','-uall'],cwd=q,text=True).splitlines()
paths=[line[3:] for line in files]
before={p:hashlib.sha256((q/p).read_bytes()).hexdigest() for p in paths}
(n/('TESTED_SOURCE_INPUTS.'+label+'.json')).write_text(json.dumps(before,indent=2)+'\n')
import qhat.analysis.sparse_response_ordering as module
assert Path(module.__file__).resolve()==q/'analysis/sparse_response_ordering.py'
cohort=['analysis/tests/test_response_resolution.py','analysis/tests/test_sparse_response_ordering.py','analysis/tests/test_sparse_ordering_updates.py',
        'analysis/tests/test_train_ordering_cli.py','numerical_audit/scalable_ordering/test_pipeline.py']
print('EXACT_SOURCE',module.__file__,flush=True)
result=subprocess.run([sys.executable,'-m','pytest',*cohort,'-q'],cwd=q)
after={p:hashlib.sha256((q/p).read_bytes()).hexdigest() for p in paths};assert before==after
(n/('TEST_RESULT.'+label+'.json')).write_text(json.dumps({'exit_code':result.returncode,'cohort':cohort,'source_inputs':str(n/('TESTED_SOURCE_INPUTS.'+label+'.json')),'source_stable':True,'interpreter':sys.executable},indent=2)+'\n')
raise SystemExit(result.returncode)
