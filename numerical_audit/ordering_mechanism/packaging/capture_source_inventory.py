#!/usr/bin/env python3
"""Post-run source inventory, explicitly not a retroactive run-time attestation."""
import hashlib,importlib.metadata,json,subprocess,sys
from pathlib import Path
from qhat.analysis import benchmark_b2_active_spaces_matrix_free
from qhat.numerical_audit import reference_tools,run_frozen_audit

root=Path('/home/bernalde/repos/qhat').resolve()
files={}
for name,mod in sorted(sys.modules.items()):
 if name=='qhat' or name.startswith('qhat.'):
  file=getattr(mod,'__file__',None)
  if not file:continue
  p=Path(file).resolve()
  if not p.is_relative_to(root):raise RuntimeError(f'Unexpected import root: {p}')
  files[str(p.relative_to(root))]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'module':name}
data={'scope':'Post-run inventory of the unchanged reviewed dependency checkout; not a retroactive import attestation. Candidate drivers independently replayed. Later drivers must enforce their import root at run time.','qhat_root':str(root),'qhat_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'qhat_status':subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True),'python':sys.version,'packages':{name:importlib.metadata.version(name) for name in ('numpy','scipy','numba','openfermion')},'files':files}
print(json.dumps(data,indent=2))
