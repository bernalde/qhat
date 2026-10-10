#!/usr/bin/env python3
"""Destructive checks operate only on a temporary fixture copy, under python -O."""
import resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
import importlib.util,json,shutil,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
checks=[]
with tempfile.TemporaryDirectory(prefix='strang-control-guards-') as directory:
 d=Path(directory)
 for name in ('run.py','independent_core.py','PROTOCOL.md','ALL_PLANS_FROZEN.json'):shutil.copyfile(root/name,d/name)
 for name in ('designs','inputs'):(d/name).mkdir()
 case='CH2O_s1.00';shutil.copyfile(root/'designs'/f'{case}.json',d/'designs'/f'{case}.json');shutil.copyfile(root/'inputs'/f'{case}.npz',d/'inputs'/f'{case}.npz')
 spec=importlib.util.spec_from_file_location('temporary_control',d/'run.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 module.verified_design(case);checks.append({'probe':'valid frozen fixture','status':'accepted'})
 for name in ('designs/'+case+'.json','inputs/'+case+'.npz','PROTOCOL.md'):
  p=d/name;original=p.read_bytes();p.write_bytes(original+b' ')
  try:module.verified_design(case)
  except ValueError as e:checks.append({'probe':'tamper '+name,'status':'rejected','message':str(e)})
  else:raise RuntimeError('Tampered fixture accepted '+name)
  p.write_bytes(original)
 module.verified_design(case)
(root/'guard_checks.json').write_text(json.dumps({'status':'passed','python_optimized':not __debug__,'source_sha256':module.sha(root/'check_guards.py'),'checks':checks},indent=2)+'\n')
print(json.dumps(checks))
