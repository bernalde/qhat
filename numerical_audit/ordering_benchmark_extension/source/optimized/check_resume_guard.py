"""Verify checkpoint recovery rejects stale inputs and source under python -O."""
from pathlib import Path
import json
import checked_prepare as checked
import common as c
import run
HERE=Path(__file__).resolve().parent;case='CH2O_s1.00';records={p.name:c.sha(p) for p in (HERE/'checkpoints').glob(case+'*')};results={}
original_load=c.load_case;original_sources=run.source_hashes

def changed_input(name):
    d,s=original_load(name);d['extension_input_provenance']['design_sha256']='0'*64;return d,s
c.load_case=changed_input
try:
    checked.prepare_checked(case,adopt=True)
except RuntimeError as e:results['changed_input_rejected']=str(e)
else:raise RuntimeError('Accepted changed input')
finally:c.load_case=original_load

def changed_sources():
    d=original_sources();d['optimized/optimize.py']='0'*64;return d
run.source_hashes=changed_sources
try:
    checked.prepare_checked(case,adopt=True)
except RuntimeError as e:results['changed_source_rejected']=str(e)
else:raise RuntimeError('Accepted changed source')
finally:run.source_hashes=original_sources
c.require(records=={p.name:c.sha(p) for p in (HERE/'checkpoints').glob(case+'*')},'Rejection modified checkpoint files')
results['checkpoints_unchanged']=True;results['python_optimized']=not __debug__;results['source_sha256']=c.sha(__file__)
c.write_atomic(HERE/'checkpoint_guard_checks.json',results);print(json.dumps(results,indent=2))
