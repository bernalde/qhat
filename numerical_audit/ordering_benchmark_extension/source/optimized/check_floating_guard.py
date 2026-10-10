"""Reproduce the pre-existing float-delta consistency guard on F2."""
from pathlib import Path
import json,sys,inspect,time
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
import optimize as opt
import frozen_candidate as frozen
n=c.normalized_case('f2_coeff1p00');bs=[{w:n['total'][w]} for w in sorted(n['total'],key=lambda w:(c.masks(w)[0],w))];xs=[c.masks(next(iter(b)))[0] for b in bs]
pairs,rho=opt.pair_setup(bs,xs,(1<<n['design']['electrons'])-1,n['denominator'],True)
# Instrument only the old guard's diagnostic message; the original algorithm,
# float operation order, candidate tests, and guard predicate stay unchanged.
src=inspect.getsource(frozen.optimize)
needle="require(abs(actual-best)<1e-10*max(value,actual,1e-20),'Insertion delta disagrees with full objective')"
replacement="require(abs(actual-best)<1e-10*max(value,actual,1e-20),json.dumps({'pass':it,'move':best_move,'old_J':value,'incremental_J':best,'exact_amplitude_J':actual,'absolute_difference':abs(actual-best),'guard_limit':1e-10*max(value,actual,1e-20)}))"
c.require(src.count(needle)==1,'Frozen source guard changed');g=dict(frozen.__dict__);exec(src.replace(needle,replacement),g);began=time.monotonic()
try:
    g['optimize'](list(range(len(bs))),pairs,rho,n['denominator']);result={'old_guard_failed':False}
except RuntimeError as e:
    result={'old_guard_failed':True,'diagnostic':json.loads(str(e))}
result.update(case='f2_coeff1p00',seconds=time.monotonic()-began,changed_algorithm=False,instrumentation='Only replaces the frozen guard diagnostic string with local scalar values.',source_sha256={'frozen_candidate.py':c.sha(HERE/'frozen_candidate.py'),'check_floating_guard.py':c.sha(__file__)})
c.require(result['old_guard_failed'],'Expected frozen guard failure did not reproduce')
c.write_atomic(HERE/'floating_guard_diagnostic.json',result);print(json.dumps(result,indent=2),flush=True)
