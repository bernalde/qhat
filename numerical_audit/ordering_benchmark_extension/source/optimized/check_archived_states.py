"""Compare r100 repetitions with their archived pre-extension state vectors."""
from pathlib import Path
import io,json,sys,zipfile
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c
checks=[]
for case in c.CASES:
    methods=('class_refined_response','class_refined_J0','class_pauli_JD') if case in c.DEV_CASES else ('class_refined_response',)
    for method in methods:
        new=HERE/'states'/case/f'{method}.b100.npz'
        if case in c.PROS_CASES:
            member=f'results/{case}.states.npz'
            with zipfile.ZipFile(c.HELDOUT) as z:blob=z.read(member)
            source=io.BytesIO(blob);key=method;source_hash=c.sha(c.HELDOUT)
        else:
            prefix=c.MECH/'diagonal_refinement/results' if method=='class_refined_response' else c.MECH/'objective_ablation/results'/method
            source=prefix/f'{case}.states.npz';key=method+'__r100__'+('development' if method=='class_refined_response' else 'ablation');source_hash=c.sha(source)
        with np.load(source) as f:old=f[key]
        with np.load(new) as f:actual=f['actual']
        residual=float(np.linalg.norm(actual-old));c.require(residual<1e-12,'Archived state replay differs beyond tolerance')
        checks.append({'case':case,'ordering':method,'source_sha256':source_hash,'new_state_sha256':c.sha(new),'bitwise_equal':bool(np.array_equal(actual,old)),'norm_difference':residual})
c.require(len(checks)==30,'Expected30 archived-state comparisons')
out={'count':len(checks),'bitwise_equal_count':sum(r['bitwise_equal'] for r in checks),'maximum_norm_difference':max(r['norm_difference'] for r in checks),'checks':checks,'source_sha256':c.sha(__file__)};c.write_atomic(HERE/'archived_state_checks.json',out);print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
